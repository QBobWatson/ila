"""
Process LaTeX in html files.

The output is...
"""

import argparse
import os
import re
import sys
from base64 import b64encode
from contextlib import suppress
from hashlib import md5
from io import StringIO
from shutil import move, copy
from subprocess import Popen, PIPE
from tempfile import TemporaryDirectory

import cairo
from lxml import html

import simpletransform
from tounicode import add_unicode_codepoints

import gi
gi.require_version('Poppler', '0.18')
from gi.repository import Poppler  # noqa: E402


BASE = os.path.dirname(__file__)
TOUNICODE = os.path.join(BASE, 'tounicode.py')

PID = os.getpid()
FONTFORGE = 'fontforge'


def log(text):
    """Log a line of output."""
    print(f"[{PID:6d}] {text}")


# Snippet to tell fontforge to delete some empty lists.
# Otherwise the Webkit CFF sanitizer balks.
# Also, FF seems to incorrectly save default values for some entries.
FIX_PRIVATE_TABLE = '''
  if(GetPrivateEntry("OtherBlues") == "[]")
     ClearPrivateEntry("OtherBlues")
  endif
  if(GetPrivateEntry("FamilyBlues") == "[]")
     ClearPrivateEntry("FamilyBlues")
  endif
  if(GetPrivateEntry("FamilyOtherBlues") == "[]")
     ClearPrivateEntry("FamilyOtherBlues")
  endif
  if(GetPrivateEntry("BlueShift") == "")
     ChangePrivateEntry("BlueShift", "7")
  endif
  if(GetPrivateEntry("BlueScale") == "")
     ChangePrivateEntry("BlueScale", ".039625")
  endif
  if(GetPrivateEntry("BlueFuzz") == "")
     ChangePrivateEntry("BlueFuzz", "1")
  endif
'''

LATEX_PREAMBLE = r'''
\documentclass[12pt,reqno]{amsart}
\usepackage[margin=0pt]{geometry}
\usepackage[charter,sfscaled,ttscaled,cal=cmcal]{mathdesign}
\renewcommand{\sfdefault}{phv}
\usepackage{textcomp}

\newwrite\boxsize
\immediate\openout\boxsize=boxsize.txt
\def\writesize#1{\write\boxsize{#1}}
\newsavebox\measurebox

\newlength\emlength

\def\postag#1{\tag*{\phantom{#1}\pdfsavepos\write\boxsize{tag:{#1},\the\pdflastypos}}}

\pagestyle{empty}

\usepackage{graphicx}
\graphicspath{{figure-images/}{.}}
'''

LATEX_BEGIN = r'''
\begin{document}%
\topskip=0pt%
\parindent=0pt%
\parskip=0pt%
\thispagestyle{empty}%
\emlength=1em\writesize{fontsize:\the\emlength}%
'''

LATEX_NEWPAGE = r'\newpage\topskip=0pt%' + '\n'

LATEX_INLINE = r'''%
\sbox{{\measurebox}}{{%
${code}$%
}}%
\vbox to 0pt{{\vss\usebox\measurebox}}%
\writesize{{inline:{{\the\wd\measurebox}}{{\the\ht\measurebox}}{{\the\dp\measurebox}}}}%
'''

LATEX_CODE_INLINE = r'''%
\sbox{{\measurebox}}{{%
{code}%
}}%
\vbox to 0pt{{\vss\usebox\measurebox}}%
\writesize{{inline:{{\the\wd\measurebox}}{{\the\ht\measurebox}}{{\the\dp\measurebox}}}}%
'''

# tounicode.py calculates the extents for displayed equations
# html is 675px ~ 7in wide
LATEX_DISPLAY = r'''%
\pdfsavepos\write\boxsize{{prepage:\the\pdflastypos}}
\begin{{minipage}}{{7in}}%
{code}%
\end{{minipage}}%
\writesize{{display:}}%
'''

PRETEX_STYLE = '''
.pretex-bind {
  display: inline-block;
}
.pretex-inline {
  display: inline-block;
}
.pretex-inline span {
  display: inline-block;
}
.pretex-inline span:last-child {
  position: relative;
}
.pretex-inline span:last-child svg.pretex {
  position: absolute;
  bottom:   0;
  height:   1em;
}
svg.pretex {
  display:      inline-block;
  overflow:     visible;
  font-variant: normal;
  font-weight:  normal;
  font-style:   normal;
}
.pretex-display {
  text-align:  center;
  margin:      1em 0;
  padding:     0;
  text-indent: 0;
  text-transform: none;
  position: relative;
}
.pretex-display svg.pretex {
  /* hack to adjust spacing */
  vertical-align: middle;
}
.pretex-display .tag {
  position:    absolute;
  right:       0;
  top:         0;
}
.pretex-display .tag > span {
  display: inline-block;
}
'''

# This is where processed images end up under build/
FIGURE_IMG_DIR = 'figure-images'


def check_proc(proc, msg='', stdin=None):
    """Run a process and die verbosely on error."""
    if stdin is not None:
        stdin = stdin.encode('ascii')
    out, err = proc.communicate(input=stdin)
    if proc.returncode != 0:
        print(msg)
        print("stdout:")
        print(out.decode())
        print("stderr:")
        print(err.decode())
        sys.exit(1)
    return out


def css_to_dict(css_str):
    """Parse a simple css string to a dict."""
    # Won't handle complicated things like semicolons in strings.
    ret = {}
    for line in css_str.split(';'):
        idx = line.find(':')
        if idx == -1:
            continue
        key = line[:idx].strip()
        val = line[idx+1:].strip()
        if val[0] == "'" or val[0] == '"':
            val = val[1:-1]
        else:
            # Assume whitespace is unimportant in unquoted values
            val = val.replace(' ', '')
        ret[key] = val
    return ret


def dict_to_css(css):
    """Convert a dict to a css string."""
    items = []
    for key, val in css.items():
        if val.find(' ') != -1:
            val = "'" + val + "'"
        items.append(key + ':' + val)
    return ';'.join(items)


def add_class(text, cls):
    """Add a CSS class to a space-separated list."""
    if text is None:
        return cls
    text = text.strip()
    if not text:
        return cls
    return ' '.join(re.split(r'\s+', text) + [cls])


def smart_float(num, decimals=5):
    """Format a float without zero-padding on the right."""
    return format(num, "." + str(decimals) + 'f').rstrip('0').rstrip('.')


# Encoding an md5 digest in base64 instead of hex reduces length from 32 to 20
def b64_hash(text):
    """Encode `text` as a base64 string."""
    if not isinstance(text, bytes):
        text = text.encode()
    return b64encode(md5(text).digest()[:15], b'-_').decode('ascii')


class CSSClasses:
    """
    Repository for css classes used as abbreviations.

    Some style attributes, like stroke-width, font-size, and font-family,
    require many characters to specify, despite there being few unique values.
    It saves space to use css classes with short names for these.  This class
    acts as a repository for such css classes.
    """

    ALPHABET = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

    def __init__(self):
        """Initialize."""
        self.class_names = {}
        self.class_vals = {}

    def _num_to_str(self, number):
        if 0 <= number < len(self.ALPHABET):
            return self.ALPHABET[number]
        text = ''
        while number != 0:
            number, i = divmod(number, len(self.ALPHABET))
            text = self.ALPHABET[i] + text
        return text

    def get(self, val):
        """Get an abbreviation class."""
        if val in self.class_vals:
            return self.class_vals[val]
        class_name = self._num_to_str(len(self.class_names))
        self.class_names[class_name] = val
        self.class_vals[val] = class_name
        return class_name

    def css(self, prefix):
        """Return the rendered css."""
        return ''.join(prefix + '.' + name + ' { ' + val + ' }\n'
                       for name, val in sorted(self.class_names.items()))


class HTMLDoc:
    """Stores all data needed to convert the LaTeX in an html file."""

    SVG_ATTRS = set(['viewBox', 'height', 'width', 'version'])

    DEFAULT_TEXT = css_to_dict('''
        writing-mode: horizontal-tb;
        fill:         #000000;
        fill-rule:    nonzero;
        fill-opacity: 1;
        stroke:       none;
    ''')
    DEFAULT_PATH = css_to_dict('''
        fill:              none;
        fill-rule:         nonzero;
        fill-opacity:      1;
        stroke:            #000000;
        stroke-linecap:    butt;
        stroke-linejoin:   miter;
        stroke-miterlimit: 10;
        stroke-dasharray:  none;
        stroke-opacity:    1;
    ''')

    def __init__(self, html_file, preamble, tmp_dir, args):
        """Initialize."""
        self.html_file = html_file
        self.output_file = os.path.join(
            args.build_dir,
            os.path.relpath(self.html_file, args.output_dir))
        with open(self.html_file, encoding='utf-8') as fobj:
            self.html_data = fobj.read()
        parser = html.HTMLParser(remove_comments=True)
        self.dom = html.parse(StringIO(self.html_data), parser=parser)
        self.to_replace = []
        self.preamble = preamble
        self.basename = b64_hash(self.html_data)
        self.base_dir = os.path.join(tmp_dir, self.basename)
        self.pdf_dir = os.path.join(self.base_dir, 'pdf')
        self.svg_dir = os.path.join(self.base_dir, 'svg')
        self.out_img_dir = os.path.join(tmp_dir, 'img')
        self.sfd_dir = os.path.join(tmp_dir, 'sfd')
        self.cache_dir = args.cache_dir
        img_dir = args.img_dir

        os.makedirs(self.base_dir, exist_ok=True)
        os.makedirs(self.pdf_dir, exist_ok=True)
        os.makedirs(self.svg_dir, exist_ok=True)
        link_dest = os.path.join(self.pdf_dir, 'figure-images')
        if not os.path.exists(link_dest):
            os.symlink(os.path.realpath(img_dir), link_dest,
                       target_is_directory=True)

        self.latex_file = os.path.join(self.pdf_dir, self.basename + '.tex')
        self.pdf_file = os.path.join(self.pdf_dir, self.basename + '.pdf')
        self.boxsize_file = os.path.join(self.pdf_dir, 'boxsize.txt')
        self.pages_extents = []
        self.num_pages = 0
        self.fonts = {}
        self.font_hashes = {}
        self.images = []
        self.contents = ''
        self.html_cache = os.path.join(self.cache_dir, self.basename + '.html')
        self.svg_cache = None
        self.path_classes = CSSClasses()
        self.tspan_classes = CSSClasses()
        self.poppler_doc = None

    @property
    def is_cached_svg(self):
        """Whether the svg files exist in the cache."""
        return os.path.exists(self.svg_cache)

    @property
    def is_cached_html(self):
        """Whether the modified html exists in the cache."""
        return os.path.exists(self.html_cache)

    def svg_file(self, num):
        """Return the svg filename for LaTeX element number `num`."""
        return os.path.join(self.svg_dir, f'out{num+1:03d}.svg')

    def make_latex(self):
        """Extract math from the html file, then make a LaTeX file."""
        self.to_replace = []
        pages = []
        for elt in self.dom.getiterator('script'):
            if not elt.attrib.get('type', '').startswith('text/x-latex-'):
                continue
            if not elt.text:
                continue
            code = elt.text.strip()
            typ = elt.attrib['type']
            if typ == 'text/x-latex-inline':
                pages.append(LATEX_INLINE.format(code=code))
                pages.append(LATEX_NEWPAGE)
            elif typ == 'text/x-latex-code-inline':
                pages.append(LATEX_CODE_INLINE.format(code=code))
                pages.append(LATEX_NEWPAGE)
            elif typ == 'text/x-latex-code-bare':
                # Use raw code
                pages.append(code)
                continue
            elif typ in ('text/x-latex-display', 'text/x-latex-code'):
                if code.find(r'\tag') != -1:
                    code = code.replace(r'\tag', r'\postag')
                pages.append(LATEX_DISPLAY.format(
                    code=code, pageno=len(pages)))
                pages.append(LATEX_NEWPAGE)
            self.to_replace.append(elt)
        if not pages:
            # Cache the unchanged file for next time
            with open(self.html_cache, 'w', encoding='utf-8') as fobj:
                fobj.write(self.html_data)
            return False
        if pages[-1] == LATEX_NEWPAGE:
            pages = pages[:-1]
        contents = ''
        with open(self.latex_file, 'w') as fobj:
            fobj.write(LATEX_PREAMBLE)
            fobj.write(self.preamble)
            fobj.write(LATEX_BEGIN)
            fobj.write(''.join(pages))
            fobj.write(r'\end{document}')
            contents += LATEX_PREAMBLE
            contents += self.preamble
            contents += LATEX_BEGIN
            contents += ''.join(pages)
            contents += r'\end{document}'
        # Now we know the hash file name
        self.contents_hash = b64_hash(contents)
        self.svg_cache = os.path.join(self.cache_dir,
                                      self.contents_hash + '.xml')
        self.contents = contents
        return True

    def latex(self):  # noqa: D402
        """Compile the file generated by self.make_latex()."""
        proc = Popen(['pdflatex', '-interaction=nonstopmode',
                      '\\input{' + os.path.basename(self.latex_file) + '}'],
                     cwd=self.pdf_dir, stdout=PIPE, stderr=PIPE)
        check_proc(proc, f'Failed to compile LaTeX in {self.html_file}\n'
                   + 'Contents of .tex file:\n'
                   + self.contents)

    def add_unicode_codepoints(self):
        """Add unicode codepoints to the pdf file."""
        add_unicode_codepoints(self.pdf_file, self.sfd_dir)

    def measure_extents(self, pageno):
        """Render a pdf page and return its extents."""
        if self.poppler_doc is None:
            self.poppler_doc = Poppler.Document.new_from_file(
                'file://' + os.path.abspath(self.pdf_file))
        page = self.poppler_doc.get_page(pageno)
        width, height = page.get_size()
        surf = cairo.RecordingSurface(cairo.Content.COLOR_ALPHA,
                                      cairo.Rectangle(0, 0, width, height))
        ctx = cairo.Context(surf)
        page.render_for_printing(ctx)
        return surf.ink_extents()

    def read_extents(self):
        """Parse boxsize.txt and populate size data."""
        self.pages_extents = []
        this_prepage = 0
        this_tags = []
        fontsize = 12
        with open(self.boxsize_file) as fobj:
            for line in fobj.readlines():
                if line.startswith('fontsize:'):
                    # Should be the first line; ends in "pt\n"
                    # Convert to big (usual) points
                    fontsize = float(line[len('fontsize:'):-3]) * 800/803
                    continue
                if line.startswith('prepage:'):
                    # Position of the top of the page (relative to the bottom)
                    this_prepage = float(line[len('prepage:'):])
                    continue
                if line.startswith('tag:'):
                    line = line[len('tag:'):]
                    match = re.match(r'{(.*)},(.*)', line)
                    if not match:
                        continue
                    contents, pos = match.groups()
                    # Position of a tag on the page
                    pos = this_prepage - float(pos)
                    # This is in in sp = 1/65536 pt
                    pos /= 65536
                    pos *= 800/803
                    this_tags.append((contents, pos))
                    continue
                match = re.match(
                    r'inline:{(.*)pt}{(.*)pt}{(.*)pt}', line)
                if match:
                    typ = 'inline'
                    width, height, depth = [float(x) for x in match.groups()]
                    # Convert to big point (everyone else's points): *= 800/803
                    width *= 800/803
                    height *= 800/803
                    depth *= 800/803
                    left = 0
                    top = -height
                elif line.startswith('display:'):
                    typ = 'display'
                    left, top, width, height \
                        = self.measure_extents(len(self.pages_extents))
                    depth = 0
                    this_tags = [(c, y - top) for c, y in this_tags]
                else:
                    continue
                # In Inkscape, 96 user units (or "px") is one inch, which is 72
                # pt.  The "width", "height", and "depth" are used to specify
                # the viewBox, which is in user units.
                page_extents = {
                    "width"  : width  * 96/72,  # noqa: E221 E203
                    "height" : height * 96/72,  # noqa: E221 E203
                    "left"   : left   * 96/72,  # noqa: E221 E203
                    "top"    : top    * 96/72,  # noqa: E221 E203
                    "depth"  : depth  * 96/72,  # noqa: E221 E203
                    "tags"   : this_tags,       # noqa: E221 E203
                    # These are used for the "width", "height", and
                    # "vertical-align" properties, which are relative to the
                    # current font size.
                    "fontsize" : fontsize,           # noqa: E221 E203
                    "widthem"  : width  / fontsize,  # noqa: E221 E203
                    "heightem" : height / fontsize,  # noqa: E221 E203
                    "depthem"  : depth  / fontsize,  # noqa: E221 E203
                    "display"  : typ == 'display',   # noqa: E221 E203
                }
                self.pages_extents.append(page_extents)
                this_tags = []
                this_prepage = 0
        self.num_pages = len(self.pages_extents)
        self.DEFAULT_TEXT['font-size'] = f"{fontsize}px"

    def inkscape_script(self):
        """Generate inkscape commands necessary to convert pdf to svg."""
        script = []
        for page_num in range(self.num_pages):
            script.append(
                'export-plain-svg:true;'
                f'open-page:{page_num+1};'  # this has to come before file-open
                f'file-open:{self.pdf_file};'
                f'export-filename:{self.svg_file(page_num)};'
                'export-do;'
                'file-close'
            )
        return '\n'.join(script) + '\n'

    def add_font(self, name, fname):
        """Add a font to the generated html file."""
        with open(fname, 'rb') as fobj:
            self.fonts[name] = fobj.read()
        self.font_hashes[name] = 'f'+b64_hash(self.fonts[name])

    def write_cache(self, style, fonts, svgs):
        """Cache the computed data in an xml file."""
        cache = html.Element('cache')
        elt = html.Element('style', id='pretex-style')
        elt.text = style
        cache.append(elt)
        elt = html.Element('style', id='pretex-fonts')
        elt.text = fonts
        cache.append(elt)
        for svg in svgs:
            svg.tail = ''
            cache.append(svg)
        with open(self.svg_cache, 'wb') as fobj:
            fobj.write(html.tostring(cache))

    def _replace_elt(self, elt, svg):
        """
        Replace an element with an svg.

        Use a binding wrapper if necessary.
        """
        # The code below is to prevent a line break occurring between an
        # equation and an adjacent piece of text, like "ith" or "(and f(x))".
        if svg.attrib['class'] == 'pretex-inline':
            head_text = ''
            tail_text = elt.tail
            need_wrap = False
            if elt.tail and not elt.tail[0].isspace():
                match = re.match(r'(\S+)(.*)', elt.tail)
                svg.tail, tail_text = match.groups()
                need_wrap = True
            # Figure out what text came right before this element
            prev = elt.getprevious()
            parent = elt.getparent()
            if prev is not None and prev.tail and not prev.tail[-1].isspace():
                match = re.match(r'(.*\s)?(\S+)', prev.tail)
                unbound, head_text = match.groups()
                prev.tail = unbound or ''
                need_wrap = True
            elif (prev is None and
                  parent.text and
                  not parent.text[-1].isspace()):
                match = re.match(r'(.*\s)?(\S+)', parent.text)
                unbound, head_text = match.groups()
                parent.text = unbound or ''
                need_wrap = True
            if need_wrap:
                # Wrap in a binding span
                wrapper = html.Element('span', {'class': 'pretex-bind'})
                wrapper.text = head_text
                wrapper.append(svg)
                wrapper.tail = tail_text
                svg = wrapper
            else:
                svg.tail = tail_text
        else:
            svg.tail = elt.tail
        parent = elt.getparent()
        parent.replace(elt, svg)

    def _rewrite_common(self, style, fonts):
        """Code common to use_cached_svg() and write_html()."""
        root = self.dom.getroot()
        try:
            root.get_element_by_id('pretex-style').text = style
        except KeyError:
            pass
        try:
            root.get_element_by_id('pretex-fonts').text = fonts
        except KeyError:
            pass
        for elt in self.dom.getiterator('script'):
            if elt.attrib.get('type', '').startswith('text/x-latex-code-bare'):
                elt.drop_tree()
        # Add base class name to all immediate children of <body>.  These are
        # the elements that are imported by the knowl mechanism.
        base_class = 'C' + self.contents_hash
        for elt in root.find('body'):
            elt.attrib['class'] = add_class(elt.attrib.get('class'),
                                            base_class)

    def use_cached_svg(self):
        """Write the cached output to the html file."""
        with open(self.svg_cache, 'rb') as fobj:
            cache = html.fromstring(fobj.read())
        style = cache[0].text
        fonts = cache[1].text
        # Replace DOM elements
        for elt in self.to_replace:
            svg = cache[2]
            self._replace_elt(elt, svg)
        self._rewrite_common(style, fonts)
        contents = html.tostring(
            self.dom, include_meta_content_type=True, encoding='utf-8')
        with open(self.output_file, 'wb') as outf:
            outf.write(contents)
        with open(self.html_cache, 'wb') as outf:
            outf.write(contents)

    def copy_over(self):
        """Write html data directly to output file."""
        copy(self.html_cache, self.output_file)

    def write_html(self):
        """Write the modified html file."""
        svgs = self.process_svgs()
        cached_elts = []
        # Replace DOM elements
        for i, elt in enumerate(self.to_replace):
            self._replace_elt(elt, svgs[i])
            cached_elts.append(svgs[i])
        style = PRETEX_STYLE
        style += r'''
svg.pretex text {{
  {}
}}
svg.pretex path {{
  {}
}}
'''.format(dict_to_css(self.DEFAULT_TEXT), dict_to_css(self.DEFAULT_PATH))
        # Add fonts
        font_style = f'\n/* pretex cache: {self.contents_hash} */\n'
        for name, data in self.fonts.items():
            name = self.font_hashes[name]
            font_style += r'''
@font-face {{
  font-family: "{name}";
  src: url(data:application/font-woff;base64,{data}) format('woff');
}}
'''.format(name=name, data=b64encode(data).decode('ascii'))
        font_style += '\n'
        # These go here so they show up in knowls too
        base_class = 'C' + self.contents_hash
        font_style += self.tspan_classes.css(
            f'.{base_class} svg.pretex tspan')
        font_style += self.path_classes.css(
            f'.{base_class} svg.pretex path')
        self._rewrite_common(style, font_style)
        content = html.tostring(
            self.dom, include_meta_content_type=True, encoding='utf-8')
        with open(self.output_file, 'wb') as outf:
            outf.write(content)
        with open(self.html_cache, 'wb') as outf:
            outf.write(content)
        self.write_cache(style, font_style, cached_elts)

    def process_svgs(self):
        """Process all generated svgs file for use in an html page."""
        svgs = []
        for page_num, page_extents in enumerate(self.pages_extents):
            with open(self.svg_file(page_num), 'rb') as fobj:
                svg = html.fromstring(fobj.read())
            # Remove extra attrs from <svg>
            for key in svg.attrib.keys():
                if key not in self.SVG_ATTRS:
                    del svg.attrib[key]
            # Auto-calculated based on height and viewBox aspect ratio
            del svg.attrib['width']
            svg.attrib['class'] = 'pretex'
            # Get rid of metadata
            metadata = svg.find('metadata')
            if metadata is not None:
                metadata.drop_tree()
            # Get rid of empty defs
            defs = svg.find('defs')
            if defs is not None and len(defs) == 0:
                defs.drop_tree()
            # Undo global page coordinate transforms
            units_in_pt = unwrap_transforms(svg)
            # Plug in actual size data
            if page_extents['display']:
                scale = 72/96 if units_in_pt else 1
                svg.attrib['viewBox'] = '{} {} {} {}'.format(
                    smart_float(scale * page_extents['left']),
                    smart_float(scale * page_extents['top']),
                    smart_float(scale * page_extents['width']),
                    smart_float(scale * page_extents['height'])
                )
                # The height is 1em.  The fonts in the pdf file are relative to
                # fontsize.
                svg.attrib['height'] = '{}em'.format(
                    smart_float(page_extents['heightem']))
            else:
                scale = 1 if units_in_pt else 96/72
                # The size of the view box doesn't matter, since the wrapper
                # and the strut take care of spacing.  Set it to a 1em square.
                svg.attrib['viewBox'] = '0 -{fs} {fs} {fs}'.format(
                    fs=smart_float(page_extents['fontsize']*(
                        1 if units_in_pt else 96/72)))
                # height is 1em; it is set in css
                del svg.attrib['height']
            # Get rid of ids
            for elt in svg.xpath('//*[@id]'):
                if not elt.xpath("ancestor::defs"):
                    del elt.attrib['id']
            # Clean up text styles
            for tspan in svg.xpath('//tspan[@style]'):
                self.process_tspan(tspan, page_extents['fontsize'])
            # Clean up path styles
            for path in svg.xpath('//path[@style]'):
                self.process_path(path)
            # Process linked images
            for img in svg.xpath('//image'):
                self.process_image(img)
            # Delete empty groups (recursively)
            todelete = svg.xpath('//g[count(*)=0]')
            while todelete:
                todelete2 = list(todelete)
                todelete = []
                for elt in todelete2:
                    parent = elt.getparent()
                    elt.drop_tree()
                    if parent.tag == 'g' and len(parent) == 0:
                        todelete.append(parent)
            if page_extents['display']:
                # Wrap displayed equations
                div = html.Element('div', {'class': 'pretex-display'})
                div.append(svg)
                svg = div
                # Add tags
                for contents, pos in page_extents['tags']:
                    tagelt = html.Element('span', {'class': 'tag'})
                    tagelt.text = '('+contents+')'
                    # This moves the tag down the calculated amount
                    htelt = html.Element('span', style='height:{}em'.format(
                        smart_float(pos / page_extents['fontsize'])))
                    tagelt.append(htelt)
                    svg.append(tagelt)
            else:
                # After much experimentation, this seems to be the most
                # reliable way to lock the origin of the svg to the baseline.
                wrapper = html.Element('span', {
                    'class': 'pretex-inline',
                    'style': 'width:{}em'.format(
                        smart_float(page_extents['widthem'])),
                })
                # make strut
                style = 'height:{}em'.format(
                    smart_float(page_extents['heightem'] +
                                page_extents['depthem']))
                if page_extents['depthem'] > 0.0:
                    style += ';vertical-align:-{}em'.format(
                        smart_float(page_extents['depthem']))
                wrapper.append(html.Element('span', style=style))
                # This last span is relatively positioned.  Its size will be
                # 0x0, so it sits right on the baseline.  The bottom of the svg
                # is then absolutely positioned to that.
                elt = html.Element('span')
                wrapper.append(elt)
                elt.append(svg)
                svg = wrapper
            svgs.append(svg)
        return svgs

    def process_tspan(self, tspan, page_font_size):
        """Simplify <tspan> tag."""
        css = css_to_dict(tspan.attrib.get('style', ''))
        # These are hard-coded into the font, but not marked as such
        css.pop('font-variant', 1)
        css.pop('font-weight', 1)
        css.pop('font-style', 1)
        css.pop('-inkscape-font-specification', 1)
        # Get rid of inherited styles
        for key in self.DEFAULT_TEXT:
            if key in css and css[key] == self.DEFAULT_TEXT[key]:
                del css[key]
        if 'writing-mode' in css and css['writing-mode'] == 'lr-tb':
            del css['writing-mode']
        # Add css class to save space
        css_val = []
        if 'font-size' in css:
            css_val.append('font-size:'+css['font-size'])
            del css['font-size']
        else:
            # Shouldn't happen
            print("WARNING: unspecified font-size in tspan")
        if 'font-family' in css:
            font_family = css['font-family'].split(',')
            if font_family and font_family[0]:
                font_family = font_family[0]
                if font_family in self.font_hashes:
                    css_val.append(
                        'font-family:'+self.font_hashes[font_family])
                else:
                    # Shouldn't happen
                    css_val.append('font-family:'+font_family)
                del css['font-family']
        else:
            # Shouldn't happen
            print("WARNING: unspecified font-family in tspan")
        tspan.attrib['style'] = dict_to_css(css)
        if not tspan.attrib['style']:
            del tspan.attrib['style']
        if css_val:
            tspan.attrib['class'] = add_class(
                tspan.attrib.get('class'),
                self.tspan_classes.get(';'.join(css_val))
            )

    def process_path(self, path):
        """Simplify <path> tag."""
        css = css_to_dict(path.attrib.get('style', ''))
        # Get rid of inherited styles
        for key in self.DEFAULT_PATH:
            if key in css and css[key] == self.DEFAULT_PATH[key]:
                del css[key]
        if css.get('fill') == '#000000':
            css['fill'] = '#000'
        # Add class to save space
        css_val = []
        if 'stroke-width' in css:
            swd = css['stroke-width']
            # Append "px" to unitless numbers
            try:
                float(swd)
            except ValueError:
                pass
            else:
                swd = swd + 'px'
            css_val.append('stroke-width:'+swd)
            del css['stroke-width']
        else:
            # The default value is 1.
            css_val.append('stroke-width:1px')
        path.attrib['style'] = dict_to_css(css)
        if not path.attrib['style']:
            del path.attrib['style']
        if css_val:
            path.attrib['class'] = add_class(
                path.attrib.get('class'),
                self.path_classes.get(';'.join(css_val))
            )

    def process_image(self, img):
        """Simplify <image> tag."""
        href = img.attrib['xlink:href']
        del img.attrib['xlink:href']
        # Inkscape has no idea where the file ended up
        fname = os.path.join(self.out_img_dir, os.path.basename(href))
        # Cache the image by a hash of its content
        with open(fname, 'rb') as fobj:
            img_hash = b64_hash(fobj.read())
        img_name = img_hash + '.png'
        self.images.append(img_name)
        img.attrib['href'] = FIGURE_IMG_DIR + '/' + img_name
        # Move to the cache directory
        move(fname, os.path.join(self.cache_dir, img_name))
        # Simplify css
        css = css_to_dict(img.get('style', ''))
        css.pop('image-rendering', 1)
        img.attrib['style'] = dict_to_css(css)
        if not img.attrib['style']:
            del img.attrib['style']


def almost_zero(num, ε=0.0001):
    """Decide if a number is close enough to zero."""
    return abs(num) < ε


def smart_round(num, decimals=8):
    """
    Round "num" to `decimals`.

    This uses the fewest decimal places possible within given precision.
    """
    # There must be a less stupid algorithm...
    if not isinstance(num, float):
        return num
    error = 1.0
    for i in range(decimals):
        error /= 10
    if num < 0:
        num *= -1
        neg = -1
    else:
        neg = 1
    for i in range(decimals):
        shift = num
        for j in range(i):
            shift *= 10
        approx1 = int(shift)
        approx2 = int(shift) + 1
        for j in range(i):
            approx1 /= 10.0
            approx2 /= 10.0
        if num - approx1 < error:
            return format(neg * approx1, "." + str(i) + "f")
        if approx2 - num < error:
            return format(neg * approx2, "." + str(i) + "f")
    return neg * num


def simplify_transforms(svg):
    """Re-format transform attributes to save characters."""
    for elt in svg.xpath('//*[@transform]'):
        mat = simpletransform.parse_transform(elt.attrib['transform'])
        # Recognize identity / translation
        if(almost_zero(mat[0][0] - 1) and
           almost_zero(mat[1][1] - 1) and
           almost_zero(mat[0][1]) and
           almost_zero(mat[1][0])):
            if almost_zero(mat[1][2]):
                if almost_zero(mat[0][2]):
                    del elt.attrib['transform']
                    continue
                elt.attrib['transform'] = 'translate({})'.format(
                    smart_round(mat[0][2]))
                continue
            elt.attrib['transform'] = 'translate({} {})'.format(
                smart_round(mat[0][2]), smart_round(mat[1][2]))
            continue
        # Recognize scale
        if(almost_zero(mat[0][1]) and
           almost_zero(mat[0][2]) and
           almost_zero(mat[1][0]) and
           almost_zero(mat[1][2])):
            if almost_zero(mat[0][0] - mat[1][1]):
                elt.attrib['transform'] = 'scale({})'.format(
                    smart_round(mat[0][0]))
                continue
            elt.attrib['transform'] = 'scale({} {})'.format(
                smart_round(mat[0][0]), smart_round(mat[1][1]))
            continue
        elt.attrib['transform'] = "matrix({},{},{},{},{},{})".format(
            smart_round(mat[0][0]), smart_round(mat[1][0]),
            smart_round(mat[0][1]), smart_round(mat[1][1]),
            smart_round(mat[0][2]), smart_round(mat[1][2]))


def unwrap_transforms(svg):
    """Undo global coordinate transformation, if there is one."""
    groups = svg.findall('g')
    if len(groups) != 1:
        return False
    group = groups[0]
    if not set(group.attrib.keys()) <= {'id', 'transform'}:
        return False
    if not group.attrib.get('transform'):
        return False
    mat = simpletransform.parse_transform(group.attrib['transform'])
    # Recognize pdf coordinate transformation
    if not (almost_zero(mat[0][0] - 4/3) and
            almost_zero(mat[1][1] + 4/3) and
            almost_zero(mat[0][1])       and  # noqa: E272
            almost_zero(mat[1][0])       and  # noqa: E272
            almost_zero(mat[0][2])):
        return False
    mat = [[1, 0, 0], [0, -1, mat[1][2]*3/4]]
    for child in group:
        child_mat = simpletransform.parse_transform(
            child.attrib.get('transform', ''), mat)
        child.attrib['transform'] = simpletransform.format_transform(child_mat)
        svg.append(child)
    group.drop_tree()
    simplify_transforms(svg)
    return True


def main():
    """Run the main routine."""
    parser = argparse.ArgumentParser(
        description='Process LaTeX in html files.')
    parser.add_argument('--preamble', default='preamble.tex', type=str,
                        help='LaTeX preamble')
    parser.add_argument('--style-path', default='', type=str,
                        help='Location of LaTeX style files')
    parser.add_argument('--cache-dir', default='pretex-cache', type=str,
                        help='Cache directory')
    parser.add_argument('--img-dir', default='figure-images', type=str,
                        help='LaTeX image include directory')
    parser.add_argument('--no-cache', action='store_true',
                        help='Ignore cache and regenerate')
    parser.add_argument('--output-dir', type=str, required=True,
                        help='HTML output directory')
    parser.add_argument('--build-dir', type=str, required=True,
                        help='Final build directory')
    parser.add_argument('htmls', type=str, nargs='+',
                        help='HTML files to process')
    args = parser.parse_args()

    with open(args.preamble) as fobj:
        preamble = fobj.read()

    if args.style_path:
        os.environ['TEXINPUTS'] = f'.:{args.style_path}:'
    os.makedirs(args.cache_dir, exist_ok=True)

    # tmpdir = "/tmp/ila"
    # if True:
    with TemporaryDirectory() as tmpdir:
        sfd_dir = os.path.join(tmpdir, 'sfd')
        os.makedirs(sfd_dir, exist_ok=True)

        html_files = [HTMLDoc(html, preamble, tmpdir, args)
                      for html in args.htmls]

        # Create pdf files
        log(f"Processing {len(html_files)} files")
        log("Extracting code and running LaTeX...")
        done = set()
        for htmlf in html_files:
            basename = os.path.basename(htmlf.html_file)
            if htmlf.is_cached_html and not args.no_cache:
                # The input html file is unchanged
                done.add(htmlf)
                htmlf.copy_over()
                continue
            if not htmlf.make_latex():
                # Nothing to TeX
                done.add(htmlf)
                htmlf.copy_over()
                continue
            if htmlf.is_cached_svg and not args.no_cache:
                log(f"Using cached svgs for {basename}")
                done.add(htmlf)
                htmlf.use_cached_svg()
                continue
            else:
                log(f"(Re)processing {basename}")
                htmlf.latex()
                htmlf.add_unicode_codepoints()
                htmlf.read_extents()
        html_files = [h for h in html_files if h not in done]
        if not html_files:
            log("Done!")
            return
        html_byhash = {h.basename: h for h in html_files}

        # Convert all fonts
        log("Converting fonts to woff format...")
        woff_dir = os.path.join(tmpdir, 'woff')
        os.makedirs(woff_dir, exist_ok=True)
        script = []
        for fname in os.listdir(sfd_dir):
            if fname[-4:] != '.sfd':
                continue
            fullpath = os.path.join(sfd_dir, fname)
            entry = ''
            entry += f'Open("{fullpath}")\n' \
                + FIX_PRIVATE_TABLE \
                + 'Generate("{}")\n'.format(
                    os.path.join(woff_dir, fname[:-4] + '.woff'))
            script.append(entry)
            # Process 1000 at a time; otherwise ff might segfault
            if len(script) == 1000:
                proc = Popen([FONTFORGE, '-lang=ff', '-script', '-'],
                             stdin=PIPE, stdout=PIPE, stderr=PIPE)
                check_proc(proc, 'Could not convert pdf fonts to woff format',
                           stdin=''.join(script))
                script = []
        proc = Popen([FONTFORGE, '-lang=ff', '-script', '-'],
                     stdin=PIPE, stdout=PIPE, stderr=PIPE)
        check_proc(proc, 'Could not convert pdf fonts to woff format',
                   stdin=''.join(script))

        # Associate the fonts with their html filse
        for fname in os.listdir(woff_dir):
            match = re.match(r'\[(.*)\](.*)\.woff', fname)
            if not match:
                continue
            hash_name, font_name = match.groups()
            # In debugging mode, all workers share a tmp directory
            with suppress(KeyError):
                html_byhash[hash_name].add_font(
                    font_name.replace('+', ' '), os.path.join(woff_dir, fname))

        # Convert all pages of all pdf files to svg files
        log("Generating svg files...")
        # inkscape exports images to the current directory
        img_dir = os.path.join(tmpdir, 'img')
        os.makedirs(img_dir, exist_ok=True)
        script = ''.join(html.inkscape_script() for html in html_files)
        proc = Popen(['inkscape', '--shell'],
                     stdout=PIPE, stderr=PIPE, stdin=PIPE,
                     cwd=img_dir)
        check_proc(proc, "SVG conversion failed", script)

        # Process svg files and write html
        log("Writing html files...")
        for htmlf in html_files:
            htmlf.write_html()
        log("Done!")


if __name__ == "__main__":
    main()
