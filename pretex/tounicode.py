#!/usr/bin/python2

# This script uses fontforge and pdfrw to add ToUnicode tables to all embedded
# fonts in the pdf files specified on the command line.  Glyphs with no
# reasonable guess get a "miscellaneous symbol" codepoint.
#
# Since we're already running python2, we also measure the extents of the
# rendered content in the output pdf.

import argparse
import os
import sys
import unicodedata

import platform
if platform.system() == 'Darwin':
    sys.path.append('/Applications/FontForge.app/Contents/Resources/opt/local/lib/python2.7/site-packages')

import cairo
import fontforge
import poppler
from pdfrw import PdfReader, PdfWriter, PdfDict

from aglfn import GLYPHS, GLYPHS_BYCP
from pdf_enc import ENCODINGS


def chunks(l, n):
    "Yield successive n-sized chunks from l."
    for i in range(0, len(l), n):
        yield l[i:i+n]

def generate_tounicode(font, pdffont):
    "Generate a ToUnicode CMAP table for the font."
    cmap_entries = []
    used_codepoints = set([GLYPHS[glyphname].codepoint
                           for glyphname in font
                           if glyphname in GLYPHS])
    pdf_encoding = ENCODINGS['NoEncoding']()
    if hasattr(pdffont, 'Encoding') and pdffont.Encoding:
        # fontforge uses the raw postscript font; it won't take into account the
        # pdf font's encoding override
        encoding = pdffont.Encoding
        if hasattr(encoding, 'BaseEncoding') and encoding.BaseEncoding:
            pdf_encoding = ENCODINGS[encoding.BaseEncoding]()
        else:
            pdf_encoding = ENCODINGS['StandardEncoding']()
        if hasattr(encoding, 'Differences') and encoding.Differences:
            pdf_encoding.modify(encoding.Differences)
    # Latin script starts here
    last_unknown = 0x0021
    for glyphname in font:
        glyph = font[glyphname]
        if glyphname in GLYPHS:
            codepoint = GLYPHS[glyphname].codepoint
        else:
            # Use the next nice symbol
            while (last_unknown in used_codepoints or
                   last_unknown not in GLYPHS_BYCP or
                   unicodedata.category(unichr(last_unknown))
                   not in ('Lu', 'Ll')):
                last_unknown += 1
            codepoint = last_unknown
            glyph.glyphname = GLYPHS_BYCP[codepoint].name
            last_unknown += 1
        glyph.unicode = codepoint
        used_codepoints.add(codepoint)
        if glyph.width == 0:
            # Safari won't display zero-width glyphs -- and the "not equal" sign
            # uses one... anyway, "1" is 1/1000 em.
            glyph.width = 1
        if glyphname in pdf_encoding:
            cmap_entries.append((pdf_encoding[glyphname], codepoint))
        else:
            cmap_entries.append((glyph.encoding, codepoint))

    out = ''
    out += ('12 dict begin\n'
            'begincmap\n'
            '1 begincodespacerange\n'
            '<00> <FF>\n'
            'endcodespacerange\n')
    for chunk in chunks(cmap_entries, 100):
        out += '{} beginbfchar\n'.format(len(chunk))
        for encoding, codepoint in chunk:
            out += '<{:02X}> <{:04X}>\n'.format(encoding, codepoint)
        out += 'endbfchar\n'
    out += ('endcmap\n'
            'end\n')
    return out

def main():
    parser = argparse.ArgumentParser(
        description='Add ToUnicode tables to PDF files.')
    parser.add_argument('--outdir', default='tmp/sfd', type=str,
                        help='Output .sfd files to this directory')
    parser.add_argument('pdfs', type=str, nargs='+',
                        help='PDF files to process')
    args = parser.parse_args()

    fontnum = 0
    for pdf in args.pdfs:
        print("Adding ToUnicode tables to PDF file {}".format(pdf))
        with open(pdf, 'rb') as fobj:
            pdfdata = fobj.read()
        doc = PdfReader(fdata=pdfdata)
        doc.read_all()
        fonts = [o for o in doc.indirect_objects.values()
                 if hasattr(o, 'Type') and o.Type == '/Font']
        fonts = {font.FontDescriptor.FontName[1:] : font
                 for font in fonts if font.FontDescriptor is not None}
        embedded_fonts = fontforge.fontsInFile(pdf)
        for fontname in embedded_fonts:
            if fontname not in fonts:
                print("WARNING: font {} not found in pdf file"
                      .format(fontname))
                continue
            print("Adding ToUnicode table to font {}".format(fontname))
            font = fontforge.open('{}({})'.format(pdf, fontname))
            fonts[fontname].ToUnicode = PdfDict()
            fonts[fontname].ToUnicode.stream = generate_tounicode(
                font, fonts[fontname])
            # Need to save the modified font because fontforge won't read
            # ToUnicode when it converts to woff later.
            font.fontname = 'pretex{:06d}'.format(fontnum)
            font.save(os.path.join(
                args.outdir, '[{}]{}.sfd'.format(
                    os.path.basename(pdf)[:-4], fontname)))
            fontnum += 1
        PdfWriter(pdf, trailer=doc).write()

        # Measure extents for displayed equations
        pdfpath = os.path.realpath(os.path.dirname(pdf))
        doc = poppler.document_new_from_file(
            'file://{}'.format(os.path.realpath(pdf)), None)
        boxsize = os.path.join(pdfpath, 'boxsize.txt')
        with open(boxsize) as fobj:
            lines = fobj.readlines()
        with open(boxsize, 'w') as fobj:
            pageno = 0
            for line in lines:
                if not (line.startswith('inline:') or
                        line.startswith('display:')):
                    fobj.write(line)
                    continue
                pageno += 1
                if not line.startswith('display:'):
                    fobj.write(line)
                    continue
                page = doc.get_page(pageno-1)
                width, height = page.get_size()
                surf = cairo.RecordingSurface(
                    cairo.Content.COLOR_ALPHA,
                    cairo.Rectangle(0, 0, width, height))
                ctx = cairo.Context(surf)
                page.render_for_printing(ctx)
                x, y, w, h = surf.ink_extents()
                fobj.write(line.strip() + '{},{},{},{}\n'
                           .format(x, y, w, h))

if __name__ == '__main__':
    main()
