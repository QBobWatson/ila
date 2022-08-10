"""
Add ToUnicode tables to all embedded fonts.

This module contains a function that uses fontforge and pdfrw to add ToUnicode
tables to all embedded fonts in the pdf files specified on the command line.
Glyphs with no reasonable guess get a "miscellaneous symbol" codepoint.
"""

import os
import sys
import unicodedata

import fontforge
from pdfrw import PdfReader, PdfWriter, PdfDict

from aglfn import GLYPHS, GLYPHS_BYCP
from pdf_enc import ENCODINGS


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i+n]


def generate_tounicode(ff_font, pdffont):
    """Generate a ToUnicode CMAP table for one font."""
    cmap_entries = []
    used_codepoints = set([GLYPHS[glyphname].codepoint
                           for glyphname in ff_font
                           if glyphname in GLYPHS])
    pdf_encoding = ENCODINGS['NoEncoding']()
    if hasattr(pdffont, 'Encoding') and pdffont.Encoding:
        # fontforge uses the raw postscript font; it won't take into account
        # the pdf font's encoding override
        encoding = pdffont.Encoding
        if hasattr(encoding, 'BaseEncoding') and encoding.BaseEncoding:
            pdf_encoding = ENCODINGS[encoding.BaseEncoding]()
        else:
            pdf_encoding = ENCODINGS['StandardEncoding']()
        if hasattr(encoding, 'Differences') and encoding.Differences:
            pdf_encoding.modify(encoding.Differences)
    # Latin script starts here
    last_unknown = 0x0021
    for glyphname in ff_font:
        glyph = ff_font[glyphname]
        if glyphname in GLYPHS:
            codepoint = GLYPHS[glyphname].codepoint
        else:
            # Use the next nice symbol
            while (last_unknown in used_codepoints or
                   last_unknown not in GLYPHS_BYCP or
                   unicodedata.category(chr(last_unknown))
                   not in ('Lu', 'Ll')):
                last_unknown += 1
            codepoint = last_unknown
            glyph.glyphname = GLYPHS_BYCP[codepoint].name
            last_unknown += 1
        glyph.unicode = codepoint
        used_codepoints.add(codepoint)
        if glyph.width == 0:
            # Safari won't display zero-width glyphs -- and the "not equal"
            # sign uses one... anyway, "1" is 1/1000 em.
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
        out += f'{len(chunk)} beginbfchar\n'
        for encoding, codepoint in chunk:
            out += f'<{encoding:02X}> <{codepoint:04X}>\n'
        out += 'endbfchar\n'
    out += ('endcmap\n'
            'end\n')
    return out


def add_unicode_codepoints(pdf, outdir):
    """
    Add ToUnicode tables to all embedded fonts in pdf files.

    This modifies the pdf file to update the embedded fonts, and saves the
    fonts in sfd format to `outdir`.
    """
    # print(f"Adding ToUnicode tables to PDF file {pdf}")
    fontnum = 0
    doc = PdfReader(pdf)
    fonts = [o for o in doc.indirect_objects.values()
             if hasattr(o, 'Type') and o.Type == '/Font']
    fonts = {font.FontDescriptor.FontName[1:]: font
             for font in fonts if font.FontDescriptor is not None}
    embedded_fonts = fontforge.fontsInFile(pdf)
    for fontname in embedded_fonts:
        if fontname not in fonts:
            print(f"WARNING: font {fontname} not found in pdf file",
                  file=sys.stderr)
            continue
        # print(f"Adding ToUnicode table to font {fontname}")
        # Low-level suppress stderr to shut up fontforge.  See
        #    https://stackoverflow.com/questions/68716139/
        try:
            devnull = open(os.devnull, 'w')
            orig_stderr = os.dup(2)
            os.dup2(devnull.fileno(), 2)
            ff_font = fontforge.open(f'{pdf}({fontname})')
        finally:
            os.dup2(orig_stderr, 2)
            os.close(orig_stderr)
            devnull.close()
        fonts[fontname].ToUnicode = PdfDict()
        fonts[fontname].ToUnicode.stream = generate_tounicode(
            ff_font, fonts[fontname])
        # Need to save the modified font because fontforge won't read
        # ToUnicode when it converts to woff later.
        ff_font.fontname = f'pretex{fontnum:06d}'
        ff_font.save(os.path.join(
            outdir, f'[{os.path.basename(pdf)[:-4]}]{fontname}.sfd'))
        fontnum += 1
    PdfWriter(pdf, trailer=doc).write()


if __name__ == "__main__":
    add_unicode_codepoints(sys.argv[1], sys.argv[2])
