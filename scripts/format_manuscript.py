"""Format the Pandoc manuscript and create its covering letter from the source text."""
from pathlib import Path

from docx import Document
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

MANUSCRIPT = Path(__file__).resolve().parents[1] / "manuscript"


def set_element(parent, name, **attributes):
    element = parent.find(qn("w:" + name))
    if element is None:
        element = OxmlElement("w:" + name)
        parent.append(element)
    for key, value in attributes.items():
        element.set(qn("w:" + key), str(value))
    return element


def page_setup(document):
    section = document.sections[0]
    section.page_height, section.page_width = Inches(11), Inches(8.5)
    section.top_margin, section.bottom_margin = Inches(0.85), Inches(0.83)
    section.left_margin, section.right_margin = Inches(1.05), Inches(0.95)
    section.footer_distance = Inches(0.38)
    document.core_properties.author = "Mohammad Javad Eslami"
    return section


def format_table(table, widths):
    table.autofit = False
    borders = set_element(table._tbl.tblPr, "tblBorders")
    for side in ("top", "bottom", "left", "right", "insideH", "insideV"):
        set_element(borders, side, val="single", sz=4, color="D9D9D9")
    for col, width in zip(table.columns, widths):
        col.width = Inches(width)
    for index, row in enumerate(table.rows):
        set_element(row._tr.get_or_add_trPr(), "cantSplit")
        if index == 0:
            set_element(row._tr.get_or_add_trPr(), "tblHeader", val="true")
        for cell, width in zip(row.cells, widths):
            cell.width = Inches(width)
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            if index == 0:
                set_element(cell._tc.get_or_add_tcPr(), "shd", fill="E9E9E9")
            for p in cell.paragraphs:
                fmt = p.paragraph_format
                fmt.keep_with_next = index < len(table.rows) - 1
                fmt.space_before, fmt.space_after = Pt(3), Pt(3)
                fmt.line_spacing = 1.12
                for run in p.runs:
                    run.bold = index == 0


def format_manuscript():
    path = MANUSCRIPT / "Manuscript.docx"
    document = Document(path)
    section = page_setup(document)
    for style in document.styles:
        if style.type == WD_STYLE_TYPE.PARAGRAPH:
            style.font.name = "Times New Roman"
            style.font.color.rgb = RGBColor(0, 0, 0)
    for name in ("Normal", "First Paragraph", "Body Text"):
        style = document.styles[name]
        style.font.size = Pt(11)
        fmt = style.paragraph_format
        fmt.space_before, fmt.space_after = Pt(0), Pt(5)
        fmt.line_spacing = 2
        fmt.widow_control = True
    for name, size in (("Title", 17), ("Heading 1", 14), ("Heading 2", 12)):
        style = document.styles[name]
        style.font.size, style.font.bold = Pt(size), True
        fmt = style.paragraph_format
        fmt.space_before, fmt.space_after = Pt(8), Pt(4)
        fmt.line_spacing, fmt.keep_with_next = 1.15, True
    title = document.styles["Title"]._element.get_or_add_pPr()
    border = title.find(qn("w:pBdr"))
    if border is not None:
        title.remove(border)

    line_numbers = set_element(section._sectPr, "lnNumType", countBy=1, restart="continuous", distance=170)
    section._sectPr.insert_element_before(line_numbers, "w:cols")
    footer = section.footer.paragraphs[0]
    footer.clear()
    set_element(footer._p.get_or_add_pPr(), "suppressLineNumbers")
    footer.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    set_element(footer._p, "fldSimple", instr="PAGE")

    in_front_matter = True
    in_references = False
    for p in document.paragraphs:
        fmt = p.paragraph_format
        if p.style.name == "Title":
            fmt.alignment = WD_ALIGN_PARAGRAPH.LEFT
            for run in p.runs:
                run.text = run.text.replace("500 nm", "500\u00a0nm")
                run.font.name = "Times New Roman"
        elif p.style.name in ("Heading 2", "Heading 3"):
            p.style = document.styles["Heading 1" if p.style.name == "Heading 2" else "Heading 2"]
            in_front_matter = False
            in_references = p.text == "References"
        elif in_front_matter and p.style.name != "Title":
            fmt.line_spacing, fmt.space_after = 1.05, Pt(2)
        elif in_references:
            fmt.keep_together = True
        elif p.text.startswith(("Table ", "Figure ")):
            fmt.keep_together, fmt.keep_with_next = True, True
            fmt.space_before, fmt.space_after = Pt(6), Pt(4)
            if p.text.startswith("Figure "):
                fmt.line_spacing = 1.5
        elif p._p.xpath(".//w:drawing"):
            fmt.line_spacing = 1
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        elif p._p.getprevious() is not None and p._p.getprevious().tag == qn("w:tbl"):
            fmt.line_spacing, fmt.space_before = 1.5, Pt(5)
            fmt.keep_together = True
    widths = ((2.0, 1.0, 1.25, 1.0, 1.25), (2.15, 0.7, 0.9, 1.15, 1.6),
              (1.05, 2.0, 0.8, 1.15, 1.5))
    for table, row_widths in zip(document.tables, widths):
        format_table(table, row_widths)
    document.save(path)


def build_cover_letter():
    document = Document()
    section = page_setup(document)
    section.top_margin = Inches(1.03)
    normal = document.styles["Normal"]
    normal.font.name, normal.font.size = "Times New Roman", Pt(11.5)
    fmt = normal.paragraph_format
    fmt.line_spacing, fmt.space_after = 1.18, Pt(9)
    blocks = (MANUSCRIPT / "Cover_letter.txt").read_text(encoding="utf-8").strip().split("\n\n")
    for block in blocks:
        document.add_paragraph(block)
    document.core_properties.title = "Cover letter to BMC Genomic Data"
    document.save(MANUSCRIPT / "Cover_letter.docx")


if __name__ == "__main__":
    format_manuscript()
    build_cover_letter()
