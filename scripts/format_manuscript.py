from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.style import WD_STYLE_TYPE
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from pathlib import Path
p=Path('manuscript/Manuscript.docx');d=Document(p)
for style in d.styles:
 if style.type==WD_STYLE_TYPE.PARAGRAPH:
  style.font.name='Times New Roman';style.font.color.rgb=RGBColor.from_string('111111');style.paragraph_format.line_spacing=2
  if style.name.startswith('Heading'):
   style.font.size=Pt(14 if style.name=='Heading 1' else 12);style.font.bold=True
d.styles['Normal'].font.size=Pt(11);d.styles['Title'].font.size=Pt(19)
for s in d.sections:
 ln=OxmlElement('w:lnNumType');ln.set(qn('w:countBy'),'1');ln.set(qn('w:restart'),'continuous');s._sectPr.insert_element_before(ln,'w:cols')
 f=s.footer.paragraphs[0];f.alignment=2
 fld=OxmlElement('w:fldSimple');fld.set(qn('w:instr'),'PAGE');f._p.append(fld)
for q in d.paragraphs:
 if q.text.startswith('Figure ') and q.runs and q.runs[0].bold:q.paragraph_format.keep_with_next=True;q.paragraph_format.keep_together=True
d.save(p)
