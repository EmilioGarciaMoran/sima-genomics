#!/usr/bin/env python3
"""Genera documents/CORRECTION_NOTICE_v2.pdf, el fichero de la versión 2 del preprint en Zenodo.

Los números de las tablas se leen de mtDNA_repro/results/*.tsv (no se teclean) y la figura es
mtDNA_repro/results/figures/figure1_mtDNA_rooted_en.png. El texto resume CORRECTION.md.

uso (desde la raíz del repositorio): python3 tools/build_correction_pdf.py
requiere: reportlab y la fuente DejaVu Sans (/usr/share/fonts/truetype/dejavu)
"""
import csv
import subprocess
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (Image, KeepTogether, Paragraph, SimpleDocTemplate, Spacer, Table,
                                TableStyle)

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "mtDNA_repro" / "results"
OUT = ROOT / "documents" / "CORRECTION_NOTICE_v2.pdf"
REPO = "https://github.com/EmilioGarciaMoran/sima-genomics"
ZENODO_V1 = "10.5281/zenodo.19389351"

FONT_DIR = Path("/usr/share/fonts/truetype/dejavu")
pdfmetrics.registerFont(TTFont("DV", FONT_DIR / "DejaVuSans.ttf"))
pdfmetrics.registerFont(TTFont("DV-B", FONT_DIR / "DejaVuSans-Bold.ttf"))
pdfmetrics.registerFont(TTFont("DV-I", FONT_DIR / "DejaVuSans-Oblique.ttf"))
pdfmetrics.registerFont(TTFont("DV-BI", FONT_DIR / "DejaVuSans-BoldOblique.ttf"))
pdfmetrics.registerFontFamily("DV", normal="DV", bold="DV-B", italic="DV-I", boldItalic="DV-BI")

ss = getSampleStyleSheet()
BODY = ParagraphStyle("body", parent=ss["Normal"], fontName="DV", fontSize=9.2, leading=12.6,
                      spaceAfter=5)
SMALL = ParagraphStyle("small", parent=BODY, fontSize=7.8, leading=10)
H1 = ParagraphStyle("h1", parent=BODY, fontName="DV-B", fontSize=15, leading=19, spaceAfter=4)
H2 = ParagraphStyle("h2", parent=BODY, fontName="DV-B", fontSize=11, leading=14, spaceBefore=8,
                    spaceAfter=4, textColor=colors.HexColor("#1F3A5F"), keepWithNext=1)
CELL = ParagraphStyle("cell", parent=BODY, fontSize=7.8, leading=9.8, spaceAfter=0)
CAP = ParagraphStyle("cap", parent=SMALL, alignment=TA_CENTER, textColor=colors.HexColor("#444444"))


def commit():
    try:
        return subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"],
                                       text=True).strip()
    except Exception:
        return "unknown"


def table(rows, widths, header_bg="#E8EEF5"):
    data = [[Paragraph(str(c), CELL) for c in r] for r in rows]
    t = Table(data, colWidths=widths, repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(header_bg)),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#B8C4D0")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 2.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
    ]))
    return t


def support_rows():
    rows = [["Alignment", "Analysis", "Sima + Denisova", "Sima + Neandertals",
             "Denisova + Neandertals"]]
    for r in csv.DictReader(open(RES / "CLADE_SUPPORT.tsv"), delimiter="\t"):
        aln, method = r["analisis"].split("_", 1)
        method = method.replace("_", " ").replace("stdboot", "standard bootstrap")
        rows.append(["full" if aln == "full" else "no D-loop", method,
                     f"<b>{r['%Sima+Denisova']}</b>", r["%Sima+Neandertales"],
                     r["%Denisova+Neandertales"]])
    return rows


def topology_rows():
    hyp = {"H0_Sima_basal": "H0 Sima basal (withdrawn claim)", "H1_Sima+Denisova": "H1 Sima + Denisova",
           "H2_Sima+Neandertales": "H2 Sima + Neandertals"}
    rows = [["Alignment", "Model", "Hypothesis", "ΔlogL", "p-AU"]]
    for r in csv.DictReader(open(RES / "TOPOLOGY_TESTS.tsv"), delimiter="\t"):
        p = float(r["p_AU"])
        rows.append(["full" if r["alineamiento"] == "full" else "no D-loop", r["modelo"],
                     hyp[r["hipotesis"]], f"{float(r['deltaL']):.2f}",
                     f"<b>{p:.3f}</b>" if p < 0.05 else f"{p:.3f}"])
    return rows


def main():
    W = A4[0] - 3.6 * cm
    s = []
    s.append(Paragraph("Correction notice", H1))
    s.append(Paragraph(
        "<b>Genomic Identity of the Sima de los Huesos Hominins: Convergent Evidence from "
        "Mitochondrial and Nuclear DNA</b> — version 2", BODY))
    s.append(Paragraph(
        f"Emilio García-Morán · 21 October 2026 · Supersedes version 1 (doi:{ZENODO_V1}, April 2026) "
        "and the October 2026 manuscript <i>Phylogenetic placement of Sima de los Huesos hominins "
        "revisited</i>.", SMALL))
    s.append(Spacer(1, 6))

    s.append(Paragraph("Summary", H2))
    s.append(table([
        ["Part", "Status"],
        ["Mitochondrial DNA", "<b>Corrected.</b> With explicit <i>Pan troglodytes</i> rooting, Sima is "
         "sister to Denisova, as reported by Meyer et al. (2014). The “basal Sima” placement is withdrawn."],
        ["Description of Meyer et al. (2014)", "<b>Corrected.</b> They used a Bayesian analysis "
         "(MrBayes, GTR+I+Γ; ~120 mtDNAs including 46 chimpanzee and bonobo genomes; D-loop removed), "
         "not an unrooted Neighbor-Joining tree."],
        ["Nuclear DNA (D-statistics, distances)", "<b>Under verification.</b> Withdrawn as results until "
         "the original ANGSD outputs are recovered and checked."],
        ["Phase 2 mtDNA analyses", "<b>Withdrawn.</b> Input data included non-hominin sequences."],
    ], [4.2 * cm, W - 4.2 * cm]))

    s.append(Paragraph("1. Mitochondrial DNA", H2))
    s.append(Paragraph(
        "<b>What was claimed.</b> Version 1 described Sima de los Huesos as the basal sister lineage to "
        "(Denisova + Neandertal), with 93% IQ-TREE bootstrap support, and interpreted the Sima–Denisova "
        "affinity reported by Meyer et al. (2014) as an artefact of an unrooted analysis.", BODY))
    s.append(Paragraph(
        "<b>What went wrong.</b> The IQ-TREE run was performed without specifying the outgroup; its own "
        "output states that the tree is unrooted and that Sima is only <i>drawn</i> at the root. The 93% "
        "value supports the split {Sima, Denisova} | {Neandertals, <i>Pan</i>}, which, rooted on "
        "<i>Pan</i>, is the Sima + Denisova clade. Reading Sima as basal was a misreading of an unrooted "
        "tree, and the published tree figures were typed by hand rather than drawn from the output. The "
        "original computations were correct: the March 2026 RAxML-NG bootstrap replicates already "
        "supported Sima + Denisova in 98.4%.", BODY))
    s.append(Paragraph(
        "<b>Corrected result.</b> Six complete mitochondrial genomes (Sima Femur XIII KF683087, "
        "Denisova 3, Vindija 33.16, Feldhofer 1, Mezmaiskaya 1, <i>Pan troglodytes</i>), re-analysed "
        "from GenBank with automatic validation of every record, MAFFT 7.526, IQ-TREE 3.0.1 and "
        "RAxML-NG 1.2.2, all rooted on <i>Pan</i>. Bootstrap support (% of 1,000 replicates):", BODY))
    s.append(table(support_rows(), [1.9 * cm, 5.6 * cm, 2.6 * cm, 2.9 * cm, W - 13.0 * cm]))
    s.append(Spacer(1, 6))
    s.append(KeepTogether([Paragraph(
        "Joint topology tests (10,000 RELL replicates) reject “Sima basal” in all four model/alignment "
        "combinations; Sima + Denisova has the highest likelihood in all four:", BODY),
        table(topology_rows(), [1.9 * cm, 2.6 * cm, 6.2 * cm, 2.2 * cm, W - 12.9 * cm])]))
    s.append(Spacer(1, 6))
    img = Image(str(RES / "figures" / "figure1_mtDNA_rooted_en.png"))
    ratio = img.imageHeight / img.imageWidth
    img.drawWidth, img.drawHeight = W, W * ratio
    s.append(KeepTogether([img, Paragraph(
        "Figure 1. Rooted mitochondrial phylogeny drawn directly from the IQ-TREE output "
        "(UFBoot values at nodes; <i>Pan</i> branch shortened for display, real length indicated).",
        CAP)]))
    s.append(Paragraph(
        "This agrees with Meyer et al. (2014) and with Posth et al. (2017), who also rooted on "
        "chimpanzee. The mito-nuclear discordance is explained by the replacement model: Sima retains "
        "the ancestral, Denisovan-like mtDNA of the early Neandertal lineage, later replaced in "
        "Neandertals by an introgressed African mtDNA between ~413 and ~268 ka (Posth et al. 2017); a "
        "parallel replacement is documented for the Y chromosome (Petr et al. 2020).", BODY))

    s.append(Paragraph("2. Phase 2 analyses (withdrawn)", H2))
    s.append(Paragraph(
        "Accession numbers in the phase 2 scripts were typed by hand and never checked against the "
        "records GenBank returned. Ten are not hominin mitochondrial DNA:", BODY))
    s.append(table([
        ["Used as", "Accession", "Actual GenBank record"],
        ["Oase 1, Tianyuan, Kostenki 14, Vindija 33, Denisova 4/8", "KX198073/74/75/78/79",
         "<i>Agaricus</i> (fungus), EF1α gene"],
        ["“Sima SH5”", "KY751401", "<i>Henneguya tucunarei</i> (fish parasite), 18S rRNA"],
        ["Ust’-Ishim", "KP295291", "<i>Phytophthora castaneae</i> (oomycete), HSP90"],
        ["Ust’-Ishim", "KC758896", "<i>Ostrinia scapulalis</i> (moth), ABC transporter"],
        ["“Denisova 2”", "NC_024829", "<i>Nephroselmis astigmatica</i> (green alga), chloroplast genome"],
        ["“Neandertal”", "EU260658", "Simian immunodeficiency virus, partial <i>tat</i> gene"],
    ], [5.6 * cm, 3.4 * cm, W - 9.0 * cm]))
    s.append(Spacer(1, 4))
    s.append(Paragraph(
        "Phase 2 results, including the “6–0” Sima–Denisova shared-site count, are withdrawn. The "
        "six genomes of the corrected analysis are unaffected; their GenBank titles and sequence hashes "
        "are now validated before alignment.", BODY))

    s.append(Paragraph("3. Nuclear DNA (under verification)", H2))
    s.append(Paragraph(
        "The nuclear section of version 1 (chromosome 1 only; D = 0.7425, Z = 30.6; pairwise "
        "distances) cannot yet be reproduced: the ANGSD commands and outputs were not archived with the "
        "code, and the figure scripts contain hand-typed values. These numbers are withdrawn as results "
        "while the original outputs are recovered. A genome-wide, fully reproducible nuclear analysis "
        "is planned as a separate project. The taxonomic proposal of version 1 is withdrawn with it.",
        BODY))

    s.append(Paragraph("4. Data, code and reproducibility", H2))
    s.append(Paragraph(
        f"Repository: <font color='#1F3A5F'>{REPO}</font> (CORRECTION.md; analysis in mtDNA_repro/). "
        "One command, <font name='DV'>bash mtDNA_repro/run_all.sh genbank</font>, reproduces the "
        "analysis from GenBank to the figure (~16 min on 2 cores). An independent clean-clone run "
        "reproduced sequences, alignments and all bootstrap frequencies bit for bit, and the author’s "
        "local run downloading from GenBank returned identical sequence hashes. Superseded scripts and "
        f"manuscripts remain in the repository with notices. This PDF was generated from commit "
        f"{commit()} by tools/build_correction_pdf.py.", BODY))

    s.append(Paragraph("References", H2))
    for ref in [
        "Meyer M. <i>et al.</i> (2014) A mitochondrial genome sequence of a hominin from Sima de los "
        "Huesos. <i>Nature</i> 505:403–406. doi:10.1038/nature12788",
        "Meyer M. <i>et al.</i> (2016) Nuclear DNA sequences from the Middle Pleistocene Sima de los "
        "Huesos hominins. <i>Nature</i> 531:504–507. doi:10.1038/nature17405",
        "Posth C. <i>et al.</i> (2017) Deeply divergent archaic mitochondrial genome provides lower "
        "time boundary for African gene flow into Neanderthals. <i>Nat. Commun.</i> 8:16046. "
        "doi:10.1038/ncomms16046",
        "Petr M. <i>et al.</i> (2020) The evolutionary history of Neanderthal and Denisovan Y "
        "chromosomes. <i>Science</i> 369:1653–1656. doi:10.1126/science.abb6460",
    ]:
        s.append(Paragraph(ref, SMALL))

    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("DV", 7)
        canvas.setFillColor(colors.HexColor("#777777"))
        canvas.drawString(1.8 * cm, 1.1 * cm, "Correction notice · Sima de los Huesos mtDNA · v2 · "
                          "21 Oct 2026")
        canvas.drawRightString(A4[0] - 1.8 * cm, 1.1 * cm, f"page {doc.page}")
        canvas.restoreState()

    doc = SimpleDocTemplate(str(OUT), pagesize=A4, leftMargin=1.8 * cm, rightMargin=1.8 * cm,
                            topMargin=1.6 * cm, bottomMargin=1.7 * cm,
                            title="Correction notice: Sima de los Huesos mtDNA (v2)",
                            author="Emilio García-Morán")
    doc.build(s, onFirstPage=footer, onLaterPages=footer)
    print("Escrito:", OUT)


if __name__ == "__main__":
    main()
