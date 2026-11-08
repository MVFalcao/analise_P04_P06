#!/usr/bin/python3
# =============================================================================
# Etapa 10 (versão Python) — Perfil das variantes a partir do MAF
#
# Faz as mesmas 7 figuras do analise_maftools.R e mais 3 figuras pensadas
# para leitura clínica (08–10).
#
# USO (no terminal, a partir da pasta do projeto):
#   /usr/bin/python3 scripts/analise_variantes.py  resultados/09_maf/coorte.maf  [dados/dados_clinicos.tsv]
#
# Figuras salvas em resultados/10_python/:
#   01_resumo_todas_variantes.png     painel-resumo de todas as variantes
#   02_oncoplot_todas_variantes.png   oncoplot (gene × paciente) de todas
#   03_resumo_variantes_raras.png     painel-resumo só das raras
#   04_oncoplot_variantes_raras.png   oncoplot das raras
#   05_oncoplot_sem_multihit.png      oncoplot das raras, 1 variante por gene/paciente
#   06_titv.png                       transições/transversões (controle de qualidade)
#   07_lollipop_TP53.png              posição das variantes raras na proteína TP53
#   08_painel_por_paciente.png        ficha por paciente: clínica + variantes raras
#   09_vaf_vs_gnomad.png              frequência no tumor (VAF) × na população (gnomAD)
#   10_seguimento.png                 tempo de seguimento e desfecho por paciente
#
# Use o Python do sistema (/usr/bin/python3, com pandas e matplotlib): o Python
# 3.15 do ambiente aula-variantes ainda não tem pandas.
# Regra usada: gnomAD AF >= 1% = polimorfismo comum; < 1% ou ausente = rara.
# Dados clínicos SIMULADOS (aula). As figuras são didáticas, não são laudo.
# =============================================================================
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle
from matplotlib.ticker import MaxNLocator

CORTE_AF = 0.01
SAIDA = "resultados/10_python"

# ---------------------------------------------------------------------------
# Cores e rótulos (cada classe tem sempre a mesma cor em todas as figuras)
# ---------------------------------------------------------------------------
CLASSES = [  # da mais para a menos deletéria
    "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
    "Translation_Start_Site", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins",
    "Missense_Mutation", "Silent",
]
COR = {
    "Missense_Mutation": "#2a78d6",
    "Nonsense_Mutation": "#eb6834",
    "Frame_Shift_Del": "#1baf7a",
    "Frame_Shift_Ins": "#eda100",
    "Splice_Site": "#e87ba4",
    "In_Frame_Del": "#008300",
    "In_Frame_Ins": "#008300",
    "Translation_Start_Site": "#e34948",
    "Nonstop_Mutation": "#e34948",
    "Multi_Hit": "#6250d6",
    "Silent": "#a3a29b",
}
NOME = {
    "Missense_Mutation": "Missense (troca de aminoácido)",
    "Nonsense_Mutation": "Nonsense (códon de parada)",
    "Frame_Shift_Del": "Frameshift (deleção)",
    "Frame_Shift_Ins": "Frameshift (inserção)",
    "Splice_Site": "Sítio de splicing",
    "In_Frame_Del": "Deleção in-frame",
    "In_Frame_Ins": "Inserção in-frame",
    "Translation_Start_Site": "Perda do códon inicial",
    "Nonstop_Mutation": "Perda do códon de parada",
    "Multi_Hit": "Multi_Hit (2+ variantes no gene)",
    "Silent": "Sinônima (silenciosa)",
}
TRUNCANTES = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"}
COR_PACIENTE = ["#2a78d6", "#eb6834", "#1baf7a"]  # até 3 pacientes por cor; acima disso, cinza
TINTA, TINTA2, TINTA3 = "#0b0b0b", "#52514e", "#8a8984"
GRADE, VAZIO = "#e6e5e0", "#f1f0ec"

# Papel dos genes (conhecimento consolidado, para orientar a leitura clínica)
GENES = {
    "MLH1": ("Reparo de DNA (MMR)",
             "Variante germinativa truncante causa síndrome de Lynch."),
    "MSH2": ("Reparo de DNA (MMR)", "Variante germinativa causa síndrome de Lynch."),
    "MSH6": ("Reparo de DNA (MMR)", "Variante germinativa causa síndrome de Lynch."),
    "PMS2": ("Reparo de DNA (MMR)", "Variante germinativa causa síndrome de Lynch."),
    "APC": ("Supressor tumoral (via WNT)",
            "Variante germinativa truncante causa polipose adenomatosa familiar (PAF)."),
    "KRAS": ("Oncogene (via RAS/MAPK)",
             "G12/G13: driver somático frequente no colorretal; prediz resistência a anti-EGFR."),
    "TP53": ("Supressor tumoral",
             "R175H é hotspot somático; se germinativa, considerar síndrome de Li-Fraumeni."),
    "BRCA1": ("Reparo por recombinação homóloga", "Germinativa patogênica: câncer de mama/ovário hereditário."),
    "BRCA2": ("Reparo por recombinação homóloga", "Germinativa patogênica: câncer de mama/ovário hereditário."),
    "EGFR": ("Oncogene (receptor tirosina-quinase)", "Alvo terapêutico em alguns tumores."),
}
SINDROME = {"MLH1": "síndrome de Lynch", "MSH2": "síndrome de Lynch", "MSH6": "síndrome de Lynch",
            "PMS2": "síndrome de Lynch", "APC": "polipose adenomatosa familiar (PAF)"}

# Proteínas para o lollipop (tamanho em aminoácidos e domínios Pfam aproximados)
PROTEINAS = {
    "TP53": (393, [("Transativação", 6, 30), ("Ligação ao DNA", 95, 288), ("Tetramerização", 319, 357)]),
}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": TINTA3,
    "axes.labelcolor": TINTA2, "xtick.color": TINTA2, "ytick.color": TINTA2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.titleweight": "bold",
    "axes.titlesize": 11, "axes.titlecolor": TINTA, "figure.facecolor": "white",
})


# ---------------------------------------------------------------------------
# Leitura e preparo
# ---------------------------------------------------------------------------
def ler_maf(caminho):
    df = pd.read_csv(caminho, sep="\t", dtype=str, na_values=[".", "", "NA"], keep_default_na=False)
    df["gnomAD_AF"] = pd.to_numeric(df["gnomAD_AF"].str.split("&").str[0], errors="coerce")
    for c in ["VAF", "t_depth", "t_ref_count", "t_alt_count", "Start_Position"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["rara"] = df["gnomAD_AF"].isna() | (df["gnomAD_AF"] < CORTE_AF)
    df["classe_snv"] = df.apply(classe_snv, axis=1)
    return df


def classe_snv(l):
    if l["Variant_Type"] != "SNP":
        return None
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    r, a = l["Reference_Allele"], l["Tumor_Seq_Allele2"]
    if r in ("A", "G"):
        r, a = comp[r], comp.get(a, a)
    return f"{r}>{a}"


def ler_clinico(caminho, amostras):
    if caminho and os.path.exists(caminho):
        clin = pd.read_csv(caminho, sep="\t", dtype=str)
    else:
        clin = pd.DataFrame({"Tumor_Sample_Barcode": amostras})
    return clin.set_index("Tumor_Sample_Barcode")


def ordem_classes(df):
    presentes = set(df["Variant_Classification"])
    return [c for c in CLASSES if c in presentes]


def salvar(fig, nome):
    caminho = os.path.join(SAIDA, nome)
    fig.savefig(caminho, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  figura:", caminho)


def rodape(fig, texto="Dados clínicos simulados (aula IEC – PUC Minas). Figura didática, não é laudo."):
    fig.text(0.01, -0.01, texto, fontsize=8, color=TINTA3, ha="left", va="top")


def barras_h(ax, rotulos, valores, cores, titulo, fmt="{:d}"):
    y = np.arange(len(rotulos))[::-1]
    ax.barh(y, valores, color=cores, height=0.6)
    ax.set_yticks(y, rotulos)
    ax.set_title(titulo, loc="left")
    m = max(valores) if len(valores) else 1
    for yi, v in zip(y, valores):
        ax.text(v + m * 0.02, yi, fmt.format(v), va="center", fontsize=9, color=TINTA2)
    ax.set_xlim(0, m * 1.18 if m else 1)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(axis="x", color=GRADE, lw=0.8)
    ax.set_axisbelow(True)


# ---------------------------------------------------------------------------
# 01 / 03 — Painel-resumo (equivalente ao plotmafSummary)
# ---------------------------------------------------------------------------
def fig_resumo(df, titulo, nome):
    classes = ordem_classes(df)
    amostras = sorted(df["Tumor_Sample_Barcode"].unique())
    fig, axs = plt.subplots(2, 3, figsize=(15, 8.5))
    fig.suptitle(titulo, x=0.01, ha="left", fontsize=14, fontweight="bold")

    cont = df["Variant_Classification"].value_counts()
    barras_h(axs[0, 0], [NOME[c] for c in classes], [int(cont[c]) for c in classes],
             [COR[c] for c in classes], "Classificação das variantes")

    tipos = df["Variant_Type"].value_counts()
    barras_h(axs[0, 1], list(tipos.index), [int(v) for v in tipos.values],
             [TINTA3] * len(tipos), "Tipo (SNP, DEL, INS)")

    snv = df["classe_snv"].dropna().value_counts().reindex(
        ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"], fill_value=0)
    barras_h(axs[0, 2], list(snv.index), [int(v) for v in snv.values],
             ["#2a78d6" if c in ("C>T", "T>C") else "#a3a29b" for c in snv.index],
             "Classe das SNVs (azul = transição)")

    # variantes por paciente, empilhadas por classe
    ax = axs[1, 0]
    tab = pd.crosstab(df["Tumor_Sample_Barcode"], df["Variant_Classification"]).reindex(
        index=amostras, columns=classes, fill_value=0)
    base = np.zeros(len(amostras))
    for c in classes:
        ax.bar(amostras, tab[c], bottom=base, color=COR[c], width=0.55, edgecolor="white", lw=1.5)
        base += tab[c].values
    med = np.median(base)
    ax.axhline(med, ls="--", color=TINTA2, lw=1)
    for i, t in enumerate(base):
        ax.text(i, t, f"{int(t)}", ha="center", va="bottom", fontsize=9, color=TINTA2)
    ax.set_title(f"Variantes por paciente (mediana {med:g}, linha tracejada)", loc="left")
    ax.set_ylabel("nº de variantes")
    ax.set_ylim(0, base.max() * 1.15)
    ax.grid(axis="y", color=GRADE, lw=0.8)
    ax.set_axisbelow(True)

    # resumo por classe: cada ponto = um paciente, traço = mediana
    ax = axs[1, 1]
    for i, c in enumerate(classes):
        v = tab[c].values
        jit = np.linspace(-0.12, 0.12, len(v)) if len(v) > 1 else [0]
        ax.scatter(np.full(len(v), i) + jit, v, color=COR[c], s=50, zorder=3, edgecolor="white", lw=1.5)
        ax.plot([i - 0.25, i + 0.25], [np.median(v)] * 2, color=TINTA, lw=2)
    ax.set_xticks(range(len(classes)), [NOME[c].split(" (")[0] for c in classes], rotation=30, ha="right")
    ax.set_title("Por classe: cada ponto = 1 paciente", loc="left")
    ax.set_ylabel("nº de variantes")
    ax.grid(axis="y", color=GRADE, lw=0.8)
    ax.set_axisbelow(True)

    # genes mais afetados
    ax = axs[1, 2]
    genes = (df.groupby("Hugo_Symbol")
               .agg(n=("Hugo_Symbol", "size"), amostras=("Tumor_Sample_Barcode", "nunique"))
               .sort_values(["amostras", "n"], ascending=False).head(10))
    y = np.arange(len(genes))[::-1]
    sub = pd.crosstab(df["Hugo_Symbol"], df["Variant_Classification"]).reindex(
        index=genes.index, columns=classes, fill_value=0)
    left = np.zeros(len(genes))
    for c in classes:
        ax.barh(y, sub[c], left=left, color=COR[c], height=0.6, edgecolor="white", lw=1.5)
        left += sub[c].values
    for yi, (g, r) in zip(y, genes.iterrows()):
        ax.text(left[list(genes.index).index(g)] + 0.15, yi,
                f"{100 * r['amostras'] / len(amostras):.0f}% dos pacientes", va="center", fontsize=9, color=TINTA2)
    ax.set_yticks(y, genes.index, fontstyle="italic")
    ax.set_xlim(0, left.max() * 1.6)
    ax.set_title("Genes mais afetados", loc="left")
    ax.set_xlabel("nº de variantes")
    ax.grid(axis="x", color=GRADE, lw=0.8)
    ax.set_axisbelow(True)

    fig.legend(handles=[Patch(color=COR[c], label=NOME[c]) for c in classes],
               loc="lower center", ncol=min(len(classes), 5), frameon=False, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    rodape(fig)
    fig.texts[-1].set_position((0.01, -0.06))
    salvar(fig, nome)


# ---------------------------------------------------------------------------
# 02 / 04 / 05 — Oncoplot (gene × paciente)
# ---------------------------------------------------------------------------
def fig_oncoplot(df, clin, titulo, nome, top=50):
    amostras_todas = sorted(clin.index) if len(clin) else sorted(df["Tumor_Sample_Barcode"].unique())
    resumo = (df.groupby("Hugo_Symbol")
                .agg(n=("Hugo_Symbol", "size"), am=("Tumor_Sample_Barcode", "nunique"))
                .sort_values(["am", "n"], ascending=False).head(top))
    genes = list(resumo.index)

    # célula = classe, ou Multi_Hit se houver 2+ variantes do gene no mesmo paciente
    celula, texto = {}, {}
    for (a, g), sub in df.groupby(["Tumor_Sample_Barcode", "Hugo_Symbol"]):
        celula[(a, g)] = "Multi_Hit" if len(sub) > 1 else sub["Variant_Classification"].iloc[0]
        texto[(a, g)] = (", ".join(sub["HGVSp_Short"].fillna("?")) if len(sub) <= 2
                         else f"{len(sub)} variantes")

    # ordena pacientes: tipo tumoral, depois padrão de mutação (como o maftools)
    def chave(a):
        tipo = clin.loc[a, "Tipo_tumoral"] if "Tipo_tumoral" in clin.columns and a in clin.index else ""
        return (tipo, a)  # mesma ordem de pacientes em todas as figuras
    amostras = sorted(amostras_todas, key=chave)

    poucos = len(amostras) <= 8  # com poucos pacientes, escreve a alteração dentro da célula
    larg_cel = 2.4 if poucos else 0.5
    fig_w = 3.5 + larg_cel * len(amostras) + 2.2
    fig_h = 2.2 + 0.42 * len(genes) + 1.3
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(3, 2, width_ratios=[larg_cel * len(amostras), 2.2],
                          height_ratios=[1.4, 0.42 * len(genes), 0.9], hspace=0.08, wspace=0.03)
    ax_top, ax = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])
    ax_dir, ax_clin = fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[2, 0])

    classes_usadas = set()
    for j, a in enumerate(amostras):
        for i, g in enumerate(genes):
            c = celula.get((a, g))
            ax.add_patch(Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, color=COR[c] if c else VAZIO, lw=0))
            if c:
                classes_usadas.add(c)
                if poucos:
                    t = texto[(a, g)]
                    t = t if len(t) <= 26 else t[:24] + "…"
                    ax.text(j + 0.5, i + 0.5, t, ha="center", va="center", fontsize=8,
                            color="white" if c != "Silent" else TINTA, fontweight="bold")
    ax.set_xlim(0, len(amostras))
    ax.set_ylim(len(genes), 0)
    ax.set_yticks(np.arange(len(genes)) + 0.5, genes, fontstyle="italic")
    ax.set_xticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)

    # barra superior: nº de variantes por paciente
    tab = pd.crosstab(df["Tumor_Sample_Barcode"], df["Variant_Classification"]).reindex(
        index=amostras, fill_value=0)
    base = np.zeros(len(amostras))
    for c in [c for c in CLASSES if c in tab.columns]:
        ax_top.bar(np.arange(len(amostras)) + 0.5, tab[c], bottom=base, width=0.6 if poucos else 0.85,
                   color=COR[c], edgecolor="white", lw=1.2)
        base += tab[c].values
    for j, t in enumerate(base):
        ax_top.text(j + 0.5, t, f"{int(t)}", ha="center", va="bottom", fontsize=9, color=TINTA2)
    ax_top.set_xlim(0, len(amostras))
    ax_top.set_ylim(0, max(base.max(), 1) * 1.25)
    ax_top.set_xticks(np.arange(len(amostras)) + 0.5, amostras, fontweight="bold")
    ax_top.xaxis.tick_top()
    ax_top.set_ylabel("variantes", fontsize=9)
    ax_top.spines["bottom"].set_visible(False)
    ax_top.set_title(titulo, loc="left", fontsize=13, pad=26)

    # barra à direita: % de pacientes com o gene alterado
    pct = 100 * resumo["am"] / len(amostras)
    ax_dir.barh(np.arange(len(genes)) + 0.5, pct.values, height=0.55, color=TINTA3)
    for i, v in enumerate(pct.values):
        ax_dir.text(v + 2, i + 0.5, f"{v:.0f}%", va="center", fontsize=9, color=TINTA2)
    ax_dir.set_ylim(len(genes), 0)
    ax_dir.set_xlim(0, 135)
    ax_dir.set_yticks([])
    ax_dir.set_xticks([])
    ax_dir.spines["left"].set_visible(False)
    ax_dir.spines["bottom"].set_visible(False)
    ax_dir.set_title("% pacientes", fontsize=9, color=TINTA2, fontweight="normal")

    # trilhas clínicas: tipo tumoral, idade, sexo
    trilhas = [c for c in ["Tipo_tumoral", "Idade", "Sexo"] if c in clin.columns]
    for i, col in enumerate(trilhas):
        for j, a in enumerate(amostras):
            v = clin.loc[a, col] if a in clin.index else "?"
            ax_clin.add_patch(Rectangle((j + 0.04, i + 0.1), 0.92, 0.8, color="#cde2fb", lw=0))
            ax_clin.text(j + 0.5, i + 0.5, f"{v} anos" if col == "Idade" else v,
                         ha="center", va="center", fontsize=8.5, color=TINTA)
    ax_clin.set_xlim(0, len(amostras))
    ax_clin.set_ylim(max(len(trilhas), 1), 0)
    ax_clin.set_yticks(np.arange(len(trilhas)) + 0.5, [t.replace("_", " ") for t in trilhas], fontsize=9)
    ax_clin.set_xticks([])
    for s in ax_clin.spines.values():
        s.set_visible(False)
    ax_clin.tick_params(length=0)

    ordem = [c for c in CLASSES + ["Multi_Hit"] if c in classes_usadas]
    fig.legend(handles=[Patch(color=COR[c], label=NOME[c]) for c in ordem], loc="upper center",
               bbox_to_anchor=(0.45, 0.0), ncol=min(len(ordem), 3), frameon=False, fontsize=9)
    rodape(fig, "Cada célula = gene alterado no paciente; o texto mostra a alteração na proteína. "
                "Dados clínicos simulados (aula).")
    fig.texts[-1].set_position((0.01, -0.09))
    salvar(fig, nome)


# ---------------------------------------------------------------------------
# 06 — Transições / transversões
# ---------------------------------------------------------------------------
def fig_titv(df, nome):
    snv = df.dropna(subset=["classe_snv"])
    classes6 = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    frac = pd.crosstab(snv["Tumor_Sample_Barcode"], snv["classe_snv"], normalize="index").reindex(
        columns=classes6, fill_value=0) * 100
    amostras = list(frac.index)
    cor6 = dict(zip(classes6, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]))
    ti = snv["classe_snv"].isin(["C>T", "T>C"])

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8), gridspec_kw={"width_ratios": [1.4, 0.8, 1.2]})
    fig.suptitle("Transições e transversões (controle de qualidade)", x=0.01, ha="left",
                 fontsize=14, fontweight="bold")

    ax = axs[0]
    for i, c in enumerate(classes6):
        v = frac[c].values
        jit = np.linspace(-0.12, 0.12, len(v)) if len(v) > 1 else [0]
        ax.scatter(np.full(len(v), i) + jit, v, color=cor6[c], s=45, zorder=3, edgecolor="white", lw=1.5)
        ax.plot([i - 0.25, i + 0.25], [np.median(v)] * 2, color=TINTA, lw=2)
    ax.set_xticks(range(6), classes6)
    ax.set_ylabel("% das SNVs do paciente")
    ax.set_title("Por classe (ponto = paciente)", loc="left")
    ax.grid(axis="y", color=GRADE, lw=0.8)
    ax.set_axisbelow(True)

    ax = axs[1]
    vals = [100 * ti.mean(), 100 * (~ti).mean()]
    ax.bar(["Transição", "Transversão"], vals, color=["#2a78d6", "#a3a29b"], width=0.55)
    for i, v in enumerate(vals):
        ax.text(i, v + 1.5, f"{v:.0f}%", ha="center", fontsize=10, color=TINTA2)
    ax.set_ylim(0, 110)
    razao = ti.sum() / max((~ti).sum(), 1)
    ax.set_title(f"Ti/Tv = {razao:.2f}".replace(".", ","), loc="left")
    ax.set_ylabel("% de todas as SNVs")
    ax.grid(axis="y", color=GRADE, lw=0.8)
    ax.set_axisbelow(True)

    ax = axs[2]
    base = np.zeros(len(amostras))
    for c in classes6:
        ax.bar(amostras, frac[c], bottom=base, color=cor6[c], width=0.55, edgecolor="white", lw=1.5, label=c)
        base += frac[c].values
    ax.set_ylabel("% das SNVs")
    ax.set_title("Composição por paciente", loc="left")
    ax.legend(frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9)
    fig.tight_layout(rect=(0, 0.02, 1, 0.93))
    rodape(fig, "Em exomas humanos, Ti/Tv perto de 2–3 indica chamadas de boa qualidade. "
                "Com poucas variantes, o valor é só indicativo.")
    salvar(fig, nome)


# ---------------------------------------------------------------------------
# 07 — Lollipop (posição na proteína)
# ---------------------------------------------------------------------------
def fig_lollipop(df, gene, n_amostras, nome):
    tamanho, dominios = PROTEINAS[gene]
    sub = df[df["Hugo_Symbol"] == gene].copy()
    sub["pos"] = sub["HGVSp_Short"].str.extract(r"p\.[A-Z*]?(\d+)")[0].astype(float)
    sub = sub.dropna(subset=["pos"])
    grupos = sub.groupby(["pos", "HGVSp_Short", "Variant_Classification"]).agg(
        n=("Tumor_Sample_Barcode", "size"), pac=("Tumor_Sample_Barcode", lambda s: ", ".join(sorted(s)))).reset_index()

    fig, ax = plt.subplots(figsize=(12, 4.8))
    cores_dom = ["#cde2fb", "#9ec5f4", "#cde2fb"]
    ax.add_patch(Rectangle((0, -0.25), tamanho, 0.5, color=GRADE, lw=0))
    for (nome_d, ini, fim), c in zip(dominios, cores_dom):
        ax.add_patch(FancyBboxPatch((ini, -0.42), fim - ini, 0.84, boxstyle="round,pad=0,rounding_size=4",
                                    color=c, lw=0, mutation_aspect=0.02))
        ax.text((ini + fim) / 2, 0, nome_d, ha="center", va="center", fontsize=9, color=TINTA)
    alt_max = max(grupos["n"].max(), 1) if len(grupos) else 1
    for _, r in grupos.iterrows():
        h = 1 + r["n"]
        ax.plot([r["pos"]] * 2, [0.42, h], color=TINTA3, lw=1.5, zorder=1)
        ax.scatter(r["pos"], h, s=160, color=COR[r["Variant_Classification"]], zorder=3,
                   edgecolor="white", lw=2)
        ax.text(r["pos"], h + 0.35, f"{r['HGVSp_Short']}\n({r['pac']})", ha="center", va="bottom",
                fontsize=9, color=TINTA)
    ax.set_xlim(-5, tamanho + 5)
    ax.set_ylim(-0.8, alt_max + 2.6)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("posição do aminoácido")
    taxa = 100 * sub["Tumor_Sample_Barcode"].nunique() / n_amostras
    ax.set_title(f"{gene}: variantes raras na proteína ({tamanho} aa) · "
                 f"{taxa:.0f}% dos pacientes", loc="left", fontsize=13)
    usadas = [c for c in CLASSES if c in set(grupos["Variant_Classification"])]
    ax.legend(handles=[Patch(color=COR[c], label=NOME[c]) for c in usadas], frameon=False,
              loc="upper right", fontsize=9)
    rodape(fig, "Domínios Pfam aproximados. Altura da haste = nº de pacientes com a variante.")
    salvar(fig, nome)


# ---------------------------------------------------------------------------
# 08 — Ficha por paciente (clínica + variantes raras + hipótese)
# ---------------------------------------------------------------------------
def tipo_curto(c):
    return {"Nonsense_Mutation": "Nonsense (truncante)", "Frame_Shift_Del": "Frameshift (truncante)",
            "Frame_Shift_Ins": "Frameshift (truncante)", "Splice_Site": "Splicing (truncante)",
            "Missense_Mutation": "Missense"}.get(c, NOME.get(c, c))


def fmt_af(v):
    if pd.isna(v):
        return "ausente"
    t = f"{100 * v:.4f}%" if v < 0.0001 else f"{100 * v:.2f}%"
    return t.replace(".", ",")


def fig_painel_pacientes(df, clin, nome):
    amostras = sorted(df["Tumor_Sample_Barcode"].unique())
    raras = df[df["rara"]]
    alturas = [2.75 + 0.62 * max(len(raras[raras["Tumor_Sample_Barcode"] == a]), 1) for a in amostras]
    fig = plt.figure(figsize=(15, sum(alturas) + 0.3 * len(amostras)))
    gs = fig.add_gridspec(len(amostras), 1, height_ratios=alturas, hspace=0.08,
                          left=0.01, right=0.99, top=0.95, bottom=0.02)
    fig.suptitle("Resumo por paciente: dados clínicos e variantes raras (gnomAD < 1%)",
                 x=0.01, ha="left", fontsize=14, fontweight="bold")
    colunas = [("Gene", 0.01), ("Alteração", 0.08), ("Tipo", 0.20), ("VAF no tumor", 0.355),
               ("gnomAD", 0.45), ("Papel do gene / relevância", 0.53)]

    for k, a in enumerate(amostras):
        ax = fig.add_subplot(gs[k])
        ax.axis("off")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.add_patch(FancyBboxPatch((0, 0), 1, 1, boxstyle="round,pad=0,rounding_size=0.015",
                                    fc="#fcfcfb", ec=GRADE, lw=1.2, transform=ax.transAxes))
        c = clin.loc[a] if a in clin.index else pd.Series(dtype=str)
        g = lambda col: c.get(col, "?") if len(c) else "?"
        da = raras[raras["Tumor_Sample_Barcode"] == a].copy()
        da["ordem"] = da["Variant_Classification"].map({v: i for i, v in enumerate(CLASSES)})
        da = da.sort_values("ordem")
        n_comuns = int((~df["rara"] & (df["Tumor_Sample_Barcode"] == a)).sum())

        h = alturas[k]
        lin = lambda pol: 1 - pol / h  # converte polegadas a partir do topo em coordenada 0–1
        ax.text(0.015, lin(0.38), f"Paciente {a}", fontsize=14, fontweight="bold", color=TINTA, va="center")
        ax.text(0.015, lin(0.75), f"{g('Idade')} anos · sexo {g('Sexo')} · {g('Diagnostico')}",
                fontsize=11, color=TINTA, va="center")
        ax.text(0.015, lin(1.05), f"História familiar: {g('Historia_familiar')}", fontsize=10, color=TINTA2, va="center")
        status = {"0": "vivo / sem evento", "1": "evento (óbito)"}.get(str(g("status")), str(g("status")))
        ax.text(0.015, lin(1.33), f"Seguimento: {g('tempo_meses')} meses · {status}   |   "
                f"{len(da)} variantes raras · {n_comuns} polimorfismos comuns filtrados",
                fontsize=10, color=TINTA2, va="center")

        y0 = 1.78
        for titulo_c, x in colunas:
            ax.text(x + 0.005, lin(y0), titulo_c, fontsize=9, color=TINTA3, fontweight="bold", va="center")
        ax.plot([0.012, 0.988], [lin(y0 + 0.18)] * 2, color=GRADE, lw=1)
        for i, (_, r) in enumerate(da.iterrows()):
            y = lin(y0 + 0.45 + 0.62 * i)
            cor = COR[r["Variant_Classification"]]
            ax.add_patch(Rectangle((0.012, y - 0.16 / h * 1.4), 0.004, 0.32 / h * 1.4, color=cor, lw=0))
            papel, nota = GENES.get(r["Hugo_Symbol"], ("", ""))
            vals = [r["Hugo_Symbol"], r["HGVSp_Short"], tipo_curto(r["Variant_Classification"]),
                    f"{100 * r['VAF']:.0f}%" if pd.notna(r["VAF"]) else "?", fmt_af(r["gnomAD_AF"])]
            for (_, x), v in zip(colunas, vals):
                ax.text(x + 0.005, y, v, fontsize=10.5, va="center", color=TINTA,
                        fontstyle="italic" if x == 0.01 else "normal",
                        fontweight="bold" if x in (0.01, 0.08) else "normal")
            ax.text(0.535, y + 0.09 / h * 1.4, papel, fontsize=9.5, va="center", color=TINTA, fontweight="bold")
            ax.text(0.535, y - 0.11 / h * 1.4, nota, fontsize=9, va="center", color=TINTA2)
        if da.empty:
            ax.text(0.015, lin(y0 + 0.45), "Nenhuma variante rara.", fontsize=10, color=TINTA2, va="center")

        # hipótese: variante truncante em gene de síndrome hereditária
        #    (lista todas; a mais próxima de VAF 50% primeiro, padrão de variante germinativa heterozigota)
        cand = da[da["Hugo_Symbol"].isin(SINDROME) & da["Variant_Classification"].isin(TRUNCANTES)].copy()
        if not cand.empty:
            cand = cand.assign(_d=(cand["VAF"] - 0.5).abs()).sort_values("_d")
            hip = " · ".join(f"{SINDROME[r['Hugo_Symbol']]} ({r['Hugo_Symbol']}, VAF {100 * r['VAF']:.0f}%)"
                             for _, r in cand.iterrows())
            ax.text(0.985, lin(0.38), f"Investigar: {hip}", ha="right", va="center",
                    fontsize=11, fontweight="bold", color="#b34514",
                    bbox=dict(boxstyle="round,pad=0.35", fc="#fdeee6", ec="#f3c3ab"))
        if not cand.empty:
            ax.text(0.015, lin(h - 0.32), "Variante truncante em gene de síndrome hereditária: confirmar origem "
                    "germinativa (sangue) e encaminhar ao aconselhamento genético.",
                    va="center", fontsize=9.5, color="#b34514")
    rodape(fig, "VAF = fração de leituras com a variante no tumor. Sequenciamento só do tumor: "
                "a VAF sozinha não define se a variante é germinativa ou somática. "
                "Dados clínicos simulados (aula); não é laudo.")
    salvar(fig, nome)


# ---------------------------------------------------------------------------
# 09 — VAF no tumor × frequência na população
# ---------------------------------------------------------------------------
def fig_vaf_gnomad(df, nome):
    amostras = sorted(df["Tumor_Sample_Barcode"].unique())
    cores = {a: (COR_PACIENTE[i] if i < len(COR_PACIENTE) else TINTA3) for i, a in enumerate(amostras)}
    marcas = ["o", "^", "s", "D", "v", "P"]
    x_aus = 1e-7  # posição usada para "ausente no gnomAD"
    fig, ax = plt.subplots(figsize=(13, 6.5))
    ax.axvspan(x_aus / 3, CORTE_AF, color="#fdeee6", lw=0, zorder=0)
    ax.axvline(CORTE_AF, color="#b34514", ls="--", lw=1.2)
    ax.text(CORTE_AF * 1.15, 1.05, "corte 1%", color="#b34514", fontsize=9, va="center")
    ax.text(x_aus / 2.5, 1.05, "RARAS (gnomAD < 1%): podem explicar o caso", color="#b34514",
            fontsize=10, fontweight="bold", va="center")
    ax.text(0.05, 1.05, "COMUNS (≥ 1%): polimorfismos, filtrados", color=TINTA2,
            fontsize=10, fontweight="bold", va="center")
    ax.axhspan(0.40, 0.60, color="#eef4fc", lw=0, zorder=0)
    ax.text(1.5, 0.5, "~50%\nheterozigota", fontsize=8.5, color=TINTA2, va="center")
    ax.text(1.5, 0.97, "~100%\nhomozigota", fontsize=8.5, color=TINTA2, va="center")

    rotulos = []
    for k, a in enumerate(amostras):
        d = df[df["Tumor_Sample_Barcode"] == a]
        x = d["gnomAD_AF"].fillna(x_aus).values * (1.08 if k % 2 else 0.93)  # leve deslocamento entre pacientes
        ax.scatter(x, d["VAF"], s=np.where(d["rara"], 130, 55), marker=marcas[k % len(marcas)],
                   color=cores[a], edgecolor="white", lw=1.5, zorder=3, alpha=np.where(d["rara"], 1, 0.55),
                   label=f"Paciente {a}")
        for (_, r), xi in zip(d[d["rara"]].iterrows(), x[d["rara"].values]):
            rotulos.append([xi, r["VAF"], r["VAF"], f"{a} · {r['Hugo_Symbol']} {r['HGVSp_Short']}"])
    # afasta rótulos verticalmente quando ficam próximos (mínimo 4,5 pontos percentuais)
    rotulos.sort(key=lambda t: t[1])
    for i in range(1, len(rotulos)):
        if rotulos[i][2] - rotulos[i - 1][2] < 0.045:
            rotulos[i][2] = rotulos[i - 1][2] + 0.045
    for xi, yi, yl, txt in rotulos:
        ax.annotate(txt, (xi, yi), xytext=(xi * 2.2, yl), textcoords="data", fontsize=9, va="center",
                    color=TINTA, arrowprops=dict(arrowstyle="-", color=TINTA3, lw=0.8, shrinkA=0, shrinkB=6))
    ax.set_xscale("log")
    ax.set_xlim(x_aus / 3, 3)
    ax.set_ylim(0, 1.12)
    ticks = [x_aus, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]
    ax.set_xticks(ticks, ["ausente", "0,0001%", "0,001%", "0,01%", "0,1%", "1%", "10%", "100%"])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1], ["0%", "25%", "50%", "75%", "100%"])
    ax.set_xlabel("Frequência na população (gnomAD exomas, escala log)")
    ax.set_ylabel("VAF no tumor (% das leituras com a variante)")
    ax.set_title("Quais variantes importam? Frequência no tumor × na população", loc="left", fontsize=13, pad=22)
    ax.grid(color=GRADE, lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower right")
    rodape(fig, "Pontos grandes e rotulados = variantes raras. Pontos claros = polimorfismos comuns "
                "(aparecem nos dois pacientes). Dados simulados (aula).")
    salvar(fig, nome)


# ---------------------------------------------------------------------------
# 10 — Seguimento (swimmer plot)
# ---------------------------------------------------------------------------
def fig_seguimento(df, clin, nome):
    if not {"tempo_meses", "status"}.issubset(clin.columns):
        print("  (sem tempo_meses/status nos dados clínicos — figura 10 não gerada)")
        return
    raras = df[df["rara"]]
    c = clin.copy()
    c["tempo"] = pd.to_numeric(c["tempo_meses"], errors="coerce")
    c = c.sort_values("tempo")
    fig, ax = plt.subplots(figsize=(12, 1.2 + 1.1 * len(c)))
    for i, (a, r) in enumerate(c.iterrows()):
        cor = "#2a78d6" if str(r["status"]) == "0" else "#eb6834"
        ax.barh(i, r["tempo"], height=0.38, color=cor)
        if str(r["status"]) == "1":
            ax.scatter(r["tempo"], i, marker="X", s=170, color=TINTA, zorder=3)
        else:
            ax.annotate("", xy=(r["tempo"] + 3, i), xytext=(r["tempo"], i),
                        arrowprops=dict(arrowstyle="-|>", color=TINTA, lw=1.5))
        genes = ", ".join(f"{g} {p}" for g, p in raras[raras["Tumor_Sample_Barcode"] == a]
                          [["Hugo_Symbol", "HGVSp_Short"]].itertuples(index=False))
        ax.text(0.5, i + 0.36, genes, fontsize=9, color=TINTA2, va="center")
        ax.text(r["tempo"] + 5, i, f"{r['tempo']:.0f} meses", fontsize=10, va="center", color=TINTA)
    ax.set_yticks(range(len(c)), [f"{a}\n{r.get('Idade', '?')} anos" for a, r in c.iterrows()],
                  fontweight="bold")
    ax.set_ylim(-0.6, len(c) - 0.1)
    ax.set_xlim(0, c["tempo"].max() * 1.25)
    ax.set_xlabel("meses desde o diagnóstico")
    ax.set_title("Seguimento por paciente", loc="left", fontsize=13)
    ax.grid(axis="x", color=GRADE, lw=0.8)
    ax.set_axisbelow(True)
    ax.legend(handles=[Patch(color="#2a78d6", label="vivo / sem evento (seta = segue em acompanhamento)"),
                       Patch(color="#eb6834", label="evento (X = óbito)")],
              frameon=False, loc="lower right", bbox_to_anchor=(1, 1.0), ncol=2, fontsize=9)
    rodape(fig, "Com poucos pacientes não há comparação estatística de sobrevida. Dados simulados (aula).")
    salvar(fig, nome)


# ---------------------------------------------------------------------------
def main():
    arquivo_maf = sys.argv[1] if len(sys.argv) > 1 else "resultados/09_maf/coorte.maf"
    arquivo_clin = sys.argv[2] if len(sys.argv) > 2 else "dados/dados_clinicos.tsv"
    os.makedirs(SAIDA, exist_ok=True)

    df = ler_maf(arquivo_maf)
    clin = ler_clinico(arquivo_clin, sorted(df["Tumor_Sample_Barcode"].unique()))
    clin = clin[clin.index.isin(df["Tumor_Sample_Barcode"])]
    n_am = df["Tumor_Sample_Barcode"].nunique()

    pd.set_option("display.width", 200, "display.max_columns", 20)
    print(f"\n========== Resumo geral ==========\n{len(df)} variantes em {n_am} pacientes")
    print("\n========== Por paciente ==========")
    print(pd.crosstab(df["Tumor_Sample_Barcode"], df["Variant_Classification"], margins=True, margins_name="total"))
    print("\n========== Por gene ==========")
    print(df.groupby("Hugo_Symbol").agg(variantes=("Hugo_Symbol", "size"),
                                        pacientes=("Tumor_Sample_Barcode", "nunique"))
            .sort_values(["pacientes", "variantes"], ascending=False))
    comuns, raras = df[~df["rara"]], df[df["rara"]]
    print(f"\n========== Variantes comuns (gnomAD_AF >= {CORTE_AF}) ==========")
    print(comuns.groupby(["Hugo_Symbol", "HGVSp_Short", "gnomAD_AF"]).size().rename("n_pacientes").reset_index())
    print(f"\n========== Variantes raras (gnomAD_AF < {CORTE_AF} ou ausentes do gnomAD) ==========")
    print(raras[["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "HGVSp_Short",
                 "HGVSc", "gnomAD_AF", "VAF"]].sort_values("Tumor_Sample_Barcode").to_string(index=False))

    print("\nGerando figuras:")
    fig_resumo(df, "Resumo de todas as variantes", "01_resumo_todas_variantes.png")
    fig_oncoplot(df, clin, "Todas as variantes (inclui polimorfismos comuns)", "02_oncoplot_todas_variantes.png")
    if raras.empty:
        print("  Nenhuma variante rara — figuras 03, 04, 05 e 07 não geradas.")
    else:
        fig_resumo(raras, "Resumo das variantes raras (gnomAD < 1%)", "03_resumo_variantes_raras.png")
        fig_oncoplot(raras, clin, "Variantes raras (gnomAD < 1%)", "04_oncoplot_variantes_raras.png")
        ordem = {c: i for i, c in enumerate(CLASSES)}
        dedup = (raras.assign(_o=raras["Variant_Classification"].map(ordem))
                      .sort_values("_o").groupby(["Tumor_Sample_Barcode", "Hugo_Symbol"]).head(1))
        fig_oncoplot(dedup, clin, "Variantes raras, a mais deletéria por gene e paciente",
                     "05_oncoplot_sem_multihit.png")
    fig_titv(df, "06_titv.png")
    for gene in PROTEINAS:
        if gene in set(raras["Hugo_Symbol"]):
            fig_lollipop(raras, gene, n_am, f"07_lollipop_{gene}.png")
        else:
            print(f"  Nenhuma variante rara em {gene} — lollipop não gerado.")
    fig_painel_pacientes(df, clin, "08_painel_por_paciente.png")
    fig_vaf_gnomad(df, "09_vaf_vs_gnomad.png")
    fig_seguimento(df, clin, "10_seguimento.png")
    print(f"\nFiguras salvas em {SAIDA}/")


if __name__ == "__main__":
    main()
