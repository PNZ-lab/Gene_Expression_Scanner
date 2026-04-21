#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%% ===========================================================================
# 1. Description
# =============================================================================
'''
This script lets the user explore how the expression of genes are correlated with other genes or survival using several datasets such as Pölönen, TCGA, and CCLE
The purpose of this script is to explore correlations between genes in different datasets.
There are several applications of this script:

Replace 'target' in section 4 and run that cell
   - This script will create a waterfall graph with genes ranked based on Pearson's R
   - Breakpoints will be identified using the Kneedle algorithm - All genes with their R- and p-values will be written to a csv together with their position relative to the breakpoints
Replace 'target' and 'target2' in section 5 and run that cell
   - Script will create a graph and calculate Pearson's R and associated p-value for the two specified genes
Replace 'gene' and 'clin_col' in cell 6 and run that cell
   - Script will create a series of boxplots for the expression levels for patients separated by unique values in a column in the clinical dataset
Replace 'gene' in section 7 and run that cell
   - Script will generate a Kaplan-Meier graph for event-free survival for that gene
Run section 8 to create a Kaplan-Meier graph for each clinical parameter
Replace 'protein_x' and 'protein_y' in section 9 and run that cell
   - Script will create a graph and calculate Pearson's R and associated p-value for the two specified proteins across triplicates in six cell lines
Replace 'gene1' and 'gene2' in section 10 and run that cell
   - Script will create a graph and calculate Pearson's R and associated p-value for the two specified genes using cancer cell line data from CCLE
Replace 'gene' and 'group_by' in section 11 and run that cell
   -  Script will create a series of boxplots for the expression levels for cancer cell lines in CCLE separated by unique values in a column in the annotation dataset
Replace the 'gene' in Section 12 and run that cell
   - Script will generate a plot that shows the distribution of gene expression levels for that gene in the Pölönen dataset
Section 13 let's you investigate frequency of mutations across expression levels of certain genes'
Section 14 calculates the relative differences in expression between a subpopulation, ranks them and highlights certain genes
Section 15 plots expression levels across many cancers using TARGET+TCGA
Section 16 plots expression levels in healthy versus tunor
Section 17 plots every possible KM plot with high and low levels of a gene (based on median) and stratified by every clinical parameter in Pölönen


Sidenote: Some of the functionalities of this script can rely on KTC_functions.py to define a set of genes. This script supports you simply defining them by yourself, but if you want this functionality (e.g. Querying MSigDB for gene sets) - you can find KTC_functions.py on this GitHub.
'''

#%% ===========================================================================
# 2. Setup and settings
# =============================================================================

#Modules
import os
import mygene
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
from kneed import KneeLocator
import matplotlib.pyplot as plt
from adjustText import adjust_text
from itertools import combinations
import matplotlib.ticker as ticker
from collections import defaultdict
from scipy.signal import find_peaks
from scipy.stats import fisher_exact
from scipy.stats import gaussian_kde
from lifelines import KaplanMeierFitter
from KTC_functions import KTC_GetGeneSet
from KTC_functions import KTC_ncbi_gene_scraper
from lifelines.statistics import logrank_test
from scipy.stats import pearsonr, mannwhitneyu
from statannotations.Annotator import Annotator
from sklearn.linear_model import LinearRegression

files_directory = '/Volumes/kachrist-1/shares/cmgg_pnlab/Kasper/Data/Interesting_Lists' #Directory where files for clinical and gene expression are stored
# files_directory = '/Users/kachrist/Desktop/Interesting_Lists_Desktop' #Directory where files for clinical and gene expression are stored
out_dir         = '/Users/kachrist/Desktop/out_dir' # Directory where files and images are written. Subdirectories for individual genes are created

#Initialization
print('--- Loading data into memory ----')
print("Loading gene expression data...")
df_gexp                  = pd.read_csv(os.path.join(files_directory, 'PeCan_gexp.csv'))
print("Loading clinical data...")
df_annot                 = pd.read_csv(os.path.join(files_directory, 'PeCan_annot.csv'))
df_M0_clinical           = pd.read_csv(os.path.join(files_directory, 'PeCan_M0_clinical.csv'))
df_M1_classifying_driver = pd.read_csv(os.path.join(files_directory, "PeCan_M1_classifying_driver.csv"))
df_M2_ETP_status         = pd.read_csv(os.path.join(files_directory, "PeCan_M2_ETP_status.csv"))
df_M3_genetic_subtype    = pd.read_csv(os.path.join(files_directory, "PeCan_M3_genetic_subtype.csv"))
df_M3_subsubtype         = pd.read_csv(os.path.join(files_directory, "PeCan_M3_subsubtype.csv"))
df_M3_subtype            = pd.read_csv(os.path.join(files_directory, "PeCan_M3_subtype.csv"))
df_M4_pathway            = pd.read_csv(os.path.join(files_directory, "PeCan_M4_pathway.csv"))
df_M5_allesions_genes    = pd.read_csv(os.path.join(files_directory, "PeCan_M5_Allesions_genes.csv"))
df_M5_allesions_variants = pd.read_csv(os.path.join(files_directory, "PeCan_M5_Allesions_variants.csv"))
df_M7_IP                 = pd.read_csv(os.path.join(files_directory, "PeCan_M7_IP.csv"))
# df_ST1                   = pd.read_excel(os.path.join(files_directory, 'Polonen_Extended_Data.xlsx'), sheet_name='ST1_Clinical_Data')
df_ST3                   = pd.read_excel(os.path.join(files_directory, 'Polonen_Extended_Data.xlsx'), sheet_name='ST3_Sample_Annotations')
df_SNV                   = pd.read_excel(os.path.join(files_directory, 'Polonen_Extended_Data.xlsx'), sheet_name='ST14_Alterations.SNVIndel')
df_CNV                   = pd.read_excel(os.path.join(files_directory, 'Polonen_Extended_Data.xlsx'), sheet_name='ST19_Alterations.CNV')
df_alt_genes             = pd.read_excel(os.path.join(files_directory, 'Polonen_Extended_Data.xlsx'), sheet_name='ST10_Alterations.genes')

print("Loading cell line MS data...")
df_cell_line_MS          = pd.read_excel(os.path.join(files_directory, 'MS_results_PRC-5607 2.xlsx'), sheet_name='S2 Quantified proteins')
print("Loading CCLE data...")
# path_CCLE_rpkm = os.path.join(files_directory, 'CCLE_RNAseq_genes_rpkm_20180929.gct')
# path_CCLE_cl   = os.path.join(files_directory, 'Cell_lines_annotations_20181226.txt')
# df_CCLE_TPM   = pd.read_csv(path_CCLE_rpkm, sep='\t', skiprows=2)
# df_cl     = pd.read_csv(path_CCLE_cl, sep='\t')


df_cl = pd.read_csv(os.path.join(files_directory, 'DepMap_CellLine_annotation_2025Q3.csv'))
df_cl['OncotreeSubtype'] = df_cl['OncotreeSubtype'].replace({
    'Undifferentiated Pleomorphic Sarcoma/Malignant Fibrous Histiocytoma/High-Grade Spindle Cell Sarcoma':
    'UPS/MFHS/Spindle Cell Sarcoma'
})

df_cl['OncotreeSubtype'] = df_cl['OncotreeSubtype'].replace({
    'Uterine Carcinosarcoma/Uterine Malignant Mixed Mullerian Tumor':
    'Uterine Carcinosarcoma/Uterine MMMT'
})




df_CCLE_TPM      = pd.read_csv(os.path.join(files_directory, 'OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv'))

print("Merging clinical data...")
for df in [df_gexp, df_M1_classifying_driver, df_M2_ETP_status, df_M3_genetic_subtype, 
           df_M3_subsubtype, df_M3_subtype, df_M4_pathway, df_M5_allesions_genes, 
           df_M5_allesions_variants, df_M7_IP]:
    if 'Unnamed: 0' in df.columns:
        df.rename(columns={'Unnamed: 0': 'Gene'}, inplace=True)

def collapse_binary_columns(df, new_column_name):
    """ Collapses multiple binary indicator columns into a single categorical column. """
    binary_columns = df.columns[df.isin([0, 1]).all()]  # Identify binary indicator columns
    
    def get_category(row):
        for col in binary_columns:
            if row[col] == 1:
                return col  # Return the name of the column that is marked as 1
        return 'Unknown'  # If no column is marked, assign 'Unknown'
    
    if len(binary_columns) > 0:
        df[new_column_name] = df.apply(get_category, axis=1)
        df.drop(columns=binary_columns, inplace=True)  # Remove the original binary columns

# Apply to all relevant dataframes
collapse_binary_columns(df_M1_classifying_driver, "Classifying Driver")
collapse_binary_columns(df_M2_ETP_status, "ETP Status")
collapse_binary_columns(df_M3_genetic_subtype, "Genetic Subtype")
collapse_binary_columns(df_M3_subsubtype, "Subsubtype")
collapse_binary_columns(df_M3_subtype, "Subtype")
collapse_binary_columns(df_M4_pathway, "Pathway")
collapse_binary_columns(df_M5_allesions_genes, "Alleles Genes")
collapse_binary_columns(df_M5_allesions_variants, "Alleles Variants")
collapse_binary_columns(df_M7_IP, "IP Status")

# List all clinical dataframes for merging
df_list = [
    df_annot, df_M0_clinical, df_M1_classifying_driver, df_M2_ETP_status, 
    df_M3_genetic_subtype, df_M3_subsubtype, df_M3_subtype, df_M4_pathway, 
    df_M5_allesions_genes, df_M5_allesions_variants, df_M7_IP, df_ST3
]

# Convert patient identifiers into index for all dataframes
df_list = [df.set_index(df.columns[0]) for df in df_list]

# Merge all dataframes while handling overlapping columns
clin_df = df_list[0]
for i, df in enumerate(df_list[1:], start=1):
    clin_df = clin_df.join(df, how='outer', lsuffix='', rsuffix=f'_df{i}')

def pct_labels_from_fraction_bins(bins: np.ndarray) -> list[str]:
    # bins are e.g. [0.0, 0.1, ..., 1.0]
    return [f"{int(round(bins[i]*100))}-{int(round(bins[i+1]*100))}%"
            for i in range(len(bins) - 1)]

def pct_labels_from_percent_bins(bins: np.ndarray) -> list[str]:
    # bins are e.g. [0, 10, ..., 100]
    return [f"{int(bins[i])}-{int(bins[i+1])}%"
            for i in range(len(bins) - 1)]

bin_specs = {
    "Age.at.Diagnosis.in.Years": {
        "bins": np.arange(0, clin_df["Age.at.Diagnosis.in.Years"].max() + 5, 5),
        "labels_fn": lambda b: [f"{int(b[i])}-{int(b[i+1])}" for i in range(len(b)-1)],
        "suffix": "bin_5y",
    },
    # fraction (0-1)
    "Blast.percentage": {
        "bins": np.linspace(0, 1.0, 11),
        "labels_fn": pct_labels_from_fraction_bins,
        "suffix": "bin_10pct",
    },
    # percent (0-100)
    "Percent.Blasts.Tumor.Sample.Diagnostic": {
        "bins": np.arange(0, 101, 10),
        "labels_fn": pct_labels_from_percent_bins,
        "suffix": "bin_10pct",
    },
    # fraction (0-1)
    "Blast_percentage_TR": {
        "bins": np.linspace(0, 1.0, 11),
        "labels_fn": pct_labels_from_fraction_bins,
        "suffix": "bin_10pct",
    },
}

for col, spec in bin_specs.items():
    bins = spec["bins"]
    labels = spec["labels_fn"](bins)
    clin_df[f"{col}_{spec['suffix']}"] = pd.cut(
        clin_df[col],
        bins=bins,
        labels=labels,
        right=False,
        include_lowest=True,
    )

# Reset index to keep patient identifiers as a column
clin_df.reset_index(inplace=True)
clin_df.rename(columns={'index': 'Patient_ID'}, inplace=True)


print("Loading in TARGET + TCGA data...")
# === Load expression + metadata ===
df_target_tcga_tpm  = pd.read_csv(os.path.join(files_directory, "tcga_target_no_normal_rsem_gene_tpm.gz"), sep="\t", index_col=0)
df_target_tcga_meta = pd.read_csv(os.path.join(files_directory, "TCGA_TARGET_phenotype"), sep="\t", index_col=0)


# === Load gene map ===
df_target_tcga_map = pd.read_csv(
    os.path.join(files_directory, "gencode.v23.annotation.gene.probemap"),
    sep="\t", skiprows=1, header=None,
    names=["ensembl", "gene", "chrom", "start", "end", "strand"]
)

# === Prepare expression data ===
df_target_tcga_tpm = df_target_tcga_tpm.transpose()
df_target_tcga_tpm.columns = df_target_tcga_tpm.columns.str.split(".").str[0]
df_target_tcga_tpm = df_target_tcga_tpm.loc[:, df_target_tcga_tpm.columns.str.startswith("ENSG")]

# Map Ensembl to symbols
df_target_tcga_map["ensembl"] = df_target_tcga_map["ensembl"].str.split(".").str[0]
ensembl_to_symbol = dict(zip(df_target_tcga_map["ensembl"], df_target_tcga_map["gene"]))
df_target_tcga_tpm = df_target_tcga_tpm.rename(columns=ensembl_to_symbol)
df_target_tcga_tpm = df_target_tcga_tpm.loc[:, df_target_tcga_tpm.columns.notnull()]
df_target_tcga_tpm = df_target_tcga_tpm.T.groupby(level=0).mean().T  # collapse duplicated symbols

# === Merge with metadata ===
df_target_tcga_meta.index.name = "sample"
df_target_tcga_merged = df_target_tcga_tpm.join(df_target_tcga_meta, how="inner")

print('Loading in TCGA healthy v tumor')
path_TCGA_expr     = os.path.join(files_directory, 'tcga_toil_merged.parquet')
df_TCGA_expr       = pd.read_parquet(path_TCGA_expr)
path_TCGA_surviv   = os.path.join(files_directory, 'tcga_survival.txt')
df_TCGA_surviv     = pd.read_csv(path_TCGA_surviv, sep='\t')

df_TCGA_expr = df_TCGA_expr.copy()
df_TCGA_expr["PATIENT"] = df_TCGA_expr.index.str.slice(0,12)

df_surv_ready = df_TCGA_expr.merge(
    df_TCGA_surviv,
    left_on="PATIENT",
    right_on="_PATIENT",   # in your survival table
    how="inner"
)

print('--- Data loaded succesfully ----')



#Options
print_corr_genes = False # Genes above and below breakpoints can be written directly to console
write_files      = True # Turns on/off the writing of csv and pngs
log_scale        = False

mg = mygene.MyGeneInfo()

#%% ===========================================================================
# 3. Main functions
# =============================================================================

def WriteFile(name):
    plt.savefig(os.path.join(out_dir, name))
    print('file created: %s' %(os.path.join(out_dir, name)))

    print(name, target)
    out_dir_target = os.path.join(out_dir, target)
    if not os.path.exists(out_dir_target):
        os.makedirs(out_dir_target)
    if name.endswith('png') or name.endswith('svg'):
        plt.savefig(os.path.join(out_dir_target, name))
        print('file created: %s' %(os.path.join(out_dir_target, name)))

# This function is called elsewhere to create pairwise correlation plots between two genes
def Grapher(
    gene1, gene2,
    split_by_subtype=False,
    subanalysis_do=False, subanalysis_col=None, subanalysis_hit=None,
    show_equation=False, set_lim_0=False, pval_scientific=False,
    top_n_residuals=0,
    min_target_expr=None,  # NEW
    show_all_points=True   # NEW: show omitted in gray
):
    values1 = df_gexp.loc[df_gexp['Gene'] == gene1].iloc[0, 1:].astype(float).to_numpy()
    values2 = df_gexp.loc[df_gexp['Gene'] == gene2].iloc[0, 1:].astype(float).to_numpy()

    if log_scale:
        values1 = values1 + 1
        values2 = values2 + 1

    pecan_samples = df_gexp.columns[1:].tolist()

    # --- NEW: joint threshold mask applies to BOTH axes ---
    keep_mask_global = np.isfinite(values1) & np.isfinite(values2)
    if min_target_expr is not None:
        keep_mask_global &= (values1 >= float(min_target_expr))
        keep_mask_global &= (values2 >= float(min_target_expr))

    n_total = values1.size
    n_kept = int(keep_mask_global.sum())
    n_omitted = n_total - n_kept

    def _fmt_p(p):
        return '{:.2e}'.format(p) if pval_scientific else f'{p:.3f}'

    # Helper to compute pearson+fit on a subset mask (mask is over full sample axis)
    def _fit_on_mask(mask):
        x = values1[mask]
        y = values2[mask]
        # guard
        if x.size < 3:
            return None
        r, p = pearsonr(x, y)

        Xs = x.reshape(-1, 1)
        Ys = y
        model = LinearRegression()
        model.fit(Xs, Ys)
        a, b = model.coef_[0], model.intercept_

        x_line = np.linspace(np.nanmin(x), np.nanmax(x), 200)
        y_line = model.predict(x_line.reshape(-1, 1))
        return dict(r=r, p=p, a=a, b=b, x_line=x_line, y_line=y_line, x=x, y=y, model=model)

    # -------------------------------------------------------------------------
    # SPLIT BY SUBTYPE branch
    # -------------------------------------------------------------------------
    if split_by_subtype:
        match = clin_df[clin_df['Patient_ID'].isin(pecan_samples)]
        unique_subtypes = match['Classifying Driver'].dropna().unique()
        sample_subtypes = {row['Patient_ID']: row['Classifying Driver'] for _, row in match.iterrows()}

        for subtype in unique_subtypes:
            subtype_mask = np.array([sample_subtypes.get(s) == subtype for s in pecan_samples], dtype=bool)

            # Apply BOTH: subtype selection AND threshold filter
            mask = subtype_mask & keep_mask_global

            fig, ax = plt.subplots(figsize=(8, 8), dpi=100)

            # Plot omitted (in this subtype) if desired
            if show_all_points:
                omitted_mask = subtype_mask & (~keep_mask_global)
                if omitted_mask.any():
                    plt.scatter(values1[omitted_mask], values2[omitted_mask],
                                color='lightgray', alpha=0.15, s=12,
                                label='Omitted (below threshold)')

            # Plot included points
            plt.scatter(values1[mask], values2[mask], color='black', alpha=0.4, s=16,
                        label='Included in correlation')

            fit = _fit_on_mask(mask)
            if fit is None:
                plt.title(f'Polonen expression: {gene1} v {gene2}\nSubtype: {subtype}\n'
                          f'Not enough samples after filtering (kept {mask.sum()})',
                          fontsize=16)
                plt.xlabel(gene1 + ' Expression (VST)', fontsize=18)
                plt.ylabel(gene2 + ' Expression (VST)', fontsize=18)
                plt.legend(fontsize=12)
                WriteFile(f'Polonen_correlation_{gene1}_v_{gene2}_{subtype.replace("/","_")}.svg')
                plt.show()
                plt.close(fig)
                continue

            label_line = f'Filtered: R={fit["r"]:.2f}, p={_fmt_p(fit["p"])}'
            if show_equation:
                label_line += f', y={fit["a"]:.2f}x + {fit["b"]:.2f}'

            plt.plot(fit["x_line"], fit["y_line"], color='red', label=label_line)

            extra = ""
            if min_target_expr is not None:
                extra = f'\nThreshold: both ≥ {min_target_expr} (kept {mask.sum()}/{subtype_mask.sum()}, omitted {(subtype_mask.sum()-mask.sum())})'

            plt.title(f'Polonen expression: {gene1} v {gene2}\nSubtype: {subtype}{extra}', fontsize=16)
            plt.xlabel(gene1 + ' Expression (VST)', fontsize=18)
            plt.ylabel(gene2 + ' Expression (VST)', fontsize=18)
            plt.tick_params(axis='both', labelsize=16)

            if set_lim_0:
                plt.ylim(0)
                plt.xlim(0)

            plt.legend(fontsize=12)
            file_name = f'Polonen_correlation_{gene1}_v_{gene2}_{subtype.replace("/","_")}.svg'
            WriteFile(file_name)
            plt.show()
            plt.close(fig)

        return  # done with split_by_subtype

    # -------------------------------------------------------------------------
    # DEFAULT (non-subtype) branch
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 8), dpi=100)

    # Base scatter: omitted + included
    if show_all_points and (~keep_mask_global).any():
        plt.scatter(values1[~keep_mask_global], values2[~keep_mask_global],
                    color='lightgray', alpha=0.15, s=12,
                    label='Omitted (below threshold)')

    plt.scatter(values1[keep_mask_global], values2[keep_mask_global],
                color='black', alpha=0.25, s=16,
                label='Included in correlation')

    # Subanalysis: highlight subgroup points + add subgroup regression line (on filtered subset)
    etp_indices = []
    if subanalysis_do:
        for i, sample in enumerate(pecan_samples):
            match = clin_df[clin_df['Patient_ID'] == sample]
            if not match.empty:
                stage = match[subanalysis_col].values[0]
                if stage == subanalysis_hit:
                    etp_indices.append(i)

        if len(etp_indices) > 1:
            etp_mask = np.zeros(n_total, dtype=bool)
            etp_mask[etp_indices] = True

            # show subgroup points (but only those included by threshold if you want to be strict)
            etp_included = etp_mask & keep_mask_global
            if etp_included.any():
                plt.scatter(values1[etp_included], values2[etp_included],
                            color='red', alpha=0.5, s=22,
                            label=f'{subanalysis_hit} (included)')

            fit_etp = _fit_on_mask(etp_included)
            if fit_etp is not None:
                label_etp = f'{subanalysis_hit}: R={fit_etp["r"]:.2f}, p={_fmt_p(fit_etp["p"])}'
                if show_equation:
                    label_etp += f', y={fit_etp["a"]:.2f}x + {fit_etp["b"]:.2f}'
                plt.plot(fit_etp["x_line"], fit_etp["y_line"], color='red', label=label_etp)

    # Main fit on filtered set
    fit_all = _fit_on_mask(keep_mask_global)
    if fit_all is None:
        plt.title(f'Polonen expression: {gene1} v {gene2}\nNot enough samples after filtering (kept {n_kept}/{n_total})',
                  fontsize=16)
        plt.xlabel(gene1 + ' Expression (VST)', fontsize=18)
        plt.ylabel(gene2 + ' Expression (VST)', fontsize=18)
        plt.legend(fontsize=12)
        WriteFile(f'Polonen_correlation_{gene1}_v_{gene2}.svg')
        plt.show()
        return

    # Optional residual highlights (closest-to-line samples), but computed on filtered set
    if top_n_residuals and top_n_residuals > 0:
        x = fit_all["x"]
        y = fit_all["y"]
        y_pred = fit_all["model"].predict(x.reshape(-1, 1))
        residuals = np.abs(y - y_pred)
        sorted_indices = np.argsort(residuals)
        k = min(int(top_n_residuals), sorted_indices.size)

        # Map back to original indices for sample IDs
        kept_indices = np.where(keep_mask_global)[0]
        highlight_global_idx = kept_indices[sorted_indices[:k]]

        plt.scatter(values1[highlight_global_idx], values2[highlight_global_idx],
                    color='blue', edgecolor='white', s=80, label='Best correlating samples')

        closest_ids = [pecan_samples[i] for i in highlight_global_idx]
        print(f"Top {k} samples closest to regression line (filtered set):")
        print("Sample IDs:", ", ".join(closest_ids))

    label_all = f'Filtered: R={fit_all["r"]:.2f}, p={_fmt_p(fit_all["p"])}'
    if show_equation:
        label_all += f', y={fit_all["a"]:.2f}x + {fit_all["b"]:.2f}'
    plt.plot(fit_all["x_line"], fit_all["y_line"], label=label_all, color='black')

    plt.xlabel(gene1 + ' Expression (VST)', fontsize=18)
    plt.ylabel(gene2 + ' Expression (VST)', fontsize=18)
    plt.tick_params(axis='both', labelsize=16)

    extra = ""
    if min_target_expr is not None:
        extra = f"\nThreshold: both ≥ {min_target_expr} (kept {n_kept}/{n_total}, omitted {n_omitted})"
    plt.title(f'Polonen expression: {gene1} v {gene2}{extra}', fontsize=18)

    if set_lim_0:
        plt.ylim(0)
        plt.xlim(0)

    plt.legend(fontsize=12)
    file_name = f'Polonen_correlation_{gene1}_v_{gene2}.svg'
    WriteFile(file_name)
    plt.show()




# This plot is called to create the waterfall plot and return genes above and below breakpoints (for csv output)
from scipy.stats import mannwhitneyu, fisher_exact
import numpy as np
import matplotlib.pyplot as plt

def WaterfallPlot(
    dictionary,
    target_gene,
    gene_set,
    label,
    n_total=None, n_kept=None, n_omitted=None,
    min_target_expr=None,
    show_breakpoint=True,
    show_gene_set_ticks=True,
    tick_alpha=0.25,
    tick_lw=0.8,
    fisher_top_n=None,   # e.g. 200 to also test enrichment in top 200 ranks
):
    """
    Waterfall of correlations to target_gene.
    Adds:
      - red vertical ticks at ranks of genes in gene_set
      - stats: Mann-Whitney (gene_set vs background), Fisher enrichment in "above elbow",
        plus optional Fisher enrichment in top N
    """

    gene_set = set(gene_set or [])

    genes, r_p_values = zip(*dictionary)
    r_values = np.array([r for r, p in r_p_values], dtype=float)

    ranks_all = np.arange(1, len(genes) + 1)

    fig, ax = plt.subplots(figsize=(8, 5), dpi=100)

    finite_mask = np.isfinite(r_values)

    # Scatter all finite
    ax.scatter(
        ranks_all[finite_mask],
        r_values[finite_mask],
        color="black",
        s=2,
        zorder=2,
        label="gene expression correlations"
    )

    ax.set_xlabel("Rank of gene to gene correlation", fontsize=16)
    ax.set_ylabel("Pearson's R", fontsize=16)

    extra = ""
    if min_target_expr is not None and n_total is not None:
        extra = f"\nFiltered: {target_gene} >= {min_target_expr} (kept {n_kept}/{n_total}, omitted {n_omitted})"
    ax.set_title(f"{target_gene} : Waterfall plot of gene correlations{extra}", fontsize=18)

    # ---- knee finding on finite prefix ----
    finite_idx = np.where(finite_mask)[0]
    if finite_idx.size == 0:
        plt.show()
        WriteFile(f"Polonen_waterfall_{target_gene}.png")
        return [], []

    last_finite = finite_idx[-1]
    ranks = np.arange(1, last_finite + 2)
    r_vals = r_values[: last_finite + 1]

    kn_positive = KneeLocator(ranks, r_vals, curve="convex", direction="decreasing")
    kn_negative = KneeLocator(ranks, r_vals, curve="concave", direction="decreasing")

    knee_pos = kn_positive.knee
    knee_neg = kn_negative.knee

    if knee_pos is None:
        knee_pos = 0
    if knee_neg is None:
        knee_neg = len(r_vals)

    knee_pos = int(knee_pos)
    knee_neg = int(knee_neg)

    genes_above_elbow = [genes[i] for i in range(0, knee_pos)]
    genes_below_elbow = [genes[i] for i in range(knee_neg, len(genes))]

    # breakpoints
    if show_breakpoint:
        if knee_pos > 0:
            ax.axvline(x=knee_pos, color="blue", linestyle="--", label="breakpoint")
        if knee_neg < len(ranks_all):
            ax.axvline(x=knee_neg, color="blue", linestyle="--")

    # ---- NEW: gene_set tick marks (red vertical lines at gene ranks) ----
    if show_gene_set_ticks and len(gene_set) > 0:
        gene_to_rank = {g: i + 1 for i, g in enumerate(genes)}  # 1-based ranks
        set_ranks = sorted([gene_to_rank[g] for g in gene_set if g in gene_to_rank])

        # draw thin red lines across y-range
        if len(set_ranks) > 0:
            ymin, ymax = np.nanmin(r_values[finite_mask]), np.nanmax(r_values[finite_mask])
            for xr in set_ranks:
                ax.vlines(xr, ymin, ymax, colors="red", linewidth=tick_lw, alpha=tick_alpha, zorder=1)

            # add a legend handle (single proxy)
            ax.plot([], [], color="red", linewidth=2, label=f"{label} positions (n={len(set_ranks)})")

    # ---- NEW: stats ----
    # 1) Mann-Whitney: do gene_set genes have different r than background?
    set_rs = []
    bg_rs = []
    gene_set_present = set()

    for g, (r, p) in dictionary:
        if not np.isfinite(r):
            continue
        if g in gene_set:
            set_rs.append(r)
            gene_set_present.add(g)
        else:
            bg_rs.append(r)

    stats_lines = []
    if len(set_rs) >= 5 and len(bg_rs) >= 20:
        U, p_mwu = mannwhitneyu(set_rs, bg_rs, alternative="two-sided")
        stats_lines.append(f"MWU ({label} vs bg): p={p_mwu:.2e}, median_r_set={np.median(set_rs):.3f}, median_r_bg={np.median(bg_rs):.3f}")
    else:
        stats_lines.append(f"MWU: insufficient data (set n={len(set_rs)}, bg n={len(bg_rs)})")

    # 2) Fisher enrichment in above elbow
    # if knee_pos > 0:
    #     # top_set = set(genes_above_elbow)
    #     # in_top = len(top_set & gene_set_present)
    #     # not_in_top = len(top_set) - in_top
    #     # in_rest = len(gene_set_present) - in_top
    #     # not_in_rest = (len([g for g in genes[: last_finite + 1] if g not in gene_set_present])) - not_in_top

    #     # build 2x2 carefully
    #     # table = np.array([[in_top, not_in_top],
    #     #                   [in_rest, max(not_in_rest, 0)]], dtype=int)

    #     # only if valid
    #     # if table.min() >= 0 and table.sum() > 0 and table.shape == (2, 2):
    #     #     OR, p_f = fisher_exact(table)
    #     #     stats_lines.append(f"Fisher in above-elbow (k={len(top_set)}): OR={OR:.2f}, p={p_f:.2e}, hits={in_top}")
    # else:
    #     stats_lines.append("Fisher above-elbow: no elbow detected")

    # 3) Optional Fisher enrichment in top N ranks
    if fisher_top_n is not None:
        N = int(fisher_top_n)
        N = max(1, min(N, last_finite + 1))
        topN = set(genes[:N])
        in_topN = len(topN & gene_set_present)
        not_in_topN = len(topN) - in_topN
        in_restN = len(gene_set_present) - in_topN
        not_in_restN = (len([g for g in genes[: last_finite + 1] if g not in gene_set_present])) - not_in_topN

        tableN = np.array([[in_topN, not_in_topN],
                           [in_restN, max(not_in_restN, 0)]], dtype=int)
        if tableN.min() >= 0 and tableN.sum() > 0 and tableN.shape == (2, 2):
            ORN, pN = fisher_exact(tableN)
            stats_lines.append(f"Fisher top{N}: OR={ORN:.2f}, p={pN:.2e}, hits={in_topN}")

    # Put stats on plot
    ax.text(
        0.02, 0.98,
        "\n".join(stats_lines),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8)
    )

    # de-dup legend
    handles, labels_ = ax.get_legend_handles_labels()
    by_label = dict(zip(labels_, handles))
    ax.legend(by_label.values(), by_label.keys(), fontsize=9, loc="lower left")

    plt.show()
    WriteFile(f"Polonen_waterfall_{target_gene}.png")

    return genes_above_elbow, genes_below_elbow




def _safe_pearsonr(a: np.ndarray, b: np.ndarray):
    """Pearson r with NaN/constant guards; returns (r, p) or (np.nan, np.nan)."""
    m = np.isfinite(a) & np.isfinite(b)
    ax, bx = a[m], b[m]
    if ax.size < 2:
        return np.nan, np.nan
    if np.std(ax) == 0 or np.std(bx) == 0:
        return np.nan, np.nan
    return pearsonr(ax, bx)

#This function calculates R and p values for all genes to one gene, returns a ranked list for the Waterfall function, and writes the csv
def top_n_comparisons(
    gene,
    gene_set,
    label,
    min_target_expr=None,
    min_n_kept=25,
    top_n=10,
    write_files=False,
    out_dir=None,
    category_mode="elbows",  # "elbows" or "top_bottom"
):
    df_num = df_gexp.copy()
    df_num = df_num.drop_duplicates(subset='Gene').set_index('Gene')
    df_num = df_num.apply(pd.to_numeric, errors='coerce')

    if gene not in df_num.index:
        raise ValueError(f"Gene '{gene}' not found in df_gexp.")

    gene_values = df_num.loc[gene].to_numpy(dtype=float)
    r_dict = {}

    for other_gene in tqdm(df_num.index, desc='Comparing genes', unit='gene'):
        if other_gene == gene:
            continue

        other_values = df_num.loc[other_gene].to_numpy(dtype=float)

        keep_mask = np.isfinite(gene_values) & np.isfinite(other_values)
        if min_target_expr is not None:
            keep_mask &= (gene_values >= min_target_expr)
            keep_mask &= (other_values >= min_target_expr)

        if keep_mask.sum() < min_n_kept:
            r_dict[other_gene] = (np.nan, np.nan)
            continue

        r, p = pearsonr(gene_values[keep_mask], other_values[keep_mask])
        r_dict[other_gene] = (r, p)

    # Sort by r desc; NaNs last
    sorted_genes = sorted(
        r_dict.items(),
        key=lambda x: (-1e9 if np.isnan(x[1][0]) else x[1][0]),
        reverse=True
    )
    
    gsea_result = correlation_gsea_plot(
    sorted_genes=sorted_genes,
    gene_set=gene_set,
    target_gene=gene,
    label=label if label is not None else "PASINI_SUZ12_TARGETS_DN",
    weight=1.0
)

    # Waterfall + elbows
    genes_above_elbow, genes_below_elbow = WaterfallPlot(
        sorted_genes,
        gene,
        gene_set,
        label,
        min_target_expr=min_target_expr
    )

    # --- NEW: optional CSV export ---
    if write_files:
        if out_dir is None:
            raise ValueError("out_dir must be provided when write_files=True")

        out_dir_target = os.path.join(out_dir, gene)
        os.makedirs(out_dir_target, exist_ok=True)

        # precompute top/bottom sets if requested
        n_half = int(round(top_n / 2))
        top_set = set([g for g, _ in sorted_genes[:n_half]])
        bot_set = set([g for g, _ in sorted_genes[-n_half:]])

        rows = []
        for g, (r, p) in sorted_genes:
            if category_mode == "top_bottom":
                if g in top_set:
                    category = "top_n"
                elif g in bot_set:
                    category = "bottom_n"
                else:
                    category = "neither"
            else:  # "elbows"
                if g in genes_above_elbow:
                    category = "above_1st_elbow"
                elif g in genes_below_elbow:
                    category = "below_2nd_elbow"
                else:
                    category = "neither"

            rows.append([g, r, p, category])

        df_result = pd.DataFrame(rows, columns=["Gene", "Pearson_r", "p_value", "Category"])
        out_path = os.path.join(out_dir_target, f"{gene}_correlation_data.csv")
        df_result.to_csv(out_path, index=False)
        print(f"Wrote: {out_path}")

    n = int(round(top_n / 2))
    return sorted_genes, sorted_genes[:n], sorted_genes[-n:]


def correlation_gsea_plot(
    sorted_genes,
    gene_set,
    target_gene,
    label=None,
    weight=1.0,
    figsize=(9, 7),
    dpi=100,
    show_rank_metric=True,
    normalize_es=False,
):
    """
    Create a GSEA-like enrichment plot from a ranked correlation list.

    Parameters
    ----------
    sorted_genes : list of tuples
        Output like [(gene, (r, p)), ...], sorted descending by r.
    gene_set : iterable
        Genes of interest.
    target_gene : str
        The gene correlations were computed against.
    label : str or None
        Label for the gene set.
    weight : float
        Weighting exponent for hits. 
        0 = unweighted (classic KS-like)
        1 = weighted by |r|
    figsize : tuple
        Figure size.
    dpi : int
        Figure DPI.
    show_rank_metric : bool
        Whether to show the bottom panel with correlation values.
    normalize_es : bool
        If True, divides ES by max possible absolute excursion for display.
        This is just visual normalization, not NES.
    
    Returns
    -------
    result : dict
        Contains ES, peak index, hit indices, running ES, genes_in_set_ranked.
    """
    gene_set = set(gene_set or [])
    label = label if label is not None else "Gene set"

    # keep only finite correlations
    genes = []
    r_vals = []
    p_vals = []
    for g, (r, p) in sorted_genes:
        if np.isfinite(r):
            genes.append(g)
            r_vals.append(r)
            p_vals.append(p)

    genes = np.array(genes, dtype=object)
    r_vals = np.array(r_vals, dtype=float)
    p_vals = np.array(p_vals, dtype=float)

    N = len(genes)
    if N == 0:
        raise ValueError("No finite correlations available.")

    hits = np.array([g in gene_set for g in genes], dtype=bool)
    Nh = hits.sum()
    Nm = N - Nh

    if Nh == 0:
        raise ValueError("None of the genes in gene_set were found in sorted_genes.")

    # Ranking weights for hits
    hit_weights = np.abs(r_vals) ** weight
    sum_hit_weights = hit_weights[hits].sum()

    # Running enrichment score
    running_es = np.zeros(N, dtype=float)
    rs = 0.0

    for i in range(N):
        if hits[i]:
            rs += hit_weights[i] / sum_hit_weights if sum_hit_weights > 0 else 0.0
        else:
            rs -= 1.0 / Nm if Nm > 0 else 0.0
        running_es[i] = rs

    # Peak: choose the larger absolute excursion
    max_idx = np.argmax(running_es)
    min_idx = np.argmin(running_es)

    if abs(running_es[max_idx]) >= abs(running_es[min_idx]):
        es = running_es[max_idx]
        peak_idx = max_idx
    else:
        es = running_es[min_idx]
        peak_idx = min_idx

    if normalize_es:
        denom = np.max(np.abs(running_es))
        if denom > 0:
            running_es = running_es / denom
            es = es / denom

    hit_indices = np.where(hits)[0]

    # MWU as a simple companion stat
    set_rs = r_vals[hits]
    bg_rs = r_vals[~hits]
    mwu_text = "MWU: insufficient data"
    if len(set_rs) >= 5 and len(bg_rs) >= 20:
        U, p_mwu = mannwhitneyu(set_rs, bg_rs, alternative="two-sided")
        mwu_text = (
            f"MWU p={p_mwu:.2e}, "
            f"median_r_set={np.median(set_rs):.3f}, "
            f"median_r_bg={np.median(bg_rs):.3f}"
        )

    # plotting
    nrows = 3 if show_rank_metric else 2
    height_ratios = [3, 0.6, 1.8] if show_rank_metric else [4, 0.8]
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=1,
        figsize=figsize,
        dpi=dpi,
        sharex=True,
        gridspec_kw={"height_ratios": height_ratios}
    )

    if show_rank_metric:
        ax_es, ax_hits, ax_metric = axes
    else:
        ax_es, ax_hits = axes

    x = np.arange(N)

    # Top panel: running ES
    ax_es.plot(x, running_es, color="green", linewidth=2)
    ax_es.axhline(0, color="black", linestyle="--", linewidth=1)
    ax_es.axvline(peak_idx, color="red", linestyle="--", linewidth=1)
    ax_es.scatter([peak_idx], [running_es[peak_idx]], color="red", s=30, zorder=3)

    ax_es.set_ylabel("Running ES", fontsize=13)
    ax_es.set_title(
        f"{target_gene}: correlation-ranked enrichment for {label}",
        fontsize=16
    )

    stats_text = (
        f"ES={es:.3f}, peak_rank={peak_idx + 1}, hits={Nh}/{N}\n"
        f"{mwu_text}"
    )
    ax_es.text(
        0.02, 0.98,
        stats_text,
        transform=ax_es.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8)
    )

    # Middle panel: hit ticks
    for idx in hit_indices:
        ax_hits.vlines(idx, 0, 1, color="black", linewidth=0.8)
    ax_hits.set_ylim(0, 1)
    ax_hits.set_yticks([])
    ax_hits.set_ylabel("Hits", fontsize=11)

    # Bottom panel: ranking metric
    if show_rank_metric:
        ax_metric.plot(x, r_vals, color="black", linewidth=1.5)
        ax_metric.axhline(0, color="black", linestyle="--", linewidth=1)
        ax_metric.set_ylabel("Pearson r", fontsize=12)
        ax_metric.set_xlabel("Rank in correlation-sorted gene list", fontsize=13)

    else:
        ax_hits.set_xlabel("Rank in correlation-sorted gene list", fontsize=13)

    plt.tight_layout()
    plt.show()

    genes_in_set_ranked = genes[hits].tolist()

    return {
        "ES": es,
        "peak_idx": int(peak_idx),
        "peak_rank": int(peak_idx + 1),
        "running_es": running_es,
        "hit_indices": hit_indices,
        "genes_in_set_ranked": genes_in_set_ranked,
        "r_values": r_vals,
        "genes": genes,
    }



def SubsetBoxplotter(
    gene, PECAN_col, do_stats=True, write_file=False, _palette='gray',
    _dotcolor='white', _fontsize=14, order=None, set_ylim_0=False,
    list_n=False, sort_median=False, do_binary=False, hit_binary=None,
    mut_show=False, mut_gene=None, mut_aa=None,
    mut_mark='s', mut_col='red', mut_mark_s=30, stat_test='Mann-Whitney',
    stat_text='star', fig_size=(8,8), dpi=200, plot_type='boxplot', plt_show=True,
    median_col='black',
    # ---- NEW ----
    show_Q1Q3=True,
    q_bar_width=0.35,
    q_bar_lw=2.0,
    q_bar_color="green",
    q_bar_alpha=0.9,
    highlight_patients=None,
    highlight_label="Highlighted",
    highlight_color="red",
    highlight_marker="s",
    highlight_size=35

    # highlight_patients = set(map(str, highlight_patients or [])),
    # data["Patient_ID"] = data["Patient_ID"].astype(str)
):
    if gene not in df_gexp['Gene'].values:
        print(f"Gene '{gene}' not found in df_gexp.")
        return
    if PECAN_col not in clin_df.columns:
        print(f"Column '{PECAN_col}' not found in clin_df.")
        return

    matching_patient_ids = set(df_gexp.columns[1:]).intersection(clin_df['Patient_ID'].unique())
    filtered_df_gexp = df_gexp.loc[df_gexp['Gene'] == gene, list(matching_patient_ids)]
    if filtered_df_gexp.empty:
        print(f"No expression data found for gene {gene}.")
        return

    values = filtered_df_gexp.iloc[0].tolist()
    clin_lookup = clin_df.set_index('Patient_ID')

    # ========= BINARY MODE =========
    if do_binary and PECAN_col and hit_binary:
        group_labels = []
        for sample in matching_patient_ids:
            val = clin_lookup.loc[sample, PECAN_col] if sample in clin_lookup.index else None
            group = hit_binary if val == hit_binary else 'Other'
            group_labels.append(group)

        data = pd.DataFrame({
            'Subtype': group_labels,
            'Expression': values,
            'Patient_ID': list(matching_patient_ids)
        })

        data = data.dropna()
        data['Subtype'] = pd.Categorical(data['Subtype'], categories=[hit_binary, 'Other'], ordered=True)
        order = [hit_binary, 'Other']
        label_order = order if not list_n else [f"{x}\n(n={data['Subtype'].value_counts().get(x, 0)})" for x in order]
        if list_n:
            code_to_label = dict(zip(range(len(label_order)), label_order))
            data['Subtype_Labeled'] = data['Subtype'].cat.codes.map(code_to_label)
        else:
            data['Subtype_Labeled'] = data['Subtype']

    # ========= MULTI-GROUP MODE =========
    else:
        subtype_dict = defaultdict(list)
        for i, sample in enumerate(matching_patient_ids):
            if sample in clin_lookup.index:
                subtype = clin_lookup.loc[sample, PECAN_col]
                subtype_dict[subtype].append((sample, values[i]))

        records = []
        for subtype, sample_expr in subtype_dict.items():
            for sample, expr in sample_expr:
                records.append({'Subtype': subtype, 'Expression': expr, 'Patient_ID': sample})

        data = pd.DataFrame(records)
        data = data[~data['Subtype'].isin(['.', 'Unevaluable'])]
        data = data.dropna(subset=['Expression', 'Subtype'])

        #Clean
        data["Subtype"] = (
            data["Subtype"]
            .astype(str)
            .str.replace(":Significant_pathways", "", regex=False)
            .str.replace("B:SAMP:genetic_subtype_", "", regex=False)
            .str.replace("_Pathway", "", regex=False)
            
            # .str.replace("B:SAMP:genetic_subtype_, "", regex=False)
        )

        if sort_median:
            group_stats = data.groupby('Subtype')['Expression'].agg(['median', 'mean'])
            group_stats = group_stats.sort_values(by=['median', 'mean'], ascending=[False, False])
            print(group_stats)
            order = group_stats.index.tolist()
        elif order:
            data['Subtype'] = pd.Categorical(data['Subtype'], categories=order, ordered=True)

        sample_counts = data['Subtype'].value_counts().to_dict()
        if list_n:
            data['Subtype_Labeled'] = data['Subtype'].map(lambda x: f"{x}\n(n={sample_counts.get(x, 0)})")
        else:
            data['Subtype_Labeled'] = data['Subtype']

        label_order = data.groupby('Subtype_Labeled')['Expression'].median().sort_values(ascending=False).index.tolist() \
            if sort_median else (
            [f"{cat}\n(n={sample_counts.get(cat, 0)})" if list_n else cat for cat in order] if order else data['Subtype_Labeled'].unique())

    # ======== Mutation Lookup ========
    mutated_patients = set()
    if mut_show and mut_gene and mut_aa:
            print(f'Mutation analysis: gene={mut_gene}, aa_change={mut_aa}')
            mut_df = df_SNV[(df_SNV['gene'] == mut_gene) & (df_SNV['aa_change'] == mut_aa)]
            mutated_patients = set(mut_df['sample'])

    # ========== Plotting ==========
    plt.figure(figsize=fig_size, dpi=dpi)
    if plot_type == 'boxplot':
        ax = sns.boxplot(
            data=data, x='Subtype_Labeled', y='Expression',
            palette=_palette, showfliers=False, order=label_order
        )
    elif plot_type == 'violinplot':
        ax = sns.violinplot(
            data=data, x='Subtype_Labeled', y='Expression',
            palette=_palette, inner=None, order=label_order,
            cut=0, linewidth=1
        )
        # Draw median lines manually
        medians = data.groupby("Subtype_Labeled")["Expression"].median()
        for i, label in enumerate(label_order):
            median_val = medians[label]
            ax.hlines(
                y=median_val,
                xmin=i - 0.2, xmax=i + 0.2,
                color=median_col, linewidth=1.5, zorder=4
            )
    else:
        print('plot_type must be either "boxplot" or "violinplot"')

    # ---- NEW: Q1 / Q3 horizontal bars (works for both boxplot + violinplot) ----
    if show_Q1Q3:
        qs = (
            data.groupby("Subtype_Labeled")["Expression"]
            .quantile([0.25, 0.75])
            .unstack()
            .reindex(label_order)  # match plotted order
        )
        for i, label in enumerate(label_order):
            if label not in qs.index:
                continue
            q1 = qs.loc[label, 0.25]
            q3 = qs.loc[label, 0.75]
            # draw Q1 and Q3 as short horizontal bars
            ax.hlines(
                y=q1,
                xmin=i - q_bar_width/2, xmax=i + q_bar_width/2,
                color=q_bar_color, linewidth=q_bar_lw,
                alpha=q_bar_alpha, zorder=6
            )
            ax.hlines(
                y=q3,
                xmin=i - q_bar_width/2, xmax=i + q_bar_width/2,
                color=q_bar_color, linewidth=q_bar_lw,
                alpha=q_bar_alpha, zorder=6
            )

    # === Always overlay stripplot ===
    sns.stripplot(
        data=data[~data['Patient_ID'].isin(mutated_patients)],
        x='Subtype_Labeled', y='Expression',
        facecolor=_dotcolor, edgecolor='black',
        linewidth=0.8, alpha=0.3, jitter=True,
        order=label_order, zorder=2
    )
    
    highlight_patients = set(map(str, highlight_patients or []))
    data["Patient_ID"] = data["Patient_ID"].astype(str)
    
    if highlight_patients:
        hi = data[data["Patient_ID"].isin(highlight_patients)]
        if not hi.empty:
            sns.scatterplot(
                data=hi, x="Subtype_Labeled", y="Expression",
                marker=highlight_marker, s=highlight_size,
                color=highlight_color,
                edgecolor="black", linewidth=0.8,
                zorder=10, label=highlight_label
            )
            plt.legend()

    if mut_show:
        mut_data = data[data['Patient_ID'].isin(mutated_patients)]
        if not mut_data.empty:
            sns.scatterplot(data=mut_data, x='Subtype_Labeled', y='Expression',
                            marker=mut_mark, s=mut_mark_s, color=mut_col,
                            edgecolor='black', linewidth=0.8, zorder=10, label=f'{mut_gene} - {mut_aa}')
        plt.legend()

    if do_stats and data['Subtype'].nunique() >= 2:
        pairs = [(label_order[i], label_order[j])
                 for i in range(len(label_order)) for j in range(i+1, len(label_order))]

        annotator = Annotator(
            ax, pairs, data=data,
            x='Subtype_Labeled',
            y='Expression',
            order=label_order
        )

        cfg = dict(test=stat_test, text_format=stat_text, loc='outside',
                   line_height=0.02, line_offset=0.02, text_offset=0.01,
                   verbose=0)
        annotator.configure(**cfg)
        annotator.hide_non_significant = False
        annotator.apply_and_annotate()

        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax + 0.1 * (ymax - ymin))

    plt.title(f'Polonen expression of {gene} grouped by: \"{PECAN_col}\"', fontsize=_fontsize)
    plt.xlabel(PECAN_col, fontsize=_fontsize)
    plt.ylabel('Expression (VST)', fontsize=_fontsize)
    plt.xticks(rotation=90 if data['Subtype'].nunique() > 4 else 0, fontsize=_fontsize)
    plt.yticks(fontsize=_fontsize)
    if set_ylim_0:
        plt.ylim(0)
    plt.tight_layout()

    if write_file:
        out_path = os.path.join(out_dir, f"{gene}_{PECAN_col}_{_palette}")
        if mut_show:
            out_path = out_path + f'_{mut_gene}_{mut_aa}'
        if do_binary:
            out_path = out_path + f'_{hit_binary}'
        out_path = out_path + '.svg'
        print(f"Saved to: {out_path}")
        plt.savefig(out_path)
    if plt_show:
        plt.show()

def Grapher_MSpr1(protein1, protein2, df_msdataset, log_scale=False, show_equation=False, pval_scientific=False):
    # List of conditions (cell lines or experiments)
    conditions = ['ALLSIL', 'DND41', 'HPBALL', 'LOUCY', 'TALL1', 'KARPAS45']
    
    # Initialize a dictionary to store the color map for each condition
    colors = {'ALLSIL': 'red', 'DND41': 'blue', 'HPBALL': 'green', 'LOUCY': 'orange', 'TALL1': 'purple', 'KARPAS45': 'brown'}
    
    # Prepare lists to hold the data points for plotting
    values1 = []
    values2 = []
    color_values = []

    def find_matching_genes(gene_name, df):
        """Search for rows containing the gene_name as a substring in the 'Gene names' column."""
        # Exact match check
        matching_rows = df[df['Gene names'].str.contains(f"^{gene_name}$", case=False, na=False)]
        
        if not matching_rows.empty:
            return matching_rows  # Return exact matches if found
        
        # Find potential alternative matches that contain the gene_name string
        alternative_matches = df[df['Gene names'].str.contains(gene_name, case=False, na=False)]
        
        if alternative_matches.empty:
            print(f"No alternative matches found for '{gene_name}' either.")
        else:
            print("Suggested alternative gene names based on string matching:")
            for index, row in alternative_matches.iterrows():
                print(f"- {row['Gene names']}")  # Print all alternative gene names

        return alternative_matches  # Return the found alternative matches if any


    # Ensure proteins exist in the dataset, if not, find similar gene names
    protein1_matches = find_matching_genes(protein1, df_msdataset)
    protein2_matches = find_matching_genes(protein2, df_msdataset)

    # Track if any gene was matched through substring search
    substring_match_protein1 = not protein1_matches.empty and protein1 != protein1_matches['Gene names'].iloc[0]
    substring_match_protein2 = not protein2_matches.empty and protein2 != protein2_matches['Gene names'].iloc[0]

    if protein1_matches.empty:
        print(f"Warning: No exact match found for {protein1}.")
        return
    if protein2_matches.empty:
        print(f"Warning: No exact match found for {protein2}.")
        return

    

    # Inform user about substring matches
    if substring_match_protein1:
        print(f"Note: Substring match used for {protein1} in gene names.")
    if substring_match_protein2:
        print(f"Note: Substring match used for {protein2} in gene names.")
    
    # Get actual gene names for use in the title and axis labels
    actual_protein1 = protein1_matches['Gene names'].iloc[0] if not protein1_matches.empty else protein1
    actual_protein2 = protein2_matches['Gene names'].iloc[0] if not protein2_matches.empty else protein2

    # Loop through conditions
    for condition in conditions:
        # Loop through replicate numbers rep1, rep2, rep3
        for replicate in ['rep1', 'rep2', 'rep3']:
            # Extracting the values for each replicate and condition
            condition_values1 = protein1_matches.loc[
                protein1_matches['Gene names'].str.contains(protein1, case=False),
                f'log2_Reporter intensity corrected {condition}_{replicate}'
            ].values.flatten()
            
            condition_values2 = protein2_matches.loc[
                protein2_matches['Gene names'].str.contains(protein2, case=False),
                f'log2_Reporter intensity corrected {condition}_{replicate}'
            ].values.flatten()

            # Append all replicate values to the lists
            values1.extend(condition_values1)
            values2.extend(condition_values2)
            color_values.extend([colors[condition]] * len(condition_values1))  # Assign color based on condition

    # Remove any NaN values from values1 and values2
    mask = ~np.isnan(values1) & ~np.isnan(values2)
    values1 = np.array(values1)[mask]
    values2 = np.array(values2)[mask]
    color_values = np.array(color_values)[mask]

    # Scatter plot with different colors for different cell lines
    fig, ax = plt.subplots(figsize=(8, 8), dpi=200)
    for condition, color in colors.items():
        condition_mask = color_values == color
        plt.scatter(values1[condition_mask], values2[condition_mask], color=color, label=condition, alpha=0.5, zorder=6)

    # Linear regression model
    model = LinearRegression()
    model.fit(values1.reshape(-1, 1), values2)
    Y_pred = model.predict(values1.reshape(-1, 1))
    
    # Calculate Pearson's correlation coefficient and p-value
    r_value, p_value = pearsonr(values1, values2)
    formatted_p_value = '{:.2e}'.format(p_value) if pval_scientific else f'{p_value:.3f}'
    
    # Add regression line to the plot with the label
    plt.plot(values1, Y_pred, color='black', label=f'R: {r_value:.2f}, p={formatted_p_value}')

    # Set plot labels and title using actual gene names
    plt.xlabel(f'{actual_protein1} (log2 intensity)', fontsize=18)
    plt.ylabel(f'{actual_protein2} (log2 intensity)', fontsize=18)
    plt.tick_params(axis='both', labelsize=16)

    # Display the legend
    plt.title(f'MS Data Correlation: {actual_protein1} vs {actual_protein2}', fontsize=22)

    plt.legend(fontsize=16)

    # Annotate if substring matches were used
    if substring_match_protein1:
        plt.text(0.1, 0.8, f'Note: Substring match used for {protein1}', transform=plt.gca().transAxes, fontsize=14, color='red', bbox=dict(facecolor='white', alpha=0.7))
    if substring_match_protein2:
        plt.text(0.1, 0.75, f'Note: Substring match used for {protein2}', transform=plt.gca().transAxes, fontsize=14, color='red', bbox=dict(facecolor='white', alpha=0.7))

    # Save the figure as a .svg file
    plt.savefig(f'MS_Correlation_{actual_protein1}_vs_{actual_protein2}.svg')
    plt.show()


KM_COLOR_MAP = {
    # single-gene
    "Low":  "#1f77b4",   # blue
    "High": "#ff7f0e",   # orange

    # quantiles / bins
    "Q1": "#1f77b4",
    "Q2": "#aec7e8",
    "Q3": "#ffbb78",
    "Q4": "#ff7f0e",

    # two-gene combinations
    "Low/Low":   "#388234",
    "High/High": "#8e3b46",  # red (stands out)
    "High/Low":  "#648C87",
    "Low/High":  "#989670",
}

def KaplanMeier(
    _gene,
    gene2=None,
    n_groups=2,                       # used only if you don't pass split_value / split_values
    split_value=None,                 # <-- for n_groups==2: numeric cutoff for gene1 (High if > cutoff)
    split_value2=None,                # <-- for n_groups==2 with gene2: numeric cutoff for gene2
    split_values=None,                # <-- for n_groups>2: explicit bin edges for gene1 (list of edges)
    split_values2=None,               # <-- for n_groups>2 with gene2: explicit bin edges for gene2
    id_col="Patient_ID",
    time_col="EFS",
    event_col="EFS.status",
    KM_ymin_0=True,
    dpi=200,
    out_dir=None,
    filename=None
):
    """
    Split behavior
    --------------
    If gene2 is None:
      - If split_value is not None: median split replaced by cutoff split (High if > split_value else Low)
      - Else if split_values is not None: bin by explicit edges into B1..Bk
      - Else: median split (n_groups==2) or qcut (n_groups>2)

    If gene2 is provided:
      - For n_groups==2: each gene split by its cutoff (split_value / split_value2) if provided,
        otherwise median; then combined into High/High, High/Low, Low/High, Low/Low.
      - For n_groups>2: each gene binned using split_values / split_values2 if provided,
        otherwise qcut; then combined into Gene1Bin/Gene2Bin groups.

    Notes
    -----
    - Cutoff rule matches your original median split semantics: High is strictly > cutoff, else Low.
    - For explicit edges, bins are right-inclusive by default: (a, b] via pandas.cut.
    """

    import numpy as np
    import pandas as pd
    import os
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test, multivariate_logrank_test

    # ---- helpers ----
    def _get_expr_series(gene, ids):
        s = df_gexp.set_index("Gene").loc[gene, ids]
        if not hasattr(s, "index"):
            s = pd.Series({i: float(s) for i in ids})
        return s.astype(float)

    def _bin_by_cutoff(series, cutoff, prefix=""):
        # High if > cutoff, else Low (keeps your original median-split behavior)
        grp = np.where(series > cutoff, f"{prefix}High", f"{prefix}Low")
        grp = pd.Series(grp, index=series.index)
        labels = [f"{prefix}Low", f"{prefix}High"]
        title = f"High vs Low (cutoff={cutoff:g})"
        return grp, labels, title

    def _bin_by_edges(series, edges, prefix=""):
        # edges must be monotonically increasing; creates len(edges)-1 bins
        edges = list(edges)
        if len(edges) < 2:
            raise ValueError("split_values must contain at least 2 edges (e.g., [0, 1, 2]).")

        # Create labels B1..B(k)
        n_bins = len(edges) - 1
        labels = [f"{prefix}B{i+1}" for i in range(n_bins)]

        grp = pd.cut(series, bins=edges, labels=labels, include_lowest=True, right=True)
        grp = grp.astype(str)

        title = f"Custom bins (edges={edges})"
        return grp, labels, title

    def _bin_expression(series, n_groups, split_value=None, split_values=None, prefix=""):
        # Priority: explicit cutoff -> explicit edges -> median/qcut
        if split_value is not None:
            return _bin_by_cutoff(series, split_value, prefix=prefix)

        if split_values is not None:
            return _bin_by_edges(series, split_values, prefix=prefix)

        # Default behavior (your old logic)
        if n_groups == 2:
            med = series.median()
            return _bin_by_cutoff(series, med, prefix=prefix[:-0] if prefix else "")
        else:
            labels = [f"{prefix}Q{i+1}" for i in range(n_groups)]
            try:
                grp = pd.qcut(series, q=n_groups, labels=labels)
                if grp.nunique() < n_groups:
                    raise ValueError("Collapsed bins in qcut; using rank-based qcut.")
            except Exception:
                ranked = series.rank(method="first")
                grp = pd.qcut(ranked, q=n_groups, labels=labels)
            title = f"{n_groups} quantiles"
            return grp.astype(str), labels, title

    # ---- 1) Align clinical and expression data ----
    matched = clin_df[clin_df[id_col].isin(df_gexp.columns)].copy()

    expr1 = _get_expr_series(_gene, matched[id_col])
    matched["Expression"] = matched[id_col].map(expr1)

    if gene2 is not None:
        expr2 = _get_expr_series(gene2, matched[id_col])
        matched["Expression2"] = matched[id_col].map(expr2)
        drop_cols = ["Expression", "Expression2", time_col, event_col]
    else:
        drop_cols = ["Expression", time_col, event_col]

    matched = matched.dropna(subset=drop_cols).copy()
    matched[event_col] = matched[event_col].astype(int)

    # ---- 2) Build groups ----
    if gene2 is None:
        # If user passes split_values, override n_groups to match edges
        if split_values is not None:
            n_groups_eff = len(split_values) - 1
        else:
            n_groups_eff = n_groups

        matched["KM_Group"], group_order, title_part = _bin_expression(
            matched["Expression"],
            n_groups=n_groups_eff,
            split_value=split_value if n_groups_eff == 2 else None,   # cutoff intended for 2-group
            split_values=split_values if n_groups_eff != 2 else None, # edges intended for multi-bin
            prefix=""
        )

    else:
        # determine effective n_groups per gene if custom edges passed
        if split_values is not None:
            n1 = len(split_values) - 1
        else:
            n1 = n_groups

        if split_values2 is not None:
            n2 = len(split_values2) - 1
        else:
            n2 = n_groups

        # For the 4-way High/Low combo you asked about, you typically want n_groups==2 on both genes.
        g1, g1_order, t1 = _bin_expression(
            matched["Expression"],
            n_groups=n1,
            split_value=split_value if n1 == 2 else None,
            split_values=split_values if n1 != 2 else None,
            prefix=""
        )
        g2, g2_order, t2 = _bin_expression(
            matched["Expression2"],
            n_groups=n2,
            split_value=split_value2 if n2 == 2 else None,
            split_values=split_values2 if n2 != 2 else None,
            prefix=""
        )

        matched["Gene1_Group"] = g1.astype(str)
        matched["Gene2_Group"] = g2.astype(str)
        matched["KM_Group"] = matched["Gene1_Group"] + "/" + matched["Gene2_Group"]

        group_order = [f"{a}/{b}" for a in g1_order for b in g2_order]

        if n1 == 2 and n2 == 2:
            title_part = f"{_gene} x {gene2}: High/Low combinations"
            if split_value is not None or split_value2 is not None:
                title_part += f" (cutoffs: {split_value if split_value is not None else 'median'} / {split_value2 if split_value2 is not None else 'median'})"
            else:
                title_part += " (median splits)"
        else:
            title_part = f"{_gene} x {gene2}: custom bin combinations"

    # If custom edges created out-of-range NaNs (pd.cut), drop them (keeps plot sane)
    matched = matched[matched["KM_Group"].notna() & (matched["KM_Group"] != "nan")].copy()

    # ---- 3) Compute p-value ----
    # Only do pairwise logrank when it's exactly Low vs High and gene2 is None
    if gene2 is None and group_order == ["Low", "High"]:
        g_low = matched.query("KM_Group == 'Low'")
        g_high = matched.query("KM_Group == 'High'")
        p_value = logrank_test(
            g_high[time_col], g_low[time_col],
            event_observed_A=g_high[event_col],
            event_observed_B=g_low[event_col]
        ).p_value
    else:
        mlr = multivariate_logrank_test(
            matched[time_col], matched["KM_Group"], matched[event_col]
        )
        p_value = float(mlr.p_value)

    # ---- 4) Plot ----
    plt.figure(figsize=(6, 5), dpi=dpi)
    ax = None
    for grp in group_order:
        sub = matched.loc[matched["KM_Group"] == grp]
        if sub.empty:
            continue
        kmf = KaplanMeierFitter()
        label = f"{grp} (n={len(sub)})"
        color = KM_COLOR_MAP.get(grp, None)  # fallback to default if unseen label
        kmf.fit(
            durations=sub[time_col],
            event_observed=sub[event_col],
            label=label
        )
        
        ax = kmf.plot_survival_function(
            ci_show=True,
            ax=ax,
            color=color
        )

    if KM_ymin_0:
        plt.ylim(0, 1)
    plt.xlabel("Days (Event-free survival)")
    plt.ylabel("Survival Probability")
    plt.title(
        f"{_gene}" + (f" + {gene2}" if gene2 is not None else "") +
        f" — Kaplan–Meier by expression: {title_part}\nlog-rank p = {p_value:.4g}"
    )
    plt.legend()
    plt.grid(True, alpha=0.3)

    # ---- 5) Save if requested ----
    if out_dir is not None:
        if filename is None:
            if gene2 is None:
                if split_value is not None:
                    suffix = f"cutoff{split_value:g}"
                elif split_values is not None:
                    suffix = "custombins"
                else:
                    suffix = "quartiles" if n_groups == 4 else (f"{n_groups}quantiles" if n_groups != 2 else "median")
                filename = f"{_gene}_KaplanMeier_{suffix}.svg"
            else:
                if (split_value is not None) or (split_value2 is not None):
                    suffix = "2gene_cutoffs"
                elif (split_values is not None) or (split_values2 is not None):
                    suffix = "2gene_custombins"
                else:
                    suffix = "2gene_defaultbins"
                filename = f"{_gene}_{gene2}_KaplanMeier_{suffix}.svg"

        out_path = os.path.join(out_dir, filename)
        try:
            WriteFile(out_path)
        except NameError:
            plt.savefig(out_path, bbox_inches="tight")

    plt.show()

    # ---- 6) Return ----
    if gene2 is None:
        return p_value, matched[[id_col, "Expression", "KM_Group"]].copy()
    else:
        return p_value, matched[[id_col, "Expression", "Expression2", "Gene1_Group", "Gene2_Group", "KM_Group"]].copy()




from lifelines.statistics import multivariate_logrank_test

def KaplanMeier_clinical(clin_column):
    # Drop NA in selected column
    df = clin_df[clin_df['Patient_ID'].isin(df_gexp.columns) & clin_df[clin_column].notna()].copy()

    # Ensure EFS and EFS.status are available and valid
    df = df[df['EFS'].notna() & df['EFS.status'].notna()]

    unique_groups = df[clin_column].unique()
    if len(unique_groups) < 2:
        print(f"Skipping {clin_column} (only one group present)")
        return

    kmf = KaplanMeierFitter()
    plt.figure(figsize=(7, 5), dpi=200)
    plt.xlabel("Days (Event-Free Survival)", fontsize=11)
    plt.ylabel("Survival Probability", fontsize=11)

    # Dynamically generate a color palette
    colors = sns.color_palette("tab20", n_colors=len(unique_groups))

    for color, group in zip(colors, unique_groups):
        mask = df[clin_column] == group
        n_patients = mask.sum()
        label = f"{group} (n={n_patients})"
        kmf.fit(df.loc[mask, 'EFS'], df.loc[mask, 'EFS.status'], label=label)
        kmf.plot_survival_function(ci_show=False, color=color)

    # Perform multivariate log-rank test
    results = multivariate_logrank_test(df['EFS'], df[clin_column], df['EFS.status'])
    pval = results.p_value

    plt.grid()
    plt.legend(title=clin_column, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.title(f"KM Survival by {clin_column} (p={pval:.4f})", fontsize=13)
    WriteFile(os.path.join(out_dir, f'{clin_column}_KM.svg'))
    plt.show()

def Grapher_CCLE(
    gene_x, gene_y,
    show_equation=False, log_scale=False, set_lim_0=False,
    filter_col=None, filter_val=None,
    label_points=False,
    df_CCLE_TPM=df_CCLE_TPM,
    df_cl=df_cl,
):
    """
    Adapted to the new data structure:
      - df_CCLE_TPM: rows = models, columns include "ModelID" and gene columns like "MEIS1 (4211)"
      - df_cl: metadata with "ModelID" and cell line name columns (e.g., StrippedCellLineName, CCLEName, CellLineName)

    Functionality preserved:
      - optional filter on df_cl column/value
      - optional log transform (log10(TPM+1)) when log_scale=True
      - Pearson r + linear regression line
      - optional axis limits from 0
      - optional point labels using adjust_text
    """

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.stats import pearsonr
    from sklearn.linear_model import LinearRegression

    # -----------------------------
    # 1) Gene column lookup (new format: "GENE (ENTREZ)")
    # -----------------------------
    def _build_gene_map(columns):
        # Map "MEIS1" -> "MEIS1 (4211)" for columns that look like that
        gm = {}
        for c in columns:
            if isinstance(c, str) and " (" in c and c.endswith(")"):
                gm[c.split(" (")[0]] = c
        return gm

    gene_map = _build_gene_map(df_CCLE_TPM.columns)

    # Allow passing exact column name too
    if gene_x in df_CCLE_TPM.columns:
        gene_col_x = gene_x
    elif gene_x in gene_map:
        gene_col_x = gene_map[gene_x]
    else:
        raise ValueError(f"Gene {gene_x} not found in expression data columns.")

    if gene_y in df_CCLE_TPM.columns:
        gene_col_y = gene_y
    elif gene_y in gene_map:
        gene_col_y = gene_map[gene_y]
    else:
        raise ValueError(f"Gene {gene_y} not found in expression data columns.")

    # -----------------------------
    # 2) Filter cell lines based on annotation (new ID: ModelID)
    # -----------------------------
    if filter_col and (filter_val is not None):
        # guard against trailing comma bug: ('OncotreeSubtype',)
        if isinstance(filter_col, (tuple, list)) and len(filter_col) == 1:
            filter_col = filter_col[0]

        filtered_cl = df_cl[df_cl[filter_col] == filter_val].copy()
    else:
        filtered_cl = df_cl.copy()

    # Model IDs to keep
    if "ModelID" not in filtered_cl.columns:
        raise ValueError("df_cl must contain 'ModelID' column.")
    model_ids = filtered_cl["ModelID"].astype(str).values

    # -----------------------------
    # 3) Subset expression (merge on ModelID)
    # -----------------------------
    if "ModelID" not in df_CCLE_TPM.columns:
        raise ValueError("df_CCLE_TPM must contain 'ModelID' column.")

    expr_sub = df_CCLE_TPM[["ModelID", gene_col_x, gene_col_y]].copy()
    expr_sub["ModelID"] = expr_sub["ModelID"].astype(str)

    # keep only filtered models
    expr_sub = expr_sub[expr_sub["ModelID"].isin(model_ids)].copy()

    # numeric + drop NA
    expr_sub[gene_col_x] = pd.to_numeric(expr_sub[gene_col_x], errors="coerce")
    expr_sub[gene_col_y] = pd.to_numeric(expr_sub[gene_col_y], errors="coerce")
    expr_sub = expr_sub.dropna(subset=[gene_col_x, gene_col_y])

    # keep a "CCLE_ID-like" index for labeling/iterrows parity
    expr_sub = expr_sub.set_index("ModelID")

    # match original naming downstream (df_subset columns should be gene symbols the user passed)
    df_subset = expr_sub[[gene_col_x, gene_col_y]].copy()
    df_subset.columns = [gene_x, gene_y]

    # -----------------------------
    # 4) Log transform if specified
    # -----------------------------
    if log_scale:
        with np.errstate(divide="ignore", invalid="ignore"):
            df_subset[gene_x] = np.log10(df_subset[gene_x] + 1)
            df_subset[gene_y] = np.log10(df_subset[gene_y] + 1)
        ylabel = "log10(TPM + 1)"
    else:
        ylabel = "TPM"

    # -----------------------------
    # 5) Regression
    # -----------------------------
    x = df_subset[gene_x].values
    y = df_subset[gene_y].values

    if len(x) > 1:
        r_value, p_value = pearsonr(x, y)
        model = LinearRegression()
        model.fit(x.reshape(-1, 1), y)
        x_range = np.linspace(np.min(x), np.max(x), 100)
        y_pred = model.predict(x_range.reshape(-1, 1))
        reg_label = f"R={r_value:.2f}, p={p_value:.2e}"
        if show_equation:
            reg_label += f"\n y = {model.coef_[0]:.2f} x + {model.intercept_:.2f}"
    else:
        r_value, p_value, reg_label = (np.nan, np.nan, "Insufficient data")

    # -----------------------------
    # 6) Plotting
    # -----------------------------
    plt.figure(figsize=(8, 8), dpi=150)
    sns.scatterplot(x=gene_x, y=gene_y, data=df_subset, alpha=0.2, edgecolor=None, color="black")

    if len(x) > 1:
        plt.plot(x_range, y_pred, color="black", label=reg_label)

    title = f"{gene_x} vs {gene_y} - CCLE"
    if filter_col and (filter_val is not None):
        title += f"\n({filter_col} = {filter_val})"

    plt.xlabel(f"{gene_x} Expression ({ylabel})", fontsize=18)
    plt.ylabel(f"{gene_y} Expression ({ylabel})", fontsize=18)
    plt.title(title, fontsize=22)
    plt.legend(fontsize=14)
    plt.tick_params(axis="both", labelsize=16)

    if set_lim_0:
        plt.xlim(left=0)
        plt.ylim(bottom=0)

    # -----------------------------
    # 7) Label points if specified
    # -----------------------------
    if label_points:
        from adjustText import adjust_text
        texts = []

        # Create a mapping from ModelID to a name (choose the most likely available name column)
        name_col_candidates = ["StrippedCellLineName", "CCLEName", "CellLineName", "ModelID"]
        name_col = next((c for c in name_col_candidates if c in df_cl.columns), "ModelID")

        id_to_name = df_cl.set_index(df_cl["ModelID"].astype(str))[name_col].to_dict()

        for i, row in df_subset.iterrows():
            label = id_to_name.get(str(i), str(i))  # fallback to ModelID
            texts.append(plt.text(row[gene_x], row[gene_y], label, fontsize=12))

        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey'), forcePoints=0.5)

    plt.tight_layout()
    plt.show()

def CCLE_Boxplotter(
    gene,
    group_by='OncotreeSubtype',
    log_scale=False,
    fig_height=8,
    fig_width=8,
    palette='gray',
    do_stats=False,
    min_samples=5,
    ylim=0,
    add_scatter=False,
    include_terms=None,       # list of specific group terms to plot
    highlight_labels=None     # list of group names to make red
):
    """
    Adapted to the new data structure:
      - df_CCLE_TPM: rows = models, columns include 'ModelID' and gene columns like 'MEIS1 (4211)'
      - df_cl: metadata with 'ModelID' and grouping columns (e.g. OncotreeSubtype, OncotreePrimaryDisease, etc.)

    Functionality preserved:
      - include_terms filter
      - min_samples filter
      - optional log10(TPM+1)
      - group ordering by median
      - optional pairwise t-test annotations (statannotations)
      - optional highlighting of specific x tick labels
    """

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from itertools import combinations

    # -----------------------------
    # 1) Gene column lookup (new format: "GENE (ENTREZ)")
    # -----------------------------
    gene_map = {
        c.split(" (")[0]: c
        for c in df_CCLE_TPM.columns
        if isinstance(c, str) and " (" in c and c.endswith(")")
    }

    # Allow passing exact column name too
    if gene in df_CCLE_TPM.columns:
        gene_col = gene
    elif gene in gene_map:
        gene_col = gene_map[gene]
    else:
        raise ValueError(f"{gene} not found in TPM columns.")

    # -----------------------------
    # 2) Construct dataframe with metadata (merge on ModelID)
    # -----------------------------
    if "ModelID" not in df_CCLE_TPM.columns:
        raise ValueError("df_CCLE_TPM must contain 'ModelID' column.")
    if "ModelID" not in df_cl.columns:
        raise ValueError("df_cl must contain 'ModelID' column.")
    if group_by not in df_cl.columns:
        raise ValueError(f"group_by '{group_by}' not found in df_cl.")

    plot_data = df_CCLE_TPM[["ModelID", gene_col]].copy()
    plot_data["ModelID"] = plot_data["ModelID"].astype(str)
    plot_data[gene_col] = pd.to_numeric(plot_data[gene_col], errors="coerce")
    plot_data = plot_data.dropna(subset=[gene_col])

    meta = df_cl[["ModelID", group_by]].copy()
    meta["ModelID"] = meta["ModelID"].astype(str)

    plot_data = plot_data.merge(meta, on="ModelID", how="left")
    plot_data = plot_data.dropna(subset=[group_by])

    # Rename to preserve downstream behavior
    plot_data = plot_data.rename(columns={gene_col: "Expression"})

    # -----------------------------
    # 3) Optional filter to only specified terms
    # -----------------------------
    if include_terms is not None:
        plot_data = plot_data[plot_data[group_by].isin(include_terms)]

    # -----------------------------
    # 4) Filter out groups with fewer than min_samples
    # -----------------------------
    group_counts = plot_data[group_by].value_counts()
    valid_groups = group_counts[group_counts >= min_samples].index
    plot_data = plot_data[plot_data[group_by].isin(valid_groups)]

    if plot_data.empty:
        raise ValueError("No data left to plot after filtering (include_terms/min_samples).")

    # -----------------------------
    # 5) Log-transform if needed
    # -----------------------------
    if log_scale:
        plot_data["Expression"] = np.log10(plot_data["Expression"] + 1)
        ylabel = "log10(TPM + 1)"
    else:
        ylabel = "TPM"

    # -----------------------------
    # 6) Order groups by median expression
    # -----------------------------
    group_order = (
        plot_data.groupby(group_by)["Expression"]
        .median()
        .sort_values(ascending=False)
        .index.tolist()
    )

    # -----------------------------
    # 7) Plot
    # -----------------------------
    plt.figure(figsize=(fig_width, fig_height), dpi=300)
    ax = sns.boxplot(
        x=group_by,
        y="Expression",
        data=plot_data,
        order=group_order,
        palette=palette
    )

    if add_scatter:
       sns.stripplot(
           x=group_by,
           y="Expression",
           data=plot_data,
           order=group_order,
           color="black",
           size=3,
           jitter=0.25,
           alpha=0.6,
           ax=ax
       )
       
    # -----------------------------
    # 8) Optional statistical annotations (unchanged behavior)
    # -----------------------------
    if do_stats and plot_data[group_by].nunique() >= 2:
        # expects statannotations to be installed and Annotator imported/available
        pairs = list(combinations(group_order, 2))
        annotator = Annotator(ax, pairs, data=plot_data, x=group_by, y="Expression")
        annotator.configure(test="t-test_ind", text_format="star", loc="inside", verbose=0)
        annotator.hide_non_significant = True
        annotator.apply_and_annotate()

    # -----------------------------
    # 9) Highlight tick labels if requested
    # -----------------------------
    if highlight_labels:
        for label in ax.get_xticklabels():
            if label.get_text() in highlight_labels:
                label.set_color("red")
                label.set_fontweight("bold")

    # -----------------------------
    # 10) Final plot settings
    # -----------------------------
    plt.title(f"{gene} Expression by {group_by}", fontsize=22)
    plt.xlabel(group_by, fontsize=18)
    plt.ylabel(ylabel, fontsize=18)
    plt.xticks(rotation=90)
    plt.tick_params(axis="both", labelsize=16)
    plt.ylim(bottom=ylim)
    plt.tight_layout()
    plt.show()

def Grapher_CCLE_density(
    gene,
    df_CCLE_TPM=df_CCLE_TPM,
    df_cl=df_cl,
    filter_col=None,
    filter_val=None,
    label_col="StrippedCellLineName",   # nicer than ModelID usually
    log2p1=True,
    kde_bw_adjust=1.0,
    label_points=True,
    figsize=(8, 5)
):

    # -----------------------------
    # 1. Gene column lookup
    # -----------------------------
    gene_map = {
        c.split(" (")[0]: c
        for c in df_CCLE_TPM.columns
        if "(" in c
    }

    if gene not in gene_map:
        raise ValueError(f"{gene} not found in TPM columns.")

    gene_col = gene_map[gene]

    # -----------------------------
    # 2. Filter metadata
    # -----------------------------
    meta = df_cl.copy()

    if filter_col and filter_val:
        meta = meta[meta[filter_col] == filter_val]

    # -----------------------------
    # 3. Merge on ModelID
    # -----------------------------
    df = (
        df_CCLE_TPM[["ModelID", gene_col]]
        .merge(meta[["ModelID", label_col]], on="ModelID")
        .dropna()
    )

    if df.empty:
        raise ValueError("No overlapping models after filtering.")

    x = df[gene_col].astype(float).values

    if log2p1:
        x_plot = np.log2(x + 1)
        xlabel = f"{gene} log2(TPM+1)"
    else:
        x_plot = x
        xlabel = f"{gene} TPM"

    # -----------------------------
    # 4. Plot
    # -----------------------------
    plt.figure(figsize=figsize, dpi=200)

    sns.kdeplot(x=x_plot, fill=True, bw_adjust=kde_bw_adjust)

    
    # --- Create beeswarm-ish y positions by binning x and stacking in y ---
    ax = plt.gca()
    ymax = ax.get_ylim()[1]
    
    # Parameters controlling vertical spread
    band_low  = 0.05 * ymax
    band_high = 0.95 * ymax
    levels_per_bin = 10          # max vertical levels used within each x-bin
    x_bins = 14                  # number of x bins across the range
    
    # Bin x
    x_min, x_max = np.min(x_plot), np.max(x_plot)
    edges = np.linspace(x_min, x_max, x_bins + 1)
    bin_id = np.clip(np.digitize(x_plot, edges) - 1, 0, x_bins - 1)
    
    # For each bin, assign y levels 0..k-1 (then wrap if >levels_per_bin)
    rng = np.random.default_rng(0)
    y_pts = np.zeros_like(x_plot, dtype=float)
    
    for b in range(x_bins):
        idx = np.where(bin_id == b)[0]
        if len(idx) == 0:
            continue
    
        # randomize order so labels don't all align
        rng.shuffle(idx)
    
        # assign levels
        k = len(idx)
        levels = np.arange(k) % levels_per_bin
    
        # map levels to band (even spacing)
        if levels_per_bin == 1:
            y_vals = np.full(k, (band_low + band_high) / 2)
        else:
            y_vals = band_low + (levels / (levels_per_bin - 1)) * (band_high - band_low)
    
        # add a tiny jitter so it looks organic
        y_vals = y_vals + rng.normal(0, 0.01 * ymax, size=k)
    
        y_pts[idx] = y_vals
    
    plt.scatter(x_plot, y_pts, s=55, alpha=0.85, zorder=5)
    plt.ylim(0, ymax)  # keep density scale
    
    plt.scatter(x_plot, y_pts, s=55, alpha=0.85, zorder=5)
    plt.ylim(0, ymax)  # keep density scale

    if label_points:
        texts = []
        for xi, yi, lab in zip(x_plot, y_pts, df[label_col]):
            texts.append(plt.text(xi, yi, str(lab), fontsize=9))

        adjust_text(
            texts,
            arrowprops=dict(arrowstyle='-', color='grey', lw=0.8)
        )

    plt.axvline(0, color="black", lw=1)

    title = f"{gene} expression distribution (CCLE)"
    if filter_col and filter_val:
        title += f"\nfiltered by {filter_col} = {filter_val}"

    plt.title(title)
    plt.xlabel(xlabel)
    plt.yticks([])
    plt.tight_layout()
    plt.show()

    return df

def Mutation_Barplotter(gene_mut, gene_split, cutoff, mut_aa=None, plot_abs=False):
    _fontsize = 18

    # Expression status
    expr_df = df_gexp[df_gexp['Gene'] == gene_split].drop('Gene', axis=1).T
    expr_df.columns = [gene_split]
    expr_df['Patient_ID'] = expr_df.index
    expr_df[f'{gene_split}_Status'] = [
        'High' if x > cutoff else 'Low' for x in expr_df[gene_split]
    ]

    global mut_df

    # Mutation filtering logic
    if mut_aa is None:
        mut_df = df_SNV[(df_SNV['gene'] == gene_mut) & (df_SNV['aa_change'] != '.')]
        mutation_label = f'mut {gene_mut}'
    else:
        mutation_string = f'p.{mut_aa}'
        mut_df = df_SNV[(df_SNV['gene'] == gene_mut) & (df_SNV['aa_change'] == mutation_string)]
        mutation_label = f'{gene_mut} $^{{{mut_aa}}}$'

    # Get unique patient IDs with mutation
    mutated_patients = set(mut_df['sample'].unique())

    expr_df['Has_Mutation'] = expr_df['Patient_ID'].apply(
        lambda x: 'Mutated' if x in mutated_patients else 'Wildtype'
    )

    # Contingency table and Fisher's test
    contingency_table = pd.crosstab(expr_df[f'{gene_split}_Status'], expr_df['Has_Mutation'])
    contingency_table = contingency_table.reindex(['Low', 'High'])

    oddsratio, pvalue = fisher_exact(contingency_table)

    # Plot data
    plot_data = contingency_table if plot_abs else contingency_table.div(contingency_table.sum(axis=1), axis=0) * 100

    # Plotting
    fig, ax = plt.subplots(figsize=(6, 6), dpi=200)
    plot_data[['Mutated', 'Wildtype']].plot(
        kind='bar',
        stacked=True,
        ax=ax,
        color=['#EF8784', '#D4D4D4'],
        edgecolor='None',
        width=0.8
    )

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    for spine in ['left', 'bottom']:
        ax.spines[spine].set_linewidth(3)
        ax.spines[spine].set_color('black')

    ax.set_ylabel('%' if not plot_abs else 'Count', fontsize=_fontsize, fontweight='bold')
    plot_data = plot_data.reindex(['Low', 'High'])

    ax.set_xticklabels([
        fr'{gene_split}$^{{low}}$', 
        fr'{gene_split}$^{{high}}$'
    ], rotation=-45, fontweight='bold', fontsize=_fontsize)
    for label in ax.get_xticklabels():
        label.set_ha('left')
    ax.tick_params(axis='y', labelsize=_fontsize)

    ax.tick_params(
        axis='both', which='major', direction='out', length=10, width=3,
        colors='black', bottom=True, top=False, left=True, right=False
    )
    ax.tick_params(
        axis='y', which='minor', direction='out', length=5, width=3,
        colors='black', bottom=True, top=False, left=True, right=False
    )
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

    # Legend
    handles, _ = ax.get_legend_handles_labels()
    handles = handles[::-1]
    custom_labels = [f'WT {gene_mut}', mutation_label]
    legend = ax.legend(
        handles=handles,
        labels=custom_labels,
        loc='center left',
        bbox_to_anchor=(1, 0.5),
        handlelength=2.5,
        handleheight=2.5,
        borderpad=0.5,
        frameon=False
    )
    plt.setp(legend.get_texts(), fontweight='bold', fontsize=_fontsize)

    ax.minorticks_on()
    if not plot_abs:
        ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(10))
    else:
        ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    plt.tight_layout()

    # Bracket and p-value
    y_vals = plot_data.sum(axis=1)
    y_max = y_vals.max()
    bracket_height = y_max + (5 if not plot_abs else 1)
    text_height = bracket_height + (2 if not plot_abs else 0.5)
    x1, x2 = 0, 1

    ax.plot([x1, x1, x2, x2], [bracket_height, bracket_height + 1, bracket_height + 1, bracket_height], color='black', lw=3)

    def pval_to_stars(p):
        if p < 0.001: return '***'
        elif p < 0.01: return '**'
        elif p < 0.05: return '*'
        else: return 'ns'

    p_text = pval_to_stars(pvalue)
    ax.text((x1 + x2) / 2, text_height, p_text, ha='center', va='bottom', fontsize=_fontsize*1.5, fontweight='bold')

    ax.set_xlabel('')
    WriteFile(f'barplot_{gene_mut}_{mut_aa}_{plot_abs}.svg')
    plt.show()

    # Print contingency table
    print(f"\nContingency Table - gene: {gene_mut}, specific mutation: {mut_aa}")
    print(contingency_table)
    print(f"Odds Ratio: {oddsratio:.3f}, P-value: {pvalue:.4g}")

def Plot_Density(gene, write_file=False):
    if gene not in df_gexp['Gene'].values:
        print(f"Gene '{gene}' not found in df_gexp.")
        return

    expr_values = df_gexp.loc[df_gexp['Gene'] == gene].iloc[0, 1:].dropna().astype(float)

    # Kernel density estimate
    kde = gaussian_kde(expr_values)
    x_grid = np.linspace(expr_values.min() - 1, expr_values.max() + 1, 1000)
    y_kde = kde(x_grid)

    # Local maxima and minima
    peak_indices, _ = find_peaks(y_kde)
    trough_indices, _ = find_peaks(-y_kde)  # local minima

    peak_x = x_grid[peak_indices]
    peak_y = y_kde[peak_indices]

    trough_x = x_grid[trough_indices]
    trough_y = y_kde[trough_indices]

    # Per-point density for rug coloring
    point_densities = kde(expr_values)

    # Normalize densities to [0,1] for color mapping
    norm = plt.Normalize(point_densities.min(), point_densities.max())
    cmap = plt.cm.plasma_r

    # Begin plot
    plt.figure(figsize=(8, 5), dpi=300)

    # KDE plot
    plt.plot(x_grid, y_kde, color='gray', linewidth=2, label='Density')
    plt.fill_between(x_grid, y_kde, color='gray', alpha=0.3)

    # Custom rug plot: colored vertical lines based on local density
    for x, d in zip(expr_values, point_densities):
        plt.axvline(x, ymin=0, ymax=0.02, color=cmap(norm(d)), linewidth=1)

    # Local maxima (peaks)
    plt.scatter(peak_x, peak_y, color='black', s=50, marker='^', label='Local Maxima', zorder=5)
    for x, y in zip(peak_x, peak_y):
        plt.text(x, y + 0.01, f'{x:.2f}', ha='center', va='bottom', fontsize=10, color='black')

    # Local minima (troughs)
    plt.scatter(trough_x, trough_y, color='black', s=50, marker='v', label='Local Minima', zorder=5)
    for x, y in zip(trough_x, trough_y):
        plt.text(x, y - 0.01, f'{x:.2f}', ha='center', va='top', fontsize=10, color='black')

    median_val = np.median(expr_values)
    plt.axvline(median_val, color='red', linestyle='--', linewidth=2, label='Median')
    plt.text(median_val*1.02, max(y_kde)*0.95, f'Median\n{median_val:.2f}', 
             ha='left', va='top', fontsize=10, color='red')

    # Labels and legend
    plt.title(f'Kernel density and rug plot for {gene} expression', fontsize=14)
    plt.xlabel('Expression (VST)', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.legend()
    plt.tight_layout()
    if write_file:
        WriteFile(os.path.join(out_dir, f'{gene}_Density.svg'))
    else:
        plt.show()


import numpy as np

def TCGA_TARGET_expression(
    gene,
    filter_diseases=None,
    figsize=(14, 6),
    point_alpha=0.2,
    point_size=10
):
    if gene not in df_target_tcga_merged.columns:
        print(f"{gene} not found.")
        return

    subset = df_target_tcga_merged[[gene, "_study", "_primary_disease"]].dropna()

    if filter_diseases is not None:
        subset = subset[subset["_primary_disease"].isin(filter_diseases)]
        if subset.empty:
            print(f"No rows left after filtering to {filter_diseases}")
            return

    # order categories by median expression
    ordered = subset.groupby("_primary_disease")[gene].median().sort_values(ascending=False).index

    # draw the boxplots (dodge by cohort)
    plt.figure(figsize=figsize, dpi=200)
    ax = sns.boxplot(
        data=subset,
        x="_primary_disease",
        y=gene,
        hue="_study",
        order=ordered,
        fliersize=0,
        linewidth=1
    )

    # --- perfectly aligned black dots ---
    # map x positions for each disease category
    x_pos = {cat: i for i, cat in enumerate(ordered)}

    # iterate per disease and cohort; choose offset:
    # - if both TCGA and TARGET present: offsets = (-0.2, +0.2)
    # - if only one cohort present: offset = 0 (centered under the single box)
    for cat in ordered:
        cat_df = subset[subset["_primary_disease"] == cat]
        cohorts_present = cat_df["_study"].unique().tolist()
        both = set(cohorts_present) == {"TCGA", "TARGET"}
        offsets = {"TCGA": -0.2, "TARGET": 0.2} if both else {"TCGA": 0.0, "TARGET": 0.0}

        for cohort in cohorts_present:
            d = cat_df[cat_df["_study"] == cohort]
            # small horizontal jitter so points don't sit on one x
            jitter = (np.random.rand(len(d)) - 0.5) * 0.12
            xs = x_pos[cat] + offsets[cohort] + jitter
            ax.scatter(
                xs, d[gene].values,
                s=point_size, alpha=point_alpha, c="black", edgecolors="none", zorder=3
            )

    # tidy up
    plt.xticks(rotation=90)
    plt.title(f"{gene} expression across TCGA and TARGET tumors")
    plt.ylabel("log2(TPM + 0.001)" if subset[gene].max() < 100 else "TPM")
    plt.xlabel("Cancer type")
    # keep the boxplot legend (cohorts). It’s already correct.
    plt.tight_layout()
    plt.show()




def TCGA_healthy_v_tumor(
    expr_annot,
    gene,                       # ENSG id or symbol (e.g. "FTO")
    filter_col="_primary_site",
    filter_value=None,
    normals="auto",             # 'auto' | 'tcga' | 'gtex' | 'both'
    add_points=True,
    add_pvalue=True,
    dpi=200,
    point_color='black',
    point_alpha=0.2
):
    # --- resolve gene ---
    display_name = gene
    if not gene.startswith("ENSG"):
        # query mygene if not an Ensembl ID
        q = mg.query(gene, scopes="symbol", fields="ensembl.gene", species="human")
        if not q["hits"]:
            raise KeyError(f"Could not resolve symbol '{gene}' via mygene.")
        ens = q["hits"][0].get("ensembl", {})
        if isinstance(ens, list):
            ens = ens[0]
        ensg = ens.get("gene")
        if ensg is None:
            raise KeyError(f"No Ensembl ID found for '{gene}'")
    else:
        ensg = gene

    if ensg not in expr_annot.columns:
        raise KeyError(f"{ensg} not in merged DataFrame columns.")

    # --- subset dataframe ---
    df = expr_annot[[ensg, "_sample_type", "_study", filter_col]].rename(columns={ensg: "expr"}).copy()

    # apply filter
    if filter_value and str(filter_value).lower() != "all":
        if isinstance(filter_value, (list, tuple, set)):
            m = False
            for v in filter_value:
                m = m | df[filter_col].astype(str).str.contains(str(v), case=False, na=False)
        else:
            m = df[filter_col].astype(str).str.contains(str(filter_value), case=False, na=False)
        df = df.loc[m]

    # define groups
    is_tumor  = (df["_study"] == "TCGA") & (df["_sample_type"] == "Primary Tumor")
    is_tcnorm = (df["_study"] == "TCGA") & (df["_sample_type"].isin(["Solid Tissue Normal", "Blood Derived Normal"]))
    is_gtnorm = (df["_study"] == "GTEX")

    if normals == "tcga":
        df = df.loc[is_tumor | is_tcnorm].copy()
        df["Group"] = np.where(is_tumor.loc[df.index], "Tumor", "Normal (TCGA)")
        order = ["Tumor", "Normal (TCGA)"]
    elif normals == "gtex":
        df = df.loc[is_tumor | is_gtnorm].copy()
        df["Group"] = np.where(is_tumor.loc[df.index], "Tumor", "Normal (GTEx)")
        order = ["Tumor", "Normal (GTEx)"]
    elif normals == "both":
        df = df.loc[is_tumor | is_tcnorm | is_gtnorm].copy()
        df["Group"] = np.where(is_tumor.loc[df.index], "Tumor",
                        np.where(is_tcnorm.loc[df.index], "Normal (TCGA)", "Normal (GTEx)"))
        order = ["Tumor", "Normal (TCGA)", "Normal (GTEx)"]
    elif normals == "auto":
        if is_tcnorm.sum() > 0:
            df = df.loc[is_tumor | is_tcnorm].copy()
            df["Group"] = np.where(is_tumor.loc[df.index], "Tumor", "Normal (TCGA)")
            order = ["Tumor", "Normal (TCGA)"]
        else:
            df = df.loc[is_tumor | is_gtnorm].copy()
            df["Group"] = np.where(is_tumor.loc[df.index], "Tumor", "Normal (GTEx)")
            order = ["Tumor", "Normal (GTEx)"]
    else:
        raise ValueError("normals must be one of {'auto','tcga','gtex','both'}")

    df = df.dropna(subset=["expr", "Group"])
    if df["Group"].nunique() < 2:
        raise ValueError("Not enough groups after filtering to draw a boxplot.")

    # --- plot ---
    plt.figure(dpi=dpi, figsize=(6,6))
    ax = sns.boxplot(data=df, x="Group", y="expr", order=order, showfliers=False)
    if add_points:
        sns.stripplot(data=df, x="Group", y="expr", order=order, alpha=point_alpha, size=3, color=point_color)

    subtitle = ""
    if add_pvalue:
        tumor_vals = df.loc[df["Group"] == "Tumor", "expr"].values
        norm_label = next((g for g in order if g != "Tumor" and g in df["Group"].unique()), None)
        if norm_label is not None:
            norm_vals = df.loc[df["Group"] == norm_label, "expr"].values
            if len(tumor_vals) > 0 and len(norm_vals) > 0:
                _, p = mannwhitneyu(tumor_vals, norm_vals, alternative="two-sided")
                subtitle = f"\nMWU p = {p:.2e} (Tumor vs {norm_label})"

    filt_txt = "All" if not filter_value or str(filter_value).lower()=="all" else str(filter_value)
    ax.set_title(f"{display_name} — {filter_col}: {filt_txt}{subtitle}")
    ax.set_ylabel("Expression (log2, DESeq2-normalized)")
    plt.tight_layout()
    plt.show()

    return df[["expr", "Group", "_sample_type", "_study", filter_col]].copy()

def TCGA_km_plot(
    df_surv_ready: pd.DataFrame,
    gene: str,
    endpoint: str = "OS",                      # one of: OS, DSS, PFI, DFI
    filter_col: str = "primary disease or tissue",
    filter_value: str = "Acute Myeloid Leukemia",
    n_groups: int = 4,                         # e.g., 2 (median), 3 tertiles, 4 quartiles
    min_per_group: int = 10,                   # ensure each KM group has enough samples
    match: str = "exact",                      # 'exact' or 'contains' (case-insensitive)
    dpi: int = 200,
    y_min_0 = True
) -> pd.DataFrame:
    """
    Returns the tidy DF used for KM (columns: expr, time, event, KM_Group).
    """
    ensg = _symbol_to_ensg(gene)
    time_col, event_col = f"{endpoint}.time", endpoint
    missing = [c for c in [ensg, time_col, event_col, filter_col] if c not in df_surv_ready.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    # --- subset to one cancer/diagnosis ---
    col = df_surv_ready[filter_col].astype(str)
    if match == "contains":
        mask = col.str.contains(str(filter_value), case=False, na=False)
    else:  # exact
        mask = col.str.casefold() == str(filter_value).casefold()
    df = df_surv_ready.loc[mask, [ensg, time_col, event_col]].copy()
    if df.empty:
        raise ValueError(f"No rows match {filter_col} == '{filter_value}' (match='{match}').")

    # --- coerce & clean ---
    df.rename(columns={ensg: "expr", time_col: "time", event_col: "event"}, inplace=True)
    df["time"]  = pd.to_numeric(df["time"], errors="coerce")
    df["event"] = pd.to_numeric(df["event"], errors="coerce")
    df = df.dropna(subset=["expr", "time", "event"])
    df = df.loc[df["time"] > 0]
    if df.empty:
        raise ValueError("No data left after cleaning (expr/time/event).")

    # --- make quantile groups robust to ties ---
    r = df["expr"].rank(method="average")
    df["KM_Group"] = pd.qcut(
        r, q=n_groups,
        labels=[f"Q{i}" for i in range(1, n_groups+1)],
        duplicates="drop"
    )
    # ensure at least 2 groups and min_per_group
    vc = df["KM_Group"].value_counts()
    # If any group too small, try falling back to median split
    if (vc.min() < min_per_group) or (df["KM_Group"].nunique() < 2):
        # fallback: median split
        med = df["expr"].median()
        df["KM_Group"] = np.where(df["expr"] > med, "High", "Low")

    if df["KM_Group"].nunique() < 2:
        raise ValueError("Not enough distinct groups after binning; try smaller n_groups or a broader filter.")

    # --- plot ---
    plt.figure(dpi=dpi)
    kmf = KaplanMeierFitter()
    order = sorted(df["KM_Group"].unique(), key=lambda x: str(x))
    for grp in order:
        gdf = df.loc[df["KM_Group"] == grp]
        kmf.fit(gdf["time"], gdf["event"], label=f"{grp} (n={len(gdf)})")
        kmf.plot_survival_function(ci_show=True)

    # --- log-rank p-value ---
    if df["KM_Group"].nunique() == 2:
        g1, g2 = order[:2]
        A, B = df[df["KM_Group"] == g1], df[df["KM_Group"] == g2]
        p = logrank_test(A["time"], B["time"], A["event"], B["event"]).p_value
    else:
        p = multivariate_logrank_test(df["time"], df["KM_Group"], df["event"]).p_value

    if y_min_0 == True:
        plt.ylim(0)

    plt.title(f"{endpoint} - {gene} ({ensg})\n{filter_col} = {filter_value}\nlog-rank p={p:.2e}")
    plt.xlabel("Time (days)")
    plt.ylabel("Survival probability")
    plt.tight_layout()
    plt.show()

    return df[["expr", "time", "event", "KM_Group"]].copy()

def _symbol_to_ensg(symbol: str) -> str:
    if symbol.upper().startswith("ENSG"):
        return symbol
    q = mg.query(symbol, scopes="symbol", fields="ensembl.gene", species="human")
    if not q.get("hits"):
        raise KeyError(f"Could not resolve symbol '{symbol}' to Ensembl ID")
    ens = q["hits"][0].get("ensembl", {})
    if isinstance(ens, list):
        ens = ens[0]
    ensg = ens.get("gene")
    if not ensg:
        raise KeyError(f"No Ensembl ID found for '{symbol}'")
    return ensg

def KaplanMeier_expression_clinical(
    _gene,
    group_col="Classifying Driver",   # choose any clinical column here
    n_groups=2,                       # still supports 2 (median), 4 (quartiles), etc.
    id_col="Patient_ID",
    time_col="EFS",
    event_col="EFS.status",
    KM_ymin_0=True,
    dpi=100,                          # your default preference
    out_dir=None,                     # if given, will save figures (one per group)
    filename_prefix=None,             # optional custom filename prefix
    min_group_size=5                  # skip groups with too few matched samples
):
    """
    Returns:
        dict_pvals : {group_value: p_value}
        dict_tables: {group_value: DataFrame[id_col, Expression, KM_Group]}
    """
    # Get which groups we will iterate over
    all_groups = (
        clin_df[group_col]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )
    all_groups = sorted(all_groups)

    dict_pvals  = {}
    dict_tables = {}

    # Prepare expression vector only once (across all samples), then subset per group
    # Align clinical IDs present in expression columns
    # (Assumes df_gexp has "Gene" rows and sample IDs as columns)
    all_match = clin_df[clin_df[id_col].isin(df_gexp.columns)].copy()
    try:
        expr_full = df_gexp.set_index("Gene").loc[_gene, all_match[id_col]].astype(float)
    except KeyError:
        raise KeyError(f"Gene '{_gene}' not found in df_gexp['Gene'].")

    # Loop groups
    for grp in all_groups:
        sub = clin_df.loc[clin_df[group_col].astype(str) == grp].copy()
        sub = sub[sub[id_col].isin(df_gexp.columns)]

        if len(sub) < min_group_size:
            # Not enough total samples for this group; skip
            continue

        # Map expression for this group
        expr_series = sub[id_col].map(expr_full)
        sub = sub.assign(Expression=expr_series).dropna(subset=["Expression", time_col, event_col]).copy()
        if len(sub) < min_group_size:
            # Not enough after NA filtering; skip
            continue

        sub[event_col] = sub[event_col].astype(int)

        # Build groups (High/Low by within-group median or quantiles)
        if n_groups == 2:
            med = sub["Expression"].median()
            sub["KM_Group"] = np.where(sub["Expression"] > med, "High", "Low")
            group_order = ["Low", "High"]
            title_part = "High vs Low (median split)"
        else:
            labels = [f"Q{i+1}" for i in range(n_groups)]
            try:
                sub["KM_Group"] = pd.qcut(sub["Expression"], q=n_groups, labels=labels)
                if sub["KM_Group"].nunique() < n_groups:
                    raise ValueError("Collapsed bins in qcut; using rank-based qcut.")
            except Exception:
                ranked = sub["Expression"].rank(method="first")
                sub["KM_Group"] = pd.qcut(ranked, q=n_groups, labels=labels)
            group_order = labels
            title_part = f"{n_groups} quantiles"

        # Compute p-value
        if n_groups == 2:
            g_low  = sub.query("KM_Group == 'Low'")
            g_high = sub.query("KM_Group == 'High'")
            if g_low.empty or g_high.empty:
                # Cannot compute log-rank with missing bins
                continue
            p_value = logrank_test(
                g_high[time_col], g_low[time_col],
                event_observed_A=g_high[event_col],
                event_observed_B=g_low[event_col]
            ).p_value
        else:
            mlr = multivariate_logrank_test(
                sub[time_col], sub["KM_Group"], sub[event_col]
            )
            p_value = float(mlr.p_value)

        # Plot
        plt.figure(figsize=(6,5), dpi=dpi)
        ax = None
        for g in group_order:
            part = sub.loc[sub["KM_Group"] == g]
            if part.empty:
                continue
            kmf = KaplanMeierFitter()
            label = f"{g} (n={len(part)})"
            kmf.fit(durations=part[time_col], event_observed=part[event_col], label=label)
            ax = kmf.plot_survival_function(ci_show=True, ax=ax)

        if KM_ymin_0:
            plt.ylim(0, 1)
        plt.xlabel("Days (Event-free survival)")
        plt.ylabel("Survival Probability")
        plt.title(f"{_gene} - {group_col}: {grp}\n{title_part}; log-rank p = {p_value:.4g}")
        plt.legend()
        plt.grid(True, alpha=0.3)

        # # Save (one file per group) if requested
        # if out_dir is not None:
        #     os.makedirs(out_dir, exist_ok=True)
        #     prefix = filename_prefix if filename_prefix else f"{_gene}"
        #     suffix = "median" if n_groups == 2 else (f"{n_groups}quantiles")
        #     safe_grp = "".join(c if c.isalnum() or c in "._- " else "_" for c in str(grp))
        #     fname = f"{prefix}_KM_{suffix}_{group_col.replace(' ','_')}_{safe_grp}.svg"
        #     out_path = os.path.join(out_dir, fname)
        #     try:
        #         WriteFile(out_path)  # your helper, if defined
        #     except NameError:
        #         plt.savefig(out_path, bbox_inches="tight")

        plt.show()

        # Store results
        dict_pvals[grp]  = p_value
        dict_tables[grp] = sub[[id_col, "Expression", "KM_Group"]].copy()

    return dict_pvals, dict_tables


def rank_median_shift_scatter(
    id_col="Patient_ID",
    clin_col="Subtype",
    clin_hit=None,                     # EXACT match only (str or list/set)
    gene_col="Gene",
    list_gene=None,
    case_insensitive=True,
    highlight_genes=None,
    eps=1e-6,
    show_plot=True,
    out_dir=None,
    stat="median",
    suffix=""
):
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    # ---------------------------
    # Expression preprocessing
    # ---------------------------
    expr = df_gexp.copy()

    if case_insensitive:
        expr[gene_col] = expr[gene_col].astype(str).str.upper()
        if list_gene is not None:
            list_gene = [str(g).upper() for g in list_gene]
        highlight_set = set([g.upper() for g in (highlight_genes or [])])
    else:
        expr[gene_col] = expr[gene_col].astype(str)
        highlight_set = set(highlight_genes or [])

    patient_cols = [c for c in expr.columns if c != gene_col]
    expr[patient_cols] = expr[patient_cols].apply(pd.to_numeric, errors="coerce")

    expr = (
        expr.groupby(gene_col, as_index=True)
            .median(numeric_only=True)
            .sort_index()
    )

    if list_gene is not None:
        keep = [g for g in list_gene if g in expr.index]
        expr = expr.loc[keep]

    # ---------------------------
    # Clinical PERFECT MATCH ONLY
    # ---------------------------
    if clin_hit is None:
        raise ValueError("clin_hit must be provided for exact matching.")

    if isinstance(clin_hit, (list, tuple, set)):
        mask = clin_df[clin_col].isin(clin_hit)
        hit_label = f"{clin_col} in {list(clin_hit)}"
    else:
        mask = clin_df[clin_col] == clin_hit
        hit_label = f"{clin_col} == {clin_hit}"

    hit_ids = (
        clin_df.loc[mask, id_col]
        .astype(str)
        .tolist()
    )

    hit_ids = [pid for pid in hit_ids if pid in expr.columns]
    if not hit_ids:
        raise ValueError("No overlap between selected clinical group and expression columns.")

    # ---------------------------
    # Statistic selection
    # ---------------------------
    label = "median"
    q = None

    if isinstance(stat, str):
        s = stat.lower()
        if s == "median":
            pass
        elif s in ("q1", "q25"):
            q = 0.25; label = "Q1"
        elif s in ("q3", "q75"):
            q = 0.75; label = "Q3"
        else:
            raise ValueError("stat must be 'median', 'q1', 'q3', or quantile in [0,1].")
    elif isinstance(stat, (float, int)):
        if 0 <= float(stat) <= 1:
            q = float(stat); label = f"Q{int(round(100*q))}"
        else:
            raise ValueError("Numeric stat must be a quantile in [0,1].")
    else:
        raise ValueError("Unsupported stat type.")

    if q is None:
        overall_stat = expr.median(axis=1, skipna=True)
        hit_stat     = expr[hit_ids].median(axis=1, skipna=True)
    else:
        overall_stat = expr.quantile(q, axis=1, interpolation="linear")
        hit_stat     = expr[hit_ids].quantile(q, axis=1, interpolation="linear")

    # ---------------------------
    # Relative shift
    # ---------------------------
    delta = hit_stat - overall_stat
    denom = np.maximum(np.abs(overall_stat), eps)
    relative_pct = 100.0 * delta / denom

    out = (
        pd.DataFrame({
            "Gene": expr.index,
            f"overall_{label}": overall_stat.values,
            f"hit_{label}": hit_stat.values,
            "delta": delta.values,
            "relative_pct": relative_pct.values,
        })
        .sort_values("relative_pct", ascending=False)
        .reset_index(drop=True)
    )
    out.insert(0, "Rank", out.index + 1)

    # ---------------------------
    # Save table
    # ---------------------------
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        fname = f"rank_shift_{label}_{suffix}.csv" if suffix else f"rank_shift_{label}.csv"
        out.to_csv(os.path.join(out_dir, fname), index=False)

    # ---------------------------
    # Plot
    # ---------------------------
    if show_plot:
        plt.figure(figsize=(6, 4.5), dpi=200)

        mask_hl = out["Gene"].isin(highlight_set)

        plt.scatter(
            out["Rank"], out["relative_pct"],
            s=12, alpha=0.8, color="lightgrey"
        )

        if mask_hl.any():
            plt.scatter(
                out.loc[mask_hl, "Rank"],
                out.loc[mask_hl, "relative_pct"],
                s=36, color="red", zorder=3
            )
            try:
                from adjustText import adjust_text
                texts = [
                    plt.text(r["Rank"], r["relative_pct"], r["Gene"], fontsize=9)
                    for _, r in out.loc[mask_hl].iterrows()
                ]
                adjust_text(texts, arrowprops=dict(arrowstyle="-", color="black"), force_explode=4)
            except Exception:
                for _, r in out.loc[mask_hl].iterrows():
                    plt.annotate(
                        r["Gene"],
                        (r["Rank"], r["relative_pct"]),
                        xytext=(0, 6),
                        textcoords="offset points",
                        fontsize=9
                    )

        plt.axhline(0, linestyle="--", linewidth=1)
        plt.xlabel(f"Rank (by relative {label} shift)")
        plt.ylabel(f"Relative {label} shift (%)")
        plt.title(f"{label} shift: {hit_label} {suffix}")
        plt.tight_layout()

        if out_dir:
            fig_name = f"rank_shift_{label}_{suffix}.svg" if suffix else f"rank_shift_{label}.svg"
            plt.savefig(os.path.join(out_dir, fig_name))

        plt.show()

    return out


# Identi
import pandas as pd
from statsmodels.stats.multitest import multipletests


#Find correlations of expressions of a gene above a certain threshold to clinical parameters
def ClinicalCorrelator(
    gene,
    cutoff=None,
    quantile=None,
    id_col="Patient_ID",
    max_categories=20,
    min_group_size=10
):
    """
    Scan clinical features enriched in gene-high population.

    Parameters
    ----------
    gene : str
        Gene symbol
    cutoff : float
        Expression cutoff to define "high"
    quantile : float (0-1)
        If cutoff not provided, use quantile threshold
    df_expr : dataframe
        Gene x sample expression matrix (Gene as index)
    clin_df : dataframe
        Clinical dataframe
    id_col : str
        Column linking samples
    """

    df_expr = df_gexp.set_index("Gene") if "Gene" in df_gexp.columns else df_gexp

    # Extract gene expression
    gene_df = (
        df_expr.loc[gene]
        .to_frame(name="expr")
        .reset_index()
        .rename(columns={"index": id_col})
    )

    df = clin_df.merge(gene_df, on=id_col, how="inner")

    # Define high group
    if cutoff is None and quantile is not None:
        cutoff = df["expr"].quantile(quantile)

    if cutoff is None:
        raise ValueError("Provide either cutoff or quantile")

    df["gene_high"] = df["expr"] > cutoff

    results_cat = []
    results_num = []

    # ---- CATEGORICAL VARIABLES ----
    for col in df.columns:

        if col in [id_col, "expr", "gene_high"]:
            continue

        if df[col].dtype == "object" or df[col].dtype.name == "category":

            if df[col].nunique() > max_categories:
                continue

            for category in df[col].dropna().unique():

                df["tmp"] = df[col] == category
                ct = pd.crosstab(df["gene_high"], df["tmp"])

                if ct.shape == (2, 2):
                    if ct.values.min() >= 1:
                        OR, p = fisher_exact(ct)
                        results_cat.append({
                            "type": "categorical",
                            "column": col,
                            "category": category,
                            "odds_ratio": OR,
                            "p_value": p
                        })

    # ---- NUMERIC VARIABLES ----
    for col in df.select_dtypes(include=["float","int"]).columns:

        if col in ["expr"]:
            continue

        high = df[df["gene_high"]][col].dropna()
        low  = df[~df["gene_high"]][col].dropna()

        if len(high) >= min_group_size and len(low) >= min_group_size:
            stat, p = mannwhitneyu(high, low)
            results_num.append({
                "type": "numeric",
                "column": col,
                "p_value": p,
                "median_high": high.median(),
                "median_low": low.median()
            })

    res_cat = pd.DataFrame(results_cat)
    res_num = pd.DataFrame(results_num)

    # Multiple testing correction
    if not res_cat.empty:
        res_cat["FDR"] = multipletests(res_cat["p_value"], method="fdr_bh")[1]
        res_cat = res_cat.sort_values("FDR")

    if not res_num.empty:
        res_num["FDR"] = multipletests(res_num["p_value"], method="fdr_bh")[1]
        res_num = res_num.sort_values("FDR")

    return res_cat



#%% ===========================================================================
# 4. Run this cell to scan for best correlations for one gene
# =============================================================================
gene            = 'MEN1'
min_target_expr = None #Usually None or 5 (threshold for null expression), change to omit expressions below this value
show_breakpoint = True
gene_set        = KTC_orthologue(KTC_GetGeneSet('WANG_MLL_TARGETS'), organism_origin='mmusculus', organism_target='hsapiens')
label           = 'WANG_MLL_TARGETS'
top_n           = 2
write_files     = True
category_mode   = 'elbows' #or "top_bottom"

# top_genes, bottom_genes = top_n_comparisons(
#     gene,
#     gene_set=gene_set,
#     label=label,
#     min_target_expr=min_target_expr,
#     top_n=top_n,
#     write_files=write_files,
#     out_dir=out_dir,
#     category_mode=category_mode
# )

sorted_genes, top_genes, bottom_genes = top_n_comparisons(
    gene=gene,
    gene_set=gene_set,
    label=label,
    min_target_expr=min_target_expr,
    top_n=top_n,
    write_files=write_files,
    out_dir=out_dir,
    category_mode='elbows'
)


for g, _ in top_genes:
    Grapher(gene, g, min_target_expr=min_target_expr)

for g, _ in bottom_genes:
    Grapher(gene, g, min_target_expr=min_target_expr)




#%% ===========================================================================
# 5. Run this cell to perform a single, specified comparison
# =============================================================================
#Overwrite 'target' and 'target2' and run this cell
#File is saved in out_dir/[target]
target           = 'MEN1' # The expression of the gene on the 1st axis
target2          = 'LMO2' # The expression of the gene on the 2nd axis
show_equation    = False
split_by_subtype = False # Instead of making one graph for all patients, make one expression graph for patients of each subtype
set_lim_0        = False
subanalysis_do   = False # Triggers the subanalysis: Make a new red line on the plot for a subset of the patients. Requires the next two folloding data.
subanalysis_col  = 'Subtype' # This column in the clinical data will be used to separate patients into two groups
subanalysis_hit  = 'ETP-like' # This value in the column above will be used to separate patients into two groups
pval_scientific  = False
top_n_residuals  = 0
min_target_expr  = 4.5

# Grapher(target, target2, split_by_subtype, subanalysis_do, subanalysis_col, subanalysis_hit,show_equation=False)
Grapher(target, target2, split_by_subtype,subanalysis_do=subanalysis_do, subanalysis_col=subanalysis_col, subanalysis_hit=subanalysis_hit,  show_equation=show_equation, set_lim_0=set_lim_0, pval_scientific=pval_scientific, top_n_residuals=top_n_residuals, min_target_expr=min_target_expr)

#%% 5b. Run all permutations for a set.
genes = ['PCNA', 'DHFR', 'NAMPT']
gene_combinations = set(itertools.combinations(genes, 2))
for target, target2 in gene_combinations:
    Grapher(target, target2, split_by_subtype, subanalysis_do=subanalysis_do, 
            subanalysis_col=subanalysis_col, subanalysis_hit=subanalysis_hit, 
            show_equation=show_equation, set_lim_0=set_lim_0)


#%% ===========================================================================
# 6. Analyze levels of expression across clinical parameters
# =============================================================================

clin_col   = 'Subtype' # Classifying Driver, ETP.STATUS, Sex, Race, CNS.Status, Insurance, Treatment.Arm, Subtype, Subsuptype, IP Status
gene       = 'CDKN2A' # The gene whose expression you want to track
palette    = 'Set2'  # The colors used in the graph. Eg. gray, Set1, Pastel1. Choose from: https://www.practicalpythonfordatascience.com/ap_seaborn_palette
dotcolor   = 'white' # The colors of the dots on top of the boxplots
fontsize   = 12 # The size of the text items
order      = None
set_ylim_0 = False # Force the 2nd axis to include 0
write_file = False # Write the graph to a file. Will be written to out_dir
list_n     = True # provide the number in each category
sort_median= True
median_col = 'black'
show_Q1Q3  = False
q_bar_color= 'black'

plot_type  = 'boxplot' # 'boxplot' or 'violinplot'

do_stats   = True # Perform a statistical analysis and include asterisks in the plot
stat_test  = 'Mann-Whitney' # 't-test_ind' or 'Mann-Whitney''
stat_text  = 'full' # 'star' for asterisks, 'full', for p-value

do_binary  = True
hit_binary = 'ETP-like'

mut_show   = False
mut_gene   = 'MYCN'
mut_aa     = 'All'
mut_col    = 'red'
mut_mark   = "."
mut_mark_s = 150


highlight_patients = []
highlight_label    = "Highlighted"
highlight_color    = "red"
highlight_marker   = "s"
highlight_size     = 35

# print(patients_loss_set)
# patients_loss_set = set(map(str, patients_loss))
# patients_mut_set  = set(map(str, patients_mut))


fig_size   = (6, 6) # width, height
dpi        = 200

SubsetBoxplotter(gene, clin_col, do_stats=do_stats, write_file=write_file, _palette=palette, _dotcolor=dotcolor, _fontsize=fontsize, set_ylim_0=set_ylim_0, list_n=list_n, sort_median=sort_median, do_binary=do_binary, hit_binary=hit_binary, order=order, mut_show=mut_show, mut_gene=mut_gene, mut_aa=mut_aa, mut_col=mut_col,mut_mark=mut_mark, mut_mark_s=mut_mark_s, stat_test=stat_test, stat_text=stat_text, fig_size=fig_size, dpi=dpi, plot_type=plot_type, show_Q1Q3=show_Q1Q3, q_bar_color=q_bar_color, median_col=median_col, highlight_patients=highlight_patients, highlight_label=highlight_label, highlight_color=highlight_color, highlight_marker=highlight_marker, highlight_size=highlight_size)



#%% =============================================================================
# 6b Polonen - Create a plot for all categories for a set of genes
# =============================================================================

for col in clin_df:
    print("'" + col + "',")

clin_cols = [
'Sex',
'Race',
'Ethnicity',
'Risk.Group',
'CNS.Status',
'Event.Type',
'Tumor.Specimen.Type',
'Germline.Specimen.Type',
'Relapse.Specimen.Type',
'ETP.STATUS',
'Standard Induction',
'Classifying Driver',
'ETP Status',
'Genetic Subtype',
'Subsubtype',
'Subtype',
'Pathway',
'IP Status',
'Driver.Gene',
'Reviewed.subtype',
'Reviewed.genetic.subtype',
'TCR.Rearrangement',
'somatic.TCR.SV',
'Classification',
'cluster.lowres.annotated',
'cluster.oncogene.annotated',
'monoallelic.expression.driver',
'ETP_status',
'Gender',
'Age.at.Diagnosis.in.Years_bin_5y',
'Blast.percentage_bin_10pct',
'Percent.Blasts.Tumor.Sample.Diagnostic_bin_10pct',
'Blast_percentage_TR_bin_10pct']




from matplotlib.backends.backend_pdf import PdfPages


gene = 'MEIS1'
out_pdf = os.path.join(out_dir, f"{gene}_clinical_covariates_violinplots.pdf")

with PdfPages(out_pdf) as pdf:
    for cc in clin_cols:
        SubsetBoxplotter(
            gene,
            cc,
            False,
            False,
            sort_median=True,
            plot_type='violinplot',
            plt_show = False,
            show_Q1Q3=True,
            q_bar_color='green',
            median_col='red'
        )
    
        # Save current figure to the PDF
        pdf.savefig(plt.gcf(), bbox_inches="tight")
        plt.close()  # important to avoid memory leaks
    print(f'pdf saved to: {out_pdf}')


#%% ===========================================================================
# 7. Polonen - Generate a Kaplan Meier plot of event-free survival for one gene
# =============================================================================
gene      = 'LMO2'
gene2     =  'MEN1' #String or None, Separates in four groups based on low/high
n_groups  = 2
KM_ymin_0 = True

split_value   = None
split_value2  = None
split_values  = None # For >2 genes
split_values2 = None # For >2 genes

KaplanMeier(gene, gene2=gene2, n_groups=n_groups, KM_ymin_0=KM_ymin_0, split_value=split_value, split_values=split_values, split_value2=split_value2, split_values2=split_values2)

#%% ===========================================================================
# 8. Polonen - Generate a Kaplan Meier plot for all clinical parameters
# =============================================================================

# Run for each clinical column
clin_cols = ['Classifying Driver', 'ETP.STATUS', 'Sex', 'Race', 'CNS.Status', 'Insurance', 
             'Treatment.Arm', 'Subtype', 'Subsubtype', 'IP Status', 'TCR.Rearrangement']

clin_cols = ['Classifying Driver']

for col in clin_cols:
    KaplanMeier_clinical(col)
    uniques = set(clin_df[col])
    uniques = list(uniques)
    uniques = [str(x) for x in uniques]
    print(col)
    print(', '.join(list(uniques)))


#%% ===========================================================================
# 9. Our own cell line MS data - Compare gene expression level trends
# =============================================================================

protein_x = 'MYCBP'
protein_y = 'ATP6AP2'
Grapher_MSpr1(protein1=protein_x, protein2=protein_y, df_msdataset=df_cell_line_MS)


#%% =============================================================================
# 10. CCLE density plot
# =============================================================================


gene          = 'KDM6B'
filter_col    = 'OncotreeSubtype' # 'DepmapModelType'
filter_val    = 'T-Cell Acute Lymphoblastic Leukemia' # 'Acute Myeloid Leukemia', T-Cell Acute Lymphoblastic Leukemia
label_col     = "StrippedCellLineName"   # nicer than ModelID usually
log2p1        = True
label_points  = True
figsize       = (8, 5)

Grapher_CCLE_density(gene=gene, filter_col=filter_col, filter_val=filter_val, log2p1=log2p1, label_points=label_points, figsize=figsize)


#%% =============================================================================
# 10b Printing the options for CCLE
# =============================================================================


print("\n[Filter Options]")
filter_columns = ['OncotreeLineage','OncotreePrimaryDisease', 'OncotreeSubtype', 'OncotreeCode']
for col in filter_columns:
    if col in df_cl.columns:
        print(f"\n\t{col}:")
        print(df_cl[col].value_counts(dropna=True).to_string())


#%% ===========================================================================
# 11. CCLE (cancer cell lines) data - Compare expression levels of two genes
# =============================================================================

Grapher_CCLE(
    gene_x        = 'MEIS1',
    gene_y        = 'GAPDH',
    show_equation = False,
    log_scale     = False,
    set_lim_0     = False,
    label_points  = True,
    filter_col    = 'OncotreeSubtype',
    filter_val    = 'B-Cell Acute Lymphoblastic Leukemia', # 'Acute Myeloid Leukemia'
    # filter_val    = 'acute_myeloid_leukaemia'
    )


#%% =============================================================================
# 12. CCLE (cancer cell lines) data - one gene across different categories
# =============================================================================

CCLE_Boxplotter(
    gene          = 'MSH6',
    group_by      = 'OncotreeCode',
    log_scale     = False,
    fig_height    = 5,
    fig_width     = 8,
    palette       = 'gray',
    do_stats      = False,
    min_samples   = 0,
    ylim          = 0,
    include_terms = ['EOV', 'HGSOC', 'MXOV', 'MOV', 'SOC', 'SCCO', 'CCOV'],
    add_scatter   = True
    )

#%% ===========================================================================
# 13. Plot the distribtution of expression levels in patients for one gene
# =============================================================================

gene='KDM6B'

Plot_Density(gene)

#%% ===========================================================================
# 14. Mutation_Barplotter (investigates mutations based on expression levels)
# =============================================================================

gene_mut   = 'MYCN' #The gene whose mutations you are interested in
mut_aa     = None #The mutation you are interested in (e.g. 'L72P')
gene_split = 'MYCN' # The gene used to split the population
cutoff     = 7.77 # THe expression level used to split the population (investigate with Plot_Density)
plot_abs   = False

Mutation_Barplotter(gene_mut, gene_split, cutoff, mut_aa, plot_abs)
Mutation_Barplotter(gene_mut, gene_split, cutoff, mut_aa, plot_abs=True)


#%% ===========================================================================
# 15. Ranked median shift 
# =============================================================================

clin_col        = 'Reviewed.subtype'
clin_hit        = 'ETP-like'
highlight_genes = ['ATP6AP2']
# highlight_genes = SPI1_targets

list_gene       = None
suffix          = ''

df_rms = rank_median_shift_scatter(clin_col=clin_col, clin_hit=clin_hit, highlight_genes=highlight_genes, suffix=suffix, list_gene=list_gene)

# df_rms.to_csv('/Users/kachrist/Desktop/out_dir/LMO2gd_like_rms.csv')

#%% =============================================================================
# 16. TARGET + TCGA expression levels across many cancers
# =============================================================================

gene            = 'ASNS'
filter_diseases = None
point_alpha     = 0.2
point_size      = 10
figsize         = (12,4)


TCGA_TARGET_expression(gene=gene, figsize=figsize, filter_diseases=filter_diseases, point_alpha=point_alpha, point_size=point_size)

#%% =============================================================================
# 16.b Report on  the data structure in TARGET:
# =============================================================================

for col in df_TCGA_expr:
    if not col.startswith("ENSG") and not col.startswith("PATIENT"):
        print(f"\n---# {col} #---")
        vc = df_TCGA_expr[col].value_counts(dropna=True)
        for val, count in vc.items():
            print(f"{val}\t{count}")

#%% =============================================================================
# 17. TCGA healthy v tumor
# =============================================================================

gene         = 'ASNS'
filter_col   = "_primary_site"
filter_value = "Ovary"
normals      = 'auto'
dpi          = 200
add_points   = True
add_pvalue   = True
point_color  = 'black'


_ = TCGA_healthy_v_tumor(df_TCGA_expr, gene=gene, filter_col=filter_col, filter_value=filter_value, normals=normals, add_points=add_points, point_color=point_color, add_pvalue=add_pvalue)

#%% ===========================================================================
# 18. Kaplan Meier plot with Polonen data using expression data AND clinical data
# =============================================================================


clin_cols = ['Classifying Driver', 'ETP.STATUS', 'Sex', 'Race', 'CNS.Status', 'Insurance', 
             'Treatment.Arm', 'Subtype', 'Subsubtype', 'IP Status']
# clin_cols = ['Classifying Driver']

gene = 'ASNS'

for col in clin_cols:
    pvals, tables = KaplanMeier_expression_clinical(
        _gene=gene,
        group_col=col,
        n_groups=2,
        out_dir="KM_plots_IGF2BP2_bySubtype"
    )

#%% =============================================================================
# 19. Correlating expressions of a gene above a certain threshold to clinical parameters
# =============================================================================

gene = 'SOX11'
cutoff = 9.55
quantile = None # Alternative to cutoff E.g. 0.9 = top 10%
max_categories = 20 #Skip a column with too many unique values
min_group_size = 10 #Fail if we have too few qualifying samples

df_clinScan = ClinicalCorrelator(gene=gene, cutoff=cutoff, quantile=quantile, max_categories=max_categories, min_group_size=min_group_size)
