#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%% ===========================================================================
# 1. Description
# =============================================================================
'''
The purpose of this script is to explore correlations between genes in different datasets.
There are several applications of this script:

1. Replace 'target' in section 4 and run that cell
   - This script will create a waterfall graph with genes ranked based on Pearson's R
   - Breakpoints will be identified using the Kneedle algorithm - All genes with their R- and p-values will be written to a csv together with their position relative to the breakpoints
3. Replace 'target' and 'target2' in section 5 and run that cell
   - Script will create a graph and calculate Pearson's R and associated p-value for the two specified genes
4. Replace 'gene' and 'clin_col' in cell 6 and run that cell
   - Script will create a series of boxplots for the expression levels for patients separated by unique values in a column in the clinical dataset
5. Replace 'gene' in section 7 and run that cell
   - Script will generate a Kaplan-Meier graph for event-free survival for that gene
6. Run section 8 to create a Kaplan-Meier graph for each clinical parameter
7. Replace 'protein_x' and 'protein_y' in section 9 and run that cell
   - Script will create a graph and calculate Pearson's R and associated p-value for the two specified proteins across triplicates in six cell lines
8. Replace 'gene1' and 'gene2' in section 10 and run that cell
   - Script will create a graph and calculate Pearson's R and associated p-value for the two specified genes using cancer cell line data from CCLE
9. Replace 'gene' and 'group_by' in section 11 and run that cell
   -  Script will create a series of boxplots for the expression levels for cancer cell lines in CCLE separated by unique values in a column in the annotation dataset
10. This section analyzes TARGET + TCGA expression levels for many cancers

Sidenote: Some of the functionalities of this script can rely on KTC_functions.py to define a set of genes. This script supports you simply defining them by yourself, but if you want this functionality (e.g. Querying MSigDB for gene sets) - you can find KTC_functions.py on this GitHub.
'''

#%% ===========================================================================
# 2. Setup and settings
# =============================================================================

#Modules
import os
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
from scipy import stats
from kneed import KneeLocator
import matplotlib.pyplot as plt
from itertools import combinations
import matplotlib.ticker as ticker
from collections import defaultdict
from scipy.signal import find_peaks
from scipy.stats import fisher_exact
from scipy.stats import gaussian_kde
from lifelines import KaplanMeierFitter
from KTC_functions import KTC_GetGeneSet
from lifelines.statistics import logrank_test
from scipy.stats import pearsonr, mannwhitneyu
from statannotations.Annotator import Annotator
from sklearn.linear_model import LinearRegression

files_directory = '/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists' #Directory where files for clinical and gene expression are stored
out_dir         = '/Users/kachrist/Desktop/out_dir' # Directory where files and images are written. Subdirectories for individual genes are created

#Initialization
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
df_wgs                   = pd.read_excel(os.path.join(files_directory, 'Polonen_Extended_Data.xlsx'), sheet_name='ST14_Alterations.SNVIndel')

print("Loading cell line MS data...")
df_cell_line_MS          = pd.read_excel(os.path.join(files_directory, 'MS_results_PRC-5607 2.xlsx'), sheet_name='S2 Quantified proteins')
print("Loading CCLE data...")
path_CCLE_rpkm = os.path.join(files_directory, 'CCLE_RNAseq_genes_rpkm_20180929.gct')
path_CCLE_cl   = os.path.join(files_directory, 'Cell_lines_annotations_20181226.txt')
df_CCLE_rpkm   = pd.read_csv(path_CCLE_rpkm, sep='\t', skiprows=2)
df_CCLE_cl     = pd.read_csv(path_CCLE_cl, sep='\t')

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
    df_M5_allesions_genes, df_M5_allesions_variants, df_M7_IP
]

# Convert patient identifiers into index for all dataframes
df_list = [df.set_index(df.columns[0]) for df in df_list]

# Merge all dataframes while handling overlapping columns
clin_df = df_list[0]
for i, df in enumerate(df_list[1:], start=1):
    clin_df = clin_df.join(df, how='outer', lsuffix='', rsuffix=f'_df{i}')

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


#Options
print_corr_genes = False # Genes above and below breakpoints can be written directly to console
write_files      = True # Turns on/off the writing of csv and pngs
log_scale        = False
show_breakpoint  = False # Identify and label breakpoints with Kneedle
#%% ===========================================================================
# 3. Main functions
# =============================================================================

def WriteFile(name):
    plt.savefig(os.path.join(out_dir, name))
    print('file created: %s' %(os.path.join(out_dir, name)))

#     print(name, target)
# 	out_dir_target = os.path.join(out_dir, target)
# 	if not os.path.exists(out_dir_target):
# 		os.makedirs(out_dir_target)
# 	if name.endswith('png') or name.endswith('svg'):
# 		plt.savefig(os.path.join(out_dir_target, name))
# 		print('file created: %s' %(os.path.join(out_dir_target, name)))

# This function is called elsewhere to create pairwise correlation plots between two genes
def Grapher(gene1, gene2, split_by_subtype=False, subanalysis_do=False, subanalysis_col=None, subanalysis_hit=None, show_equation=False, set_lim_0=False, pval_scientific=False, top_n_residuals=0):
    values1 = df_gexp.loc[df_gexp['Gene'] == gene1].iloc[0, 1:].tolist()
    values2 = df_gexp.loc[df_gexp['Gene'] == gene2].iloc[0, 1:].tolist()
    
    print(len(values1), len(values2))
    
    if log_scale:
        values1 = [value + 1 for value in values1]
        values2 = [value + 1 for value in values2]
    
    pecan_samples = df_gexp.columns[1:].tolist()
    sample_colors = ['black' for _ in pecan_samples]
    
    etp_indices = []
    if subanalysis_do:
        for i, sample in enumerate(pecan_samples):
            match = clin_df[clin_df['Patient_ID'] == sample]  # Changed to Patient_ID
            if not match.empty:
                stage = match[subanalysis_col].values[0]
                if stage == subanalysis_hit:
                    sample_colors[i] = 'red'
                    etp_indices.append(i)

    values1 = np.array(values1)
    values2 = np.array(values2)
    sample_colors = np.array(sample_colors)

    if split_by_subtype:
        match = clin_df[clin_df['Patient_ID'].isin(pecan_samples)]  # Changed to Patient_ID
        unique_subtypes = match['Classifying Driver'].dropna().unique()
        sample_subtypes = {row['Patient_ID']: row['Classifying Driver'] for _, row in match.iterrows()}  # Changed to Patient_ID
        
        for subtype in unique_subtypes:
            indices = [i for i, sample in enumerate(pecan_samples) if sample_subtypes.get(sample) == subtype]
            values1_sub = values1[indices]
            values2_sub = values2[indices]
            
            if len(values1_sub) > 1:
                fig, ax = plt.subplots(figsize=(8, 8), dpi=200)
                plt.scatter(values1_sub, values2_sub, alpha=0.5)
                
                model_sub = LinearRegression()
                model_sub.fit(values1_sub.reshape(-1, 1), values2_sub)
                Y_pred_sub = model_sub.predict(values1_sub.reshape(-1, 1))
                
                r_sub, p_sub = pearsonr(values1_sub, values2_sub)
                formatted_p_sub = '{:.2e}'.format(p_sub) if pval_scientific else f'{p_sub:.3f}'
                
                plt.plot(values1_sub, Y_pred_sub, label=f'R={r_sub:.2f}, p={formatted_p_sub}' + (f', y={model_sub.coef_[0]:.2f}x + {model_sub.intercept_:.2f}' if show_equation else ''), color='red')
                plt.xlabel(gene1 + ' Expression (VST)', fontsize=18)
                plt.ylabel(gene2 + ' Expression (VST)', fontsize=18)
                plt.tick_params(axis='both', labelsize=16)
                plt.title('PECAN expression: %s v %s\nSubset of patients with subtype: (%s)' % (gene1, gene2, subtype), fontsize=22)
                plt.legend(fontsize=16)
                file_name = 'PECAN_correlation_%s_v_%s_%s.svg' % (gene1, gene2, subtype.replace('/','_'))
                if set_lim_0 == True:
                    plt.ylim(0)
                    plt.xlim(0)
                WriteFile(file_name)
                plt.show()
                plt.close(fig)
    else:
        fig, ax = plt.subplots(figsize=(8, 8), dpi=200)
        plt.scatter(values1, values2, color=sample_colors, alpha=0.2)
        
        if subanalysis_do and len(etp_indices) > 1:
            values1_etp = values1[etp_indices]
            values2_etp = values2[etp_indices]
            
            r_etp, p_etp = pearsonr(values1_etp, values2_etp)
            formatted_p_etp = '{:.2e}'.format(p_etp) if pval_scientific else f'{p_etp:.3f}'

            
            model_etp = LinearRegression()
            model_etp.fit(values1_etp.reshape(-1, 1), values2_etp)
            Y_pred_etp = model_etp.predict(values1_etp.reshape(-1, 1))
            
            plt.plot(values1_etp, Y_pred_etp, label=f'{subanalysis_hit}: R={r_etp:.2f}, p={formatted_p_etp}' + (f', y={model_etp.coef_[0]:.2f}x + {model_etp.intercept_:.2f}' if show_equation else ''), color='red')
        
        if log_scale:
            log_values1 = np.log10(values1)
            log_values2 = np.log10(values2)
            r_value, p_value = pearsonr(log_values1, log_values2)
            Xs = log_values1.reshape(-1, 1)
            Ys = log_values2
        else:
            r_value, p_value = pearsonr(values1, values2)
            Xs = values1.reshape(-1, 1)
            Ys = values2
        
        formatted_p_value = '{:.2e}'.format(p_value) if pval_scientific else p_value

        
        model = LinearRegression()
        model.fit(Xs, Ys)
        Y_pred = model.predict(Xs)
        a, b = model.coef_[0], model.intercept_
        
        if top_n_residuals!=0:
            # --- Residual analysis ---
            residuals = np.abs(Ys - Y_pred)
            sorted_indices = np.argsort(residuals)
            highlight_indices = sorted_indices[:top_n_residuals]
            highlight_x = values1[highlight_indices]
            highlight_y = values2[highlight_indices]
            plt.scatter(highlight_x, highlight_y, color='blue', edgecolor='white', s=80, label='Best correlating samples')
            
            print(f"Top {top_n_residuals} samples closest to regression line:")
            print(f"Residual range: {residuals[sorted_indices[0]]:.4f} to {residuals[sorted_indices[top_n_residuals - 1]]:.4f}")
            closest_ids = [pecan_samples[i] for i in sorted_indices[:top_n_residuals]]
            print("Sample IDs:", ", ".join(closest_ids))
        
        plt.plot(values1, Y_pred, label='All: R=%.2f, p=%s%s' % (r_value, formatted_p_value, (', y=%.2fx + %.2f' % (a, b) if show_equation else '')), color='black')
        plt.xlabel(gene1 + ' Expression (VST)', fontsize=18)
        plt.ylabel(gene2 + ' Expression (VST)', fontsize=18)
        plt.tick_params(axis='both', labelsize=16)
        plt.title('PECAN expression: %s v %s' % (gene1, gene2), fontsize=22)
        plt.legend(fontsize=16)
        file_name = 'PECAN_correlation_%s_v_%s.svg' % (gene1, gene2)
        if set_lim_0 == True:
            plt.ylim(0)
            plt.xlim(0)
        WriteFile(file_name)
        plt.show()

# This plot is called to create the waterfall plot and return genes above and below breakpoints (for csv output)
def WaterfallPlot(dictionary, target_gene, gene_set, label):
   genes, r_p_values = zip(*dictionary)
   r_values = [r for r, p in r_p_values]
   ranks = np.arange(1, len(genes) + 1)
   fig, ax = plt.subplots(figsize=(8,5), dpi=200)
   plt.scatter(ranks, r_values, color='black', s=2, zorder=2,label='gene expression correlations')

   #
   highlight_genes = gene_set
   if len(highlight_genes) != 0:
      highlight_indices = [i for i, gene in enumerate(genes) if gene in highlight_genes]
      highlight_ranks = np.array(ranks)[highlight_indices]
      highlight_r_values = np.array(r_values)[highlight_indices]
      background_r_values = [r for r in r_values if r not in highlight_r_values]
      stat, p_value = mannwhitneyu(highlight_r_values, background_r_values, alternative='two-sided')
      for rank, r_value in zip(highlight_ranks, highlight_r_values):
         plt.vlines(x=rank, ymin=min(r_values), ymax=max(r_values), color='red', linewidth=0.2, zorder=1, label= '%s (%.0f, p=%f)' %(label, stat, p_value))

   plt.xlabel('Rank of gene to gene correlation', fontsize=16)
   plt.ylabel('Pearson\'s R', fontsize=16)
   plt.title('%s : Waterfall plot of gene correlations' %(target_gene), fontsize=22)
   kn_positive = KneeLocator(ranks, r_values, curve='convex', direction='decreasing')
   kn_negative = KneeLocator(ranks, r_values, curve='concave', direction='decreasing')

   genes_above_elbow = [genes[i] for i in range(kn_positive.knee)]
   genes_below_elbow = [genes[i] for i in range(kn_negative.knee, len(genes))]
#    ax.text(400, 0, ' genes above: %i' %(len(genes_above_elbow)), fontsize=12)
#    ax.text(14000, 0, 'genes below: %i ►' %(len(genes_below_elbow)), fontsize=12)
   if show_breakpoint:
      plt.axvline(x=kn_positive.knee, color='blue', linestyle='--', label='breakpoint')
      plt.axvline(x=kn_negative.knee, color='blue', linestyle='--')
   plt.legend(fontsize=12)

   # This section trims the legend down to only one per unique item
   handles, labels = plt.gca().get_legend_handles_labels()
   by_label = dict(zip(labels, handles))
   plt.legend(by_label.values(), by_label.keys(), fontsize=12)
   plt.show()

   file_name = 'PECAN_waterfall_%s.png' %(target_gene)
   WriteFile(file_name)

   if print_corr_genes:
      for gene in genes_above_elbow:
         print(gene)
      print()
      for gene in genes_below_elbow:
         print(gene)

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
def top_n_comparisons(gene, gene_set, label):
    # --- 1) Build a numeric matrix indexed by Gene once ---
    # expects df_gexp with column 'Gene' and sample columns after it
    df_num = df_gexp.copy()
    if 'Gene' not in df_num.columns:
        raise ValueError("df_gexp must contain a 'Gene' column.")
    df_num = df_num.drop_duplicates(subset='Gene').set_index('Gene')
    df_num = df_num.apply(pd.to_numeric, errors='coerce')  # coerce all to float

    if gene not in df_num.index:
        raise ValueError(f"Gene '{gene}' not found in df_gexp.")

    gene_values = df_num.loc[gene].to_numpy(dtype=float)

    r_dict = {}
    all_genes = df_num.index.tolist()

    # --- 2) Compare against every other gene safely ---
    for other_gene in tqdm(all_genes, desc='Comparing genes', unit='gene'):
        if other_gene == gene:
            continue
        other_gene_values = df_num.loc[other_gene].to_numpy(dtype=float)
        r_value, pearson_p = _safe_pearsonr(gene_values, other_gene_values)
        r_dict[other_gene] = (r_value, pearson_p)

    # --- 3) Sort by r (descending), putting NaNs at the end ---
    def sort_key(item):
        r = item[1][0]
        return (-1e9 if np.isnan(r) else r)

    sorted_genes = sorted(r_dict.items(), key=sort_key, reverse=True)

    # --- 4) WaterfallPlot + export (uses your existing globals) ---
    genes_above_elbow, genes_below_elbow = WaterfallPlot(sorted_genes, gene, gene_set, label)

    if write_files:
        out_dir_target = os.path.join(out_dir, target)
        os.makedirs(out_dir_target, exist_ok=True)
        rows = []
        for gene_name, (r_value, p_value) in sorted_genes:
            if gene_name in genes_above_elbow:
                category = 'above_1st_elbow'
            elif gene_name in genes_below_elbow:
                category = 'below_2nd_elbow'
            else:
                category = 'neither'
            rows.append([gene_name, r_value, p_value, category])
        df_result = pd.DataFrame(rows, columns=['Gene', 'Pearson_r', 'p_value', 'Category'])
        df_result.to_csv(os.path.join(out_dir_target, f'{gene}_correlation_data.csv'), index=False)
        print(f'file created: {out_dir_target}')

    # --- 5) Return top/bottom N (uses your global top_n) ---
    n = int(round(top_n/2))
    return (sorted_genes[:n], sorted_genes[-n:])


def SubsetBoxplotter(
    gene, PECAN_col, do_stats=True, write_file=False, _palette='gray',
    _dotcolor='white', _fontsize=14, order=None, set_ylim_0=False,
    list_n=False, sort_median=False, do_binary=False, hit_binary=None,
    mut_show=False, mut_gene=None, mut_aa=None,
    mut_mark='s', mut_col='red', mut_mark_s=30, stat_test='Mann-Whitney',
    stat_text='star', fig_size=(8,8), dpi=200, plot_type='boxplot'
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

        if sort_median:
            # Compute both median and mean
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
            mut_df = df_wgs[(df_wgs['gene'] == mut_gene) & (df_wgs['aa_change'] == mut_aa)]
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
                color="black", linewidth=1.5, zorder=4
            )
    else:
        print('plot_type must be either "boxplot" or "violinplot"')

    # === Always overlay stripplot ===
    sns.stripplot(
        data=data[~data['Patient_ID'].isin(mutated_patients)],
        x='Subtype_Labeled', y='Expression',
        facecolor=_dotcolor, edgecolor='black',
        linewidth=0.8, alpha=0.3, jitter=True,
        order=label_order, zorder=2
    )

    if mut_show:
        mut_data = data[data['Patient_ID'].isin(mutated_patients)]
        if not mut_data.empty:
            sns.scatterplot(data=mut_data, x='Subtype_Labeled', y='Expression',
                            marker=mut_mark, s=mut_mark_s, color=mut_col,
                            edgecolor='black', linewidth=0.8, zorder=10, label=f'{mut_gene} - {mut_aa}')
        plt.legend()

    if do_stats and data['Subtype'].nunique() >= 2:
        # 1) Use the same order as the plot
        pairs = [(label_order[i], label_order[j])
                 for i in range(len(label_order)) for j in range(i+1, len(label_order))]

        annotator = Annotator(
            ax, pairs, data=data,
            x='Subtype_Labeled',
            y='Expression',
            order=label_order
        )

        # 2) Put annotations above the violins (cleanest) and give offsets
        cfg = dict(test=stat_test, text_format=stat_text, loc='outside',
                   line_height=0.02, line_offset=0.02, text_offset=0.01,
                   verbose=0)
        annotator.configure(**cfg)
        annotator.hide_non_significant = True
        annotator.apply_and_annotate()

        # 3) Ensure there's headroom for the stars
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax + 0.1 * (ymax - ymin))

    plt.title(f'PeCan expression of {gene} grouped by: \"{PECAN_col}\"', fontsize=_fontsize)
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

def KaplanMeier(
    _gene,
    n_groups=2,                     # 2 = median split; 4 = quartiles; any int >=2 works
    id_col="Patient_ID",
    time_col="EFS",
    event_col="EFS.status",
    dpi=200,                        # your default preference
    out_dir=None,                   # if given, will save a figure
    filename=None                   # custom filename if desired
):
    # 1) Align clinical and expression data
    matched = clin_df[clin_df[id_col].isin(df_gexp.columns)].copy()
    expr_series = df_gexp.set_index("Gene").loc[_gene, matched[id_col]].astype(float)
    matched["Expression"] = matched[id_col].map(expr_series)

    # Drop rows with missing essentials
    matched = matched.dropna(subset=["Expression", time_col, event_col]).copy()
    matched[event_col] = matched[event_col].astype(int)

    # 2) Build groups
    if n_groups == 2:
        med = matched["Expression"].median()
        matched["KM_Group"] = np.where(matched["Expression"] > med, "High", "Low")
        group_order = ["Low", "High"]
        title_part = "High vs Low (median split)"
    else:
        # Try qcut directly; fall back to ranking if bin edges collide
        labels = [f"Q{i+1}" for i in range(n_groups)]
        try:
            matched["KM_Group"] = pd.qcut(matched["Expression"], q=n_groups, labels=labels)
            # Ensure we got the requested number of groups; if not, fallback to rank-based
            if matched["KM_Group"].nunique() < n_groups:
                raise ValueError("Collapsed bins in qcut; using rank-based qcut.")
        except Exception:
            ranked = matched["Expression"].rank(method="first")
            matched["KM_Group"] = pd.qcut(ranked, q=n_groups, labels=labels)

        group_order = labels  # Q1(lowest) ... Qn(highest)
        title_part = f"{n_groups} quantiles"

    # 3) Compute p-value (pairwise or multivariate log-rank)
    if n_groups == 2:
        g_low  = matched.query("KM_Group == 'Low'")
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

    # 4) Plot
    plt.figure(figsize=(6,5), dpi=dpi)
    ax = None
    for grp in group_order:
        sub = matched.loc[matched["KM_Group"] == grp]
        if sub.empty:
            continue
        kmf = KaplanMeierFitter()
        label = f"{grp} (n={len(sub)})"
        kmf.fit(durations=sub[time_col], event_observed=sub[event_col], label=label)
        ax = kmf.plot_survival_function(ci_show=True, ax=ax)

    plt.ylim(0, 1)
    plt.xlabel("Days (Event-free survival)")
    plt.ylabel("Survival Probability")
    plt.title(f"{_gene} — Kaplan–Meier by expression: {title_part}\nlog-rank p = {p_value:.4g}")
    plt.legend()
    plt.grid(True, alpha=0.3)

    # 5) Save if requested (compatible with your WriteFile helper if present)
    if out_dir is not None:
        if filename is None:
            suffix = "quartiles" if n_groups == 4 else (f"{n_groups}quantiles" if n_groups != 2 else "median")
            filename = f"{_gene}_KaplanMeier_{suffix}.svg"
        out_path = os.path.join(out_dir, filename)
        try:
            WriteFile(out_path)  # your helper, if defined
        except NameError:
            plt.savefig(out_path, bbox_inches="tight")

    plt.show()
    return p_value, matched[[id_col, "Expression", "KM_Group"]].copy()


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

def Grapher_CCLR(gene_x, gene_y, show_equation=False, log_scale=False, set_lim_0=False,
                 filter_col=None, filter_val=None, label_points=False):
    # --- Map gene symbols to ENSG IDs --- 
    gene_map = df_CCLE_rpkm[['Name', 'Description']].set_index('Description')['Name'].to_dict()
    if gene_x not in gene_map or gene_y not in gene_map:
        raise ValueError(f"Gene {gene_x} or {gene_y} not found in expression data.")
    
    ensg_x = gene_map[gene_x]
    ensg_y = gene_map[gene_y]

    # --- Filter cell lines based on annotation --- 
    if filter_col and filter_val:
        filtered_cl = df_CCLE_cl[df_CCLE_cl[filter_col] == filter_val]
    else:
        filtered_cl = df_CCLE_cl.copy()
    cell_lines = filtered_cl['CCLE_ID'].values

    # --- Extract and transpose expression matrix --- 
    df_expr = df_CCLE_rpkm.set_index('Name').drop(columns='Description')
    df_expr_T = df_expr.T
    df_expr_T.index.name = 'CCLE_ID'

    # --- Subset for the genes and cell lines of interest --- 
    available_lines = df_expr_T.index.intersection(cell_lines)
    df_subset = df_expr_T.loc[available_lines, [ensg_x, ensg_y]].dropna()
    df_subset.columns = [gene_x, gene_y]

    # --- Log transform if specified --- 
    if log_scale:
        with np.errstate(divide='ignore', invalid='ignore'):
            df_subset[gene_x] = np.log10(df_subset[gene_x] + 1)
            df_subset[gene_y] = np.log10(df_subset[gene_y] + 1)
        ylabel = 'log10(RPKM + 1)'
    else:
        ylabel = 'RPKM'

    # --- Perform regression --- 
    x = df_subset[gene_x].values
    y = df_subset[gene_y].values
    if len(x) > 1:
        r_value, p_value = pearsonr(x, y)
        model = LinearRegression()
        model.fit(x.reshape(-1, 1), y)
        x_range = np.linspace(min(x), max(x), 100)
        y_pred = model.predict(x_range.reshape(-1, 1))
        reg_label = f'R={r_value:.2f}, p={p_value:.2e}'
        if show_equation:
            reg_label += f'\n y = {model.coef_[0]:.2f} x + {model.intercept_:.2f}'
    else:
        r_value, p_value, reg_label = (np.nan, np.nan, "Insufficient data")

    # --- Plotting --- 
    plt.figure(figsize=(8, 8), dpi=150)
    sns.scatterplot(x=gene_x, y=gene_y, data=df_subset, alpha=0.2, edgecolor=None, color='black')

    if len(x) > 1:
        plt.plot(x_range, y_pred, color='black', label=reg_label)

    title = f'{gene_x} vs {gene_y} - CCLE'
    if filter_col and filter_val:
        title += f'\n({filter_col} = {filter_val})'

    plt.xlabel(f'{gene_x} Expression ({ylabel})', fontsize=18)
    plt.ylabel(f'{gene_y} Expression ({ylabel})', fontsize=18)
    plt.title(title, fontsize=22)
    plt.legend(fontsize=14)
    plt.tick_params(axis='both', labelsize=16)

    if set_lim_0:
        plt.xlim(left=0)
        plt.ylim(bottom=0)

    # --- Label points if specified --- 
    if label_points:
            from adjustText import adjust_text
            texts = []
    
            # Create a mapping from CCLE_ID to Name
            id_to_name = df_CCLE_cl.set_index('CCLE_ID')['Name'].to_dict()
    
            for i, row in df_subset.iterrows():
                label = id_to_name.get(i, i)  # fallback to IDs if Name not found
                texts.append(plt.text(row[gene_x], row[gene_y], label, fontsize=12))
    
            adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey'), forcePoints=0.5)

    plt.tight_layout()
    plt.show()

def CCLE_Boxplotter(
    gene,
    group_by='Histology',
    log_scale=False,
    fig_height=8,
    fig_width=8,
    palette='gray',
    do_stats=False,
    min_samples=5,
    ylim=0,
    include_terms=None,       # list of specific group terms to plot
    highlight_labels=None     # NEW: list of group names to make red
):
    # --- Step 1: Get expression data for the gene ---
    gene_expr = df_CCLE_rpkm.loc[df_CCLE_rpkm['Description'] == gene] \
                            .drop(['Name', 'Description'], axis=1) \
                            .squeeze()  # Series: index = CCLE_IDs

    # --- Step 2: Construct dataframe with metadata ---
    plot_data = pd.DataFrame({
        'CCLE_ID': gene_expr.index,
        'Expression': gene_expr.values
    })

    # Merge with annotation dataframe
    plot_data = plot_data.merge(df_CCLE_cl[['CCLE_ID', group_by]], on='CCLE_ID', how='left')
    plot_data = plot_data.dropna(subset=[group_by])

    # --- Step 3: Optional filter to only specified terms ---
    if include_terms is not None:
        plot_data = plot_data[plot_data[group_by].isin(include_terms)]

    # --- Step 4: Filter out groups with fewer than min_samples ---
    group_counts = plot_data[group_by].value_counts()
    valid_groups = group_counts[group_counts >= min_samples].index
    plot_data = plot_data[plot_data[group_by].isin(valid_groups)]

    # --- Step 5: Log-transform if needed ---
    if log_scale:
        plot_data['Expression'] = np.log10(plot_data['Expression'] + 1)
        ylabel = 'log10(RPKM + 1)'
    else:
        ylabel = 'RPKM'

    # --- Step 6: Order groups by median expression ---
    group_order = plot_data.groupby(group_by)['Expression'] \
                           .median() \
                           .sort_values(ascending=False) \
                           .index.tolist()

    # --- Step 7: Plot ---
    plt.figure(figsize=(fig_width, fig_height), dpi=300)
    ax = sns.boxplot(
        x=group_by,
        y='Expression',
        data=plot_data,
        order=group_order,
        palette=palette
    )

    # --- Step 8: Optional statistical annotations ---
    if do_stats and plot_data[group_by].nunique() >= 2:
        pairs = list(combinations(group_order, 2))
        annotator = Annotator(ax, pairs, data=plot_data, x=group_by, y='Expression')
        annotator.configure(test='t-test_ind', text_format='star', loc='inside', verbose=0)
        annotator.hide_non_significant = True
        annotator.apply_and_annotate()

    # --- Step 9: Highlight tick labels if requested ---
    if highlight_labels:
        for label in ax.get_xticklabels():
            if label.get_text() in highlight_labels:
                label.set_color('red')
                label.set_fontweight('bold')

    # --- Step 10: Final plot settings ---
    plt.title(f'{gene} Expression by {group_by}', fontsize=22)
    plt.xlabel(group_by, fontsize=18)
    plt.ylabel(ylabel, fontsize=18)
    plt.xticks(rotation=90)
    plt.tick_params(axis='both', labelsize=16)
    plt.ylim(bottom=ylim)
    plt.tight_layout()
    plt.show()

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
        mut_df = df_wgs[(df_wgs['gene'] == gene_mut) & (df_wgs['aa_change'] != '.')]
        mutation_label = f'mut {gene_mut}'
    else:
        mutation_string = f'p.{mut_aa}'
        mut_df = df_wgs[(df_wgs['gene'] == gene_mut) & (df_wgs['aa_change'] == mutation_string)]
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

def Plot_Density(gene):
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
    WriteFile(os.path.join(out_dir, f'{gene}_Density.svg'))
    plt.show()


# === Plotting function ===
def TCGA_TARGET_expression(gene):
    if gene not in df_target_tcga_merged.columns:
        print(f"{gene} not found.")
        return

    subset = df_target_tcga_merged[[gene, "_study", "_primary_disease"]].dropna()

    # Duplicate T-ALL rows and relabel _primary_disease and _study
    # tall_labeled = subset[subset["TALL_status"]].copy()
    # tall_labeled["_primary_disease"] = "T-ALL (only)"
    # tall_labeled["_study"] = "T-ALL-only"  # New hue category
    # print(tall_labeled)

    # Append to existing data
    # subset_combined = pd.concat([subset, tall_labeled])

    # Create x-axis order (including new label)
    ordered = subset.groupby("_primary_disease")[gene].median().sort_values(ascending=False).index

    # Plot
    plt.figure(figsize=(14, 6), dpi=200)
    sns.boxplot(data=subset, x="_primary_disease", y=gene, hue="_study", order=ordered, fliersize=0, linewidth=1)

    plt.xticks(rotation=90)
    plt.title(f"{gene} expression across TCGA and TARGET tumors)")
    plt.ylabel("log2(TPM + 0.001)" if subset[gene].max() < 100 else "TPM")
    plt.xlabel("Cancer type")
    plt.legend(title="Cohort", loc='upper right')
    plt.tight_layout()
    plt.show()


#%% ===========================================================================
# 4. Run this cell to scan for best correlations for one gene
# =============================================================================

# # Gene sets whose genes will be highlighted on the Waterfall Plot
gene_set = KTC_GetGeneSet('splicing_factors')
label = ', '.join(KTC_GetGeneSet('splicing_factors'))

# Every gene in the list below ('targets') will have all gene-gene correlations to it calculated and plotted.
targets = ['PCNA']

top_n            = 10 # E.g. 10 will generate graphs for the 5 most positive and most negatively correlated genes
for target in targets:
    bottom_genes, top_genes = top_n_comparisons(target, gene_set, label)
    if top_n > 0:
        for bottom_gene, _ in bottom_genes:
            Grapher(target, bottom_gene)
        for top_gene, _ in top_genes:
            Grapher(target, top_gene)


#%% ===========================================================================
# 5. Run this cell to perform a single, specified comparison
# =============================================================================
#Overwrite 'target' and 'target2' and run this cell
#File is saved in out_dir/[target]
target  = 'MYC' # The expression of the gene on the 1st axis
target2 = 'TONSL' # The expression of the gene on the 2nd axis
show_equation    = False
split_by_subtype = False # Instead of making one graph for all patients, make one expression graph for patients of each subtype
set_lim_0        = False
subanalysis_do   = False # Triggers the subanalysis: Make a new red line on the plot for a subset of the patients. Requires the next two folloding data.
subanalysis_col  = 'Subtype' # This column in the clinical data will be used to separate patients into two groups
subanalysis_hit  = 'SPI1' # This value in the column above will be used to separate patients into two groups
pval_scientific  = True
top_n_residuals  = 0

# Grapher(target, target2, split_by_subtype, subanalysis_do, subanalysis_col, subanalysis_hit,show_equation=False)
Grapher(target, target2, split_by_subtype,subanalysis_do=subanalysis_do, subanalysis_col=subanalysis_col, subanalysis_hit=subanalysis_hit,  show_equation=show_equation, set_lim_0=set_lim_0, pval_scientific=pval_scientific, top_n_residuals=top_n_residuals)

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

clin_col   = 'ETP.STATUS' #Classifying Driver, ETP.STATUS, Sex, Race, CNS.Status, Insurance, Treatment.Arm, Subtype, Subsuptype, IP Status
gene       = 'TYMS' # The gene whose expression you want to track
palette    = 'gray'  # The colors used in the graph. Choose from: https://www.practicalpythonfordatascience.com/ap_seaborn_palette
dotcolor   = 'white' # The colors of the dots on top of the boxplots
fontsize   = 12 # The size of the text items
#order      = ['ETP', 'Near-ETP', 'Non-ETP'] # Specify the order. Set to None or make sure the items are represented in the clin_col
order      = None
set_ylim_0 = False # Force the 2nd axis to include 0
write_file = True # Write the graph to a file. Will be written to out_dir
list_n     = True # provide the number in each category
sort_median= True

plot_type  = 'boxplot' # 'boxplot' or 'violinplot'

do_stats   = True # Perform a statistical analysis and include asterisks in the plot
stat_test  = 'Mann-Whitney' # 't-test_ind' or 'Mann-Whitney''
stat_text  = 'star' # 'star' for asterisks, 'full', for p-value

do_binary  = True
hit_binary = 'ETP'

mut_show   = False
mut_gene   = 'MYCN'
mut_aa     = 'p.P44L'
mut_col    = 'red'
mut_mark   = "."
mut_mark_s = 150

fig_size   = (4.5, 6) # width, height
dpi        = 200

SubsetBoxplotter(gene, clin_col, do_stats=do_stats, write_file=write_file, _palette=palette, _dotcolor=dotcolor, _fontsize=fontsize, set_ylim_0=set_ylim_0, list_n=list_n, sort_median=sort_median, do_binary=do_binary, hit_binary=hit_binary, order=order, mut_show=mut_show, mut_gene=mut_gene, mut_aa=mut_aa, mut_col=mut_col,mut_mark=mut_mark, mut_mark_s=mut_mark_s, stat_test=stat_test, stat_text=stat_text, fig_size=fig_size, dpi=dpi, plot_type=plot_type)



#%% =============================================================================
# 6b Polonen - Create a plot for all categories for a set of genes
# =============================================================================

clin_cols = ['Classifying Driver', 'ETP.STATUS', 'Sex', 'Race', 'CNS.Status', 'Insurance', 'Treatment.Arm', 'Subtype', 'Subsuptype', 'IP Status']
# clin_cols = ['ETP.STATUS', ]
genes     = ['KDM6B']
for cc in clin_cols:
    for gene in genes:
        # SubsetBoxplotter(gene, cc,True,True)
        SubsetBoxplotter(gene, cc, False, False)


#%% ===========================================================================
# 7. Polonen - Generate a Kaplan Meier plot of event-free survival for one gene
# =============================================================================

gene = 'TONSL'
KaplanMeier(gene, n_groups=2)

#%% ===========================================================================
# 8. Polonen - Generate a Kaplan Meier plot for all clinical parameters
# =============================================================================

# Run for each clinical column
clin_cols = ['Classifying Driver', 'ETP.STATUS', 'Sex', 'Race', 'CNS.Status', 'Insurance', 
             'Treatment.Arm', 'Subtype', 'Subsubtype', 'IP Status']

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

protein_x = 'NOTCH1'
protein_y = 'NAMPT'
Grapher_MSpr1(protein1=protein_x, protein2=protein_y, df_msdataset=df_cell_line_MS)


#%% ===========================================================================
# 10. CCLE (cancer cell lines) data - Compare expression levels of two genes
# =============================================================================

Grapher_CCLR(
    gene_x        = 'NOTCH1',
    gene_y        = 'NAMPT',
    show_equation = False,
    log_scale     = False,
    set_lim_0     = False,
    label_points  = True,
    filter_col    = 'Hist_Subtype1',
    filter_val    = 'acute_lymphoblastic_T_cell_leukaemia'
    )

print("\n[Filter Options]")
filter_columns = ['Site_Primary','Histology', 'Hist_Subtype1', 'Hist_Subtype2', 'Gender', 'Life_Stage', 'Race']
for col in filter_columns:
    if col in df_CCLE_cl.columns:
        vals = sorted(df_CCLE_cl[col].dropna().unique())
        print(f"\n\t{col}:\n{vals}")


#%% =============================================================================
# 11. CCLE (cancer cell lines) data - one gene across different categories
# =============================================================================

for gene in KTC_GetGeneSet('Laura'):
    
    try:
        CCLE_Boxplotter(
            gene          = gene,
            group_by      = 'Hist_Subtype1',
            log_scale     = False,
            fig_height    = 5,
            fig_width     = 6,
            palette       = 'gray',
            do_stats      = False,
            min_samples   = 5,
            ylim          = 0,
            include_terms = ['acute_lymphoblastic_B_cell_leukaemia', 'acute_lymphoblastic_T_cell_leukaemia', 'acute_myeloid_leukaemia'],
            )
        
        CCLE_Boxplotter(
            gene          = gene,
            group_by      = 'Hist_Subtype1',
            log_scale     = False,
            fig_height    = 5,
            fig_width     = 15,
            palette       = 'gray',
            do_stats      = False,
            min_samples   = 5,
            ylim          = 0,
            highlight_labels=['acute_lymphoblastic_B_cell_leukaemia', 'acute_lymphoblastic_T_cell_leukaemia', 'acute_myeloid_leukaemia']
            )
    except:
        print(f'failed for {gene}')


#%% ===========================================================================
# 12. Plot the distribtution of expression levels in patients for one gene
# =============================================================================

Plot_Density('TONSL')

#%% ===========================================================================
# 13. Mutation_Barplotter (investigates mutations based on expression levels)
# =============================================================================

gene_mut   = 'NOTCH1' #The gene whose mutations you are interested in
mut_aa     = None #The mutation you are interested in (e.g. 'L72P')
gene_split = 'NOTCH1' # The gene used to split the population
cutoff     = 13.36 # THe expression level used to split the population (investigate with Plot_Density)
plot_abs   = False

Mutation_Barplotter(gene_mut, gene_split, cutoff, mut_aa, plot_abs)
Mutation_Barplotter(gene_mut, gene_split, cutoff, mut_aa, plot_abs=True)

#%% =============================================================================
# Plotting gene expression levels across many cancers (TARGET + TCGA)
# =============================================================================

TCGA_TARGET_expression("TONSL")
