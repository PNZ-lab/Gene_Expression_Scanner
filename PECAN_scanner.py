#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%% ===========================================================================
# 1. Description
# =============================================================================
'''
The purpose of this script is to explore correlations between genes in the PECAN dataset.
There are four applications of this script:
	1. Replace 'target' in section 4 and run that cell
		- This script will create a waterfall graph with genes ranked based on Pearson's R
		- Breakpoints will be identified using the Kneedle algorithm
		- All genes with their R- and p-values will be written to a csv together with their position relative to the breakpoints
	2. Replace 'target' and 'target2' in section 5 and run that cell
		- Script will create a graph and calculate Pearson's R and associated p-value for the two specified genes
	3. Replace 'gene' and 'clin_col' in cell 6 and run that cell
		Script will create a series of boxplots for the expression levels for patients separated by unique values in a column in the clinical dataset
	4. Replace 'gene' in section 7 and run that cell.
		Script will generate a Kaplan-Meier graph for event-free survival for that gene
'''

#%% ===========================================================================
# 2. Setup and settings
# =============================================================================

#Modules
import os
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
from scipy import stats
from kneed import KneeLocator
import matplotlib.pyplot as plt
from itertools import combinations
from collections import defaultdict
from KTC_functions import KTC_GetGeneSet
from scipy.stats import pearsonr, mannwhitneyu
from statannotations.Annotator import Annotator
from sklearn.linear_model import LinearRegression

files_directory = '/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists' #Directory where files for clinical and gene expression are stored
out_dir         = '/Users/kasperthorhaugechristensen/Desktop/Dumpbox/Cristina' # Directory where files and images are written. Subdirectories for individual genes are created

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
print("Loading cell line MS data")
df_cell_line_MS          = pd.read_excel(os.path.join(files_directory, 'MS_results_PRC-5607 2.xlsx'), sheet_name='S2 Quantified proteins')

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


#Options
top_n            = 10 # E.g. 10 will generate graphs for the 5 most positive and most negatively correlated genes
print_corr_genes = False # Genes above and below breakpoints can be written directly to console
write_files      = True # Turns on/off the writing of csv and pngs
log_scale        = False
show_breakpoint  = False # Identify and label breakpoints with Kneedle
#%% ===========================================================================
# 3. Main functions: WriteFile: (), Grapher (expression correlation of two proteins),
# =============================================================================

def WriteFile(name):
	print(name, target)
	out_dir_target = os.path.join(out_dir, target)
	if not os.path.exists(out_dir_target):
		os.makedirs(out_dir_target)
	if name.endswith('png') or name.endswith('svg'):
		plt.savefig(os.path.join(out_dir_target, name))
		print('file created: %s' %(os.path.join(out_dir_target, name)))

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

   file_name = 'PECAN_waterfall_%s.png' %(target_gene)
   WriteFile(file_name)

   if print_corr_genes:
      for gene in genes_above_elbow:
         print(gene)
      print()
      for gene in genes_below_elbow:
         print(gene)

   return genes_above_elbow, genes_below_elbow

#This function calculates R and p values for all genes to one gene, returns a ranked list for the Waterfall function, and writes the csv
def top_n_comparisons(gene, gene_set, label):
   gene_values = df_gexp.loc[df_gexp['Gene'] == gene].iloc[0, 1:].values
   r_dict = {}
   all_genes = df_gexp['Gene']

   for other_gene in tqdm(all_genes, desc='Comparing genes', unit='gene'):
      if other_gene != gene:
         other_gene_values = df_gexp.loc[df_gexp['Gene'] == other_gene].iloc[0, 1:].values
         r_value, pearson_p = pearsonr(gene_values, other_gene_values)
         r_dict[other_gene] = (r_value, pearson_p)

   sorted_genes = sorted(r_dict.items(), key=lambda item: item[1][0], reverse=True)
   genes_above_elbow, genes_below_elbow = WaterfallPlot(sorted_genes, gene, gene_set, label)

   if write_files:
      out_dir_target = os.path.join(out_dir, target)
      if not os.path.exists(out_dir_target):
         os.makedirs(out_dir_target)
      data = []
      for gene_name, (r_value, p_value) in sorted_genes:
         if gene_name in genes_above_elbow:
            category = 'above_1st_elbow'
         elif gene_name in genes_below_elbow:
            category = 'below_2nd_elbow'
         else:
            category = 'neither'
         data.append([gene_name, r_value, p_value, category])
      df_result = pd.DataFrame(data, columns=['Gene', 'Pearson_r', 'p_value', 'Category'])
      df_result.to_csv(os.path.join(out_dir_target, '%s_correlation_data.csv' %(gene)), index=False)
      print('file created: %s' %(out_dir_target))

   return(sorted_genes[:round(top_n/2)], sorted_genes[-round(top_n/2):])

#%% ===========================================================================
# 4. Run this cell to scan for best correlations for one gene
# =============================================================================

# # Gene sets whose genes will be highlighted on the Waterfall Plot
#Use KTC_GetGene
# name_gene_set = 'HALLMARK_MYC_TARGETS_V1'
# name_gene_set = 'HALLMARK_ESTROGEN_RESPONSE_EARLY'
label = ', '.join(KTC_GetGeneSet('Freya'))
# label = 'm6a_readers'
#KTC_GetGeneSet can take: a list of gene names, a single gene name, or the name of a gene set from Msigdb
gene_set = KTC_GetGeneSet('Freya')

#Targets is a list of gene names who will have correlations for all other genes calculated (top hits will be )
targets = ['NAMPT']

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
#Overwrite 'target' and 'target2' abd run this cell
#File is saved in out_dir/[target]
target  = 'MYC' # The expression of the gene on the 1st axis
target2 = 'IGF2BP2' # The expression of the gene on the 2nd axis
show_equation    = False
split_by_subtype = True # Instead of making one graph for all patients, make one expression graph for patients of each subtype
set_lim_0        = False
subanalysis_do   = False # Triggers the subanalysis: Make a new red line on the plot for a subset of the patients. Requires the next two folloding data.
subanalysis_col  = 'CNS.Status' # This column in the clinical data will be used to separate patients into two groups
subanalysis_hit  = 'CNS 3c' # This value in the column above will be used to separate patients into two groups
pval_scientific  = True
top_n_residuals  = 0

# Grapher(target, target2, split_by_subtype, subanalysis_do, subanalysis_col, subanalysis_hit,show_equation=False)
Grapher(target, target2, split_by_subtype,subanalysis_do=subanalysis_do, subanalysis_col=subanalysis_col, subanalysis_hit=subanalysis_hit,  show_equation=show_equation, set_lim_0=set_lim_0, pval_scientific=pval_scientific, top_n_residuals=top_n_residuals)

#%% 5b. Run all permutations for a set.
import itertools
genes = ['MTOR', 'EIF']
gene_combinations = set(itertools.combinations(genes, 2))
for target, target2 in gene_combinations:
    Grapher(target, target2, split_by_subtype, subanalysis_do=subanalysis_do, 
            subanalysis_col=subanalysis_col, subanalysis_hit=subanalysis_hit, 
            show_equation=show_equation, set_lim_0=set_lim_0)


#%% ===========================================================================
# 6. Analyze levels of expression across clinical parameters
# =============================================================================

from itertools import combinations
from collections import defaultdict
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from statannotations.Annotator import Annotator

def SubsetBoxplotter(gene, PECAN_col, do_stats=True, write_file=False, _palette='gray',
                     _dotcolor='white', _fontsize=14, order=None, set_ylim_0=False,
                     list_n=False, sort_mean=False, do_binary=False, col_binary=None, hit_binary=None):

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

    # ======================= #
    # Binary grouping enabled #
    # ======================= #
    if do_binary and col_binary and hit_binary:
        group_labels = []
        for i, sample in enumerate(matching_patient_ids):
            val = clin_lookup.loc[sample, col_binary] if sample in clin_lookup.index else None
            group = hit_binary if val == hit_binary else 'Other'
            group_labels.append(group)

        data = pd.DataFrame({
            'Subtype': group_labels,
            'Expression': values
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
    # ========================= #
    # Default (multi-group) mode #
    # ========================= #
    else:
        subtype_dict = defaultdict(list)
        for i, sample in enumerate(matching_patient_ids):
            if sample in clin_lookup.index:
                subtype = clin_lookup.loc[sample, PECAN_col]
                subtype_dict[subtype].append(values[i])

        data = pd.DataFrame([
            {'Subtype': subtype, 'Expression': expression}
            for subtype, expressions in subtype_dict.items()
            for expression in expressions
        ])

        data = data[~data['Subtype'].isin(['.', 'Unevaluable'])]
        data = data.dropna(subset=['Expression', 'Subtype'])
        data['Subtype'] = data['Subtype'].astype(str)

        if sort_mean:
            mean_order = data.groupby('Subtype')['Expression'].mean().sort_values(ascending=False).index.tolist()
            order = mean_order
        elif order:
            data['Subtype'] = pd.Categorical(data['Subtype'], categories=order, ordered=True)

        sample_counts = data['Subtype'].value_counts().to_dict()
        if list_n:
            data['Subtype_Labeled'] = data['Subtype'].map(lambda x: f"{x}\n(n={sample_counts.get(x, 0)})")
        else:
            data['Subtype_Labeled'] = data['Subtype']

        label_order = data.groupby('Subtype_Labeled')['Expression'].mean().sort_values(ascending=False).index.tolist() \
            if sort_mean else (
            [f"{cat}\n(n={sample_counts.get(cat, 0)})" if list_n else cat for cat in order] if order else data['Subtype_Labeled'].unique())

    # ========== #
    # Plotting   #
    # ========== #
    plt.figure(figsize=(6, 6), dpi=200)
    sns.boxplot(data=data, x='Subtype_Labeled', y='Expression', palette=_palette,
                showfliers=False, order=label_order)
    sns.stripplot(data=data, x='Subtype_Labeled', y='Expression', facecolor=_dotcolor,
                  edgecolor='black', linewidth=0.8, alpha=0.2, jitter=True, order=label_order)

    if do_stats and data['Subtype'].nunique() >= 2:
        pairs = list(combinations(data['Subtype'].dropna().unique(), 2))
        annotator = Annotator(plt.gca(), pairs, data=data, x='Subtype', y='Expression')
        annotator.configure(test='t-test_ind', text_format='star', loc='inside', verbose=2)
        annotator.hide_non_significant = True
        annotator.apply_and_annotate()

    plt.title(f'PeCan expression of {gene} grouped by: \"{PECAN_col}\"', fontsize=_fontsize)
    plt.xlabel(PECAN_col, fontsize=_fontsize)
    plt.ylabel('Expression (VST)', fontsize=_fontsize)
    plt.xticks(rotation=90 if data['Subtype'].nunique() > 4 else 0, fontsize=_fontsize)
    plt.yticks(fontsize=_fontsize)
    if set_ylim_0:
        plt.ylim(0)
    plt.tight_layout()

    if write_file:
        out_path = os.path.join(out_dir, f"{gene}_{PECAN_col}_{_palette}.svg")
        print(f"Saved to: {out_path}")
        plt.savefig(out_path)

    plt.show()



# Example usage with custom order
clin_col   = 'Subtype' #Classifying Driver, ETP.STATUS, Sex, Race, CNS.Status, Insurance, Treatment.Arm, Subtype, Subsuptype, IP Status
gene       = 'KDM6B' # The gene whose expression you want to track
palette    = 'pastel'  # The colors used in the graph. Choose from: https://www.practicalpythonfordatascience.com/ap_seaborn_palette
dotcolor   = 'white' # The colors of the dots on top of the boxplots
fontsize   = 16 # The size of the text items
# order      = ['ETP', 'Near-ETP', 'Non-ETP'] # Specify the order. Set to None or make sure the items are represented in the clin_col
set_ylim_0 = False # Force the 2nd axis to include 0
write_file = False # Write the graph to a file. Will be written to out_dir
do_stats   = False # Perform a statistical analysis and include asterisks in the plot
list_n     = False # provide the number in each category
sort_mean  = True
# do_binary  = True
col_binary = 'Subtype'
hit_binary = 'ETP-like'
# col_binary = 'CNS.Status',
# hit_binary = 'CNS 1'

SubsetBoxplotter(gene, clin_col, do_stats=do_stats, write_file=write_file, _palette=palette, _dotcolor=dotcolor, _fontsize=16, set_ylim_0=set_ylim_0, list_n=list_n, sort_mean=sort_mean, do_binary=True, col_binary=col_binary, hit_binary=hit_binary)
# SubsetBoxplotter(gene, clin_col, do_stats=do_stats, write_file=write_file, _palette=palette, _dotcolor=dotcolor, _fontsize=16, order=order, set_ylim_0=set_ylim_0)
# SubsetBoxplotter(gene, clin_col, do_stats=False, write_file=True, order=['ETP', 'Near-ETP', 'Non-ETP', 'Unknown'], list_n=list_n)

# SubsetBoxplotter('AKT1',clin_col, do_stats=True, write_file=True)

#%% 6b Creating a plot for all categories for a set of genes:
clin_cols = ['Classifying Driver', 'ETP.STATUS', 'Sex', 'Race', 'CNS.Status', 'Insurance', 'Treatment.Arm', 'Subtype', 'Subsuptype', 'IP Status']
# clin_cols = ['ETP.STATUS', ]
genes     = ['KDM6B']
for cc in clin_cols:
    for gene in genes:
        # SubsetBoxplotter(gene, cc,True,True)
        SubsetBoxplotter(gene, cc, False, False)


#%% ===========================================================================
#  7. Run this to create a Kaplan Meier plot of event-free survival for one gene
# =============================================================================
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

def KaplanMeier(_gene):
    # Ensure matching Patient_IDs
    matched_clin_df = clin_df[clin_df['Patient_ID'].isin(df_gexp.columns)]
    
    # Extract expression values for the selected gene
    gene_expression = df_gexp.set_index("Gene").loc[_gene, matched_clin_df['Patient_ID']]
    
    # Compute the median expression
    median_expression = gene_expression.median()
    
    # Split into High and Low expression groups
    matched_clin_df["Expression_Group"] = ["High" if gene_expression[pid] > median_expression else "Low" for pid in matched_clin_df["Patient_ID"]]
    
    # Extract survival data
    time_high = matched_clin_df.loc[matched_clin_df["Expression_Group"] == "High", "EFS"]
    event_high = matched_clin_df.loc[matched_clin_df["Expression_Group"] == "High", "EFS.status"]
    time_low = matched_clin_df.loc[matched_clin_df["Expression_Group"] == "Low", "EFS"]
    event_low = matched_clin_df.loc[matched_clin_df["Expression_Group"] == "Low", "EFS.status"]
    
    # Initialize Kaplan-Meier fitters
    kmf_high = KaplanMeierFitter()
    kmf_low = KaplanMeierFitter()
    
    # Fit data for both groups
    kmf_high.fit(time_high, event_high, label="High Expression")
    kmf_low.fit(time_low, event_low, label="Low Expression")
    
    # Perform log-rank test
    logrank_p = logrank_test(time_high, time_low, event_high, event_low).p_value
    
    # Plot survival curves
    plt.figure(figsize=(6,5),dpi=200)
    plt.ylim(0)
    kmf_high.plot_survival_function()
    kmf_low.plot_survival_function()
    plt.title(f"Kaplan-Meier Survival by {_gene} Expression (p value: %.4f)" %(logrank_p))
    plt.xlabel("Days (Event-free survival)")
    plt.ylabel("Survival Probability")
    plt.legend()
    plt.grid()
    WriteFile(os.path.join(out_dir, '%s_KaplanMeier.svg' %(_gene)))
    plt.show()


# Define the gene of interest
gene = "NAMPT"  # Change this to any gene from df_gexp
KaplanMeier(gene)

#%% ===========================================================================
# 8. Kaplan Meier plot for all clinical parameters
# =============================================================================
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
import matplotlib.pyplot as plt

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
    plt.title(f"Kaplan-Meier Survival by {clin_column}", fontsize=13)
    plt.xlabel("Days (Event-Free Survival)", fontsize=11)
    plt.ylabel("Survival Probability", fontsize=11)

    for group in unique_groups:
        mask = df[clin_column] == group
        n_patients = mask.sum()
        label = f"{group} (n={n_patients})"
        kmf.fit(df.loc[mask, 'EFS'], df.loc[mask, 'EFS.status'], label=label)
        kmf.plot_survival_function(ci_show=False)

    # Perform multivariate log-rank test
    results = multivariate_logrank_test(df['EFS'], df[clin_column], df['EFS.status'])
    pval = results.p_value

    plt.grid()
    plt.legend(title=clin_column, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.ylim(0)
    plt.tight_layout()
    plt.title(f"KM Survival by {clin_column} (p={pval:.4f})")
    WriteFile(os.path.join(out_dir, f'{clin_column}_KM.svg'))
    plt.show()



# Run for each clinical column
clin_cols = ['Classifying Driver', 'ETP.STATUS', 'Sex', 'Race', 'CNS.Status', 'Insurance', 
             'Treatment.Arm', 'Subtype', 'Subsubtype', 'IP Status']

for col in clin_cols:
    KaplanMeier_clinical(col)


#%% ===========================================================================
# 9. Plotting levels in cell lines based on MS
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr

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
            print(f"Suggested alternative gene names based on string matching:")
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

# Example usage:
protein_x = 'IGF2BP2'
protein_y = 'TLX3'
Grapher_MSpr1(protein1=protein_x, protein2=protein_y, df_msdataset=df_cell_line_MS)



