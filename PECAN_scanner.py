#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%% ===========================================================================
# 1. Description
# =============================================================================
'''
This purpose of this script is to explore correlations between genes in the PECAN dataset.
There are two applications of this script:
	1. Replace 'target' in section 4 and run that cell
		- This script will create a waterfall graph with genes ranked based on Pearson's R
		- Breakpoints will be identified using the Kneedle algorithm
		- All genes with their R- and p-values will be written to a csv together with their position relative to the breakpoints
	2. Replace 'target' and 'target2' in section 5 and run that cell
		- Script will create a graph and calculate Pearson's R and associated p-value for the two specified genes
	3. Replace 'gene' and 'clin_col' in cell 6 and run that cell
		Script will create a series of boxplots for the expression levels for patients separated by unique values in a column in the clinical dataset

Optionally, the script can perform a separate coloring and calculation of an input column and hit using the PECAN_clinical_reference.tsv
	- e.g. (group==TAL1) or (ETP status==ETP)
	- Important: The calculation for the subset is in addition to the normal data - that is the black trendline and R includes the subset

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



#Initialization
#PECAN_in  = '/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/Expression_Pecan.txt' # PECAN FPKM data
#clin_data = '/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/PECAN_clinical_reference.tsv'
print("Loading gene expression data...")
df_gexp                  = pd.read_csv('/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/PeCan_gexp.csv')
print("Loading clinical data...")
df_annot                 = pd.read_csv('/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/PeCan_annot.csv')
df_M0_clinical           = pd.read_csv("/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/PeCan_M0_clinical.csv")
df_M1_classifying_driver = pd.read_csv("/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/PeCan_M1_classifying_driver.csv")
df_M2_ETP_status         = pd.read_csv("/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/PeCan_M2_ETP_status.csv")
df_M3_genetic_subtype    = pd.read_csv("/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/PeCan_M3_genetic_subtype.csv")
df_M3_subsubtype         = pd.read_csv("/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/PeCan_M3_subsubtype.csv")
df_M3_subtype            = pd.read_csv("/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/PeCan_M3_subtype.csv")
df_M4_pathway            = pd.read_csv("/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/PeCan_M4_pathway.csv")
df_M5_allesions_genes    = pd.read_csv("/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/PeCan_M5_Allesions_genes.csv")
df_M5_allesions_variants = pd.read_csv("/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/PeCan_M5_Allesions_variants.csv")
df_M7_IP                 = pd.read_csv("/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/PeCan_M7_IP.csv")

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

out_dir   = '/Users/kasperthorhaugechristensen/Desktop/Dumpbox' # Directory where files and images are written. Subdirectories for individual genes are created


#%%

#df_gexp  = pd.read_csv(PECAN_in, sep='\t')
#clin_df   = pd.read_csv(clin_data, sep='\t')
top_n     = 2 # E.g. 10 will generate graphs for the 5 most positive and most negatively correlated genes

#%% Optional
print_corr_genes = False # Genes above and below breakpoints can be written directly to console
write_files      = True # Turns on/off the writing of csv and pngs
log_scale        = False
show_breakpoint  = False # Identify and label breakpoints with Kneedle
#%% ===========================================================================
# 3. Functions
# =============================================================================

def WriteFile(name):
	print(name)
	out_dir_target = os.path.join(out_dir, target)
	if not os.path.exists(out_dir_target):
		os.makedirs(out_dir_target)
	if name.endswith('png') or name.endswith('svg'):
		plt.savefig(os.path.join(out_dir_target, name))
		print('file created: %s' %(os.path.join(out_dir_target, name)))

# This function is called elsewhere to create pairwise correlation plots between two genes
def Grapher(gene1, gene2, split_by_subtype=False, subanalysis_do=False, subanalysis_col=None, subanalysis_hit=None, show_equation=False):
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
                
                plt.plot(values1_sub, Y_pred_sub, label='%s%s' % (subtype, (', y=%.2fx + %.2f' % (model_sub.coef_[0], model_sub.intercept_) if show_equation else '')), linestyle='dashed')
                plt.xlabel(gene1 + ' (FPKM)', fontsize=18)
                plt.ylabel(gene2 + ' (FPKM)', fontsize=18)
                plt.tick_params(axis='both', labelsize=16)
                plt.title('PECAN expression: %s v %s (%s)' % (gene1, gene2, subtype), fontsize=22)
                plt.legend(fontsize=16)
                file_name = 'PECAN_correlation_%s_v_%s_%s.svg' % (gene1, gene2, subtype.replace('/','_'))
                WriteFile(file_name)
                plt.show()
                plt.close(fig)
    else:
        fig, ax = plt.subplots(figsize=(8, 8), dpi=200)
        plt.scatter(values1, values2, color=sample_colors, alpha=0.2)
        
        if subanalysis_do and len(etp_indices) > 1:
            values1_etp = values1[etp_indices]
            values2_etp = values2[etp_indices]
            
            model_etp = LinearRegression()
            model_etp.fit(values1_etp.reshape(-1, 1), values2_etp)
            Y_pred_etp = model_etp.predict(values1_etp.reshape(-1, 1))
            
            plt.plot(values1_etp, Y_pred_etp, label='%s%s' % (subanalysis_hit, ', y=%.2fx + %.2f' % (model_etp.coef_[0], model_etp.intercept_) if show_equation else ''), color='red')
        
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
        
        model = LinearRegression()
        model.fit(Xs, Ys)
        Y_pred = model.predict(Xs)
        a, b = model.coef_[0], model.intercept_
        
        plt.plot(values1, Y_pred, label='All: R=%.2f, p=%f%s' % (r_value, p_value, (', y=%.2fx + %.2f' % (a, b) if show_equation else '')), color='black')
        
        plt.xlabel(gene1 + ' (FPKM)', fontsize=18)
        plt.ylabel(gene2 + ' (FPKM)', fontsize=18)
        plt.tick_params(axis='both', labelsize=16)
        plt.title('PECAN expression: %s v %s' % (gene1, gene2), fontsize=22)
        plt.legend(fontsize=16)
        file_name = 'PECAN_correlation_%s_v_%s.svg' % (gene1, gene2)
        WriteFile(file_name)




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

   return(sorted_genes[:round(top_n/2)], sorted_genes[-round(top_n/2):])

#%% ===========================================================================
# 4. Scan for best correlations for one gene
# =============================================================================

# # Gene sets whose genes will be highlighted on the Waterfall Plot
#Use KTC_GetGene
# name_gene_set = 'HALLMARK_MYC_TARGETS_V1'
# name_gene_set = 'HALLMARK_ESTROGEN_RESPONSE_EARLY'
label = 'Neuroblastoma Breakpoint Family'
#KTC_GetGeneSet can take: a list of gene names, a single gene name, or the name of a gene set from Msigdb
gene_set = KTC_GetGeneSet(['CPT2','ACAT1','OXCT1','BDH1','SLC2A1','SLC16A1','UCP2','SLC25A20','HMGC2','HMGCL', 'OXCT1','CDK6'])

#Targets is a list of gene names who will have correlations for all other genes calculated (top hits will be )
targets = ['KDM6B']

for target in targets:
    bottom_genes, top_genes = top_n_comparisons(target, gene_set, label)
    if top_n > 0:
        for bottom_gene, _ in bottom_genes:
            Grapher(target, bottom_gene)
        for top_gene, _ in top_genes:
            Grapher(target, top_gene)


#%% ===========================================================================
# 5. Perform a single, specified comparison
# =============================================================================
#Overwrite 'target' and 'target2' abd run this cell
#File is saved in out_dir/[target]
#DHFR, NAMPT, IDO1, NAPRT1
target  = 'KDM6B'
target2 = 'TMEM88'
split_by_subtype = False
subanalysis_do = False
show_equation = True
subanalysis_col = 'ETP.STATUS'
subanalysis_hit = 'ETP'
# Grapher(target, target2, split_by_subtype, subanalysis_do, subanalysis_col, subanalysis_hit,show_equation=False)
Grapher(target, target2, split_by_subtype,subanalysis_do=subanalysis_do, subanalysis_col=subanalysis_col, subanalysis_hit=subanalysis_hit,  show_equation=show_equation)

#%% ===========================================================================
# 6. Analyze levels of expression across clinical parameters
# =============================================================================

def SubsetBoxplotter(gene, PECAN_col, perform_statistics=True, write_file=False, _palette='gray', _dotcolor='white', _fontsize=14, order=None):
    if gene not in df_gexp['Gene'].values:
        print(f"Gene '{gene}' not found in df_gexp.")
        return
    
    if PECAN_col not in clin_df.columns:
        print(f"Column '{PECAN_col}' not found in clin_df.")
        return

    # Use the matching patient IDs based on previous intersection
    matching_patient_ids = set(df_gexp.columns[1:]).intersection(clin_df['Patient_ID'].unique())
    
    # Filter df_gexp columns to only those that match the patient IDs
    filtered_df_gexp = df_gexp.loc[df_gexp['Gene'] == gene, list(matching_patient_ids)]

    # Now extract the values for the matching samples
    values = filtered_df_gexp.iloc[0].tolist()

    subtype_dict = defaultdict(list)

    # Update the clin_lookup to use 'Patient_ID' instead of 'RNAseq_id_D'
    clin_lookup = clin_df.set_index('Patient_ID')
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
    data = data.dropna(subset=['Expression'])

    if data.empty:
        print(f"No matching data found for gene {gene}.")
        return

    # Set the category order if provided
    if order:
        data['Subtype'] = pd.Categorical(data['Subtype'], categories=order, ordered=True)

    # Plot the boxplot and stripplot
    plt.figure(figsize=(8, 8), dpi=200)
    sns.boxplot(data=data, x='Subtype', y='Expression', palette=_palette, showfliers=False, order=order)
    sns.stripplot(
        data=data,
        x='Subtype',
        y='Expression',
        facecolor=_dotcolor,
        edgecolor='black',
        linewidth=0.8,
        alpha=0.2,
        jitter=True,
        order=order
    )

    if perform_statistics:
        # Specify the pairs to compare
        pairs = list(combinations(data['Subtype'].unique(), 2))

        # Create Annotator and calculate p-values using statannotations
        annotator = Annotator(plt.gca(), pairs, data=data, x='Subtype', y='Expression')
        annotator.configure(test='t-test_ind', text_format='star', loc='inside', verbose=2)
        annotator.hide_non_significant = True

        # Apply annotations for significant pairs
        _, results = annotator.apply_and_annotate()

    # Set the plot title and labels
    plt.title(f'PeCan expression of {gene} in patients grouped by column: \"{PECAN_col}\"', fontsize=_fontsize)
    plt.xlabel(PECAN_col, fontsize=_fontsize)
    plt.ylabel('Expression (FPKM)', fontsize=_fontsize)
    plt.xticks(rotation=45, fontsize=_fontsize)
    plt.yticks(fontsize=_fontsize)

    plt.tight_layout()
    if write_file:
        out_path = os.path.join(out_dir, f"{gene}_{PECAN_col}.svg")
        plt.savefig(out_path)
    plt.show()

# Example usage with custom order
clin_col = 'ETP.STATUS' #Classifying Driver, ETP.STATUS, Sex, Race, CNS.Status, Insurance, Treatment.Arm, Subtype, Subsuptype, IP Status
# Choose from: 'Maturation stage', 'group', 'Gender', 'Race', 'CNS_at_Dx', 'ETP status'
gene = 'CEBPZ'
colors = 'pastel'  # Choose from: https://www.practicalpythonfordatascience.com/ap_seaborn_palette
SubsetBoxplotter(gene, clin_col, perform_statistics=False, write_file=True, _palette=colors, _dotcolor='white', _fontsize=16, order=['ETP', 'Near-ETP', 'Non-ETP', 'Unknown'])
# SubsetBoxplotter(gene, clin_col, perform_statistics=False, write_file=True, _palette=colors, order=['ETP', 'Near-ETP', 'Non-ETP', 'Unknown'])

# SubsetBoxplotter(gene,clin_col, perform_statistics=True, write_file=True,_palette=colors, order=None)

#%% 6c Creating a plot for all category for a set of genes:
clin_cols = ['Maturation stage', 'group', 'Gender', 'Race', 'CNS_at_Dx', 'ETP status']
genes     = ['METTL3', 'METTL14', 'FTO', 'HNRNPC', 'MYC', 'HMGCS1', 'FDFT1', 'DHCR7']
for cc in clin_cols:
    for gene in genes:
        SubsetBoxplotter(gene, cc, False, True)

