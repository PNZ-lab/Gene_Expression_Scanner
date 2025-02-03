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
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from tqdm import tqdm
from kneed import KneeLocator
from scipy.stats import pearsonr, mannwhitneyu
import seaborn as sns
import os
from KTC_functions import KTC_GetGeneSet


#Initialization
PECAN_in  = '/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/Expression_Pecan.txt' # PECAN FPKM data
PECAN_in  = '/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/Expression_Pecan.txt' # PECAN FPKM data
clin_data = '/Users/kasperthorhaugechristensen/Library/CloudStorage/OneDrive-UGent/PNG lab/Kasper/PECAN_scanner/PECAN_clinical_reference.tsv'
out_dir   = '/Users/kasperthorhaugechristensen/Desktop/Dumpbox' # Directory where files and images are written. Subdirectories for individual genes are created
df_PECAN  = pd.read_csv(PECAN_in, sep='\t')
clin_df   = pd.read_csv(clin_data, sep='\t')
top_n     = 2 # E.g. 10 will generate graphs for the 5 most positive and most negatively correlated genes

#%% Optional
print_corr_genes = False # Genes above and below breakpoints can be written directly to console
write_files      = True # Turns on/off the writing of csv and pngs
log_scale        = False
show_breakpoint  = False # Identify and label breakpoints with Kneedle
subanalysis_do   = False # Turns on/off the coloring and separate Pearson's R calculation of data based on a clinical parameter
subanalysis_col  = 'CNS_at_Dx' 
subanalysis_hit  = 'CNS 3' 
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
def Grapher(gene1, gene2, split_by_subtype=False, show_equation=False):
    values1 = df_PECAN.loc[df_PECAN['Gene'] == gene1].iloc[0, 1:].tolist()
    values2 = df_PECAN.loc[df_PECAN['Gene'] == gene2].iloc[0, 1:].tolist()
    if log_scale:
        values1 = [value + 1 for value in values1]
        values2 = [value + 1 for value in values2]

    pecan_samples = df_PECAN.columns[1:].tolist()
    sample_colors = ['black' for _ in pecan_samples]

    if split_by_subtype:
        match = clin_df[clin_df['RNAseq_id_D'].isin(pecan_samples)]
        unique_subtypes = match['group'].dropna().unique()
        sample_subtypes = {row['RNAseq_id_D']: row['group'] for _, row in match.iterrows()}
        
        for subtype in unique_subtypes:
            indices = [i for i, sample in enumerate(pecan_samples) if sample_subtypes.get(sample) == subtype]
            values1_sub = np.array(values1)[indices]
            values2_sub = np.array(values2)[indices]
            
            if len(values1_sub) > 1:
                fig, ax = plt.subplots(figsize=(8, 8), dpi=200)
                
                if log_scale:
                    log_values1_sub = np.log10(values1_sub)
                    log_values2_sub = np.log10(values2_sub)
                    r_value, p_value = pearsonr(log_values1_sub, log_values2_sub)
                    Xs = log_values1_sub.reshape(-1, 1)
                    Ys = log_values2_sub
                else:
                    r_value, p_value = pearsonr(values1_sub, values2_sub)
                    Xs = values1_sub.reshape(-1, 1)
                    Ys = values2_sub
                
                model = LinearRegression()
                model.fit(Xs, Ys)
                Y_pred = model.predict(Xs)
                a, b = model.coef_[0], model.intercept_
                
                plt.plot(values1_sub, Y_pred, label='%s: R=%.2f, p=%f%s' % (subtype, r_value, p_value, (', y=%.2fx + %.2f' % (a, b) if show_equation else '')), color='black')
                plt.scatter(values1_sub, values2_sub, color='black', alpha=0.3)
                
                if log_scale:
                    plt.yscale('log')
                    plt.xscale('log')
                    plt.xlabel(gene1 + ' (FPKM+1)', fontsize=18)
                    plt.ylabel(gene2 + ' (FPKM+1)', fontsize=18)
                else:
                    plt.xlabel(gene1 + ' (FPKM)', fontsize=18)
                    plt.ylabel(gene2 + ' (FPKM)', fontsize=18)
                plt.tick_params(axis='both', labelsize=16)
                plt.title('PECAN expression: %s v %s (%s)' % (gene1, gene2, subtype), fontsize=22)
                plt.legend(fontsize=16)
                file_name = 'PECAN_correlation_%s_v_%s_%s.svg' % (gene1, gene2, subtype.replace('/','_'))
                WriteFile(file_name)
    
    fig, ax = plt.subplots(figsize=(8, 8), dpi=200)
    if log_scale:
        log_values1 = np.log10(values1)
        log_values2 = np.log10(values2)
        r_value, p_value = pearsonr(log_values1, log_values2)
        Xs = log_values1.reshape(-1, 1)
        Ys = log_values2
    else:
        r_value, p_value = pearsonr(values1, values2)
        Xs = np.array(values1).reshape(-1, 1)
        Ys = np.array(values2)
    
    model = LinearRegression()
    model.fit(Xs, Ys)
    Y_pred = model.predict(Xs)
    a, b = model.coef_[0], model.intercept_
    
    plt.plot(values1, Y_pred, label='All: R=%.2f, p=%f%s' % (r_value, p_value, (', y=%.2fx + %.2f' % (a, b) if show_equation else '')), color='black')
    plt.scatter(values1, values2, color=sample_colors, alpha=0.3)
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
   gene_values = df_PECAN.loc[df_PECAN['Gene'] == gene].iloc[0, 1:].values
   r_dict = {}
   all_genes = df_PECAN['Gene']

   for other_gene in tqdm(all_genes, desc='Comparing genes', unit='gene'):
      if other_gene != gene:
         other_gene_values = df_PECAN.loc[df_PECAN['Gene'] == other_gene].iloc[0, 1:].values
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
targets = ['IGF2BP2']

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
target  = 'HNRNPC'
target2 = 'COPS4'
Grapher(target, target2)

#%% ===========================================================================
# 6. Analyze levels of expression across clinical parameters
# =============================================================================
from statannotations.Annotator import Annotator
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
from collections import defaultdict
from scipy import stats

def SubsetBoxplotter(gene, PECAN_col, perform_statistics=True, write_file=False, _palette='gray'):
    if gene not in df_PECAN['Gene'].values:
        print(f"Gene '{gene}' not found in df_PECAN.")
        return
    
    if PECAN_col not in clin_df.columns:
        print(f"Column '{PECAN_col}' not found in clin_df.")
        return

    pecan_samples = df_PECAN.columns[1:].tolist()
    values = df_PECAN.loc[df_PECAN['Gene'] == gene].iloc[0, 1:].tolist()

    subtype_dict = defaultdict(list)

    clin_lookup = clin_df.set_index('RNAseq_id_D')
    for i, sample in enumerate(pecan_samples):
        if sample in clin_lookup.index:
            subtype = clin_lookup.loc[sample, PECAN_col]
            subtype_dict[subtype].append(values[i])

    data = pd.DataFrame([
        {'Subtype': subtype, 'Expression': expression}
        for subtype, expressions in subtype_dict.items()
        for expression in expressions
    ])

    data = data[~data['Subtype'].isin(['.', 'Unevaluable'])]


    if data.empty:
        print(f"No matching data found for gene {gene}.")
        return

    # # Perform the Shapiro-Wilk test for normality on each group
    # for subtype in data['Subtype'].unique():
    #     group_data = data[data['Subtype'] == subtype]['Expression']
    #     stat, p_value = stats.shapiro(group_data)

    #     print(f"Shapiro-Wilk test for {subtype}:")
    #     print(f"  Test statistic: {stat}, p-value: {p_value}")
    #     if p_value > 0.05:
    #         print(f"  The data for {subtype} is normally distributed.")
    #     else:
    #         print(f"  The data for {subtype} is NOT normally distributed.")

    # Plot the boxplot and stripplot
    plt.figure(figsize=(8, 8), dpi=200)
    sns.boxplot(data=data, x='Subtype', y='Expression', palette=_palette, showfliers=False)
    sns.stripplot(
        data=data,
        x='Subtype',
        y='Expression',
        facecolor='white',
        edgecolor='black',
        linewidth=0.8,
        alpha=0.6,
        jitter=True
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
    
        # Filter out insignificant pairs based on p-value threshold (e.g., p < 0.05)
        significant_pairs = []
        for idx, res in enumerate(results):
            p_value = res.data.pvalue  # Extract p-value from the result
            if p_value < 0.05:  # Only include significant p-values
                significant_pairs.append(pairs[idx])  # Use the pair index to get the corresponding pair

    # Set the plot title and labels
    plt.title(f'PeCan expression of {gene} in patients grouped by column: \"{PECAN_col}\"', fontsize=14)
    plt.xlabel(PECAN_col, fontsize=12)
    plt.ylabel('Expression (FPKM)', fontsize=12)
    plt.xticks(rotation=45, fontsize=10)
    plt.yticks(fontsize=10)

    plt.tight_layout()
    if write_file:
        out_path = os.path.join(out_dir, '%s_%s.svg' %(gene, PECAN_col))
        plt.savefig(out_path)
    plt.show()


# Example usage
clin_col = 'ETP status' # Choose from: 'Maturation stage', 'group', 'Gender', 'Race', 'CNS_at_Dx, 'ETP status'
gene     = 'KDM6B'
SubsetBoxplotter(gene, clin_col, True, False, 'pastel')

#ETP near, not


#%% 6c Creating a plot for all category for a set of genes:
clin_cols = ['Maturation stage', 'group', 'Gender', 'Race', 'CNS_at_Dx', 'ETP status']
genes     = ['METTL3', 'METTL14', 'FTO', 'HNRNPC', 'MYC', 'HMGCS1', 'FDFT1', 'DHCR7']
for cc in clin_cols:
    for gene in genes:
        SubsetBoxplotter(gene, cc, False, True)
