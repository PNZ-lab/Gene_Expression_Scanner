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
subanalysis_do   = False # Turns on/off the coloring and separate Pearson's R calculation of data based on a clinical parameter
log_scale        = False
show_breakpoint  = False # Identify and label breakpoints with Kneedle
subanalysis_col  = 'CNS_at_Dx' 
subanalysis_hit  = 'CNS 3' 
#%% ===========================================================================
# 3. Functions
# =============================================================================

def WriteFile(name):
	out_dir_target = os.path.join(out_dir, target)
	if not os.path.exists(out_dir_target):
		os.makedirs(out_dir_target)
	if name.endswith('png'):
		plt.savefig(os.path.join(out_dir_target, name))
		print('file created: %s' %(os.path.join(out_dir_target, name)))

# This function is called elsewhere to create pairwise correlation plots between two genes
def Grapher(gene1, gene2):
    values1 = df_PECAN.loc[df_PECAN['Gene'] == gene1].iloc[0, 1:].tolist()
    values2 = df_PECAN.loc[df_PECAN['Gene'] == gene2].iloc[0, 1:].tolist()
    if log_scale:
        values1 = [value + 1 for value in values1]
        values2 = [value + 1 for value in values2]

    pecan_samples = df_PECAN.columns[1:].tolist()

    if subanalysis_do:
        etp_indices = []
        sample_colors = []
        for i, sample in enumerate(pecan_samples):
            # Check if the sample exists in the clinical data
            match = clin_df[clin_df['RNAseq_id_D'] == sample]
            if not match.empty:
                stage = match[subanalysis_col].values[0]
                if stage == subanalysis_hit:
                    sample_colors.append('red')
                    etp_indices.append(i)
                else:
                    sample_colors.append('black')
            else:
                sample_colors.append('black')
        values1_etp = np.array(values1)[etp_indices]
        values2_etp = np.array(values2)[etp_indices]
        if len(values1_etp) > 1:
            if log_scale:
                # Perform regression in log-transformed space
                log_values1_etp = np.log10(values1_etp)
                log_values2_etp = np.log10(values2_etp)
                r_value_etp, p_value_etp = pearsonr(log_values1_etp, log_values2_etp)
                Xs_etp = log_values1_etp.reshape(-1, 1)
                Ys_etp = log_values2_etp
                model_etp = LinearRegression()
                model_etp.fit(Xs_etp, Ys_etp)
                Y_pred_etp = model_etp.predict(Xs_etp)
            else:
                r_value_etp, p_value_etp = pearsonr(values1_etp, values2_etp)
                Xs_etp = values1_etp.reshape(-1, 1)
                Ys_etp = values2_etp
                model_etp = LinearRegression()
                model_etp.fit(Xs_etp, Ys_etp)
                Y_pred_etp = model_etp.predict(Xs_etp)
    else:
        sample_colors = ['black' for sample in pecan_samples]

    if log_scale:
        # Perform regression in log-transformed space
        log_values1 = np.log10(values1)
        log_values2 = np.log10(values2)
        r_value, p_value = pearsonr(log_values1, log_values2)
        Xs = log_values1.reshape(-1, 1)
        Ys = log_values2
        model = LinearRegression()
        model.fit(Xs, Ys)
        Y_pred = model.predict(Xs)

        # Back-transform Xs for plotting
        sorted_indices = np.argsort(log_values1)
        sorted_log_values1 = log_values1[sorted_indices]
        sorted_Y_pred = Y_pred[sorted_indices]
        sorted_original_values1 = 10**sorted_log_values1
        sorted_Y_pred_back = 10**sorted_Y_pred
    else:
        Xs = np.array(values1).reshape(-1, 1)
        Ys = np.array(values2)
        r_value, p_value = pearsonr(values1, values2)
        model = LinearRegression()
        model.fit(Xs, Ys)
        Y_pred = model.predict(Xs)
        sorted_indices = np.argsort(values1)
        sorted_original_values1 = np.array(values1)[sorted_indices]
        sorted_Y_pred_back = Y_pred[sorted_indices]

    ax, fig = plt.subplots(figsize=(8, 8), dpi=200)
    plt.plot(sorted_original_values1, sorted_Y_pred_back, color='black', label='All: R=%.2f, p=%f' % (r_value, p_value))
    if subanalysis_do:
        if log_scale:
            plt.plot(10**log_values1_etp, 10**Y_pred_etp, color='red', label='%s: R=%.2f, p=%f' % (subanalysis_hit, r_value_etp, p_value_etp))
        else:
            plt.plot(values1_etp, Y_pred_etp, color='red', label='%s: R=%.2f, p=%f' % (subanalysis_hit, r_value_etp, p_value_etp))
    plt.scatter(values1, values2, color=sample_colors, alpha=0.3)
    if log_scale:
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel(gene1 + ' (FPKM+1)', fontsize=18)
        plt.ylabel(gene2 + ' (FPKM+1)', fontsize=18)
    else:
        plt.xlabel(gene1 + ' (FPKM)', fontsize=18)
        plt.ylabel(gene2 + ' (FPKM)', fontsize=18)
    plt.tick_params(axis='both', labelsize=16)
    plt.title('PECAN expression: %s v %s' % (gene1, gene2), fontsize=22)
    plt.legend(fontsize=16)
    file_name = 'PECAN_correlation_%s_v_%s.png' % (gene1, gene2)
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
label = 'WP_KETOGENESIS_AND_KETOLYSIS_custom'
#KTC_GetGeneSet can take: a list of gene names, a single gene name, or the name of a gene set from Msigdb
gene_set = KTC_GetGeneSet(['CPT2','ACAT1','OXCT1','BDH1','SLC2A1','SLC16A1','UCP2','SLC25A20','HMGC2','HMGCL', 'OXCT1','CDK6'])


#Targets is a list of gene names who will have correlations for all other genes calculated (top hits will be )
targets = ['CDK4']

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
target  = 'DHFR'
target2 = 'NAMPT'
Grapher(target, target2)
