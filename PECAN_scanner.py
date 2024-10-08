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

'''

'/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/contrast_Taz_v_ctrl.tsv'
'/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/contrast_ATAC_Taz_v_ctrl.tsv'
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

#Initialization
PECAN_in  = '/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/Expression_Pecan.txt' # PECAN FPKM data
clin_data = '/Users/kasperthorhaugechristensen/Library/CloudStorage/OneDrive-UGent/PNG lab/Kasper/PECAN_scanner/PECAN_clinical_reference.tsv'
out_dir   = '/Users/kasperthorhaugechristensen/Desktop/Dumpbox' # Directory where files and images are written. Subdirectories for individual genes are created
df_PECAN  = pd.read_csv(PECAN_in, sep='\t')
clin_df   = pd.read_csv(clin_data, sep='\t')
top_n     = 10 # E.g. 10 will generate graphs for the 5 most positive and most negatively correlated genes

#%%
#Optional
print_corr_genes = False # Genes above and below breakpoints can be written directly to console
write_files      = True # Turns on/off the writing of csv and pngs
subanalysis_do   = True # Turns on/off the coloring and separate R calculation of data based on a clinical parameter
subanalysis_col  = 'group'
subanalysis_hit  = 'NKX2_1'
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
	values1 = df_PECAN.loc[df_PECAN['Gene'] == gene1].iloc[0,1:].tolist()
	values2 = df_PECAN.loc[df_PECAN['Gene'] == gene2].iloc[0,1:].tolist()

	pecan_samples = df_PECAN.columns[1:].tolist()
# 	print(pecan_samples)
# 	print(clin_df['RNAseq_id_D'].tolist())
	if subanalysis_do:
		etp_indices   = []
		sample_colors = []
		for i, sample in enumerate(pecan_samples):
			# Check if the sample exists in the clinical data
			match = clin_df[clin_df['RNAseq_id_D'] == sample]
			if not match.empty:
# 				print('match found')
				stage = match[subanalysis_col].values[0]
# 				print(stage)
				if stage == subanalysis_hit:
# 					print('\tSubanalysis hit')
					sample_colors.append('red')
					etp_indices.append(i)
				else:
					sample_colors.append('black')
			else:
				# If the sample doesn't exist, default to black
				sample_colors.append('black')
		values1_etp = np.array(values1)[etp_indices]
		values2_etp = np.array(values2)[etp_indices]
		if len(values1_etp) > 1:
			r_value_etp, p_value_etp = pearsonr(values1_etp, values2_etp)
			Xs_etp = values1_etp.reshape(-1, 1)
			Ys_etp = values2_etp
			model_etp = LinearRegression()
			model_etp.fit(Xs_etp, Ys_etp)
			Y_pred_etp = model_etp.predict(Xs_etp)
	else:
		sample_colors = ['black' for sample in pecan_samples]

	Xs = np.array(values1).reshape(-1,1)
	Ys = np.array(values2)

	r_value, p_value = pearsonr(values1, values2)
	model = LinearRegression()
	model.fit(Xs, Ys)
	Y_pred = model.predict(Xs)

	ax, fig = plt.subplots(figsize=(8,8), dpi=200)
	plt.plot(values1, Y_pred, color='black', label='R=%f, p=%f' %(r_value, p_value))
	if subanalysis_do:
		plt.plot(values1_etp, Y_pred_etp, color='red', label='%s: R=%f, p=%f' % (subanalysis_hit, r_value_etp, p_value_etp))
	plt.scatter(values1, values2, color=sample_colors, alpha=0.3)
	plt.xlabel(gene1 + ' (FPKM)', fontsize=18)
	plt.ylabel(gene2 + ' (FPKM)', fontsize=18)
	plt.tick_params(axis='both', labelsize=16)
	plt.title('PECAN: %s v %s' %(gene1, gene2), fontsize=22)
	plt.legend(fontsize=16)
	file_name = 'PECAN_correlation_%s_v_%s.png' %(gene1, gene2)
	WriteFile(file_name)

# This plot is called to create the waterfall plot and return genes above and below breakpoints (for csv output)
def WaterfallPlot(dictionary, target_gene, gene_set):
	genes, r_p_values = zip(*dictionary)
	r_values = [r for r, p in r_p_values]
	ranks = np.arange(1, len(genes) + 1)
	fig, ax = plt.subplots(figsize=(8,5), dpi=200)
	plt.scatter(ranks, r_values, color='black', s=2, zorder=2,label='gene expression correlations')

	#
	highlight_genes = gene_sets[gene_set]
	if len(highlight_genes) != 0:
		highlight_indices = [i for i, gene in enumerate(genes) if gene in highlight_genes]
		highlight_ranks = np.array(ranks)[highlight_indices]
		highlight_r_values = np.array(r_values)[highlight_indices]
		background_r_values = [r for r in r_values if r not in highlight_r_values]
		print(len(highlight_r_values), len(background_r_values))
		stat, p_value = mannwhitneyu(highlight_r_values, background_r_values, alternative='two-sided')
		for rank, r_value in zip(highlight_ranks, highlight_r_values):
			plt.vlines(x=rank, ymin=min(r_values), ymax=max(r_values), color='red', linewidth=0.2, zorder=1, label= '%s (%.0f, p=%f)' %(gene_set, stat, p_value))

	plt.xlabel('Rank of gene to gene correlation', fontsize=16)
	plt.ylabel('Pearson\'s R', fontsize=16)
	plt.title('%s : Waterfall plot of gene correlations' %(target_gene), fontsize=22)
	kn_positive = KneeLocator(ranks, r_values, curve='convex', direction='decreasing')
	kn_negative = KneeLocator(ranks, r_values, curve='concave', direction='decreasing')

	genes_above_elbow = [genes[i] for i in range(kn_positive.knee)]
	genes_below_elbow = [genes[i] for i in range(kn_negative.knee, len(genes))]
# 	ax.text(400, 0, ' genes above: %i' %(len(genes_above_elbow)), fontsize=12)
# 	ax.text(14000, 0, 'genes below: %i ►' %(len(genes_below_elbow)), fontsize=12)
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
def top_n_comparisons(gene, gene_set):
	gene_values = df_PECAN.loc[df_PECAN['Gene'] == gene].iloc[0, 1:].values
	r_dict = {}
	all_genes = df_PECAN['Gene']

	for other_gene in tqdm(all_genes, desc='Comparing genes', unit='gene'):
		if other_gene != gene:
			other_gene_values = df_PECAN.loc[df_PECAN['Gene'] == other_gene].iloc[0, 1:].values
			r_value, pearson_p = pearsonr(gene_values, other_gene_values)
			r_dict[other_gene] = (r_value, pearson_p)

	sorted_genes = sorted(r_dict.items(), key=lambda item: item[1][0], reverse=True)
	genes_above_elbow, genes_below_elbow = WaterfallPlot(sorted_genes, gene, gene_set)

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

# Gene sets whose genes will be highlighted on the Waterfall Plot
# Leave gene_set below as 'None' to ignore
gene_sets = {
	'None' : [],
	'Test' : ['ZCCHC8'],
	'm6a readers' : ['HNRNPC', 'YTHDF1', 'YTHDF2', 'YTHDF3', 'YTHDC1', 'YTHDC2', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3', 'HNRNPA2B1'],
	'splicing factors' : ['SFRS7', 'CUGBP1', 'DAZAP1', 'CUGBP2', 'FMR1', 'A2BP1', 'RBFOX2', 'HNRNPA0', 'HNRNPA1', 'HNRNPA2B1', 'HNRNPC', 'HNRNPC', 'HNRNPC', 'HNRNPD', 'HNRNPD', 'HNRPDL', 'PCBP1', 'PCBP2', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPH3', 'PTBP1', 'HNRNPK', 'HNRNPK', 'HNRNPL', 'HNRPLL', 'HNRNPM', 'FUS', 'HNRNPU', 'TRA2A', 'TRA2B', 'ELAVL2', 'ELAVL4', 'ELAVL1', 'KHSRP', 'MBNL1', 'NOVA1', 'NOVA2', 'PTBP2', 'SFPQ', 'RBM25', 'RBM4', 'KHDRBS1', 'SF3B1', 'SFRS2', 'SF1', 'SFRS1', 'KHDRBS2', 'KHDRBS3', 'SFRS3', 'SFRS9', 'SFRS13A', 'SFRS5', 'SFRS11', 'SFRS6', 'SFRS4', 'TARDBP', 'TIA1', 'TIAL1', 'YBX1', 'ZRANB2', 'ELAVL3', 'RBM5', 'SYNCRIP', 'HNRNPA3', 'QKI', 'RBMX', 'SRRM1', 'ESRP1', 'ESRP2'], # From SpliceAid
	'EpiFactors' : ['A1CF', 'ACINU', 'ACTB', 'ACTL6A', 'ACTL6B', 'ACTR3B', 'ACTR5', 'ACTR6', 'ACTR8', 'ADNP', 'AEBP2', 'AICDA', 'AIRE', 'ALKBH1', 'ALKBH1', 'ALKBH4', 'ALKBH5', 'ANKRD32', 'ANP32A', 'ANP32B', 'ANP32E', 'APBB1', 'APEX1', 'APOBEC1', 'APOBEC2', 'APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D', 'APOBEC3F', 'APOBEC3G', 'APOBEC3H', 'ARID1A', 'ARID1B', 'ARID2', 'ARID4A', 'ARID4B', 'ARNTL', 'ARRB1', 'ASF1A', 'ASF1B', 'ASH1L', 'ASH2L', 'ASXL1', 'ASXL2', 'ASXL3', 'ATAD2', 'ATAD2B', 'ATF2', 'ATF7IP', 'ATM', 'ATN1', 'ATR', 'ATRX', 'ATXN7', 'ATXN7L3', 'AURKA', 'AURKB', 'AURKC', 'BABAM1', 'BAHD1', 'BANP', 'BAP1', 'BARD1', 'BAZ1A', 'BAZ1B', 'BAZ2A', 'BAZ2B', 'BCOR', 'BCORL1', 'BMI1', 'BPTF', 'BRCA1', 'BRCA2', 'BRCC3', 'BRD1', 'BRD2', 'BRD3', 'BRD4', 'BRD7', 'BRD8', 'BRD9', 'BRDT', 'BRE', 'BRMS1', 'BRMS1L', 'BRPF1', 'BRPF3', 'BRWD1', 'BRWD3', 'BUB1', 'C11orf30', 'C14orf169', 'C17orf49', 'CARM1', 'CBLL1', 'CBX1', 'CBX2', 'CBX3', 'CBX4', 'CBX5', 'CBX6', 'CBX7', 'CBX8', 'CCDC101', 'CDC6', 'CDC73', 'CDK1', 'CDK17', 'CDK2', 'CDK3', 'CDK5', 'CDK7', 'CDK9', 'CDY1', 'CDY1B', 'CDY2A', 'CDY2B', 'CDYL', 'CDYL2', 'CECR2', 'CELF1', 'CELF2', 'CELF3', 'CELF4', 'CELF5', 'CELF6', 'CENPC', 'CHAF1A', 'CHAF1B', 'CHD1', 'CHD1L', 'CHD2', 'CHD3', 'CHD4', 'CHD5', 'CHD6', 'CHD7', 'CHD8', 'CHD9', 'CHEK1', 'CHRAC1', 'CHTOP', 'CHUK', 'CIR1', 'CIT', 'CLNS1A', 'CLOCK', 'CRB2', 'CREBBP', 'CSNK2A1', 'CSRP2BP', 'CTBP1', 'CTBP2', 'CTCF', 'CTCFL', 'CTR9', 'CUL1', 'CUL2', 'CUL3', 'CUL4A', 'CUL4B', 'CUL5', 'CXXC1', 'DAPK3', 'DAXX', 'DDB1', 'DDB2', 'DDX17', 'DDX21', 'DDX5', 'DDX50', 'DEK', 'DHX9', 'DMAP1', 'DNAJC1', 'DNAJC2', 'DND1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DNMT3L', 'DNTTIP2', 'DOT1L', 'DPF1', 'DPF2', 'DPF3', 'DPPA3', 'DPY30', 'DR1', 'DTX3L', 'DZIP3', 'E2F6', 'EED', 'EEF1AKMT3', 'EEF1AKMT4', 'EEF1AKNMT', 'EHMT1', 'EHMT2', 'EID1', 'EID2', 'EID2B', 'EIF4A3', 'ELP2', 'ELP3', 'ELP4', 'ELP5', 'ELP6', 'ENY2', 'EP300', 'EP400', 'EPC1', 'EPC2', 'ERBB4', 'ERCC6', 'EXOSC1', 'EXOSC2', 'EXOSC3', 'EXOSC4', 'EXOSC5', 'EXOSC6', 'EXOSC7', 'EXOSC8', 'EXOSC9', 'EYA1', 'EYA2', 'EYA3', 'EYA4', 'EZH1', 'EZH2', 'FAM175A', 'FAM175B', 'FBL', 'FBRS', 'FBRSL1', 'FOXA1', 'FOXO1', 'FOXP1', 'FOXP2', 'FOXP3', 'FOXP4', 'FTO', 'GADD45A', 'GADD45B', 'GADD45G', 'GATAD1', 'GATAD2A', 'GATAD2B', 'GFI1', 'GFI1B', 'GLYR1', 'GSE1', 'GSG2', 'GTF2I', 'GTF3C4', 'HAT1', 'HCFC1', 'HCFC2', 'HDAC1', 'HDAC10', 'HDAC11', 'HDAC2', 'HDAC3', 'HDAC4', 'HDAC5', 'HDAC6', 'HDAC7', 'HDAC8', 'HDAC9', 'HDGF', 'HDGFL2', 'HELLS', 'HIF1AN', 'HINFP', 'HIRA', 'HIRIP3', 'HJURP', 'HLCS', 'HLTF', 'HMG20A', 'HMG20B', 'HMGB1', 'HMGN1', 'HMGN2', 'HMGN3', 'HMGN4', 'HMGN5', 'HNRNPU', 'HNRPL', 'HNRPM', 'HP1BP3', 'HR', 'HSFX3', 'HSPA1A', 'HSPA1A', 'HSPA1B', 'HSPA1B', 'HUWE1', 'IKBKAP', 'IKZF1', 'IKZF3', 'ING1', 'ING2', 'ING3', 'ING4', 'ING5', 'INO80', 'INO80B', 'INO80C', 'INO80D', 'INO80E', 'JADE1', 'JADE2', 'JADE3', 'JAK2', 'JARID2', 'JDP2', 'JMJD1C', 'JMJD6', 'KANSL1', 'KANSL2', 'KANSL3', 'KAT2A', 'KAT2B', 'KAT5', 'KAT6A', 'KAT6B', 'KAT7', 'KAT8', 'KDM1A', 'KDM1B', 'KDM2A', 'KDM2B', 'KDM3A', 'KDM3B', 'KDM4A', 'KDM4B', 'KDM4C', 'KDM4D', 'KDM4E', 'KDM5A', 'KDM5B', 'KDM5C', 'KDM5D', 'KDM6A', 'KDM6B', 'KDM7A', 'KDM8', 'KEAP1', 'KHDRBS1', 'KLF18', 'KMT2A', 'KMT2B', 'KMT2C', 'KMT2D', 'KMT2E', 'L3MBTL1', 'L3MBTL2', 'L3MBTL3', 'L3MBTL4', 'LAS1L', 'LBR', 'LEO1', 'LRWD1', 'MAGOH', 'MAP3K7', 'MAPKAPK3', 'MASTL', 'MAX', 'MAZ', 'MBD1', 'MBD2', 'MBD3', 'MBD4', 'MBD5', 'MBD6', 'MBIP', 'MBNL1', 'MBNL3', 'MBTD1', 'MCRS1', 'MDC1', 'MEAF6', 'MECP2', 'MEN1', 'METTL11B', 'METTL14', 'METTL16', 'METTL21A', 'METTL3', 'METTL4', 'MGA', 'MGEA5', 'MINA', 'MIS18A', 'MIS18BP1', 'MLLT1', 'MLLT10', 'MLLT6', 'MORF4L1', 'MORF4L2', 'MOV10', 'MPHOSPH8', 'MPND', 'MRGBP', 'MSH6', 'MSL1', 'MSL2', 'MSL3', 'MST1', 'MTA1', 'MTA2', 'MTA3', 'MTF2', 'MUM1', 'MYBBP1A', 'MYO1C', 'MYSM1', 'NAA60', 'NAP1L1', 'NAP1L2', 'NAP1L4', 'NASP', 'NAT10', 'NAT10', 'NBN', 'NCL', 'NCOA1', 'NCOA2', 'NCOA3', 'NCOA6', 'NCOR1', 'NCOR2', 'NEK6', 'NEK9', 'NFRKB', 'NFYB', 'NFYC', 'NIPBL', 'NOC2L', 'NPAS2', 'NPM1', 'NPM2', 'NSD1', 'NSL1', 'NSRP1', 'NSUN2', 'NSUN6', 'NTMT1', 'NUP98', 'OGT', 'OIP5', 'PADI1', 'PADI2', 'PADI3', 'PADI4', 'PAF1', 'PAGR1', 'PAK2', 'PARG', 'PARP1', 'PARP2', 'PARP3', 'PAXIP1', 'PBK', 'PBRM1', 'PCGF1', 'PCGF2', 'PCGF3', 'PCGF5', 'PCGF6', 'PCNA', 'PDP1', 'PELP1', 'PHC1', 'PHC2', 'PHC3', 'PHF1', 'PHF10', 'PHF12', 'PHF13', 'PHF14', 'PHF19', 'PHF2', 'PHF20', 'PHF20L1', 'PHF21A', 'PHF8', 'PHIP', 'PIWIL4', 'PKM', 'PKN1', 'POGZ', 'POLE3', 'PPARGC1A', 'PPM1G', 'PPP2CA', 'PPP4C', 'PPP4R2', 'PQBP1', 'PRDM1', 'PRDM11', 'PRDM12', 'PRDM13', 'PRDM14', 'PRDM16', 'PRDM2', 'PRDM4', 'PRDM5', 'PRDM6', 'PRDM7', 'PRDM8', 'PRDM9', 'PRKAA1', 'PRKAA2', 'PRKAB1', 'PRKAB2', 'PRKAG1', 'PRKAG2', 'PRKAG3', 'PRKCA', 'PRKCB', 'PRKCD', 'PRKDC', 'PRMT1', 'PRMT2', 'PRMT5', 'PRMT6', 'PRMT7', 'PRMT8', 'PRMT9', 'PRPF31', 'PRR14', 'PSIP1', 'PTBP1', 'PTBP1', 'PUF60', 'RAD51', 'RAD54B', 'RAD54L', 'RAD54L2', 'RAG1', 'RAG2', 'RAI1', 'RARA', 'RB1', 'RBBP4', 'RBBP5', 'RBBP7', 'RBFOX1', 'RBM11', 'RBM15', 'RBM15B', 'RBM17', 'RBM24', 'RBM25', 'RBM4', 'RBM5', 'RBM7', 'RBM8A', 'RBMY1A1', 'RBX1', 'RCC1', 'RCOR1', 'RCOR3', 'REST', 'RFOX1', 'RING1', 'RLIM', 'RMI1', 'RNF168', 'RNF2', 'RNF20', 'RNF40', 'RNF8', 'RNPS1', 'RPS6KA3', 'RPS6KA4', 'RPS6KA5', 'RPUSD3', 'RRP8', 'RSF1', 'RSRC1', 'RUVBL1', 'RUVBL2', 'RYBP', 'SAFB', 'SAP130', 'SAP18', 'SAP25', 'SAP30', 'SAP30L', 'SATB1', 'SATB2', 'SCMH1', 'SCML2', 'SCML4', 'SENP1', 'SENP3', 'SET', 'SETD1A', 'SETD1B', 'SETD2', 'SETD3', 'SETD5', 'SETD6', 'SETD7', 'SETD8', 'SETDB1', 'SETDB2', 'SETMAR', 'SF3B1', 'SF3B3', 'SFMBT1', 'SFMBT2', 'SFPQ', 'SFSWAP', 'SHPRH', 'SIN3A', 'SIN3B', 'SIRT1', 'SIRT2', 'SIRT6', 'SIRT7', 'SKP1', 'SLU7', 'SMARCA1', 'SMARCA2', 'SMARCA4', 'SMARCA5', 'SMARCAD1', 'SMARCAL1', 'SMARCB1', 'SMARCC1', 'SMARCC2', 'SMARCD1', 'SMARCD2', 'SMARCD3', 'SMARCE1', 'SMEK1', 'SMEK2', 'SMYD1', 'SMYD2', 'SMYD3', 'SMYD4', 'SNAI2', 'SP1', 'SP100', 'SP140', 'SPEN', 'SPOP', 'SRCAP', 'SRRM4', 'SRSF1', 'SRSF10', 'SRSF12', 'SRSF3', 'SRSF6', 'SS18L1', 'SS18L2', 'SSRP1', 'STK4', 'SUDS3', 'SUPT16H', 'SUPT3H', 'SUPT6H', 'SUPT7L', 'SUV39H1', 'SUV39H2', 'SUV420H1', 'SUV420H2', 'SUZ12', 'SYNCRIP', 'TADA1', 'TADA2A', 'TADA2B', 'TADA3', 'TAF1', 'TAF10', 'TAF12', 'TAF1L', 'TAF2', 'TAF3', 'TAF4', 'TAF5', 'TAF5L', 'TAF6', 'TAF6L', 'TAF7', 'TAF8', 'TAF9', 'TAF9B', 'TBL1XR1', 'TDG', 'TDRD3', 'TDRD7', 'TDRKH', 'TET1', 'TET2', 'TET3', 'TEX10', 'TFDP1', 'TFPT', 'THRAP3', 'TLE1', 'TLE2', 'TLE4', 'TLK1', 'TLK2', 'TNP1', 'TNP2', 'TONSL', 'TOP2A', 'TOP2B', 'TP53', 'TP53BP1', 'TRA2B', 'TRIM16', 'TRIM24', 'TRIM27', 'TRIM28', 'TRIM33', 'TRRAP', 'TRUB2', 'TSSK6', 'TTK', 'TYW5', 'U2AF2', 'UBE2A', 'UBE2B', 'UBE2D1', 'UBE2D3', 'UBE2E1', 'UBE2H', 'UBE2N', 'UBE2T', 'UBN1', 'UBR2', 'UBR5', 'UBR7', 'UCHL5', 'UHRF1', 'UHRF2', 'UIMC1', 'USP11', 'USP12', 'USP15', 'USP16', 'USP17L2', 'USP21', 'USP22', 'USP3', 'USP36', 'USP44', 'USP46', 'USP49', 'USP7', 'UTY', 'VDR', 'VIRMA', 'VPS72', 'VRK1', 'WAC', 'WDR5', 'WDR77', 'WDR82', 'WHSC1', 'WHSC1L1', 'WSB2', 'WTAP', 'YAF2', 'YEATS2', 'YEATS4', 'YTHDC1', 'YWHAB', 'YWHAE', 'YWHAZ', 'YY1', 'ZBTB16', 'ZBTB33', 'ZBTB7A', 'ZBTB7C', 'ZC3H13', 'ZCWPW1', 'ZFP57', 'ZGPAT', 'ZHX1', 'ZMYM2', 'ZMYM3', 'ZMYND11', 'ZMYND8', 'ZNF217', 'ZNF516', 'ZNF532', 'ZNF541', 'ZNF592', 'ZNF687', 'ZNF711', 'ZNHIT1', 'ZRANB3', 'ZZZ3']
	}
gene_set = 'splicing factors'

#Targets is a list of gene names who will have correlations for all other genes calculated (top hits will be )
targets = [
	'KDM6B'
	]

for target in targets:
	bottom_genes, top_genes = top_n_comparisons(target, gene_set)
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
target  = 'KDM6B'
target2 = 'SPI1'
Grapher(target, target2)
