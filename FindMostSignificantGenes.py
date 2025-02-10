from collections import defaultdict
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from scipy.stats import ttest_ind
import os
from statsmodels.stats.multitest import multipletests

def FindMostSignificantGenes(PECAN_col, group1, group2, _palette='gray', write_file=False, top_n=5, csv_output='gene_statistics.csv'):
    if PECAN_col not in clin_df.columns:
        print(f"Column '{PECAN_col}' not found in clin_df.")
        return
    
    pecan_samples = df_PECAN.columns[1:].tolist()
    clin_lookup = clin_df.set_index('RNAseq_id_D')
    
    results = []
    total_genes = len(df_PECAN['Gene'].values)
    
    for idx, gene in enumerate(df_PECAN['Gene'].values):
        if idx % 100 == 0:
            print(f"Processing gene {idx+1}/{total_genes}...")
        
        values = df_PECAN.loc[df_PECAN['Gene'] == gene].iloc[0, 1:].tolist()
        subtype_dict = defaultdict(list)
        
        for i, sample in enumerate(pecan_samples):
            if sample in clin_lookup.index:
                subtype = clin_lookup.loc[sample, PECAN_col]
                if subtype in [group1, group2]:
                    subtype_dict[subtype].append(values[i])
        
        if group1 in subtype_dict and group2 in subtype_dict:
            mean_group1 = pd.Series(subtype_dict[group1]).mean()
            mean_group2 = pd.Series(subtype_dict[group2]).mean()
            stat, p_value = ttest_ind(subtype_dict[group1], subtype_dict[group2], equal_var=False)
            
            if pd.notna(p_value) and p_value != float('inf') and p_value != -float('inf'):
                results.append((gene, mean_group1, mean_group2, p_value))
    
    if not results:
        print("No genes with sufficient data found.")
        return
    
    results_df = pd.DataFrame(results, columns=['Gene', f'Mean_{group1}', f'Mean_{group2}', 'p_value'])
    
    if len(results_df) == 0:
        print("No valid p-values to adjust.")
        return
    
    _, adj_p_values, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
    results_df['adj_p_value'] = adj_p_values
    results_df.sort_values(by='adj_p_value', inplace=True)
    results_df.to_csv(csv_output, index=False)
    print(f"Results saved to {csv_output}")
    
    top_genes = results_df.head(top_n)['Gene'].tolist()
    for gene in top_genes:
        print(f"Plotting {gene}...")
        values = df_PECAN.loc[df_PECAN['Gene'] == gene].iloc[0, 1:].tolist()
        subtype_dict = defaultdict(list)
        for i, sample in enumerate(pecan_samples):
            if sample in clin_lookup.index:
                subtype = clin_lookup.loc[sample, PECAN_col]
                if subtype in [group1, group2]:
                    subtype_dict[subtype].append(values[i])
        
        data = pd.DataFrame([
            {'Subtype': subtype, 'Expression': expression}
            for subtype, expressions in subtype_dict.items()
            for expression in expressions
        ])
        
        plt.figure(figsize=(8, 8), dpi=200)
        sns.boxplot(data=data, x='Subtype', y='Expression', palette=_palette, showfliers=False, order=[group1, group2])
        sns.stripplot(
            data=data, x='Subtype', y='Expression', facecolor='white', edgecolor='black', 
            linewidth=0.8, alpha=0.6, jitter=True, order=[group1, group2]
        )
        
        pairs = [(group1, group2)]
        annotator = Annotator(plt.gca(), pairs, data=data, x='Subtype', y='Expression')
        annotator.configure(test='t-test_ind', text_format='star', loc='inside', verbose=2)
        annotator.apply_and_annotate()
        
        plt.title(f'PeCan Expression: {gene} (adj_p={results_df[results_df.Gene == gene].adj_p_value.values[0]:.4e})', fontsize=14)
        plt.xlabel(PECAN_col, fontsize=12)
        plt.ylabel('Expression (FPKM)', fontsize=12)
        plt.xticks(rotation=45, fontsize=10)
        plt.yticks(fontsize=10)
        plt.tight_layout()
        
        if write_file:
            out_path = os.path.join(out_dir, f"{gene}_{PECAN_col}.svg")
            plt.savefig(out_path)
        plt.show()

FindMostSignificantGenes('CNS_at_Dx', 'CNS 1', 'CNS 3')
# FindMostSignificantGenes('Gender', 'Male', 'Female')
