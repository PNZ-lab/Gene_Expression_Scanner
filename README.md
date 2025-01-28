# PECAN_scanner.py:
This purpose of this script is to explore correlations between genes in the PECAN dataset.
There are three applications of this script:
1. Replace 'target' in section 4 and run that cell
    - This script will create a waterfall graph with genes ranked based on Pearson's R
    - Breakpoints will be identified using the Kneedle algorithm
    - All genes with their R- and p-values will be written to a csv together with their position relative to the breakpoints
2. Replace 'target' and 'target2' in section 5 and run that cell
    - Script will create a graph and calculate Pearson's R and associated p-value for the two specified genes
    - Optionally, this script can perform a separate analysis in the same graph for a subset of the cohort based on clinical data.
3. Replace 'gene' and 'clin_col' in cell 6 and run that cell 
    - Script will create a series of boxplots for the expression levels for patients separated by unique values in a column in the clinical dataset

Sidenote: On this github you can find PECAN_CNS_Scanner.py, specifically written to compare levels of expression across different levels of invasion of the central nervous system of patients in the PeCan dataset.

## Graph produced with Section 4 (with accompanying .csv):
- Input one gene and scan all other genes for correlations (Pearson's R) between their expression levels to the target gene. <br>
- Statistic and p-value calculated using Mannwhitney U test.

<img width="450" alt="image" src="https://github.com/user-attachments/assets/7e027f59-77fb-42e4-96fd-a19409ab2db7">
<br>
<img width="450" alt="image" src="https://github.com/user-attachments/assets/c2f8dde5-1570-4936-a130-4129f171e2da">

## Graph produced for a single comparison (Section 6):
- R and p-values are calculated using Pearson's R.

<img width="450" alt="image" src="https://github.com/user-attachments/assets/93c5ce2d-4f2b-4afa-900b-89d564917d85">

## Graph produced for a single comparison using subanalysis_do=True:
Here:
```python
subanalysis_do=True
subanalysis_col='ETP status'
subanalysis_hit='ETP'
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/02d20bd9-f1f3-405b-af4b-19a53773a81c">

## Graphs produced with Section 6:
```
gene     = 'DHFR'
clin_col = 'ETP status'
SubsetBoxplotter(gene, clin_col)
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/9c77e48e-1b53-468f-9cd6-a367b6c34db9">

```
gene     = 'METTL3'
clin_col = 'group'
SubsetBoxplotter(gene, clin_col, True, False, 'pastel')
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/57f2e61d-ed08-4f5d-8139-7bb871dc08ea">





