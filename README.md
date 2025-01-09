# PECAN_scanner.py:
This purpose of this script is to explore correlations between genes in the PECAN dataset.
There are two applications of this script:
1. Replace 'target' in section 4 and run that cell
    - This script will create a waterfall graph with genes ranked based on Pearson's R
    - Breakpoints will be identified using the Kneedle algorithm
    - All genes with their R- and p-values will be written to a csv together with their position relative to the breakpoints
2. Replace 'target' and 'target2' in section 5 and run that cell
    - Script will create a graph and calculate Pearson's R and associated p-value for the two specified genes


## Graph produced with Section 4 (with accompanying .csv):
- Input one gene and scan all other genes for correlations (Pearson's R) between their expression levels to the target gene. <br>
- Statistic and p-value calculated using Mannwhitney U test.

<img width="450" alt="image" src="https://github.com/user-attachments/assets/7e027f59-77fb-42e4-96fd-a19409ab2db7">
<br>
<img width="450" alt="image" src="https://github.com/user-attachments/assets/c2f8dde5-1570-4936-a130-4129f171e2da">

## Graph produced for a single comparison:
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
