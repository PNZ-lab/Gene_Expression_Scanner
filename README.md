# PECAN_scanner.py:
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

Sidenote: On this github you can find PECAN_CNS_Scanner.py, specifically written to compare levels of expression across different levels of invasion of the central nervous system of patients in the PeCan dataset. Also - some of the functionalities of this script may rely on KTC_functions.py which you can find on this GitHub.

## Graph produced with Section 4 (with accompanying .csv):
- Input one gene and scan all other genes for correlations (Pearson's R) between their expression levels to the target gene. <br>
- Statistic and p-value calculated using Mannwhitney U test.

<img width="450" alt="image" src="https://github.com/user-attachments/assets/7e027f59-77fb-42e4-96fd-a19409ab2db7">
<br>
<img width="450" alt="image" src="https://github.com/user-attachments/assets/c2f8dde5-1570-4936-a130-4129f171e2da">

## Graphs produced for a single comparison (Section 5):
- R and p-values are calculated using Pearson's R.

### Basic comparison of two genes for all patients

```python
target  = 'HNRNPC'
target2 = 'COPS4'
Grapher(target, target2)
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/eec534b5-cbfa-4c65-90f0-6dc93f1323e1">

### split_by_subtypes==True produces plots per subtype, show_equation==True displays slope and intersect

```python
target  = 'HNRNPC'
target2 = 'MYC'
Grapher(target, target2, split_by_subtype=True, show_equation=True)
```


<img width="450" alt="image" src="https://github.com/user-attachments/assets/7ebd4cee-fc17-4941-886c-41cf66b2aa00">


## Graph produced for a single comparison using subanalysis_do=True:
Here:
```python
subanalysis_do=True
subanalysis_col='ETP status'
subanalysis_hit='ETP'
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/02d20bd9-f1f3-405b-af4b-19a53773a81c">

## Expression levels with patients split based on a clinical parameter (Section 6):

### Basic boxplot with gene expression split by ETP status
```
gene     = 'DHFR'
clin_col = 'ETP status'
SubsetBoxplotter(gene, clin_col)
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/9c77e48e-1b53-468f-9cd6-a367b6c34db9">

### perform_statistics and write_file toggle showing significance and saving the plot, _palette changes color scheme, order specificies the order of the boxplots

```
gene     = 'METTL3'
clin_col = 'group'
SubsetBoxplotter(gene, clin_col, perform_statistics=False, write_file=False, _palette='pastel', order=['ETP', 'nearETP', 'notETP'])
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/3f1c704b-4935-4675-8176-2bb70f4e9b3a">

## Kaplan-Meier plots of event-free survival based on gene expression levels (Section 7):

```
gene = "IKZF1"
KaplanMeier(gene)
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/7a87e812-5be0-4191-8c4f-8ddee8b61b0a">




