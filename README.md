# PECAN_scanner.py:
The purpose of this script is to explore correlations between genes in the PECAN dataset (and six cell lines). There are five applications of this script:
1. Replace 'target' in section 4 and run that cell
   - This script will create a waterfall graph with genes ranked based on Pearson's R
   - Breakpoints will be identified using the Kneedle algorithm - All genes with their R- and p-values will be written to a csv together with their position relative to the breakpoints
3. Replace 'target' and 'target2' in section 5 and run that cell
   - Script will create a graph and calculate Pearson's R and associated p-value for the two specified genes
4. Replace 'gene' and 'clin_col' in cell 6 and run that cell
   - Script will create a series of boxplots for the expression levels for patients separated by unique values in a column in the clinical dataset
5. Replace 'gene' in section 7 and run that cell
   - Script will generate a Kaplan-Meier graph for event-free survival for that gene
6. Replace 'protein_x' and 'protein_y' in section 8 and run that cell
   - Script will create a graph and calculate Pearson's R and associated p-value for the two specified proteins across triplicates in six cell lines

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
show_equation    = False
split_by_subtype = False
set_lim_0        = False
subanalysis_do   = False
subanalysis_col  = 'Subtype'
subanalysis_hit  = 'ETP-like'
pval_scientific  = True
top_n_residuals  = 0
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/ccd18204-55b8-4cc4-b245-cded156ad4e5">

### top_n_residuals>0 indicates the patients which lie the closest to the regression line and prints them in the terminal

```python
target  = 'HNRNPC'
target2 = 'COPS4'
show_equation    = False
split_by_subtype = False
set_lim_0        = False
subanalysis_do   = False
subanalysis_col  = 'Subtype'
subanalysis_hit  = 'ETP-like'
pval_scientific  = True
top_n_residuals  = 100
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/c53123a6-a5dc-45ad-b2b2-dadc7412fb4e">

```
Top 100 samples closest to regression line:
Residual range: 0.0003 to 0.0255
Sample IDs: PAVSHE, PAVXPT, PAUGLY, PARLIB, PATZYE, PAVXBD, PAVHVN, PAURFC, PAVKZV, PARUCV, PASMHF, PARMJA, PASJIY, PATXNK, PAVESV, PARIHY, PAWLPU, PAWAWJ, PAWAYT, PASVJM, PAVKWV, PAUGAA, PAWEEX, PAUXKE, PAUKWJ, PAUFZC, PATSDS, PAUXZN, PAWIXS, PASPJR, PASXUL, PAUEIL, PASKLK, PAVGSF, PATKZH, PASRAF, PATRAB, PASWNU, PAUYIK, PAVFWF, PAUAIZ, PAUDPH, PAVMGW, PATGVX, PAWCWY, PAWBMY, PAUCIB, PASILW, PAWKKF, PAUSEY, PATEFF, PAVPEN, PARJFH, PAUEYF, PATWHB, PAVFIC, PASZTU, PAWGLE, PATEBX, PATGZA, PAWCWL, PAUIFR, PAUMLV, PATBYK, PAVLMX, PAVCRK, PARFLD, PATDRC, PAWLSM, PARFWJ, PATYYJ, PAVYTZ, PAUXAI, PAWFBT, PAVVHL, PATZWM, PAWIBU, PAVCEA, PATKRF, PATTTE, PASKIC, PASYCN, PATIHV, PAVMPB, PAUHDN, PAWGDG, PAWGIT, PAUIRF, PATIFI, PATKVR, PASYKJ, PATGLV, PAVRHT, PARXTW, PAUCDC, PAWAIN, PAVFBV, PAWIGC, PAUHAF, PAWALP
```

### split_by_subtypes=True produces plots per subtype, show_equation=True displays slope and intersect, set_lim_0=True forces origo inclusion

```python
target  = 'HNRNPC'
target2 = 'COPS4'
show_equation    = True
split_by_subtype = True
set_lim_0        = True
subanalysis_do   = True
subanalysis_col  = 'Subtype'
subanalysis_hit  = 'ETP-like'
pval_scientific  = True
top_n_residuals  = 0
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/c0e7a459-fa7a-4360-ac35-5ef5de930904">


## subanalysis_do=True only shows includes data points for a subset of patients that satisfy criteria in the clinical dataframe:

```python
subanalysis_do=True
subanalysis_col='ETP status'
subanalysis_hit='ETP'
```
<img width="450" alt="image" src="https://github.com/user-attachments/assets/02d20bd9-f1f3-405b-af4b-19a53773a81c">


Available options for subanlysis_col and subanalysis_hit:
```
Classifying Driver
   MLLT10, TLX3, LYL1, HOXA13, HOXA9, BCL11B, ZFP36L2, KMT2A, Unknown, LMO1, NKX2-5, NUP98, NUP214, MYB, ETV6, SPI1, TLX1, TAL2, NKX2-1, MED12, TAL1, LMO2
ETP.STATUS
   Near-ETP, Unknown, ETP, Non-ETP
Sex
   Female, Male
Race
   Native Hawaiian or other Pacific Islander, Reported Unknown, Black or African American, Asian, American Indian or Alaska Native, White
CNS.Status
   CNS 3c, CNS 3, CNS 3a, CNS 2b, CNS 2a, CNS 2c, CNS 3b, Unknown, CNS 2, CNS 1
Insurance
   PRIVATE INSURANCE, UNKNOWN, MEDICAID, MILITARY SPONSORED (INCLUDING CHAMPUS & TRICARE), MEDICARE, nan, NO MEANS OF PAYMENT (NO INSURANCE), SELF PAY (NO INSURANCE), OTHER, MEDICAID AND MEDICARE, MEDICARE AND PRIVATE INSURANCE
Treatment.Arm
   Arm C, Arm B, nan, Arm D, Standard Induction, Arm A
Subtype
   MLLT10, TLX3, STAG2&LMO2, KMT2A, NKX2-1, HOXA9 TCR, TME-enriched, NUP98, SPI1, NKX2-5, TAL1 DP-like, NUP214, TAL1 αβ-like, TLX1, LMO2 γδ-like, ETP-like, BCL11B
Subsubtype
   TLX3 DP-like, TLX3 Immature, HOXA9 ETP-like, TAL1 DP-like&RPL10, NKX2-1 Other, NKX2-1 TCR, Unknown, HOXA13 ETP-like, Other ETP-like
IP Status
   IP_LR_ETP-like, IP_LR_Myeloid-like, IP_LR_DP-like, Unknown, IP_LR_αβ-like
```

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

## Kaplan-Meier plots of all clinical parameters (Section 8):
e.g.
<img width="450" alt="image" src="https://github.com/user-attachments/assets/2e6fdf5a-64d1-47fb-937b-062c20a192e9">

## Correlation plot of protein levels in six cell lines (Section 9):

```
protein_x = "EZH2"
protein_y = "IGF2BP2"
Grapher_MSpr1(protein1=protein_x, protein2=protein_y, df_msdataset=df_cell_line_MS)
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/27774327-2e35-463a-bbb2-db82d8d0a644">



