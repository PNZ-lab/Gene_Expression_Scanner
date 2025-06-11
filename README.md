# Gene_Expression_Scanner.py:
The purpose of this script is to explore correlations between genes in the Polonen and CCLE dataset (and our own proteomics on six cell lines). There are several applications of this script:
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

## Section 1 - Description
This is for the human reading the script itself to parse the ambition and use behind it.

## Section 2 - Setup and Settings
This section reads and formats all the necessary data. The directory that contains all the necessary files are at:
```
cmgg_pnlab/Kasper/Data/Interesting_Lists
```
But in_dir (and potentially the out_dir where files are written to) will need to be modified to match the relative path to that directory for the user launching the script.
Running Section 2 (and Section 3) is mandatory. The analyses will not run without the data loaded into memory (Section 2) - or the functions to analyze the data (Section 3).

## Section 3 - Main Functions
This section contains all the functions for every type of analysis performed by this script. Do not worry, while running this cell is mandatory, no analysis is performed until the functions are called in the sections below. As such running this cell just loads the functions into memory - and is thus instantaneous.

## Section 4 - Graph produced with Section 4
- Input one gene and scan all other genes for correlations (Pearson's R) between their expression levels to the target gene. <br>
- Statistic and p-value calculated using Mannwhitney U test.

<img width="450" alt="image" src="https://github.com/user-attachments/assets/7e027f59-77fb-42e4-96fd-a19409ab2db7">
<br>
<img width="450" alt="image" src="https://github.com/user-attachments/assets/c2f8dde5-1570-4936-a130-4129f171e2da">

## Section 5 - Graphs produced for a single gene-to-gene comparison
- R and p-values are calculated using Pearson's R.

#### Basic comparison of two genes for all patients

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

#### top_n_residuals>0 indicates the patients which lie the closest to the regression line and prints them in the terminal

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

#### split_by_subtypes=True produces plots per subtype, show_equation=True displays slope and intersect, set_lim_0=True forces origo inclusion

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


#### subanalysis_do=True separates data points for a subset of patients that satisfy criteria in the clinical dataframe:

```python
subanalysis_do=True
subanalysis_col='ETP status'
subanalysis_hit='ETP'
```
<img width="450" alt="image" src="https://github.com/user-attachments/assets/02d20bd9-f1f3-405b-af4b-19a53773a81c">


Available options for subanalysis_col and subanalysis_hit:
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

## Section 6 - Expression levels with patients split based on a clinical parameter

#### Basic boxplot with gene expression split by ETP status, order variable should be set to None or a predefined order to appear in graph
```
clin_col   = 'ETP.STATUS'
gene       = 'KDM6B'
palette    = 'pastel'
dotcolor   = 'white'
fontsize   = 16
order      = ['ETP', 'Near-ETP', 'Non-ETP']
set_ylim_0 = False
write_file = False
do_stats   = False
list_n     = False
sort_mean  = False
do_binary  = False
hit_binary = 'Near-ETP'
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/4c209ad5-5832-4f28-9c04-892821e857c0">


#### palette variable changes color scheme, set_ylim_0=True forces 0 to be included on the y-axis, do_stats=True performs a t-test and shows asterisks, list_n=True shows the number of samples, sort_mean=True orders the groups by their means, do_binary=True isolates all samples with hit_binary in clin_col and compares them to the rest

```
clin_col   = 'ETP.STATUS'
gene       = 'KDM6B'
palette    = 'Set2'
dotcolor   = 'white'
fontsize   = 16
order      = None
set_ylim_0 = True
write_file = False
do_stats   = True
list_n     = True
sort_mean  = True
do_binary  = True
hit_binary = 'Near-ETP'
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/2cbf4806-ec97-4175-94d5-070549ba876e">


## Section 7 - Kaplan-Meier plots of event-free survival based on gene expression levels

```
gene = "IKZF1"
KaplanMeier(gene)
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/7a87e812-5be0-4191-8c4f-8ddee8b61b0a">

## Section 8 - Kaplan-Meier plots of all clinical parameters
e.g.
<img width="450" alt="image" src="https://github.com/user-attachments/assets/2e6fdf5a-64d1-47fb-937b-062c20a192e9">

## Section 9 - Correlation plot of protein levels in six cell lines

```
protein_x = "EZH2"
protein_y = "IGF2BP2"
Grapher_MSpr1(protein1=protein_x, protein2=protein_y, df_msdataset=df_cell_line_MS)
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/27774327-2e35-463a-bbb2-db82d8d0a644">


## Section 10 - Correlation plot of gene-to-gene expression levels in CCLE data (cancer cell lines)

```
gene1         = 'RECQL4',
gene2         = 'TONSL',
show_equation = False,
log_scale     = False,
set_lim_0     = False,
filter_col    = 'Hist_Subtype1',
filter_val    = 'acute_lymphoblastic_T_cell_leukaemia'
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/6168d3a8-d762-4f68-9f68-02b21416d246">




## Section 11 - Expression levels of a gene split by clinical parameters

```
CCLE_Boxplotter(
    gene       = 'ZFX',
    group_by   = 'Gender',
    log_scale  = False,
    fig_height = 5,
    fig_width  = 5,
    palette    = 'gray',
    do_stats   = True
    )
```

<img width="450" alt="image" src="https://github.com/user-attachments/assets/466aab4b-c665-48e5-8bdb-3c02306b5ce1">





