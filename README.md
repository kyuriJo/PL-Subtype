# PL-Subtype
* Prerequisite : R 2.15.* (limma), Python 3.* (Scipy).
* To run PL-Subtype, You need 1) **two input files** and 2) **pathway of interest**.
* Description of the input files is in **Input Files**.
* If you have your pathways of interest, see **Case1**.
* If you would like to analyze on recommended pathways from DEG analysis, see **Case2**.
* After generating two result files, follow instructions in **Subtype analysis**.

## Input Files
1. Gene expression file
   * Gene expression file should be tab-delimited.
   * Example:
```
Sample1  Sample2	...	SampleN
Gene1	8	11	...	19
...	19	21	...	11
Gene2	14	20	...	10
```
2. Sample information file
   * Sample information file should be tab-delimited.
   * First column: Sample name, second column: subtype name
   * Normal samples should be included
```
Sample1	Normal
Sample2	Normal
...	...
SampleN-1	Basal
SampleN	Basal
```

## Case 1.
* Select your pathway of interest in the list below:
* **Csf1, Egf, Hgf, Ifnab, Ifng, Igf1, IL1, IL2, IL4, IL6, IL12, Ins, Ngf, Pdg, Tgfb, Tnf, Vegf, Cd40, Serum, Lps, PolyIC, TLR9, Adriamycin, Aniso, Bleomycin, Etoposide, Hydroxyurea, NCS, PMA, Sorbitol, UV, IRad**
* Description of the pathways can be found here: http://pl.csl.sri.com/stm7-guide.html

```console
Rscript run_limma.R [expression_file] [sample_file]
python subtype_analysis.py [pathway name]
```
## Case 2. 
```console
Rscript run_limma.R [expression_file] [sample_file]
python enrichment.py [sample_file]
```
enrichment.py generates `res/Enrichment-test_results.txt` where recommended pathways are shown in the order of increasing P-value from DEG enrichment test.
Select one of the pathways and continue.

```console
python subtype_analysis.py [subtype name] [pathway name]
```

## Result files (in `res` folder)
run_limma.R will generate two types of files:
1. `Limma_Normal-[subtype].txt`
2. `Limma_[subtype]-Others.txt`

subtype_analysis.py will generate two files:
1. `DEGs_in_[subtype]-[Pathway].txt`
   - This file includes the list of DEGs in the user-selected pathway(s).

2. `Subtype-specific_TFs_in_[subtype]-[Pathway].txt`
   - For each TF in user-selected pathway, TF analysis is performed using hypergeometric test to find whether the target genes of TF enrich in subtype DEGs. This file includes the hypergeometric test result of each TF.
   
## Subtype analysis
Install Java and download PLA Online, and proceed analysis with result files.
Java : https://www.oracle.com/technetwork/java/javase/downloads/index.html
PLA Online : http://pl.csl.sri.com/online.html

1. Launch PLA Online and select your pathway of interest.
2. Select Occurences tab in the lower right panel and repeat the following:
   1. Find Occurences that match with DEGs in your `DEGs_in_[subtype]-[Pathway].txt` file
   2. Click the node in the graph panel and this will open the context tab in the lower right panel.
   3. Click on the 'make occ as a goal' button and the node will turn green.
3. Press the button 'Subnet' in the tool bar. This will show you the subtype-specific paths to DEGs.
4. Among the proteins/genes in the subtype-specific path, find TFs that passed the enrichment test in `Subtype-specific_TFs_in_[subtype]-[Pathway].txt`.

For Step 2, please refer to 'Subnets and Pathnets' section of http://pl.csl.sri.com/stm7-guide.html

