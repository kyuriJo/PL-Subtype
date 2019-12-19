# PL-Subtype
* Prerequisite : R 2.15.* (limma), Python 3.* (Scipy).
* To run PL-Subtype, You need 1) **two input files** and 2) **pathway of interest**.
* Description of the input files is in **Input Files**.
* If you have your pathways of interest, see **Case1**.
* If you would like to analyze on recommended pathways from DEG analysis, see **Case2**.

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
python enrichment.py [enrichment_res_file]
```
In [enrichment_res_file], recommended pathways are shown in the order of increasing P-value from DEG enrichment test.
Select one of the pathways and continue.

```console
python subtype_analysis.py [pathway name]
```

## Result files
subtype_analysis.py will generate two files:
1. `DEGs_Normal-[subtype]_in_[Pathway].txt`
   - DEG analysis is performed by Limma comparing Normal vs. each subtype. This file includes the list of DEGs in the user-selected pathway(s).

2. `Subtype-specific_TFs_in_[Pathway].txt`
   - For each TF in user-selected pathway, TF analysis is performed using hypergeometric test to find whether the target genes of TF enrich in subtype DEGs. This file includes the hypergeometric test result of each TF.
