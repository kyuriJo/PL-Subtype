# PL-Subtype
* Prerequisite : R 2.15.* (limma), Python 3.* (Scipy).
* To run PL-Subtype, You need 1) **two input files** and 2) **pathway of interest**.
1. If you have your pathways of interest, see **Case1**.
2. If you would like to analyze on recommended pathways from DEG analysis, see **Case2**.

## Case 1.
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
1. DEGs_Normal-[subtype]_in_[Pathway].txt
   - DEG analysis is performed by Limma comparing Normal vs. each subtype. This file includes the list of DEGs in the user-selected pathway(s).

2. Subtype-specific_TFs_in_[Pathway].txt
   - For each TF in user-selected pathway, TF analysis is performed using hypergeometric test to find whether the target genes of TF enrich in subtype DEGs. This file includes the hypergeometric test result of each TF.
