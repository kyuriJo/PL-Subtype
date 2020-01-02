library(edgeR)
library(limma)

args = commandArgs(trailingOnly=TRUE)
expr <- read.table(args[1], sep='\t', stringsAsFactors=FALSE)
targets <- readTargets(args[2], sep='\t', row.names=NULL)

dir.create(file.path('res'), showWarnings = FALSE)
subs = unique(targets$Subtype)
for (sub in subs) {
  if (sub == 'Normal') { next; }
  # Run Limma for Normal vs. Subtype
  filter = which(targets$Subtype=="Normal"|targets$Subtype==sub)
  group = targets$Subtype[filter]
  expr_f = expr[,filter]
  dge = DGEList(counts=expr_f, group=group)
  design = model.matrix(~factor(group))
  v = voom(dge, design)
  fit = lmFit(v, design)
  fit = eBayes(fit)
  result = topTable(fit, coef=ncol(design), number=nrow(fit), adjust.method="fdr")
  write.table(result, file=paste0("res/Limma_Normal-", sub, ".txt"), sep="\t")

  # Run Limma for Subtype vs. Others
  filter = which(targets$Subtype!="Normal")
  group = targets$Subtype[filter]
  group[group!=sub]="Others"
  expr_f = expr[,filter]
  dge = DGEList(counts=expr_f, group=group)
  design = model.matrix(~factor(group))
  v = voom(dge, design)
  fit = lmFit(v, design)
  fit = eBayes(fit)
  result = topTable(fit, coef=ncol(design), number=nrow(fit), adjust.method="fdr")
  write.table(result, file=paste0("res/Limma_", sub, "-Others.txt"), sep="\t")
}
