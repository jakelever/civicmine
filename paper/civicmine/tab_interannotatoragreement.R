source('civicmine/dependencies.R')

biomarker1 <- read.table('civicmine/interannotatoragreement_biomarker1.tsv',sep='\t',header=T,row.names = 1)
biomarker2 <- read.table('civicmine/interannotatoragreement_biomarker2.tsv',sep='\t',header=T,row.names = 1)
variant <- read.table('civicmine/interannotatoragreement_variant.tsv',sep='\t',header=T,row.names = 1)

tab_biomarker1 <- tableGrob(biomarker1,cols=c("Annotator 2","Annotator 3"))
tab_biomarker2 <- tableGrob(biomarker2,cols=c("Annotator 2","Annotator 3"))
tab_variant <- tableGrob(variant,cols=c("Annotator 2","Annotator 3"))

tab_biomarker1 <- arrangeGrob(tab_biomarker1,top="Biomarker Group 1")
tab_biomarker2 <- arrangeGrob(tab_biomarker2,top="Biomarker Group 2")
tab_variant <- arrangeGrob(tab_variant,top="Variant Group")

tab_interannotatoragreement <- arrangeGrob(tab_biomarker1,tab_biomarker2,tab_variant,ncol=2)

grid.arrange(tab_interannotatoragreement)


