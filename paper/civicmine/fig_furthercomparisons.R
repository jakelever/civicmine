

source('civicmine/dependencies.R')

#setwd(".")
civicmineFilename <- 'civicmine/civicmine_collated.tsv'
civicmineFilename <- normalizePath(civicmineFilename)
civicmine <- read.table(civicmineFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")

civicmine$variant_withcivic <- tolower(str_replace_all(civicmine$variant_withsub, fixed(' (substitution)'), ''))

cgiPredictive <- read.table('civicmine/cgi_predictive.tsv',sep='\t',header=F,comment.char='',encoding="UTF-8")
colnames(cgiPredictive) <- c('gene','variant','cancer_id','cancer','drug_id','drug')
cgiPredictive <- cgiPredictive[cgiPredictive$cancer_id!='*',]

tmpCivicmine <- civicmine[civicmine$evidencetype=='Predictive',]
tmpCivicmine$combined <- paste(tmpCivicmine$gene_normalized,tmpCivicmine$cancer_id,tmpCivicmine$drug_id,tmpCivicmine$variant_withcivic,sep='_')
cgiPredictive$combined <- paste(cgiPredictive$gene,cgiPredictive$cancer_id,cgiPredictive$drug_id,str_to_lower(cgiPredictive$variant),sep='_')

tmpCivicmine$combinedNoVariant <- paste(tmpCivicmine$gene_normalized,tmpCivicmine$cancer_id,tmpCivicmine$drug_id,sep='_')
cgiPredictive$combinedNoVariant <- paste(cgiPredictive$gene,cgiPredictive$cancer_id,cgiPredictive$drug_id,sep='_')


cgiPredictivePlot <- venn.diagram(
  x = list("CGI"=unique(cgiPredictive$combined) , "CIViCmine"=unique(tmpCivicmine$combined) ),
  scaled=F,
  fill = c("grey", "white"),
  cat.fontface = 2,
  cat.dist = vennDist,
  cat.pos = c(vennOffset,-vennOffset),
  filename=NULL)
cgiPredictivePlot <- gTree(children=cgiPredictivePlot)
grid.arrange(cgiPredictivePlot)



cgiPredisposing <- read.table('civicmine/cgi_predisposing.tsv',sep='\t',header=F,comment.char='',encoding="UTF-8")
colnames(cgiPredisposing) <- c('gene','variant','cancer_id','cancer')
cgiPredisposing <- cgiPredisposing[cgiPredisposing$cancer_id!='*',]

tmpCivicmine <- civicmine[civicmine$evidencetype=='Predisposing',]
tmpCivicmine$combined <- paste(tmpCivicmine$gene_normalized,tmpCivicmine$cancer_id,tmpCivicmine$variant_withcivic,sep='_')
tmpCivicmine$combinedNoVariant <- paste(tmpCivicmine$gene_normalized,tmpCivicmine$cancer_id,sep='_')
cgiPredisposing$combined <- paste(cgiPredisposing$gene,cgiPredisposing$cancer_id,str_to_lower(cgiPredisposing$variant),sep='_')
cgiPredisposing$combinedNoVariant <- paste(cgiPredisposing$gene,cgiPredisposing$cancer_id,sep='_')

#uniqueToCGIPredisposing <- cgiPredisposing[!cgiPredisposing$combined %in% tmpCivicmine$combined,]
#write.table(uniqueToCGIPredisposing,'uniqueToCGIPredisposing.tsv',sep='\t',quote=F,row.names=F)

cgiPredisposingPlot <- venn.diagram(
  x = list("CGI"=unique(cgiPredisposing$combinedNoVariant) , "CIViCmine"=unique(tmpCivicmine$combinedNoVariant) ),
  scaled=F,
  fill = c("grey", "white"),
  cat.fontface = 2,
  cat.dist = vennDist,
  cat.pos = c(vennOffset,-vennOffset),
  filename=NULL)
cgiPredisposingPlot <- gTree(children=cgiPredisposingPlot)
#cgiPredisposingPlot <- arrangeGrob(cgiPredisposingPlot, top='Predisposing')
grid.arrange(cgiPredisposingPlot)



oncokbPredictive <- read.table('civicmine/oncoKB_predictive.tsv',sep='\t',header=F,comment.char='',encoding="UTF-8")
colnames(oncokbPredictive) <- c('gene_id','gene','variant','cancer_id','cancer','drug_id','drug')
oncokbPredictive <- oncokbPredictive[oncokbPredictive$cancer_id!='*',]

tmpCivicmine <- civicmine[civicmine$evidencetype=='Predictive',]
tmpCivicmine$combined <- paste(tmpCivicmine$gene_entrez_id,tmpCivicmine$cancer_id,tmpCivicmine$drug_id,tmpCivicmine$variant_withcivic,sep='_')
tmpCivicmine$combinedNoVariant <- paste(tmpCivicmine$gene_entrez_id,tmpCivicmine$cancer_id,tmpCivicmine$drug_id,sep='_')
oncokbPredictive$combined <- paste(oncokbPredictive$gene_id,oncokbPredictive$cancer_id,oncokbPredictive$drug_id,str_to_lower(oncokbPredictive$variant),sep='_')
oncokbPredictive$combinedNoVariant <- paste(oncokbPredictive$gene_id,oncokbPredictive$cancer_id,oncokbPredictive$drug_id,sep='_')

overlap <- oncokbPredictive[oncokbPredictive$combined %in% tmpCivicmine$combined,]

#uniqueToOncoKB <- oncokbPredictive[!oncokbPredictive$combined %in% tmpCivicmine$combined,]
#uniqueToOncoKB$drugInCIViCmine <- uniqueToOncoKB$drug_id %in% civicmine$drug_id
#write.table(uniqueToOncoKB,'uniqueToOncoKB.tsv',sep='\t',quote=F,row.names=F)

#egfrSubset <- civicmine[civicmine$evidencetype=='Predictive' & civicmine$gene_normalized=='EGFR',]

oncoKBPredictivePlot <- venn.diagram(
  x = list("OncoKB"=unique(oncokbPredictive$combined) , "CIViCmine"=unique(tmpCivicmine$combined) ),
  scaled=F,
  fill = c("grey", "white"),
  cat.fontface = 2,
  cat.dist = vennDist,
  cat.pos = c(vennOffset,-vennOffset),
  filename=NULL)
oncoKBPredictivePlot <- gTree(children=oncoKBPredictivePlot)
grid.arrange(oncoKBPredictivePlot)

vennjunk <- dir(path=".", pattern="VennDiagram*") # ?dir
removed <- file.remove(vennjunk) # ?file.remove
