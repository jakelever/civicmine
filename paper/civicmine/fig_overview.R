
source('civicmine/dependencies.R')

#setwd(".")
civicmineFilename <- 'civicmine/civicmine_conservative.tsv'
civicmineFilename <- normalizePath(civicmineFilename)

civicmine <- read.table(civicmineFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")

civicmine <- civicmine[grep(";",civicmine$cancer_normalized,fixed=T,invert=T),]
civicmine <- civicmine[grep(";",civicmine$gene_normalized,fixed=T,invert=T),]
civicmine <- civicmine[grep(";",civicmine$drug_normalized,fixed=T,invert=T),]
civicmine <- civicmine[grep("|",civicmine$gene_normalized,fixed=T,invert=T),]

# Remove the entity location columns and unique the rows
nonLocationColumns <- grep("(start|end)",colnames(civicmine),invert=T)
civicmine <- civicmine[,nonLocationColumns]
civicmine <- civicmine[!duplicated(civicmine),]

civicmineQuadsWithPublications <- civicmine[,c('pmid','evidencetype','cancer_normalized','gene_normalized','drug_normalized','variant_normalized')]
civicmineQuadsWithPublications <- civicmineQuadsWithPublications[!duplicated(civicmineQuadsWithPublications),]

civicmineQuadsWithCounts <- plyr::count(civicmineQuadsWithPublications[,c('evidencetype','cancer_normalized','gene_normalized','drug_normalized','variant_normalized')])

paper.prognostic_erbb2_breastcancer <- civicmineQuadsWithCounts[civicmineQuadsWithCounts$evidencetype=='Prognostic' & civicmineQuadsWithCounts$cancer_normalized=='breast cancer' & civicmineQuadsWithCounts$gene_normalized=='ERBB2' & civicmineQuadsWithCounts$variant_normalized=='overexpression','freq']
paper.predisposing_brca1_breastcancer_mutation <- civicmineQuadsWithCounts[civicmineQuadsWithCounts$evidencetype=='Predisposing' & civicmineQuadsWithCounts$cancer_normalized=='breast cancer' & civicmineQuadsWithCounts$gene_normalized=='BRCA1' & civicmineQuadsWithCounts$variant_normalized=='mutation','freq']

paper.mentionCount <- nrow(civicmine)
paper.paperCount <- length(unique(civicmine$pmid))
paper.abstractCount <- length(unique(civicmine[civicmine$section=='abstract','pmid']))
paper.articleCount <- length(unique(civicmine[civicmine$section=='article','pmid']))
paper.biomarkerCount <- nrow(civicmineQuadsWithCounts)
paper.cancerCount <- length(unique(civicmineQuadsWithPublications$cancer_normalized))
paper.geneCount <- length(unique(civicmineQuadsWithPublications$gene_normalized))
paper.drugCount <- length(unique(civicmineQuadsWithPublications$drug_normalized))
paper.sentenceCount <- length(unique(civicmine$sentence))
paper.multiCitationCount <- sum(civicmineQuadsWithCounts$freq==1)
paper.multiCitationPerc <- round(100*paper.multiCitationCount/paper.biomarkerCount,1)

paper.mentionCount <- prettyNum(paper.mentionCount,big.mark=",")
paper.paperCount <- prettyNum(paper.paperCount,big.mark=",")
paper.abstractCount <- prettyNum(paper.abstractCount,big.mark=",")
paper.articleCount <- prettyNum(paper.articleCount,big.mark=",")
paper.biomarkerCount <- prettyNum(paper.biomarkerCount,big.mark=",")
paper.cancerCount <- prettyNum(paper.cancerCount,big.mark=",")
paper.geneCount <- prettyNum(paper.geneCount,big.mark=",")
paper.drugCount <- prettyNum(paper.drugCount,big.mark=",")
paper.sentenceCount <- prettyNum(paper.sentenceCount,big.mark=",")
paper.multiCitationCount <- prettyNum(paper.multiCitationCount,big.mark=",")


