
source('civicmine/dependencies.R')

#setwd(".")
civicmineFilename <- 'civicmine/civicmine_collated.tsv'
civicmineFilename <- normalizePath(civicmineFilename)
civicmine <- read.table(civicmineFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")

civicmineSentencesFilename <- 'civicmine/civicmine_sentences.tsv'
civicmineSentencesFilename <- normalizePath(civicmineSentencesFilename)
civicmineSentences <- read.table(civicmineSentencesFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")



paper.prognostic_erbb2_breastcancer <- civicmine[civicmine$evidencetype=='Prognostic' & civicmine$cancer_normalized=='breast cancer' & civicmine$gene_normalized=='ERBB2' & civicmine$variant_withsub=='overexpression','citation_count']
paper.predisposing_brca1_breastcancer_mutation <- civicmine[civicmine$evidencetype=='Predisposing' & civicmine$cancer_normalized=='breast cancer' & civicmine$gene_normalized=='BRCA1' & civicmine$variant_group=='mutation','citation_count']

paper.mentionCount <- nrow(civicmine)
paper.paperCount <- length(unique(civicmineSentences$pmid))
paper.abstractCount <- length(unique(civicmineSentences[civicmineSentences$section=='abstract','pmid']))
paper.articleCount <- length(unique(civicmineSentences[civicmineSentences$section=='article','pmid']))
paper.biomarkerCount <- nrow(civicmine)
paper.cancerCount <- length(unique(civicmine$cancer_normalized))
paper.geneCount <- length(unique(civicmine$gene_normalized))
paper.drugCount <- length(unique(civicmine$drug_normalized))
paper.sentenceCount <- length(unique(civicmineSentences$sentence))
paper.multiCitationCount <- sum(civicmine$citation_count>1)
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


