

source('civicmine/dependencies.R')

#setwd(".")
civicmineFilename <- 'civicmine/civicmine_collated.tsv'
civicmineFilename <- normalizePath(civicmineFilename)
civicmine <- read.table(civicmineFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")

civicmineSentencesFilename <- 'civicmine/civicmine_sentences.tsv'
civicmineSentencesFilename <- normalizePath(civicmineSentencesFilename)
civicmineSentences <- read.table(civicmineSentencesFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")


civicmine$variant_withcivic <- tolower(str_replace_all(civicmine$variant_withsub, fixed(' (substitution)'), ''))
civicmine$combined <- paste(civicmine$evidencetype,civicmine$gene_entrez_id,gsub("DOID:","",civicmine$cancer_id),civicmine$drug_normalized,civicmine$variant_withcivic,sep='_')


################
# Load CIViCdb #
################
civicdbfilename <- 'civicmine/nightly-ClinicalEvidenceSummaries.tsv'
civicdb <- read.table(civicdbfilename,header=T,sep='\t',quote='',comment.char='',encoding='utf8')

civicdbfileInfo <- file.info(civicdbfilename)
paper.civicdbDate <- format(as.Date(civicdbfileInfo$mtime), format="%d %B %Y")

civicdb <- civicdb[,c('pubmed_id','entrez_id','doid','drugs','evidence_type','variant')]
civicdb$drugs <- tolower(civicdb$drugs)
civicdb$variant <- tolower(civicdb$variant)

#looksLikeVariant <- grepl("[a-z][0-9]+.*",civicdb$variant,perl=T)
#civicdb[looksLikeVariant,'variant'] <- 'mutation'
#looksLikeMutation <- grep("mutation", civicdb$variant)
#civicdb[looksLikeMutation,'variant'] <- 'mutation'
#civicdb[civicdb$variant=='deleterious mutation','variant'] <- 'mutation'
#unique(civicdb$variant)

civicdb$combined <- paste(civicdb$evidence_type,civicdb$entrez_id,civicdb$doid,civicdb$drugs,civicdb$variant,sep='_')

paper.civicCount <- length(unique(civicdb$combined))
paper.civicCount <- prettyNum(paper.civicCount,big.mark=",")

biomarkerComparisonPlot <- venn.diagram(
  x = list(CIViC=unique(civicdb$combined) , CIViCmine=unique(civicmine$combined) ),
  scaled=F,
  fill = c("grey", "white"),
  cat.fontface = 2,
  cat.pos = 0,
  filename=NULL)
biomarkerComparisonPlot <- gTree(children=biomarkerComparisonPlot)

pmidComparisonPlot <- venn.diagram(
  x = list(CIViC=unique(as.character(civicdb$pubmed_id)) , CIViCmine=unique(as.character(civicmineSentences$pmid)) ),
  scaled=F,
  fill = c("grey", "white"),
  cat.fontface = 2,
  cat.pos = 0,
  filename=NULL)
pmidComparisonPlot <- gTree(children=pmidComparisonPlot)

paper.pmidsInCIViC <- length(unique(as.character(civicdb$pubmed_id)))
paper.pmidsInBoth <- length(intersect(unique(as.character(civicdb$pubmed_id)),unique(as.character(civicmineSentences$pmid))))

paper.pmidsInCIViC <- prettyNum(paper.pmidsInCIViC,big.mark=",")
paper.pmidsInBoth <- prettyNum(paper.pmidsInBoth,big.mark=",")


#listInput <- list(CIViC=unique(civicdb$combined),CIViCmine=)
#upset(fromList(listInput), order.by = "freq")
#grid.edit('arrange',name='arrange2')
#biomarkerComparisonPlot <- grid.grab()
biomarkerComparisonPlot <- grid.arrange(biomarkerComparisonPlot,top='(a)')



#listInput <- list(CIViC=unique(as.character(civicdb$pubmed_id)),CIViCmine=unique(as.character(civicmine$pmid)))
#upset(fromList(listInput), order.by = "freq")
#grid.edit('arrange',name='arrange2')
#pmidComparisonPlot <- grid.grab()
pmidComparisonPlot <- grid.arrange(pmidComparisonPlot,top='(b)')

fig_comparisons <- arrangeGrob(biomarkerComparisonPlot,pmidComparisonPlot,nrow=1)
grid.arrange(fig_comparisons)

paper.novelGeneCount <- length(setdiff(unique(civicmine$gene_entrez_id),unique(civicdb$entrez_id)))
paper.novelCancerCount <- length(setdiff(unique(gsub("DOID:","",civicmine$cancer_id)),unique(civicdb$doid)))
paper.novelDrugCount <- length(setdiff(unique(civicmine$drug_normalized),unique(civicdb$drugs)))

paper.novelGeneCount <- prettyNum(paper.novelGeneCount,big.mark=",")
paper.novelCancerCount <- prettyNum(paper.novelCancerCount,big.mark=",")
paper.novelDrugCount <- prettyNum(paper.novelDrugCount,big.mark=",")


