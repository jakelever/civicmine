

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

###################
# Extract variant #
###################
civicmine$variant_normalized <- as.character(civicmine$variant_normalized)
civicmine[civicmine$variant_normalized=='','variant_normalized'] <- '[unknown]'

civicmine$variant_withsub <- as.character(civicmine$variant_normalized)
civicmine[civicmine$variant_normalized=='substitution','variant_withsub'] <- as.character(as.matrix(as.data.frame(civicmine[civicmine$variant_normalized=='substitution','variant_id'])))
civicmine$variant_withsub <- as.factor(civicmine$variant_withsub)

civicmine$variant_withcivic <- tolower(str_replace_all(civicmine$variant_withsub, fixed('substitution|'), ''))

subs <- civicmine[civicmine$variant_normalized=='substitution',]
specificSub <- data.frame(do.call('rbind', strsplit(as.character(subs$variant_id),'|',fixed=TRUE)))
civicmine$variant_sub <- ''
civicmine[civicmine$variant_normalized=='substitution','variant_sub'] <- as.character(specificSub$X2)


################
# Map Gene IDs #
################
genemappingFilename <- 'civicmine/gene_with_protein_product.txt'
genemapping <- read.table(genemappingFilename,header=T,sep='\t',quote='',comment.char='')
rownames(genemapping) <- as.character(genemapping$hgnc_id)
genemapping <- genemapping[,'entrez_id',drop=F]

civicmine$gene_entrez_id <- genemapping[as.character(civicmine$gene_id),]

#civicmine[civicmine$variant_normalized=='substitution','variant_withcivic'] <- 'mutation'
civicmine$combined <- paste(civicmine$evidencetype,civicmine$gene_entrez_id,gsub("DOID:","",civicmine$cancer_id),civicmine$drug_normalized,civicmine$variant_withcivic,sep='_')


civicmineEntrezMapped <- civicmine[!is.na(civicmine$gene_entrez_id),]

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
paper.biomarkerCountEntrezMapped <- length(unique(civicmineEntrezMapped$combined))

paper.civicCount <- prettyNum(paper.civicCount,big.mark=",")
paper.biomarkerCountEntrezMapped <- prettyNum(paper.biomarkerCountEntrezMapped,big.mark=",")



biomarkerComparisonPlot <- venn.diagram(
  x = list(CIViC=unique(civicdb$combined) , CIViCmine=unique(civicmineEntrezMapped$combined) ),
  scaled=F,
  fill = c("grey", "white"),
  cat.fontface = 2,
  cat.pos = 0,
  filename=NULL)
biomarkerComparisonPlot <- gTree(children=biomarkerComparisonPlot)

pmidComparisonPlot <- venn.diagram(
  x = list(CIViC=unique(as.character(civicdb$pubmed_id)) , CIViCmine=unique(as.character(civicmine$pmid)) ),
  scaled=F,
  fill = c("grey", "white"),
  cat.fontface = 2,
  cat.pos = 0,
  filename=NULL)
pmidComparisonPlot <- gTree(children=pmidComparisonPlot)

paper.pmidsInCIViC <- length(unique(as.character(civicdb$pubmed_id)))
paper.pmidsInBoth <- length(intersect(unique(as.character(civicdb$pubmed_id)),unique(as.character(civicmine$pmid))))

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

paper.novelGeneCount <- length(setdiff(unique(civicmineEntrezMapped$gene_entrez_id),unique(civicdb$entrez_id)))
paper.novelCancerCount <- length(setdiff(unique(gsub("DOID:","",civicmine$cancer_id)),unique(civicdb$doid)))
paper.novelDrugCount <- length(setdiff(unique(civicmine$drug_normalized),unique(civicdb$drugs)))

paper.novelGeneCount <- prettyNum(paper.novelGeneCount,big.mark=",")
paper.novelCancerCount <- prettyNum(paper.novelCancerCount,big.mark=",")
paper.novelDrugCount <- prettyNum(paper.novelDrugCount,big.mark=",")


