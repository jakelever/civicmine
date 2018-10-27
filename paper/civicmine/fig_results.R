
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

evidenceTypeCounts <- plyr::count(civicmineQuadsWithCounts[,'evidencetype',drop=F])

fig_evidenceTypes <- barchart(freq ~ evidencetype, 
                              evidenceTypeCounts, 
                              col="black",
                              xlab="Evidence Type",
                              ylab="# of biomarkers",
                              ylim=c(0,1.1*max(evidenceTypeCounts$freq)))


paper.topCount <- 20

myColours <- brewer.pal(4,"Set2")
my.settings <- list(
  superpose.polygon=list(col=myColours),
  #strip.background=list(col=myColours),
  superpose.line=list(col=myColours),
  strip.border=list(col="black")
)

typeCounts <- plyr::count(civicmineQuadsWithCounts[,c('evidencetype'),drop=F])
typeCounts <- typeCounts[order(typeCounts$freq,decreasing=T),]

geneCounts <- plyr::count(civicmineQuadsWithCounts[,c('gene_normalized'),drop=F])
geneCounts <- geneCounts[order(geneCounts$freq,decreasing=T),]
geneCounts <- geneCounts[1:paper.topCount,]
geneAndTypeCounts <- plyr::count(civicmineQuadsWithCounts[,c('evidencetype','gene_normalized')])
geneAndTypeCounts <- geneAndTypeCounts[geneAndTypeCounts$gene_normalized %in% geneCounts$gene_normalized,]
geneAndTypeCounts$gene_normalized <- factor(as.character(geneAndTypeCounts$gene_normalized),levels=as.character(geneCounts$gene_normalized))
geneAndTypeCounts$evidencetype <- factor(as.character(geneAndTypeCounts$evidencetype),levels=as.character(typeCounts$evidencetype))
topGenesPlot <- barchart(freq ~ gene_normalized,
                         geneAndTypeCounts,
                         groups=evidencetype,
                         stack=T,
                         ylim=c(0,1.1*max(geneCounts$freq)),
                         xlab="Gene",
                         ylab="# of biomarkers",
                         par.settings = my.settings,
                         auto.key=list(space="top", columns=2, 
                                       points=FALSE, rectangles=TRUE),
                         scales=list(x=list(rot=45)))

cancerCounts <- plyr::count(civicmineQuadsWithCounts[,c('cancer_normalized'),drop=F])
cancerCounts <- cancerCounts[order(cancerCounts$freq,decreasing=T),]
cancerCounts <- cancerCounts[1:paper.topCount,]
cancerAndTypeCounts <- plyr::count(civicmineQuadsWithCounts[,c('evidencetype','cancer_normalized')])
cancerAndTypeCounts <- cancerAndTypeCounts[cancerAndTypeCounts$cancer_normalized %in% cancerCounts$cancer_normalized,]
cancerAndTypeCounts$cancer_normalized <- factor(as.character(cancerAndTypeCounts$cancer_normalized),levels=as.character(cancerCounts$cancer_normalized))
cancerAndTypeCounts$evidencetype <- factor(as.character(cancerAndTypeCounts$evidencetype),levels=as.character(typeCounts$evidencetype))
topCancersPlot <- barchart(freq ~ cancer_normalized,
                           cancerAndTypeCounts,
                         groups=evidencetype,
                         stack=T,
                         ylim=c(0,1.1*max(cancerCounts$freq)),
                         xlab="Cancer Type",
                         ylab="# of biomarkers",
                         par.settings = my.settings,
                         scales=list(x=list(rot=45)))

drugCounts <- plyr::count(civicmineQuadsWithCounts[,c('drug_normalized'),drop=F])
drugCounts <- drugCounts[drugCounts$drug_normalized != '',]
drugCounts <- drugCounts[order(drugCounts$freq,decreasing=T),]
drugCounts <- drugCounts[1:paper.topCount,]
drugAndTypeCounts <- plyr::count(civicmineQuadsWithCounts[,c('evidencetype','drug_normalized')])
drugAndTypeCounts <- drugAndTypeCounts[drugAndTypeCounts$drug_normalized %in% drugCounts$drug_normalized,]
drugAndTypeCounts$drug_normalized <- factor(as.character(drugAndTypeCounts$drug_normalized),levels=as.character(drugCounts$drug_normalized))
drugAndTypeCounts$evidencetype <- factor(as.character(drugAndTypeCounts$evidencetype),levels=as.character(typeCounts$evidencetype))
topDrugsPlot <- barchart(freq ~ drug_normalized,
                           drugAndTypeCounts,
                           groups=evidencetype,
                           stack=T,
                           ylim=c(0,1.1*max(drugCounts$freq)),
                           xlab="Drug",
                           ylab="# of biomarkers",
                           par.settings = my.settings,
                           scales=list(x=list(rot=45)))

variantCounts <- plyr::count(civicmineQuadsWithCounts[,c('variant_normalized'),drop=F])
variantCounts$variant_normalized <- as.character(variantCounts$variant_normalized)
variantCounts[variantCounts$variant_normalized == '','variant_normalized'] <- '[unknown]'
variantCounts <- variantCounts[order(variantCounts$freq,decreasing=T),]
variantCounts <- variantCounts[1:paper.topCount,]
variantAndTypeCounts <- plyr::count(civicmineQuadsWithCounts[,c('evidencetype','variant_normalized')])
variantAndTypeCounts$variant_normalized <- as.character(variantAndTypeCounts$variant_normalized)
variantAndTypeCounts[variantAndTypeCounts$variant_normalized == '','variant_normalized'] <- '[unknown]'
variantAndTypeCounts <- variantAndTypeCounts[variantAndTypeCounts$variant_normalized %in% variantCounts$variant_normalized,]
variantAndTypeCounts$variant_normalized <- factor(as.character(variantAndTypeCounts$variant_normalized),levels=as.character(variantCounts$variant_normalized))
variantAndTypeCounts$evidencetype <- factor(as.character(variantAndTypeCounts$evidencetype),levels=as.character(typeCounts$evidencetype))
topVariantsPlot <- barchart(freq ~ variant_normalized,
                            variantAndTypeCounts,
                         groups=evidencetype,
                         stack=T,
                         ylim=c(0,1.1*max(variantCounts$freq)),
                         xlab="Variant Type",
                         ylab="# of biomarkers",
                         par.settings = my.settings,
                         scales=list(x=list(rot=45)))

topGenesPlot <- arrangeGrob(topGenesPlot,top='(a)')
topCancersPlot <- arrangeGrob(topCancersPlot,top='(b)')
topDrugsPlot <- arrangeGrob(topDrugsPlot,top='(c)')
topVariantsPlot <- arrangeGrob(topVariantsPlot,top='(d)')

fig_entitycounts <- arrangeGrob(topGenesPlot,topCancersPlot,topDrugsPlot,topVariantsPlot)
grid.arrange(fig_entitycounts)

civicmine$journalShort <- strtrim(civicmine$journal,35)
journalCounts <- plyr::count(civicmine[,c('journalShort'),drop=F])
journalCounts <- journalCounts[order(journalCounts$freq,decreasing=T),]
journalCounts <- journalCounts[1:paper.topCount,]
journalAndTypeCounts <- plyr::count(civicmine[,c('evidencetype','journalShort')])
journalAndTypeCounts <- journalAndTypeCounts[journalAndTypeCounts$journal %in% journalCounts$journal,]
journalAndTypeCounts$journalShort <- factor(as.character(journalAndTypeCounts$journalShort),levels=as.character(journalCounts$journalShort))
journalAndTypeCounts$evidencetype <- factor(as.character(journalAndTypeCounts$evidencetype),levels=as.character(typeCounts$evidencetype))
topJournalsPlot <- barchart(freq ~ journalShort,
                            journalAndTypeCounts,
                            groups=evidencetype,
                            stack=T,
                            ylim=c(0,1.1*max(variantCounts$freq)),
                            xlab="Journal",
                            ylab="# of biomarkers",
                            par.settings = my.settings,
                            scales=list(x=list(rot=45)))
grid.arrange(topJournalsPlot)

quads <- civicmine[order(civicmine$year,civicmine$month),]
quads$novel <- !duplicated(quads[,c('evidencetype','cancer_normalized','gene_normalized','drug_normalized','variant_normalized')])

yearCounts <- plyr::count(quads[,c('year')])
yearAndNovelCounts <- plyr::count(quads[,c('year','novel')])
yearAndNovelCounts <- yearAndNovelCounts[yearAndNovelCounts$year<2018,]
yearAndNovelCounts$novelTxt <- 'Not Novel'
yearAndNovelCounts$novelTxt[yearAndNovelCounts$novel] <- 'Novel'
yearAndNovelCounts$yearTxt <- as.character(yearAndNovelCounts$year)
labels <- sort(unique(yearAndNovelCounts$year))
labels[labels %% 5 != 0] <- ''
barchart(freq ~ yearTxt, 
         groups=novelTxt,
         stack=T,
         ylim=c(0,1.1*max(yearCounts$freq)),
         par.settings = my.settings,
         scales=list(x=list(rot=45,labels=labels)),
         auto.key=list(space="top", columns=2, 
                       points=FALSE, rectangles=TRUE),
         yearAndNovelCounts)
head(civicmine,1)

#fig_examplesentences <- tableGrob(exampleSentences,row=NULL,cols=c("Evidence Type","Pubmed ID","Sentence"))
#grid.arrange(fig_examplesentences)

bySentences <- plyr::count(civicmine[,'sentence',drop=F])
paper.sentenceCount <- nrow(bySentences)
paper.multiplePerSentence <- sum(bySentences$freq>1)
paper.percSentencesWithMultiple <- round(100*paper.multiplePerSentence/paper.sentenceCount,1)
  
paper.multiplePerSentence <- prettyNum(paper.multiplePerSentence,big.mark=",")
paper.sentenceCount <- prettyNum(paper.sentenceCount,big.mark=",")

expressionRelated <- civicmineQuadsWithCounts[grep("expression",civicmineQuadsWithCounts$variant_normalized),]
nonexpressionRelated <- civicmineQuadsWithCounts[grep("expression",civicmineQuadsWithCounts$variant_normalized,invert=T),]
nonexpressionRelated <- nonexpressionRelated[nonexpressionRelated$variant_normalized != '',]

paper.expressionPrognosticPerc <- round(100.0*sum(expressionRelated$evidencetype=="Prognostic") / nrow(expressionRelated),1)
paper.nonexpressionPrognosticPerc <- round(100.0*sum(nonexpressionRelated$evidencetype=="Prognostic") / nrow(nonexpressionRelated),1)
