

source('civicmine/dependencies.R')

germlineOrSomatic <- read.table('civicmine/civicmine_germlineOrSomatic.tsv',header=F,sep='\t',quote='',comment.char='',encoding="UTF-8")
colnames(germlineOrSomatic) <- c('evidencetype','gene_entrez_id','gene','variant','inDBSNP','inCOSMIC')

germlineOrSomatic$inDBSNP <-  germlineOrSomatic$inDBSNP == 'True'
germlineOrSomatic$inCOSMIC <-  germlineOrSomatic$inCOSMIC == 'True'

typeCounts <- plyr::count(germlineOrSomatic[,c('evidencetype'), drop=F])
counts <- plyr::count(germlineOrSomatic[,c('evidencetype','inDBSNP','inCOSMIC')])

merged <- merge(counts,typeCounts,by='evidencetype')
colnames(merged) <- c('evidencetype','inDBSNP','inCOSMIC','count','total')
#merged$desc <- paste(merged$inDBSNP,merged$inCOSMIC)
merged$desc <- ''
merged$desc[merged$inDBSNP==F & merged$inCOSMIC==F] <- 'Neither'
merged$desc[merged$inDBSNP==T & merged$inCOSMIC==F] <- 'dbSNP'
merged$desc[merged$inDBSNP==F & merged$inCOSMIC==T] <- 'COSMIC'
merged$desc[merged$inDBSNP==T & merged$inCOSMIC==T] <- 'dbSNP & COSMIC'
merged$desc <- factor(as.character(merged$desc),levels=c('Neither','dbSNP','COSMIC','dbSNP & COSMIC'))
merged$perc <- 100 * merged$count / merged$total


myColours <- brewer.pal(4,"Set3")
my.settings <- list(
  superpose.polygon=list(col=myColours),
  #strip.background=list(col=myColours),
  superpose.line=list(col=myColours),
  strip.border=list(col="black")
)

germlineOrSomaticPlot <- barchart(perc ~ evidencetype, 
         merged,
         stack=T,
         auto.key=list(columns=2),
         xlab='Evidence Type',
         ylab='Percentage',
         par.settings = my.settings,
         groups=desc)
