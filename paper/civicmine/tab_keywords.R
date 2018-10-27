
source('civicmine/dependencies.R')

data <- read.table('civicmine/civic_keywords.table.tsv',sep='\t')
colnames(data) <- c('General','Diagnostic','Predictive','Predisposing','Prognostic')

tab_keywords <- tableGrob(data,rows=NULL)

grid.arrange(tab_keywords)
