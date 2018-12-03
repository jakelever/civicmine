
library(lattice)
library(dplyr)
library(plyr)
library(reshape2)

data <- read.table('civicmine/evaluation_results.txt',header=T,sep='\t')
data <- data[data$correct!='N/A',c('evidencetype','correct','usable','needed')]
colnames(data) <- c('evidencetype','Correct','Usable','Needed')

data$Correct[data$Correct==''] <- 'no'
data$Usable[data$Usable==''] <- 'no'
data$Needed[data$Needed==''] <- 'no'
data$Usable[data$Usable=='?'] <- 'maybe'
data$Usable[data$Usable=='probably'] <- 'maybe'
data$Needed[data$Needed=='probably'] <- 'maybe'

data$Correct <- revalue(data$Correct,c('no'='No','yes'='Yes','maybe'='Maybe'))
data$Usable <- revalue(data$Usable,c('no'='No','yes'='Yes','maybe'='Maybe'))
data$Needed <- revalue(data$Needed,c('no'='No','yes'='Yes','maybe'='Maybe'))

#barchart(freq ~ Needed | evidencetype, summary)

summary <- plyr::count(data[,c('evidencetype','Correct')])
wide <- dcast(summary, evidencetype ~ Correct, value.var='freq', fill=0)
wide$total <- wide$No + wide$Maybe + wide$Yes
wide$No <- 100*wide$No / wide$total
wide$Maybe <- 100*wide$Maybe / wide$total
wide$Yes <- 100*wide$Yes / wide$total
wide$metric <- 'Correct'
Correct <- wide


summary <- plyr::count(data[,c('evidencetype','Usable')])
wide <- dcast(summary, evidencetype ~ Usable, value.var='freq', fill=0)
wide$total <- wide$No + wide$Maybe + wide$Yes
wide$No <- 100*wide$No / wide$total
wide$Maybe <- 100*wide$Maybe / wide$total
wide$Yes <- 100*wide$Yes / wide$total
wide$metric <- 'Usable'
Usable <- wide

summary <- plyr::count(data[,c('evidencetype','Needed')])
wide <- dcast(summary, evidencetype ~ Needed, value.var='freq', fill=0)
wide$total <- wide$No + wide$Maybe + wide$Yes
wide$No <- 100*wide$No / wide$total
wide$Maybe <- 100*wide$Maybe / wide$total
wide$Yes <- 100*wide$Yes / wide$total
wide$metric <- 'Needed'
Needed <- wide

combined <- rbind(Correct,Usable,Needed)

combined$metric <- factor(combined$metric,c('Correct','Usable','Needed'))
combined$evidencetype <- factor(combined$evidencetype,c('Predisposing','Prognostic','Diagnostic','Predictive'))

combined$Intermediate <- combined$Maybe

fig_evaluationresults <- barchart(No + Intermediate + Yes ~ metric | evidencetype, combined, 
         stack=T, 
         auto.key=list(columns=3), 
         horizontal=F, 
         xlab="Evaluation metric",
         ylab="Percentage of evaluated evidence items")

grid.arrange(fig_evaluationresults)
