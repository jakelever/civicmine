
library(lattice)
library(dplyr)
library(plyr)
library(reshape2)

data <- read.table('civicmine/evaluation_results.txt',header=T,sep='\t')

data <- data[data$correct!='N/A',c('evidencetype','correct','usable','needed')]
data$correct[data$correct==''] <- 'no'
data$usable[data$usable==''] <- 'no'
data$needed[data$needed==''] <- 'no'
data$usable[data$usable=='?'] <- 'maybe'
data$usable[data$usable=='probably'] <- 'maybe'
data$needed[data$needed=='probably'] <- 'maybe'

#barchart(freq ~ needed | evidencetype, summary)

summary <- plyr::count(data[,c('evidencetype','correct')])
wide <- dcast(summary, evidencetype ~ correct, value.var='freq', fill=0)
wide$total <- wide$no + wide$maybe + wide$yes
wide$no <- 100*wide$no / wide$total
wide$maybe <- 100*wide$maybe / wide$total
wide$yes <- 100*wide$yes / wide$total
wide$metric <- 'correct'
correct <- wide


summary <- plyr::count(data[,c('evidencetype','usable')])
wide <- dcast(summary, evidencetype ~ usable, value.var='freq', fill=0)
wide$total <- wide$no + wide$maybe + wide$yes
wide$no <- 100*wide$no / wide$total
wide$maybe <- 100*wide$maybe / wide$total
wide$yes <- 100*wide$yes / wide$total
wide$metric <- 'usable'
usable <- wide

summary <- plyr::count(data[,c('evidencetype','needed')])
wide <- dcast(summary, evidencetype ~ needed, value.var='freq', fill=0)
wide$total <- wide$no + wide$maybe + wide$yes
wide$no <- 100*wide$no / wide$total
wide$maybe <- 100*wide$maybe / wide$total
wide$yes <- 100*wide$yes / wide$total
wide$metric <- 'needed'
needed <- wide

combined <- rbind(correct,usable,needed)

combined$metric <- factor(combined$metric,c('correct','usable','needed'))
combined$evidencetype <- factor(combined$evidencetype,c('Predisposing','Prognostic','Diagnostic','Predictive'))

combined$intermediate <- combined$maybe

fig_evaluationresults <- barchart(no + intermediate + yes ~ metric | evidencetype, combined, 
         stack=T, 
         auto.key=list(columns=3), 
         horizontal=F, 
         xlab="Evaluation metric",
         ylab="Percentage of evaluated evidence items")

grid.arrange(fig_evaluationresults)
