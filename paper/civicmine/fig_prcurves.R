
source('civicmine/dependencies.R')

AssociatedVariant <- read.table('civicmine/prCurves/AssociatedVariant.txt',header=T)
colnames(AssociatedVariant) <- c('threshold','precision','recall')
AssociatedVariant$reltype <- 'AssociatedVariant'

Diagnostic <- read.table('civicmine/prCurves/Diagnostic.txt',header=T)
colnames(Diagnostic) <- c('threshold','precision','recall')
Diagnostic$reltype <- 'Diagnostic'

Predictive <- read.table('civicmine/prCurves/Predictive.txt',header=T)
colnames(Predictive) <- c('threshold','precision','recall')
Predictive$reltype <- 'Predictive'

Predisposing <- read.table('civicmine/prCurves/Predisposing.txt',header=T)
colnames(Predisposing) <- c('threshold','precision','recall')
Predisposing$reltype <- 'Predisposing'

Prognostic <- read.table('civicmine/prCurves/Prognostic.txt',header=T)
colnames(Prognostic) <- c('threshold','precision','recall')
Prognostic$reltype <- 'Prognostic'

data <- rbind(AssociatedVariant,Diagnostic,Predictive,Predisposing,Prognostic)
data <- data[order(data$precision,decreasing=T),]
data <- data[order(data$recall),]

data$reltype <- factor(as.character(data$reltype),levels=c("Predisposing","Prognostic","AssociatedVariant","Diagnostic","Predictive"))
prcurvesPlot <- xyplot(precision ~ recall | reltype, 
       xlab="Recall", ylab="Precision",
       #xlim=c(0,1),ylim=c(0,1),
       data, 
       lwd=3,
       type="l")
prcurvesPlot <- arrangeGrob(prcurvesPlot,top="(a)")

data <- data[order(data$threshold),]

myColours <- c(brewer.pal(3,"Dark2"),"#000000")
my.settings <- list(
  superpose.polygon=list(col=myColours),
  #strip.background=list(col=myColours),
  superpose.line=list(col=myColours),
  strip.border=list(col="black")
)


thresholdPlot <- xyplot(precision + recall ~ threshold | reltype, 
       xlab="Threshold", ylab="Precision / Recall",
       #auto.key=T,
       par.settings = my.settings,
       auto.key=list(space="top", columns=2, 
                     points=FALSE, rectangles=TRUE),
       data, type="l")
thresholdPlot <- arrangeGrob(thresholdPlot,top="(b)")

fig_prcurves <- arrangeGrob(prcurvesPlot,thresholdPlot,ncol=1)

grid.arrange(fig_prcurves)

targetMatching <- data
targetMatching$target <- 0.9
targetMatching[targetMatching$reltype=='AssociatedVariant','target'] <- 0.94
targetMatching$closeToTarget <- (targetMatching$precision-targetMatching$target)^2
targetMatching <- targetMatching[order(targetMatching$closeToTarget),]

thresholdChoice <- targetMatching[!duplicated(targetMatching$reltype),]
thresholdChoice <- thresholdChoice[order(as.character(thresholdChoice$reltype)),]
thresholdChoice$precision <- round(thresholdChoice$precision,3)
thresholdChoice$recall <- round(thresholdChoice$recall,3)

paper.performanceTable <- thresholdChoice[,c('reltype','threshold','precision','recall')]
#fig_performance <- tableGrob(thresholdChoice[,c('reltype','threshold','precision','recall')], 
#                             rows=NULL, 
#                             cols=c("Extracted Relation","Threshold","Precision","Recall"))
#grid.arrange(fig_performance)
          