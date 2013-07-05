#
# Pre-processing Data
#
#exprData <- read.delim("pnas_expression.txt")
# remove final "len"
exprData <- exprData[,-9]
row.names(exprData)<-exprData$ensembl_ID
exprData<-exprData[,-1]
# store gene with valid expression data under ≥ 1 conditions
atLeastOne<-apply(exprData,1,function(row) any (row!=0))
exprData.fil<-exprData[atLeastOne,]
colnames(exprData.fil)<-c(paste("Control",1:4,sep="_"),paste("DHT",1:3,sep="_"))

#
# Basic
#
exprData.fil$ControlMean<-apply(exprData.fil,1,function(x){mean(x[1:4])})
exprData.fil$DHTMean<-apply(exprData.fil,1,function(x){mean(x[5:7])})
effect<-log2(exprData.fil$DHTMean/exprData.fil$ControlMean)
exprData.fil$Effect<-effect
# t.test cannot operate when
# 1) zero variance
# 2) observation number is less than 1
# solution 1
my.t.test.p.value<-function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
#effect.sig<-apply(exprData.fil,1,function(row)my.t.test.p.value(row[1:4],row[5:7]))
#
# solution 2
effect.sig<-apply(exprData.fil,1,function(row)tryCatch(t.test(row[1:4],row[5:7])$p.value,error=function(row) NA))
#
#wilcox.test is another statistical test
#effect.sig<-apply(exprData.fil,1,function(row)wilcox.test(row[1:4],row[5:7])$p.value)
exprData.fil$Effect.sig<- (-log10(effect.sig))
exprData.fil$Mean<- log2(exprData.fil$ControlMean * exprData.fil$DHTMean)/2
#
# MA-plot
#
library(ggplot2)
p<-ggplot(exprData.fil,aes(Mean,Effect))+geom_point(aes(colour=Effect.sig),alpha=0.7)+scale_colour_gradient(low="red",high="green")
p+xlab(expression(1/2 * log[2](bar(Ctrl) %*% bar(DHT))))+ylab(expression(log[2](bar(DHT) / bar(Ctrl))))
#
# Volcano Plot
#
v<-ggplot(exprData.fil,aes(Effect,Effect.sig))+geom_point(aes(colour=Effect.sig),alpha=0.7)+scale_colour_gradient(low="red",high="green")
v+xlab(expression(log[2](bar(DHT) / bar(Ctrl))))+ylab(expression(-log[10]("p.value")))