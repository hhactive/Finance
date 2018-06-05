##Load the package
library(plyr)
## import the data
trans <- read.csv("trans.asc", sep = ";")
#Assign False to all 'inRange'
trans[,'inRange']=FALSE
##Assign all trans in [9500,9999] with TRUE
trans[which((trans[,6]>=9500)&(trans[,6]<=9999)),'inRange']<-TRUE
##Based on 'account_id do summarize'
transByAccount <- ddply(trans, "account_id", summarize, fractionInRange  = mean(inRange))
##'acc' take all account which have more than 10% in range.
acc<-transByAccount[which(transByAccount[,2]>=0.1),]
##This function summarise for absolute trans records.
valuesInRange<-function(arr){
  try<-ddply(arr[arr[,'inRange'],], "account_id", summarize, cc=length(unique(amount)))
  return<-try
}
vv<-as.vector(valuesInRange(trans))
##Filter all trans accounts more than 5 records in range.
accff<-vv[which(vv[,'cc']>=5),][,'account_id']
##Get the answer qualified with 2 conditions.
answer<-intersect(accff,acc[,'account_id'])
#Output
print("We suggest these accounts for investigation of structuring.")
print(answer)

