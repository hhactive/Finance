##1-3 import the data
trans <- read.csv("trans.asc", sep = ";")
##4 try sample for get columns 
transid<-trans['trans_id'] 
##5 get the first digits, save them
amm<-trans['amount']
ammc<-c(1:1056320)
for (i in 1: 1056320)
{
  ammc[i]=substr(amm[i,1],1,1)
}

##6  use table function to do frquency 
ammfrq<-as.data.frame(table(as.numeric(ammc)))

##7 Calculate the fraction from 1 to 9
ammtfreq<-c(1:9)
asum<-sum(ammfrq[2:10,2])
for (i in 1:9)
{
  ammtfreq[i]<-ammfrq[i+1,2]/asum
}

##8  import the benford as benchmark
benford<-c(0.301,0.176,0.125,0.097,0.079,0.067,0.058,0.051,0.046)
check<-ammtfreq-benford
##9 compare and make conclusion

##bo<-TRUE
##for (i in 1:9)
##  if (abs(check[i])>0.01) 
##  {
##    bo<-FALSE
##  }
if (max(abs(check))>0.01) {
  print("This dataset looks fake according to Benford's law")
}else
{ 
  print("This dataset looks real according to Benford's law")
}
