library(plyr)
library(rjson)

dat = read.csv('BalMaterialSum.csv')

flipFall = function(falltype, stimtype){
  if(stimtype=='N'){return(falltype)}
  if(falltype=='B'){return('B')}
  if(falltype=='L'){return('R')}
  if(falltype=='R'){return('L')}
  return(NULL)
}

makeCondDF = function(cnm) {
  clist = fromJSON(file=cnm)
  sspl = strsplit(clist,'_')
  rdf = data.frame(Trial = sapply(sspl, function(x){paste(x[1],x[2],x[3],x[4],x[5],sep='_')}),
                   Rev = sapply(sspl, function(x){x[6]}))
  
  falls = as.character(sapply(as.character(rdf$Trial),
                              function(tnm){return(dat$Falls[which(dat$Trial==tnm)])}))
  
  rdf$Falls = mapply(flipFall,falls,as.character(rdf$Rev))
  rdf$CondType = sapply(sspl, function(x){x[1]})
  
  return(rdf)
}

cond1 = makeCondDF('Cond1.json')
cond1$Condition = 'C1'
cond2 = makeCondDF('Cond2.json')
cond2$Condition = 'C2'
cond3 = makeCondDF('Cond3.json')
cond3$Condition = 'C3'
cond4 = makeCondDF('Cond4.json')
cond4$Condition = 'C4'

conditions = rbind(cond1,cond2,cond3,cond4)
with(conditions,table(Condition,Falls))
with(conditions,table(Condition,CondType))

numobs = function(trnm) {
  no = 0
  if(trnm %in% cond1$Trial) {no = no + 1}
  if(trnm %in% cond2$Trial) {no = no + 1}
  if(trnm %in% cond3$Trial) {no = no + 1}
  if(trnm %in% cond4$Trial) {no = no + 1}
  return(no)
}
dat$NSeen = sapply(as.character(dat$Trial),numobs)
table(dat$NSeen)
with(dat,table(Classification,NSeen))
