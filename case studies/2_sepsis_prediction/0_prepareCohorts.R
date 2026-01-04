


dat=read.csv("./raw_data/ProgimpComputeMCMC.csv")

#remove character and human-created features
remFeat=c("X_biliFlag", "X_fio2Flag", "X_sofaliverFlag", "encounternumber", "c_abe_grp", "c_ph_grp", "c_pco2_grp", 
          "c_leuko_grp", "c_temp_grp", "c_blutz_grp")
dat=dat[,!is.element(colnames(dat), remFeat)]

#check the similarity over imputations
datList=lapply(1:100, function(x)dat[dat$X_Imputation_==x, ])
dim(datList[[1]])
rr=lapply(1:62, function(x){
  m=sapply(1:100, function(y)datList[[y]][,x])
  return(cor(m))
})
sapply(rr, function(x) min(x))
# [1]        NA 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
# [12] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000        NA 1.0000000
# [23] 1.0000000 0.9716977 0.9601065 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 0.9636509 1.0000000 1.0000000
# [34] 1.0000000 1.0000000 1.0000000 0.9777564 0.9630765 0.9697142 1.0000000 0.9602967 0.9216965 0.9907746 1.0000000
# [45] 0.9822059 0.9558862 0.9781459 0.9687485 0.9696126 1.0000000 0.9690628 0.9722369 0.9861927 0.9819349 0.9600723
# [56] 0.9566153 0.9971666 1.0000000 1.0000000 1.0000000 0.9343108 1.0000000

#average the data over imputations
datImp=matrix(0, nrow=401, ncol=62)
for(i in 1:100){
  datImp=datImp + datList[[i]]
}
datImp=datImp/100
datImp$X_Imputation_ = NULL

rownames(datImp) = datImp$encid
datImp$encid=NULL

saveRDS(file="0_datImp.rds", datImp)


dat_cohort1 = datImp[datImp$cohort==1, ]
dat_cohort2 = datImp[datImp$cohort==2, ]
dat_cohort1$cohort=NULL
dat_cohort2$cohort=NULL

saveRDS(file="0_data_cohort1.rds", dat_cohort1)
saveRDS(file="0_data_cohort2.rds", dat_cohort2)

#read onset data, only keep the outcomes
dat=read.csv("./raw_data/DiffimpCompute_model1MCMC.csv")
cols = c("X_Imputation_", "stratum", "encid", "event", "n_krea", "n_harn", "n_laktat")
datOutcomes = dat[,cols]

#average the data over imputations
imputMats = list()
for(i in 1:100){
  print(i)
  subDat = datOutcomes[datOutcomes$X_Imputation_==i, ]
  imputMats[[i]] = subDat
}
outcomes = imputMats[[1]]
for(i in 2:100){
  print(i)
  outcomes = outcomes + imputMats[[i]]
}
outcomes=outcomes/100
outcomes$X_Imputation_=NULL

#remove healthy times points data in cases
outcomesCases = outcomes[outcomes$event==1,]
outcomesContr = outcomes[outcomes$event!=1,]
outcomesContr = outcomesContr[!is.element(outcomesContr$encid, outcomesCases$encid), ]

#average healthy time points data for controls
contIdx = unique(outcomesContr$encid)
outcomesContr2 = vector()
for (i in contIdx ){
  subDat = outcomesContr[outcomesContr$encid==i, , drop=F]
  outcomesContr2 = rbind(outcomesContr2, subDat[median(subDat$stratum), , drop=F])
}
outcomesContr2 = outcomesContr2[,-1]
outcomesCases = outcomesCases[,-1]

final.outcomes = rbind(outcomesCases, outcomesContr2)
saveRDS(final.outcomes, file="0_onset_outcomes.rds")


#match adminsion and onset data
final.outcomes = readRDS("0_onset_outcomes.rds")
dat_cohort1 = readRDS("0_data_cohort1.rds")
dat_cohort2 = readRDS("0_data_cohort2.rds")

out_cohort1 = final.outcomes[match(rownames(dat_cohort1), as.character(final.outcomes[,1])), ]
out_cohort2 = final.outcomes[match(rownames(dat_cohort2), as.character(final.outcomes[,1])), ]
dat_cohort1$event==out_cohort1[,2]
dat_cohort2$event==out_cohort2[,2]
saveRDS(file="0_out_cohort1.rds", out_cohort1)
saveRDS(file="0_out_cohort2.rds", out_cohort2)






