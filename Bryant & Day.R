# BRYANT & DAY

####### Bryant & Day modif simulation patients pour AI #####

#N1=28 ; CR1=7 ; CT1=20 ; N2=60 ; CR2=18 ; CT2=46 #seuils trouvés avec la fonction
N1=28 ; CR1=6 ; CT1=19 ; N2=60 ; CR2=17 ; CT2=45 #seuils -1

ntrials=10000
npat=N2

res_etape1=res_etape2=rep(NA,ntrials)
result1=result2=rep(NA,ntrials)
nb_eff=nb_ech=nb_tox=nb_notox=0
sample1=sample2=0
time=rep(NA,ntrials)
t_obs=4

param_eff=c(0.2,0.4,0.4,0.2,0.4,0.2,0.4,0.3)
param_tox=c(0.35,0.15,0.05,0.05,0.45,0.15,0.35,0.2)
scenario=paste("Scenario",1:8)
grid=cbind(data.frame(scenario,param_eff,param_tox))
colnames(grid)=c("Scenario","Peff","Ptox")


for (i in 1:nrow(grid)) {
  
  print(paste("Scenario",i))
  
  peff<-grid$Peff[i]
  ptox<-grid$Ptox[i]
  print(paste("peff=",peff,"ptox=",ptox))
  
  pefftox=peff*ptox
  peffnotox=peff-pefftox
  pnoefftox=ptox-pefftox
  pnoeffnotox=(1-peff)-pnoefftox
  
  res_etape1=res_etape2=rep(NA,ntrials)
  result1=result2=rep(NA,ntrials)
  nb_eff=nb_ech=nb_tox=nb_notox=0
  sample1=sample2=0
  time=rep(NA,ntrials)
  check_stop1=rep(NA,ntrials)
  
  
  set.seed(3003)
  for (n in 1:ntrials) {
    
    #print(paste("essai",n))
    #simulation des patients
    #patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))),cumsum(c(0,rep(2,npat-1)))) #temps fixe: cumsum(c(0,rep(2,npat-1))) temps aléatoire: cumsum(c(0,rpois(npat-1,2)))
    patient=t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox)))
    pois=c(0,rpois(N2-1,2))
    patient_7=pois[7]
    patient_26=pois[26]
    patient_29=pois[29]
    
    if (patient_7==0) {
      while(patient_7==0) {
        patient_7<-rpois(1,2)
      }
      pois[7]<-patient_7
    }
    if (patient_26==0) {
      while(patient_26==0) {
        patient_26<-rpois(1,2)
      }
      pois[26]<-patient_26
    }
    if (patient_29==0) {
      while(patient_29==0) {
        patient_29<-rpois(1,2)
      }
      pois[29]<-patient_29
    }
    patient=cbind(patient,cumsum(pois))
    
    #ETAPE 1
    stage1=patient[1:N1,]
    eff1=sum(stage1[,1])+sum(stage1[,2])
    ech1=sum(stage1[,3])+sum(stage1[,4])
    tox1=sum(stage1[,2])+sum(stage1[,4])
    notox1=sum(stage1[,1])+sum(stage1[,3])
    
    if (eff1 > CR1 && notox1 > CT1) {res_etape1[n]<-"succes"} else {res_etape1[n] <- "echec"}
    sample1=sample1+N1
    nb_eff=nb_eff+eff1
    nb_ech=nb_ech+ech1
    nb_tox=nb_tox+tox1
    nb_notox=nb_notox+notox1
    
    #raison de l'arrêt précoce
    if (eff1<=CR1 && notox1<=CT1) {result1[n]<-"tox:excessive ; rep:inadequate"}
    if (eff1 >CR1 && notox1<=CT1) {result1[n]<-"tox:excessive ; rep:insuf data"}
    if (eff1<=CR1 && notox1 >CT1) {result1[n]<-"tox:insuf data ; rep:inadequate"}
    
    #durée de l'essai
    if (res_etape1[n]=="echec") {time[n]<-patient[N1,5]+t_obs}
    #print(paste(N1,"patients",eff1,"eff",notox1,"notox"))
    #print(paste("Résultat:",res_etape1[n],result1[n]))
    #print(paste("durée",time[n]))
    
    #vérification inclusion AI
    if (patient[N1,5]==patient[N1+1,5]) {check_stop1[n] <- "STOP"} else {"OK"}
    
    #ETAPE 2
    if (res_etape1[n]=="succes") {
      stage2=patient[(N1+1):N2,]
      eff2=sum(stage2[,1])+sum(stage2[,2])
      ech2=sum(stage2[,3])+sum(stage2[,4])
      tox2=sum(stage2[,2])+sum(stage2[,4])
      notox2=sum(stage2[,1])+sum(stage2[,3])
      
      if (eff1+eff2 > CR2 && notox1+notox2 > CT2) {res_etape2[n]<-"succes"} else {res_etape2[n] <- "echec"}
      sample1=sample1+N2-N1
      nb_eff=nb_eff+eff2
      nb_ech=nb_ech+ech2
      nb_tox=nb_tox+tox2
      nb_notox=nb_notox+notox2
      
      #raison de l'échec 
      if (eff1+eff2<=CR2 && notox1+notox2<=CT2) {result2[n]<-"tox:excessive ; rep:inadequate"}
      if (eff1+eff2 >CR2 && notox1+notox2<=CT2) {result2[n]<-"tox:excessive ; rep:OK"}
      if (eff1+eff2<=CR2 && notox1+notox2 >CT2) {result2[n]<-"tox:OK ; rep:inadequate"}
      
      #durée de l'essai
      if (sum(patient[1:N1-1,1])+sum(patient[1:N1-1,2])==CR1 | sum(patient[1:N1-1,1])+sum(patient[1:N1-1,3])==CT1-1) {time[n]<-patient[N2,5]+2*t_obs}
      if (sum(patient[1:N1-1,1])+sum(patient[1:N1-1,2])!=CR1 && sum(patient[1:N1-1,1])+sum(patient[1:N1-1,3])!=CT1-1) {time[n]<-patient[N2,5]+t_obs}
      #print(paste(N2-N1,"patients",eff2,"eff",notox2,"notox"))
      #print(paste("total eff:",eff1+eff2,"total notox:",notox1+notox2))
      #print(paste("Résultat:",res_etape2[n],result2[n]))
      #print(paste("durée",time[n]))
      
      
    }
    
  }
  
  #effectif attendu
  grid$ess[i] <- (sample1+sample2)/ntrials
  
  #nb d'efficacités 
  grid$moy_eff[i] <- nb_eff/ntrials
  
  #nb d'échecs
  grid$moy_ech[i] <- nb_ech/ntrials
  
  #nb de toxicités 
  grid$moy_tox[i] <- nb_tox/ntrials
  
  #durée de l'essai
  grid$moy_time[i] <- mean(time)
  
  #% conclu d'efficacité
  if (all(is.na(res_etape2[res_etape2=="succes"]))) {
    grid$succes[i] <- 0
  } else {
    grid$succes[i] <- ((table(res_etape2[res_etape2=="succes"]))*100)/ntrials
  }
  
  #% arrêt précoce
  if (all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$early_ter[i] <- 0
  } else {
    grid$early_ter[i] <- (table(res_etape1[res_etape1=="echec"])*100)/ntrials
  }
  
  #nb de succès et échecs à chaque étape
  if (all(is.na(res_etape1[res_etape1=="succes"]))) {
    grid$res1_succes[i] <- 0
  } else {
    grid$res1_succes[i] <- table(res_etape1[res_etape1=="succes"])
  }
  if (all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$res1_echecs[i] <- 0
  } else {
    grid$res1_echecs[i] <- table(res_etape1[res_etape1=="echec"])
  }
  
  if (all(is.na(res_etape2[res_etape2=="succes"]))) {
    grid$res2_succes[i] <- 0
  } else {
    grid$res2_succes[i] <- table(res_etape2[res_etape2=="succes"])
  }
  if (all(is.na(res_etape2[res_etape2=="echec"]))) {
    grid$res2_echecs[i] <- 0
  } else {
    grid$res2_echecs[i] <- table(res_etape2[res_etape2=="echec"])
  }
  
  #raison de l'arrêt à l'étape 2
  if (all(is.na(result2[result2=="tox:excessive ; rep:inadequate"]))) {
    grid$res2_toxE_repI[i] <- 0
  } else {
    grid$res2_toxE_repI[i] <- table(result2[result2=="tox:excessive ; rep:inadequate"])
  }
  if (all(is.na(result2[result2=="tox:excessive ; rep:OK"]))) {
    grid$res2_toxE_repO[i] <- 0
  } else {
    grid$res2_toxE_repO[i] <- table(result2[result2=="tox:excessive ; rep:OK"])
  }
  if (all(is.na(result2[result2=="tox:OK ; rep:inadequate"]))) {
    grid$res2_toxO_repI[i] <- 0
  } else {
    grid$res2_toxO_repI[i] <- table(result2[result2=="tox:OK ; rep:inadequate"])
  }
  #raison de l'arrêt à l'étape 1
  if (all(is.na(result1[result1=="tox:excessive ; rep:inadequate"]))) {
    grid$res1_toxE_repI[i] <- 0
  } else {
    grid$res1_toxE_repI[i] <- table(result1[result1=="tox:excessive ; rep:inadequate"])
  }
  if (all(is.na(result1[result1=="tox:excessive ; rep:OK"]))) {
    grid$res1_toxE_repO[i] <- 0
  } else {
    grid$res1_toxE_repO[i] <- table(result1[result1=="tox:excessive ; rep:OK"])
  }
  if (all(is.na(result1[result1=="tox:OK ; rep:inadequate"]))) {
    grid$res1_toxO_repI[i] <- 0
  } else {
    grid$res1_toxO_repI[i] <- table(result1[result1=="tox:OK ; rep:inadequate"])
  }
  
  #% arrêt pour tox
  grid$stop_tox[i] <- grid$res1_toxE_repI[i]+grid$res1_toxE_repO[i]+grid$res2_toxE_repI[i]+grid$res2_toxE_repO[i]
  
  #vérification s'il y a pb lors des AI (patients inclus en même temps)
  if (all(is.na(check_stop1[check_stop1=="STOP"]))) {
    grid$check_AI1[i] <- 0
  } else {
    grid$check_AI1[i] <- table(check_stop1[check_stop1=="STOP"])
  }
  
  
}

writexl::write_xlsx(grid, "C:\\util\\marie\\2022-2023\\stage\\R\\résultats\\Bryant&Day_tpspoisAI.xlsx")




####### Bryant & Day avant modif pour AI ######
#N1=28 ; CR1=7 ; CT1=20 ; N2=60 ; CR2=18 ; CT2=46 #seuils trouvés avec la fonction
N1=28 ; CR1=6 ; CT1=19 ; N2=60 ; CR2=17 ; CT2=45 #seuils -1

ntrials=10000
npat=N2

res_etape1=res_etape2=rep(NA,ntrials)
result1=result2=rep(NA,ntrials)
nb_eff=nb_ech=nb_tox=nb_notox=0
sample1=sample2=0
time=rep(NA,ntrials)
t_obs=4

param_eff=c(0.2,0.4,0.4,0.2,0.4,0.2,0.4,0.3)
param_tox=c(0.35,0.15,0.05,0.05,0.45,0.15,0.35,0.2)
scenario=paste("Scenario",1:8)
grid=cbind(data.frame(scenario,param_eff,param_tox))
colnames(grid)=c("Scenario","Peff","Ptox")


for (i in 1:nrow(grid)) {
  
  print(paste("Scenario",i))
  
  peff<-grid$Peff[i]
  ptox<-grid$Ptox[i]
  print(paste("peff=",peff,"ptox=",ptox))
  
  pefftox=peff*ptox
  peffnotox=peff-pefftox
  pnoefftox=ptox-pefftox
  pnoeffnotox=(1-peff)-pnoefftox
  
  res_etape1=res_etape2=rep(NA,ntrials)
  result1=result2=rep(NA,ntrials)
  nb_eff=nb_ech=nb_tox=nb_notox=0
  sample1=sample2=0
  time=rep(NA,ntrials)
  check_stop1=rep(NA,ntrials)
  
  
  set.seed(3003)
  for (n in 1:ntrials) {
    
    #print(paste("essai",n))
    #simulation des patients
    #patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))),cumsum(c(0,rep(2,npat-1)))) #temps fixe: cumsum(c(0,rep(2,npat-1))) temps aléatoire: cumsum(c(0,rpois(npat-1,2)))
    patient=t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox)))
    #fixe=cumsum(c(0,rep(2,npat-1)))
    pois=cumsum(c(0,rpois(npat-1,2)))
    patient=cbind(patient,pois) #patient=cbind(patient,pois)
    
    #ETAPE 1
    stage1=patient[1:N1,]
    eff1=sum(stage1[,1])+sum(stage1[,2])
    ech1=sum(stage1[,3])+sum(stage1[,4])
    tox1=sum(stage1[,2])+sum(stage1[,4])
    notox1=sum(stage1[,1])+sum(stage1[,3])
    
    if (eff1 > CR1 && notox1 > CT1) {res_etape1[n]<-"succes"} else {res_etape1[n] <- "echec"}
    sample1=sample1+N1
    nb_eff=nb_eff+eff1
    nb_ech=nb_ech+ech1
    nb_tox=nb_tox+tox1
    nb_notox=nb_notox+notox1
    
    #raison de l'arrêt précoce
    if (eff1<=CR1 && notox1<=CT1) {result1[n]<-"tox:excessive ; rep:inadequate"}
    if (eff1 >CR1 && notox1<=CT1) {result1[n]<-"tox:excessive ; rep:insuf data"}
    if (eff1<=CR1 && notox1 >CT1) {result1[n]<-"tox:insuf data ; rep:inadequate"}
    
    #durée de l'essai
    if (res_etape1[n]=="echec") {time[n]<-patient[N1,5]+t_obs}
    #print(paste(N1,"patients",eff1,"eff",notox1,"notox"))
    #print(paste("Résultat:",res_etape1[n],result1[n]))
    #print(paste("durée",time[n]))
    
    #vérification inclusion AI
    if (patient[N1,5]==patient[N1+1,5]) {check_stop1[n] <- "STOP"} else {"OK"}
    
    #ETAPE 2
    if (res_etape1[n]=="succes") {
      stage2=patient[(N1+1):N2,]
      eff2=sum(stage2[,1])+sum(stage2[,2])
      ech2=sum(stage2[,3])+sum(stage2[,4])
      tox2=sum(stage2[,2])+sum(stage2[,4])
      notox2=sum(stage2[,1])+sum(stage2[,3])
      
      if (eff1+eff2 > CR2 && notox1+notox2 > CT2) {res_etape2[n]<-"succes"} else {res_etape2[n] <- "echec"}
      sample1=sample1+N2-N1
      nb_eff=nb_eff+eff2
      nb_ech=nb_ech+ech2
      nb_tox=nb_tox+tox2
      nb_notox=nb_notox+notox2
      
      #raison de l'échec 
      if (eff1+eff2<=CR2 && notox1+notox2<=CT2) {result2[n]<-"tox:excessive ; rep:inadequate"}
      if (eff1+eff2 >CR2 && notox1+notox2<=CT2) {result2[n]<-"tox:excessive ; rep:OK"}
      if (eff1+eff2<=CR2 && notox1+notox2 >CT2) {result2[n]<-"tox:OK ; rep:inadequate"}
      
      #durée de l'essai
      if (sum(patient[1:N1-1,1])+sum(patient[1:N1-1,2])==CR1 | sum(patient[1:N1-1,1])+sum(patient[1:N1-1,3])==CT1-1) {time[n]<-patient[N2,5]+2*t_obs}
      if (sum(patient[1:N1-1,1])+sum(patient[1:N1-1,2])!=CR1 && sum(patient[1:N1-1,1])+sum(patient[1:N1-1,3])!=CT1-1) {time[n]<-patient[N2,5]+t_obs}
      #print(paste(N2-N1,"patients",eff2,"eff",notox2,"notox"))
      #print(paste("total eff:",eff1+eff2,"total notox:",notox1+notox2))
      #print(paste("Résultat:",res_etape2[n],result2[n]))
      #print(paste("durée",time[n]))
      
      
    }
    
  }
  
  #effectif attendu
  grid$ess[i] <- (sample1+sample2)/ntrials
  
  #nb d'efficacités 
  grid$moy_eff[i] <- nb_eff/ntrials
  
  #nb d'échecs
  grid$moy_ech[i] <- nb_ech/ntrials
  
  #nb de toxicités 
  grid$moy_tox[i] <- nb_tox/ntrials
  
  #durée de l'essai
  grid$moy_time[i] <- mean(time)
  
  #% conclu d'efficacité
  if (all(is.na(res_etape2[res_etape2=="succes"]))) {
    grid$succes[i] <- 0
  } else {
    grid$succes[i] <- ((table(res_etape2[res_etape2=="succes"]))*100)/ntrials
  }
  
  #% arrêt précoce
  if (all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$early_ter[i] <- 0
  } else {
    grid$early_ter[i] <- (table(res_etape1[res_etape1=="echec"])*100)/ntrials
  }
  
  #nb de succès et échecs à chaque étape
  if (all(is.na(res_etape1[res_etape1=="succes"]))) {
    grid$res1_succes[i] <- 0
  } else {
    grid$res1_succes[i] <- table(res_etape1[res_etape1=="succes"])
  }
  if (all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$res1_echecs[i] <- 0
  } else {
    grid$res1_echecs[i] <- table(res_etape1[res_etape1=="echec"])
  }
  
  if (all(is.na(res_etape2[res_etape2=="succes"]))) {
    grid$res2_succes[i] <- 0
  } else {
    grid$res2_succes[i] <- table(res_etape2[res_etape2=="succes"])
  }
  if (all(is.na(res_etape2[res_etape2=="echec"]))) {
    grid$res2_echecs[i] <- 0
  } else {
    grid$res2_echecs[i] <- table(res_etape2[res_etape2=="echec"])
  }
  
  #raison de l'arrêt à l'étape 2
  if (all(is.na(result2[result2=="tox:excessive ; rep:inadequate"]))) {
    grid$res2_toxE_repI[i] <- 0
  } else {
    grid$res2_toxE_repI[i] <- table(result2[result2=="tox:excessive ; rep:inadequate"])
  }
  if (all(is.na(result2[result2=="tox:excessive ; rep:OK"]))) {
    grid$res2_toxE_repO[i] <- 0
  } else {
    grid$res2_toxE_repO[i] <- table(result2[result2=="tox:excessive ; rep:OK"])
  }
  if (all(is.na(result2[result2=="tox:OK ; rep:inadequate"]))) {
    grid$res2_toxO_repI[i] <- 0
  } else {
    grid$res2_toxO_repI[i] <- table(result2[result2=="tox:OK ; rep:inadequate"])
  }
  #raison de l'arrêt à l'étape 1
  if (all(is.na(result1[result1=="tox:excessive ; rep:inadequate"]))) {
    grid$res1_toxE_repI[i] <- 0
  } else {
    grid$res1_toxE_repI[i] <- table(result1[result1=="tox:excessive ; rep:inadequate"])
  }
  if (all(is.na(result1[result1=="tox:excessive ; rep:OK"]))) {
    grid$res1_toxE_repO[i] <- 0
  } else {
    grid$res1_toxE_repO[i] <- table(result1[result1=="tox:excessive ; rep:OK"])
  }
  if (all(is.na(result1[result1=="tox:OK ; rep:inadequate"]))) {
    grid$res1_toxO_repI[i] <- 0
  } else {
    grid$res1_toxO_repI[i] <- table(result1[result1=="tox:OK ; rep:inadequate"])
  }
  
  #% arrêt pour tox
  grid$stop_tox[i] <- grid$res1_toxE_repI[i]+grid$res1_toxE_repO[i]+grid$res2_toxE_repI[i]+grid$res2_toxE_repO[i]
  
  #vérification s'il y a pb lors des AI (patients inclus en même temps)
  if (all(is.na(check_stop1[check_stop1=="STOP"]))) {
    grid$check_AI1[i] <- 0
  } else {
    grid$check_AI1[i] <- table(check_stop1[check_stop1=="STOP"])
  }
  
  
}

writexl::write_xlsx(grid, "C:\\util\\marie\\2022-2023\\stage\\R\\résultats\\Bryant&Day_tpspois2.xlsx")

table(res_etape1)
table(res_etape2)

table(result1)
table(result2)

#nb de patients
(sample1+sample2)/ntrials

#nb d'efficacités
nb_eff/ntrials

#nb d'échecs
nb_ech/ntrials

#nb de toxicités
nb_tox/ntrials

#durée de l'essai
mean(time)
