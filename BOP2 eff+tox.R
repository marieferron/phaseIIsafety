# BOP2 efficacy + toxicity

########### BOP2 eff + tox modif simulation patients AI ###########


#alpha=5%  
n1=6 ; r1=0 ; to1=3 ; n2=25 ; r2=4 ; to2=9 ; N=60 ; R=14 ; TO=19

ntrials=10000
npat=N

res_int1=res_int2=res_int3=rep(NA,ntrials) #résultat de l'essai à chaque étape (succès ou échec)
result1=result2=result3=rep(NA,ntrials) #résultat de toxicité à chaque étape 
sample1=sample2=sample3=0 #comptage du nombre de patients à chaque étape
time=rep(NA,ntrials)
t_obs=4 #mesure de l'outcome après 28j (4 semaines)


param_eff=c(0.2,0.2,0.4,0.4,0.4,0.2,0.4,0.3)
param_tox=c(0.35,0.15,0.35,0.15,0.05,0.05,0.45,0.25)
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
  
  res_int1=res_int2=res_int3=rep(NA,ntrials) #résultat de l'essai à chaque étape (succès ou échec)
  result1=result2=result3=rep(NA,ntrials) #résultat de toxicité à chaque étape 
  nb_eff=0 #nombre de réponses (succès) à chaque étape
  nb_ech=0 #nombre de non-réponses (échecs) à chaque étape
  nb_tox=0 #nombre de toxicités à chaque étape
  sample=sample1=sample2=sample3=0 #comptage du nombre de patients à chaque étape
  time=rep(NA,ntrials)
  p_eff=p_ech=p_tox=rep(NA,ntrials)
  check_stop1=check_stop2=rep(NA,ntrials)
  
  
  set.seed(3003)
  for (n in 1:ntrials) {
    #print(paste("essai",n))
    
    #simulation des patients et du moment d'inclusion
    #patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))), cumsum(c(0,rpois(N-1,2)))) #temps aléatoire: cumsum(c(0,rpois(N-1,2))) temps fixe: cumsum(c(0,rep(2,N-1)))
    patient=t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox)))
    pois=c(0,rpois(N-1,2))
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
    
    nbeff1=nbeff2=nbeff3=0
    nbech1=nbech2=nbech3=0
    nbtox1=nbtox2=nbtox3=0
    sample1=sample2=sample3=0
    
    
    #1ère analyse intermédiaire (6 patients)
    stage1=patient[1:n1,]
    eff1=sum(stage1[,1])+sum(stage1[,2])
    tox1=sum(stage1[,2])+sum(stage1[,4])
    ech1=sum(stage1[,3])+sum(stage1[,4])
    
    sample=sample+n1
    sample1=n1
    nb_eff=nb_eff+eff1
    nb_tox=nb_tox+tox1
    nb_ech=nb_ech+ech1
    nbeff1=eff1
    nbech1=ech1
    nbtox1=tox1
    
    if (eff1 > r1 && tox1 < to1) {res_int1[n] <- "succes"} else {res_int1[n] <- "echec"}
    #print(paste(eff1,"eff",tox1,"tox",ech1,"échecs"))
    #print(paste("résultat 1:", res_int1[n]))
    
    if (res_int1[n] == "echec") {time[n] <- patient[n1,5]+t_obs}
    
    if (eff1 > r1 && tox1 < to1) {result1[n]<-"tox:OK eff:OK"}
    if (eff1 > r1 && tox1 >= to1) {result1[n]<-"tox:excessive eff:OK"}
    if (eff1 <= r1 && tox1 < to1) {result1[n]<-"tox:OK eff:inadequate"}
    if (eff1 <= r1 && tox1 >= to1) {result1[n]<-"tox:excessive eff:inadequate"}
    
    #vérification inclusion AI
    if (patient[n1,5]==patient[n1+1,5]) {check_stop1[n] <- "STOP"} else {"OK"}
    if (patient[n2,5]==patient[n2+1,5]) {check_stop2[n] <- "STOP"} else {"OK"}
    
    #2ème analyse intermédiaire (19 patients)
    if (res_int1[n] == "succes") {
      stage2=patient[n1+1:(n2-n1),]
      eff2=sum(stage2[,1])+sum(stage2[,2])
      tox2=sum(stage2[,2])+sum(stage2[,4])
      ech2=sum(stage2[,3])+sum(stage2[,4])
      
      sample=sample+n2-n1
      sample2=n2-n1
      nb_eff=nb_eff+eff2
      nb_tox=nb_tox+tox2
      nb_ech=nb_ech+ech2
      nbeff2=eff2
      nbech2=ech2
      nbtox2=tox2
      
      if (eff1+eff2 > r2 && tox1+tox2 < to2) {res_int2[n] <- "succes"} else {res_int2[n] <- "echec"}
      #print(paste(eff2,"eff",tox2,"tox",ech2,"échecs"))
      #print(paste("résultat 2:", res_int2[n]))
      
      if (res_int2[n] == "echec") {
        if (sum(patient[1:5,1])+sum(patient[1:5,2])==r1 | sum(patient[1:5,2])+sum(patient[1:5,4])==to1-1) {
          if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,1])+sum(patient[7:24,2]) != to2-1) {time[n] <- patient[n2,5]+t_obs} #ici t_obs pour l'attente lors de la première analyse intermédiaire
          if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,1])+sum(patient[7:24,2]) == to2-1) {time[n] <- patient[n2,5]+2*t_obs}
        }
        if (sum(patient[1:5,1])+sum(patient[1:5,2])!=r1 && sum(patient[1:5,2])+sum(patient[1:5,4])!=to1-1) {
          if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,1])+sum(patient[7:24,2]) != to2-1) {time[n] <- patient[n2,5]} 
          if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,1])+sum(patient[7:24,2]) == to2-1) {time[n] <- patient[n2,5]+t_obs}
        }
        #print(paste("Durée (2):",time[n]))
      }
      
      if (eff1+eff2 > r2 && tox1+tox2 < to2) {result2[n]<-"tox:OK eff:OK"}
      if (eff1+eff2 > r2 && tox1+tox2 >= to2) {result2[n]<-"tox:excessive eff:OK"}
      if (eff1+eff2 <= r2 && tox1+tox2 < to2) {result2[n]<-"tox:OK eff:inadequate"}
      if (eff1+eff2 <= r2 && tox1+tox2 >= to2) {result2[n]<-"tox:excessive eff:inadequate"}
      
      
      #3ème analyse intermédiaire (35 patients)
      if (res_int2[n] == "succes") { 
        eff3=sum(patient[n2+1:(N-n2),1])+sum(patient[n2+1:(N-n2),2])
        tox3=sum(patient[n2+1:(N-n2),2])+sum(patient[n2+1:(N-n2),4])
        ech3=sum(patient[n2+1:(N-n2),3])+sum(patient[n2+1:(N-n2),4])
        
        sample=sample+N-n2
        sample3=N-n2
        nb_eff=nb_eff+eff3
        nb_tox=nb_tox+tox3
        nb_ech=nb_ech+ech3
        nbeff3=eff3
        nbech3=ech3
        nbtox3=tox3
        
        if (eff1+eff2+eff3 > R && tox1+tox2+tox3 < TO) {res_int3[n] <- "succes"} else {res_int3[n] <- "echec"}
        #print(paste(eff3,"eff",tox3,"tox",ech3,"échecs"))
        #print(paste("résultat 3:", res_int3[n]))
        
        
        if (sum(patient[1:5,1])+sum(patient[1:5,2])==r1 | sum(patient[1:5,2])+sum(patient[1:5,4])==to1-1) {
          if (eff1+eff2+sum(patient[26:59,1])+sum(patient[26:59,2]) != R && tox1+tox2+sum(patient[26:59,2])+sum(patient[26:59,4]) != TO-1) {
            if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) == to2-1) {time[n] <- patient[N,5]+3*t_obs}
            if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) != to2-1) {time[n] <- patient[N,5]+2*t_obs}
          }
          if (eff1+eff2+sum(patient[26:59,1])+sum(patient[26:59,2]) == R | tox1+tox2+sum(patient[26:59,2])+sum(patient[26:59,4]) == TO-1) {
            if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) == to2-1) {time[n] <- patient[N,5]+3*t_obs}
            if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) != to2-1) {time[n] <- patient[N,5]+2*t_obs}
          }
          
        }
        if (sum(patient[1:5,1])+sum(patient[1:5,2])!=r1 && sum(patient[1:5,2])+sum(patient[1:5,4])!=to1-1) {
          if (eff1+eff2+sum(patient[26:59,1])+sum(patient[26:59,2]) != R | tox1+tox2+sum(patient[26:59,2])+sum(patient[26:59,4]) != TO-1) {
            if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) == to2-1) {time[n] <- patient[N,5]+2*t_obs}
            if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) != to2-1) {time[n] <- patient[N,5]+t_obs}
          }
          if (eff1+eff2+sum(patient[26:59,1])+sum(patient[26:59,2]) == R | tox1+tox2+sum(patient[26:59,2])+sum(patient[26:59,4]) == TO-1) {
            if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) == to2-1) {time[n] <- patient[N,5]+2*t_obs}
            if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) != to2-1) {time[n] <- patient[N,5]+t_obs}
          }
        }
        #print(paste("Durée (3):",time[n]))
        
        if (eff1+eff2+eff3 > R && tox1+tox2+tox3 < TO) {result3[n]<-"tox:OK eff:OK"}
        if (eff1+eff2+eff3 > R && tox1+tox2+tox3 >= TO) {result3[n]<-"tox:excessive eff:OK"}
        if (eff1+eff2+eff3 <= R && tox1+tox2+tox3 < TO) {result3[n]<-"tox:OK eff:inadequate"}
        if (eff1+eff2+eff3 <= R && tox1+tox2+tox3 >= TO) {result3[n]<-"tox:excessive eff:inadequate"}
        
      }
      
    }
    
    p_eff[n] <- (nbeff1+nbeff2+nbeff3) / (sample1+sample2+sample3) 
    p_ech[n] <- (nbech1+nbech2+nbech3) / (sample1+sample2+sample3) 
    p_tox[n] <- (nbtox1+nbtox2+nbtox3) / (sample1+sample2+sample3)
  
  }
  
  
  #effectif attendu
  grid$ess[i] <- sample/ntrials
  
  #nb d'efficacités 
  grid$moy_eff[i] <- nb_eff/ntrials
  grid$prob_eff[i] <- nb_eff/sample
  grid$prob_eff2[i] <- mean(p_eff)
  
  #nb d'échecs
  grid$moy_ech[i] <- nb_ech/ntrials
  grid$prob_ech[i] <- nb_ech/sample
  grid$prob_ech2[i] <- mean(p_ech)
  
  #nb de toxicités 
  grid$moy_tox[i] <- nb_tox/ntrials
  grid$prob_tox[i] <- nb_tox/sample
  grid$prob_tox2[i] <- mean(p_tox)
  
  #durée de l'essai
  grid$moy_time[i] <- mean(time)
  
  
  #% conclu d'efficacité
  if (all(is.na(res_int3[res_int3=="succes"]))) {
    grid$succes[i] <- 0
  } else {
    grid$succes[i] <- ((table(res_int3[res_int3=="succes"]))*100)/ntrials
  }
  
  #% arrêt précoce
  if (all(is.na(res_int1[res_int1=="echec"])) && all(is.na(res_int2[res_int2=="echec"]))) {
    grid$early_ter[i] <- 0
  } else {
    grid$early_ter[i] <- ((table(res_int1[res_int1=="echec"])+table(res_int2[res_int2=="echec"]))*100)/ntrials
  }
  
  #nb de succès et échecs à chaque étape
  if (all(is.na(res_int1[res_int1=="succes"]))) {
    grid$res1_succes[i] <- 0
  } else {
    grid$res1_succes[i] <- table(res_int1[res_int1=="succes"])
  }
  if (all(is.na(res_int1[res_int1=="echec"]))) {
    grid$res1_echecs[i] <- 0
  } else {
    grid$res1_echecs[i] <- table(res_int1[res_int1=="echec"])
  }
  
  if (all(is.na(res_int2[res_int2=="succes"]))) {
    grid$res2_succes[i] <- 0
  } else {
    grid$res2_succes[i] <- table(res_int2[res_int2=="succes"])
  }
  if (all(is.na(res_int2[res_int2=="echec"]))) {
    grid$res2_echecs[i] <- 0
  } else {
    grid$res2_echecs[i] <- table(res_int2[res_int2=="echec"])
  }
  
  if (all(is.na(res_int3[res_int3=="succes"]))) {
    grid$res3_succes[i] <- 0
  } else {
    grid$res3_succes[i] <- table(res_int3[res_int3=="succes"])
  }
  if (all(is.na(res_int3[res_int3=="echec"]))) {
    grid$res3_echecs[i] <- 0
  } else {
    grid$res3_echecs[i] <- table(res_int3[res_int3=="echec"])
  }
  
  #raison de l'arrêt à l'étape 1
  if (all(is.na(result1[result1=="tox:excessive eff:inadequate"]))) {
    grid$res1_toxE_repI[i] <- 0
  } else {
    grid$res1_toxE_repI[i] <- table(result1[result1=="tox:excessive eff:inadequate"])
  }
  if (all(is.na(result1[result1=="tox:excessive eff:OK"]))) {
    grid$res1_toxE_repO[i] <- 0
  } else {
    grid$res1_toxE_repO[i] <- table(result1[result1=="tox:excessive eff:OK"])
  }
  if (all(is.na(result1[result1=="tox:OK eff:inadequate"]))) {
    grid$res1_toxO_repI[i] <- 0
  } else {
    grid$res1_toxO_repI[i] <- table(result1[result1=="tox:OK eff:inadequate"])
  }
  #raison de l'arrêt à l'étape 2
  if (all(is.na(result2[result2=="tox:excessive eff:inadequate"]))) {
    grid$res2_toxE_repI[i] <- 0
  } else {
    grid$res2_toxE_repI[i] <- table(result2[result2=="tox:excessive eff:inadequate"])
  }
  if (all(is.na(result2[result2=="tox:excessive eff:OK"]))) {
    grid$res2_toxE_repO[i] <- 0
  } else {
    grid$res2_toxE_repO[i] <- table(result2[result2=="tox:excessive eff:OK"])
  }
  if (all(is.na(result2[result2=="tox:OK eff:inadequate"]))) {
    grid$res2_toxO_repI[i] <- 0
  } else {
    grid$res2_toxO_repI[i] <- table(result2[result2=="tox:OK eff:inadequate"])
  }
  #raison de l'arrêt à l'étape 2
  if (all(is.na(result3[result3=="tox:excessive eff:inadequate"]))) {
    grid$res3_toxE_repI[i] <- 0
  } else {
    grid$res3_toxE_repI[i] <- table(result3[result3=="tox:excessive eff:inadequate"])
  }
  if (all(is.na(result3[result3=="tox:excessive eff:OK"]))) {
    grid$res3_toxE_repO[i] <- 0
  } else {
    grid$res3_toxE_repO[i] <- table(result3[result3=="tox:excessive eff:OK"])
  }
  if (all(is.na(result3[result3=="tox:OK eff:inadequate"]))) {
    grid$res3_toxO_repI[i] <- 0
  } else {
    grid$res3_toxO_repI[i] <- table(result2[result2=="tox:OK eff:inadequate"])
  }
  
  
  #% arrêt pour tox
  grid$stop_tox[i] <- grid$res1_toxE_repI[i]+grid$res1_toxE_repO[i]+grid$res2_toxE_repI[i]+grid$res2_toxE_repO[i]+grid$res3_toxE_repI[i]+grid$res3_toxE_repO[i]
  
  #vérification s'il y a pb lors des AI (patients inclus en même temps)
  if (all(is.na(check_stop1[check_stop1=="STOP"]))) {
    grid$check_AI1[i] <- 0
  } else {
    grid$check_AI1[i] <- table(check_stop1[check_stop1=="STOP"])
  }
  if (all(is.na(check_stop2[check_stop2=="STOP"]))) {
    grid$check_AI2[i] <- 0
  } else {
    grid$check_AI2[i] <- table(check_stop2[check_stop2=="STOP"])
  }
  
}

writexl::write_xlsx(grid, "C:\\util\\marie\\2022-2023\\stage\\R\\résultats\\BOP2efftox_tpspoisAI+prob.xlsx")





##### test #### 
n1=6 ; r1=0 ; to1=3 ; n2=25 ; r2=4 ; to2=9 ; N=60 ; R=14 ; TO=19

ntrials=10000
npat=N

res_int1=res_int2=res_int3=rep(NA,ntrials) #résultat de l'essai à chaque étape (succès ou échec)
result1=result2=result3=rep(NA,ntrials) #résultat de toxicité à chaque étape 
sample1=sample2=sample3=0 #comptage du nombre de patients à chaque étape
time=rep(NA,ntrials)
t_obs=4 #mesure de l'outcome après 28j (4 semaines)

peff=0.2
ptox=0.35
pefftox=peff*ptox
peffnotox=peff-pefftox
pnoefftox=ptox-pefftox
pnoeffnotox=(1-peff)-pnoefftox

res_int1=res_int2=res_int3=rep(NA,ntrials) #résultat de l'essai à chaque étape (succès ou échec)
result1=result2=result3=rep(NA,ntrials) #résultat de toxicité à chaque étape 
nb_eff=0 #nombre de réponses (succès) à chaque étape
nb_ech=0 #nombre de non-réponses (échecs) à chaque étape
nb_tox=0 #nombre de toxicités à chaque étape
sample=sample1=sample2=sample3=0 #comptage du nombre de patients à chaque étape
time=rep(NA,ntrials)
check_stop1=check_stop2=rep(NA,ntrials)
p_eff=p_ech=p_tox=rep(NA,ntrials)


set.seed(3003)
for (n in 1:ntrials) {
  #print(paste("essai",n))
  
  #simulation des patients et du moment d'inclusion
  #patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))), cumsum(c(0,rpois(N-1,2)))) #temps aléatoire: cumsum(c(0,rpois(N-1,2))) temps fixe: cumsum(c(0,rep(2,N-1)))
  patient=t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox)))
  pois=c(0,rpois(N-1,2))
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
  
  nbeff1=nbeff2=nbeff3=0
  nbech1=nbech2=nbech3=0
  nbtox1=nbtox2=nbtox3=0
  sample1=sample2=sample3=0
  
  
  #1ère analyse intermédiaire (6 patients)
  stage1=patient[1:n1,]
  eff1=sum(stage1[,1])+sum(stage1[,2])
  tox1=sum(stage1[,2])+sum(stage1[,4])
  ech1=sum(stage1[,3])+sum(stage1[,4])
  
  sample=sample+n1
  sample1=n1
  nb_eff=nb_eff+eff1
  nb_tox=nb_tox+tox1
  nb_ech=nb_ech+ech1
  nbeff1=eff1
  nbech1=ech1
  nbtox1=tox1
  
  if (eff1 > r1 && tox1 < to1) {res_int1[n] <- "succes"} else {res_int1[n] <- "echec"}
  #print(paste(eff1,"eff",tox1,"tox",ech1,"échecs"))
  #print(paste("résultat 1:", res_int1[n]))
  
  if (res_int1[n] == "echec") {time[n] <- patient[n1,5]+t_obs}
  
  if (eff1 > r1 && tox1 < to1) {result1[n]<-"tox:OK eff:OK"}
  if (eff1 > r1 && tox1 >= to1) {result1[n]<-"tox:excessive eff:OK"}
  if (eff1 <= r1 && tox1 < to1) {result1[n]<-"tox:OK eff:inadequate"}
  if (eff1 <= r1 && tox1 >= to1) {result1[n]<-"tox:excessive eff:inadequate"}
  
  #vérification inclusion AI
  if (patient[n1,5]==patient[n1+1,5]) {check_stop1[n] <- "STOP"} else {"OK"}
  if (patient[n2,5]==patient[n2+1,5]) {check_stop2[n] <- "STOP"} else {"OK"}
  
  #2ème analyse intermédiaire (19 patients)
  if (res_int1[n] == "succes") {
    stage2=patient[n1+1:(n2-n1),]
    eff2=sum(stage2[,1])+sum(stage2[,2])
    tox2=sum(stage2[,2])+sum(stage2[,4])
    ech2=sum(stage2[,3])+sum(stage2[,4])
    
    sample=sample+n2-n1
    sample2=n2-n1
    nb_eff=nb_eff+eff2
    nb_tox=nb_tox+tox2
    nb_ech=nb_ech+ech2
    nbeff2=eff2
    nbech2=ech2
    nbtox2=tox2
    
    if (eff1+eff2 > r2 && tox1+tox2 < to2) {res_int2[n] <- "succes"} else {res_int2[n] <- "echec"}
    #print(paste(eff2,"eff",tox2,"tox",ech2,"échecs"))
    #print(paste("résultat 2:", res_int2[n]))
    
    if (res_int2[n] == "echec") {
      if (sum(patient[1:5,1])+sum(patient[1:5,2])==r1 | sum(patient[1:5,2])+sum(patient[1:5,4])==to1-1) {
        if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,1])+sum(patient[7:24,2]) != to2-1) {time[n] <- patient[n2,5]+t_obs} #ici t_obs pour l'attente lors de la première analyse intermédiaire
        if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,1])+sum(patient[7:24,2]) == to2-1) {time[n] <- patient[n2,5]+2*t_obs}
      }
      if (sum(patient[1:5,1])+sum(patient[1:5,2])!=r1 && sum(patient[1:5,2])+sum(patient[1:5,4])!=to1-1) {
        if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,1])+sum(patient[7:24,2]) != to2-1) {time[n] <- patient[n2,5]} 
        if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,1])+sum(patient[7:24,2]) == to2-1) {time[n] <- patient[n2,5]+t_obs}
      }
      #print(paste("Durée (2):",time[n]))
    }
    
    if (eff1+eff2 > r2 && tox1+tox2 < to2) {result2[n]<-"tox:OK eff:OK"}
    if (eff1+eff2 > r2 && tox1+tox2 >= to2) {result2[n]<-"tox:excessive eff:OK"}
    if (eff1+eff2 <= r2 && tox1+tox2 < to2) {result2[n]<-"tox:OK eff:inadequate"}
    if (eff1+eff2 <= r2 && tox1+tox2 >= to2) {result2[n]<-"tox:excessive eff:inadequate"}
    
    
    #3ème analyse intermédiaire (35 patients)
    if (res_int2[n] == "succes") { 
      eff3=sum(patient[n2+1:(N-n2),1])+sum(patient[n2+1:(N-n2),2])
      tox3=sum(patient[n2+1:(N-n2),2])+sum(patient[n2+1:(N-n2),4])
      ech3=sum(patient[n2+1:(N-n2),3])+sum(patient[n2+1:(N-n2),4])
      
      sample=sample+N-n2
      sample3=N-n2
      nb_eff=nb_eff+eff3
      nb_tox=nb_tox+tox3
      nb_ech=nb_ech+ech3
      nbeff3=eff3
      nbech3=ech3
      nbtox3=tox3
      
      if (eff1+eff2+eff3 > R && tox1+tox2+tox3 < TO) {res_int3[n] <- "succes"} else {res_int3[n] <- "echec"}
      #print(paste(eff3,"eff",tox3,"tox",ech3,"échecs"))
      #print(paste("résultat 3:", res_int3[n]))
      
      
      if (sum(patient[1:5,1])+sum(patient[1:5,2])==r1 | sum(patient[1:5,2])+sum(patient[1:5,4])==to1-1) {
        if (eff1+eff2+sum(patient[26:59,1])+sum(patient[26:59,2]) != R && tox1+tox2+sum(patient[26:59,2])+sum(patient[26:59,4]) != TO-1) {
          if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) == to2-1) {time[n] <- patient[N,5]+3*t_obs}
          if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) != to2-1) {time[n] <- patient[N,5]+2*t_obs}
        }
        if (eff1+eff2+sum(patient[26:59,1])+sum(patient[26:59,2]) == R | tox1+tox2+sum(patient[26:59,2])+sum(patient[26:59,4]) == TO-1) {
          if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) == to2-1) {time[n] <- patient[N,5]+3*t_obs}
          if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) != to2-1) {time[n] <- patient[N,5]+2*t_obs}
        }
        
      }
      if (sum(patient[1:5,1])+sum(patient[1:5,2])!=r1 && sum(patient[1:5,2])+sum(patient[1:5,4])!=to1-1) {
        if (eff1+eff2+sum(patient[26:59,1])+sum(patient[26:59,2]) != R | tox1+tox2+sum(patient[26:59,2])+sum(patient[26:59,4]) != TO-1) {
          if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) == to2-1) {time[n] <- patient[N,5]+2*t_obs}
          if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) != to2-1) {time[n] <- patient[N,5]+t_obs}
        }
        if (eff1+eff2+sum(patient[26:59,1])+sum(patient[26:59,2]) == R | tox1+tox2+sum(patient[26:59,2])+sum(patient[26:59,4]) == TO-1) {
          if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) == to2-1) {time[n] <- patient[N,5]+2*t_obs}
          if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) != to2-1) {time[n] <- patient[N,5]+t_obs}
        }
      }
      #print(paste("Durée (3):",time[n]))
      
      if (eff1+eff2+eff3 > R && tox1+tox2+tox3 < TO) {result3[n]<-"tox:OK eff:OK"}
      if (eff1+eff2+eff3 > R && tox1+tox2+tox3 >= TO) {result3[n]<-"tox:excessive eff:OK"}
      if (eff1+eff2+eff3 <= R && tox1+tox2+tox3 < TO) {result3[n]<-"tox:OK eff:inadequate"}
      if (eff1+eff2+eff3 <= R && tox1+tox2+tox3 >= TO) {result3[n]<-"tox:excessive eff:inadequate"}
      
    }
    
  }
  
  p_eff[n] <- (nbeff1+nbeff2+nbeff3) / (sample1+sample2+sample3) 
  p_ech[n] <- (nbech1+nbech2+nbech3) / (sample1+sample2+sample3) 
  p_tox[n] <- (nbtox1+nbtox2+nbtox3) / (sample1+sample2+sample3)
  
}

#100 simulations (peff=0.4 et ptox=0.15)
mean(p_eff) #=0.3935
nb_eff/sample #=0.3952709
#100 simulations (peff=0.2 et ptox=0.35)
mean(p_eff) #=0.1664333
nb_eff/sample #=0.1898907

#10000 simulations (peff=0.4 et ptox=0.15)
mean(p_eff) #=0.4009767
nb_eff/sample #=0.4007756
#10000 simulations (peff=0.2 et ptox=0.35)
mean(p_eff) #=0.1675233
nb_eff/sample #=0.2012923

########### BOP2 eff + tox avant modif pour AI ###########
#code pour tous les scenarios

#POUR PTOX1=0.2
#alpha=5% et ptox1=0,2
#n1=6 ; r1=0 ; to1=3 ; n2=25 ; r2=4 ; to2=9 ; N=60 ; R=14 ; TO=19

#POUR PTOX1=0.15  
#alpha=5%  
n1=6 ; r1=0 ; to1=3 ; n2=25 ; r2=4 ; to2=9 ; N=60 ; R=14 ; TO=19

ntrials=10000
npat=N

res_int1=res_int2=res_int3=rep(NA,ntrials) #résultat de l'essai à chaque étape (succès ou échec)
result1=result2=result3=rep(NA,ntrials) #résultat de toxicité à chaque étape 
sample1=sample2=sample3=0 #comptage du nombre de patients à chaque étape
time=rep(NA,ntrials)
t_obs=4 #mesure de l'outcome après 28j (4 semaines)


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
  
  res_int1=res_int2=res_int3=rep(NA,ntrials) #résultat de l'essai à chaque étape (succès ou échec)
  result1=result2=result3=rep(NA,ntrials) #résultat de toxicité à chaque étape 
  nb_eff=0 #nombre de réponses (succès) à chaque étape
  nb_ech=0 #nombre de non-réponses (échecs) à chaque étape
  nb_tox=0 #nombre de toxicités à chaque étape
  sample1=sample2=sample3=0 #comptage du nombre de patients à chaque étape
  time=rep(NA,ntrials)
  check_stop1=check_stop2=rep(NA,ntrials)
  
  
  set.seed(3003)
  for (n in 1:ntrials) {
    #print(paste("essai",n))
    
    #simulation des patients et du moment d'inclusion
    patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))), cumsum(c(0,rpois(N-1,2)))) #temps aléatoire: cumsum(c(0,rpois(N-1,2))) temps fixe: cumsum(c(0,rep(2,N-1)))
    #patient=t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox)))
    #fixe=cumsum(c(0,rep(2,npat-1)))
    #pois=cumsum(c(0,rpois(npat-1,2)))
    #patient=cbind(patient,fixe) #patient=cbind(patient,pois)
    
    
    #1ère analyse intermédiaire (6 patients)
    stage1=patient[1:n1,]
    eff1=sum(stage1[,1])+sum(stage1[,2])
    tox1=sum(stage1[,2])+sum(stage1[,4])
    ech1=sum(stage1[,3])+sum(stage1[,4])
    
    sample1=sample1+n1
    nb_eff=nb_eff+eff1
    nb_tox=nb_tox+tox1
    nb_ech=nb_ech+ech1
    
    if (eff1 > r1 && tox1 < to1) {res_int1[n] <- "succes"} else {res_int1[n] <- "echec"}
    #print(paste(eff1,"eff",tox1,"tox",ech1,"échecs"))
    #print(paste("résultat 1:", res_int1[n]))
    
    if (res_int1[n] == "echec") {time[n] <- patient[n1,5]+t_obs}
    
    if (eff1 > r1 && tox1 < to1) {result1[n]<-"tox:OK eff:OK"}
    if (eff1 > r1 && tox1 >= to1) {result1[n]<-"tox:excessive eff:OK"}
    if (eff1 <= r1 && tox1 < to1) {result1[n]<-"tox:OK eff:inadequate"}
    if (eff1 <= r1 && tox1 >= to1) {result1[n]<-"tox:excessive eff:inadequate"}
    
    #vérification inclusion AI
    if (patient[n1,5]==patient[n1+1,5]) {check_stop1[n] <- "STOP"} else {"OK"}
    if (patient[n2,5]==patient[n2+1,5]) {check_stop2[n] <- "STOP"} else {"OK"}
    
    #2ème analyse intermédiaire (19 patients)
    if (res_int1[n] == "succes") {
      stage2=patient[n1+1:(n2-n1),]
      eff2=sum(stage2[,1])+sum(stage2[,2])
      tox2=sum(stage2[,2])+sum(stage2[,4])
      ech2=sum(stage2[,3])+sum(stage2[,4])
      
      sample2=sample2+n2-n1
      nb_eff=nb_eff+eff2
      nb_tox=nb_tox+tox2
      nb_ech=nb_ech+ech2
      
      if (eff1+eff2 > r2 && tox1+tox2 < to2) {res_int2[n] <- "succes"} else {res_int2[n] <- "echec"}
      #print(paste(eff2,"eff",tox2,"tox",ech2,"échecs"))
      #print(paste("résultat 2:", res_int2[n]))
      
      if (res_int2[n] == "echec") {
        if (sum(patient[1:5,1])+sum(patient[1:5,2])==r1 | sum(patient[1:5,2])+sum(patient[1:5,4])==to1-1) {
          if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,1])+sum(patient[7:24,2]) != to2-1) {time[n] <- patient[n2,5]+t_obs} #ici t_obs pour l'attente lors de la première analyse intermédiaire
          if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,1])+sum(patient[7:24,2]) == to2-1) {time[n] <- patient[n2,5]+2*t_obs}
        }
        if (sum(patient[1:5,1])+sum(patient[1:5,2])!=r1 && sum(patient[1:5,2])+sum(patient[1:5,4])!=to1-1) {
          if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,1])+sum(patient[7:24,2]) != to2-1) {time[n] <- patient[n2,5]} 
          if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,1])+sum(patient[7:24,2]) == to2-1) {time[n] <- patient[n2,5]+t_obs}
        }
        #print(paste("Durée (2):",time[n]))
      }
      
      if (eff1+eff2 > r2 && tox1+tox2 < to2) {result2[n]<-"tox:OK eff:OK"}
      if (eff1+eff2 > r2 && tox1+tox2 >= to2) {result2[n]<-"tox:excessive eff:OK"}
      if (eff1+eff2 <= r2 && tox1+tox2 < to2) {result2[n]<-"tox:OK eff:inadequate"}
      if (eff1+eff2 <= r2 && tox1+tox2 >= to2) {result2[n]<-"tox:excessive eff:inadequate"}
      
      
      #3ème analyse intermédiaire (35 patients)
      if (res_int2[n] == "succes") { 
        eff3=sum(patient[n2+1:(N-n2),1])+sum(patient[n2+1:(N-n2),2])
        tox3=sum(patient[n2+1:(N-n2),2])+sum(patient[n2+1:(N-n2),4])
        ech3=sum(patient[n2+1:(N-n2),3])+sum(patient[n2+1:(N-n2),4])
        
        sample3=sample3+N-n2
        nb_eff=nb_eff+eff3
        nb_tox=nb_tox+tox3
        nb_ech=nb_ech+ech3
        
        if (eff1+eff2+eff3 > R && tox1+tox2+tox3 < TO) {res_int3[n] <- "succes"} else {res_int3[n] <- "echec"}
        #print(paste(eff3,"eff",tox3,"tox",ech3,"échecs"))
        #print(paste("résultat 3:", res_int3[n]))
        
        
        if (sum(patient[1:5,1])+sum(patient[1:5,2])==r1 | sum(patient[1:5,2])+sum(patient[1:5,4])==to1-1) {
          if (eff1+eff2+sum(patient[26:59,1])+sum(patient[26:59,2]) != R && tox1+tox2+sum(patient[26:59,2])+sum(patient[26:59,4]) != TO-1) {
            if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) == to2-1) {time[n] <- patient[N,5]+3*t_obs}
            if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) != to2-1) {time[n] <- patient[N,5]+2*t_obs}
          }
          if (eff1+eff2+sum(patient[26:59,1])+sum(patient[26:59,2]) == R | tox1+tox2+sum(patient[26:59,2])+sum(patient[26:59,4]) == TO-1) {
            if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) == to2-1) {time[n] <- patient[N,5]+3*t_obs}
            if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) != to2-1) {time[n] <- patient[N,5]+2*t_obs}
          }
         
        }
        if (sum(patient[1:5,1])+sum(patient[1:5,2])!=r1 && sum(patient[1:5,2])+sum(patient[1:5,4])!=to1-1) {
          if (eff1+eff2+sum(patient[26:59,1])+sum(patient[26:59,2]) != R | tox1+tox2+sum(patient[26:59,2])+sum(patient[26:59,4]) != TO-1) {
            if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) == to2-1) {time[n] <- patient[N,5]+2*t_obs}
            if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) != to2-1) {time[n] <- patient[N,5]+t_obs}
          }
          if (eff1+eff2+sum(patient[26:59,1])+sum(patient[26:59,2]) == R | tox1+tox2+sum(patient[26:59,2])+sum(patient[26:59,4]) == TO-1) {
            if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) == r2 | tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) == to2-1) {time[n] <- patient[N,5]+2*t_obs}
            if (eff1+sum(patient[7:24,1])+sum(patient[7:24,2]) != r2 && tox1+sum(patient[7:24,2])+sum(patient[7:24,4]) != to2-1) {time[n] <- patient[N,5]+t_obs}
          }
        }
        #print(paste("Durée (3):",time[n]))
        
        if (eff1+eff2+eff3 > R && tox1+tox2+tox3 < TO) {result3[n]<-"tox:OK eff:OK"}
        if (eff1+eff2+eff3 > R && tox1+tox2+tox3 >= TO) {result3[n]<-"tox:excessive eff:OK"}
        if (eff1+eff2+eff3 <= R && tox1+tox2+tox3 < TO) {result3[n]<-"tox:OK eff:inadequate"}
        if (eff1+eff2+eff3 <= R && tox1+tox2+tox3 >= TO) {result3[n]<-"tox:excessive eff:inadequate"}
        
      }
      
    }
    
  }
  
  
  #effectif attendu
  grid$ess[i] <- (sample1+sample2+sample3)/ntrials
  
  #nb d'efficacités 
  grid$moy_eff[i] <- nb_eff/ntrials
  
  #nb d'échecs
  grid$moy_ech[i] <- nb_ech/ntrials
  
  #nb de toxicités 
  grid$moy_tox[i] <- nb_tox/ntrials
  
  #durée de l'essai
  grid$moy_time[i] <- mean(time)
  
  
  #% conclu d'efficacité
  if (all(is.na(res_int3[res_int3=="succes"]))) {
    grid$succes[i] <- 0
  } else {
    grid$succes[i] <- ((table(res_int3[res_int3=="succes"]))*100)/ntrials
  }
  
  #% arrêt précoce
  if (all(is.na(res_int1[res_int1=="echec"])) && all(is.na(res_int2[res_int2=="echec"]))) {
    grid$early_ter[i] <- 0
  } else {
    grid$early_ter[i] <- ((table(res_int1[res_int1=="echec"])+table(res_int2[res_int2=="echec"]))*100)/ntrials
  }
  
  #nb de succès et échecs à chaque étape
  if (all(is.na(res_int1[res_int1=="succes"]))) {
    grid$res1_succes[i] <- 0
  } else {
    grid$res1_succes[i] <- table(res_int1[res_int1=="succes"])
  }
  if (all(is.na(res_int1[res_int1=="echec"]))) {
    grid$res1_echecs[i] <- 0
  } else {
    grid$res1_echecs[i] <- table(res_int1[res_int1=="echec"])
  }
  
  if (all(is.na(res_int2[res_int2=="succes"]))) {
    grid$res2_succes[i] <- 0
  } else {
    grid$res2_succes[i] <- table(res_int2[res_int2=="succes"])
  }
  if (all(is.na(res_int2[res_int2=="echec"]))) {
    grid$res2_echecs[i] <- 0
  } else {
    grid$res2_echecs[i] <- table(res_int2[res_int2=="echec"])
  }
  
  if (all(is.na(res_int3[res_int3=="succes"]))) {
    grid$res3_succes[i] <- 0
  } else {
    grid$res3_succes[i] <- table(res_int3[res_int3=="succes"])
  }
  if (all(is.na(res_int3[res_int3=="echec"]))) {
    grid$res3_echecs[i] <- 0
  } else {
    grid$res3_echecs[i] <- table(res_int3[res_int3=="echec"])
  }
  
  #raison de l'arrêt à l'étape 1
  if (all(is.na(result1[result1=="tox:excessive eff:inadequate"]))) {
    grid$res1_toxE_repI[i] <- 0
  } else {
    grid$res1_toxE_repI[i] <- table(result1[result1=="tox:excessive eff:inadequate"])
  }
  if (all(is.na(result1[result1=="tox:excessive eff:OK"]))) {
    grid$res1_toxE_repO[i] <- 0
  } else {
    grid$res1_toxE_repO[i] <- table(result1[result1=="tox:excessive eff:OK"])
  }
  if (all(is.na(result1[result1=="tox:OK eff:inadequate"]))) {
    grid$res1_toxO_repI[i] <- 0
  } else {
    grid$res1_toxO_repI[i] <- table(result1[result1=="tox:OK eff:inadequate"])
  }
  #raison de l'arrêt à l'étape 2
  if (all(is.na(result2[result2=="tox:excessive eff:inadequate"]))) {
    grid$res2_toxE_repI[i] <- 0
  } else {
    grid$res2_toxE_repI[i] <- table(result2[result2=="tox:excessive eff:inadequate"])
  }
  if (all(is.na(result2[result2=="tox:excessive eff:OK"]))) {
    grid$res2_toxE_repO[i] <- 0
  } else {
    grid$res2_toxE_repO[i] <- table(result2[result2=="tox:excessive eff:OK"])
  }
  if (all(is.na(result2[result2=="tox:OK eff:inadequate"]))) {
    grid$res2_toxO_repI[i] <- 0
  } else {
    grid$res2_toxO_repI[i] <- table(result2[result2=="tox:OK eff:inadequate"])
  }
  #raison de l'arrêt à l'étape 2
  if (all(is.na(result3[result3=="tox:excessive eff:inadequate"]))) {
    grid$res3_toxE_repI[i] <- 0
  } else {
    grid$res3_toxE_repI[i] <- table(result3[result3=="tox:excessive eff:inadequate"])
  }
  if (all(is.na(result3[result3=="tox:excessive eff:OK"]))) {
    grid$res3_toxE_repO[i] <- 0
  } else {
    grid$res3_toxE_repO[i] <- table(result3[result3=="tox:excessive eff:OK"])
  }
  if (all(is.na(result3[result3=="tox:OK eff:inadequate"]))) {
    grid$res3_toxO_repI[i] <- 0
  } else {
    grid$res3_toxO_repI[i] <- table(result2[result2=="tox:OK eff:inadequate"])
  }
  
  
  #% arrêt pour tox
  grid$stop_tox[i] <- grid$res1_toxE_repI[i]+grid$res1_toxE_repO[i]+grid$res2_toxE_repI[i]+grid$res2_toxE_repO[i]+grid$res3_toxE_repI[i]+grid$res3_toxE_repO[i]
  
  #vérification s'il y a pb lors des AI (patients inclus en même temps)
  if (all(is.na(check_stop1[check_stop1=="STOP"]))) {
    grid$check_AI1[i] <- 0
  } else {
    grid$check_AI1[i] <- table(check_stop1[check_stop1=="STOP"])
  }
  if (all(is.na(check_stop2[check_stop2=="STOP"]))) {
    grid$check_AI2[i] <- 0
  } else {
    grid$check_AI2[i] <- table(check_stop2[check_stop2=="STOP"])
  }
  
}

writexl::write_xlsx(grid, "C:\\util\\marie\\2022-2023\\stage\\R\\résultats\\BOP2efftox_tpsfixe.xlsx")

table(res_int1)
table(res_int2)
table(res_int3)

table(result1)
table(result2)
table(result3)

(sample1+sample2+sample3)/ntrials

#nb d'efficacités
(nb_rep1+nb_rep2+nb_rep3)/ntrials

#nb d'échecs
(nb_ech1+nb_ech2+nb_ech3)/ntrials

#nb de toxicités
(nb_tox1+nb_tox2+nb_tox3)/ntrials

#durée de l'essai
mean(time)

