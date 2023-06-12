# SIMON + TOX 

######### Simon + TOX modif simulation patients #########

#simulation des patients -> lors des AI, les patients ne peuvent pas être inclus en même temps 

param_eff=c(0.2,0.4,0.4,0.2,0.4,0.2,0.4,0.3)
param_tox=c(0.35,0.15,0.05,0.05,0.45,0.15,0.35,0.2)
scenario=paste("Scenario",1:8)
grid=cbind(data.frame(scenario,param_eff,param_tox))
colnames(grid)=c("Scenario","Peff","Ptox")

n1=6 ; n2=25 ; N=60
r1=4 ; R=15 #mêmes règles que Simon + safety
ct1=3 ; ct2=8 ; CT=16


ntrials=10000
npat=N

nb_eff=nb_ech=nb_tox=0
res_etape1=res_etape2=res_etape3=rep(NA,ntrials)
sample1=sample2=sample3=0
result1=result2=result3=rep(NA,ntrials) #pour compter le nb d'échecs pour toxicité 
time=rep(NA,ntrials)
t_obs=4
theta0=0.2
tau=0.9



for (i in 1:nrow(grid))  {
  
  print(paste("LIGNE",i))
  
  peff<-grid$Peff[i]
  ptox<-grid$Ptox[i]
  print(paste("peff=",peff,"ptox=",ptox))
  
  pefftox=peff*ptox
  peffnotox=peff-pefftox
  pnoefftox=ptox-pefftox
  pnoeffnotox=(1-peff)-pnoefftox
  
  nb_eff=nb_ech=nb_tox=0
  res_etape1=res_etape2=res_etape3=rep(NA,10000)
  sample1=sample2=sample3=0
  result1=result2=result3=rep(NA,10000) #pour compter le nb d'échecs pour toxicité 
  check_stop1=check_stop2=rep(NA,ntrials)
  
  set.seed(3003)
  for (n in 1:ntrials) {
    
    #Simulation des patients
    #patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))),cumsum(c(0,rpois(N-1,2))))
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
    
    #ETAPE 1 (évaluation de la toxicité uniquement)
    stage1=patient[1:n1,]
    eff1=sum(stage1[,1])+sum(stage1[,2])
    ech1=sum(stage1[,3])+sum(stage1[,4])
    tox1=sum(stage1[,2])+sum(stage1[,4])
    
    if (tox1<ct1) {res_etape1[n]<-"succes"} else {res_etape1[n]<-"echec"}
    #if (pbeta(theta0, a0+tox1, b0+(n1-tox1), lower.tail=F) < tau) {res_etape1[n]<-"succes"} else {res_etape1[n]<-"echec"}
    sample1=sample1+n1
    nb_eff=nb_eff+eff1
    nb_ech=nb_ech+ech1
    nb_tox=nb_tox+tox1
    
    #raison de l'arrêt
    if (res_etape1[n]=="echec") {result1[n]<-"tox:excessive"}
    
    #durée de l'essai
    if (res_etape1[n]=="echec") {time[n]<-patient[n1,5]+t_obs}
    
    #vérification inclusion AI
    if (patient[n1,5]==patient[n1+1,5]) {check_stop1[n] <- "STOP"} else {"OK"}
    if (patient[n2,5]==patient[n2+1,5]) {check_stop2[n] <- "STOP"} else {"OK"}
    
    
    #ETAPE 2 (évaluation de l'efficacité et de la toxicité)
    if (res_etape1[n]=="succes") {
      stage2=patient[(n1+1):n2,]
      eff2=sum(stage2[,1])+sum(stage2[,2])
      ech2=sum(stage2[,3])+sum(stage2[,4])
      tox2=sum(stage2[,2])+sum(stage2[,4])
      
      if (eff2>r1 && tox1+tox2<ct2) {res_etape2[n]<-"succes"} else {res_etape2[n]<-"echec"}
      #if (pbeta(theta0, a0+tox2+tox1, b0+(n2-tox2-tox1), lower.tail=F) < tau && eff2 > r1) {res_etape2[n]<-"succes"} else {res_etape2[n]<-"echec"}
      sample2=sample2+n2-n1
      nb_eff=nb_eff+eff2
      nb_ech=nb_ech+ech2
      nb_tox=nb_tox+tox2
      
      #raison de l'arrêt
      if (res_etape2[n]=="echec") {
        if (eff2<=r1 && tox1+tox2>=ct2) {result2[n]<-"tox:excessive ; rep:inadequate"}
        if (eff2 >r1 && tox1+tox2>=ct2) {result2[n]<-"tox:excessive ; rep:OK"}
        if (eff2<=r1 && tox1+tox2 <ct2) {result2[n]<-"tox:OK ; rep:inadequate"}
      }
      
      #durée de l'essai
      if (res_etape2[n]=="echec") {
        if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])==ct1-1) {time[n]<-patient[n2,5]+2*t_obs}
        if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])!=ct1-1) {time[n]<-patient[n2,5]+t_obs}
      }
      
      
      #ETAPE 3 (évaluation de l'efficacité et de la toxicité)
      if (res_etape2[n]=="succes") {
        stage3=patient[(n2+1):N,]
        eff3=sum(stage3[,1])+sum(stage3[,2])
        ech3=sum(stage3[,3])+sum(stage3[,4])
        tox3=sum(stage3[,2])+sum(stage3[,4])
        
        if (eff2+eff3>R && tox1+tox2+tox3<CT) {res_etape3[n]<-"succes"} else {res_etape3[n]<-"echec"}
        #if (pbeta(theta0, a0+tox3+tox2+tox1, b0+(N-tox3-tox2-tox1), lower.tail=F) < tau && eff2+eff3 > R) {res_etape3[n]<-"succes"} else {res_etape3[n]<-"echec"}
        sample3=sample3+N-n2
        nb_eff=nb_eff+eff3
        nb_ech=nb_ech+ech3
        nb_tox=nb_tox+tox3
        
        #raison de l'arrêt
        if (res_etape3[n]=="echec") {
          if (eff2+eff3<=R && tox1+tox2+tox3>=CT) {result3[n]<-"tox:excessive ; rep:inadequate"}
          if (eff2+eff3 >R && tox1+tox2+tox3>=CT) {result3[n]<-"tox:excessive ; rep:OK"}
          if (eff2+eff3<=R && tox1+tox2+tox3 <CT) {result3[n]<-"tox:OK ; rep:inadequate"}
        }
        
        #durée de l'essai 
        if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])==ct1-1) {
          if (sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])==r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])==ct2-1) {time[n]<-patient[N,5]+3*t_obs}
          if (sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])!=r1 && tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])!=ct2-1) {time[n]<-patient[N,5]+2*t_obs}
        }
        if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])!=ct1-1) {
          if (sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])==r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])==ct2-1) {time[n]<-patient[N,5]+2*t_obs}
          if (sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])!=r1 && tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])!=ct2-1) {time[n]<-patient[N,5]+t_obs}
        }
        
      } 
      
    }
    
  }
  
  #% de succès total
  if (all(is.na(res_etape3[res_etape3=="succes"]))) {
    grid$H11[i] <- 0
  } else {
    grid$H11[i] <- table(res_etape3[res_etape3=="succes"])/100
  }
  print(paste(grid$H11[i],"% de succès"))
  
  #% arrêt précoce
  if (all(is.na(res_etape2[res_etape2=="echec"])) && all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$early_ter[i] <- 0
  } else {
    grid$early_ter[i] <- ((table(res_etape2[res_etape2=="echec"])+table(res_etape1[res_etape1=="echec"]))*100)/ntrials
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
  
  
  #nb d'échecs et de succès à chaque étape 
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
  
  if (all(is.na(res_etape3[res_etape3=="succes"]))) {
    grid$res3_succes[i] <- 0
  } else {
    grid$res3_succes[i] <- table(res_etape3[res_etape3=="succes"])
  }
  if (all(is.na(res_etape3[res_etape3=="echec"]))) {
    grid$res3_echecs[i] <- 0
  } else {
    grid$res3_echecs[i] <- table(res_etape3[res_etape3=="echec"])
  }
  
  
  #raison de l'arrêt à l'étape 1
  if (all(is.na(result1[result1=="tox:excessive"]))) {
    grid$res1_toxE[i] <- 0
  } else {
    grid$res1_toxE[i] <- table(result1[result1=="tox:excessive"])
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
  #raison de l'arrêt à l'étape 3
  if (all(is.na(result3[result3=="tox:excessive ; rep:inadequate"]))) {
    grid$res3_toxE_repI[i] <- 0
  } else {
    grid$res3_toxE_repI[i] <- table(result3[result3=="tox:excessive ; rep:inadequate"])
  }
  if (all(is.na(result3[result3=="tox:excessive ; rep:OK"]))) {
    grid$res3_toxE_repO[i] <- 0
  } else {
    grid$res3_toxE_repO[i] <- table(result3[result3=="tox:excessive ; rep:OK"])
  }
  if (all(is.na(result3[result3=="tox:OK ; rep:inadequate"]))) {
    grid$res3_toxO_repI[i] <- 0
  } else {
    grid$res3_toxO_repI[i] <- table(result3[result3=="tox:OK ; rep:inadequate"])
  }
  
  #% arrêt pour tox
  grid$stop_tox[i] <- grid$res1_toxE[i]+grid$res2_toxE_repI[i]+grid$res2_toxE_repO[i]+grid$res3_toxE_repI[i]+grid$res3_toxE_repO[i]
  
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

writexl::write_xlsx(grid, "C:\\util\\marie\\2022-2023\\stage\\R\\résultats\\Simon+TOX_tpspoisAI.xlsx")

######### Simon + TOX (grid search) #########
#code avec matrice : résultats des 8 scenarios pour les 56 valeurs possibles de theta0 et tau

paramA=seq(0.2,0.4,0.05)
paramB=seq(0.85,0.98,0.01)
grid=expand.grid(paramA,paramB)
colnames(grid)=c("theta0","tau")

a0=1
b0=1

peff=0.4
ptox=0.35

pefftox=peff*ptox
peffnotox=peff-pefftox
pnoefftox=ptox-pefftox
pnoeffnotox=(1-peff)-pnoefftox
pefftox+peffnotox+pnoefftox+pnoeffnotox==1

n1=6 ; n2=25 ; N=60
r1=4 ; R=15 #mêmes règles que Simon + safety


ntrials=10000
npat=N

nb_eff=nb_ech=nb_tox=0
res_etape1=res_etape2=res_etape3=rep(NA,ntrials)
sample1=sample2=sample3=0
result1=result2=result3=rep(NA,ntrials) #pour compter le nb d'échecs pour toxicité 
time=rep(NA,ntrials)
t_obs=4


set.seed(3003)
for (i in 1:nrow(grid))  {
  
  print(paste("MATRICE LIGNE",i))
  
  theta0 <- grid$theta0[i]
  tau <- grid$tau[i]
  print(paste("theta0=",theta0,"tau=",tau))
  
  #remise des compteurs à 0 pour chaque couple (theta0,tau)
  nb_eff=nb_ech=nb_tox=0
  res_etape1=res_etape2=res_etape3=rep(NA,10000)
  sample1=sample2=sample3=0
  result1=result2=result3=rep(NA,10000) #pour compter le nb d'échecs pour toxicité 
  
  set.seed(3003)
  for (n in 1:ntrials) {
    
    #Simulation des patients
    #patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))),cumsum(c(0,rep(2,N-1))))
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
    
    #ETAPE 1 (évaluation de la toxicité uniquement)
    stage1=patient[1:n1,]
    eff1=sum(stage1[,1])+sum(stage1[,2])
    ech1=sum(stage1[,3])+sum(stage1[,4])
    tox1=sum(stage1[,2])+sum(stage1[,4])
    
    #if (tox1<ct1) {res_etape1[n]<-"succes"} else {res_etape1[n]<-"echec"}
    if (pbeta(theta0, a0+tox1, b0+(n1-tox1), lower.tail=F) < tau) {res_etape1[n]<-"succes"} else {res_etape1[n]<-"echec"}
    sample1=sample1+n1
    nb_eff=nb_eff+eff1
    nb_ech=nb_ech+ech1
    nb_tox=nb_tox+tox1
    
    #raison de l'arrêt
    if (res_etape1[n]=="echec") {result1[n]<-"tox:excessive"}
    
    #durée de l'essai
    if (res_etape1[n]=="echec") {time[n]<-patient[n1,5]+t_obs}
    # if (res_etape1[n]=="echec") {
    #   if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])==ct1) {time[n]<-patient[n1,5]+t_obs}
    #   if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])!=ct1) {time[n]<-patient[n1-1,5]+t_obs}
    # }
    
    
    #ETAPE 2 (évaluation de l'efficacité et de la toxicité)
    if (res_etape1[n]=="succes") {
      stage2=patient[(n1+1):n2,]
      eff2=sum(stage2[,1])+sum(stage2[,2])
      ech2=sum(stage2[,3])+sum(stage2[,4])
      tox2=sum(stage2[,2])+sum(stage2[,4])
      
      #if (eff1+eff2>r1 && tox1+tox2<ct2) {res_etape2[n]<-"succes"} else {res_etape2[n]<-"echec"}
      if (pbeta(theta0, a0+tox2+tox1, b0+(n2-tox2-tox1), lower.tail=F) < tau && eff2 > r1) {res_etape2[n]<-"succes"} else {res_etape2[n]<-"echec"}
      sample2=sample2+n2-n1
      nb_eff=nb_eff+eff2
      nb_ech=nb_ech+ech2
      nb_tox=nb_tox+tox2
      
      #raison de l'arrêt
      if (res_etape2[n]=="echec") {
        if (eff1+eff2<=r1 && pbeta(theta0, a0+tox2+tox1, b0+(n2-tox2-tox1), lower.tail=F)>=tau) {result2[n]<-"tox:excessive ; rep:inadequate"}
        if (eff1+eff2 >r1 && pbeta(theta0, a0+tox2+tox1, b0+(n2-tox2-tox1), lower.tail=F)>=tau) {result2[n]<-"tox:excessive ; rep:OK"}
        if (eff1+eff2<=r1 && pbeta(theta0, a0+tox2+tox1, b0+(n2-tox2-tox1), lower.tail=F) <tau) {result2[n]<-"tox:OK ; rep:inadequate"}
      }
      
      #durée de l'essai
      if (res_etape2[n]=="echec") {time[n]<-patient[n2,5]+t_obs}
      # if (res_etape2[n]=="echec") {
      #   if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])==ct1) { #{time[n]<-patient[n1,5]+t_obs}
      #     if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])==r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])==ct2) {time[n]<-patient[n2,5]+2*t_obs}
      #     if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])!=r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])!=ct2) {time[n]<-patient[n2-1,5]+2*t_obs}
      #   } 
      #   if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])!=ct1) { #{time[n]<-patient[n1-1,5]+t_obs}
      #     if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])==r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])==ct2) {time[n]<-patient[n2,5]+t_obs}
      #     if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])!=r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])!=ct2) {time[n]<-patient[n2-1,5]+t_obs}
      #   } 
      # }
      
      
      #ETAPE 3 (évaluation de l'efficacité et de la toxicité)
      if (res_etape2[n]=="succes") {
        stage3=patient[(n2+1):N,]
        eff3=sum(stage3[,1])+sum(stage3[,2])
        ech3=sum(stage3[,3])+sum(stage3[,4])
        tox3=sum(stage3[,2])+sum(stage3[,4])
        
        #if (eff1+eff2+eff3>R && tox1+tox2+tox3<TO) {res_etape3[n]<-"succes"} else {res_etape3[n]<-"echec"}
        if (pbeta(theta0, a0+tox3+tox2+tox1, b0+(N-tox3-tox2-tox1), lower.tail=F) < tau && eff2+eff3 > R) {res_etape3[n]<-"succes"} else {res_etape3[n]<-"echec"}
        sample3=sample3+N-n2
        nb_eff=nb_eff+eff3
        nb_ech=nb_ech+ech3
        nb_tox=nb_tox+tox3
        
        #raison de l'arrêt
        if (res_etape3[n]=="echec") {
          if (eff1+eff2+eff3<=R && pbeta(theta0, a0+tox3+tox2+tox1, b0+(N-tox3-tox2-tox1), lower.tail=F)>=tau) {result3[n]<-"tox:excessive ; rep:inadequate"}
          if (eff1+eff2+eff3 >R && pbeta(theta0, a0+tox3+tox2+tox1, b0+(N-tox3-tox2-tox1), lower.tail=F)>=tau) {result3[n]<-"tox:excessive ; rep:OK"}
          if (eff1+eff2+eff3<=R && pbeta(theta0, a0+tox3+tox2+tox1, b0+(N-tox3-tox2-tox1), lower.tail=F) <tau) {result3[n]<-"tox:OK ; rep:inadequate"}
        }
        
        #durée de l'essai AJOUTER LES CONDITIONS POUR L'ETAPE 3
        if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])==r1) {time[n]<-patient[N,5]+2*t_obs}
        if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])!=r1) {time[n]<-patient[N,5]+t_obs}
        # if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])==ct1) { 
        #   if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])==r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])==ct2) {
        #     if (eff1+eff2+sum(stage3[1:(N-n2-n1-1),1])+sum(stage3[1:(N-n2-n1-1),2])==R | tox1+tox2+sum(stage3[1:(N-n2-n1-1),2])+sum(stage3[1:(N-n2-n1-1),4])==TO) {time[n]<-patient[N,5]+3*t_obs}
        #     if (eff1+eff2+sum(stage3[1:(N-n2-n1-1),1])+sum(stage3[1:(N-n2-n1-1),2])!=R | tox1+tox2+sum(stage3[1:(N-n2-n1-1),2])+sum(stage3[1:(N-n2-n1-1),4])==TO) {time[n]<-patient[N-1,5]+3*t_obs}
        #   }
        #   if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])!=r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])!=ct2) {
        #     if (eff1+eff2+sum(stage3[1:(N-n2-n1-1),1])+sum(stage3[1:(N-n2-n1-1),2])==R | tox1+tox2+sum(stage3[1:(N-n2-n1-1),2])+sum(stage3[1:(N-n2-n1-1),4])==TO) {time[n]<-patient[N,5]+2*t_obs}
        #     if (eff1+eff2+sum(stage3[1:(N-n2-n1-1),1])+sum(stage3[1:(N-n2-n1-1),2])!=R | tox1+tox2+sum(stage3[1:(N-n2-n1-1),2])+sum(stage3[1:(N-n2-n1-1),4])==TO) {time[n]<-patient[N-1,5]+2*t_obs}
        #   }
        # } 
        # if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])!=ct1) { 
        #   if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])==r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])==ct2) {#{time[n]<-patient[n2,5]+t_obs}
        #     if (eff1+eff2+sum(stage3[1:(N-n2-n1-1),1])+sum(stage3[1:(N-n2-n1-1),2])==R | tox1+tox2+sum(stage3[1:(N-n2-n1-1),2])+sum(stage3[1:(N-n2-n1-1),4])==TO) {time[n]<-patient[N,5]+2*t_obs}
        #     if (eff1+eff2+sum(stage3[1:(N-n2-n1-1),1])+sum(stage3[1:(N-n2-n1-1),2])!=R | tox1+tox2+sum(stage3[1:(N-n2-n1-1),2])+sum(stage3[1:(N-n2-n1-1),4])==TO) {time[n]<-patient[N-1,5]+2*t_obs}
        #   }
        #   if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])!=r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])!=ct2) {#{time[n]<-patient[n2-1,5]+t_obs}
        #     if (eff1+eff2+sum(stage3[1:(N-n2-n1-1),1])+sum(stage3[1:(N-n2-n1-1),2])==R | tox1+tox2+sum(stage3[1:(N-n2-n1-1),2])+sum(stage3[1:(N-n2-n1-1),4])==TO) {time[n]<-patient[N,5]+t_obs}
        #     if (eff1+eff2+sum(stage3[1:(N-n2-n1-1),1])+sum(stage3[1:(N-n2-n1-1),2])!=R | tox1+tox2+sum(stage3[1:(N-n2-n1-1),2])+sum(stage3[1:(N-n2-n1-1),4])==TO) {time[n]<-patient[N-1,5]+t_obs}
        #   }
        #   
        # }
        
      } 
      
    }
    
  }
  
  #% de succès total
  if (all(is.na(res_etape3[res_etape3=="succes"]))) {
    grid$H11[i] <- 0
  } else {
    grid$H11[i] <- table(res_etape3[res_etape3=="succes"])/100
  }
  print(paste(grid$H11[i],"% de succès"))
  
  #nb d'échecs et de succès à chaque étape 
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
  
  if (all(is.na(res_etape3[res_etape3=="succes"]))) {
    grid$res3_succes[i] <- 0
  } else {
    grid$res3_succes[i] <- table(res_etape3[res_etape3=="succes"])
  }
  if (all(is.na(res_etape3[res_etape3=="echec"]))) {
    grid$res3_echecs[i] <- 0
  } else {
    grid$res3_echecs[i] <- table(res_etape3[res_etape3=="echec"])
  }
  
  
  #raison de l'arrêt à l'étape 1
  if (all(is.na(result1[result1=="tox:excessive"]))) {
    grid$res1_toxE[i] <- 0
  } else {
    grid$res1_toxE[i] <- table(result1[result1=="tox:excessive"])
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
  #raison de l'arrêt à l'étape 3
  if (all(is.na(result3[result3=="tox:excessive ; rep:inadequate"]))) {
    grid$res3_toxE_repI[i] <- 0
  } else {
    grid$res3_toxE_repI[i] <- table(result3[result3=="tox:excessive ; rep:inadequate"])
  }
  if (all(is.na(result3[result3=="tox:excessive ; rep:OK"]))) {
    grid$res3_toxE_repO[i] <- 0
  } else {
    grid$res3_toxE_repO[i] <- table(result3[result3=="tox:excessive ; rep:OK"])
  }
  if (all(is.na(result3[result3=="tox:OK ; rep:inadequate"]))) {
    grid$res3_toxO_repI[i] <- 0
  } else {
    grid$res3_toxO_repI[i] <- table(result3[result3=="tox:OK ; rep:inadequate"])
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
  #grid$moy_time[i] <- mean(time)

}

writexl::write_xlsx(grid, "C:\\util\\marie\\2022-2023\\stage\\R\\résultats Simon + tox\\H10.xlsx")

table(res_etape1)
table(res_etape2)
table(res_etape3)

table(result1)
table(result2)
table(result3)

#effecifs
(sample1+sample2+sample3)/ntrials

#nb d'efficacités 
nb_eff/ntrials

#nb d'échecs
nb_ech/ntrials

#nb de toxicités 
nb_tox/ntrials

#durée de l'essai
mean(time)



a0=1
b0=1
#étape 1
pbeta(0.35, a0+0, b0+6, lower.tail = F) # =0.05
pbeta(0.35, a0+1, b0+5, lower.tail = F) # =0.23
pbeta(0.35, a0+2, b0+4, lower.tail = F) # =0.53
pbeta(0.35, a0+3, b0+3, lower.tail = F) # =0.80
pbeta(0.35, a0+4, b0+2, lower.tail = F) # =0.94
pbeta(0.35, a0+5, b0+1, lower.tail = F) # =0.99
pbeta(0.35, a0+6, b0+0, lower.tail = F) # =0.99
# à partir de 4 toxicités -> échec 

a0=1
b0=1
N=seq(1,60,1)
npat=seq(1,60,1)
distrib=rep(NA,60)
valeurs=expand.grid(npat,N)

for (n in 1:60) {
  
  print(paste(n,"patients"))
  N <- n
  distrib=rep(NA,60)
  
  for (i in 1:60) {
    
    distrib[i] <- pbeta(0.09, a0+i, b0+(N-i), lower.tail=F)
    print(paste("tox:",a0+i,"notox:",b0+(N-i),"distrib:",distrib[i]))
    seuil=table(distrib>=0.95)
    
    valeurs$nb_pat[i] <- N
    valeurs$distrib[i] <- seuil
    valeurs$tox[i] <- a0+i
    valeurs$notox[i] <- b0+(N-i)
    
  }
  
}


table(distrib>0.95)

######### Simon + TOX (avec tau et theta0) #############
# code pour scenarios (avec pbeta())

param_eff=c(0.2,0.4,0.4,0.2,0.4,0.2,0.4,0.3)
param_tox=c(0.35,0.15,0.05,0.05,0.45,0.15,0.35,0.2)
scenario=paste("Scenario",1:8)
grid=cbind(data.frame(scenario,param_eff,param_tox))
colnames(grid)=c("Scenario","Peff","Ptox")

a0=1
b0=1

peff=0.3
ptox=0.2

pefftox=peff*ptox
peffnotox=peff-pefftox
pnoefftox=ptox-pefftox
pnoeffnotox=(1-peff)-pnoefftox
pefftox+peffnotox+pnoefftox+pnoeffnotox==1

n1=6 ; n2=25 ; N=60
r1=4 ; R=15 #mêmes règles que Simon + safety


ntrials=10000
npat=N

nb_eff=nb_ech=nb_tox=0
res_etape1=res_etape2=res_etape3=rep(NA,ntrials)
sample1=sample2=sample3=0
result1=result2=result3=rep(NA,ntrials) #pour compter le nb d'échecs pour toxicité 
#time=rep(NA,ntrials)
#t_obs=4
theta0=0.2
tau=0.9

set.seed(3003)
for (i in 1:nrow(grid))  {
  
  print(paste("LIGNE",i))
  
  peff<-grid$Peff[i]
  ptox<-grid$Ptox[i]
  print(paste("peff=",peff,"ptox=",ptox))
  
  pefftox=peff*ptox
  peffnotox=peff-pefftox
  pnoefftox=ptox-pefftox
  pnoeffnotox=(1-peff)-pnoefftox
  
  #remise des compteurs à 0 pour chaque couple (theta0,tau)
  nb_eff=nb_ech=nb_tox=0
  res_etape1=res_etape2=res_etape3=rep(NA,10000)
  sample1=sample2=sample3=0
  result1=result2=result3=rep(NA,10000) #pour compter le nb d'échecs pour toxicité 
  
  set.seed(3003)
  for (n in 1:ntrials) {
    
    #Simulation des patients
    patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))),cumsum(c(0,rep(2,N-1))))
    
    #ETAPE 1 (évaluation de la toxicité uniquement)
    stage1=patient[1:n1,]
    eff1=sum(stage1[,1])+sum(stage1[,2])
    ech1=sum(stage1[,3])+sum(stage1[,4])
    tox1=sum(stage1[,2])+sum(stage1[,4])
    notox1=sum(stage1[,1])+sum(stage1[,3])
    
    #if (tox1<ct1) {res_etape1[n]<-"succes"} else {res_etape1[n]<-"echec"}
    if (pbeta(theta0, a0+tox1, b0+notox1, lower.tail=F) < tau) {res_etape1[n]<-"succes"} else {res_etape1[n]<-"echec"}
    sample1=sample1+n1
    nb_eff=nb_eff+eff1
    nb_ech=nb_ech+ech1
    nb_tox=nb_tox+tox1
    
    #raison de l'arrêt
    if (res_etape1[n]=="echec") {result1[n]<-"tox:excessive"}
    
    #durée de l'essai
    # if (res_etape1[n]=="echec") {
    #   if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])==ct1) {time[n]<-patient[n1,5]+t_obs}
    #   if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])!=ct1) {time[n]<-patient[n1-1,5]+t_obs}
    # }
    
    
    #ETAPE 2 (évaluation de l'efficacité et de la toxicité)
    if (res_etape1[n]=="succes") {
      stage2=patient[(n1+1):n2,]
      eff2=sum(stage2[,1])+sum(stage2[,2])
      ech2=sum(stage2[,3])+sum(stage2[,4])
      tox2=sum(stage2[,2])+sum(stage2[,4])
      notox2=sum(stage2[,1])+sum(stage2[,3])
      
      #if (eff1+eff2>r1 && tox1+tox2<ct2) {res_etape2[n]<-"succes"} else {res_etape2[n]<-"echec"}
      if (pbeta(theta0, a0+tox2+tox1, b0+notox2+notox1, lower.tail=F) < tau && eff2 > r1) {res_etape2[n]<-"succes"} else {res_etape2[n]<-"echec"}
      sample2=sample2+n2-n1
      nb_eff=nb_eff+eff2
      nb_ech=nb_ech+ech2
      nb_tox=nb_tox+tox2
      
      #raison de l'arrêt
      if (res_etape2[n]=="echec") {
        if (eff1+eff2<=r1 && pbeta(theta0, a0+tox2+tox1, b0+(n2-tox2-tox1), lower.tail=F)>=tau) {result2[n]<-"tox:excessive ; rep:inadequate"}
        if (eff1+eff2 >r1 && pbeta(theta0, a0+tox2+tox1, b0+(n2-tox2-tox1), lower.tail=F)>=tau) {result2[n]<-"tox:excessive ; rep:OK"}
        if (eff1+eff2<=r1 && pbeta(theta0, a0+tox2+tox1, b0+(n2-tox2-tox1), lower.tail=F) <tau) {result2[n]<-"tox:OK ; rep:inadequate"}
      }
      
      #durée de l'essai
      # if (res_etape2[n]=="echec") {
      #   if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])==ct1) { #{time[n]<-patient[n1,5]+t_obs}
      #     if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])==r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])==ct2) {time[n]<-patient[n2,5]+2*t_obs}
      #     if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])!=r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])!=ct2) {time[n]<-patient[n2-1,5]+2*t_obs}
      #   } 
      #   if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])!=ct1) { #{time[n]<-patient[n1-1,5]+t_obs}
      #     if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])==r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])==ct2) {time[n]<-patient[n2,5]+t_obs}
      #     if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])!=r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])!=ct2) {time[n]<-patient[n2-1,5]+t_obs}
      #   } 
      # }
      
      
      #ETAPE 3 (évaluation de l'efficacité et de la toxicité)
      if (res_etape2[n]=="succes") {
        stage3=patient[(n2+1):N,]
        eff3=sum(stage3[,1])+sum(stage3[,2])
        ech3=sum(stage3[,3])+sum(stage3[,4])
        tox3=sum(stage3[,2])+sum(stage3[,4])
        notox3=sum(stage3[,1])+sum(stage3[,3])
        
        #if (eff1+eff2+eff3>R && tox1+tox2+tox3<TO) {res_etape3[n]<-"succes"} else {res_etape3[n]<-"echec"}
        if (pbeta(theta0, a0+tox3+tox2+tox1, b0+notox3+notox2+notox1, lower.tail=F) < tau && eff2+eff3 > R) {res_etape3[n]<-"succes"} else {res_etape3[n]<-"echec"}
        sample3=sample3+N-n2
        nb_eff=nb_eff+eff3
        nb_ech=nb_ech+ech3
        nb_tox=nb_tox+tox3
        
        #raison de l'arrêt
        if (res_etape3[n]=="echec") {
          if (eff1+eff2+eff3<=R && pbeta(theta0, a0+tox3+tox2+tox1, b0+(N-tox3-tox2-tox1), lower.tail=F)>=tau) {result3[n]<-"tox:excessive ; rep:inadequate"}
          if (eff1+eff2+eff3 >R && pbeta(theta0, a0+tox3+tox2+tox1, b0+(N-tox3-tox2-tox1), lower.tail=F)>=tau) {result3[n]<-"tox:excessive ; rep:OK"}
          if (eff1+eff2+eff3<=R && pbeta(theta0, a0+tox3+tox2+tox1, b0+(N-tox3-tox2-tox1), lower.tail=F) <tau) {result3[n]<-"tox:OK ; rep:inadequate"}
        }
        
        #durée de l'essai AJOUTER LES CONDITIONS POUR L'ETAPE 3
        # if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])==ct1) { 
        #   if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])==r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])==ct2) {
        #     if (eff1+eff2+sum(stage3[1:(N-n2-n1-1),1])+sum(stage3[1:(N-n2-n1-1),2])==R | tox1+tox2+sum(stage3[1:(N-n2-n1-1),2])+sum(stage3[1:(N-n2-n1-1),4])==TO) {time[n]<-patient[N,5]+3*t_obs}
        #     if (eff1+eff2+sum(stage3[1:(N-n2-n1-1),1])+sum(stage3[1:(N-n2-n1-1),2])!=R | tox1+tox2+sum(stage3[1:(N-n2-n1-1),2])+sum(stage3[1:(N-n2-n1-1),4])==TO) {time[n]<-patient[N-1,5]+3*t_obs}
        #   }
        #   if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])!=r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])!=ct2) {
        #     if (eff1+eff2+sum(stage3[1:(N-n2-n1-1),1])+sum(stage3[1:(N-n2-n1-1),2])==R | tox1+tox2+sum(stage3[1:(N-n2-n1-1),2])+sum(stage3[1:(N-n2-n1-1),4])==TO) {time[n]<-patient[N,5]+2*t_obs}
        #     if (eff1+eff2+sum(stage3[1:(N-n2-n1-1),1])+sum(stage3[1:(N-n2-n1-1),2])!=R | tox1+tox2+sum(stage3[1:(N-n2-n1-1),2])+sum(stage3[1:(N-n2-n1-1),4])==TO) {time[n]<-patient[N-1,5]+2*t_obs}
        #   }
        # } 
        # if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])!=ct1) { 
        #   if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])==r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])==ct2) {#{time[n]<-patient[n2,5]+t_obs}
        #     if (eff1+eff2+sum(stage3[1:(N-n2-n1-1),1])+sum(stage3[1:(N-n2-n1-1),2])==R | tox1+tox2+sum(stage3[1:(N-n2-n1-1),2])+sum(stage3[1:(N-n2-n1-1),4])==TO) {time[n]<-patient[N,5]+2*t_obs}
        #     if (eff1+eff2+sum(stage3[1:(N-n2-n1-1),1])+sum(stage3[1:(N-n2-n1-1),2])!=R | tox1+tox2+sum(stage3[1:(N-n2-n1-1),2])+sum(stage3[1:(N-n2-n1-1),4])==TO) {time[n]<-patient[N-1,5]+2*t_obs}
        #   }
        #   if (eff1+sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])!=r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])!=ct2) {#{time[n]<-patient[n2-1,5]+t_obs}
        #     if (eff1+eff2+sum(stage3[1:(N-n2-n1-1),1])+sum(stage3[1:(N-n2-n1-1),2])==R | tox1+tox2+sum(stage3[1:(N-n2-n1-1),2])+sum(stage3[1:(N-n2-n1-1),4])==TO) {time[n]<-patient[N,5]+t_obs}
        #     if (eff1+eff2+sum(stage3[1:(N-n2-n1-1),1])+sum(stage3[1:(N-n2-n1-1),2])!=R | tox1+tox2+sum(stage3[1:(N-n2-n1-1),2])+sum(stage3[1:(N-n2-n1-1),4])==TO) {time[n]<-patient[N-1,5]+t_obs}
        #   }
        #   
        # }
        
      } 
      
    }
    
  }
  
  #% de succès total
  if (all(is.na(res_etape3[res_etape3=="succes"]))) {
    grid$H11[i] <- 0
  } else {
    grid$H11[i] <- table(res_etape3[res_etape3=="succes"])/100
  }
  print(paste(grid$H11[i],"% de succès"))
  
  #% arrêt précoce
  if (all(is.na(res_etape2[res_etape2=="echec"])) && all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$early_ter[i] <- 0
  } else {
    grid$early_ter[i] <- ((table(res_etape2[res_etape2=="echec"])+table(res_etape1[res_etape1=="echec"]))*100)/ntrials
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
  #grid$moy_time[i] <- mean(time)
  
  
  #nb d'échecs et de succès à chaque étape 
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
  
  if (all(is.na(res_etape3[res_etape3=="succes"]))) {
    grid$res3_succes[i] <- 0
  } else {
    grid$res3_succes[i] <- table(res_etape3[res_etape3=="succes"])
  }
  if (all(is.na(res_etape3[res_etape3=="echec"]))) {
    grid$res3_echecs[i] <- 0
  } else {
    grid$res3_echecs[i] <- table(res_etape3[res_etape3=="echec"])
  }
  
  
  # #raison de l'arrêt à l'étape 1
  # if (all(is.na(result1[result1=="tox:excessive"]))) {
  #   grid$res1_toxE[i] <- 0
  # } else {
  #   grid$res1_toxE[i] <- table(result1[result1=="tox:excessive"])
  # }
  # #raison de l'arrêt à l'étape 2
  # if (all(is.na(result2[result2=="tox:excessive ; rep:inadequate"]))) {
  #   grid$res2_toxE_repI[i] <- 0
  # } else {
  #   grid$res2_toxE_repI[i] <- table(result2[result2=="tox:excessive ; rep:inadequate"])
  # }
  # if (all(is.na(result2[result2=="tox:excessive ; rep:OK"]))) {
  #   grid$res2_toxE_repO[i] <- 0
  # } else {
  #   grid$res2_toxE_repO[i] <- table(result2[result2=="tox:excessive ; rep:OK"])
  # }
  # if (all(is.na(result2[result2=="tox:OK ; rep:inadequate"]))) {
  #   grid$res2_toxO_repI[i] <- 0
  # } else {
  #   grid$res2_toxO_repI[i] <- table(result2[result2=="tox:OK ; rep:inadequate"])
  # }
  # #raison de l'arrêt à l'étape 3
  # if (all(is.na(result3[result3=="tox:excessive ; rep:inadequate"]))) {
  #   grid$res3_toxE_repI[i] <- 0
  # } else {
  #   grid$res3_toxE_repI[i] <- table(result3[result3=="tox:excessive ; rep:inadequate"])
  # }
  # if (all(is.na(result3[result3=="tox:excessive ; rep:OK"]))) {
  #   grid$res3_toxE_repO[i] <- 0
  # } else {
  #   grid$res3_toxE_repO[i] <- table(result3[result3=="tox:excessive ; rep:OK"])
  # }
  # if (all(is.na(result3[result3=="tox:OK ; rep:inadequate"]))) {
  #   grid$res3_toxO_repI[i] <- 0
  # } else {
  #   grid$res3_toxO_repI[i] <- table(result3[result3=="tox:OK ; rep:inadequate"])
  # }
  
  #% arrêt pour tox
  #grid$stop_tox[i] <- grid$res1_toxE[i]+grid$res2_toxE_repI[i]+grid$res2_toxE_repO[i]+res3_toxE_repI[i]+res3_toxE_repO[i]
  
}

writexl::write_xlsx(grid, "C:\\util\\marie\\2022-2023\\stage\\R\\résultats\\Simon+TOX_codesansseuil.xlsx")





######### Simon + TOX (avec seuils trouvés) #############
# code pour scenarios (avec seuils trouvés)

param_eff=c(0.2,0.4,0.4,0.2,0.4,0.2,0.4,0.3)
param_tox=c(0.35,0.15,0.05,0.05,0.45,0.15,0.35,0.2)
scenario=paste("Scenario",1:8)
grid=cbind(data.frame(scenario,param_eff,param_tox))
colnames(grid)=c("Scenario","Peff","Ptox")

n1=6 ; n2=25 ; N=60
r1=4 ; R=15 #mêmes règles que Simon + safety
ct1=3 ; ct2=8 ; CT=16


ntrials=10000
npat=N

nb_eff=nb_ech=nb_tox=0
res_etape1=res_etape2=res_etape3=rep(NA,ntrials)
sample1=sample2=sample3=0
result1=result2=result3=rep(NA,ntrials) #pour compter le nb d'échecs pour toxicité 
time=rep(NA,ntrials)
t_obs=4
theta0=0.2
tau=0.9



for (i in 1:nrow(grid))  {
  
  print(paste("LIGNE",i))
  
  peff<-grid$Peff[i]
  ptox<-grid$Ptox[i]
  print(paste("peff=",peff,"ptox=",ptox))
  
  pefftox=peff*ptox
  peffnotox=peff-pefftox
  pnoefftox=ptox-pefftox
  pnoeffnotox=(1-peff)-pnoefftox
  
  nb_eff=nb_ech=nb_tox=0
  res_etape1=res_etape2=res_etape3=rep(NA,10000)
  sample1=sample2=sample3=0
  result1=result2=result3=rep(NA,10000) #pour compter le nb d'échecs pour toxicité 
  check_stop1=check_stop2=rep(NA,ntrials)
  
  set.seed(3003)
  for (n in 1:ntrials) {
    
    #Simulation des patients
    patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))),cumsum(c(0,rpois(N-1,2))))
    #patient=t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox)))
    #fixe=cumsum(c(0,rep(2,N-1)))
    #pois=cumsum(c(0,rpois(N-1,2)))
    #patient=cbind(patient,pois)
    
    #ETAPE 1 (évaluation de la toxicité uniquement)
    stage1=patient[1:n1,]
    eff1=sum(stage1[,1])+sum(stage1[,2])
    ech1=sum(stage1[,3])+sum(stage1[,4])
    tox1=sum(stage1[,2])+sum(stage1[,4])
    
    if (tox1<ct1) {res_etape1[n]<-"succes"} else {res_etape1[n]<-"echec"}
    #if (pbeta(theta0, a0+tox1, b0+(n1-tox1), lower.tail=F) < tau) {res_etape1[n]<-"succes"} else {res_etape1[n]<-"echec"}
    sample1=sample1+n1
    nb_eff=nb_eff+eff1
    nb_ech=nb_ech+ech1
    nb_tox=nb_tox+tox1
    
    #raison de l'arrêt
    if (res_etape1[n]=="echec") {result1[n]<-"tox:excessive"}
    
    #durée de l'essai
    if (res_etape1[n]=="echec") {time[n]<-patient[n1,5]+t_obs}
    
    #vérification inclusion AI
    if (patient[n1,5]==patient[n1+1,5]) {check_stop1[n] <- "STOP"} else {"OK"}
    if (patient[n2,5]==patient[n2+1,5]) {check_stop2[n] <- "STOP"} else {"OK"}
    
    
    #ETAPE 2 (évaluation de l'efficacité et de la toxicité)
    if (res_etape1[n]=="succes") {
      stage2=patient[(n1+1):n2,]
      eff2=sum(stage2[,1])+sum(stage2[,2])
      ech2=sum(stage2[,3])+sum(stage2[,4])
      tox2=sum(stage2[,2])+sum(stage2[,4])
      
      if (eff2>r1 && tox1+tox2<ct2) {res_etape2[n]<-"succes"} else {res_etape2[n]<-"echec"}
      #if (pbeta(theta0, a0+tox2+tox1, b0+(n2-tox2-tox1), lower.tail=F) < tau && eff2 > r1) {res_etape2[n]<-"succes"} else {res_etape2[n]<-"echec"}
      sample2=sample2+n2-n1
      nb_eff=nb_eff+eff2
      nb_ech=nb_ech+ech2
      nb_tox=nb_tox+tox2
      
      #raison de l'arrêt
      if (res_etape2[n]=="echec") {
        if (eff2<=r1 && tox1+tox2>=ct2) {result2[n]<-"tox:excessive ; rep:inadequate"}
        if (eff2 >r1 && tox1+tox2>=ct2) {result2[n]<-"tox:excessive ; rep:OK"}
        if (eff2<=r1 && tox1+tox2 <ct2) {result2[n]<-"tox:OK ; rep:inadequate"}
      }
      
      #durée de l'essai
      if (res_etape2[n]=="echec") {
        if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])==ct1-1) {time[n]<-patient[n2,5]+2*t_obs}
        if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])!=ct1-1) {time[n]<-patient[n2,5]+t_obs}
      }
      
      
      #ETAPE 3 (évaluation de l'efficacité et de la toxicité)
      if (res_etape2[n]=="succes") {
        stage3=patient[(n2+1):N,]
        eff3=sum(stage3[,1])+sum(stage3[,2])
        ech3=sum(stage3[,3])+sum(stage3[,4])
        tox3=sum(stage3[,2])+sum(stage3[,4])
        
        if (eff2+eff3>R && tox1+tox2+tox3<CT) {res_etape3[n]<-"succes"} else {res_etape3[n]<-"echec"}
        #if (pbeta(theta0, a0+tox3+tox2+tox1, b0+(N-tox3-tox2-tox1), lower.tail=F) < tau && eff2+eff3 > R) {res_etape3[n]<-"succes"} else {res_etape3[n]<-"echec"}
        sample3=sample3+N-n2
        nb_eff=nb_eff+eff3
        nb_ech=nb_ech+ech3
        nb_tox=nb_tox+tox3
        
        #raison de l'arrêt
        if (res_etape3[n]=="echec") {
          if (eff2+eff3<=R && tox1+tox2+tox3>=CT) {result3[n]<-"tox:excessive ; rep:inadequate"}
          if (eff2+eff3 >R && tox1+tox2+tox3>=CT) {result3[n]<-"tox:excessive ; rep:OK"}
          if (eff2+eff3<=R && tox1+tox2+tox3 <CT) {result3[n]<-"tox:OK ; rep:inadequate"}
        }
        
        #durée de l'essai 
        if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])==ct1-1) {
          if (sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])==r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])==ct2-1) {time[n]<-patient[N,5]+3*t_obs}
          if (sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])!=r1 && tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])!=ct2-1) {time[n]<-patient[N,5]+2*t_obs}
        }
        if (sum(stage1[1:(n1-1),2])+sum(stage1[1:(n1-1),4])!=ct1-1) {
          if (sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])==r1 | tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])==ct2-1) {time[n]<-patient[N,5]+2*t_obs}
          if (sum(stage2[1:(n2-n1-1),1])+sum(stage2[1:(n2-n1-1),2])!=r1 && tox1+sum(stage2[1:(n2-n1-1),2])+sum(stage2[1:(n2-n1-1),4])!=ct2-1) {time[n]<-patient[N,5]+t_obs}
        }
        
      } 
      
    }
    
  }
  
  #% de succès total
  if (all(is.na(res_etape3[res_etape3=="succes"]))) {
    grid$H11[i] <- 0
  } else {
    grid$H11[i] <- table(res_etape3[res_etape3=="succes"])/100
  }
  print(paste(grid$H11[i],"% de succès"))
  
  #% arrêt précoce
  if (all(is.na(res_etape2[res_etape2=="echec"])) && all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$early_ter[i] <- 0
  } else {
    grid$early_ter[i] <- ((table(res_etape2[res_etape2=="echec"])+table(res_etape1[res_etape1=="echec"]))*100)/ntrials
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
  
  
  #nb d'échecs et de succès à chaque étape 
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
  
  if (all(is.na(res_etape3[res_etape3=="succes"]))) {
    grid$res3_succes[i] <- 0
  } else {
    grid$res3_succes[i] <- table(res_etape3[res_etape3=="succes"])
  }
  if (all(is.na(res_etape3[res_etape3=="echec"]))) {
    grid$res3_echecs[i] <- 0
  } else {
    grid$res3_echecs[i] <- table(res_etape3[res_etape3=="echec"])
  }
  
  
  #raison de l'arrêt à l'étape 1
  if (all(is.na(result1[result1=="tox:excessive"]))) {
    grid$res1_toxE[i] <- 0
  } else {
    grid$res1_toxE[i] <- table(result1[result1=="tox:excessive"])
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
  #raison de l'arrêt à l'étape 3
  if (all(is.na(result3[result3=="tox:excessive ; rep:inadequate"]))) {
    grid$res3_toxE_repI[i] <- 0
  } else {
    grid$res3_toxE_repI[i] <- table(result3[result3=="tox:excessive ; rep:inadequate"])
  }
  if (all(is.na(result3[result3=="tox:excessive ; rep:OK"]))) {
    grid$res3_toxE_repO[i] <- 0
  } else {
    grid$res3_toxE_repO[i] <- table(result3[result3=="tox:excessive ; rep:OK"])
  }
  if (all(is.na(result3[result3=="tox:OK ; rep:inadequate"]))) {
    grid$res3_toxO_repI[i] <- 0
  } else {
    grid$res3_toxO_repI[i] <- table(result3[result3=="tox:OK ; rep:inadequate"])
  }
  
  #% arrêt pour tox
  grid$stop_tox[i] <- grid$res1_toxE[i]+grid$res2_toxE_repI[i]+grid$res2_toxE_repO[i]+grid$res3_toxE_repI[i]+grid$res3_toxE_repO[i]
  
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

writexl::write_xlsx(grid, "C:\\util\\marie\\2022-2023\\stage\\R\\résultats\\Simon+TOX_tpspois2.xlsx")


######### Simon + TOX en continu grid search ##############
# code : monitoring tox en continu 
#grid search



paramA=seq(0.25,0.4,0.05)
paramB=seq(0.85,0.98,0.01)
grid=expand.grid(paramA,paramB)
colnames(grid)=c("theta0","tau")

a0=1
b0=1

peff=0.4
ptox=0.35

pefftox=peff*ptox
peffnotox=peff-pefftox
pnoefftox=ptox-pefftox
pnoeffnotox=(1-peff)-pnoefftox
pefftox+peffnotox+pnoefftox+pnoeffnotox==1

n1=6 ; n2=25 ; N=60
r1=4 ; R=15 #mêmes règles que Simon + safety

tox=notox=eff=ech=0

ntrials=10000
npat=N

res_tox1=res_etape1=res_etape2=rep(NA,ntrials)
stop_trial=rep(NA,ntrials) #raison de l'échec 
monito_tox=rep(NA,ntrials) #nb de patients inclus avant arrêt pour tox
sample=0
nb_tox=nb_ech=nb_eff=nb_notox=0
time=eff_1=rep(NA,ntrials)
t_obs=4



for (g in 1:nrow(grid)) {
  
  print(paste("MATRICE LIGNE",g))
  
  theta0 <- grid$theta0[g]
  tau <- grid$tau[g]
  print(paste("theta0=",theta0,"tau=",tau))
  
  res_tox1=res_etape1=res_etape2=rep(NA,ntrials)
  stop_trial=rep(NA,ntrials) #raison de l'échec 
  monito_tox=rep(NA,ntrials) #nb de patients inclus avant arrêt pour tox
  sample=0
  nb_tox=nb_ech=nb_eff=nb_notox=0
  time=eff_1=rep(NA,ntrials)
  
  
  set.seed(3003)
  for (n in 1:ntrials) {
    
    #print(paste("ESSAI",n))
    tox=notox=eff=ech=0
    
    #simulation des patients et outcomes
    patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))),cumsum(c(0,rep(2,N-1))))
    
    
    for (i in 1:npat) {
      
      eff = eff + patient[,1][i]+patient[,2][i]
      ech = ech + patient[,3][i]+patient[,4][i]
      tox = tox + patient[,2][i]+patient[,4][i]
      notox = notox + patient[,1][i]+patient[,3][i]
      #print(paste("patient n°",i))
      
      if (i == 6) {
        if (pbeta(theta0, a0+tox, b0+notox, lower.tail=F) < tau) {res_tox1[n] <-"succes"} else {res_tox1[n]<-"echec"}
        #print(paste("res AI1:",res_tox1[n]))
        #print(paste("tox",tox,"eff",eff))
      }
      
      if(i == n2-1) {eff_1[n] <- eff}
      
      if (i == n2) {
        if (eff > r1 && pbeta(theta0, a0+tox, b0+notox, lower.tail=F) < tau) {res_etape1[n]<-"succes"} else {
          res_etape1[n]<-"echec"
          if (eff <= r1) {stop_trial[n] <- "arrêt futilité"}
          if (pbeta(theta0, a0+tox, b0+notox, lower.tail=F) >= tau) {stop_trial[n] <- "arrêt tox"}
          if (eff <= r1 && pbeta(theta0, a0+tox, b0+notox, lower.tail=F) >= tau) {stop_trial[n] <- "arrêt tox et futilité"}
        }
        #print(paste("res AI2:",res_etape1[n]))
        #print(paste("tox",tox,"eff",eff))
        if (eff <= r1 | pbeta(theta0, a0+tox, b0+notox, lower.tail=F) >= tau) {break()} 
      }
      
      if (res_etape1[n]=="succes" && i == N) {
        if (eff > R && pbeta(theta0, a0+tox, b0+notox, lower.tail=F) < tau) {res_etape2[n]<-"succes"} else {
          res_etape2[n]<-"echec"
          if (eff <= R) {stop_trial[n] <- "arrêt futilité"}
          if (pbeta(theta0, a0+tox, b0+notox, lower.tail=F) >= tau) {stop_trial[n] <- "arrêt tox"}
          if (eff <= R && pbeta(theta0, a0+tox, b0+notox, lower.tail=F) >= tau) {stop_trial[n] <- "arrêt tox et futilité"}
          }
        #print(paste("res AI3:",res_etape2[n]))
        #print(paste("tox",tox,"eff",eff))
      }
      
      if (pbeta(theta0, a0+tox, b0+notox, lower.tail=F) >= tau) {
        monito_tox[n] <- i
        stop_trial[n] <- "arrêt tox"
      } 
      
      if (pbeta(theta0, a0+tox, b0+notox, lower.tail=F) >= tau) {break()}
    }
    
    nb_eff=nb_eff+eff
    nb_ech=nb_ech+ech
    nb_tox=nb_tox+tox
    nb_notox=nb_notox+notox
    
    sample=sample+i
  
    if (is.na(eff_1[n]) | eff_1[n]!=r1) {time[n]<-patient[i,5]+t_obs
    } else {
      if (eff_1[n]==r1 && is.na(res_etape1[n])) {time[n]<-patient[i,5]+t_obs
      } else {
        if (eff_1[n]==r1 && res_etape1[n]!="succes") {time[n]<-patient[i,5]+t_obs} 
        if (eff_1[n]==r1 && res_etape1[n]=="succes") {time[n]<-patient[i,5]+2*t_obs}
      } 
    }
    #print(paste("durée:",time[n]))
    
  }
  
  #% de succès
  if (all(is.na(res_etape2[res_etape2=="succes"]))) {
    grid$H11[g] <- 0
  } else {
    grid$H11[g] <- table(res_etape2[res_etape2=="succes"])/100
  }
  print(paste(grid$H11[g],"% de succès"))
  
  #effectif attendu
  grid$ess[g] <- sample/ntrials
  
  #nb efficacités
  grid$moy_eff[g] <- nb_eff/ntrials
  
  #nb d'échecs
  grid$moy_ech[g] <- nb_ech/ntrials
  
  #nb de toxicités 
  grid$moy_tox[g] <- nb_tox/ntrials
  
  #durée de l'essai
  grid$moy_time[g] <- mean(time)
  
  #étape 1 : nb d'échecs et succès 
  if (all(is.na(res_etape1[res_etape1=="succes"]))) {
    grid$res1_succes[g] <- 0
  } else {
    grid$res1_succes[g] <- table(res_etape1[res_etape1=="succes"])
  }
  if (all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$res1_echecs[g] <- 0
  } else {
    grid$res1_echecs[g] <- table(res_etape1[res_etape1=="echec"])
  }
  
  #étape 2 : nb d'échecs et succès
  if (all(is.na(res_etape2[res_etape2=="succes"]))) {
    grid$res2_succes[g] <- 0
  } else {
    grid$res2_succes[g] <- table(res_etape2[res_etape2=="succes"])
  }
  if (all(is.na(res_etape2[res_etape2=="echec"]))) {
    grid$res2_echecs[g] <- 0
  } else {
    grid$res2_echecs[g] <- table(res_etape2[res_etape2=="echec"])
  }
  
  #raison de l'arrêt 
  if (all(is.na(stop_trial[stop_trial=="arrêt futilité"]))) {
    grid$toxO_repI[g] <- 0
  } else {
    grid$toxO_repI[g] <- table(stop_trial[stop_trial=="arrêt futilité"])
  }
  
  if (all(is.na(stop_trial[stop_trial=="arrêt tox"]))) {
    grid$toxE_repO[g] <- 0
  } else {
    grid$toxE_repO[g] <- table(stop_trial[stop_trial=="arrêt tox"])
  }
  
  if (all(is.na(stop_trial[stop_trial=="arrêt tox et futilité"]))) {
    grid$toxE_repI[g] <- 0
  } else {
    grid$toxE_repI[g] <- table(stop_trial[stop_trial=="arrêt tox et futilité"])
  }
  
  
}


writexl::write_xlsx(grid, "C:\\util\\marie\\2022-2023\\stage\\R\\résultats Simon + tox\\monitoring continu\\scenario6_H10.xlsx")

table(res_tox1)
table(res_etape1)
table(res_etape2)

sample/ntrials

#raison de l'échec
table(stop_trial)

#durée
table(is.na(time))
mean(time)

#pour 5 essais (peff=0.2 et ptox=0.35)
#essai 1: 32 patients
#essai 2: 25 patients
#essai 3: 25 patients 
#essai 4:  2 patients
#essai 5:  2 patients
#86 patients au total


#pour 5 essais (avec code corrigé)
#essai 1: 21 patients
#essai 2: 25 patients
#essai 3: 11 patients 
#essai 4: 60 patients
#essai 5: 60 patients
#86 patients au total





######### Simon + TOX en continu ############
# code monitoring continu 
#scenarios



param_eff=c(0.2,0.4,0.4,0.2,0.4,0.2,0.4,0.3)
param_tox=c(0.35,0.15,0.05,0.05,0.45,0.15,0.35,0.2)
scenario=paste("Scenario",1:8)
grid=cbind(data.frame(scenario,param_eff,param_tox))
colnames(grid)=c("Scenario","Peff","Ptox")

a0=1
b0=1

theta0=0.35
tau=0.95

n1=6 ; n2=25 ; N=60
r1=4 ; R=15 #mêmes règles que Simon + safety

tox=notox=eff=ech=0

ntrials=10000
npat=N


for (g in 1:nrow(grid)) {
  
  print(paste("MATRICE LIGNE",g))
  
  peff <- grid$Peff[g]
  ptox <- grid$Ptox[g]
  print(paste("peff=",peff,"ptox=",ptox))
  
  pefftox=peff*ptox
  peffnotox=peff-pefftox
  pnoefftox=ptox-pefftox
  pnoeffnotox=(1-peff)-pnoefftox
  
  res_tox1=res_etape1=res_etape2=rep(NA,ntrials)
  stop_trial=rep(NA,ntrials) #raison de l'échec 
  monito_tox=rep(NA,ntrials) #nb de patients inclus avant arrêt pour tox
  sample=0
  nb_tox=nb_ech=nb_eff=nb_notox=0
  time=eff_1=rep(NA,ntrials)
  t_obs=4
  
  
  set.seed(3003)
  for (n in 1:ntrials) {
    
    #print(paste("ESSAI",n))
    tox=notox=eff=ech=0
    
    #simulation des patients et outcomes
    patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))),cumsum(c(0,rep(2,N-1))))
    
    
    for (i in 1:npat) {
      
      eff = eff + patient[,1][i]+patient[,2][i]
      ech = ech + patient[,3][i]+patient[,4][i]
      tox = tox + patient[,2][i]+patient[,4][i]
      notox = notox + patient[,1][i]+patient[,3][i]
      #print(paste("patient n°",i))
      
      if (i == 6) {
        if (pbeta(theta0, a0+tox, b0+notox, lower.tail=F) < tau) {res_tox1[n] <-"succes"} else {res_tox1[n]<-"echec"}
        #print(paste("res AI1:",res_tox1[n]))
        #print(paste("tox",tox,"eff",eff))
      }
      
      if(i == n2-1) {eff_1[n] <- eff}
      
      if (i == n2) {
        if (eff > r1 && pbeta(theta0, a0+tox, b0+notox, lower.tail=F) < tau) {res_etape1[n]<-"succes"} else {
          res_etape1[n]<-"echec"
          if (eff <= r1) {stop_trial[n] <- "arrêt futilité"}
          if (pbeta(theta0, a0+tox, b0+notox, lower.tail=F) >= tau) {stop_trial[n] <- "arrêt tox"}
          if (eff <= r1 && pbeta(theta0, a0+tox, b0+notox, lower.tail=F) >= tau) {stop_trial[n] <- "arrêt tox et futilité"}
        }
        #print(paste("res AI2:",res_etape1[n]))
        #print(paste("tox",tox,"eff",eff))
        if (eff <= r1 | pbeta(theta0, a0+tox, b0+notox, lower.tail=F) >= tau) {break()} 
      }
      
      if (res_etape1[n]=="succes" && i == N) {
        if (eff > R && pbeta(theta0, a0+tox, b0+notox, lower.tail=F) < tau) {res_etape2[n]<-"succes"} else {
          res_etape2[n]<-"echec"
          if (eff <= R) {stop_trial[n] <- "arrêt futilité"}
          if (pbeta(theta0, a0+tox, b0+notox, lower.tail=F) >= tau) {stop_trial[n] <- "arrêt tox"}
          if (eff <= R && pbeta(theta0, a0+tox, b0+notox, lower.tail=F) >= tau) {stop_trial[n] <- "arrêt tox et futilité"}
        }
        #print(paste("res AI3:",res_etape2[n]))
        #print(paste("tox",tox,"eff",eff))
      }
      
      if (pbeta(theta0, a0+tox, b0+notox, lower.tail=F) >= tau) {
        monito_tox[n] <- i
        stop_trial[n] <- "arrêt tox"
      } 
      
      if (pbeta(theta0, a0+tox, b0+notox, lower.tail=F) >= tau) {break()}
    }
    
    nb_eff=nb_eff+eff
    nb_ech=nb_ech+ech
    nb_tox=nb_tox+tox
    nb_notox=nb_notox+notox
    
    sample=sample+i
    
    if (is.na(eff_1[n]) | eff_1[n]!=r1) {time[n]<-patient[i,5]+t_obs
    } else {
      if (eff_1[n]==r1 && is.na(res_etape1[n])) {time[n]<-patient[i,5]+t_obs
      } else {
        if (eff_1[n]==r1 && res_etape1[n]!="succes") {time[n]<-patient[i,5]+t_obs} 
        if (eff_1[n]==r1 && res_etape1[n]=="succes") {time[n]<-patient[i,5]+2*t_obs}
      } 
    }
    #print(paste("durée:",time[n]))
    
  }
  
  #% de succès
  if (all(is.na(res_etape2[res_etape2=="succes"]))) {
    grid$H11[g] <- 0
  } else {
    grid$H11[g] <- table(res_etape2[res_etape2=="succes"])/100
  }
  print(paste(grid$H11[g],"% de succès"))
  
  #effectif attendu
  grid$ess[g] <- sample/ntrials
  
  #nb efficacités
  grid$moy_eff[g] <- nb_eff/ntrials
  
  #nb d'échecs
  grid$moy_ech[g] <- nb_ech/ntrials
  
  #nb de toxicités 
  grid$moy_tox[g] <- nb_tox/ntrials
  
  #durée de l'essai
  grid$moy_time[g] <- mean(time)
  
  #étape 1 : nb d'échecs et succès 
  if (all(is.na(res_etape1[res_etape1=="succes"]))) {
    grid$res1_succes[g] <- 0
  } else {
    grid$res1_succes[g] <- table(res_etape1[res_etape1=="succes"])
  }
  if (all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$res1_echecs[g] <- 0
  } else {
    grid$res1_echecs[g] <- table(res_etape1[res_etape1=="echec"])
  }
  
  #étape 2 : nb d'échecs et succès
  if (all(is.na(res_etape2[res_etape2=="succes"]))) {
    grid$res2_succes[g] <- 0
  } else {
    grid$res2_succes[g] <- table(res_etape2[res_etape2=="succes"])
  }
  if (all(is.na(res_etape2[res_etape2=="echec"]))) {
    grid$res2_echecs[g] <- 0
  } else {
    grid$res2_echecs[g] <- table(res_etape2[res_etape2=="echec"])
  }
  
  #raison de l'arrêt 
  if (all(is.na(stop_trial[stop_trial=="arrêt futilité"]))) {
    grid$toxO_repI[g] <- 0
  } else {
    grid$toxO_repI[g] <- table(stop_trial[stop_trial=="arrêt futilité"])
  }
  
  if (all(is.na(stop_trial[stop_trial=="arrêt tox"]))) {
    grid$toxE_repO[g] <- 0
  } else {
    grid$toxE_repO[g] <- table(stop_trial[stop_trial=="arrêt tox"])
  }
  
  if (all(is.na(stop_trial[stop_trial=="arrêt tox et futilité"]))) {
    grid$toxE_repI[g] <- 0
  } else {
    grid$toxE_repI[g] <- table(stop_trial[stop_trial=="arrêt tox et futilité"])
  }
  
  
}

writexl::write_xlsx(grid, "C:\\util\\marie\\2022-2023\\stage\\R\\résultats\\Simon+TOX continu.xlsx")




######### recherche des valeurs seuils #################
# Simon + tox 
# déterminer les valeurs seuils 

#1ère analyse intermédiaire : 6 patients

a0=b0=1
theta0=0.2
tau=0.9

tox=0:6
pbeta(theta0, a0+tox, b0+(6-tox), lower.tail=F) > tau


#2ème analyse intermédiare
tox2=0:25

pbeta(theta0, a0+tox2, b0+(25-tox2), lower.tail=F) > tau


#3ème analyse intermédiaire
tox3=0:60

pbeta(theta0, a0+tox3, b0+(60-tox3), lower.tail=F) > tau


#fonction 
a0=b0=1
n=6 #ou 25 ou 60
tau=0.95
theta0=0.3

seuil <- function(theta0,tau,a0,b0,n) {
  
  tox <- 1:n
  
  for (i in tox) {
    if (pbeta(theta0, a0+i, b0+(n-i), lower.tail=F) > tau) {
      return(i)
    }
  }
  return(NA)
}

seuil(theta0 = 0.2, tau = 0.9, a0=1, b0=1, n=6)
#1ère AI : 3
#2ème AI : 8
#3ème AI : 16


#grid search des seuils (monito séquentiel)

paramA=seq(0.25,0.4,0.05)
paramB=seq(0.85,0.98,0.01)
grid=expand.grid(paramA,paramB)
colnames(grid)=c("theta0","tau")

a0=b0=1

n1=6 ; n2=25 ; N=60


for (i in 1:nrow(grid)) {
  
  theta0 <- grid$theta0[i]
  tau <- grid$tau[i]
  
  #1ère analyse intermédiaire
  grid$ai1[i] <- seuil(theta0 = theta0, tau = tau, a0=1, b0=1, n=n1)
  
  #2ème analyse intermédiaire
  grid$ai2[i] <- seuil(theta0 = theta0, tau = tau, a0=1, b0=1, n=n2)
  
  #3ème analyse intermédiaire
  grid$ai3[i] <- seuil(theta0 = theta0, tau = tau, a0=1, b0=1, n=N)
  
}

writexl::write_xlsx(grid, "C:\\util\\marie\\2022-2023\\stage\\R\\résultats Simon + tox\\grille_seuils.xlsx")



#test avec résultats article (Ivanova)
a0=1=b0=1
theta0=0.2
tau=0.98

n=1:60

result=c()



n=1:60
for (i in n) {
  
  tox <- 0:n

  for (i in tox) {
    if (pbeta(theta0, a0+i, b0+(n-i), lower.tail=F) > tau) {
      return(i)
    }
  }
  return(NA)
}
