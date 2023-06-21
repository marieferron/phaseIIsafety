# SAFETY LEAD-IN + SIMON

############## Safety (1 dose, 2 sem) + Simon modif temps AI #########

#alpha = 5%
r1 = 4 ; n1 = 19 ; R = 15 ; N = 54

#alpha = 2,5%
#r1 = 5 ; n1 = 22 ; R = 19 ; N = 66

ntrials=10000
npat=60
n_sl=6

res_safety=res_etape1=res_etape2=rep(NA,ntrials) #nombre d'échecs et succès pour chaque essai
i=0
sample_safety=sample_etape1=sample_etape2=0
nb_eff=0 #nombre de réponses (succès) à chaque étape
nb_ech=0 #nombre de non-réponses (échecs) à chaque étape
nb_tox=0 #nombre de toxicités à chaque étape
time=rep(NA,ntrials)
t_obs=4 #4 semaines (28 jours)

#param_eff=c(0.2,0.4,0.4,0.2,0.4,0.2,0.4,0.3)
#param_tox=c(0.35,0.15,0.05,0.05,0.45,0.15,0.35,0.2)
param_eff=c(0.2,0.2,0.4,0.4,0.4,0.2,0.4,0.3)
param_tox=c(0.35,0.15,0.35,0.15,0.05,0.05,0.45,0.25)
scenario=paste("Scenario",1:8)
grid=cbind(data.frame(scenario,param_eff,param_tox))
colnames(grid)=c("Scenario","Peff","Ptox")


for (g in 1:nrow(grid)) {
  
  print(paste("Scenario",g))
  
  peff<-grid$Peff[g]
  ptox<-grid$Ptox[g]
  print(paste("peff=",peff,"ptox=",ptox))
  
  pefftox=peff*ptox
  peffnotox=peff-pefftox
  pnoefftox=ptox-pefftox
  pnoeffnotox=(1-peff)-pnoefftox
  
  res_safety=res_etape1=res_etape2=rep(NA,ntrials) #nombre d'échecs et succès pour chaque essai
  i=0
  sample=sample_s=sample1=sample2=0
  nb_eff=0 #nombre de réponses (succès) à chaque étape
  nb_ech=0 #nombre de non-réponses (échecs) à chaque étape
  nb_tox=0 #nombre de toxicités à chaque étape
  time=rep(NA,ntrials)
  check_stop1=check_stop2=rep(NA,ntrials)
  p_eff=p_ech=p_tox=rep(NA,ntrials)
  
  #SAFETY LEAD-IN (1 DOSE) + SIMON
  set.seed(3003)
  for (n in 1:ntrials) {
    
    #Simulation des patients
    #patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))),cumsum(c(0,rpois(N+n_sl-1,2)))) 
    patient=t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox)))
    pois=c(0,rpois(npat-1,2))
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
    
    limite=0
    rep_safe=rep_echec=0
    
    nbeff1=nbeff2=nbeff3=0
    nbech1=nbech2=nbech3=0
    nbtox1=nbtox2=nbtox3=0
    sample_s=sample1=sample2=0
    
    #SAFETY LEAD-IN
    for (i in 1:n_sl) {
      limite=limite + patient[,2][i] + patient[,4][i]
      rep_safe=rep_safe + patient[,1][i] + patient[,2][i]
      rep_echec=rep_echec + patient[,3][i] + patient[,4][i]
      if (limite==2) {break()}
    }
    if (limite==2) {res_safety[n]<-"echec"} else {res_safety[n]<-"succes_s"} 
    if (res_safety[n]=="echec") {
      sample=sample+i+1
      sample_s=i+1
      nb_eff=nb_eff+rep_safe+patient[,1][i+1] + patient[,2][i+1]
      nb_ech=nb_ech+rep_echec+patient[,3][i+1] + patient[,4][i+1]
      nb_tox=nb_tox+limite+patient[,2][i+1] + patient[,4][i+1]
      nbeff1=rep_safe+patient[,1][i+1] + patient[,2][i+1]
      nbech1=rep_echec+patient[,3][i+1] + patient[,4][i+1]
      nbtox1=limite+patient[,2][i+1] + patient[,4][i+1]
    } else {
      sample=sample+i
      sample_s=i
      nb_eff=nb_eff+rep_safe
      nb_ech=nb_ech+rep_echec
      nb_tox=nb_tox+limite
      nbeff1=rep_safe
      nbech1=rep_echec
      nbtox1=limite
    }
    
    #print(paste("essai",n))
    #print(paste(i,"patients",limite,"toxicity",rep_safe,"efficacités"))
    #print(paste("Résultat s_l:",res_safety[n]))
    
    if (res_safety[n]=="echec") {
      time[n]<-patient[i+1,5]+t_obs
      #cat("Durée de l'essai(saf):",time[n])
    }
    
    #vérification inclusion AI
    if (patient[n_sl,5]==patient[n_sl+1,5]) {check_stop1[n] <- "STOP"} else {"OK"}
    if (patient[n_sl+n1,5]==patient[n_sl+n1+1,5]) {check_stop2[n] <- "STOP"} else {"OK"}
    
    #SIMON ETAPE 1
    if (res_safety[n]=="succes_s") {
      
      eff1=sum(patient[n_sl+1:n1,1]) + sum(patient[n_sl+1:n1,2])
      ech1=sum(patient[n_sl+1:n1,3]) + sum(patient[n_sl+1:n1,4])
      tox1=sum(patient[n_sl+1:n1,2]) + sum(patient[n_sl+1:n1,4])
      
      if (eff1>r1) {res_etape1[n]<-"succes1"} else {res_etape1[n]<-"echec"} #0=échec ; 1=succès (passage à l'étape 2)
      sample=sample+n1
      sample1=n1
      nb_eff=nb_eff+eff1
      nb_ech=nb_ech+ech1
      nb_tox=nb_tox+tox1
      nbeff2=eff1
      nbech2=ech1
      nbtox2=tox1
      #print(paste(n1,"patients",eff1,"efficacités",tox1,"toxicités"))
      #print(paste("Résultat étape1:",res_etape1[n]))
      
      if (res_etape1[n]=="echec") {
        if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) == limite-1) {time[n]<-patient[n1+n_sl,5]+2*t_obs}
        else {time[n]<-patient[n1+n_sl,5]+t_obs}
        #cat("Durée de l'essai(1):",time[n])
      }
      
      
      #SIMON ETAPE 2
      if (res_etape1[n]=="succes1") { 
        
        eff2=sum(patient[n_sl+n1+1:(N-n1),1]) + sum(patient[n_sl+n1+1:(N-n1),2])
        ech2=sum(patient[n_sl+n1+1:(N-n1),3]) + sum(patient[n_sl+n1+1:(N-n1),4])
        tox2=sum(patient[n_sl+n1+1:(N-n1),2]) + sum(patient[n_sl+n1+1:(N-n1),4])
        
        if (eff1+eff2>R) {res_etape2[n]<-"succes2"} else {res_etape2[n]<-"echec"} #0=échec ; 1=succès (conclu d'efficacité)
        sample=sample+N-n1
        sample2=N-n1
        nb_eff=nb_eff+eff2
        nb_ech=nb_ech+ech2
        nb_tox=nb_tox+tox2
        nbeff3=eff2
        nbech3=ech2
        nbtox3=tox2
        #print(paste(N-n1,"patients",eff2,"efficacités",tox2,"toxicités"))
        #print(paste("Résultat étape2:",res_etape2[n]))
        
        
        if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) == limite-1) {
          if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) == r1) {time[n] <- patient[N+n_sl,5] + 3*t_obs}
          if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) != r1) {time[n] <- patient[N+n_sl,5] + 2*t_obs}
        }
        if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) != limite-1) {
          if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) == r1) {time[n] <- patient[N+n_sl,5] + 2*t_obs}
          if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) != r1) {time[n] <- patient[N+n_sl,5] + t_obs}
        }
        #cat("Durée de l'essai(2):",time[n])
        
        
      }
      
    }
    
    p_eff[n] <- (nbeff1+nbeff2+nbeff3) / (sample_s+sample1+sample2) 
    p_ech[n] <- (nbech1+nbech2+nbech3) / (sample_s+sample1+sample2) 
    p_tox[n] <- (nbtox1+nbtox2+nbtox3) / (sample_s+sample1+sample2)
    
  }
  
  #effectif attendu
  grid$ess[g] <- sample/ntrials
  
  #nb d'efficacités 
  grid$moy_eff[g] <- nb_eff/ntrials
  grid$prob_eff[g] <- nb_eff/sample
  grid$prob_eff2[g] <- mean(p_eff)
  
  #nb d'échecs
  grid$moy_ech[g] <- nb_ech/ntrials
  grid$prob_ech[g] <- nb_ech/sample
  grid$prob_ech2[g] <- mean(p_ech)
  
  #nb de toxicités 
  grid$moy_tox[g] <- nb_tox/ntrials
  grid$prob_tox[g] <- nb_tox/sample
  grid$prob_tox2[g] <- mean(p_tox)
  
  #durée de l'essai
  grid$moy_time[g] <- mean(time)
  
  #% conclu d'efficacité
  if (all(is.na(res_etape2[res_etape2=="succes2"]))) {
    grid$succes[g] <- 0
  } else {
    grid$succes[g] <- ((table(res_etape2[res_etape2=="succes2"]))*100)/ntrials
  }
  
  #% arrêt précoce
  if (all(is.na(res_safety[res_safety=="echec"])) && all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$early_ter[g] <- 0
  } else {
    grid$early_ter[g] <- ((table(res_safety[res_safety=="echec"])+table(res_etape1[res_etape1=="echec"]))*100)/ntrials
  }
  
  #% arrêt pour tox
  if (all(is.na(res_safety[res_safety=="echec"]))) {
    grid$stop_tox[g] <- 0
  } else {
    grid$stop_tox[g] <- (table(res_safety[res_safety=="echec"])*100)/ntrials
  }
  
  
  #nb d'échecs et de succès à chaque étape 
  if (all(is.na(res_safety[res_safety=="succes_s"]))) {
    grid$ressafe_succes[g] <- 0
  } else {
    grid$ressafe_succes[g] <- table(res_safety[res_safety=="succes_s"])
  }
  if (all(is.na(res_safety[res_safety=="echec"]))) {
    grid$ressafe_echecs[g] <- 0
  } else {
    grid$ressafe_echecs[g] <- table(res_safety[res_safety=="echec"])
  }
  
  if (all(is.na(res_etape1[res_etape1=="succes1"]))) {
    grid$res1_succes[g] <- 0
  } else {
    grid$res1_succes[g] <- table(res_etape1[res_etape1=="succes1"])
  }
  if (all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$res1_echecs[g] <- 0
  } else {
    grid$res1_echecs[g] <- table(res_etape1[res_etape1=="echec"])
  }
  
  if (all(is.na(res_etape2[res_etape2=="succes2"]))) {
    grid$res2_succes[g] <- 0
  } else {
    grid$res2_succes[g] <- table(res_etape2[res_etape2=="succes2"])
  }
  if (all(is.na(res_etape2[res_etape2=="echec"]))) {
    grid$res2_echecs[g] <- 0
  } else {
    grid$res2_echecs[g] <- table(res_etape2[res_etape2=="echec"])
  }
  
  #vérification inclusion AI
  if (all(is.na(check_stop1[check_stop1=="STOP"]))) {
    grid$check_AI1[g] <- 0
  } else {
    grid$check_AI1[g] <- table(check_stop1[check_stop1=="STOP"])
  }
  if (all(is.na(check_stop2[check_stop2=="STOP"]))) {
    grid$check_AI2[g] <- 0
  } else {
    grid$check_AI2[g] <- table(check_stop2[check_stop2=="STOP"])
  }
  
}



writexl::write_xlsx(grid,"C:\\util\\marie\\2022-2023\\stage\\R\\résultats\\safety1dose_Simon_tpspoisAI+prob.xlsx")





############## Safety (2 doses, 2 sem) + Simon modif temps AI #########

ntrials=1000
npat=60

#alpha = 5%
n_sl=6 ; r1 = 4 ; n1 = 19 ; R = 15 ; N = 54

res_safety=res_etape1=res_etape2=rep(NA,ntrials) #nombre d'échecs et succès pour chaque essai
i=l=0
sample_s=sample1_d1=sample2_d1=0
sample_s2=sample1_d2=sample2_d2=0
nb_eff=0 #nombre de réponses (succès) à chaque étape
nb_ech=0 #nombre de non-réponses (échecs) à chaque étape
nb_tox=0 #nombre de toxicités à chaque étape
time=rep(NA,ntrials)
t_obs=4 #mesure de l'outcome après 28j (4 semaines)
timesafe=rep(NA,ntrials) #stocke la durée du 1er safety lead-in (utile si désescalade de dose)

# param_eff=c(0.2,0.4,0.2,0.4,0.4,0.2,0.2,0.3)
# param_tox=c(0.2,0.2,0.05,0.05,0.35,0.35,0.1,0.2)
# param_eff2=c(0.4,0.5,0.4,0.5,0.5,0.4,0.3,0.4)
# param_tox2=c(0.35,0.35,0.2,0.2,0.45,0.45,0.2,0.3)
param_eff=c(0.2, 0.2, 0.4, 0.4, 0.2, 0.4, 0.2, 0.1, 0.4, 0.4, 0.2, 0.3, 0.1)
param_tox=c(0.35, 0.05, 0.15, 0.05, 0.35, 0.45, 0.15, 0.05, 0.35, 0.15, 0.15, 0.2, 0.15)
param_eff2=c(0.4, 0.4, 0.5, 0.5, 0.3, 0.6, 0.4, 0.2, 0.5, 0.6, 0.3, 0.4, 0.2)
param_tox2=c(0.5, 0.15, 0.35, 0.15, 0.45, 0.5, 0.35, 0.15, 0.45, 0.45, 0.2, 0.35, 0.35)
scenario=paste("Scenario",9:21)
grid=cbind(data.frame(scenario,param_eff,param_tox,param_eff2,param_tox2))
colnames(grid)=c("Scenario","Peff","Ptox","Peff2","Ptox2")



for (g in 1:nrow(grid)) {
  print(paste("Scenario",g))
  
  #pour DL1
  peff <- grid$Peff2[g]
  ptox <- grid$Ptox2[g]
  
  pefftox=peff*ptox
  peffnotox=peff-pefftox
  pnoefftox=ptox-pefftox
  pnoeffnotox=(1-peff)-pnoefftox
  
  #pour DL-1
  peff2 <- grid$Peff[g]
  ptox2 <- grid$Ptox[g]
  
  pefftox2=peff2*ptox2 
  peffnotox2=peff2-pefftox2
  pnoefftox2=ptox2-pefftox2
  pnoeffnotox2=(1-peff2)-pnoefftox2
  
  
  res_safety=res_etape1=res_etape2=rep(NA,ntrials) #nombre d'échecs et succès pour chaque essai
  i=l=0
  sample_s=sample1_d1=sample2_d1=0
  sample_s2=sample1_d2=sample2_d2=0
  nb_eff=0 #nombre de réponses (succès) à chaque étape
  nb_ech=0 #nombre de non-réponses (échecs) à chaque étape
  nb_tox=0 #nombre de toxicités à chaque étape
  time=rep(NA,ntrials)
  timesafe=rep(NA,ntrials) #stocke la durée du 1er safety lead-in (utile quand désescalade de dose)
  
  
  set.seed(3003)
  for (n in 1:ntrials) {
    #patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))),cumsum(c(0,rep(2,n_sl+N-1)))) #temps fixe: cumsum(c(0,rep(2,n_sl+N-1)
    patient=t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox)))
    pois=c(0,rpois(npat-1,2))
    patient_7=pois[7]
    patient_26=pois[26]
    #patient_29=pois[29]
    
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
    # if (patient_29==0) {
    #   while(patient_29==0) {
    #     patient_29<-rpois(1,2)
    #   }
    #   pois[29]<-patient_29
    # }
    patient=cbind(patient,cumsum(pois))
    
    limite=limite2=0
    rep_safe=rep_echec=0
    
    #SAFETY LEAD-IN (DL1)
    for (i in 1:n_sl) {
      limite=limite + patient[,2][i] + patient[,4][i]
      rep_safe=rep_safe + patient[,1][i] + patient[,2][i]  
      rep_echec=rep_echec + patient[,3][i] + patient[,4][i]
      if (limite==2) {break()}
    }
    if (limite==2) {res_safety[n]<-"echec"} else {res_safety[n]<-"succes DL1"} #0=échec ; 1=succès (passage à l'étape 1)
    if (res_safety[n]=="echec") {
      sample_s=sample_s+i+1
      nb_eff=nb_eff+rep_safe+patient[,1][i+1] + patient[,2][i+1]
      nb_ech=nb_ech+rep_echec+patient[,3][i+1] + patient[,4][i+1]
      nb_tox=nb_tox+limite+patient[,2][i+1] + patient[,4][i+1]
    } else {
      sample_s=sample_s+i
      nb_eff=nb_eff+rep_safe
      nb_ech=nb_ech+rep_echec
      nb_tox=nb_tox+limite
    }
    #print(paste("essai",n))
    #print(paste(i,"patients",limite,"toxicity",rep_safe,"efficacités"))
    #print(paste("Résultat s_l:",res_safety[n]))
    
    if (res_safety[n]=="echec") {
      time[n]<-patient[i+1,5]+t_obs
      timesafe[n]<-patient[i+1,5]+t_obs
      #cat("Durée de l'essai(saf):",time[n])
    }
    
    #SIMON ETAPE 1 (DL1)
    if (res_safety[n]=="succes DL1") {
      eff1= sum(patient[n_sl+1:n1,1]) + sum(patient[n_sl+1:n1,2])
      ech1= sum(patient[n_sl+1:n1,3]) + sum(patient[n_sl+1:n1,4])
      tox1= sum(patient[n_sl+1:n1,2]) + sum(patient[n_sl+1:n1,4])
      
      if (eff1>r1) {res_etape1[n]<-"succes1 DL1"} else {res_etape1[n]<-"echec"} #0=échec ; 1=succès (passage à l'étape 2)
      sample1_d1=sample1_d1+n1
      nb_eff=nb_eff+eff1
      nb_ech=nb_ech+ech1
      nb_tox=nb_tox+tox1
      #print(paste(n1,"patients",eff1,"efficacités",tox1,"toxicités"))
      #print(paste("Résultat étape1:",res_etape1[n]))
      
      if (res_etape1[n]=="echec") {
        if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) == 1) {time[n]<-patient[n1+n_sl,5]+2*t_obs}
        if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) != 1) {time[n]<-patient[n1+n_sl,5]+t_obs}
        #cat("Durée de l'essai(1):",time[n])
      }
      
      
      #SIMON ETAPE 2 (DL1)
      if (res_etape1[n]=="succes1 DL1") {
        
        eff2=sum(patient[n_sl+n1+1:(N-n1),1]) + sum(patient[n_sl+n1+1:(N-n1),2])
        ech2=sum(patient[n_sl+n1+1:(N-n1),3]) + sum(patient[n_sl+n1+1:(N-n1),4])
        tox2=sum(patient[n_sl+n1+1:(N-n1),2]) + sum(patient[n_sl+n1+1:(N-n1),4])
        
        if (eff1+eff2>R) {res_etape2[n]<-"succes2 DL1"} else {res_etape2[n]<-"echec"} #0=échec ; 1=succès (conclu d'efficacité)
        sample2_d1=sample2_d1+N-n1
        nb_eff=nb_eff+eff2
        nb_ech=nb_ech+ech2
        nb_tox=nb_tox+tox2
        #print(paste(N-n1,"patients",eff2,"efficacités",tox2,"toxicités"))
        #print(paste("Résultat étape2:",res_etape2[n]))
        
        if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) == 1) {
          if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) == r1) {time[n] <- patient[N+n_sl,5] + 3*t_obs}
          if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) != r1) {time[n] <- patient[N+n_sl,5] + 2*t_obs}
        }
        if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) != 1) {
          if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) == r1) {time[n] <- patient[N+n_sl,5] + 2*t_obs}
          if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) != r1) {time[n] <- patient[N+n_sl,5] + t_obs}
        }
      }
      
    }
    
    if (res_safety[n]=="echec") {
      
      #SAFETY LEAD-IN (DL-1)
      #patient2=cbind(t(rmultinom(npat,1,c(peffnotox2,pefftox2,pnoeffnotox2,pnoefftox2))),cumsum(c(0,rep(2,n_sl+N-1)))) #temps fixe:cumsum(c(0,rep(2,n_sl+N-1)
      patient2=t(rmultinom(npat,1,c(peffnotox2,pefftox2,pnoeffnotox2,pnoefftox2)))
      pois2=c(0,rpois(npat-1,2))
      patient_7_2=pois2[7]
      patient_26_2=pois2[26]
      #patient_29_2=pois2[29]
      
      if (patient_7_2==0) {
        while(patient_7_2==0) {
          patient_7_2<-rpois(1,2)
        }
        pois2[7]<-patient_7_2
      }
      if (patient_26_2==0) {
        while(patient_26_2==0) {
          patient_26_2<-rpois(1,2)
        }
        pois2[26]<-patient_26_2
      }
      # if (patient_29_2==0) {
      #   while(patient_29_2==0) {
      #     patient_29_2<-rpois(1,2)
      #   }
      #   pois2[29]<-patient_29_2
      # }
      patient2=cbind(patient2,cumsum(pois2))
      
      limite2=0
      rep_safe=rep_echec=0
      
      for (l in 1:n_sl) {
        limite2 = limite2 + patient2[,2][l] + patient2[,4][l]
        rep_safe=rep_safe + patient2[,1][l] + patient2[,2][l]  
        rep_echec=rep_echec + patient2[,3][l] + patient2[,4][l]
        if (limite2==2) {break()}
      }
      if (limite2==2) {res_safety[n]<-"echec"} else {res_safety[n]<-"succes DL-1"}
      if (res_safety[n]=="echec") {
        sample_s2=sample_s2+l+1
        nb_eff=nb_eff+rep_safe+patient2[,1][l+1] + patient2[,2][l+1]
        nb_ech=nb_ech+rep_echec+patient2[,3][l+1] + patient2[,4][l+1]
        nb_tox=nb_tox+limite2+patient2[,2][l+1] + patient2[,4][l+1]
      } else {
        sample_s2=sample_s2+l
        nb_eff=nb_eff+rep_safe
        nb_ech=nb_ech+rep_echec
        nb_tox=nb_tox+limite2
      }
      #print(paste("essai",n))
      #print(paste(l,"patients",limite,"toxicity",rep_safe,"efficacités"))
      #print(paste("Résultat s_l:",res_safety[n]))
      
      if (res_safety[n]=="echec") {
        time[n]<-timesafe[n] + patient2[l+1,5]+t_obs
        #cat("Durée de l'essai(saf):",time[n])
      }
      
      #SIMON ETAPE 1 (DL-1)
      if (res_safety[n]=="succes DL-1") {
        eff1=sum(patient2[n_sl+1:n1,1]) + sum(patient2[n_sl+1:n1,2])
        ech1=sum(patient2[n_sl+1:n1,3]) + sum(patient2[n_sl+1:n1,4])
        tox1=sum(patient2[n_sl+1:n1,2]) + sum(patient2[n_sl+1:n1,4])
        
        if (eff1>r1) {res_etape1[n]<-"succes1 DL-1"} else {res_etape1[n]<-"echec"} #0=échec ; 1=succès (passage à l'étape 2)
        sample1_d2=sample1_d2+n1
        nb_eff=nb_eff+eff1
        nb_ech=nb_ech+ech1
        nb_tox=nb_tox+tox1
        #print(paste(n1,"patients",eff1,"efficacités",tox1,"toxicités"))
        #print(paste("Résultat étape1:",res_etape1[n]))
        
        if (res_etape1[n]=="echec") {
          if (sum(patient2[1:(n_sl-1),2])+sum(patient2[1:(n_sl-1),4]) == 1) {time[n]<-timesafe[n]+patient2[n1+n_sl,5]+2*t_obs}
          if (sum(patient2[1:(n_sl-1),2])+sum(patient2[1:(n_sl-1),4]) != 1) {time[n]<-timesafe[n]+patient2[n1+n_sl,5]+t_obs}
          #cat("Durée de l'essai(1):",time[n])
        }
        
        
        #SIMON ETAPE 2 (DL-1)
        if (res_etape1[n]=="succes1 DL-1") {
          eff2=sum(patient2[n_sl+n1+1:(N-n1),1]) + sum(patient2[n_sl+n1+1:(N-n1),2])
          ech2=sum(patient2[n_sl+n1+1:(N-n1),3]) + sum(patient2[n_sl+n1+1:(N-n1),4])
          tox2=sum(patient2[n_sl+n1+1:(N-n1),2]) + sum(patient2[n_sl+n1+1:(N-n1),4])
          
          if (eff1+eff2>R) {res_etape2[n]<-"succes2 DL-1"} else {res_etape2[n]<-"echec"} #0=échec ; 1=succès (conclu d'efficacité)
          sample2_d2=sample2_d2+N-n1
          nb_eff=nb_eff+eff2
          nb_ech=nb_ech+ech2
          nb_tox=nb_tox+tox2
          #print(paste(N-n1,"patients",eff2,"efficacités",tox2,"toxicités"))
          #print(paste("Résultat étape2:",res_etape2[n]))
          
          if (sum(patient2[1:(n_sl-1),2])+sum(patient2[1:(n_sl-1),4]) == 1) {
            if (sum(patient2[n_sl+1:(n1-1),1])+sum(patient2[n_sl+1:(n1-1),2]) == r1) {time[n] <- timesafe[n]+patient2[N+n_sl,5] + 3*t_obs}
            if (sum(patient2[n_sl+1:(n1-1),1])+sum(patient2[n_sl+1:(n1-1),2]) != r1) {time[n] <- timesafe[n]+patient2[N+n_sl,5] + 2*t_obs}
          }
          if (sum(patient2[1:(n_sl-1),2])+sum(patient2[1:(n_sl-1),4]) != 1) {
            if (sum(patient2[n_sl+1:(n1-1),1])+sum(patient2[n_sl+1:(n1-1),2]) == r1) {time[n] <- timesafe[n]+patient2[N+n_sl,5] + 2*t_obs}
            if (sum(patient2[n_sl+1:(n1-1),1])+sum(patient2[n_sl+1:(n1-1),2]) != r1) {time[n] <- timesafe[n]+patient2[N+n_sl,5] + t_obs}
          }
          
        }
      }
    }
  }
  
  
  #effectif attendu
  grid$ess[g] <- (sample_s+sample_s2+sample1_d1+sample2_d1+sample1_d2+sample2_d2)/ntrials
  
  #nb d'efficacités 
  grid$moy_eff[g] <- nb_eff/ntrials
  
  #nb d'échecs
  grid$moy_ech[g] <- nb_ech/ntrials
  
  #nb de toxicités 
  grid$moy_tox[g] <- nb_tox/ntrials
  
  #durée de l'essai
  grid$moy_time[g] <- mean(time)
  
  
  #% conclu d'efficacité
  if (all(is.na(res_etape2[res_etape2=="succes2 DL-1"]))) {
    grid$succes[g] <- ((table(res_etape2[res_etape2=="succes2 DL1"]))*100)/ntrials
  } else {
    if (all(is.na(res_etape2[res_etape2=="succes2 DL1"]))) {
      grid$succes[g] <- ((table(res_etape2[res_etape2=="succes2 DL-1"]))*100)/ntrials
    } else {
      if (all(is.na(res_etape2[res_etape2=="succes2 DL1"])) && all(is.na(res_etape2[res_etape2=="succes2 DL-1"]))) {
        grid$succes[g] <- 0
      } else {
        grid$succes[g] <- ((table(res_etape2[res_etape2=="succes2 DL1"]) + table(res_etape2[res_etape2=="succes2 DL-1"]))*100)/ntrials
      }
    }
  }
  
  #% arrêt précoce
  if (all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$early_ter[g] <- ((table(res_safety[res_safety=="echec"]))*100)/ntrials
  } else {
    if (all(is.na(res_safety[res_safety=="echec"]))) {
      grid$early_ter[g] <- ((table(res_etape1[res_etape1=="echec"]))*100)/ntrials
    } else {
      if (all(is.na(res_etape1[res_etape1=="echec"])) && all(is.na(res_safety[res_safety=="echec"]))) {
        grid$early_ter[g] <- 0
      } else {
        grid$early_ter[g] <- ((table(res_etape1[res_etape1=="echec"])+table(res_safety[res_safety=="echec"]))*100)/ntrials
      }
    }
  }
  
  #% arrêt pour tox
  if (all(is.na(res_safety[res_safety=="echec"]))) {
    grid$stop_tox[g] <- 0
  } else {
    grid$stop_tox[g] <- (table(res_safety[res_safety=="echec"])*100)/ntrials
  }
  
  
  #nb de succès et échecs à chaque étape
  if (all(is.na(res_safety[res_safety=="succes DL-1"])) && all(is.na(res_safety[res_safety=="succes DL1"]))) {
    grid$res_safe_succes[g] <- 0
  } else {
    grid$res_safe_succes[g] <- table(res_safety[res_safety=="succes DL-1"])+table(res_safety[res_safety=="succes DL1"])
  }
  if (all(is.na(res_safety[res_safety=="echec"]))) {
    grid$res_safe_echecs[g] <- 0
  } else {
    grid$res_safe_echecs[g] <- table(res_safety[res_safety=="echec"])
  }
  
  if (all(is.na(res_etape1[res_etape1=="succes1 DL-1"]))) {
    grid$res1_succes[g] <- table(res_etape1[res_etape1=="succes1 DL1"])
  } else {
    if (all(is.na(res_etape1[res_etape1=="succes1 DL1"]))) {
      grid$res1_succes[g] <- table(res_etape1[res_etape1=="succes1 DL-1"])
    } else {
      if (all(is.na(res_etape1[res_etape1=="succes1 DL1"])) && all(is.na(res_etape1[res_etape1=="succes1 DL-1"]))) {
        grid$res1_succes[g] <- 0
      } else {
        grid$res1_succes[g] <- table(res_etape1[res_etape1=="succes1 DL1"]) + table(res_etape1[res_etape1=="succes1 DL-1"])
      }
    }
  }
  if (all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$res1_echecs[g] <- 0
  } else {
    grid$res1_echecs[g] <- table(res_etape1[res_etape1=="echec"])
  }
  
  if (all(is.na(res_etape2[res_etape2=="succes2 DL-1"]))) {
    grid$res2_succes[g] <- table(res_etape2[res_etape2=="succes2 DL1"])
  } else {
    if (all(is.na(res_etape2[res_etape2=="succes2 DL1"]))) {
      grid$res2_succes[g] <- table(res_etape2[res_etape2=="succes2 DL-1"])
    } else {
      if (all(is.na(res_etape2[res_etape2=="succes2 DL1"])) && all(is.na(res_etape2[res_etape2=="succes2 DL-1"]))) {
        grid$res2_succes[g] <- 0
      } else {
        grid$res2_succes[g] <- table(res_etape2[res_etape2=="succes2 DL1"]) + table(res_etape2[res_etape2=="succes2 DL-1"])
      }
    }
  }
  if (all(is.na(res_etape2[res_etape2=="echec"]))) {
    grid$res2_echecs[g] <- 0
  } else {
    grid$res2_echecs[g] <- table(res_etape2[res_etape2=="echec"])
  }
  
  
  
}


writexl::write_xlsx(grid, "C:\\util\\marie\\2022-2023\\stage\\R\\résultats\\safety2doses_Simon_tpspoisAI.xlsx")


############## Safety (1 dose, 4 sem) + Simon ###################
# Safety lead-in (1 dose) + Simon
# --> CORRECTION OK
# test : code avec matrice -> fait tous les scenarios à la suite 

#alpha = 5%
r1 = 4 ; n1 = 19 ; R = 15 ; N = 54

#alpha = 2,5%
#r1 = 5 ; n1 = 22 ; R = 19 ; N = 66

ntrials=10000
npat=60
n_sl=6

res_safety=res_etape1=res_etape2=rep(NA,ntrials) #nombre d'échecs et succès pour chaque essai
i=0
sample_safety=sample_etape1=sample_etape2=0
nb_eff=0 #nombre de réponses (succès) à chaque étape
nb_ech=0 #nombre de non-réponses (échecs) à chaque étape
nb_tox=0 #nombre de toxicités à chaque étape
time=rep(NA,ntrials)
t_obs=4 #4 semaines (28 jours)

param_eff=c(0.2,0.4,0.4,0.2,0.4,0.2,0.4,0.3)
param_tox=c(0.35,0.15,0.05,0.05,0.45,0.15,0.35,0.2)
scenario=paste("Scenario",1:8)
grid=cbind(data.frame(scenario,param_eff,param_tox))
colnames(grid)=c("Scenario","Peff","Ptox")


for (g in 1:nrow(grid)) {
  
  print(paste("Scenario",g))
  
  peff<-grid$Peff[g]
  ptox<-grid$Ptox[g]
  print(paste("peff=",peff,"ptox=",ptox))
  
  pefftox=peff*ptox
  peffnotox=peff-pefftox
  pnoefftox=ptox-pefftox
  pnoeffnotox=(1-peff)-pnoefftox
  
  res_safety=res_etape1=res_etape2=rep(NA,ntrials) #nombre d'échecs et succès pour chaque essai
  i=0
  sample_safety=sample_etape1=sample_etape2=0
  nb_eff=0 #nombre de réponses (succès) à chaque étape
  nb_ech=0 #nombre de non-réponses (échecs) à chaque étape
  nb_tox=0 #nombre de toxicités à chaque étape
  time=rep(NA,ntrials)
  
  #SAFETY LEAD-IN (1 DOSE) + SIMON
  set.seed(3003)
  for (n in 1:ntrials) {
    
    #Simulation des patients
    patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))),cumsum(c(0,rep(4,n_sl),rpois(N-1,2)))) #temps aléatoire: cumsum(c(0,rep(4,n_sl),rpois(N-1,2))) ; temps fixe: cumsum(c(0,rep(4,n_sl),rep(2,N-1)))
    #patient=t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox)))
    #fixe=cumsum(c(0,rep(4,n_sl),rep(2,N-1)))
    #pois=cumsum(c(0,rep(4,n_sl),rpois(N-1,2)))
    #patient=cbind(patient,pois)
    
    limite=0
    rep_safe=rep_echec=0
    
    #SAFETY LEAD-IN
    for (i in 1:n_sl) {
      limite=limite + patient[,2][i] + patient[,4][i]
      rep_safe=rep_safe + patient[,1][i] + patient[,2][i]
      rep_echec=rep_echec + patient[,3][i] + patient[,4][i]
      if (limite==2) {break()}
    }
    if (limite==2) {res_safety[n]<-"echec"} else {res_safety[n]<-"succes_s"} #0=échec ; 1=succès (passage à l'étape 1)
    sample_safety=sample_safety+i
    nb_eff=nb_eff+rep_safe
    nb_ech=nb_ech+rep_echec
    nb_tox=nb_tox+limite
    #print(paste("essai",n))
    #print(paste(i,"patients",limite,"toxicity",rep_safe,"efficacités"))
    #print(paste("Résultat s_l:",res_safety[n]))
    
    if (res_safety[n]=="echec") {
      time[n]<-patient[i,5]+t_obs
      #cat("Durée de l'essai(saf):",time[n])
    }
    
    #SIMON ETAPE 1
    if (res_safety[n]=="succes_s") {
      
      eff1=sum(patient[n_sl+1:n1,1]) + sum(patient[n_sl+1:n1,2])
      ech1=sum(patient[n_sl+1:n1,3]) + sum(patient[n_sl+1:n1,4])
      tox1=sum(patient[n_sl+1:n1,2]) + sum(patient[n_sl+1:n1,4])
      
      if (eff1>r1) {res_etape1[n]<-"succes1"} else {res_etape1[n]<-"echec"} #0=échec ; 1=succès (passage à l'étape 2)
      sample_etape1=sample_etape1+n1
      nb_eff=nb_eff+eff1
      nb_ech=nb_ech+ech1
      nb_tox=nb_tox+tox1
      #print(paste(n1,"patients",eff1,"efficacités",tox1,"toxicités"))
      #print(paste("Résultat étape1:",res_etape1[n]))
      
      if (res_etape1[n]=="echec") {time[n]<-patient[n1+n_sl,5]+t_obs}
      
      #SIMON ETAPE 2
      if (res_etape1[n]=="succes1") { 
        
        eff2=sum(patient[n_sl+n1+1:(N-n1),1]) + sum(patient[n_sl+n1+1:(N-n1),2])
        ech2=sum(patient[n_sl+n1+1:(N-n1),3]) + sum(patient[n_sl+n1+1:(N-n1),4])
        tox2=sum(patient[n_sl+n1+1:(N-n1),2]) + sum(patient[n_sl+n1+1:(N-n1),4])
        
        if (eff1+eff2>R) {res_etape2[n]<-"succes2"} else {res_etape2[n]<-"echec"} #0=échec ; 1=succès (conclu d'efficacité)
        sample_etape2=sample_etape2+N-n1
        nb_eff=nb_eff+eff2
        nb_ech=nb_ech+ech2
        nb_tox=nb_tox+tox2
        #print(paste(N-n1,"patients",eff2,"efficacités",tox2,"toxicités"))
        #print(paste("Résultat étape2:",res_etape2[n]))
        
        if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) == r1) {time[n] <- patient[N+n_sl,5] + 2*t_obs}
        if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) != r1) {time[n] <- patient[N+n_sl,5] + t_obs}

      }
      
    }
    
  }
  
  #effectif attendu
  grid$ess[g] <- (sample_safety+sample_etape1+sample_etape2)/ntrials
  
  #nb d'efficacités 
  grid$moy_eff[g] <- nb_eff/ntrials
  
  #nb d'échecs
  grid$moy_ech[g] <- nb_ech/ntrials
  
  #nb de toxicités 
  grid$moy_tox[g] <- nb_tox/ntrials
  
  #durée de l'essai
  grid$moy_time[g] <- mean(time)
  
  #% conclu d'efficacité
  if (all(is.na(res_etape2[res_etape2=="succes2"]))) {
    grid$succes[g] <- 0
  } else {
    grid$succes[g] <- ((table(res_etape2[res_etape2=="succes2"]))*100)/ntrials
  }
  
  #% arrêt précoce
  if (all(is.na(res_safety[res_safety=="echec"])) && all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$early_ter[g] <- 0
  } else {
    grid$early_ter[g] <- ((table(res_safety[res_safety=="echec"])+table(res_etape1[res_etape1=="echec"]))*100)/ntrials
  }
  
  #% arrêt pour tox
  if (all(is.na(res_safety[res_safety=="echec"]))) {
    grid$stop_tox[g] <- 0
  } else {
    grid$stop_tox[g] <- (table(res_safety[res_safety=="echec"])*100)/ntrials
  }
  
  
  #nb d'échecs et de succès à chaque étape 
  if (all(is.na(res_safety[res_safety=="succes_s"]))) {
    grid$ressafe_succes[g] <- 0
  } else {
    grid$ressafe_succes[g] <- table(res_safety[res_safety=="succes_s"])
  }
  if (all(is.na(res_safety[res_safety=="echec"]))) {
    grid$ressafe_echecs[g] <- 0
  } else {
    grid$ressafe_echecs[g] <- table(res_safety[res_safety=="echec"])
  }
  
  if (all(is.na(res_etape1[res_etape1=="succes1"]))) {
    grid$res1_succes[g] <- 0
  } else {
    grid$res1_succes[g] <- table(res_etape1[res_etape1=="succes1"])
  }
  if (all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$res1_echecs[g] <- 0
  } else {
    grid$res1_echecs[g] <- table(res_etape1[res_etape1=="echec"])
  }
  
  if (all(is.na(res_etape2[res_etape2=="succes2"]))) {
    grid$res2_succes[g] <- 0
  } else {
    grid$res2_succes[g] <- table(res_etape2[res_etape2=="succes2"])
  }
  if (all(is.na(res_etape2[res_etape2=="echec"]))) {
    grid$res2_echecs[g] <- 0
  } else {
    grid$res2_echecs[g] <- table(res_etape2[res_etape2=="echec"])
  }
  
}

writexl::write_xlsx(grid, "C:\\util\\marie\\2022-2023\\stage\\R\\résultats\\safety1dose_Simon_tpspois.xlsx")





############## Safety (2 doses, 4 sem) + Simon ####################
# Safety lead-in (2 doses) + Simon


ntrials=10000
npat=60

#alpha = 5%
r1 = 4 ; n1 = 19 ; R = 15 ; N = 54

#alpha = 2,5%
#r1 = 5 ; n1 = 22 ; R = 19 ; N = 66

n_sl=6

res_safety=res_etape1=res_etape2=rep(NA,ntrials) #nombre d'échecs et succès pour chaque essai
i=l=0
sample_s=sample1_d1=sample2_d1=0
sample_s2=sample1_d2=sample2_d2=0
nb_eff=0 #nombre de réponses (succès) à chaque étape
nb_ech=0 #nombre de non-réponses (échecs) à chaque étape
nb_tox=0 #nombre de toxicités à chaque étape
time=rep(NA,ntrials)
t_obs=4 #mesure de l'outcome après 28j (4 semaines)
timesafe=rep(NA,ntrials) #stocke la durée du 1er safety lead-in (utile si désescalade de dose)


# param_eff=c(0.2,0.4,0.2,0.4,0.4,0.2,0.2,0.3)
# param_tox=c(0.2,0.2,0.05,0.05,0.35,0.35,0.1,0.2)
# param_eff2=c(0.4,0.5,0.4,0.5,0.5,0.4,0.3,0.4)
# param_tox2=c(0.35,0.35,0.2,0.2,0.45,0.45,0.2,0.3)
param_eff=c(0.2, 0.2, 0.4, 0.4, 0.2, 0.4, 0.2, 0.1, 0.4, 0.4, 0.2, 0.3, 0.1)
param_tox=c(0.35, 0.05, 0.15, 0.05, 0.35, 0.45, 0.15, 0.05, 0.35, 0.15, 0.15, 0.2, 0.15)
param_eff2=c(0.4, 0.4, 0.5, 0.5, 0.3, 0.6, 0.4, 0.2, 0.5, 0.6, 0.3, 0.4, 0.2)
param_tox2=c(0.5, 0.15, 0.35, 0.15, 0.45, 0.5, 0.35, 0.15, 0.45, 0.45, 0.2, 0.35, 0.35)
scenario=paste("Scenario",1:13)
grid=cbind(data.frame(scenario,param_eff,param_tox,param_eff2,param_tox2))
colnames(grid)=c("Scenario","Peff","Ptox","Peff2","Ptox2")



for (g in 1:nrow(grid)) {
  print(paste("Scenario",g))
  
  #pour DL1
  peff <- grid$Peff2[g]
  ptox <- grid$Ptox2[g]
  
  pefftox=peff*ptox
  peffnotox=peff-pefftox
  pnoefftox=ptox-pefftox
  pnoeffnotox=(1-peff)-pnoefftox
  
  #pour DL-1
  peff2 <- grid$Peff[g]
  ptox2 <- grid$Ptox[g]
  
  pefftox2=peff2*ptox2 
  peffnotox2=peff2-pefftox2
  pnoefftox2=ptox2-pefftox2
  pnoeffnotox2=(1-peff2)-pnoefftox2
  
  
  res_safety=res_etape1=res_etape2=rep(NA,ntrials) #nombre d'échecs et succès pour chaque essai
  i=l=0
  sample_s=sample1_d1=sample2_d1=0
  sample_s2=sample1_d2=sample2_d2=0
  nb_eff=0 #nombre de réponses (succès) à chaque étape
  nb_ech=0 #nombre de non-réponses (échecs) à chaque étape
  nb_tox=0 #nombre de toxicités à chaque étape
  time=rep(NA,ntrials)
  timesafe=rep(NA,ntrials) #stocke la durée du 1er safety lead-in (utile quand désescalade de dose)
  
  
  set.seed(3003)
  for (n in 1:ntrials) {
    patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))),cumsum(c(0,rep(4,n_sl),rpois(N-1,2)))) #temps aléatoire:  cumsum(c(0,rep(4,n_sl),rpois(N-1,2))) temps fixe: cumsum(c(0,rep(4,n_sl),rep(2,N-1)))
    #patient1=t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox)))
    #fixe=cumsum(c(0,rep(4,n_sl),rep(2,N-1)))
    #pois=cumsum(c(0,rep(4,n_sl),rpois(N-1,2)))
    #patient=cbind(patient1,pois)
    
    limite=limite2=0
    rep_safe=rep_echec=0
    
    #SAFETY LEAD-IN (DL1)
    for (i in 1:n_sl) {
      limite=limite + patient[,2][i] + patient[,4][i]
      rep_safe=rep_safe + patient[,1][i] + patient[,2][i]  
      rep_echec=rep_echec + patient[,3][i] + patient[,4][i]
      if (limite==2) {break()}
    }
    if (limite==2) {res_safety[n]<-"echec"} else {res_safety[n]<-"succes DL1"} #0=échec ; 1=succès (passage à l'étape 1)
    sample_s=sample_s+i
    nb_eff=nb_eff+rep_safe
    nb_ech=nb_ech+rep_echec
    nb_tox=nb_tox+limite
    #print(paste("essai",n))
    #print(paste(i,"patients",limite,"toxicity",rep_safe,"efficacités"))
    #print(paste("Résultat s_l:",res_safety[n]))
    
    if (res_safety[n]=="echec") {
      time[n]<-patient[i,5]+t_obs
      timesafe[n]<-patient[i,5]+t_obs
      #cat("Durée de l'essai(saf):",time[n])
    }
    
    #SIMON ETAPE 1 (DL1)
    if (res_safety[n]=="succes DL1") {
      eff1= sum(patient[n_sl+1:n1,1]) + sum(patient[n_sl+1:n1,2])
      ech1= sum(patient[n_sl+1:n1,3]) + sum(patient[n_sl+1:n1,4])
      tox1= sum(patient[n_sl+1:n1,2]) + sum(patient[n_sl+1:n1,4])
      
      if (eff1>r1) {res_etape1[n]<-"succes1 DL1"} else {res_etape1[n]<-"echec"} #0=échec ; 1=succès (passage à l'étape 2)
      sample1_d1=sample1_d1+n1
      nb_eff=nb_eff+eff1
      nb_ech=nb_ech+ech1
      nb_tox=nb_tox+tox1
      #print(paste(n1,"patients",eff1,"efficacités",tox1,"toxicités"))
      #print(paste("Résultat étape1:",res_etape1[n]))
      
      if (res_etape1[n]=="echec") {time[n]<-patient[n1+n_sl,5]+t_obs}
      
      #SIMON ETAPE 2 (DL1)
      if (res_etape1[n]=="succes1 DL1") {
        
        eff2=sum(patient[n_sl+n1+1:(N-n1),1]) + sum(patient[n_sl+n1+1:(N-n1),2])
        ech2=sum(patient[n_sl+n1+1:(N-n1),3]) + sum(patient[n_sl+n1+1:(N-n1),4])
        tox2=sum(patient[n_sl+n1+1:(N-n1),2]) + sum(patient[n_sl+n1+1:(N-n1),4])
        
        if (eff1+eff2>R) {res_etape2[n]<-"succes2 DL1"} else {res_etape2[n]<-"echec"} #0=échec ; 1=succès (conclu d'efficacité)
        sample2_d1=sample2_d1+N-n1
        nb_eff=nb_eff+eff2
        nb_ech=nb_ech+ech2
        nb_tox=nb_tox+tox2
        #print(paste(N-n1,"patients",eff2,"efficacités",tox2,"toxicités"))
        #print(paste("Résultat étape2:",res_etape2[n]))
        
        if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) == r1) {time[n]<-patient[N+n_sl,5]+2*t_obs}
        if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) != r1) {time[n]<-patient[N+n_sl,5]+t_obs}
        
      }
      
    }
    
    if (res_safety[n]=="echec") {
      
      #SAFETY LEAD-IN (DL-1)
      patient2=cbind(t(rmultinom(npat,1,c(peffnotox2,pefftox2,pnoeffnotox2,pnoefftox2))),cumsum(c(0,rep(4,n_sl),rpois(N-1,2))))
      #patient_2=t(rmultinom(npat,1,c(peffnotox2,pefftox2,pnoeffnotox2,pnoefftox2)))
      #fixe=cumsum(c(0,rep(4,n_sl),rep(2,N-1)))
      #pois2=cumsum(c(0,rep(4,n_sl),rpois(N-1,2)))
      #patient2=cbind(patient_2,pois2)
      
      limite2=0
      rep_safe=rep_echec=0
      
      for (l in 1:n_sl) {
        limite2 = limite2 + patient2[,2][l] + patient2[,4][l]
        rep_safe=rep_safe + patient2[,1][l] + patient2[,2][l]  
        rep_echec=rep_echec + patient2[,3][l] + patient2[,4][l]
        if (limite2==2) {break()}
      }
      if (limite2==2) {res_safety[n]<-"echec"} else {res_safety[n]<-"succes DL-1"}
      sample_s2=sample_s2+l
      nb_eff=nb_eff+rep_safe
      nb_ech=nb_ech+rep_echec
      nb_tox=nb_tox+limite2
      #print(paste("essai",n))
      #print(paste(l,"patients",limite,"toxicity",rep_safe,"efficacités"))
      #print(paste("Résultat s_l:",res_safety[n]))
      
      if (res_safety[n]=="echec") {
        time[n]<-timesafe[n] + patient2[l,5]+t_obs
        #cat("Durée de l'essai(saf):",time[n])
      }
      
      #SIMON ETAPE 1 (DL-1)
      if (res_safety[n]=="succes DL-1") {
        eff1=sum(patient2[n_sl+1:n1,1]) + sum(patient2[n_sl+1:n1,2])
        ech1=sum(patient2[n_sl+1:n1,3]) + sum(patient2[n_sl+1:n1,4])
        tox1=sum(patient2[n_sl+1:n1,2]) + sum(patient2[n_sl+1:n1,4])
        
        if (eff1>r1) {res_etape1[n]<-"succes1 DL-1"} else {res_etape1[n]<-"echec"} #0=échec ; 1=succès (passage à l'étape 2)
        sample1_d2=sample1_d2+n1
        nb_eff=nb_eff+eff1
        nb_ech=nb_ech+ech1
        nb_tox=nb_tox+tox1
        #print(paste(n1,"patients",eff1,"efficacités",tox1,"toxicités"))
        #print(paste("Résultat étape1:",res_etape1[n]))
        
        if (res_etape1[n]=="echec") {patient2[n1+n_sl,5]+t_obs + timesafe[n]}
      
        
        #SIMON ETAPE 2 (DL-1)
        if (res_etape1[n]=="succes1 DL-1") {
          eff2=sum(patient2[n_sl+n1+1:(N-n1),1]) + sum(patient2[n_sl+n1+1:(N-n1),2])
          ech2=sum(patient2[n_sl+n1+1:(N-n1),3]) + sum(patient2[n_sl+n1+1:(N-n1),4])
          tox2=sum(patient2[n_sl+n1+1:(N-n1),2]) + sum(patient2[n_sl+n1+1:(N-n1),4])
          
          if (eff1+eff2>R) {res_etape2[n]<-"succes2 DL-1"} else {res_etape2[n]<-"echec"} #0=échec ; 1=succès (conclu d'efficacité)
          sample2_d2=sample2_d2+N-n1
          nb_eff=nb_eff+eff2
          nb_ech=nb_ech+ech2
          nb_tox=nb_tox+tox2
          #print(paste(N-n1,"patients",eff2,"efficacités",tox2,"toxicités"))
          #print(paste("Résultat étape2:",res_etape2[n]))
          
          if (sum(patient2[n_sl+1:(n1-1),1])+sum(patient2[n_sl+1:(n1-1),2]) == r1) {time[n]<-timesafe[n]+patient2[N+n_sl,5]+2*t_obs}
          if (sum(patient2[n_sl+1:(n1-1),1])+sum(patient2[n_sl+1:(n1-1),2]) != r1) {time[n]<-timesafe[n]+patient2[N+n_sl,5]+t_obs}
          
        }
      }
    }
  }
  
  
  #effectif attendu
  grid$ess[g] <- (sample_s+sample_s2+sample1_d1+sample2_d1+sample1_d2+sample2_d2)/ntrials
  
  #nb d'efficacités 
  grid$moy_eff[g] <- nb_eff/ntrials
  
  #nb d'échecs
  grid$moy_ech[g] <- nb_ech/ntrials
  
  #nb de toxicités 
  grid$moy_tox[g] <- nb_tox/ntrials
  
  #durée de l'essai
  grid$moy_time[g] <- mean(time)
  
  
  #% conclu d'efficacité
  if (all(is.na(res_etape2[res_etape2=="succes2 DL-1"]))) {
    grid$succes[g] <- ((table(res_etape2[res_etape2=="succes2 DL1"]))*100)/ntrials
  } else {
    if (all(is.na(res_etape2[res_etape2=="succes2 DL1"]))) {
      grid$succes[g] <- ((table(res_etape2[res_etape2=="succes2 DL-1"]))*100)/ntrials
    } else {
      if (all(is.na(res_etape2[res_etape2=="succes2 DL1"])) && all(is.na(res_etape2[res_etape2=="succes2 DL-1"]))) {
        grid$succes[g] <- 0
      } else {
        grid$succes[g] <- ((table(res_etape2[res_etape2=="succes2 DL1"]) + table(res_etape2[res_etape2=="succes2 DL-1"]))*100)/ntrials
      }
    }
  }
  
  #% arrêt précoce
  if (all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$early_ter[g] <- ((table(res_safety[res_safety=="echec"]))*100)/ntrials
  } else {
    if (all(is.na(res_safety[res_safety=="echec"]))) {
      grid$early_ter[g] <- ((table(res_etape1[res_etape1=="echec"]))*100)/ntrials
    } else {
      if (all(is.na(res_etape1[res_etape1=="echec"])) && all(is.na(res_safety[res_safety=="echec"]))) {
        grid$early_ter[g] <- 0
      } else {
        grid$early_ter[g] <- ((table(res_etape1[res_etape1=="echec"])+table(res_safety[res_safety=="echec"]))*100)/ntrials
      }
    }
  }
  
  #% arrêt pour tox
  if (all(is.na(res_safety[res_safety=="echec"]))) {
    grid$stop_tox[g] <- 0
  } else {
    grid$stop_tox[g] <- (table(res_safety[res_safety=="echec"])*100)/ntrials
  }
  
  
  #nb de succès et échecs à chaque étape
  if (all(is.na(res_safety[res_safety=="succes DL-1"])) && all(is.na(res_safety[res_safety=="succes DL1"]))) {
    grid$res_safe_succes[g] <- 0
  } else {
    grid$res_safe_succes[g] <- table(res_safety[res_safety=="succes DL-1"])+table(res_safety[res_safety=="succes DL1"])
  }
  if (all(is.na(res_safety[res_safety=="echec"]))) {
    grid$res_safe_echecs[g] <- 0
  } else {
    grid$res_safe_echecs[g] <- table(res_safety[res_safety=="echec"])
  }
  
  if (all(is.na(res_etape1[res_etape1=="succes1 DL-1"]))) {
    grid$res1_succes[g] <- table(res_etape1[res_etape1=="succes1 DL1"])
  } else {
    if (all(is.na(res_etape1[res_etape1=="succes1 DL1"]))) {
      grid$res1_succes[g] <- table(res_etape1[res_etape1=="succes1 DL-1"])
    } else {
      if (all(is.na(res_etape1[res_etape1=="succes1 DL1"])) && all(is.na(res_etape1[res_etape1=="succes1 DL-1"]))) {
        grid$res1_succes[g] <- 0
      } else {
        grid$res1_succes[g] <- table(res_etape1[res_etape1=="succes1 DL1"]) + table(res_etape1[res_etape1=="succes1 DL-1"])
      }
    }
  }
  if (all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$res1_echecs[g] <- 0
  } else {
    grid$res1_echecs[g] <- table(res_etape1[res_etape1=="echec"])
  }
  
  if (all(is.na(res_etape2[res_etape2=="succes2 DL-1"]))) {
    grid$res2_succes[g] <- table(res_etape2[res_etape2=="succes2 DL1"])
  } else {
    if (all(is.na(res_etape2[res_etape2=="succes2 DL1"]))) {
      grid$res2_succes[g] <- table(res_etape2[res_etape2=="succes2 DL-1"])
    } else {
      if (all(is.na(res_etape2[res_etape2=="succes2 DL1"])) && all(is.na(res_etape2[res_etape2=="succes2 DL-1"]))) {
        grid$res2_succes[g] <- 0
      } else {
        grid$res2_succes[g] <- table(res_etape2[res_etape2=="succes2 DL1"]) + table(res_etape2[res_etape2=="succes2 DL-1"])
      }
    }
  }
  if (all(is.na(res_etape2[res_etape2=="echec"]))) {
    grid$res2_echecs[g] <- 0
  } else {
    grid$res2_echecs[g] <- table(res_etape2[res_etape2=="echec"])
  }
  
  
  
}


writexl::write_xlsx(grid, "C:\\util\\marie\\2022-2023\\stage\\R\\résultats\\safety2doses_Simon_tpspois.xlsx")


table(res_safety)
table(res_etape1)
table(res_etape2)

(sample_s+sample_s2+
    sample1_d1+sample2_d1+
    sample1_d2+sample2_d2)/ntrials

#nombre de réponses/efficacités
nb_eff/ntrials

#nombre de non réponses/non-efficacités
nb_ech/ntrials

#nombre de toxicités
nb_tox/ntrials

#durée de l'essai
mean(time)
table(time)





############## code initial (sans grid) #################
# Safety + Simon

#modification du code pour le safety lead-in 
#inclusion des patients toutes les 2 semaines (comme pour le reste des designs)

#alpha = 5%
r1 = 4 ; n1 = 19 ; R = 15 ; N = 54

#alpha = 2,5%
#r1 = 5 ; n1 = 22 ; R = 19 ; N = 66


peff=0.2
ptox=0.35

pefftox=peff*ptox
peffnotox=peff-pefftox
pnoefftox=ptox-pefftox
pnoeffnotox=(1-peff)-pnoefftox

ntrials=10000
npat=60
n_sl=6

res_safety=res_etape1=res_etape2=rep(NA,ntrials) #nombre d'échecs et succès pour chaque essai
i=0
sample_safety=sample_etape1=sample_etape2=0
nb_eff=0 #nombre de réponses (succès) à chaque étape
nb_ech=0 #nombre de non-réponses (échecs) à chaque étape
nb_tox=0 #nombre de toxicités à chaque étape
time=rep(NA,ntrials)
t_obs=4 #4 semaines (28 jours)


set.seed(3003)
for (n in 1:ntrials) {
  
  #Simulation des patients
  patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))),cumsum(c(0,rep(2,N+n_sl-1)))) #temps aléatoire: cumsum(c(0,rep(4,n_sl),rpois(N-1,2))) ; temps fixe: cumsum(c(0,rep(4,n_sl),rep(2,N-1)))
  #patient=t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox)))
  #fixe=cumsum(c(0,rep(4,n_sl),rep(2,N-1)))
  #pois=cumsum(c(0,rep(4,n_sl),rpois(N-1,2)))
  #patient=cbind(patient,pois)
  
  limite=0
  rep_safe=rep_echec=0
  
  #SAFETY LEAD-IN
  for (i in 1:n_sl) {
    limite=limite + patient[,2][i] + patient[,4][i]
    rep_safe=rep_safe + patient[,1][i] + patient[,2][i]
    rep_echec=rep_echec + patient[,3][i] + patient[,4][i]
    if (limite==2) {break()}
  }
  if (limite==2) {res_safety[n]<-"echec"} else {res_safety[n]<-"succes_s"} 
  if (res_safety[n]=="echec") {
    sample_safety=sample_safety+i+1
    nb_eff=nb_eff+rep_safe+patient[,1][i+1] + patient[,2][i+1]
    nb_ech=nb_ech+rep_echec+patient[,3][i+1] + patient[,4][i+1]
    nb_tox=nb_tox+limite+patient[,2][i+1] + patient[,4][i+1]
  } else {
      sample_safety=sample_safety+i
      nb_eff=nb_eff+rep_safe
      nb_ech=nb_ech+rep_echec
      nb_tox=nb_tox+limite
      }
  
  #print(paste("essai",n))
  #print(paste(i,"patients",limite,"toxicity",rep_safe,"efficacités"))
  #print(paste("Résultat s_l:",res_safety[n]))
  
  if (res_safety[n]=="echec") {
    time[n]<-patient[i+1,5]+t_obs
    #cat("Durée de l'essai(saf):",time[n])
  }
  
  #SIMON ETAPE 1
  if (res_safety[n]=="succes_s") {
    
    eff1=sum(patient[n_sl+1:n1,1]) + sum(patient[n_sl+1:n1,2])
    ech1=sum(patient[n_sl+1:n1,3]) + sum(patient[n_sl+1:n1,4])
    tox1=sum(patient[n_sl+1:n1,2]) + sum(patient[n_sl+1:n1,4])
    
    if (eff1>r1) {res_etape1[n]<-"succes1"} else {res_etape1[n]<-"echec"} #0=échec ; 1=succès (passage à l'étape 2)
    sample_etape1=sample_etape1+n1
    nb_eff=nb_eff+eff1
    nb_ech=nb_ech+ech1
    nb_tox=nb_tox+tox1
    #print(paste(n1,"patients",eff1,"efficacités",tox1,"toxicités"))
    #print(paste("Résultat étape1:",res_etape1[n]))
    
    if (res_etape1[n]=="echec") {
      if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) == limite-1) {time[n]<-patient[n1+n_sl,5]+2*t_obs}
      else {time[n]<-patient[n1+n_sl,5]+t_obs}
      #cat("Durée de l'essai(1):",time[n])
    }
    
    #SIMON ETAPE 2
    if (res_etape1[n]=="succes1") { 
      
      eff2=sum(patient[n_sl+n1+1:(N-n1),1]) + sum(patient[n_sl+n1+1:(N-n1),2])
      ech2=sum(patient[n_sl+n1+1:(N-n1),3]) + sum(patient[n_sl+n1+1:(N-n1),4])
      tox2=sum(patient[n_sl+n1+1:(N-n1),2]) + sum(patient[n_sl+n1+1:(N-n1),4])
      
      if (eff1+eff2>R) {res_etape2[n]<-"succes2"} else {res_etape2[n]<-"echec"} #0=échec ; 1=succès (conclu d'efficacité)
      sample_etape2=sample_etape2+N-n1
      nb_eff=nb_eff+eff2
      nb_ech=nb_ech+ech2
      nb_tox=nb_tox+tox2
      #print(paste(N-n1,"patients",eff2,"efficacités",tox2,"toxicités"))
      #print(paste("Résultat étape2:",res_etape2[n]))
      
      if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) == limite-1) {
        if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) == r1) {time[n] <- patient[N+n_sl,5] + 3*t_obs}
        if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) != r1) {time[n] <- patient[N+n_sl,5] + 2*t_obs}
      }
      if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) != limite-1) {
        if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) == r1) {time[n] <- patient[N+n_sl,5] + 2*t_obs}
        if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) != r1) {time[n] <- patient[N+n_sl,5] + t_obs}
      }
    
    }
    
  }
  
}

table(res_safety)
table(res_etape1)
table(res_etape2)

(sample_safety+sample_etape1+sample_etape2)/ntrials

nb_eff/ntrials

nb_ech/ntrials

nb_tox/ntrials

mean(time)



############## Safety (1 dose, 2 sem) + Simon ################
# code pour tous les scenarios 

#alpha = 5%
r1 = 4 ; n1 = 19 ; R = 15 ; N = 54

#alpha = 2,5%
#r1 = 5 ; n1 = 22 ; R = 19 ; N = 66

ntrials=10000
npat=60
n_sl=6

res_safety=res_etape1=res_etape2=rep(NA,ntrials) #nombre d'échecs et succès pour chaque essai
i=0
sample_safety=sample_etape1=sample_etape2=0
nb_eff=0 #nombre de réponses (succès) à chaque étape
nb_ech=0 #nombre de non-réponses (échecs) à chaque étape
nb_tox=0 #nombre de toxicités à chaque étape
time=rep(NA,ntrials)
t_obs=4 #4 semaines (28 jours)

param_eff=c(0.2,0.4,0.4,0.2,0.4,0.2,0.4,0.3)
param_tox=c(0.35,0.15,0.05,0.05,0.45,0.15,0.35,0.2)
scenario=paste("Scenario",1:8)
grid=cbind(data.frame(scenario,param_eff,param_tox))
colnames(grid)=c("Scenario","Peff","Ptox")


for (g in 1:nrow(grid)) {
  
  print(paste("Scenario",g))
  
  peff<-grid$Peff[g]
  ptox<-grid$Ptox[g]
  print(paste("peff=",peff,"ptox=",ptox))
  
  pefftox=peff*ptox
  peffnotox=peff-pefftox
  pnoefftox=ptox-pefftox
  pnoeffnotox=(1-peff)-pnoefftox
  
  res_safety=res_etape1=res_etape2=rep(NA,ntrials) #nombre d'échecs et succès pour chaque essai
  i=0
  sample_safety=sample_etape1=sample_etape2=0
  nb_eff=0 #nombre de réponses (succès) à chaque étape
  nb_ech=0 #nombre de non-réponses (échecs) à chaque étape
  nb_tox=0 #nombre de toxicités à chaque étape
  time=rep(NA,ntrials)
  check_stop1=check_stop2=rep(NA,ntrials)
  
  #SAFETY LEAD-IN (1 DOSE) + SIMON
  set.seed(3003)
  for (n in 1:ntrials) {
    
    #Simulation des patients
    patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))),cumsum(c(0,rpois(N+n_sl-1,2)))) 
    #patient=t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox)))
    #fixe=cumsum(c(0,rep(4,n_sl),rep(2,N-1)))
    #pois=cumsum(c(0,rep(4,n_sl),rpois(N-1,2)))
    #patient=cbind(patient,pois)
    
    limite=0
    rep_safe=rep_echec=0
    
    #SAFETY LEAD-IN
    for (i in 1:n_sl) {
      limite=limite + patient[,2][i] + patient[,4][i]
      rep_safe=rep_safe + patient[,1][i] + patient[,2][i]
      rep_echec=rep_echec + patient[,3][i] + patient[,4][i]
      if (limite==2) {break()}
    }
    if (limite==2) {res_safety[n]<-"echec"} else {res_safety[n]<-"succes_s"} 
    if (res_safety[n]=="echec") {
      sample_safety=sample_safety+i+1
      nb_eff=nb_eff+rep_safe+patient[,1][i+1] + patient[,2][i+1]
      nb_ech=nb_ech+rep_echec+patient[,3][i+1] + patient[,4][i+1]
      nb_tox=nb_tox+limite+patient[,2][i+1] + patient[,4][i+1]
    } else {
      sample_safety=sample_safety+i
      nb_eff=nb_eff+rep_safe
      nb_ech=nb_ech+rep_echec
      nb_tox=nb_tox+limite
    }
    
    #print(paste("essai",n))
    #print(paste(i,"patients",limite,"toxicity",rep_safe,"efficacités"))
    #print(paste("Résultat s_l:",res_safety[n]))
    
    if (res_safety[n]=="echec") {
      time[n]<-patient[i+1,5]+t_obs
      #cat("Durée de l'essai(saf):",time[n])
    }
    
    #vérification inclusion AI
    if (patient[n_sl,5]==patient[n_sl+1,5]) {check_stop1[n] <- "STOP"} else {"OK"}
    if (patient[n_sl+n1,5]==patient[n_sl+n1+1,5]) {check_stop2[n] <- "STOP"} else {"OK"}
    
    #SIMON ETAPE 1
    if (res_safety[n]=="succes_s") {
      
      eff1=sum(patient[n_sl+1:n1,1]) + sum(patient[n_sl+1:n1,2])
      ech1=sum(patient[n_sl+1:n1,3]) + sum(patient[n_sl+1:n1,4])
      tox1=sum(patient[n_sl+1:n1,2]) + sum(patient[n_sl+1:n1,4])
      
      if (eff1>r1) {res_etape1[n]<-"succes1"} else {res_etape1[n]<-"echec"} #0=échec ; 1=succès (passage à l'étape 2)
      sample_etape1=sample_etape1+n1
      nb_eff=nb_eff+eff1
      nb_ech=nb_ech+ech1
      nb_tox=nb_tox+tox1
      #print(paste(n1,"patients",eff1,"efficacités",tox1,"toxicités"))
      #print(paste("Résultat étape1:",res_etape1[n]))
      
      if (res_etape1[n]=="echec") {
        if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) == limite-1) {time[n]<-patient[n1+n_sl,5]+2*t_obs}
        else {time[n]<-patient[n1+n_sl,5]+t_obs}
        #cat("Durée de l'essai(1):",time[n])
      }
      
      
      #SIMON ETAPE 2
      if (res_etape1[n]=="succes1") { 
        
        eff2=sum(patient[n_sl+n1+1:(N-n1),1]) + sum(patient[n_sl+n1+1:(N-n1),2])
        ech2=sum(patient[n_sl+n1+1:(N-n1),3]) + sum(patient[n_sl+n1+1:(N-n1),4])
        tox2=sum(patient[n_sl+n1+1:(N-n1),2]) + sum(patient[n_sl+n1+1:(N-n1),4])
        
        if (eff1+eff2>R) {res_etape2[n]<-"succes2"} else {res_etape2[n]<-"echec"} #0=échec ; 1=succès (conclu d'efficacité)
        sample_etape2=sample_etape2+N-n1
        nb_eff=nb_eff+eff2
        nb_ech=nb_ech+ech2
        nb_tox=nb_tox+tox2
        #print(paste(N-n1,"patients",eff2,"efficacités",tox2,"toxicités"))
        #print(paste("Résultat étape2:",res_etape2[n]))
        
        
        if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) == limite-1) {
          if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) == r1) {time[n] <- patient[N+n_sl,5] + 3*t_obs}
          if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) != r1) {time[n] <- patient[N+n_sl,5] + 2*t_obs}
        }
        if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) != limite-1) {
          if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) == r1) {time[n] <- patient[N+n_sl,5] + 2*t_obs}
          if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) != r1) {time[n] <- patient[N+n_sl,5] + t_obs}
        }
        #cat("Durée de l'essai(2):",time[n])
        
        
      }
      
    }
    
  }
  
  #effectif attendu
  grid$ess[g] <- (sample_safety+sample_etape1+sample_etape2)/ntrials
  
  #nb d'efficacités 
  grid$moy_eff[g] <- nb_eff/ntrials
  
  #nb d'échecs
  grid$moy_ech[g] <- nb_ech/ntrials
  
  #nb de toxicités 
  grid$moy_tox[g] <- nb_tox/ntrials
  
  #durée de l'essai
  grid$moy_time[g] <- mean(time)
  
  #% conclu d'efficacité
  if (all(is.na(res_etape2[res_etape2=="succes2"]))) {
    grid$succes[g] <- 0
  } else {
    grid$succes[g] <- ((table(res_etape2[res_etape2=="succes2"]))*100)/ntrials
  }
  
  #% arrêt précoce
  if (all(is.na(res_safety[res_safety=="echec"])) && all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$early_ter[g] <- 0
  } else {
    grid$early_ter[g] <- ((table(res_safety[res_safety=="echec"])+table(res_etape1[res_etape1=="echec"]))*100)/ntrials
  }
  
  #% arrêt pour tox
  if (all(is.na(res_safety[res_safety=="echec"]))) {
    grid$stop_tox[g] <- 0
  } else {
    grid$stop_tox[g] <- (table(res_safety[res_safety=="echec"])*100)/ntrials
  }
  
  
  #nb d'échecs et de succès à chaque étape 
  if (all(is.na(res_safety[res_safety=="succes_s"]))) {
    grid$ressafe_succes[g] <- 0
  } else {
    grid$ressafe_succes[g] <- table(res_safety[res_safety=="succes_s"])
  }
  if (all(is.na(res_safety[res_safety=="echec"]))) {
    grid$ressafe_echecs[g] <- 0
  } else {
    grid$ressafe_echecs[g] <- table(res_safety[res_safety=="echec"])
  }
  
  if (all(is.na(res_etape1[res_etape1=="succes1"]))) {
    grid$res1_succes[g] <- 0
  } else {
    grid$res1_succes[g] <- table(res_etape1[res_etape1=="succes1"])
  }
  if (all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$res1_echecs[g] <- 0
  } else {
    grid$res1_echecs[g] <- table(res_etape1[res_etape1=="echec"])
  }
  
  if (all(is.na(res_etape2[res_etape2=="succes2"]))) {
    grid$res2_succes[g] <- 0
  } else {
    grid$res2_succes[g] <- table(res_etape2[res_etape2=="succes2"])
  }
  if (all(is.na(res_etape2[res_etape2=="echec"]))) {
    grid$res2_echecs[g] <- 0
  } else {
    grid$res2_echecs[g] <- table(res_etape2[res_etape2=="echec"])
  }
  
  #vérification inclusion AI
  if (all(is.na(check_stop1[check_stop1=="STOP"]))) {
    grid$check_AI1[g] <- 0
  } else {
    grid$check_AI1[g] <- table(check_stop1[check_stop1=="STOP"])
  }
  if (all(is.na(check_stop2[check_stop2=="STOP"]))) {
    grid$check_AI2[g] <- 0
  } else {
    grid$check_AI2[g] <- table(check_stop2[check_stop2=="STOP"])
  }
  
}



writexl::write_xlsx(grid,"C:\\util\\marie\\2022-2023\\stage\\R\\résultats\\safety1dose_Simon_tpsfixe_2sem.xlsx")


















############## Safety (2 doses, 2 sem) + Simon ################
# Safety lead-in (2 doses) + Simon

#temps d'inclusion = 2 semaines partout


ntrials=10000
npat=60

#alpha = 5%
r1 = 4 ; n1 = 19 ; R = 15 ; N = 54

#alpha = 2,5%
#r1 = 5 ; n1 = 22 ; R = 19 ; N = 66

n_sl=6

res_safety=res_etape1=res_etape2=rep(NA,ntrials) #nombre d'échecs et succès pour chaque essai
i=l=0
sample_s=sample1_d1=sample2_d1=0
sample_s2=sample1_d2=sample2_d2=0
nb_eff=0 #nombre de réponses (succès) à chaque étape
nb_ech=0 #nombre de non-réponses (échecs) à chaque étape
nb_tox=0 #nombre de toxicités à chaque étape
time=rep(NA,ntrials)
t_obs=4 #mesure de l'outcome après 28j (4 semaines)
timesafe=rep(NA,ntrials) #stocke la durée du 1er safety lead-in (utile si désescalade de dose)


# param_eff=c(0.2,0.4,0.2,0.4,0.4,0.2,0.2,0.3)
# param_tox=c(0.2,0.2,0.05,0.05,0.35,0.35,0.1,0.2)
# param_eff2=c(0.4,0.5,0.4,0.5,0.5,0.4,0.3,0.4)
# param_tox2=c(0.35,0.35,0.2,0.2,0.45,0.45,0.2,0.3)
param_eff=c(0.2, 0.2, 0.4, 0.4, 0.2, 0.4, 0.2, 0.1, 0.4, 0.4, 0.2, 0.3, 0.1)
param_tox=c(0.35, 0.05, 0.15, 0.05, 0.35, 0.45, 0.15, 0.05, 0.35, 0.15, 0.15, 0.2, 0.15)
param_eff2=c(0.4, 0.4, 0.5, 0.5, 0.3, 0.6, 0.4, 0.2, 0.5, 0.6, 0.3, 0.4, 0.2)
param_tox2=c(0.5, 0.15, 0.35, 0.15, 0.45, 0.5, 0.35, 0.15, 0.45, 0.45, 0.2, 0.35, 0.35)
scenario=paste("Scenario",1:13)
grid=cbind(data.frame(scenario,param_eff,param_tox,param_eff2,param_tox2))
colnames(grid)=c("Scenario","Peff","Ptox","Peff2","Ptox2")
 


for (g in 1:nrow(grid)) {
  print(paste("Scenario",g))
  
  #pour DL1
  peff <- grid$Peff2[g]
  ptox <- grid$Ptox2[g]
  
  pefftox=peff*ptox
  peffnotox=peff-pefftox
  pnoefftox=ptox-pefftox
  pnoeffnotox=(1-peff)-pnoefftox
  
  #pour DL-1
  peff2 <- grid$Peff[g]
  ptox2 <- grid$Ptox[g]
  
  pefftox2=peff2*ptox2 
  peffnotox2=peff2-pefftox2
  pnoefftox2=ptox2-pefftox2
  pnoeffnotox2=(1-peff2)-pnoefftox2
  
  
  res_safety=res_etape1=res_etape2=rep(NA,ntrials) #nombre d'échecs et succès pour chaque essai
  i=l=0
  sample_s=sample1_d1=sample2_d1=0
  sample_s2=sample1_d2=sample2_d2=0
  nb_eff=0 #nombre de réponses (succès) à chaque étape
  nb_ech=0 #nombre de non-réponses (échecs) à chaque étape
  nb_tox=0 #nombre de toxicités à chaque étape
  time=rep(NA,ntrials)
  timesafe=rep(NA,ntrials) #stocke la durée du 1er safety lead-in (utile quand désescalade de dose)
  
  
  set.seed(3003)
  for (n in 1:ntrials) {
    patient=cbind(t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox))),cumsum(c(0,rep(2,n_sl+N-1)))) #temps fixe: cumsum(c(0,rep(2,n_sl+N-1)
    #patient1=t(rmultinom(npat,1,c(peffnotox,pefftox,pnoeffnotox,pnoefftox)))
    #fixe=cumsum(c(0,rep(4,n_sl),rep(2,N-1)))
    #pois=cumsum(c(0,rep(4,n_sl),rpois(N-1,2)))
    #patient=cbind(patient1,pois)
    
    limite=limite2=0
    rep_safe=rep_echec=0
    
    #SAFETY LEAD-IN (DL1)
    for (i in 1:n_sl) {
      limite=limite + patient[,2][i] + patient[,4][i]
      rep_safe=rep_safe + patient[,1][i] + patient[,2][i]  
      rep_echec=rep_echec + patient[,3][i] + patient[,4][i]
      if (limite==2) {break()}
    }
    if (limite==2) {res_safety[n]<-"echec"} else {res_safety[n]<-"succes DL1"} #0=échec ; 1=succès (passage à l'étape 1)
    if (res_safety[n]=="echec") {
      sample_s=sample_s+i+1
      nb_eff=nb_eff+rep_safe+patient[,1][i+1] + patient[,2][i+1]
      nb_ech=nb_ech+rep_echec+patient[,3][i+1] + patient[,4][i+1]
      nb_tox=nb_tox+limite+patient[,2][i+1] + patient[,4][i+1]
    } else {
      sample_s=sample_s+i
      nb_eff=nb_eff+rep_safe
      nb_ech=nb_ech+rep_echec
      nb_tox=nb_tox+limite
    }
    #print(paste("essai",n))
    #print(paste(i,"patients",limite,"toxicity",rep_safe,"efficacités"))
    #print(paste("Résultat s_l:",res_safety[n]))
    
    if (res_safety[n]=="echec") {
      time[n]<-patient[i+1,5]+t_obs
      timesafe[n]<-patient[i+1,5]+t_obs
      #cat("Durée de l'essai(saf):",time[n])
    }
    
    #SIMON ETAPE 1 (DL1)
    if (res_safety[n]=="succes DL1") {
      eff1= sum(patient[n_sl+1:n1,1]) + sum(patient[n_sl+1:n1,2])
      ech1= sum(patient[n_sl+1:n1,3]) + sum(patient[n_sl+1:n1,4])
      tox1= sum(patient[n_sl+1:n1,2]) + sum(patient[n_sl+1:n1,4])
      
      if (eff1>r1) {res_etape1[n]<-"succes1 DL1"} else {res_etape1[n]<-"echec"} #0=échec ; 1=succès (passage à l'étape 2)
      sample1_d1=sample1_d1+n1
      nb_eff=nb_eff+eff1
      nb_ech=nb_ech+ech1
      nb_tox=nb_tox+tox1
      #print(paste(n1,"patients",eff1,"efficacités",tox1,"toxicités"))
      #print(paste("Résultat étape1:",res_etape1[n]))
      
      if (res_etape1[n]=="echec") {
        if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) == limite-1) {time[n]<-patient[n1+n_sl,5]+2*t_obs}
        else {time[n]<-patient[n1+n_sl,5]+t_obs}
        #cat("Durée de l'essai(1):",time[n])
      }
      
      
      #SIMON ETAPE 2 (DL1)
      if (res_etape1[n]=="succes1 DL1") {
        
        eff2=sum(patient[n_sl+n1+1:(N-n1),1]) + sum(patient[n_sl+n1+1:(N-n1),2])
        ech2=sum(patient[n_sl+n1+1:(N-n1),3]) + sum(patient[n_sl+n1+1:(N-n1),4])
        tox2=sum(patient[n_sl+n1+1:(N-n1),2]) + sum(patient[n_sl+n1+1:(N-n1),4])
        
        if (eff1+eff2>R) {res_etape2[n]<-"succes2 DL1"} else {res_etape2[n]<-"echec"} #0=échec ; 1=succès (conclu d'efficacité)
        sample2_d1=sample2_d1+N-n1
        nb_eff=nb_eff+eff2
        nb_ech=nb_ech+ech2
        nb_tox=nb_tox+tox2
        #print(paste(N-n1,"patients",eff2,"efficacités",tox2,"toxicités"))
        #print(paste("Résultat étape2:",res_etape2[n]))
        
        if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) == limite-1) {
          if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) == r1) {time[n] <- patient[N+n_sl,5] + 3*t_obs}
          if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) != r1) {time[n] <- patient[N+n_sl,5] + 2*t_obs}
        }
        if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) != limite-1) {
          if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) == r1) {time[n] <- patient[N+n_sl,5] + 2*t_obs}
          if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) != r1) {time[n] <- patient[N+n_sl,5] + t_obs}
        }
      }
      
    }
    
    if (res_safety[n]=="echec") {
      
      #SAFETY LEAD-IN (DL-1)
      patient2=cbind(t(rmultinom(npat,1,c(peffnotox2,pefftox2,pnoeffnotox2,pnoefftox2))),cumsum(c(0,rep(2,n_sl+N-1)))) #temps fixe:cumsum(c(0,rep(2,n_sl+N-1)
      #patient_2=t(rmultinom(npat,1,c(peffnotox2,pefftox2,pnoeffnotox2,pnoefftox2)))
      #fixe=cumsum(c(0,rep(4,n_sl),rep(2,N-1)))
      #pois2=cumsum(c(0,rep(4,n_sl),rpois(N-1,2)))
      #patient2=cbind(patient_2,pois2)
      
      limite2=0
      rep_safe=rep_echec=0
      
      for (l in 1:n_sl) {
        limite2 = limite2 + patient2[,2][l] + patient2[,4][l]
        rep_safe=rep_safe + patient2[,1][l] + patient2[,2][l]  
        rep_echec=rep_echec + patient2[,3][l] + patient2[,4][l]
        if (limite2==2) {break()}
      }
      if (limite2==2) {res_safety[n]<-"echec"} else {res_safety[n]<-"succes DL-1"}
      if (res_safety[n]=="echec") {
        sample_s2=sample_s2+l+1
        nb_eff=nb_eff+rep_safe+patient[,1][l+1] + patient[,2][l+1]
        nb_ech=nb_ech+rep_echec+patient[,3][l+1] + patient[,4][l+1]
        nb_tox=nb_tox+limite+patient[,2][l+1] + patient[,4][l+1]
      } else {
        sample_s2=sample_s2+l
        nb_eff=nb_eff+rep_safe
        nb_ech=nb_ech+rep_echec
        nb_tox=nb_tox+limite
      }
      #print(paste("essai",n))
      #print(paste(l,"patients",limite,"toxicity",rep_safe,"efficacités"))
      #print(paste("Résultat s_l:",res_safety[n]))
      
      if (res_safety[n]=="echec") {
        time[n]<-timesafe[n] + patient2[l+1,5]+t_obs
        #cat("Durée de l'essai(saf):",time[n])
      }
      
      #SIMON ETAPE 1 (DL-1)
      if (res_safety[n]=="succes DL-1") {
        eff1=sum(patient2[n_sl+1:n1,1]) + sum(patient2[n_sl+1:n1,2])
        ech1=sum(patient2[n_sl+1:n1,3]) + sum(patient2[n_sl+1:n1,4])
        tox1=sum(patient2[n_sl+1:n1,2]) + sum(patient2[n_sl+1:n1,4])
        
        if (eff1>r1) {res_etape1[n]<-"succes1 DL-1"} else {res_etape1[n]<-"echec"} #0=échec ; 1=succès (passage à l'étape 2)
        sample1_d2=sample1_d2+n1
        nb_eff=nb_eff+eff1
        nb_ech=nb_ech+ech1
        nb_tox=nb_tox+tox1
        #print(paste(n1,"patients",eff1,"efficacités",tox1,"toxicités"))
        #print(paste("Résultat étape1:",res_etape1[n]))
        
        if (res_etape1[n]=="echec") {
          if (sum(patient2[1:(n_sl-1),2])+sum(patient2[1:(n_sl-1),4]) == limite-1) {time[n]<-timesafe[n]+patient2[n1+n_sl,5]+2*t_obs}
          else {time[n]<-timesafe[n]+patient2[n1+n_sl,5]+t_obs}
          #cat("Durée de l'essai(1):",time[n])
        }
        
        
        #SIMON ETAPE 2 (DL-1)
        if (res_etape1[n]=="succes1 DL-1") {
          eff2=sum(patient2[n_sl+n1+1:(N-n1),1]) + sum(patient2[n_sl+n1+1:(N-n1),2])
          ech2=sum(patient2[n_sl+n1+1:(N-n1),3]) + sum(patient2[n_sl+n1+1:(N-n1),4])
          tox2=sum(patient2[n_sl+n1+1:(N-n1),2]) + sum(patient2[n_sl+n1+1:(N-n1),4])
          
          if (eff1+eff2>R) {res_etape2[n]<-"succes2 DL-1"} else {res_etape2[n]<-"echec"} #0=échec ; 1=succès (conclu d'efficacité)
          sample2_d2=sample2_d2+N-n1
          nb_eff=nb_eff+eff2
          nb_ech=nb_ech+ech2
          nb_tox=nb_tox+tox2
          #print(paste(N-n1,"patients",eff2,"efficacités",tox2,"toxicités"))
          #print(paste("Résultat étape2:",res_etape2[n]))
          
          if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) == limite-1) {
            if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) == r1) {time[n] <- timesafe[n]+patient[N+n_sl,5] + 3*t_obs}
            if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) != r1) {time[n] <- timesafe[n]+patient[N+n_sl,5] + 2*t_obs}
          }
          if (sum(patient[1:(n_sl-1),2])+sum(patient[1:(n_sl-1),4]) != limite-1) {
            if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) == r1) {time[n] <- timesafe[n]+patient[N+n_sl,5] + 2*t_obs}
            if (sum(patient[n_sl+1:(n1-1),1])+sum(patient[n_sl+1:(n1-1),2]) != r1) {time[n] <- timesafe[n]+patient[N+n_sl,5] + t_obs}
          }
          
        }
      }
    }
  }
  
  
  #effectif attendu
  grid$ess[g] <- (sample_s+sample_s2+sample1_d1+sample2_d1+sample1_d2+sample2_d2)/ntrials
  
  #nb d'efficacités 
  grid$moy_eff[g] <- nb_eff/ntrials
  
  #nb d'échecs
  grid$moy_ech[g] <- nb_ech/ntrials
  
  #nb de toxicités 
  grid$moy_tox[g] <- nb_tox/ntrials
  
  #durée de l'essai
  grid$moy_time[g] <- mean(time)
  
  
  #% conclu d'efficacité
  if (all(is.na(res_etape2[res_etape2=="succes2 DL-1"]))) {
    grid$succes[g] <- ((table(res_etape2[res_etape2=="succes2 DL1"]))*100)/ntrials
  } else {
    if (all(is.na(res_etape2[res_etape2=="succes2 DL1"]))) {
      grid$succes[g] <- ((table(res_etape2[res_etape2=="succes2 DL-1"]))*100)/ntrials
    } else {
      if (all(is.na(res_etape2[res_etape2=="succes2 DL1"])) && all(is.na(res_etape2[res_etape2=="succes2 DL-1"]))) {
        grid$succes[g] <- 0
      } else {
        grid$succes[g] <- ((table(res_etape2[res_etape2=="succes2 DL1"]) + table(res_etape2[res_etape2=="succes2 DL-1"]))*100)/ntrials
      }
    }
  }
  
  #% arrêt précoce
  if (all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$early_ter[g] <- ((table(res_safety[res_safety=="echec"]))*100)/ntrials
  } else {
    if (all(is.na(res_safety[res_safety=="echec"]))) {
      grid$early_ter[g] <- ((table(res_etape1[res_etape1=="echec"]))*100)/ntrials
    } else {
      if (all(is.na(res_etape1[res_etape1=="echec"])) && all(is.na(res_safety[res_safety=="echec"]))) {
        grid$early_ter[g] <- 0
      } else {
        grid$early_ter[g] <- ((table(res_etape1[res_etape1=="echec"])+table(res_safety[res_safety=="echec"]))*100)/ntrials
      }
    }
  }
  
  #% arrêt pour tox
  if (all(is.na(res_safety[res_safety=="echec"]))) {
    grid$stop_tox[g] <- 0
  } else {
    grid$stop_tox[g] <- (table(res_safety[res_safety=="echec"])*100)/ntrials
  }
  
  
  #nb de succès et échecs à chaque étape
  if (all(is.na(res_safety[res_safety=="succes DL-1"])) && all(is.na(res_safety[res_safety=="succes DL1"]))) {
    grid$res_safe_succes[g] <- 0
  } else {
    grid$res_safe_succes[g] <- table(res_safety[res_safety=="succes DL-1"])+table(res_safety[res_safety=="succes DL1"])
  }
  if (all(is.na(res_safety[res_safety=="echec"]))) {
    grid$res_safe_echecs[g] <- 0
  } else {
    grid$res_safe_echecs[g] <- table(res_safety[res_safety=="echec"])
  }
  
  if (all(is.na(res_etape1[res_etape1=="succes1 DL-1"]))) {
    grid$res1_succes[g] <- table(res_etape1[res_etape1=="succes1 DL1"])
  } else {
    if (all(is.na(res_etape1[res_etape1=="succes1 DL1"]))) {
      grid$res1_succes[g] <- table(res_etape1[res_etape1=="succes1 DL-1"])
    } else {
      if (all(is.na(res_etape1[res_etape1=="succes1 DL1"])) && all(is.na(res_etape1[res_etape1=="succes1 DL-1"]))) {
        grid$res1_succes[g] <- 0
      } else {
        grid$res1_succes[g] <- table(res_etape1[res_etape1=="succes1 DL1"]) + table(res_etape1[res_etape1=="succes1 DL-1"])
      }
    }
  }
  if (all(is.na(res_etape1[res_etape1=="echec"]))) {
    grid$res1_echecs[g] <- 0
  } else {
    grid$res1_echecs[g] <- table(res_etape1[res_etape1=="echec"])
  }
  
  if (all(is.na(res_etape2[res_etape2=="succes2 DL-1"]))) {
    grid$res2_succes[g] <- table(res_etape2[res_etape2=="succes2 DL1"])
  } else {
    if (all(is.na(res_etape2[res_etape2=="succes2 DL1"]))) {
      grid$res2_succes[g] <- table(res_etape2[res_etape2=="succes2 DL-1"])
    } else {
      if (all(is.na(res_etape2[res_etape2=="succes2 DL1"])) && all(is.na(res_etape2[res_etape2=="succes2 DL-1"]))) {
        grid$res2_succes[g] <- 0
      } else {
        grid$res2_succes[g] <- table(res_etape2[res_etape2=="succes2 DL1"]) + table(res_etape2[res_etape2=="succes2 DL-1"])
      }
    }
  }
  if (all(is.na(res_etape2[res_etape2=="echec"]))) {
    grid$res2_echecs[g] <- 0
  } else {
    grid$res2_echecs[g] <- table(res_etape2[res_etape2=="echec"])
  }
  
  
  
}


writexl::write_xlsx(grid, "C:\\util\\marie\\2022-2023\\stage\\R\\résultats\\safety2doses_Simon_tpsfixe_2sem.xlsx")


table(res_safety)
table(res_etape1)
table(res_etape2)

(sample_s+sample_s2+
    sample1_d1+sample2_d1+
    sample1_d2+sample2_d2)/ntrials

#nombre de réponses/efficacités
nb_eff/ntrials

#nombre de non réponses/non-efficacités
nb_ech/ntrials

#nombre de toxicités
nb_tox/ntrials

#durée de l'essai
mean(time)
table(time)