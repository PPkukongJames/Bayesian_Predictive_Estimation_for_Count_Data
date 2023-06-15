library("DWreg")
library("cubature")
library("spatstat")
comb <- function(n, x) {
  return(factorial(n) / factorial(n-x) / factorial(x))
}
MLE_Poi <- function(X,data){
  sum_d <- 0
  for (i in 1:length(X)) {
    sum_d<-sum_d+(X[i]*data[i])
  }
  sum_d<-sum_d/sum(data)
  return(sum_d)
}
MLE_Geo <- function(X,data){
  sum_d<-0
  for (i in 1:length(X)){## 0-6##1-3397 1-1065 1-369 1-116 1-41 1-8 1-4
      sum_d<-sum_d+(X[i]*data[i])
  }
  n <- sum(data)
  return(n/(sum_d+n))
}

MLE_NB <- function(X,data,r){
  sum_d<-0
  for (i in 1:length(X)){## 0-6##1-3397 1-1065 1-369 1-116 1-41 1-8 1-4
    sum_d<-sum_d+(X[i]*data[i])
  }
  return((sum(data)*r)/(sum_d+(sum(data)*r)))
}

GOF <- function(X,data,Seta,r,D){
  E<-rep(c(0),each=length(X))
  P<-rep(c(0),each=length(X))
  for (i in 1:length(X)-1){
    if(D=="Poi"){
      P[i]<- exp(-1*Seta)*(Seta^X[i])/factorial(X[i])
    }else if(D=="Geo"){
      P[i]<- Seta*(1-Seta)^X[i]
    }else if(D=="NB"){
      P[i]<-comb(r+X[i]-1,X[i])*(Seta^r)*(1-Seta)^(X[i])
    }else if(D=="WB"){#Seta = q , r = B
      P[i]<- (Seta^(X[i]^r))-(Seta^((X[i]+1)^r))
    }
  }
  P[length(X)] = 1-(sum(P))
  for (i in 1:length(X)){
    E[i]<-sum(data)*P[i]
  }
  Chi<- 0
  for (i in 1:length(X)){
    Chi<-Chi+((data[i]-E[i])^2)/E[i]
  }
  return(Chi)
}


X <- c(0,1,2,3,4,5)
data <- c(15,2,1,1,2,1)

ALL_parameter <- rep(c(0),each=4)
ALL_test <- rep(c(0),each=4)
ALL_NB <- rep(c(0),each=min(data))
###Poison
p_Poi <- MLE_Poi(X,data)
ALL_parameter[1] <- GOF(X,data,p_Poi,0,"Poi")
ALL_test[1] <- qchisq(p=0.95,df=length(X)-1)
 
###Geometrix
p_Geo <- MLE_Geo(X,data)
p_Geo
ALL_parameter[2] <- GOF(X,data,p_Geo,0,"Geo")
ALL_test[2] <- qchisq(p=0.95,df=length(X)-1)

###NB

r=min(data)
p_NB <- MLE_NB(X,data,r)
GOFN<-GOF(X,data,p_NB,r,"NB")


###wibul
for_WB <-rep(c(0),each=length(sum(data)))
k=1
for(i in 1:length(X)){
  for(j in 1:data[i]){
    for_WB[k]<-X[i]
    k<-k+1
  }
}
length(for_WB)
Weibull <- dw.parest(for_WB)
q <- Weibull$q
B <- Weibull$beta
ALL_parameter[3] <- GOF(X,data,q,B,"WB")
ALL_test[3] <- qchisq(p=0.95,df=length(X)-2)
#====================================================
#====================================================
# ALL_parameter [Poisol,Geometric,Wibul]
n <- sum(data)
sumdx<-0
for (i in 1:length(X) ){
  sumdx <- sumdx+(X[i]*data[i])
}

#############################Poisson#######################
#prior 1 y^-c
c_increase <- 1
c_end <- 500
len_c <- c_end/c_increase
c_data <- rep(c(0),each=len_c)
Poi_prior_1 <- rep(c(0),each=len_c)
for(i in 1:len_c){
  c_data[i]<-i*c_increase
}
# c variable is rows
Poi_1_c <-rep(c(0),each=len_c)
Poi_1_alpha <- rep(c(0),each=len_c)
Poi_1_beta <-rep(c(0),each=len_c)
for (i in 1:len_c){
  Poi_1_c[i] <- c_data[i]
  alpha <- sumdx-Poi_1_c[i]+1
  Poi_1_alpha[i] <- alpha
  beta <- n
  Poi_1_beta[i] <- beta
  Poi_prior_1[i] <- pgamma(p_Poi,alpha,beta)
}
data_Poi_1 <- data.frame (
  a = Poi_1_alpha,
  b = Poi_1_beta,
  c = Poi_1_c,
  p = Poi_prior_1
)
round_Poi_1 <- len_c
Poi_prior_1
#prior 2 gamma

increase <- 0.2#increase 0.2 per round (gamma)
loop <-1
while (TRUE){
  if(p_Poi/(increase*loop) < increase){
    loop+1
    break
  }
  loop <- loop+1
} 
loop
a<-rep(c(0),each=loop)
b<-rep(c(0),each=loop)
Poi_prior_2 <- rep(c(0),each=loop)
for (i in 1:loop){
  a[i]<-increase*i
  b[i]<-p_Poi/a[i]
  Poi_prior_2[i]<- pgamma(p_Poi,sumdx+a[i],n+b[i])
}

round_Poi_2 <- loop
data_Poi_2 <- data.frame (
  a = a,
  b = b,
  p = Poi_prior_2
)

####################Weibull##############################
dataWei<-data.frame (
  a = c(0),
  b = c(0),
  p = c(0)
)

for(i in 1:length(d_X)){
  wb=i
  wa<-B/wb
  wa
  f.wuibull_1<-function(p){
    exp(sum( log((p[1]^(for_WB^p[2]))-(p[1]^((for_WB+1)^p[2]))) ) )* dgamma(p[2],wa,wb)
  }
  constwb_1 <- hcubature(f.wuibull_1,lower=c(0,0),upper=c(1,Inf))$integral
  
  dataWei <- rbind( dataWei, c(wa,wb,exp(sum( log((q^(for_WB^B))-(q^((for_WB+1)^B))) ) )* dgamma(B,wa,wb)/constwb_1))
  dataWei <- subset(dataWei,p<=1)
}
tempWei <- dataWei[dataWei$p == max(dataWei$p),]
wa <- tempWei$a
wb <- tempWei$b
wa
wb
#############################Geometric########################
#prior 1 Beta
startBeta <- sumdx
i<-1
if(startBeta<0){
  startBeta <- 1
}
while(n/(10**i)>10){
  i<-i+1
}
increase <- c((10**i),(10**i)/2,(10**i)/5)

dataG_1<-data.frame (
  a = c(0),
  b = c(0),
  p = c(0)
)
dataG_2<-data.frame (
  a = c(0),
  b = c(0),
  p = c(0)
)
dataG_3<-data.frame (
  a = c(0),
  b = c(0),
  p = c(0)
)
b_1 <- startBeta
b_2 <- startBeta
b_3 <- startBeta
for(i in 1:100) {
  dataG_1 <- rbind( dataG_1, c( ((p_Geo*b_1)/(1-p_Geo)) , b_1 , pbeta( p_Geo , ((p_Geo*b_1)/(1-p_Geo)) , b_1 ) ) )
  dataG_2 <- rbind( dataG_2, c( ((p_Geo*b_2)/(1-p_Geo)) , b_2 , pbeta( p_Geo , ((p_Geo*b_2)/(1-p_Geo)) , b_2 ) ) )
  dataG_3 <- rbind( dataG_3, c( ((p_Geo*b_3)/(1-p_Geo)) , b_3 , pbeta( p_Geo , ((p_Geo*b_3)/(1-p_Geo)) , b_3 ) ) )
  b_1 <- b_1+(increase[1]*i)
  b_2 <- b_2+(increase[2]*i)
  b_3 <- b_3+(increase[3]*i)
}

dataG_1 = subset(dataG_1,a>sum(data))
dataG_2 = subset(dataG_2,a>sum(data))
dataG_3 = subset(dataG_3,a>sum(data))

###################Negative Binomial
startBeta <- sumdx
ALL_NB_1 =  data.frame(
  a = c(0),
  b = c(0),
  p = c(0),
  r = c(0)
)
ALL_NB_2 =  data.frame(
  a = c(0),
  b = c(0),
  p = c(0),
  r = c(0)
)
ALL_NB_3 =  data.frame(
  a = c(0),
  b = c(0),
  p = c(0),
  r = c(0)
)
i<-1
if(startBeta<0){
  startBeta <- 1
}
while(n/(10**i)>10){
  i<-i+1
}

increase <- c((10**i),(10**i)/2,(10**i)/5)
b_1 <- startBeta
b_2 <- startBeta
b_3 <- startBeta
for(i in 1:100) {
  ALL_NB_1 <- rbind( ALL_NB_1, c( ((p_NB*b_1)/(1-p_NB)) , b_1 , dbeta( p_NB , ((p_NB*b_1)/(1-p_NB)) , b_1 ) , r) )
  ALL_NB_2 <- rbind( ALL_NB_2, c( ((p_NB*b_2)/(1-p_NB)) , b_2 , dbeta( p_NB , ((p_NB*b_2)/(1-p_NB)) , b_2 ) , r) )
  ALL_NB_3 <- rbind( ALL_NB_3, c( ((p_NB*b_3)/(1-p_NB)) , b_3 , dbeta( p_NB , ((p_NB*b_3)/(1-p_NB)) , b_3 ) , r) )
  b_1 <- b_1+(increase[1]*i)
  b_2 <- b_2+(increase[2]*i)
  b_3 <- b_3+(increase[3]*i)
}
ALL_NB_1 = subset(ALL_NB_1,a>n*r)
ALL_NB_2 = subset(ALL_NB_2,a>n*r)
ALL_NB_3 = subset(ALL_NB_3,a>n*r)

#prior 2 Uniform

NB_prior_2 <- pbeta(p_NB,(r*n)+1,sumdx+1)
NB_prior_2




a<-c(0,0,0)
b<-c(0,0,0)
for(i in 2:length(dataG_1$p)){
  if(as.integer(dataG_1$p[i]*10000)==as.integer(dataG_1$p[i-1]*10000)){
    a[1]<-dataG_1$a[i-1]
    b[1]<-dataG_1$b[i-1]
    
    break
  }
}
for(i in 2:length(dataG_2$p)){
  if(as.integer(dataG_2$p[i]*10000)==as.integer(dataG_2$p[i-1]*10000)){
    a[2]<-dataG_2$a[i-1]
    b[2]<-dataG_2$b[i-1]
    break
  }
}
for(i in 2:length(dataG_3$p)){
  if(as.integer(dataG_3$p[i]*10000)==as.integer(dataG_3$p[i-1]*10000)){
    a[3]<-dataG_3$a[i-1]
    b[3]<-dataG_3$b[i-1]
    break
  }
}
ga <- a[1]
gb <- b[1]
for (i in 2:3){
  if(b[i]<gb){
    gb<-b[i]
    ga<-a[i]
  }
}
gb
ga


a<-c(0,0,0)
b<-c(0,0,0)
for(i in 2:length(ALL_NB_1$p)){
  if(as.integer(ALL_NB_1$p[i]*10000)==as.integer(ALL_NB_1$p[i-1]*10000)){
    a[1]<-ALL_NB_1$a[i-1]
    b[1]<-ALL_NB_1$b[i-1]
    
    break
  }
}
for(i in 2:length(ALL_NB_2$p)){
  if(as.integer(ALL_NB_2$p[i]*10000)==as.integer(ALL_NB_2$p[i-1]*10000)){
    a[2]<-ALL_NB_2$a[i-1]
    b[2]<-ALL_NB_2$b[i-1]
    break
  }
}
for(i in 2:length(ALL_NB_3$p)){
  if(as.integer(ALL_NB_3$p[i]*10000)==as.integer(ALL_NB_3$p[i-1]*10000)){
    a[3]<-ALL_NB_3$a[i-1]
    b[3]<-ALL_NB_3$b[i-1]
    break
  }
}
NBa <- a[1]
NBb <- b[1]
for (i in 2:3){
  if(b[i]<NBa){
    NBb<-b[i]
    NBa<-a[i]
  }
}
NBb
NBa


#===========Poisson distribution===========
m=0
Poi_1<-rep(c(0),each=length(X)+m)
Poi_2<-rep(c(0),each=length(X)+m)
Geo_1<-rep(c(0),each=length(X)+m)
Geo_2<-rep(c(0),each=length(X)+m)
NB_1<-rep(c(0),each=length(X)+m)
NB_2<-rep(c(0),each=length(X)+m)
Weibull_1<-rep(c(0),each=length(X)+m)
Weibull_2<-rep(c(0),each=length(X)+m)

for (i in 1:length(X)+m){
  if(X[1]==0){
    k=i-1
  }
  data_Poi_2
  c=data_Poi_1$c[3]
  pa = data_Poi_2$a[3]
  pb = data_Poi_2$b[3]
  f.Poi_pri_1<-function(ramda){
    dpois(k,ramda)*pgamma(ramda,sumdx-c+1,n)
  }
  f.Poi_pri_2<-function(ramda){
    dpois(k,ramda)*pgamma(ramda,sumdx+pa,n+pb)
  }
  Poi_1[i]<-hcubature(f.Poi_pri_1,lower=0,upper=Inf)$integral
  Poi_2[i]<-hcubature(f.Poi_pri_2,lower=0,upper=Inf)$integral
  #===========Geometrix distribution===========
  f.Geo_pri_1<-function(p){
    dgeom(k,p)*dbeta(p,n+ga,sumdx+gb)
  }
  f.Geo_pri_2<-function(p){
      dgeom(k,p)*dbeta(p,n,sumdx)
  }
  Geo_1[i]<-hcubature(f.Geo_pri_1,lower=0,upper=1)$integral
  Geo_2[i]<-hcubature(f.Geo_pri_2,lower=0,upper=1)$integral
  #===========Negative Binomial distribution===========
  r=min(data)
  f.NB_pri_1<-function(p){
    dnbinom(k,r,p)*pbeta(p,r*n+NBa,sumdx+NBb)
  }
  f.NB_pri_2<-function(p){
    dnbinom(k,r,p)*pbeta(p,n*r + 1,sumdx + 1)
  }
  NB_1[i]<-hcubature(f.NB_pri_1,lower=0,upper=1)$integral
  NB_2[i]<-hcubature(f.NB_pri_2,lower=0,upper=1)$integral
  #===========Discrete Weibull distribution===========
  f.wuibull_1<-function(p){
    exp(sum(  log((p[1]^(for_WB^p[2]))-(p[1]^((for_WB+1)^p[2]))   )))* dgamma(p[2],wa,wb)
  }
  constwb_1 <- hcubature(f.wuibull_1,lower=c(0,0),upper=c(1,Inf))$integral
  f.wuibull_1_1<-function(p){
    ( (p[1]^(k^p[2]))-(p[1]^((k+1)^p[2])) ) * exp(sum( log((p[1]^(for_WB^p[2]))-(p[1]^((for_WB+1)^p[2])) ))) * dgamma(p[2],wa,wb) /constwb_1
  }
  Weibull_1[i] <- hcubature(f.wuibull_1_1,lower=c(0,0),upper=c(1,Inf))$integral
  ###############################################
  f.wuibull_2<-function(p){
    exp(sum( log((p[1]^(for_WB^p[2]))-(p[1]^((for_WB+1)^p[2])) ))) * (1/sqrt(p[1]*(1-p[1])*p[2]))
  }
  constwb_2 <- hcubature(f.wuibull_2,lower=c(0,0),upper=c(1,Inf))$integral
  f.wuibull_2_1<-function(p){
    ( (p[1]^(k^p[2]))-(p[1]^((k+1)^p[2])) )* exp(sum( log((p[1]^(for_WB^p[2]))-(p[1]^((for_WB+1)^p[2])) ))) * (1/sqrt(p[1]*(1-p[1])*p[2]))/constwb_2
  }
  Weibull_2[i] <- hcubature(f.wuibull_2_1,lower=c(0,0),upper=c(1,Inf))$integral
}

Poi_1
Poi_2
Geo_1
Geo_2
NB_1
NB_2
Weibull_1
Weibull_2

#####################################MOVOE############################
m=0
Poi_1_re<-rep(c(0),each=length(X)+m)
Poi_2_re<-rep(c(0),each=length(X)+m)
Geo_1_re<-rep(c(0),each=length(X)+m)
Geo_2_re<-rep(c(0),each=length(X)+m)
NB_1_re<-rep(c(0),each=length(X)+m)
NB_2_re<-rep(c(0),each=length(X)+m)
Weibull_1_re<-rep(c(0),each=length(X)+m)
Weibull_2_re<-rep(c(0),each=length(X)+m)

for (i in 1:length(X)+m){
  if(X[1]==0){
    k=i-1
  }
  for_WB_re <- for_WB
  for(j in 1:length(for_WB_re)){
    if(for_WB_re[j]==X[i]){
      for_WB_re<-for_WB_re[-j]
      break
    }
  }
  sumdx_re <- sum(for_WB_re)
  
  c=data_Poi_1$c[3]
  pa = data_Poi_2$a[3]
  pb = data_Poi_2$b[3]
  f.Poi_pri_1<-function(ramda){
    dpois(k,ramda)*pgamma(ramda,sumdx_re-c+1,n)
  }
  f.Poi_pri_2<-function(ramda){
    dpois(k,ramda)*pgamma(ramda,sumdx_re+pa,n+pb)
  }
  Poi_1_re[i]<-hcubature(f.Poi_pri_1,lower=0,upper=Inf)$integral
  Poi_2_re[i]<-hcubature(f.Poi_pri_2,lower=0,upper=Inf)$integral
  #===========Geometrix distribution===========
  f.Geo_pri_1<-function(p){
    dgeom(k,p)*dbeta(p,n+ga,sumdx_re+gb)
  }
  f.Geo_pri_2<-function(p){
    dgeom(k,p)*dbeta(p,n,sumdx_re)
  }
  Geo_1_re[i]<-hcubature(f.Geo_pri_1,lower=0,upper=1)$integral
  Geo_2_re[i]<-hcubature(f.Geo_pri_2,lower=0,upper=1)$integral
  #===========Negative Binomial distribution===========
  r=min(data)
  f.NB_pri_1<-function(p){
    dnbinom(k,r,p)*pbeta(p,r*n+NBa,sumdx_re+NBb)
  }
  f.NB_pri_2<-function(p){
    dnbinom(k,r,p)*pbeta(p,n*r + 1,sumdx_re + 1)
  }
  NB_1_re[i]<-hcubature(f.NB_pri_1,lower=0,upper=1)$integral
  NB_2_re[i]<-hcubature(f.NB_pri_2,lower=0,upper=1)$integral
  #===========Discrete Weibull distribution===========
  f.wuibull_1<-function(p){
    exp(sum( log((p[1]^(for_WB_re^p[2]))-(p[1]^((for_WB_re+1)^p[2]))) ) )* dgamma(p[2],wa,wb)
  }
  constwb_1 <- hcubature(f.wuibull_1,lower=c(0,0),upper=c(1,Inf))$integral
  f.wuibull_1_1<-function(p){
    ( (p[1]^(k^p[2]))-(p[1]^((k+1)^p[2])) ) * exp(sum( log((p[1]^(for_WB_re^p[2]))-(p[1]^((for_WB_re+1)^p[2])) ))) * dgamma(p[2],wa,wb) /constwb_1
  }
  Weibull_1_re[i] <- hcubature(f.wuibull_1_1,lower=c(0,0),upper=c(1,Inf))$integral

  ###############################################
  f.wuibull_2<-function(p){
    exp(sum( log((p[1]^(for_WB_re^p[2]))-(p[1]^((for_WB_re+1)^p[2]))) )) * (1/sqrt(p[1]*(1-p[1])*p[2]))
  }
  constwb_2 <- hcubature(f.wuibull_2,lower=c(0,0),upper=c(1,Inf))$integral
  f.wuibull_2_1<-function(p){
    ( (p[1]^(k^p[2]))-(p[1]^((k+1)^p[2])) )* exp(sum( log((p[1]^(for_WB_re^p[2]))-(p[1]^((for_WB_re+1)^p[2]))) )) * (1/sqrt(p[1]*(1-p[1])*p[2]))/constwb_2
  }
  Weibull_2_re[i] <- hcubature(f.wuibull_2_1,lower=c(0,0),upper=c(1,Inf))$integral
}

d_X <-for_WB

MAVOE_mean<-function(data,d_X,p){
  total <- 0
  for(i in 1:length(d_X)){
    temp_mean<-0
    for(j in 1:length(d_X)){
      if (i!=j){
        temp_p<-0
        if(d_X[1]==0){
          temp_p<-d_X[j]+1
        }else{
          temp_p<-d_X[j]
        }
        temp_mean<-temp_mean+(p[temp_p]*d_X[j])
      }
    }
    total <- total+abs(d_X[i]-temp_mean)
  }
  total/length(d_X)
}

MAVOE_mean(data,d_X,Poi_1_re)
MAVOE_mean(data,d_X,Poi_2_re)
MAVOE_mean(data,d_X,Geo_1_re)
MAVOE_mean(data,d_X,Geo_2_re)
MAVOE_mean(data,d_X,NB_1_re)
MAVOE_mean(data,d_X,Weibull_1_re)
MAVOE_mean(data,d_X,Weibull_2_re)

median_Geo_1 <- 0
median_Geo_2 <- 0
median_Wei_1 <- 0
median_Wei_2 <- 0
for(i in 1:length(d_X)){
  Geo_1_re_median<-rep(c(0),each=length(X))
  Geo_2_re_median<-rep(c(0),each=length(X))
  Weibull_1_re_median <-rep(c(0),each=length(X))
  Weibull_2_re_median <-rep(c(0),each=length(X))
  X_re <- d_X[-i]
  temp_dx <- d_X[i]
  n=length(X_re)
  sumX_median <- sum(X_re)
  for(j in 1:length(X)){
    k<-X[j]
  ####################################
    #===========Geometrix distribution===========
    f.Geo_pri_1<-function(p){
      dgeom(k,p)*dbeta(p,n+ga,sumX_median+gb)
    }
    f.Geo_pri_2<-function(p){
      dgeom(k,p)*dbeta(p,n,sumX_median)
    }
    #===========Discrete Weibull distribution===========
    f.wuibull_1<-function(p){
      exp(sum( log((p[1]^(X_re^p[2]))-(p[1]^((X_re+1)^p[2]))) ) )* dgamma(p[2],wa,wb)
    }
    constwb_1 <- hcubature(f.wuibull_1,lower=c(0,0),upper=c(1,Inf))$integral
    f.wuibull_1_1<-function(p){
      ( (p[1]^(k^p[2]))-(p[1]^((k+1)^p[2])) ) * exp(sum( log((p[1]^(X_re^p[2]))-(p[1]^((X_re+1)^p[2])) ))) * dgamma(p[2],wa,wb) /constwb_1
    }
    ###############################################
    f.wuibull_2<-function(p){
      exp(sum( log((p[1]^(X_re^p[2]))-(p[1]^((X_re+1)^p[2]))) )) * (1/sqrt(p[1]*(1-p[1])*p[2]))
    }
    constwb_2 <- hcubature(f.wuibull_2,lower=c(0,0),upper=c(1,Inf))$integral
    f.wuibull_2_1<-function(p){
      ( (p[1]^(k^p[2]))-(p[1]^((k+1)^p[2])) )* exp(sum( log((p[1]^(X_re^p[2]))-(p[1]^((X_re+1)^p[2]))) )) * (1/sqrt(p[1]*(1-p[1])*p[2]))/constwb_2
    }
    Geo_1_re_median[j]<-hcubature(f.Geo_pri_1,lower=0,upper=1)$integral
    Geo_2_re_median[j]<-hcubature(f.Geo_pri_2,lower=0,upper=1)$integral
    Weibull_1_re_median[j] <- hcubature(f.wuibull_1_1,lower=c(0,0),upper=c(1,Inf))$integral
    Weibull_2_re_median[j] <- hcubature(f.wuibull_2_1,lower=c(0,0),upper=c(1,Inf))$integral
  }
  median_Geo_1 <- median_Geo_1+abs(temp_dx - weighted.median(X,Geo_1_re_median))
  median_Geo_2 <- median_Geo_2+abs(temp_dx - weighted.median(X,Geo_2_re_median))
  median_Wei_1 <- median_Wei_1+abs(temp_dx - weighted.median(X,Weibull_1_re_median))
  median_Wei_2 <- median_Wei_2+abs(temp_dx - weighted.median(X,Weibull_2_re_median))
}
median_Geo_1 <- median_Geo_1/length(d_X)
median_Geo_2 <- median_Geo_2/length(d_X)
median_Wei_1 <- median_Wei_1/length(d_X)
median_Wei_2 <- median_Wei_2/length(d_X)

median_Geo_1
median_Geo_2
median_Wei_1
median_Wei_2
