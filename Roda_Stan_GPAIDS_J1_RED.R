##############################################################
# Dados de AIDS
#setwd("C:/Users/Public/2022/Code_generalized Poisson/GP_BAYES")
dados =read.csv("sobrAIDS.csv", sep = '\t',  header = T,na.strings = "NA")
#dados1=dados[-c(47,211),]
##############################################################################
y=as.numeric(dados$tempo/365.25)
x1=as.numeric(dados$idade)

dados$raca=as.numeric(dados$raca)
x2=factor(ifelse(dados$raca==0,0,1))
x3 = factor(dados$sexo)
x4 = factor(dados$prevoi)
status=dados$status
X=model.matrix(~1+x1+x4)
p=ncol(X)
n=nrow(X)
library(rstan)
first.dso <-stan_model('GP_GCRATE_Exp.stan',model_name = 'First Model')
dat <- list(
  n = n,
 p = p,
  y = y,
  status = status,
 X=X
)
first.stanfit<-sampling(first.dso, dat,chains = 4, iter = 10000,thin=5)
            
print(first.stanfit,, pars=c("beta", "lambda1", "gama"),digits=3)
     print(first.stanfit, pars=c("beta", "lambda1", "gama"), probs=c(.1,.5,.9))
       traceplot(first.stanfit, pars=c("lambda1", "gama", "beta")) 
pairs(first.stanfit, pars=c("beta", "lambda1", "gama"))
stan_plot(first.stanfit)
stan_plot(first.stanfit, point_est = "mean", show_density = TRUE, fill_color = "maroon")
Dad=extract(first.stanfit,pars=c("lambda1", "gama", "beta"))
dd=cbind(Dad$lambda1,Dad$gama,Dad$beta)
dim(dd)

####################################################################
## Log-verossimilhan�a
####################################################################
library(LambertW)
library(lamW)

fLogexpp =function(vp){
  
  R=nrow(vp)
  flog = matrix(0,nrow=n,ncol = R)
  invf=flog
  for (i in 1:n) {
    for(a in 1:R) {
      vpar=vp[a,]
      lam = vpar[1]        # dist exp por partes
      theta =vpar[2]        # theta (GP)
      beta = vpar[-c(1:2)]   # covariáveis
      vf = exp(lam)*exp(-exp(lam)*y[i])
      vF = 1-exp(-exp(lam)*y[i])
      linear=c(X[i,]%*%beta)
      mu=exp(linear)
      lambda=(1-theta)*mu
      ax = lambertW0(-theta*(1-vF)*exp(-theta))
      Spop <- exp(-(lambda/theta)*(ax+theta))     
      fpop <- -(lambda/theta)*(vf/(1-vF))*Spop*(ax/(1+ax)) 
      flog[i,a]=  status[i]*log(fpop)+(1-status[i])*log(Spop)
      invf[i,a]=1/exp(flog[i,a])    
    }
    
  }
  salida=list(ll=flog, iv=invf)
  return(salida)
}

Sal=fLogexpp(dd)
CPO=1/apply(Sal$iv, 1, mean)
LPML=sum(log(CPO))

#########################################################################
IHPD=function (x, alpha=0.05)
{
  n <- length(x)
  m <- max(1, ceiling(alpha * n))
  y <- sort(x)
  a <- y[1:m]
  b <- y[(n - m + 1):n]
  i <- order(b - a)[1]
  structure(c(a[i], b[i]), names = c("Lower Bound", "Upper Bound"))
}

Resumo <- function(x) {
  s <- c(mean(x, na.rm = T),median(x, na.rm = T),
         sqrt(var(x, na.rm = T)),IHPD(x))
  names(s)<-c("media","mediana" , "desvio Padrão", c("LI","LS"))
  return(s)
}
colnames(dd)=c("lambda", "gamma",paste("beta_", 0:(p-1), sep = ""))
Resultado=apply(dd,2,Resumo)
colnames(Resultado)=colnames(dd)
print(t(Resultado),3)
dJ1=t(Resultado)
library(xtable)
xtable(t(Resultado), digits = 3)
LPML
dc=rbind(t(dJ1),t(d47),t(d211),t(d47211))
xtable(dcj1digits = 3)
##############################################################
## RR
RR1=exp(dd[,5])
RR2=exp(dd[,4])
Resumo(RR1)
hist(RR1,prob=T, main="", ,xlab="Hazard ratio",cex.lab=1.5,cex.axis=1.5)
lines(density(RR1))
hist(RR2,prob=T, main="", ,xlab="Hazard ratio",cex.lab=1.5,cex.axis=1.5)
lines(density(RR2))
#############################################
###############################################
# GP Survival function
###############################################
Sgp=function(y,X,vpar) {
  lam = vpar[1]        # dist exp por partes
  theta =vpar[2]        # theta (GP)
  beta = vpar[-c(1:2)]   # covariáveis
  vF = 1-exp(-exp(lam)*y)
  linear=c(X%*%beta)
  mu=exp(linear)
  lambda=(1-theta)*mu
  ax = lambertW0(-theta*(1-vF)*exp(-theta))
  Spop <- exp(-(lambda/theta)*(ax+theta))     
  
}
vp=apply(dd,2,mean)
curve(-log(Sgp(x,p2,vp)),0,5,lty=1, col=1,
      lwd=3, xlab="Time (years)", ylab="Survival function")
curve(-log(Sgp(x,p1,vp)),0,5,add=T,lty=2,lwd=3,col=2)
curve(-log(Sgp(x,p4,vp)),0,20,add=F, lty=3, lwd=3)
curve(-log(Sgp(x,p3,vp)),0,20,add=T, lty=4, lwd=3)


#################################################
#    Residual
#################################################
Residual=function(vpar,X,y,status,Sp,nrep){
  spop=Sp(y,X,vpar)
  #nrep = 5
  set.seed(1234)
  mqresid = NULL
  u = status * (1 - spop) + (1 - status) * runif(length(y), 1 - spop)
  for (i in 1:nrep) {
    qresid = sort(qnorm(runif(u)))
    mqresid = cbind(mqresid, qresid)
  }
  qresid = apply(mqresid, 1, median)
  qresid
}
EB=apply(dd,2,mean)
Res_GP=Residual(EMVb,X,y,status,Sgp,5)
qqnorm(Res_GP, pch = 20, main = "", xlab = "N(0, 1) quantiles",
       ylab = "Normalized randomized quantile residuals", cex.lab = 1.5,
       cex.axis = 1.5, ylim = c(-3.5,3.5))
abline(0, 1, lty = 2, lwd = 2)
shapiro.test(Res_GP)
plot(density(Res_GP))
hist(Res_GP,prob=T,main="",xlab="Normalized randomized quantile residuals", cex.axis=1.5, cex.lab=1.5)
lines(density(Res_GP))

#############################################################

########################################################
# Diagnostic Analysis
N=length(y)
R=length(dd[,1])
K_L = numeric(N)
M_J = numeric(N)
L1= numeric(N)
CH= numeric(N)

for(i in 1:N){
  kl=0
  jk=0
  ll=0
  ch=0
  for (k in 1:R) {
    theta=dados[k,]
    z=CPO[i]/exp(Sal$ll[i,k])
    kl <- kl- log(z) ;
    jk=jk+log(z)*(z-1)
    ll=ll+0.5*abs(z-1)
    #ch=ch+(z-1)*(z-1)
    ch=ch+z*(1/z-1)*(1/z-1)
  }
  K_L[i]=kl/R
  M_J[i]=jk/R
  L1[i]=ll/R
  CH[i]=ch/R
}
###
p=0.75
cex.axis=1.25
cex.lab=1.25
par(mai=c(0.7,0.8,0.7,0.3))
par(mfrow=c(2,2))
par(mar=c(5,6,3,3))

Index=1:length(K_L)
h=(-log(2*p)- log(2*(1-p)) ) /2
plot(Index,(K_L), ylim=c(0,max((max(K_L)+max(K_L)), h*1.2)) ,
     xlab="Index",
     ylab=expression(paste("K-L divergence")),type="p", cex.axis=cex.axis, cex.lab=cex.lab,lwd=2)
abline(h=h,lty=2)
identify(1:length(K_L),K_L)

h=( (2*p-1)*log(2*p)+(2*(1-p)-1)*log(2*(1-p)) )*0.5
plot(Index,M_J,type="p", ylim=c(0,max((max(M_J)+.1*max(M_J)), h*1.2)) ,xlab="Index",
     ylab=expression(paste("J-distance")), cex.axis=cex.axis, cex.lab=cex.lab,lwd=2)
abline(h=h,lty=2)
identify(1:length(M_J),M_J)

h=(abs(2*p-1)+abs(2*(1-p)-1))/4
plot(Index,L1,type="p", ylim=c(0,max((max(L1)+.1*max(L1)),  h*1.2)) ,xlab="Index",
     ylab=expression(paste(L[1],"-distance")), cex.axis=cex.axis, cex.lab=cex.lab,lwd=2)
abline(h=h,lty=2)
identify(1:length(L1),L1)


#h=((2*p-1)^2+(2*(1-p)-1)^2)/2
h=(2*p*(1/(2*p)-1)^2+2*(1-p)*(1/(2*(1-p))-1)^2)/2
plot(Index,CH,type="p", ylim=c(0,max((max(CH)+.1*max(CH)),  h*1.2)) ,xlab="Index",
     ylab=expression(paste(chi^2,"-divergence")), cex.axis=cex.axis, cex.lab=cex.lab,lwd=2)
abline(h=h,lty=2)
identify(1:length(CH),CH)
Ind=cbind(y[c(47,211)],status[c(47,211)],X[c(47,211),-1],K_L[c(47,211)],
          M_J[c(47,211)],L1[c(47,211)],CH[c(47,211)])


xtable(Ind,digits = 3)
####

####### RR 

###############################################
# Probability
###############################################
Prob=function(y,X,vpar) {
  lam = vpar[1]        # dist exp por partes
  theta =vpar[2]        # theta (GP)
  beta = vpar[-c(1:2)]   # covariáveis
  vS = exp(-exp(lam)*y)
  linear=c(X%*%beta)
  mu=exp(linear)
  lambda=(1-theta)*mu
  ax = lambertW0(-theta*exp(-theta)*vS)
  Spop <- exp(-(lambda/theta)*(ax+theta))     
  p0 <- exp(-lambda)
  prob=p0/Spop
  prob=p0/Spop
  saida=c(p0,prob)
  saida
  
}
aa=Prob(2,XX[1,], dd[400,])
################################################
###############################################
# Probability
###############################################
P0=function(x,theta,b0,b1,b2) {
    mu0=exp(b0+b1*x)
    mu1=exp(b0+b2+b1*x)
  lambda0=(1-theta)*mu0
  lambda1=(1-theta)*mu1
  p0 <- exp(-lambda0)
  p1 <- exp(-lambda1)
  p=list(p0=p0,p1=p1)
 p
}

curve(P0(16,EB[2],EB[3],EB[4],EB[5])$p1,16,30)
curve(P0(x,EB[2],EB[3],EB[4],EB[5])$p0,20,60,ylim=c(0,1),add=F, 
      xlab="Age (years)", 
      ylab=" Probability of not dying",lwd=2, lty=1)
curve(P0(x,EB[2],EB[3],EB[4],EB[5])$p1,20,60,add=T,lwd=2,lty=2)
leg=c("no","yes")
legend("bottomleft",leg, lty=1:2 , title = "POI",bty="n",lwd=2)

Pcr=function(dados,y,X){
  R=length(dados[,1])
  Sal=matrix(0,R,2)
  for(k in 1:R){
    Sal[k,]=Prob(y,X,dados[k,]) 
  }
  colnames(Sal)=c("p0","Prob")
  Sal
}
AA=Pcr(dd,2,XX[1,])
AA1=Pcr(dd,2,XX[2,])
AA2=Pcr(dd,2,XX[3,])
AA3=Pcr(dd,2,XX[4,])
plot(density(AA[,1]),xlim=c(0,1), main="", 
     xlab="probability of not dying ",cex.lab=1.5,
     cex.axis=1.5, col=1,lwd=2)
text(locator(1),"A")
lines(density(AA1[,1]), prob=T,col=2,lwd=3)
text(locator(1),"B")
lines(density(AA2[,1]), prob=T,col=3,lwd=2)
text(locator(1),"C")
lines(density(AA3[,1]), prob=T,col=4,lwd=2)
text(locator(1),"D")

#################################################
# Patients
#################################################
p1=c(1,24,0)    # >65,no, CH I 
p2=c(1,24,1)    # =<65,no, CH I
p3=c(1,52,0)    # >65,yes, CH I
p4=c(1,52,1)    # =<65,yes, CH I
Resumo1 <- function(x) {
  s <- c(mean(x, na.rm = T), IHPD(x))
  names(s)<-c("m�dia", c("LI","LS"))
  return(s)
}
XX=rbind(p1,p2,p3,p4) 
nn=nrow(XX)
pp=NULL
for(i in 1:nrow(XX)){
  A=Pcr(dd,yA,XX[i,]) 
  Ss=c(apply(A,2,Resumo1))
  #names(Ss)=c("p0","Lp","Up","Pr","Lpr","Upr")
  pp=rbind(pp,Ss)
  colnames(pp)=c("p0","Lp","Up","Pr","Lpr","Upr")
}
S1=cbind(XX[,-1],pp)
print(pp,digits=3)
print(xtable(S1,digits = 3))
