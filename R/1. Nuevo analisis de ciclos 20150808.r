#GUION PARA ANALISIS DE DATOS DE ...
#CREADO 20130813 FL+EN
#MODIFICACIONES ...

#----------------------------------------------
# Directorios y librerias

rm(list=ls())
setwd ("~/felix.leiva@ulagos.cl/Laboratorio.Poblaciones.Marinas/Publicaciones/Leiva et al JEB/R")
getwd() #Verifica directorio
#library(Rcmdr)
today<-Sys.Date()
library(nlme);library(car);library(Hmisc);library(mgcv)
library(synchrony);
#----------------------------------------------
# Leo y exploro datos
#----------------------------------------------
datos1 <-read.table("datos/tidalneun.csv",header=TRUE, sep=",", dec=".", strip.white=TRUE)
names(datos1)
datos1$animal<-as.factor(datos1$animal)
datos1<-datos1[order(datos1$exper1,datos1$time2),]
datos1$lectura<-rep(1:129,each=9)
par(mfrow=c(1,1))
#Mostrar variabilidad en un ejemplo
with(subset(datos1, exper1=="E2" & animal=="7"),plot(mo2gs~altura))
# Mostrar variabilidad por experimento
with(datos1,boxplot(mo2gs~as.factor(dias2),at=unique(datos1$dias2)))
#Elimina datos extremos
datos2<-na.omit(subset(datos1,mo2gs<8 & !(animal==2 & time==4.859375)))
datos2<-subset(datos2,time2<2 & !(time2<0.45 & exper1=="E5"))
datos2<-subset(datos2,!(ID==258))
#----------------------------------------------
# estandariza datos de metabolismo y marea
#----------------------------------------------
datos2$mo2gs.c<-datos2$altura.c<-NA
for (e in unique(datos2$exper1))
{
datos2$altura.c[datos2$exper1==e]<-scale(datos2$altura[datos2$exper1==e])
datos2$mo2gs.c[datos2$exper1==e]<-scale(datos2$mo2gs[datos2$exper1==e])
}
#----------------------------------------------
# H0: modelo nulo
#----------------------------------------------
lme.0<-lme(mo2gs.c~1,random=~1|animal,cor=corAR1(form=~time2|animal/exper1),data=datos2)
AIC(lme.0);logLik(lme.0)  
#----------------------------------------------
# H1: ciclo mareal continuo
#----------------------------------------------
lme.1<-lme(mo2gs.c~altura.c*exper1,random=~1|animal,cor=corAR1(form=~time2|animal/exper1),data=datos2)
Anova(lme.1,type="II");AIC(lme.1);logLik(lme.1)  
#----------------------------------------------
# H2: ciclo mareal discreto
#----------------------------------------------
lme.2<-lme(mo2gs.c~(altura.c>0)*exper1,random=~1|animal,cor=corAR1(form=~time2|animal/exper1),data=datos2)
Anova(lme.2,type="II");AIC(lme.2);logLik(lme.2)  
#----------------------------------------------
# H3: ciclo mareal continuo con desfase
#----------------------------------------------
#Predictor de altura mareal por experimento
gam.a.lst<-lapply(unique(datos2$exper1),function(x)
 {
   gam(altura.c~s(time2),data=subset(datos2,exper1==x))
 })
 names(gam.a.lst)<-unique(datos2$exper1)
#  marea1<-function(ee,tt)
#  {
#    tt.df<-data.frame(time2=tt,pred=NA)
#   sapply(1:dim(tt.df)[1],function(x) predict(gam.a.lst[[ee[x]]],newdata=tt.df[x,]))
# }
#funcion predictora
marea2<-function(ee,tt)
{
  de<-ifelse(ee=="E1",0.0689,
             ifelse(ee=="E2",-0.2316,
                    ifelse(ee=="E3",0.1237,
                           ifelse(ee=="E4",-0.1164,0.1555))))
  1.5*sin((1.93*(tt+de)*360-90)*pi/180)
}
#modelo no lineal (un set de parametros por experimento)
nls1<-nlme(mo2gs.c~(b1*marea2(ee=exper1,tt=(time2+b2))),
           fixed = b1+b2 ~ exper1,
           random =b1 ~ 1|animal,
           start=c(b1=rep(0.5,5),b2=rep(0,5)),data=subset(datos2,mo2gs.c<3.8))
anova(nls1,type="marginal")
AIC(nls1);logLik(nls1)  
summary(nls1)$tTable
(desfase=summary(nls1)$tTable[2,1]*24)
(desfase.ee=(summary(nls1)$tTable[2,2]^2*24^2)^0.5)
#modelo no lineal (un set de parametros comun)
nls2<-nlme(mo2gs.c~(b1*marea2(ee=exper1,tt=(time2+b2))),
           fixed = b1+b2 ~ 1,
           random =b1 ~ 1|animal,
           start=c(b1=rep(0.5,1),b2=rep(0,1)),data=subset(datos2,mo2gs.c<3.8))
anova(nls2,type="marginal")
AIC(nls2);logLik(nls2)  
summary(nls2)$tTable


#----------------------------------------------
# H4: ciclo diurno
#----------------------------------------------
#Maximo a mediodia
datos2$d1<-sin((datos2$time2*360-90)*pi/180)
lme.2<-lme(mo2gs.c~d1*exper1,random=~1|animal,cor=corAR1(form=~time2|animal/exper1),data=datos2)
Anova(lme.2,type="II");AIC(lme.2);logLik(lme.2)  
#----------------------------------------------
# H5: ciclo semi-diurno
#----------------------------------------------
datos2$d2<-sin((2*datos2$time2*360-90)*pi/180)
lme.3<-lme(mo2gs.c~d2*exper1,random=~1|animal,cor=corAR1(form=~time2|animal/exper1),data=datos2)
Anova(lme.3,type="II");AIC(lme.3);logLik(lme.3)  
#----------------------------------------------
# H6: sinusoidal endógeno
#----------------------------------------------
#Buscando otro ciclo (b1=media,b2=numero de ciclos por dia)
ciclo<-function(b1,b2,time2) (b1)*(sin((b2*time2*360-90)*pi/180))
# un ciclo por experimento
nls3<-nlme(mo2gs.c~ciclo(b1,b2,time2),
           fixed = b1+b2 ~ exper1,
           random =b1 ~ 1|animal,
           start=c(b1=rep(2,5),b2=rep(1.93,5)),data=subset(datos2,mo2gs.c<3.8))
anova(nls3);
logLik(nls3);AIC(nls3)
summary(nls3)$tTable
# un ciclo comun para todos los experimentos
nls4<-nlme(mo2gs.c~ciclo(b1,b2,time2),
           fixed = b1+b2 ~ 1,
           random =b1 ~ 1|animal,
           start=c(b1=rep(4),b2=rep(2)),data=subset(datos2,mo2gs.c<3.8))
# compara ambos modelos
anova(nls1,nls2,nls3,nls4)
# observa parametros
summary(nls1)$tTable
# observa residuales
par(mfrow=c(1,1))
plot(nls3);qqPlot(resid(nls3,type="pearson"))

# Estima periodo de los ciclos por experimento
ss1<-summary(nls3)$tTable
resumen1<-data.frame()
for (x in 1:1000){
s1<-apply(ss1,1,function(y) rnorm(1,y[1],y[2]))
resumen0<-data.frame(days=c(1,8,15,22,29),
                    period=round(((asin(1)*180/pi)+90)/(360*c(s1[6],(s1[6]+s1[7:10])))*2*24,1),
                    amplitude=round(abs(c(s1[1],(s1[1]+s1[2:5]))),2))
resumen1=rbind(resumen1,resumen0)
}
resumen.mean<-aggregate(cbind(period,amplitude)~days,data=resumen1,mean)
resumen.ee<-aggregate(cbind(period,amplitude)~days,data=resumen1,sd)
resumen<-merge(resumen.mean,resumen.ee,by="days")
resumen

resumen$letra<-c("A)","B)","C)","D)","E)")
# Grafica
jpeg(paste("1.1. mo2gs crudos + gam + marea (",today,").jpg",sep=""),width=12,height=20,units="cm",res=300,quality=100)
par(mfrow=c(5,1),mar=c(0,2,1,1),oma=c(4.2,4,0,4));
for (i in 1:5)
{
  subset1<-subset(datos2,dias2==(resumen$days[i]-1))
  #estandariza metabolismo
  subset1$mo2gs.c<-scale(subset1$mo2gs)
  y.all<-mean(subset1$mo2gs);s.all<-sd(subset1$mo2gs)
  
  #Crea matriz predictora 
  plot(mo2gs.c~time2,subset1,cex=0.5,ylim=c(-3,4),xlim=c(0.4,1.52),xaxt="n",las=1)
  xmin<-min(subset1$time2);xmax<-max(subset1$time2)
  lines(seq(xmin,xmax,0.01),
        predict(nls3,newdata=data.frame(exper1=unique(subset1$exper1),time2=seq(xmin,xmax,0.01)),level=0),
        col="red")
  #Grafica
  axis(4,at=(seq(0,7,2)-y.all)/s.all,seq(0,7,2),las=1)
  if(i==5) {axis(1,at=seq(0.2,1.6,0.083),round(seq(0.2,1.6,0.083)*24,0),las=1)}else
    {axis(1,at=seq(0.2,1.6,0.083),labels=FALSE)}
  #text(0.38,4,paste("Days ex-situ=",resumen$days[i],sep=""),adj=c(0,1),cex=0.8)
  text(0.38,4,paste(resumen$letra[i],sep=""),adj=c(0,1),cex=0.8,font=2)
  text(0.38,3,paste("T=",round(resumen$period.x[i],1)," h",sep=""),adj=c(0,1),cex=0.8)
  text(0.38,2.4,paste("RA=",round(resumen$amplitude.x[i],1),sep=""),adj=c(0,1),cex=0.8)
  #Marea (alturas suavizadas por gam)
  altura<-predict(gam.a.lst[[i]],newdata=data.frame(time2=seq(xmin,xmax,0.01)))
  lines(seq(xmin,xmax,0.01),altura,cex=0.5,lty=2,lwd=2,col="blue")
}
#mtext(side=1,line=3,outer=T,"Time (hours after first midnight)")
mtext(side=1,line=3,outer=T,"Tiempo (horas)")
#mtext(side=2,line=1,outer=T,"Standardized OCR")
mtext(side=2,line=1,outer=T,"OCR estandarizado")
mtext(side=4,line=2,outer=T,bquote("OCR"~paste("(",mu,"mol ",.(O[2]~h^-1),g^-1,")",sep="")))
dev.off()

#-------------------------------------------------------------------------------
# Evalúa hipotesis de sincronia entre individuos
# F-value es el cuociente entre la varianza entre mediciones v/s la varianza dentro
# de las mediciones. El test de F es válido para evaluar si este cuociente es 
# mayor que uno.
#-------------------------------------------------------------------------------
datos3<-datos2
datos3$time3=0.01*round(datos3$time2/0.01,0)
lme.1<-lme(mo2gs.c~as.factor(time3),random=~1|animal,cor=corAR1(form=~time3|animal),data=subset(datos3,exper1=="E1"))
anova(lme.1,type="sequential")      
lme.2<-lme(mo2gs.c~as.factor(time3),random=~1|animal,cor=corAR1(form=~time3|animal),data=subset(datos3,exper1=="E2"))
anova(lme.2,type="sequential")      
lme.3<-lme(mo2gs.c~as.factor(time3),random=~1|animal,cor=corAR1(form=~time3|animal),data=subset(datos3,exper1=="E3"))
anova(lme.3,type="sequential")      
lme.4<-lme(mo2gs.c~as.factor(time3),random=~1|animal,cor=corAR1(form=~time3|animal),data=subset(datos3,exper1=="E4"))
anova(lme.4,type="sequential")      
lme.5<-lme(mo2gs.c~as.factor(time3),random=~1|animal,cor=corAR1(form=~time3|animal),data=subset(datos3,exper1=="E5"))
anova(lme.5,type="sequential")      
#-------------------------------------------------------------------------------
# Utilizando ajustes no lineales (gam) de los promedios por lectura en función del
# tiempo evaluamos si los patrones de respuesta son similares o diferentes entre
# experimentos: formulamos y comparamos diversas hipotesis mediante AIC
#-------------------------------------------------------------------------------
# Hipotesis 1: no existen diferencias en patrones entre experimentos
gam.h1<-gamm(mo2gs.c~s(time2,k=12),data=subset(datos2),random=list(animal=~1,exper1=~1),method="ML")
aic.h1<-AIC(gam.h1$lme)
# Hipotesis 2: Todos los patrones son distintos entre experimentos
# gam.h2<-gam(lm.pred~s(time2,k=13,by=exper),data=subset(datos4),cor=corAR1(form=~time2))
# aic.h2<-AIC(gam.h2)
gam.h21<-gamm(mo2gs.c~s(time2,k=12),data=subset(datos3,exper1=="E1"),random=list(animal=~1,exper1=~1),method="ML")
gam.h22<-gamm(mo2gs.c~s(time2,k=12),data=subset(datos3,exper1=="E2"),random=list(animal=~1,exper1=~1),method="ML")
gam.h23<-gamm(mo2gs.c~s(time2,k=12),data=subset(datos3,exper1=="E3"),random=list(animal=~1,exper1=~1),method="ML")
gam.h24<-gamm(mo2gs.c~s(time2,k=12),data=subset(datos3,exper1=="E4"),random=list(animal=~1,exper1=~1),method="ML")
gam.h25<-gamm(mo2gs.c~s(time2,k=12),data=subset(datos3,exper1=="E5"),random=list(animal=~1,exper1=~1),method="ML")
aic.h2<-sum(AIC(gam.h21$lme),AIC(gam.h22$lme),AIC(gam.h23$lme),AIC(gam.h24$lme),AIC(gam.h25$lme))
# Hipotesis 3: Los experimentos 1 y 2 poseen el mismo patron
gam.h3<-gamm(mo2gs.c~s(time2,k=12),data=subset(datos3,exper1 %in% c("E1","E2")),random=list(animal=~1,exper1=~1),method="ML")
aic.h3<-sum(AIC(gam.h3$lme),AIC(gam.h23$lme),AIC(gam.h24$lme),AIC(gam.h25$lme))
# Hipotesis 4: Los experimentos 1, 2 y 3 poseen el mismo patron
gam.h4<-gamm(mo2gs.c~s(time2,k=12),data=subset(datos3,exper1 %in% c("E1","E2","E3")),random=list(animal=~1,exper1=~1),method="ML")
aic.h4<-sum(AIC(gam.h4$lme),AIC(gam.h24$lme),AIC(gam.h25$lme))
# Hipotesis 5: Los experimentos 1, 2, 3 y 4 poseen el mismo patron
gam.h5<-gamm(mo2gs.c~s(time2,k=12),data=subset(datos3,exper1 %in% c("E1","E2","E3","E4")),random=list(animal=~1,exper1=~1),method="ML")
aic.h5<-sum(AIC(gam.h5$lme),AIC(gam.h25$lme))
#compara
aic.h1;aic.h2;aic.h3;aic.h4;aic.h5
#-----------------------------------------------------------------------
# Conclusión: la hipotesis más informativa fue aquella que consideró  un patrón
# común para los experimentos 1 y 2 y patrones distintos para los experimentos 3-5

# #-----------------------------------------------------------------------
# #Analizo utilizando synchronize
# #------------------------------------------------------------------------------------
# E="E3"
# t1=cbind(datos3.1[,c("time2","altura.c")])
# t2=cbind(datos3.1[,c("time2","altura.c")])
# ps1<-phase.sync(t1[,2],t2[,2],nrands=100)
# ps1$pval
# ## Plot distribution of phase difference
# hist(ps1$deltaphase$phasediff)
# ps1$Q.obs;ps1$pval
# 
# ps2<-peaks(t1,t2,nrands=100)
# ps2$pval
