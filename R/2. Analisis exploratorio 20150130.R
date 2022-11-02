#GUION PARA ANALISIS DE DATOS DE BIORITMOS EN NEOTRYPAEA UNCINATA
#CREADO 20140930 FL+EN
#MODIFICACIONES ...

#----------------------------------------------
# Directorios y librerias
#----------------------------------------------
rm(list=ls())
setwd ("C:/Users/IMAR/DRIVE ULAGOS/Presentaciones a Congresos y Publicaciones/Publicaciones/Leiva et al JEB/R/datos")
getwd() #Verifica directorio
today<-Sys.Date()
library(nlme);library(car);library(multcomp);library(plotrix)

#----------------------------------------------
# Leo y exploro datos
#----------------------------------------------
datos1 <-read.table("labterreno2.csv",header=TRUE, sep=",", dec=".", strip.white=TRUE)
summary(datos1)
datos1$ANIMAL=as.factor(datos1$ANIMAL)
#datos1$EXPERIMENTO=as.factor(datos1$EXPERIMENTO)
str(datos1)

#1. Formulo modelo completo, considerando todas las variables como efectos fijos
lme1<-lme(MO2_DW~MAREA*OXIGENO*ACLIMATACION,data=datos1,random=~1|ANIMAL)
Anova(lme1)

#2. Evaluo Residuales...
res1<-resid(lme1,type="normalized")
#3. Normalidad
qqPlot(res1)
shapiro.test(res1)
#4. Homocedasticidad
boxplot(res1~MAREA+OXIGENO+ACLIMATACION,data=datos1); 
leveneTest(res1~MAREA*OXIGENO*ACLIMATACION,data=datos1)

lme2<-lme(MO2_DW~MAREA*OXIGENO+MAREA*ACLIMATACION+OXIGENO*ACLIMATACION,data=datos1,random=~1|ANIMAL)
Anova(lme2)
lme3<-lme(MO2_DW~MAREA*OXIGENO+OXIGENO*ACLIMATACION,data=datos1,random=~1|ANIMAL)
Anova(lme3)

datos1$trt=with(datos1,factor(paste(MAREA,OXIGENO,ACLIMATACION,sep=".")))

lme4<-lme(MO2_DW~trt,data=datos1,random=~1|ANIMAL)
glht1<-glht(lme4,linfct=mcp(trt="Tukey"))#test a posteriori
summary(glht1)
cld(glht1)#asiga letras de diferencias significativas
pred.mat=expand.grid(MAREA=unique(datos1$MAREA),OXIGENO=unique(datos1$OXIGENO),
                     ACLIMATACION=unique(datos1$ACLIMATACION))
pred.mat$pred=predict(lme3,newdata=pred.mat,level=0)



#Figuras para panel CCM 2016

jpeg(paste("2.1. Figure 1 (",today,").jpeg"),
     height = 15,width =35,units="cm",quality=100,res=300) 
par(mfrow=c(1,2),oma=c(1.5,0,0,0),mar=c(2.5,5.5,1,1));#mar fija margenes alrededor de figuras

b1<-boxplot(MO2_DW~MAREA+OXIGENO,data=subset(datos1,ACLIMATACION=="TERRENO"),
            boxwex=0.6,#ancho de la caja del boxplot
            ylab=bquote("OCR"~paste("(",mu,"mol ",.(O[2]~h^-1),g^-1,")",sep="")),
            cex.lab=1.7,ylim=c(0,10.1),las=1,
            xaxt="n",
            cex.lab=1,cex.axis=1.4,
            col=c("steelblue1","steelblue1","darkgoldenrod3","darkgoldenrod3"))#colores de las cajas
axis(1,at=c(1:4),
     labels=c("Alta","Baja","Alta","Baja"),cex.axis=1.7)
#mtext("Hipoxia                Normoxia",cex=1.7,
      side=1,outer=TRUE,line=0.35,adj=0.229)
text(0.6,c(10),"A)",cex=2) 
dif<-text(c(1:4),(b1$stats[5,]+0.4),c("ab","ab","b","c"),cex=1.4)


b2<-boxplot(MO2_DW~MAREA+OXIGENO,data=subset(datos1,ACLIMATACION=="LABORATORIO"),
            boxwex=0.6,#ancho de la caja del boxplot
            ylab=bquote("OCR"~paste("(",mu,"mol ",.(O[2]~h^-1),g^-1,")",sep="")),
            cex.lab=1.7,ylim=c(0,10.1),las=1,
            xaxt="n",
            cex.lab=1,cex.axis=1.4,
            col=c("tomato1","tomato1","olivedrab3","olivedrab3"))#colores de las cajas
axis(1,at=c(1:4),
     labels=c("Alta","Baja","Alta","Baja"),cex.axis=1.7)
#mtext("Hipoxia                Normoxia",cex=1.7,
side=1,outer=TRUE,line=0.35,adj=0.229)
text(0.6,c(10),"B)",cex=2) 
dif<-text(c(1:4),(b2$stats[5,]+0.4),c("a","a","ab","b"),cex=1.4) 

# Dibujo rectángulos sobre boxplot rec(xleft,ybottom,xright,ytop)
rect((1-0.3),b2$stats[2,1],(1+0.3),b2$stats[4,1],angle=45,density=15,lwd=0.5) 
rect((2-0.3),b2$stats[2,2],(2+0.3),b2$stats[4,2],angle=45,density=15,lwd=0.5) 
rect((3-0.3),b2$stats[2,3],(3+0.3),b2$stats[4,3],angle=45,density=15,lwd=0.5) 
rect((4-0.3),b2$stats[2,4],(4+0.3),b2$stats[4,4],angle=45,density=15,lwd=0.5) 

dev.off()