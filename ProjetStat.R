#Préréglages

install.packages("tidyverse") #Installation de la librairie tidyverse, on utilise  notamment les packages dplyr et ggplot2
library(tidyverse) 
install.packages('VaRES') #La logCauchy est pas fournie avec le package de base de R, on a besoin de cette libraire pour l'avoir. L'expression de la log cauchy est fausse
library(VaRES)
setwd("D:/R/datasets") #définition du repository

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#Importation du fichier , extraction des valeurs importantes et calcul de grandeurs statistiques

Earthquakes=read.csv("Earthquakes.csv", header = TRUE, sep = ";", quote = "\"", dec = ".",
                     fill = TRUE, comment.char = "") #importation du fichier
Eff=read.csv("Eff.csv", header = TRUE, sep = ";", quote = "\"", dec = ".",
             fill = TRUE, comment.char = "") #Pour des histo
Eff2=read.csv("Eff2.csv", header = TRUE, sep = ";", quote = "\"", dec = ".",
             fill = TRUE, comment.char = "") #Pour des barplots
MyDataSet=Earthquakes%>%
  filter(YEAR>=1950,DAMAGE_MILLIONS_DOLLARS>0)%>%
  select(I_D,YEAR,DAMAGE_MILLIONS_DOLLARS)%>%
  arrange(DAMAGE_MILLIONS_DOLLARS)                   #Création du jeu exploité, on prend seulement les colonnes DAMAGE_MILLIONS_DOLLARS, year et I_D. Seulement les séismes après 1950 sont considérées
MyDataSetSummarized=MyDataSet%>%
  summarize(Moyenne=mean(DAMAGE_MILLIONS_DOLLARS),
            Mediane=median(DAMAGE_MILLIONS_DOLLARS),
            Quartile1=quantile(DAMAGE_MILLIONS_DOLLARS,0.25),
            Quartile3=quantile(DAMAGE_MILLIONS_DOLLARS,0.75),
            Variance=var(DAMAGE_MILLIONS_DOLLARS),
            EcartType=sd(DAMAGE_MILLIONS_DOLLARS),
            mu=log(median(MyDataSet$DAMAGE_MILLIONS_DOLLARS)),
            sigma=(log(quantile(MyDataSet$DAMAGE_MILLIONS_DOLLARS,0.75)/quantile(MyDataSet$DAMAGE_MILLIONS_DOLLARS,0.25)))/2) #Calcul de la moyenne, médiane, 1er et 3ème quartile,variance et écart  et stockage dans un tab

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#Calcul des paramètres de la Cauchy et log Cauchy, normale, exponentielle et Pareto

x0=median(MyDataSet$DAMAGE_MILLIONS_DOLLARS) #Paramètre de position de la Cauchy
gamma=(quantile(MyDataSet$DAMAGE_MILLIONS_DOLLARS,0.75)-quantile(MyDataSet$DAMAGE_MILLIONS_DOLLARS,0.25))/2 #Paramètre d'échelle de la Cauchy
logx0=log(median(MyDataSet$DAMAGE_MILLIONS_DOLLARS))#Paramètre de position de la logCauchy
loggamma=(log(quantile(MyDataSet$DAMAGE_MILLIONS_DOLLARS,0.75)/quantile(MyDataSet$DAMAGE_MILLIONS_DOLLARS,0.25)))/2 #Paramètre d'échelle de la logCauchy
moy=mean(MyDataSet$DAMAGE_MILLIONS_DOLLARS)
sd=sd(MyDataSet$DAMAGE_MILLIONS_DOLLARS)
  lambda=1/moy
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#Visualisation des différents tableaux
View(Earthquakes) #visualisation du tableau
View(MyDataSet) #visualisation du jeu de données
View(MyDataSetSummarized) #visualisation de la moyenne,médiane,quartiles, variance, écart type
View(effT)
View(Eff)
View(Eff2)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#Graphiques

x=seq(-100 ,1000,by =1)



HistLg=tibble(effExp,Damage_millions_dollars=seq(0,1000,by=50))


View(HistLg)

plot(x,dcauchy(x,location=x0,scale=gamma)) #graphe d'une cauchy de paramètre x0 et gamma calculé lignes 19&20
xlog=seq(1,1000,by=0.05)
plot(xlog,dlogcauchy(xlog,mu=logx0,sigma=loggamma)) #graphe d'une logcauchy de paramètre logx0 et loggamma calc ligne 21&22
xnorm=seq(-20000,40000, by =15)
plot(xnorm,dnorm(xnorm,moy,sd)) #graphe d'une loi normale de paramètre moyenne & écart type du jeu de données
hist(MyDataSet$DAMAGE_MILLIONS_DOLLARS,breaks=seq(from=0,to=1001,by=50))
Eff%>%
  ggplot(aes(x=DAMAGE.MillionUSD.Exp))+geom_histogram(binwidth = 50)
Eff%>%
  ggplot(aes(x=DAMAGE.MillionUSD.TH,color="red"))+geom_histogram(binwidth = 50)
MyDataSet%>%
  ggplot(aes(x=DAMAGE_MILLIONS_DOLLARS))+geom_density(kernel="gaussian") #marche pas
Eff2%>%
  ggplot(aes(DAMAGE.MillionUSD,color=ExpOuTh))+geom_bar()
Eff2%>%
  ggplot(aes(x=DAMAGE.MillionUSD,fill=ExpOuTh))+geom_histogram(binwidth = 50,stat="count",alpha=1,position=position_dodge())+
  theme_classic()+labs(title="Modélisation des coûts matériels lié à un séisme",x="Dégats en Millions de dollars américains",
  y="Compteur",fill="Type de données")+
  scale_fill_discrete(labels=c("Mesuré","Modelisé"))
Eff2%>%
  ggplot(aes(x=DAMAGE.MillionUSD,fill=ExpOuTh))+geom_bar(binwidth = 50,stat="count",alpha=1,position=position_dodge())+
  theme_classic()+labs(title="Modélisation des coûts matériels lié à un séisme",x="Dégats en Millions de dollars américains",
                       y="Compteur",fill="Type de données")+
  scale_fill_discrete(labels=c("Mesuré","Modelisé"))
MyDataSet%>%
  ggplot(aes(x=DAMAGE_MILLIONS_DOLLARS))+geom_area()
#hist(Eff2$DAMAGE.MillionUSD,breaks=seq(froù=0,to=1001,by=50))
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#tests du khi2 degré de liberté 20-1-1-1=17

effTCauchy=rep(0,21)
for(i in 1:20)
  {
  effTCauchy[i]=(pcauchy(i*50,location=x0,scale=gamma)-pcauchy((i-1)*50,location=x0,scale=gamma))*480
}
effTCauchy[21]=pcauchy(1000,location=x0,scale=gamma,lower.tail = TRUE)*480

effTLogCauchy=rep(0,21)
for(i in 1:20)
{
  print(i)
  effTLogCauchy[i]=(plogCauchy(i*50,logx0,loggamma)-plogCauchy((i-1)*50,logx0,loggamma))*480
}
effTLogCauchy[21]=(1-plogCauchy(1000,logx0,loggamma))*480

effTNormale=rep(0,21)
effTNormale[1]=pnorm(0,moy,sd)*480
for(i in 2:20)
{
  print(i)
  effTNormale[i]=(pnorm(i*50,moy,sd)-pnorm((i-1)*50,moy,sd))*480
}
effTNormale[21]=(1-pnorm(1000,moy,sd))*480

effTExp=rep(0,21)
for(i in 1:20)
{
  print(i)
  effTExp[i]=(pexp(i*50,lambda)-pexp((i-1)*50,lambda))*480
}
effTExp[21]=(1-pexp(1000,lambda))*480

print(effTCauchy)
print(effTLogCauchy)
print(effTNormale)
print(effTExp)
effExp=c(298,39,14,10,10,10,3,4,4,6,6,1,0,3,2,3,1,1,0,2,63) #Les valeurs 
print(effExp)
#chisq.test(effExp,p=effTCauchy)
monchi2=function(t,e)
{
  S=0
  r=0
  for(i in 1:21)
  {
    r=((t[i]-e[i])^2)/t[i]
    #print(r)
    S=S+r
  }
  return(S)
}
masomme=function(tab)
  {
  S=0
  for(i in 1:21)
    {
    S=S+tab[i]
  }
  return(S)
}
SEffExp=masomme(effExp)
SEffCauchy=masomme(effTCauchy)
SEffLogCauchy=masomme(effTLogCauchy)
SEffNormale=masomme(effTNormale)
SEffTExp=masomme(effTExp)
K2Cauchy=monchi2(effTCauchy,effExp)
K2LogCauchy=monchi2(effTLogCauchy,effExp)
K2Normale=monchi2(effTNormale,effExp)
K2Exp=monchi2(effTExp,effExp)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#Autres
pcauchy(x0,location=x0,scale=gamma) #calcul de P(X<x0)
dlogCauchy=function(x,u,s)
  {
  ret=0
  if(x>0)
    {
      ret=((1/(x*pi))*(s)/(((log(x)-u)^2)+s^2))
    }
    return(ret)
}
plogCauchy=function(x,u,s)
{
  ret=0
  if(x>0)
  {
    ret=((1/pi)*atan((log(x)-u)/s)+1/2)
  }
  return(ret)
}
plogCauchy(75,logx0,loggamma)
plot(xlog,dlogCauchy(xlog,logx0,loggamma)) #graphe d'une logcauchy de paramètre logx0 et loggamma calc ligne 21&22
