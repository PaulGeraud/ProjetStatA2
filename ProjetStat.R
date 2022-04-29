#Préréglages

 
library(tidyverse) #Installation de la librairie tidyverse, on utilise  notamment les packages dplyr et ggplot2
library(VaRES) #La logCauchy & Pareto n'est pas fournie avec le package de base de R, on a besoin de cette libraire pour l'avoir. L'expression de la log cauchy est fausse
setwd("D:/R/datasets") #définition du repository, à modifier selon l'installation

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#Importation du fichier , extraction des valeurs importantes et calcul de grandeurs statistiques

Earthquakes=read.csv("Earthquakes.csv", header = TRUE, sep = ";", quote = "\"", dec = ".",
                     fill = TRUE, comment.char = "") #importation du fichier

MyDataSet=Earthquakes%>%
  filter(YEAR>=1950,DAMAGE_MILLIONS_DOLLARS>0)%>%
  select(YEAR,DAMAGE_MILLIONS_DOLLARS)%>%
  arrange(DAMAGE_MILLIONS_DOLLARS)%>%
  mutate(logDMG=log(DAMAGE_MILLIONS_DOLLARS))#Création du jeu exploité, on prend seulement les colonnes DAMAGE_MILLIONS_DOLLARS, year et I_D. Seulement les séismes après 1950 sont considérées
MyDataSetSummarized=MyDataSet%>%
  summarize(Moyenne=mean(DAMAGE_MILLIONS_DOLLARS),
            Mediane=median(DAMAGE_MILLIONS_DOLLARS),
            Quartile1=quantile(DAMAGE_MILLIONS_DOLLARS,0.25),
            Quartile3=quantile(DAMAGE_MILLIONS_DOLLARS,0.75),
            Variance=var(DAMAGE_MILLIONS_DOLLARS),
            EcartType=sd(DAMAGE_MILLIONS_DOLLARS),
            x0=log(median(MyDataSet$DAMAGE_MILLIONS_DOLLARS)),
            gamma=(log(quantile(MyDataSet$DAMAGE_MILLIONS_DOLLARS,0.75)/quantile(MyDataSet$DAMAGE_MILLIONS_DOLLARS,0.25)))/2) #Calcul de la moyenne, médiane, 1er et 3ème quartile,variance et écart  et stockage dans un tab
Eff2=read.csv("Eff2.csv", header = TRUE, sep = ";", quote = "\"", dec = ".",
              fill = TRUE, comment.char = "") #Pour des barplots
plage=read.csv("plage.csv",sep=";") #Juste des valeurs allant de -10 à 10, c'est utilisé pour faire le graphique qui compare loi normale et loi de Cauchy
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#Calcul des paramètres de la Cauchy et log Cauchy, normale, exponentielle et Pareto

x0=median(MyDataSet$DAMAGE_MILLIONS_DOLLARS) #Paramètre de position de la Cauchy
gamma=(quantile(MyDataSet$DAMAGE_MILLIONS_DOLLARS,0.75)-quantile(MyDataSet$DAMAGE_MILLIONS_DOLLARS,0.25))/2 #Paramètre d'échelle de la Cauchy
logx0=log(median(MyDataSet$DAMAGE_MILLIONS_DOLLARS))#Paramètre de position de la logCauchy
loggamma=(log(quantile(MyDataSet$DAMAGE_MILLIONS_DOLLARS,0.75)/quantile(MyDataSet$DAMAGE_MILLIONS_DOLLARS,0.25)))/2 #Paramètre d'échelle de la logCauchy
moy=mean(MyDataSet$DAMAGE_MILLIONS_DOLLARS) #moyenne empirique
sd=sd(MyDataSet$DAMAGE_MILLIONS_DOLLARS) #Ecart type non biaisé
lambda=1/moy #Paramètre d'une loi exponentielle
xm=min(MyDataSet$DAMAGE_MILLIONS_DOLLARS) #Scale parameter Pareto
theta=480/(sum(MyDataSet$logDMG)-480*xm) #Shape parameter Pareto
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#Visualisation des différents tableaux
View(Earthquakes) #visualisation du tableau
View(MyDataSet) #visualisation du jeu de données
View(MyDataSetSummarized) #visualisation de la moyenne,médiane,quartiles, variance, écart type
View(Eff2)
View(plage)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#Graphiques

Eff2%>%
  ggplot(aes(x=DAMAGE.MillionUSD,fill=Model))+geom_bar(stat="count",alpha=1,position=position_dodge())+
  theme_minimal()+labs(title="Modélisation des coûts matériels lié à un séisme",x="Dégats en Millions de dollars américains",y="Compteur",fill="Type de données")+
  scale_x_continuous(breaks=seq(0,1000, by =50))+
  scale_y_continuous(breaks=seq(0,450, by =50))+
  scale_fill_discrete(labels=c("Exponentielle","LogCauchy","Mesuré","Pareto"))+
  annotate(geom="text",x =1040, y=-10,label="Coût >=1000",color="black",size=3)
  #Diagramme bâton
Eff2%>%
  ggplot(aes(x=DAMAGE.MillionUSD))+
  stat_function(fun=dlogcauchy,n=100,args=list(mu=logx0,sigma=loggamma),colour="red")+annotate(geom="text",x =1000, y=0.0055,label="Log-Cauchy",color="red")+
  stat_function(fun=dexp,n=100,args=list(rate=lambda),colour="blue")+annotate(geom="text",x =1000, y=0.005,label="Exponentielle",color="blue")+
  stat_function(fun=dpareto,n=100,args=list(K=xm,c=theta),colour="green4")+annotate(geom="text",x =1000, y=0.0045,label="Pareto",color="green4")+
  theme_classic()+
  labs(title="Fonction de densité de différentes distributions",x="x",y="y") #Représentation graphiques de densités,
Eff2%>%
  ggplot(aes(x=DAMAGE.MillionUSD))+
  stat_function(fun=plogcauchy,n=100,args=list(mu=logx0,sigma=loggamma),colour="red")+annotate(geom="text",x =1000, y=0.95,label="Log-Cauchy",color="red")+
  stat_function(fun=pexp,n=100,args=list(rate=lambda),colour="blue")+annotate(geom="text",x =1000, y=0.90,label="Exponentielle",color="blue")+
  stat_function(fun=ppareto,n=100,args=list(K=xm,c=theta),colour="green4")+annotate(geom="text",x =1000, y=0.85,label="Pareto",color="green4")+
  theme_classic()+
  labs(title="Fonction de densité de différentes distributions",x="x",y="y") #Représentation graphiques des fonctions de répartitions,

plage%>%
  ggplot(aes(x=Plage))+stat_function(fun=dnorm,n=100)+stat_function(fun=dcauchy,n=100,colour="red")+
  scale_y_continuous(breaks = c(0,0.32,0.4))+scale_x_continuous(breaks=seq(-10,10, by =10))+theme_classic()+
  annotate(geom="text",x =8.3, y=0.5,label="N(0,1)",color="black")+
  annotate(geom="text",x =9, y=0.45,label="Cauchy(0,1)",color="red")+
  labs(title="Comparaison graphique entre la loi de Cauchy(0,1) et la loi Normale(0,1)",x="x",y="y")
 #représentation graphique d'une densité d'une loi N(0,1) et d'une Cauchy(0,1)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#Calculs des effectifs théoriques qui seraient attendus si cela suivait une loi donnée (pour préparer le khi 2)
#tests du khi2 degré de liberté 21-1-1-1=18 (car 2 paramètres estimés) 19 pour l'exponentielle car 1 seul paramètre
#Les effectifs sont organisés par tranche de 50 MUSD, sauf la dernière tranche qui est >1000 MUSD

#Suivant une Cauchy, ne marche pas
effTCauchy=rep(0,21)
effTCauchy[1]=pcauchy(0,location=x0,scale=gamma)*480
for(i in 1:20)
  {
  effTCauchy[i]=(pcauchy(i*50,location=x0,scale=gamma)-pcauchy((i-1)*50,location=x0,scale=gamma))*480
}
effTCauchy[21]=pcauchy(1000,location=x0,scale=gamma,lower.tail = TRUE)*480

#Suivant une log Cauchy, je n'utilise pas la fonction du package VaRES car pour une raison que j'ignore, elle retourne des valeurs erronées pour le 21ème groupe
plogCauchy=function(x,u,s)
{
  ret=0
  if(x>0)
  {
    ret=((1/pi)*atan((log(x)-u)/s)+1/2)
  }
  return(ret)
}
effTLogCauchy=rep(0,21)
for(i in 1:20)
{
  
  effTLogCauchy[i]=(plogCauchy(i*50,logx0,loggamma)-plogCauchy((i-1)*50,logx0,loggamma))*480
}
effTLogCauchy[21]=(1-plogCauchy(1000,logx0,loggamma))*480

#Suivant une Pareto
effTPareto=rep(0,21)
effTPareto[1]=ppareto(50,K=xm,c=theta)*480
for(i in 2:20)
{

  effTPareto[i]=(ppareto(i*50,K=xm,c=theta)-ppareto((i-1)*50,K=xm,c=theta))*480
}
effTPareto[21]=(1-ppareto(1000,K=xm,c=theta))*480

#Suivant une normale
effTNormale=rep(0,21)
effTNormale[1]=pnorm(0,moy,sd)*480
for(i in 2:20)
{
  print(i)
  effTNormale[i]=(pnorm(i*50,moy,sd)-pnorm((i-1)*50,moy,sd))*480
}
effTNormale[21]=(1-pnorm(1000,moy,sd))*480

#Suivant une exponentielle
effTExp=rep(0,21)
for(i in 1:20)
{
  print(i)
  effTExp[i]=(pexp(i*50,lambda)-pexp((i-1)*50,lambda))*480
}
effTExp[21]=(1-pexp(1000,lambda))*480

#Affichage de ces effectifs théoriques
print(effTCauchy)
print(effTLogCauchy)
print(effTNormale)
print(effTExp)
print(effTPareto)

#Les effectifs observées dans le tableau, comptés à la main car je ne savais pas comment faire autrement
effObs=c(298,39,14,10,10,10,3,4,4,6,6,1,0,3,2,3,1,1,0,2,63) 
print(effObs)
#chisq.test(effExp,p=effTCauchy) Je n'utilise pas cette méthode

#Fonction fait main qui retourne le coeff du khi 2, t correspond aux effectifs théoriques, e aux effectifs mesurées
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

#Somme tous les effectifs, pour être sûr que il y a bien 480 valeurs
SEffExp=masomme(effExp)
SEffCauchy=masomme(effTCauchy)
SEffLogCauchy=masomme(effTLogCauchy)
SEffNormale=masomme(effTNormale)
SEffTExp=masomme(effTExp)
SEffTPareto=masomme(effTPareto)

#Calcul des coefficients du khi2
K2Cauchy=monchi2(effTCauchy,effObs) #Ne pas prendre en compte
K2LogCauchy=monchi2(effTLogCauchy,effObs)
K2Normale=monchi2(effTNormale,effObs)
K2Exp=monchi2(effTExp,effObs)
K2Pareto=monchi2(effTPareto,effObs)

