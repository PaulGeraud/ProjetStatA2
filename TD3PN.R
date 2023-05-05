rm(list=ls())#libere l'espace memoire des variables précedemment utilisées
setwd("D:/R/datasets")
library(ggforce)
library(tidyverse)
library(plot3D)
library(rayshader)
library(rgl)
library(plot3Drgl)
library(matlib)
library(mvtnorm)
#Ex1
g=function(x){return(exp(x))}
n=10000
X=runif(n)
Estimation=mean(g(X))
Estimation
I=exp(1)-1
#valeur cohérente
#var théorique et estimation
(-exp(2)+4*exp(1)-3)/2
var(g(X))

#Question 7: création de Sn
X=runif(n)
Y=1-X
EstimationAnti=mean((g(X)+g(Y))/2)
EstimationAnti
#la variance, donc la vitesse de convergence a été réduite
#-----------------------------------------------------------------------------------------
#Ex2
rm(list=ls())

