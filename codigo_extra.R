############## Clase práctica de GWAS
############## MATERIAL SUPLEMENTARIO
############## Dr. Martín Nahuel García
##############

#########Previo a correr el código: generar el objeto nsimout
####
####simular familias de hermanos completos
####

library(simulMGF)
#library(factoMineR)
library(factoextra)

simulN(20, 1000, 50, 0.5, 30, 0.5)
geno.madres = nsimout$geno[1:10,]
geno.padres = nsimout$geno[11:20,]

#generar matriz de genotipos de 10 familias de 50 hermanos completos
simulFS(geno.madres, geno.padres, 50)
#[1] "simulatedFS was generated"

scaled_data <- scale(simulatedFS, center = TRUE, scale = TRUE)
pca_result <- prcomp(scaled_data, graph = FALSE)

#rotulos de familias
familias = c(rep("fam1", 50),rep("fam2", 50),rep("fam3", 50),rep("fam4", 50),rep("fam5", 50),
		rep("fam6", 50),rep("fam7", 50),rep("fam8", 50),rep("fam9", 50),rep("fam10", 50))
str(familias)
# chr [1:500] "fam1" "fam1" "fam1" "fam1" "fam1" "fam1" "fam1" "fam1" "fam1" ...

p <- fviz_pca_ind(pca_result, label="none", habillage=familias,
             addEllipses=TRUE, ellipse.level=0.95)
print(p)

##########
##########
##########

# Para generar un fenotipo usamos la función simPheno:
# simPheno(x, Nqtl, Esigma, Pmean, Perror)
# donde x = matriz de genotipos
# Nqtl = número de QTLs (marcadores asociados)
# Esigma = desvío estándar de la distribución Normal de donde se muestrean los efectos
# Pmean = media del fenotipo
# Perror = desvío estándar del fenotipo

simPheno(simulatedFS, 50, 0.5, 30, 0.5)
#[1] "simP was generated"
str(simP)
#List of 3
# $ pheno   : num [1:500, 1] 31.8 28.5 29.3 ...
# $ QTN     : int [1:50] 516 250 349 821 ...
# $ Meffects: num [1:50] -0.2742 1.2245 -0.2029 0.2292 ...

###########
###########
###########

fam10FS = cbind(familias, simP$pheno[1:500])
fam10FS = as.data.frame(fam10FS)
colnames(fam10FS) = c("familias", "feno")
fam10FS$feno = as.numeric(fam10FS$feno)
head(fam10FS)
#  familias     feno
#1     fam1 31.78020
#2     fam1 28.53431
#3     fam1 30.42450

hist(fam10FS$feno)

library(dplyr)
library(ggplot2)

medias_por_grupo <- fam10FS %>%
  group_by(familias) %>%
  summarise(media_grupo = mean(feno))

ggplot(medias_por_grupo, aes(x = factor(familias), y = media_grupo)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Media por Grupo",
       x = "Familia",
       y = "Media")

boxplot(feno ~ familias, data = fam10FS)

 



#######################
########################
setwd("E:/MGYGV/2025_clase_practica")
ped = read.csv("pedigree.csv", header = T)
head(ped)
tail(ped)

library(AGHmatrix)

A = Amatrix(ped, ploidy=2)
A[1:30,1:30]

library(fields)
image.plot(A, 1:520, 1:520)

### matriz G: matriz de relaciones realizadas

Gmatrix <- Gmatrix(simulatedFS, method="VanRaden", ploidy=2, ratio=FALSE)  
colnames(Gmatrix) = paste0("I", c(1:500))
rownames(Gmatrix) = paste0("I", c(1:500))
Gmatrix[1:10, 1:10]
#           I1        I2        I3        I4        I5        I6        I7        I8        I9       I10
#I1  1.0206027 0.6396616 0.6499071 0.6721293 0.6742828 0.6721459 0.6469005 0.6685512 0.6923513 0.6776041
#I2  0.6396616 1.0124361 0.6458238 0.6576928 0.6412104 0.6659919 0.6635237 0.6686092 0.6385724 0.6652382
#I3  0.6499071 0.6458238 1.0122208 0.6803622 0.6452440 0.6700255 0.7006877 0.6602189 0.6591712 0.6423534
#I4  0.6721293 0.6576928 0.6803622 1.0111109 0.6653956 0.6736119 0.6566491 0.6638053 0.6586163 0.6832114
#I5  0.6742828 0.6412104 0.6452440 0.6653956 0.9802168 0.6633415 0.6650146 0.6742414 0.6276394 0.6356693
#I6  0.6721459 0.6659919 0.6700255 0.6736119 0.6633415 1.0504863 0.6960080 0.6638219 0.6834806 0.6997932
#I7  0.6469005 0.6635237 0.7006877 0.6566491 0.6650146 0.6960080 1.0538325 0.6924134 0.6789417 0.6869717
#I8  0.6685512 0.6686092 0.6602189 0.6638053 0.6742414 0.6638219 0.6924134 1.0205199 0.6343317 0.6713507
#I9  0.6923513 0.6385724 0.6591712 0.6586163 0.6276394 0.6834806 0.6789417 0.6343317 1.0225657 0.7013627
#I10 0.6776041 0.6652382 0.6423534 0.6832114 0.6356693 0.6997932 0.6869717 0.6713507 0.7013627 1.0158485

image(1:500, 1:500, Gmatrix, main = "matriz G")

Amatrix = A[-c(1:20), -c(1:20)]
 dim(Amatrix)
#[1] 500 500
Amatrix[1:10,1:10]
#     I1  I2  I3  I4  I5  I6  I7  I8  I9 I10
#I1  1.0 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
#I2  0.5 1.0 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
#I3  0.5 0.5 1.0 0.5 0.5 0.5 0.5 0.5 0.5 0.5
#I4  0.5 0.5 0.5 1.0 0.5 0.5 0.5 0.5 0.5 0.5
#I5  0.5 0.5 0.5 0.5 1.0 0.5 0.5 0.5 0.5 0.5
#I6  0.5 0.5 0.5 0.5 0.5 1.0 0.5 0.5 0.5 0.5
#I7  0.5 0.5 0.5 0.5 0.5 0.5 1.0 0.5 0.5 0.5
#I8  0.5 0.5 0.5 0.5 0.5 0.5 0.5 1.0 0.5 0.5
#I9  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 1.0 0.5
#I10 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 1.0

image(1:500, 1:500, Amatrix, main = "matriz A")

#### correlación entre A y G
cor(as.vector(Amatrix), as.vector(Gmatrix))
[1] 0.9881251
plot(as.vector(Amatrix), as.vector(Gmatrix))


pedG = cbind(as.vector(Amatrix), as.vector(Gmatrix))
colnames(pedG) = c("A", "G")
head(pedG)
boxplot(G ~ A, data = pedG)


######################								######################
###################### CALCULO DE HEREDABILIDAD DEL CARACTER SIMULADO 	###################### 
######################								######################
str(simP)
#List of 3
 #$ pheno   : num [1:500, 1] 31.8 28.5 30.4 30.9 29.3 ...
 #$ QTN     : int [1:50] 516 250 349 185 218 480 936 177 821 389 ...
 #$ Meffects: num [1:50] -0.2742 1.2245 -0.2029 0.0898 0.2292 ...

# nuestra simulación está basada en:  y = mu + bX + e

y = simP$pheno
X = simulatedFS[, simP$QTN]
b = simP$Meffects

# Calcular valores genéticos
g <- X %*% b

# Calcular varianzas y heredabilidad
sigma2_g <- var(g)
sigma2_p <- var(y)
h2 <- sigma2_g / sigma2_p

print(paste("Heredabilidad estimada:", round(h2, 3)))