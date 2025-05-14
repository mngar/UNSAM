##############################################
## Clase práctica de GWAS                   ##
## Docente invitado: Dr. Martin N. Garcia   ##
## Fecha: 26/05/2025                        ##
## Universidad Nacional de San Martin       ##
##############################################

#CARGA DE PAQUETES
library(simulMGF)
library(tidyverse)
library(nortest)
library(rMVP)
library(data.table)
library(CJAMP)


#setear la semilla (para que los resultados puedan reproducirse exactamente)
set.seed(1234)


#SIMULACION DE DATOS
Nind <- 10000                   #numero de individuos
Nmarkers <- 10000               #numero de marcadores
Nqtl <- 50                      #numero de QTLs
Esigma <- .5                    #desvio estandard en la distribucion del caracter
Pmean <- 25                     #media del caracter
Perror <- .25                   #error debido al "ambiente"

simulN(Nind, Nmarkers, Nqtl, Esigma, Pmean, Perror)         #funcion para simular un caracter con distribucion normal
str(nsimout)                                                #objeto con la salida de la simulacion

nsimout$geno[1:10, 1:10]
nsimout$pheno[1:10]
nsimout$QTN
nsimout$Meffects

#asignamos nombres a los marcadores (M1, M2, M3,...,M10000) y a los individuos (I1, I2, I3,...,I10000)
pop <- as.data.frame(nsimout$geno)
colnames(pop) <- paste0("M",c(1:10000))
rownames(pop) <- paste0("I",c(1:10000))

pop[1:10, 1:10]


#matrices de datos genomicos y fenotipicos
#3 poblaciones de estudio (200, 1000 y 5000 individuos)
geno1 <- pop[1:200,]
geno2 <- pop[1:1000,]
geno3 <- pop[1:5000,]

feno1 <- nsimout$pheno[1:200]
feno2 <- nsimout$pheno[1:1000]
feno3 <- nsimout$pheno[1:5000]

feno1 <- data.frame(IND = rownames(geno1), Pheno = feno1)
feno2 <- data.frame(IND = rownames(geno2), Pheno = feno2)
feno3 <- data.frame(IND = rownames(geno3), Pheno = feno3)


#1 poblacion de validacion (5000 individuos, distintos de los de las poblaciones de estudio)
genoV <- pop[5001:10000,]
fenoV <- nsimout$pheno[5001:10000]
fenoV <- data.frame(IND = rownames(genoV), Pheno = fenoV)


#tablas con los QTLs (SNPs asociados)
QTL <- cbind(nsimout$QTN, nsimout$Meffects)
QTL <- cbind(QTL,abs(nsimout$Meffects))
colnames(QTL) <- c("marker", "effect", "effabs")
QTL <- as.data.frame(QTL)
QTL <- QTL[order(-QTL$effabs),]
QTL$SNP = paste0("M",QTL$marker)

head(QTL)
tail(QTL)


#mapa (ubicacion de los marcadores en 5 cromosoma)
map <- data.frame (SNP = colnames(pop),
                   Chromosome = c(rep(1,(Nmarkers/5)),
                                  rep(2,(Nmarkers/5)),
                                  rep(3,(Nmarkers/5)),
                                  rep(4,(Nmarkers/5)),
                                  rep(5,(Nmarkers/5))),
                   Position = c(1:(Nmarkers/5),
                                1:(Nmarkers/5),
                                1:(Nmarkers/5),
                                1:(Nmarkers/5),
                                1:(Nmarkers/5)))

head(map)
tail(map)

#PRUEBAS DE NORMALIDAD
#CHECKEO VISUAL

# HISTOGRAMA Y CURVA NORMAL: Consiste en representar los datos mediante un histograma
#y superponer la curva de una distribucion normal con la misma media y desviacion estandar
#que muestran los datos.

ggplot(data = feno1, aes(x = Pheno)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 2, color = "black", fill = "lightgrey") +
  scale_fill_gradient(low = "#DCDCDC", high = "#7C7C7C") +
  stat_function(fun = dnorm, colour = "firebrick",
                args = list(mean = mean(feno1$Pheno),
                            sd = sd(feno1$Pheno))) +
  ggtitle("Histograma + curva normal teorica feno1") +
  theme_bw()


ggplot(data = feno2, aes(x = Pheno)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 2, color = "black", fill = "lightgrey") +
  scale_fill_gradient(low = "#DCDCDC", high = "#7C7C7C") +
  stat_function(fun = dnorm, colour = "firebrick",
                args = list(mean = mean(feno2$Pheno),
                            sd = sd(feno2$Pheno))) +
  ggtitle("Histograma + curva normal teorica feno2") +
  theme_bw()


ggplot(data = feno3, aes(x = Pheno)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 2, color = "black", fill = "lightgrey") +
  scale_fill_gradient(low = "#DCDCDC", high = "#7C7C7C") +
  stat_function(fun = dnorm, colour = "firebrick",
                args = list(mean = mean(feno3$Pheno),
                            sd = sd(feno3$Pheno))) +
  ggtitle("Histograma + curva normal teorica feno3") +
  theme_bw()


ggplot(data = fenoV, aes(x = Pheno)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 2, color = "black", fill = "lightgrey") +
  scale_fill_gradient(low = "#DCDCDC", high = "#7C7C7C") +
  stat_function(fun = dnorm, colour = "firebrick",
                args = list(mean = mean(fenoV$Pheno),
                            sd = sd(fenoV$Pheno))) +
  ggtitle("Histograma + curva normal teorica fenoV") +
  theme_bw()



#Grafico de cuantiles teoricos (QQplot)
#Consiste en comparar los cuantiles de la distribucion observada con los cuantiles
#teoricos de una distribucion normal con la misma media y desvio estandar que los datos.
#Cuanto mas se aproximen los datos a una normal, mas alineados estan los puntos entorno a la
#recta.
qqnorm(feno1$Pheno, pch = 19, col = "gray50", main = "QQplot - feno1")
qqline(feno1$Pheno)

qqnorm(feno2$Pheno, pch = 19, col = "gray50", main = "QQplot - feno2")
qqline(feno2$Pheno)

qqnorm(feno3$Pheno, pch = 19, col = "gray50", main = "QQplot - feno3")
qqline(feno3$Pheno)

qqnorm(fenoV$Pheno, pch = 19, col = "gray50", main = "QQplot - fenoV")
qqline(fenoV$Pheno)


#CHECKEO ANALITICO

#test  Lilliefors: modificacion del  Kolmogorov-Smirnov.
#El test Lilliefors asume que la media y varianza son desconocidas, estando especialmente desarrollado
#para contrastar la normalidad. Es la alternativa al test de Shapiro-Wilk cuando el numero de
#observaciones es mayor de 50 (***). 

#La funcion lillie.test() del paquete nortest permite aplicarlo.   ***
lillie.test(x = feno1$Pheno)
lillie.test(x = feno2$Pheno)
lillie.test(x = feno3$Pheno)
lillie.test(x = fenoV$Pheno)

#Si el valor de p es mayor a 0.05 por lo que podemos decir que nuestros datos siguen una distribución normal.

# ***
# Si bien se puede recomendar el test de Lilliefors para tamaños de muestra mayores a 50 debido a suponer una media y varianza desconocidas,
# vale la pena señalar que típicamente exhibe una menor potencia en comparación con el test de Shapiro-Wilk. También se ha demostrado que
# el test de Shapiro-Wilk es el más potente, independientemente de la distribución y el tamaño de la muestra (Mendes & Pala, 2003). 
# Inicialmente, el test de Shapiro-Wilk, introducido por Shapiro y Wilk (1965), tenía una limitación en cuanto al tamaño de muestra, 
# siendo típicamente aplicable a muestras menores de 50. Sin embargo, la función shapiro.test() en R incorpora avances. Utiliza un 
# algoritmo mejorado, que incluye aproximaciones refinadas de pesos y el algoritmo AS R94 desarrollado por Royston en 1995 (Royston, 1995). 
# Esta mejora permite que el test se aplique en un espectro más amplio de tamaños de muestra, que va de 3 a 5000, como se detalla en la 
# documentación de R para la función:
#  https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/shapiro.test.

shapiro.test(x = feno1$Pheno)
shapiro.test(x = feno2$Pheno)
shapiro.test(x = feno3$Pheno)
shapiro.test(x = fenoV$Pheno)

#Nuevamente, si el valor de p es mayor a 0.05 por lo que podemos decir que nuestros datos siguen una distribución normal.


##########
## GWAS ##
##########

###################################
### Poblacion de estudio "pop1" ###
###################################

dir.create(file.path("E:/pp/", "pop1"), showWarnings = FALSE)       #creamos un directorio de trabajo, adonde iran los archivos de resultados
setwd("E:/pp/pop1")                                                 #seteamos el directorio de trabajo que generamos

imMVP1 <- MVP(                                                      #funcion para GWAS
  phe=feno1,          #Permite los datos perdidos en fenotipo (NA)  #precisa que coincidan las dimensiones con el nro. de individuos de geno
  geno=as.big.matrix(geno1),			                                  #solo matriz numerica y tiene que ser clase big.matrix
  map=map,                #mapa
  nPC.GLM=5,              #Los datos que hemos simulado no tienen estructura genetica, sin embargo
  nPC.MLM=3,              #  los autores del paquete de R recomiendan utilizar 3 o 5, para el 
  nPC.FarmCPU=3,          #  parametro nPC.
  maxLine=10000,          #Mientras menor sea el valor de maxLine (numero de marcadores evaluados a la vez) menor sera el costo en memoria.
  vc.method="BRENT",      #Metodo para la estimacion de los componentes de varianza (solo funciona para MLM)
  method.bin="static",    # "FaST-LMM", "static" (#only works for FarmCPU)
  threshold=0.05,         #Umbral para la correccion de Bonferroni.
  method=c("GLM", "MLM", "FarmCPU"),
  file.output=c("pmap", "pmap.signal", "plot", "log")
)

#MVP.Report(imMVP1,plot.type="q",col=c("dodgerblue1", "olivedrab3", "darkgoldenrod1"),threshold=1e6,
#           signal.pch=19,signal.cex=1.5,signal.col="red", # conf.int.col = "grey",
#           box=FALSE,multracks=TRUE, file.type="jpg",memo="",dpi=300)

###################################
### Poblacion de estudio "pop2" ###
###################################
dir.create(file.path("E:/pp/", "pop2"), showWarnings = FALSE)
setwd("E:/pp/pop2")
imMVP2 <- MVP(
  phe=feno2,      
  geno=as.big.matrix(geno2),		
  map=map,
  nPC.GLM=5,              
  nPC.MLM=3,              
  nPC.FarmCPU=3,
  maxLine=10000,          
  #ncpus=10,
  vc.method="BRENT",      
  method.bin="static",    
  threshold=0.05,
  method=c("GLM", "MLM", "FarmCPU"),
  file.output=c("pmap", "pmap.signal", "plot", "log")
)

###################################
### Poblacion de estudio "pop3" ###
###################################

dir.create(file.path("E:/pp/", "pop3"), showWarnings = FALSE)
setwd("E:/pp/pop3")
imMVP3 <- MVP(
  phe=feno3,      
  geno=as.big.matrix(geno3),		
  map=map,
  nPC.GLM=5,             
  nPC.MLM=3,             
  nPC.FarmCPU=3,
  maxLine=10000,         
  #ncpus=10,
  vc.method="BRENT",     
  method.bin="static",   
  threshold=0.05,
  method=c("GLM", "MLM", "FarmCPU"),
  file.output=c("pmap", "pmap.signal", "plot", "log")
)
#################################
### MODELOS PREDICTIVOS     #####
#################################

#################################
#                               #
#             GLM               #
#     MODELO LINEAL GENERAL     #
#################################

# pop1
# GLM (ver archivo signals): "Pheno.GLM_signals.CSV"
setwd("E:/pp/pop1")
glm.signal1 <- read.csv("Pheno.GLM_signals.CSV")
glm.signal1
#SNP Chromosome Position    MAF    Effect        SE    Pheno.GLM
#1 M6087          4       87 0.5000 -1.052039 0.2152592 2.143759e-06
#2 M7894          4     1894 0.4550  1.288934 0.2073545 3.081871e-09
#3 M9846          5     1846 0.4975  1.020128 0.2103834 2.548702e-06  
  


# cuanto explican los marcadores de la variacion del fenotipo?
glm.res1 = as.data.frame(imMVP1$glm.results)
caus.glm1 = c(rep(FALSE, 10000))
caus.glm1[c(6087, 7894, 9846)] = TRUE

pev.glm1 = compute_expl_var(genodata = geno1, phenodata = feno1$Pheno,
                           type = c("Rsquared_unadj", "Rsquared_adj"),
                           causal_idx = caus.glm1, effect_causal = glm.res1$Effect)

pev.glm1
#$Rsquared_unadj
#[1] 0.3313
#$Rsquared_adj
#[1] 0.3210648


# prediccion en popV con glm de pop1
e.glm1 = glm.signal1$Effect
e.glm1 = as.vector(e.glm1)
e.glm1

glm.geno.causal1 = as.matrix(genoV[,glm.signal1$SNP])
head(glm.geno.causal1)


yPred.glm1 = glm.geno.causal1%*%e.glm1


plot(fenoV$Pheno, yPred.glm1)


yPred.glm1 = glm.geno.causal1%*%e.glm1+mean(feno1$Pheno)
plot(fenoV$Pheno, yPred.glm1)



pdf("observado_vs_predicho_GLM_pop1.pdf")
plot(fenoV$Pheno, yPred.glm1)
dev.off()






# pop2
# GLM (ver archivo signals): "Pheno.GLM_signals.CSV"
setwd("E:/pp/pop2")
glm.signal2 <- read.csv("Pheno.GLM_signals.CSV")
glm.signal2
#     SNP Chromosome Position    MAF     Effect         SE    Pheno.GLM
#1   M469          1      469 0.4635 -0.4468979 0.09541521 3.208747e-06
#2  M1272          1     1272 0.4785  0.6595759 0.09037504 5.964966e-13
#3  M1922          1     1922 0.4940  0.5207612 0.09123793 1.511035e-08
#4  M4011          3       11 0.4925 -0.4737641 0.09343551 4.732314e-07
#5  M4253          3      253 0.4930 -0.9056890 0.08982566 7.908194e-23
#6  M5630          3     1630 0.4800  0.6437930 0.09133178 3.370590e-12
#7  M5946          3     1946 0.4915 -0.5062866 0.09314626 6.881072e-08
#8  M6087          4       87 0.4970 -0.6067330 0.09268383 9.442424e-11
#9  M6251          4      251 0.4820 -0.4582316 0.09214247 7.761080e-07
#10 M6252          4      252 0.4830  0.6835924 0.09042183 9.148296e-14
#11 M7894          4     1894 0.4925  1.0874446 0.08779682 7.246433e-33
#12 M8271          5      271 0.4805 -0.5890695 0.09242924 2.829479e-10
#13 M8544          5      544 0.4895  0.4751016 0.09351943 4.501563e-07
#14 M8550          5      550 0.4875  0.6058024 0.09208485 7.663113e-11
#15 M9452          5     1452 0.4805  0.4507586 0.09250547 1.280202e-06
#16 M9650          5     1650 0.4935 -0.8652863 0.09003830 5.721880e-21
#17 M9750          5     1750 0.4920  0.5819160 0.09292465 5.641308e-10
#18 M9846          5     1846 0.4905  0.7207329 0.08943801 2.205626e-15



# cuanto explican los marcadores de la variacion del fenotipo?
glm.res2 = as.data.frame(imMVP2$glm.results)
caus.glm2 = c(rep(FALSE, 10000))
pre.caus.glm2 <- as.integer(substr(glm.signal2$SNP, start = 2, stop = 5))

caus.glm2[pre.caus.glm2] = TRUE

pev.glm2 = compute_expl_var(genodata = geno2, phenodata = feno2$Pheno,
                            type = c("Rsquared_unadj", "Rsquared_adj"),
                            causal_idx = caus.glm2, effect_causal = glm.res2$Effect)

pev.glm2
#$Rsquared_unadj
#[1] 0.7829765
#$Rsquared_adj
#[1] 0.7789944


# prediccion en popV con glm de pop2
e.glm2 = glm.signal2$Effect
e.glm2 = as.vector(e.glm2)
e.glm2

glm.geno.causal2 = as.matrix(genoV[,glm.signal2$SNP])
head(glm.geno.causal2)


yPred.glm2 = glm.geno.causal2%*%e.glm2+mean(feno2$Pheno)


plot(fenoV$Pheno, yPred.glm2)

pdf("observado_vs_predicho_GLM_pop2.pdf")
plot(fenoV$Pheno, yPred.glm2+mean(feno2$Pheno))
dev.off()








# pop3
# GLM (ver archivo signals): "Pheno.GLM_signals.CSV"
setwd("E:/pp/pop3")
glm.signal3 <- read.csv("Pheno.GLM_signals.CSV")
glm.signal3
#     SNP Chromosome Position    MAF     Effect         SE     Pheno.GLM
#1   M469          1      469 0.4857 -0.5079588 0.04234254  1.040945e-32
#2   M830          1      830 0.4967 -0.2972692 0.04239335  2.659341e-12
#3  M1211          1     1211 0.4995 -0.2766921 0.04208237  5.362124e-11
#4  M1272          1     1272 0.4928  0.6916407 0.04140813  5.487073e-61
#5  M1922          1     1922 0.4813  0.4392673 0.04176075  1.310822e-25
#6  M2939          2      939 0.4982  0.4820710 0.04215266  6.472213e-30
#7  M3235          2     1235 0.4991 -0.3740134 0.04251368  1.897247e-18
#8  M3698          2     1698 0.4920  0.3803445 0.04240164  4.117332e-19
#9  M4011          3       11 0.4981 -0.5073903 0.04235899  1.286947e-32
#10 M4253          3      253 0.4986 -0.9496993 0.04045899 1.149763e-115
#11 M4920          3      920 0.4983 -0.3808042 0.04272285  6.821836e-19
#12 M4983          3      983 0.4966  0.3020916 0.04236178  1.136989e-12
#13 M5352          3     1352 0.4885  0.2600892 0.04251439  1.021907e-09
#14 M5630          3     1630 0.4983  0.5711802 0.04184869  1.140597e-41
#15 M5702          3     1702 0.4944  0.3815839 0.04226456  2.439171e-19
#16 M5934          3     1934 0.4975  0.3129915 0.04223205  1.461681e-13
#17 M5946          3     1946 0.4977 -0.5107098 0.04177292  6.873444e-34
#18 M5985          3     1985 0.4976  0.2573808 0.04252009  1.524095e-09
#19 M6087          4       87 0.4956 -0.4777585 0.04224546  2.680015e-29
#20 M6251          4      251 0.4952 -0.3486503 0.04234162  2.285624e-16
#21 M6252          4      252 0.4927  0.7994052 0.04094801  7.454093e-82
#22 M7209          4     1209 0.4930 -0.2311973 0.04279486  6.878571e-08
#23 M7465          4     1465 0.4989  0.3334116 0.04218695  3.320029e-15
#24 M7894          4     1894 0.4944  1.1075419 0.03985735 3.838147e-158
#25 M8048          5       48 0.4964  0.2874223 0.04218258  1.063268e-11
#26 M8271          5      271 0.4991 -0.5460361 0.04184730  2.731654e-38
#27 M8535          5      535 0.4992 -0.1967981 0.04267233  4.090855e-06
#28 M8550          5      550 0.4928  0.5150039 0.04221293  9.342822e-34
#29 M8591          5      591 0.4952  0.2330355 0.04215835  3.410928e-08
#30 M9080          5     1080 0.4981 -0.5507267 0.04202528  1.334156e-38
#31 M9452          5     1452 0.4840  0.3802989 0.04219305  2.794834e-19
#32 M9650          5     1650 0.4968 -0.8507723 0.04087727  2.492441e-92
#33 M9750          5     1750 0.4938  0.6255445 0.04174176  1.075125e-49
#34 M9846          5     1846 0.4960  0.5435948 0.04177919  4.375084e-38




# cuanto explican los marcadores de la variacion del fenotipo?
glm.res3 = as.data.frame(imMVP3$glm.results)
caus.glm3 = c(rep(FALSE, 10000))
pre.caus.glm3 <- as.integer(substr(glm.signal3$SNP, start = 2, stop = 5))

caus.glm3[pre.caus.glm3] = TRUE

pev.glm3 = compute_expl_var(genodata = geno3, phenodata = feno3$Pheno,
                            type = c("Rsquared_unadj", "Rsquared_adj"),
                            causal_idx = caus.glm3, effect_causal = glm.res3$Effect)

pev.glm3
#$Rsquared_unadj
#[1] 0.9772662
#$Rsquared_adj
#[1] 0.9771105


# prediccion en popV con glm de pop2
e.glm3 = glm.signal3$Effect
e.glm3 = as.vector(e.glm3)
e.glm3

glm.geno.causal3 = as.matrix(genoV[,glm.signal3$SNP])
head(glm.geno.causal3)


yPred.glm3 = glm.geno.causal3%*%e.glm3+mean(feno3$Pheno)


plot(fenoV$Pheno, yPred.glm3)

pdf("observado_vs_predicho_GLM_pop3.pdf")
plot(fenoV$Pheno, yPred.glm3)
dev.off()


par(mfrow=c(1,3))
plot(fenoV$Pheno, yPred.glm1)
plot(fenoV$Pheno, yPred.glm2)
plot(fenoV$Pheno, yPred.glm3)

cor(fenoV$Pheno, yPred.glm1)
cor(fenoV$Pheno, yPred.glm2)
cor(fenoV$Pheno, yPred.glm3)


resultados.glm = cbind(fenoV$Pheno, yPred.glm1)
resultados.glm = cbind(resultados.glm, yPred.glm2)
resultados.glm = cbind(resultados.glm, yPred.glm3)

colnames(resultados.glm) = c("validacion", "glm1", "glm2", "glm3")
resultados.glm = as.data.frame(resultados.glm)

#si seleccionaramos el 25% mayor de acuerdo a cada modelo predictivo (por encima de linea azul)

par(mfrow=c(1,3))
plot(resultados.glm[,c(1,2)])
abline(h=quantile(resultados.glm$glm1,0.75),v=quantile(resultados.glm$validacion,0.75), col=c("blue","red"))
plot(resultados.glm[,c(1,3)])
abline(h=quantile(resultados.glm$glm2,0.75),v=quantile(resultados.glm$validacion,0.75), col=c("blue","red"))
plot(resultados.glm[,c(1,4)])
abline(h=quantile(resultados.glm$glm3,0.75),v=quantile(resultados.glm$validacion,0.75), col=c("blue","red"))

#Eficiencia de selección

ind.sup = rownames(resultados.glm[which(resultados.glm$validacion>quantile(resultados.glm$validacion,0.75)),])
sel.glm1 = rownames(resultados.glm[which(resultados.glm$glm1>quantile(resultados.glm$glm1,0.75)),])
sel.glm2 = rownames(resultados.glm[which(resultados.glm$glm2>quantile(resultados.glm$glm2,0.75)),])
sel.glm3 = rownames(resultados.glm[which(resultados.glm$glm3>quantile(resultados.glm$glm3,0.75)),])
  
ef.glm1 = length(intersect(ind.sup,sel.glm1))/length(sel.glm1)
ef.glm2 = length(intersect(ind.sup,sel.glm2))/length(sel.glm2)
ef.glm3 = length(intersect(ind.sup,sel.glm3))/length(sel.glm3)

ef.glm1
#[1] 0.4390
ef.glm2
#[1] 0.7144
ef.glm3
#[1] 0.9096


#poder estadistico
poder.glm1 = length(unique(glm.signal1$SNP,QTL$SNP))/Nqtl
poder.glm2 = length(unique(glm.signal2$SNP,QTL$SNP))/Nqtl
poder.glm3 = length(unique(glm.signal3$SNP,QTL$SNP))/Nqtl

poder.glm1
poder.glm2
poder.glm3

#falsos positivos
falpos.glm1 = abs(length(unique(glm.signal1$SNP,QTL$SNP))-length(glm.signal1$SNP))
falpos.glm2 = abs(length(unique(glm.signal2$SNP,QTL$SNP))-length(glm.signal2$SNP))
falpos.glm3 = abs(length(unique(glm.signal3$SNP,QTL$SNP))-length(glm.signal3$SNP))

falpos.glm1
falpos.glm2
falpos.glm3



###################################################################### 
### Cual es el minimo tamaño de efecto detectado por GLM en cada caso
#200 ind
mergexsnp_glm1_QTL = merge(glm.signal1,QTL, by = "SNP")
head(mergexsnp_glm1_QTL)
min(mergexsnp_glm1_QTL$effabs)
#[1] 0.4257824
#1000 ind
mergexsnp_glm2_QTL = merge(glm.signal2,QTL, by = "SNP")
head(mergexsnp_glm2_QTL)
min(mergexsnp_glm2_QTL$effabs)
#[1] 0.3326874
#5000 ind
mergexsnp_glm3_QTL = merge(glm.signal3,QTL, by = "SNP")
head(mergexsnp_glm3_QTL)
min(mergexsnp_glm3_QTL$effabs)
#[1] 0.1557778
#######################################################################




#################################
#                               #
#             MLM               #
#      MODELO LINEAL MIXTO      #
#################################
#
# pop1
# MLM (ver archivo signals): "Pheno.MLM_signals.CSV"
setwd("E:/pp/pop1")
mlm.signal1 <- read.csv("Pheno.MLM_signals.CSV")
mlm.signal1

# cuanto explican los marcadores de la variacion del fenotipo?
mlm.res1 = as.data.frame(imMVP1$mlm.results)
caus.mlm1 = c(rep(FALSE, 10000))
pre.caus.mlm1 <- as.integer(substr(mlm.signal1$SNP, start = 2, stop = 5))
caus.mlm1[pre.caus.mlm1] = TRUE


pev.mlm1 = compute_expl_var(genodata = geno1, phenodata = feno1$Pheno,
                            type = c("Rsquared_unadj", "Rsquared_adj"),
                            causal_idx = caus.mlm1, effect_causal = mlm.res1$Effect)

pev.mlm1

# prediccion en popV con glm de pop1
e.mlm1 = mlm.signal1$Effect
e.mlm1 = as.vector(e.mlm1)
e.mlm1

mlm.geno.causal1 = as.matrix(genoV[,mlm.signal1$SNP])
head(mlm.geno.causal1)


yPred.mlm1 = mlm.geno.causal1%*%e.mlm1


plot(fenoV$Pheno, yPred.mlm1)


yPred.mlm1 = mlm.geno.causal1%*%e.mlm1+mean(feno1$Pheno)
plot(fenoV$Pheno, yPred.mlm1)



pdf("observado_vs_predicho_MLM_pop1.pdf")
plot(fenoV$Pheno, yPred.mlm1)
dev.off()






# pop2
# MLM (ver archivo signals): "Pheno.MLM_signals.CSV"
setwd("E:/pp/pop2")
mlm.signal2 <- read.csv("Pheno.MLM_signals.CSV")
mlm.signal2

# cuanto explican los marcadores de la variacion del fenotipo?
mlm.res2 = as.data.frame(imMVP2$mlm.results)
caus.mlm2 = c(rep(FALSE, 10000))
pre.caus.mlm2 <- as.integer(substr(mlm.signal2$SNP, start = 2, stop = 5))

caus.mlm2[pre.caus.mlm2] = TRUE

pev.mlm2 = compute_expl_var(genodata = geno2, phenodata = feno2$Pheno,
                            type = c("Rsquared_unadj", "Rsquared_adj"),
                            causal_idx = caus.mlm2, effect_causal = mlm.res2$Effect)

pev.mlm2

# prediccion en popV con glm de pop2
e.mlm2 = mlm.signal2$Effect
e.mlm2 = as.vector(e.mlm2)
e.mlm2

mlm.geno.causal2 = as.matrix(genoV[,mlm.signal2$SNP])
head(mlm.geno.causal2)


yPred.mlm2 = mlm.geno.causal2%*%e.mlm2+mean(feno2$Pheno)


plot(fenoV$Pheno, yPred.mlm2)

pdf("observado_vs_predicho_MLM_pop2.pdf")
plot(fenoV$Pheno, yPred.mlm2+mean(feno2$Pheno))
dev.off()

# pop3
# MLM (ver archivo signals): "Pheno.MLM_signals.CSV"
setwd("E:/pp/pop3")
mlm.signal3 <- read.csv("Pheno.MLM_signals.CSV")
mlm.signal3

# cuanto explican los marcadores de la variacion del fenotipo?
mlm.res3 = as.data.frame(imMVP3$mlm.results)
caus.mlm3 = c(rep(FALSE, 10000))
pre.caus.mlm3 <- as.integer(substr(mlm.signal3$SNP, start = 2, stop = 5))

caus.mlm3[pre.caus.mlm3] = TRUE

pev.mlm3 = compute_expl_var(genodata = geno3, phenodata = feno3$Pheno,
                            type = c("Rsquared_unadj", "Rsquared_adj"),
                            causal_idx = caus.mlm3, effect_causal = mlm.res3$Effect)

pev.mlm3

# prediccion en popV con glm de pop2
e.mlm3 = mlm.signal3$Effect
e.mlm3 = as.vector(e.mlm3)
e.mlm3

mlm.geno.causal3 = as.matrix(genoV[,mlm.signal3$SNP])
head(mlm.geno.causal3)


yPred.mlm3 = mlm.geno.causal3%*%e.mlm3+mean(feno3$Pheno)


plot(fenoV$Pheno, yPred.mlm3)

pdf("observado_vs_predicho_MLM_pop3.pdf")
plot(fenoV$Pheno, yPred.mlm3)
dev.off()


par(mfrow=c(1,3))
plot(fenoV$Pheno, yPred.mlm1)
plot(fenoV$Pheno, yPred.mlm2)
plot(fenoV$Pheno, yPred.mlm3)

cor(fenoV$Pheno, yPred.mlm1)
cor(fenoV$Pheno, yPred.mlm2)
cor(fenoV$Pheno, yPred.mlm3)


resultados.mlm = cbind(fenoV$Pheno, yPred.mlm1)
resultados.mlm = cbind(resultados.mlm, yPred.mlm2)
resultados.mlm = cbind(resultados.mlm, yPred.mlm3)

colnames(resultados.mlm) = c("validacion", "mlm1", "mlm2", "mlm3")
resultados.mlm = as.data.frame(resultados.mlm)

#si seleccionaramos el 25% mayor de acuerdo a cada modelo predictivo (por encima de linea azul)

par(mfrow=c(1,3))
plot(resultados.mlm[,c(1,2)])
abline(h=quantile(resultados.mlm$mlm1,0.75),v=quantile(resultados.mlm$validacion,0.75), col=c("blue","red"))
plot(resultados.mlm[,c(1,3)])
abline(h=quantile(resultados.mlm$mlm2,0.75),v=quantile(resultados.mlm$validacion,0.75), col=c("blue","red"))
plot(resultados.mlm[,c(1,4)])
abline(h=quantile(resultados.mlm$mlm3,0.75),v=quantile(resultados.mlm$validacion,0.75), col=c("blue","red"))

#Eficiencia de selección

ind.sup = rownames(resultados.mlm[which(resultados.mlm$validacion>quantile(resultados.mlm$validacion,0.75)),])
sel.mlm1 = rownames(resultados.mlm[which(resultados.mlm$mlm1>quantile(resultados.mlm$mlm1,0.75)),])
sel.mlm2 = rownames(resultados.mlm[which(resultados.mlm$mlm2>quantile(resultados.mlm$mlm2,0.75)),])
sel.mlm3 = rownames(resultados.mlm[which(resultados.mlm$mlm3>quantile(resultados.mlm$mlm3,0.75)),])

ef.mlm1 = length(intersect(ind.sup,sel.mlm1))/length(sel.mlm1)
ef.mlm2 = length(intersect(ind.sup,sel.mlm2))/length(sel.mlm2)
ef.mlm3 = length(intersect(ind.sup,sel.mlm3))/length(sel.mlm3)

ef.mlm1
#[1] NaN -> en este caso ~33% de la poblacion esta en el puesto 1, debido a que se utilizo un unico marcador para generar el modelo
ef.mlm2
ef.mlm3


#poder estadistico
poder.mlm1 = length(unique(mlm.signal1$SNP,QTL$SNP))/Nqtl
poder.mlm2 = length(unique(mlm.signal2$SNP,QTL$SNP))/Nqtl
poder.mlm3 = length(unique(mlm.signal3$SNP,QTL$SNP))/Nqtl

poder.mlm1
poder.mlm2
poder.mlm3

#falsos positivos
falpos.mlm1 = abs(length(unique(mlm.signal1$SNP,QTL$SNP))-length(mlm.signal1$SNP))
falpos.mlm2 = abs(length(unique(mlm.signal2$SNP,QTL$SNP))-length(mlm.signal2$SNP))
falpos.mlm3 = abs(length(unique(mlm.signal3$SNP,QTL$SNP))-length(mlm.signal3$SNP))

falpos.mlm1
falpos.mlm2
falpos.mlm3



###################################################################### 
### Cual es el minimo tamaño de efecto detectado por MLM en cada caso
#200 ind
mergexsnp_mlm1_QTL = merge(mlm.signal1,QTL, by = "SNP")
head(mergexsnp_mlm1_QTL)
min(mergexsnp_mlm1_QTL$effabs)

#1000 ind
mergexsnp_mlm2_QTL = merge(mlm.signal2,QTL, by = "SNP")
head(mergexsnp_mlm2_QTL)
min(mergexsnp_mlm2_QTL$effabs)

#5000 ind
mergexsnp_mlm3_QTL = merge(mlm.signal3,QTL, by = "SNP")
head(mergexsnp_mlm3_QTL)
min(mergexsnp_mlm3_QTL$effabs)

#######################################################################








#################################
#                               #
#          FarmCPU              #
#                               #
#################################
# pop1
# FarmCPU (ver archivo signals): "Pheno.FarmCPU_signals.CSV"
setwd("E:/pp/pop1")
FarmCPU.signal1 <- read.csv("Pheno.FarmCPU_signals.CSV")
FarmCPU.signal1

# cuanto explican los marcadores de la variacion del fenotipo?
FarmCPU.res1 = as.data.frame(imMVP1$FarmCPU.results)
caus.FarmCPU1 = c(rep(FALSE, 10000))
pre.caus.FarmCPU1 <- as.integer(substr(FarmCPU.signal1$SNP, start = 2, stop = 5))

caus.FarmCPU1[pre.caus.FarmCPU1] = TRUE

pev.FarmCPU1 = compute_expl_var(genodata = geno1, phenodata = feno1$Pheno,
                            type = c("Rsquared_unadj", "Rsquared_adj"),
                            causal_idx = caus.FarmCPU1, effect_causal = FarmCPU.res1$Effect)

pev.FarmCPU1

# prediccion en popV con glm de pop1
e.FarmCPU1 = FarmCPU.signal1$Effect
e.FarmCPU1 = as.vector(e.FarmCPU1)
e.FarmCPU1

FarmCPU.geno.causal1 = as.matrix(genoV[,FarmCPU.signal1$SNP])
head(FarmCPU.geno.causal1)


yPred.FarmCPU1 = FarmCPU.geno.causal1%*%e.FarmCPU1

par(mfrow=c(1,1))
plot(fenoV$Pheno, yPred.FarmCPU1)


yPred.FarmCPU1 = FarmCPU.geno.causal1%*%e.FarmCPU1+mean(feno1$Pheno)
plot(fenoV$Pheno, yPred.FarmCPU1)



pdf("observado_vs_predicho_FarmCPU_pop1.pdf")
plot(fenoV$Pheno, yPred.FarmCPU1)
dev.off()






# pop2
# FarmCPU (ver archivo signals): "Pheno.FarmCPU_signals.CSV"
setwd("E:/pp/pop2")
FarmCPU.signal2 <- read.csv("Pheno.FarmCPU_signals.CSV")
FarmCPU.signal2

# cuanto explican los marcadores de la variacion del fenotipo?
FarmCPU.res2 = as.data.frame(imMVP2$FarmCPU.results)
caus.FarmCPU2 = c(rep(FALSE, 10000))
pre.caus.FarmCPU2 <- as.integer(substr(FarmCPU.signal2$SNP, start = 2, stop = 5))

caus.FarmCPU2[pre.caus.FarmCPU2] = TRUE

pev.FarmCPU2 = compute_expl_var(genodata = geno2, phenodata = feno2$Pheno,
                            type = c("Rsquared_unadj", "Rsquared_adj"),
                            causal_idx = caus.FarmCPU2, effect_causal = FarmCPU.res2$Effect)

pev.FarmCPU2

# prediccion en popV con glm de pop2
e.FarmCPU2 = FarmCPU.signal2$Effect
e.FarmCPU2 = as.vector(e.FarmCPU2)
e.FarmCPU2

FarmCPU.geno.causal2 = as.matrix(genoV[,FarmCPU.signal2$SNP])
head(FarmCPU.geno.causal2)


yPred.FarmCPU2 = FarmCPU.geno.causal2%*%e.FarmCPU2+mean(feno2$Pheno)


plot(fenoV$Pheno, yPred.FarmCPU2)

pdf("observado_vs_predicho_FarmCPU_pop2.pdf")
plot(fenoV$Pheno, yPred.FarmCPU2+mean(feno2$Pheno))
dev.off()








# pop3
# FarmCPU (ver archivo signals): "Pheno.FarmCPU_signals.CSV"
setwd("E:/pp/pop3")
FarmCPU.signal3 <- read.csv("Pheno.FarmCPU_signals.CSV")
FarmCPU.signal3

# cuanto explican los marcadores de la variacion del fenotipo?
FarmCPU.res3 = as.data.frame(imMVP3$FarmCPU.results)
caus.FarmCPU3 = c(rep(FALSE, 10000))
pre.caus.FarmCPU3 <- as.integer(substr(FarmCPU.signal3$SNP, start = 2, stop = 5))

caus.FarmCPU3[pre.caus.FarmCPU3] = TRUE

pev.FarmCPU3 = compute_expl_var(genodata = geno3, phenodata = feno3$Pheno,
                            type = c("Rsquared_unadj", "Rsquared_adj"),
                            causal_idx = caus.FarmCPU3, effect_causal = FarmCPU.res3$Effect)

pev.FarmCPU3

# prediccion en popV con glm de pop2
e.FarmCPU3 = FarmCPU.signal3$Effect
e.FarmCPU3 = as.vector(e.FarmCPU3)
e.FarmCPU3

FarmCPU.geno.causal3 = as.matrix(genoV[,FarmCPU.signal3$SNP])
head(FarmCPU.geno.causal3)


yPred.FarmCPU3 = FarmCPU.geno.causal3%*%e.FarmCPU3+mean(feno3$Pheno)


plot(fenoV$Pheno, yPred.FarmCPU3)

pdf("observado_vs_predicho_FarmCPU_pop3.pdf")
plot(fenoV$Pheno, yPred.FarmCPU3)
dev.off()


par(mfrow=c(1,3))
plot(fenoV$Pheno, yPred.FarmCPU1)
plot(fenoV$Pheno, yPred.FarmCPU2)
plot(fenoV$Pheno, yPred.FarmCPU3)

cor(fenoV$Pheno, yPred.FarmCPU1)
cor(fenoV$Pheno, yPred.FarmCPU2)
cor(fenoV$Pheno, yPred.FarmCPU3)


resultados.FarmCPU = cbind(fenoV$Pheno, yPred.FarmCPU1)
resultados.FarmCPU = cbind(resultados.FarmCPU, yPred.FarmCPU2)
resultados.FarmCPU = cbind(resultados.FarmCPU, yPred.FarmCPU3)

colnames(resultados.FarmCPU) = c("validacion", "FarmCPU1", "FarmCPU2", "FarmCPU3")
resultados.FarmCPU = as.data.frame(resultados.FarmCPU)

#si seleccionaramos el 25% mayor de acuerdo a cada modelo predictivo (por encima de linea azul)

par(mfrow=c(1,3))
plot(resultados.FarmCPU[,c(1,2)])
abline(h=quantile(resultados.FarmCPU$FarmCPU1,0.75),v=quantile(resultados.FarmCPU$validacion,0.75), col=c("blue","red"))
plot(resultados.FarmCPU[,c(1,3)])
abline(h=quantile(resultados.FarmCPU$FarmCPU2,0.75),v=quantile(resultados.FarmCPU$validacion,0.75), col=c("blue","red"))
plot(resultados.FarmCPU[,c(1,4)])
abline(h=quantile(resultados.FarmCPU$FarmCPU3,0.75),v=quantile(resultados.FarmCPU$validacion,0.75), col=c("blue","red"))

#Eficiencia de selección

ind.sup = rownames(resultados.FarmCPU[which(resultados.FarmCPU$validacion>quantile(resultados.FarmCPU$validacion,0.75)),])
sel.FarmCPU1 = rownames(resultados.FarmCPU[which(resultados.FarmCPU$FarmCPU1>quantile(resultados.FarmCPU$FarmCPU1,0.75)),])
sel.FarmCPU2 = rownames(resultados.FarmCPU[which(resultados.FarmCPU$FarmCPU2>quantile(resultados.FarmCPU$FarmCPU2,0.75)),])
sel.FarmCPU3 = rownames(resultados.FarmCPU[which(resultados.FarmCPU$FarmCPU3>quantile(resultados.FarmCPU$FarmCPU3,0.75)),])

ef.FarmCPU1 = length(intersect(ind.sup,sel.FarmCPU1))/length(sel.FarmCPU1)
ef.FarmCPU2 = length(intersect(ind.sup,sel.FarmCPU2))/length(sel.FarmCPU2)
ef.FarmCPU3 = length(intersect(ind.sup,sel.FarmCPU3))/length(sel.FarmCPU3)

ef.FarmCPU1
ef.FarmCPU2
ef.FarmCPU3

#poder estadistico
poder.FarmCPU1 = length(unique(FarmCPU.signal1$SNP,QTL$SNP))/Nqtl
poder.FarmCPU2 = length(unique(FarmCPU.signal2$SNP,QTL$SNP))/Nqtl
poder.FarmCPU3 = length(unique(FarmCPU.signal3$SNP,QTL$SNP))/Nqtl

poder.FarmCPU1
poder.FarmCPU2
poder.FarmCPU3

#falsos positivos
falpos.FarmCPU1 = abs(length(unique(FarmCPU.signal1$SNP,QTL$SNP))-length(FarmCPU.signal1$SNP))
falpos.FarmCPU2 = abs(length(unique(FarmCPU.signal2$SNP,QTL$SNP))-length(FarmCPU.signal2$SNP))
falpos.FarmCPU3 = abs(length(unique(FarmCPU.signal3$SNP,QTL$SNP))-length(FarmCPU.signal3$SNP))

falpos.FarmCPU1
falpos.FarmCPU2
falpos.FarmCPU3



###################################################################### 
### Cual es el minimo tamaño de efecto detectado por FarmCPU en cada caso
#200 ind
mergexsnp_FarmCPU1_QTL = merge(FarmCPU.signal1,QTL, by = "SNP")
head(mergexsnp_FarmCPU1_QTL)
min(mergexsnp_FarmCPU1_QTL$effabs)

#1000 ind
mergexsnp_FarmCPU2_QTL = merge(FarmCPU.signal2,QTL, by = "SNP")
head(mergexsnp_FarmCPU2_QTL)
min(mergexsnp_FarmCPU2_QTL$effabs)

#5000 ind
mergexsnp_FarmCPU3_QTL = merge(FarmCPU.signal3,QTL, by = "SNP")
head(mergexsnp_FarmCPU3_QTL)
min(mergexsnp_FarmCPU3_QTL$effabs)

#######################################################################





#BIBLIOGRAFIA RECOMENDADA
#  Tibbs Cortes, L., Zhang, Z., & Yu, J. (2021). Status and prospects of genome‐wide association studies in plants.
#  The plant genome, 14(1), e20077.

