#Librarys a carregar
library(knitr)
library(tidyverse)
library(kableExtra)
library(VIM)
library(gridExtra)
library(nortest)
library(corrplot)
library(reshape2)
library(car)
library(MASS)
library(vcd)
library(pROC)

# Importació del dataset
download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/00423/hcc-survival.zip",
              "../data/hcc_data.zip")
unzip("../data/hcc_data.zip",exdir = "../data")
hcc<-read.csv("../data/hcc-survival/hcc-data.txt",header = F, dec = ".",
              stringsAsFactors = F, na.strings = "?")
write.csv(hcc, file = "../data/hcc_data.csv")

# Noms de les variables
hcc_header<-c("Gender", "Symptoms", "Alcohol", "Hepatitis B Surface Antigen",
              "Hepatitis B e Antigen", "Hepatitis B Core Antibody",
              "Hepatitis C Virus Antibody",  "Cirrhosis", 
              "Endemic Countries","Smoking",  "Diabetes", "Obesity",
              "Hemochromatosis", "Arterial Hypertension", 
              "Chronic Renal Insufficiency", "Human Immunodeficiency Virus",
              "Nonalcoholic Steatohepatitis","Esophageal Varices",
              "Splenomegaly",  "Portal Hypertension", 
              "Portal Vein Thrombosis", "Liver Metastasis", 
              "Radiological Hallmark", "Age at diagnosis", 
              "Grams of Alcohol per day", "Packs of cigarets per year",
              "Performance Status", "Encephalopathy degree", "Ascites degree",
              "International Normalised Ratio", "Alpha-Fetoprotein (ng/mL)",
              "Haemoglobin (g/dL)", "Mean Corpuscular Volume",
              "Leukocytes(G/L)", "Platelets", "Albumin (mg/dL)",
              "Total Bilirubin(mg/dL)",  "Alanine transaminase (U/L)"
              , "Aspartate transaminase (U/L)",
              "Gamma glutamyl transferase (U/L)", 
              "Alkaline phosphatase (U/L)","Total Proteins (g/dL)",
              "Creatinine (mg/dL)","Number of Nodules",
              "Major dimension of nodule (cm)", 
              "Direct Bilirubin (mg/dL)", "Iron", "Oxygen Saturation (%)",
              "Ferritin (ng/mL)", "Class Attribute"   )

colnames(hcc)<-hcc_header

#Mirar de quin tipus son les variables carregades
tipVar <- c()
for (i in 1:ncol(hcc))  tipVar <- c(tipVar,is(hcc[,i])[1])
tipVar <- table(tipVar)

kable(x = tipVar, format = "latex", caption = "Tipus de dades carregades",
      booktabs = TRUE, linesep = '') %>% 
  kable_styling(latex_options = c("HOLD_position")) %>% 
  row_spec(row=0, bold = TRUE)

# Definició de tipus de dades per columna
hcc_factor <-c(1:23,50)
hcc_order <- c(27:29)
hcc_factorT<-c(1:23,27:29,50)
hcc_num<-c(24:26,30:49)

# Factorització de les columnes categòriques
hcc <- hcc %>% mutate_at(vars(c(1:23,50)), as.factor)
hcc <- hcc %>% mutate_at(vars(c(27:29)), as.factor)

#Valors nuls per pacient
tauNuls <- table(apply(hcc, 1, function(x) sum(is.na(x))))
kable(x = tauNuls, format = "latex", caption = "Valors nuls per pacient",
      booktabs = TRUE, linesep = '', col.names = c("#NA per pacient", "Freq")) %>% 
  kable_styling(latex_options = c("HOLD_position")) %>% 
  row_spec(row=0, bold = TRUE)

#Funció que calcula percentatge
funNA <- function(a, n){
  a = round(100*a/n,1)
}

#Calcul dels valors NA per variable
totNA <- hcc %>% dplyr::select(everything()) %>%
  summarise_all(funs(sum(is.na(.)))) 
perNA <- totNA %>% mutate_all(funNA, n= nrow(hcc))
tauNA <- totNA %>% bind_rows(perNA)
tauNA <- as_tibble(t(tauNA), rownames = "Variable") %>% 
  rename(`total NA` = V1, `%NA` = V2) %>% 
  arrange(-`total NA`)

#Taula en dos columnes
tau <- cbind(tauNA[1:25,],tauNA[26:50,])

#Taula amb resultats
kable(x = tau, format = "latex", caption = "Variables amb NA",
      booktabs = TRUE, linesep = '') %>% 
  kable_styling(latex_options = c("HOLD_position", "scale_down")) %>% 
  row_spec(row=0, bold = TRUE)
# Valoració de la distribució dels NA a les variables amb respecte la variable classe 
cero  <- hcc %>% filter(`Class Attribute` == 0)
uno  <- hcc %>% filter(`Class Attribute` == 1)

probTest <- tibble()
for (i in names(hcc[,hcc_factorT])) {
  if (sum(is.na(hcc[,i]))>0){
    casos<-c(sum(is.na(cero[,i])),sum(is.na(uno[,i])))
    long<-c(length(cero[,i]),length(uno[,i]))
    test<-prop.test(x=casos,n=long)
    probTest <- probTest %>%
      bind_rows(c(Clase = i, p_value = test$p.value,
                  prob_0=(casos/long)[1], prob_1=(casos/long)[2],
                  numNA_0=casos[1], numNA_1=casos[2]))
  }
}

#Filtrar aquelles que rebutjen l'hip nul·la
testSig <- probTest %>% filter(p_value <= 0.05) %>% 
  mutate_at(.vars = c("p_value", "prob_0","prob_1","numNA_0", "numNA_1"), as.numeric)

#Presentació de resultats
kable(x = testSig, format = "latex",
      caption = "Variables amb una diferència significativa de NA entre els grups de la classe",
      booktabs = TRUE, digits = 4, linesep = '') %>% 
  kable_styling(latex_options = c("HOLD_position")) %>% 
  row_spec(row=0, bold = TRUE)

#computem tots els valors NA de variables factor
factorNames<- colnames(hcc[hcc_factorT,])

# Computar per kNN, abm valors standars, k = 5
hcc <- kNN(hcc, variable = factorNames) %>% subset(select = Gender:`Class Attribute`)

# Correcció valors leucocitosi i plaquetes +
# assignació mediana als valors NA de variables numeriques
hcc <- hcc %>% 
  mutate(`Leukocytes(G/L)` = ifelse(`Leukocytes(G/L)` > 100,`Leukocytes(G/L)`/100, `Leukocytes(G/L)`),
         Platelets = ifelse(Platelets < 1000,Platelets*1000, Platelets))

for (i in hcc_num){
  hcc[,i] <- ifelse(is.na(hcc[,i]),median(hcc[,i], na.rm = T),hcc[,i])
}

# Exportació de les dades netes en .csv
hcc_clean <- hcc
write.csv(hcc_clean, file = "../data/hcc_data_clean.csv")

# Estudi distribució variables Numeriques
ind_CA <- which(colnames(hcc) == "Class Attribute")
figCapsNum <- names(hcc)[hcc_num]

# grafic 1 - histograma
# grafic 2 - boxplot
for (i in 1:length(hcc_num)) {
  data <- hcc[,c(hcc_num[i],ind_CA )]
  name_var <- names(data)
  
  a1<-data %>%
    ggplot(aes(x=data[,1], fill=`Class Attribute`))+
    geom_histogram() +
    labs(fill="Clase", y = "Frequencia", x =name_var[1] ) +
    theme(legend.position = "bottom")
  
  a2<-data %>% 
    ggplot(aes(x=`Class Attribute`,y=data[,1])) +
    geom_boxplot() +
    labs(x = "Clase",  y = name_var[1])
  
  grid.arrange(a1,a2,nrow=1)
  cat('\n\n')
}

# Transformació logarítmica de variables numèriques
ind_CA <- which(colnames(hcc) == "Class Attribute")
hcc_log<-c(31,37:43,46)
figCapsLog <- names(hcc)[hcc_log]

# grafic 1 - histograma
# grafic 2 - boxplot
for (i in 1:length(hcc_log)){
  #transformacio logaritimica
  hcc[,hcc_log[i]]<- log(hcc[,hcc_log[i]])
  colnames(hcc)[hcc_log[i]] <- paste("log_", colnames(hcc)[hcc_log[i]], sep = '')
  
  data <- hcc[,c(hcc_log[i],ind_CA )]
  name_var <- names(data) 
  
  a1<-data %>%
    ggplot(aes(x=data[,1], fill=`Class Attribute`))+
    geom_histogram() +
    labs(fill="Clase", y = "Frequencia", x =name_var[1] ) +
    theme(legend.position = "bottom")
  
  a2<-data %>% 
    ggplot(aes(x=`Class Attribute`,y=data[,1])) +
    geom_boxplot() +
    labs(x = "Clase",  y = name_var[1])
  
  grid.arrange(a1,a2,nrow=1)
  cat('\n\n')
}

# Distribució variables quantitatives
hcc_factorT<-c(1:23,27:29)
figCapsFac <- names(hcc)[hcc_factorT]

# grafic 1 - diagrama barres
# grafic 2 - diagrama barres proporció
for (i in 1:length(hcc_factorT) ) {
  data <- hcc[,c(hcc_factorT[i],ind_CA )]
  name_var <- names(data)
  
  a1<-data %>% 
    ggplot(aes(x=data[,1],fill=`Class Attribute`))+ 
    geom_bar() +
    labs(fill="Clase", x = name_var[1], y = "Frequencia") + 
    theme(legend.position = "bottom",
          plot.caption = element_text(hjust = 0.5))
  
  a2<- data %>% 
    ggplot(aes(x=data[,1],fill=`Class Attribute`))+ 
    geom_bar(position = "fill") +
    labs(fill="Clase", x = name_var[1], y = "Proporcio") +
    theme(legend.position = "bottom",
          plot.caption = element_text(hjust = 0.5))
  
  grid.arrange(a1,a2, nrow= 1)
  cat('\n\n')
}

#Test Levene per homogeneïtat de la variància
tVar <- tibble()

for (i in var_normales){
  p_val = leveneTest(group=hcc$`Class Attribute`, y=hcc[,i])$`Pr(>F)`[1]
  tVar <- tVar %>% bind_rows(c("Variable" = i,"p_value" = p_val))
}

tVar <- tVar %>%mutate_at(.vars = c("p_value"), as.numeric)

#presentacio de resultats
kable(x = tVar, format = "latex", caption = "Homogeneitat de la variància a partir del test de Levene.",
      booktabs = TRUE, digits = 4, linesep = '') %>% 
  kable_styling(latex_options = c("HOLD_position")) %>% 
  row_spec(row=0, bold = TRUE)

#Test chi-quadrat per a variables qualitatives
testChi <- tibble()
for (i in hcc_factorT) {
  tau=table(hcc[,i], hcc$`Class Attribute`)
  chi=chisq.test(tau)
  testChi <- testChi %>% bind_rows(c(Variable=names(hcc[i]),p_value = chi$p.value))
}

#Filtrar les que rebutgen hipotesis nul·la
tau <- testChi %>% filter(p_value < 0.1) %>% 
  mutate_at(.vars = c("p_value"), as.numeric)

#Guardar vairables significatives
varSig <- tau

#Presentacio de resultats
kable(x = tau, format = "latex", caption = "Variables categòriques amb p<0.10 entre la classe",
      booktabs = TRUE, digits = 4, linesep = '') %>% 
  kable_styling(latex_options = c("HOLD_position")) %>% 
  row_spec(row=0, bold = TRUE)

#Test Mann-Whitney-Wilcoxon per a variables quantitatives
testWil <- tibble()
for (i in hcc_num) {
  wil=wilcox.test(hcc[hcc$`Class Attribute`==1,i],
                  hcc[hcc$`Class Attribute`==0,i],
                  mu = 0,paired = FALSE, conf.int = 0.95)
  testWil <- testWil %>%
    bind_rows(c(Variable=names(hcc[i]),p_value = wil$p.value))
}

#Filtrar les que rebutgen hipotesis nul·la
tau <- testWil %>% filter(p_value < 0.1) %>% 
  mutate_at(.vars = c("p_value"), as.numeric)

#Guardar vairables significatives
varSig <- varSig %>% bind_rows(tau)

#Resultats en dos columnes
tau <- cbind(tau[1:6,],tau[7:12,])

#Presentacio de resultats
kable(x = tau, format = "latex", caption = "Variables numèriques amb p<0.10 entre la classe",
      booktabs = TRUE, digits = 4, linesep = '') %>% 
  kable_styling(latex_options = c("HOLD_position", "scale_down")) %>% 
  row_spec(row=0, bold = TRUE)

# Correlació variables independents.
hcc_norm <- hcc
atr<-names(hcc) 

#normalitzacio variables quantitatives
hcc_norm[,hcc_num]<-scale(hcc_norm[,hcc_num])

#variables quantitatives tipus numeric
hcc2<-as.data.frame(lapply(hcc_norm,as.numeric))
names(hcc2)<-atr

#correlacio Spearman + filtrar alta correlació
hcc_cor<-cor(hcc2[,varSig$Variable], method = "spearman")
col.names <- colnames(hcc_cor)
hcc_corTop <- as.tibble(hcc_cor) %>% melt()
hcc_corTop$variable2 <- col.names
hcc_corTop <- hcc_corTop %>%
  dplyr::select(variable, variable2, value) %>%
  rename(corr = value) %>% 
  filter(variable != variable2) %>% 
  filter(abs(corr) >= 0.8) %>% 
  arrange(-corr) %>% 
  group_by(corr) %>% 
  mutate(Variables = paste(variable, collapse = '-')) %>% 
  dplyr::select(Variables, corr) %>% 
  unique()

#Presentació variables altament correlades
kable(x = hcc_corTop, format = "latex",
      caption = "Correlació entre variables, valors entre -1 i -08 o 0-8 i 1.",
      booktabs = TRUE, digits = 4, linesep = '') %>% 
  kable_styling(latex_options = c("HOLD_position")) %>% 
  row_spec(row=0, bold = TRUE)

#Grafica de totes les correlacions
corrplot(hcc_cor, cl.pos='n',tl.srt = 45,
         tl.cex = 1.25, type="upper",method = "circle",
         order="FPC", diag=FALSE)

# MELD calculat vs clase
hcc <- hcc %>%
  mutate(MELD = 3.78*`log_Total Bilirubin(mg/dL)`+
           11.2*log(`International Normalised Ratio`)+
           9.57*hcc$`log_Creatinine (mg/dL)`+
           6.43 )

#Boxplot per veure la distribució repsecte clase
hcc %>%
  ggplot(aes(x=`Class Attribute`,y=MELD)) +
  geom_boxplot() +
  labs(x = "Clase",  y = "MELD")

#Test Mann-Whitney-Wilcoxon
x = hcc %>% filter(`Class Attribute`==1) %>% pull(MELD)
y = hcc %>% filter(`Class Attribute`==0) %>% pull(MELD)
wilcox.test(x,y, mu = 0,paired = FALSE, conf.int = 0.95)

#Model logit
varElim <- c("Oxygen Saturation (%)","log_Total Bilirubin(mg/dL)",
             "International Normalised Ratio","log_Creatinine (mg/dL)")

#Variables seleccionades
var_selecc<- setdiff(varSig$Variable, varElim)

#Afegir les que falten
hcc_sel <- hcc[c(var_selecc, "MELD", "Class Attribute")] %>% 
  mutate_at(vars("Ascites degree","Encephalopathy degree"), as.numeric)

#Model logistic
modelo <- glm(`Class Attribute` ~ ., data = hcc_sel, family = "binomial")
summary(modelo)

#Prediccions de resultats + matriu de confusió
pred1 <- ifelse(test = modelo$fitted.values > 0.50, yes = 1, no = 0)
matConf <- table(hcc_sel$`Class Attribute`, pred1,
                 dnn = c("observacions", "prediccions"))

kable(x = matConf, format = "latex",
      caption = "Matriu de confusió model general, prediccions per columna, observacions per files",
      booktabs = TRUE, linesep = '') %>% 
  kable_styling(latex_options = c("HOLD_position")) %>% 
  row_spec(row=0, bold = TRUE)

#Grafic de matriu de confusió
mosaic(matConf, shade = T, colorize = T,
       gp = gpar(fill = matrix(c("green3", "red2", "red2", "green3"), 2, 2)))

#Èxit del model 
presPred <- 100*(matConf[1]+matConf[4])/nrow(hcc_sel)

#Predir els valors del model
prob <- predict(modelo,type = "response")
hcc$prob <- prob

#Pintar corba ROC, legacy.axes = TRUE per tal de pintar 1- Specifity, si no pinta Specifity.
plot.roc(`Class Attribute` ~ prob,data = hcc,legacy.axes = TRUE)

#Area sota la corba
auc(roc(`Class Attribute` ~ prob,data = hcc))

#Model logit mètode backward
modback <- stepAIC(modelo, trace=FALSE, direction="backward")
summary(modback)

#prediccio resultat + matriu de confusió
pred2 <- ifelse(test = modback$fitted.values > 0.50, yes = 1, no = 0)
matConf<- table(hcc_sel$`Class Attribute`, pred2,
                dnn = c("observacions", "prediccions"))

kable(x = matConf, format = "latex",
      caption = "Matriu de confusió model reduit, prediccions per columna, observacions per files",
      booktabs = TRUE, linesep = '') %>% 
  kable_styling(latex_options = c("HOLD_position")) %>% 
  row_spec(row=0, bold = TRUE)

#grafica matriu de confusió
mosaic(matConf, shade = T, colorize = T,
       gp = gpar(fill = matrix(c("green3", "red2", "red2", "green3"), 2, 2)))

#Èxit del model 
presPred2 <- 100*(matConf[1]+matConf[4])/nrow(hcc_sel)

#predicció variables
prob2 <- predict(modback,type = "response")
hcc$prob2 <- prob2

#Pintar corba ROC, legacy.axes = TRUE per tal de pintar 1- Specifity, si no pinta Specifity.
plot.roc(`Class Attribute` ~ prob2,data = hcc,legacy.axes = TRUE)

#Area sota la corba
auc(roc(`Class Attribute` ~ prob2,data = hcc))


tauContribucions <- tibble(Contribucions = c("Investigació previa", 
                                             "Redacció de les respostes", 
                                             "Desenvolupament del codi"),
                           Firma = c("A.D - D-L","A.D - D-L","A.D - D-L"))

kable(x = tauContribucions, format = "latex", caption = "Taula de contribucions",
      booktabs = TRUE, linesep = '') %>% 
  kable_styling(latex_options = c("HOLD_position")) %>% 
  row_spec(row=0, bold = TRUE)

