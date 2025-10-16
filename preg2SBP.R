#############################################
##### TREBALL 1 - MICROBIOTA I OBESITAT #####
#############################################

#--------------------------------------------------
# Carregar paquets necessaris
#--------------------------------------------------

#Instal·lar 
install.packages("ggplot2") 
install.packages("dplyr")
install.packages("compositions")
install.packages("coda.base")
install.packages("tidyr")

#Carregar
library(ggplot2)
library(dplyr)
library(compositions)
library(coda.base)
library(tidyr)


#--------------------------------------------------
# Carregar l'arxiu de dades
#--------------------------------------------------

load("MICROBIOMA6.Rdata")
ls() #objecte carregat, anomenat 'data'

#--------------------------------------------------
# Neteja i preparació de les dades
#--------------------------------------------------
str(data) #per veure els tipus de dades i primers valors
head(data) #per visualitzar les primeres files
dim(data) #files i columnes 

# 1. Comprovació de valors nulls

colSums(is.na(data))
#Les variables fat_to i waist_ci tenen 5 i 3 valors N/A respectivament, decidim imputar-los, és a dir, substituir el seu valor per la mitjana
#de la resta de valors, per no perdre cap pacient i no alterar molt la distribució.
data$waist_ci[is.na(data$waist_ci)] <- mean(data$waist_ci, na.rm = TRUE)
data$fat_to[is.na(data$fat_to)] <- mean(data$fat_to, na.rm = TRUE)

colSums(is.na(data))#tornem a comptar els valors nulls, comprovem que ara n'hi ha 0 en tot el nostre set de dades

# 2. Comprovació de duplicats pel codi del pacient
sum(duplicated(data$codeIM))  
#Comprovat, no hi ha cap duplicat.

#--------------------------------------------------
# Càlcul i comprovació de les variables derivades
#--------------------------------------------------
#Verifiquem que la variable BMI està ben calculada
data$BMI_check <- data$weight / (data$height/100)^2
summary(data$bmi - data$BMI_check)

#Comprovem la coherència entre "bmi" i "obesity"
table(data$obesity, data$bmi >= 30)
#La classificació d'obesitat és coherent amb el criteri BMI>=30, per tant, no cal recalcular cap variable derivada

#----------------------------------------------------------------
# PREGUNTA 1: Descripció de les variables "obesity" i "groupIM"
#----------------------------------------------------------------

# ----- Variable obesity -----
#Taula de freqüències absolutes
table(data$obesity)

#Taula de freqüències relatives (percentatges)
prop.table(table(data$obesity))*100


#Afegim noms més clars als nivells per facilitar la lectura
data$obesity <- factor(data$obesity, 
                       levels = c("0", "1"),
                       labels = c("No obesitat", "Obesitat"))

#Resum descriptiu
summary(data$obesity)

#Gràfic de barres
ggplot(data, aes(x = obesity, fill = obesity)) +
  geom_bar() +
  labs(title = "Distribució de la variable Obesitat",
       x = "Categoria d'obesitat", 
       y = "Nombre de pacients") +
  theme_minimal() +
  theme(legend.position = "none")

pie(table(data$obesity))

#Interpretació: El gràfic i les taules mostren quants pacients del total (n=120) presenten obesitat
#segons el criteri aplicat (BIM>=30). Observem que hi ha 70 obesos (58.3%) i 50 no obesos (41.6%).

# ----- Variable groupIM -----

#Reetiquetem els nivells per claredat
data$groupIM <- factor(data$groupIM,
                       levels = c("1.0","2.0","3.0","4.0","5.0","6.0"),
                       labels = c("Women pre-menopause with obesity",
                                  "Women post-menopause with obesity",
                                  "Men with obesity",
                                  "Women pre-menopause without obesity",
                                  "Women post-menopause without obesity",
                                  "Men without obesity"))

#Taula de freqüències absolutes
table(data$groupIM)

#Taula de percentatges
prop.table(table(data$groupIM))*100

#Gràfic de barres
ggplot(data, aes(x = groupIM, fill = groupIM)) +
  geom_bar() +
  labs(title = "Distribució de la variable groupIM",
       x = "Grup de pacients", 
       y = "Nombre de pacients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")

#Diagrama de sectors
pie(table(data$groupIM))



#--------------------------------------------------------------------------------
# PREGUNTA 2: Anàlisi descriptiva multivariant de la composició de la microbiota
#--------------------------------------------------------------------------------
#Seleccionem les variables de microbiota
microbiota <- data[, 11:22]

#Comprovem dimensions i resum
dim(microbiota)
head(microbiota)
summary(microbiota)

microbiota_long <- microbiota %>%
  tidyr::pivot_longer(cols = everything(), names_to = "Taxon", values_to = "Abundància")

ggplot(microbiota_long, aes(x = Taxon, y = Abundància)) +
  geom_boxplot(fill = "skyblue", alpha = 0.7) +
  labs(title = "Distribució de les abundàncies relatives dels components de la microbiota",
       x = "Taxó",
       y = "Abundància relativa") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Comprovar si cada fila suma 1 i normalitzar
row_sums <- rowSums(microbiota)
microbiota_prop <- microbiota/row_sums
summary(rowSums(microbiota_prop)) #comprovat: totes les files sumen 1

#Còpia de seguretat
microbiota_pc <- microbiota_prop

#Transformaió CLR (per poder fer anàlisi composicional)
microbiota_clr <- clr(as.matrix(microbiota_pc))
microbiota_clr_df <- as.data.frame(microbiota_clr)

#PCA sobre coordenades CLR
pca_clr <- prcomp(microbiota_clr_df, center=TRUE, scale.=FALSE)
summary(pca_clr)

#Visualització multivariant de la composició
scores <- as.data.frame(pca_clr$x)
scores$obesity <- data$obesity
ggplot(scores, aes(x = PC1, y = PC2, color = obesity)) +
  geom_point(size = 2) +
  labs(title = "PCA sobre coordenades CLR (microbiota)",
       x = paste0("PC1 (", round(100*summary(pca_clr)$importance[2,1],1), "%)"),
       y = paste0("PC2 (", round(100*summary(pca_clr)$importance[2,2],1), "%)")) +
  theme_minimal()


#Interpretació de la PCA (ja calculada anteriorment)
summary(pca_clr)

# Mostrem les càrregues (loadings) dels components
loadings <- as.data.frame(pca_clr$rotation)
head(loadings)

# Taxons més importants en els dos primers components
loadings_PC1 <- sort(abs(loadings[,1]), decreasing = TRUE)
loadings_PC2 <- sort(abs(loadings[,2]), decreasing = TRUE)

cat("Taxons més importants en PC1:\n")
print(names(loadings_PC1[1:5]))
cat("Taxons més importants en PC2:\n")
print(names(loadings_PC2[1:5]))

# Gràfic de les càrregues del primer component
ggplot(loadings, aes(x = reorder(rownames(loadings), PC1), y = PC1)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Contribució dels taxons al primer component principal (PC1)",
       x = "Taxó", y = "Càrrega (loading)")

#Definició d'una SBP (Sequential Binary Partition)
# Creem una SBP equilibrada (divideix successivament els 12 components)
# Això és una decisió neutral quan no hi ha criteri biològic concret
n_taxa <- ncol(microbiota_prop)
SBP <- matrix(0, nrow = n_taxa, ncol = n_taxa - 1)
for (i in 1:(n_taxa - 1)) {
  SBP[1:i, i] <- 1
  SBP[(i+1):n_taxa, i] <- -1
}

# Comprovem la SBP
SBP

library(coda.base)

microbiota_olr <- coordinates(
  as.matrix(microbiota_prop),
  basis = sbp_basis(SBP)
)

microbiota_olr_df <- as.data.frame(microbiota_olr)
head(microbiota_olr_df)


