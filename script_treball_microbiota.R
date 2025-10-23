#############################################
##### TREBALL 1 - MICROBIOTA I OBESITAT #####
#############################################

#--------------------------------------------------
# Carregar paquets necessaris
#--------------------------------------------------

#Instal·lar 
#install.packages("ggplot2") 
#install.packages("dplyr")
#install.packages("compositions")
#install.packages("coda.base")
#install.packages("tidyr")

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

#--------------------------------------------------
# Busquem quina és la millor SBP
#--------------------------------------------------
#--------------------------------------------------
# SBP jerarquica 
#--------------------------------------------------
# Transformació CLR i clustering
n_taxa <- ncol(microbiota_prop)
taxa_names <- colnames(microbiota_prop)
microbiota_prop_nz <- tryCatch({
  cmultRepl(microbiota_prop, method = "CZM")
}, error = function(e){
  microbiota_prop
})
microbiota_clr_cols <- clr(as.matrix(microbiota_prop_nz))
taxa_dist <- dist(t(microbiota_clr_cols), method = "euclidean")
hc_taxa <- hclust(taxa_dist, method = "ward.D2")

# Dibuixar dendrograma (per inspecció)
plot(hc_taxa, main = "Dendrograma de tàxons (basat en CLR)")

# Construcció SBP jeràrquic
merges <- hc_taxa$merge
get_leaves <- function(node_index) {
  if (node_index < 0) return(-node_index)
  left <- merges[node_index, 1]
  right <- merges[node_index, 2]
  c(get_leaves(left), get_leaves(right))
}

SBP_hc <- matrix(0, nrow = n_taxa, ncol = n_taxa - 1,
                 dimnames = list(taxa_names, paste0("hc_b", 1:(n_taxa - 1))))

for (i in 1:(n_taxa - 1)) {
  left_node <- merges[i, 1]
  right_node <- merges[i, 2]
  
  # Obtenim les fulles associades de manera segura
  left_leaves <- taxa_names[get_leaves(left_node)]
  right_leaves <- taxa_names[get_leaves(right_node)]
  
  SBP_hc[left_leaves, i] <- 1
  SBP_hc[right_leaves, i] <- -1
}

# Coordenades OLR jeràrquiques
basis_hc <- sbp_basis(SBP_hc)
microbiota_olr_hc <- coordinates(as.matrix(microbiota_prop_nz), basis = basis_hc)
microbiota_olr_hc_df <- as.data.frame(microbiota_olr_hc)
colnames(microbiota_olr_hc_df) <- paste0("hc_b", 1:ncol(microbiota_olr_hc_df))
microbiota_olr_hc_df$obesity <- data$obesity

# T-test pel primer balanç jeràrquic
tt_hc_b1 <- t.test(hc_b1 ~ obesity, data = microbiota_olr_hc_df)

# Boxplot
p_hc <- ggplot(microbiota_olr_hc_df, aes(x = obesity, y = hc_b1, fill = obesity)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "SBP Jeràrquica — Balanç 1 vs Obesitat", y = "OLR (hc_b1)", x = "Obesitat") +
  theme_minimal() +
  theme(legend.position = "none")

print(p_hc)

# Resultats resum
results_hc_all <- data.frame(
  Balance = colnames(microbiota_olr_hc_df)[-ncol(microbiota_olr_hc_df)],
  p_value = sapply(colnames(microbiota_olr_hc_df)[-ncol(microbiota_olr_hc_df)], function(b){
    t.test(microbiota_olr_hc_df[[b]] ~ microbiota_olr_hc_df$obesity)$p.value
  })
)
results_hc_all
results_hc_all$SBP <- "Jeràrquica"






# ------------------------------------------------------------------
# PAS 1: Definir el balanç d'interès
# ------------------------------------------------------------------
# Pels nostres resultats, el balanç més significatiu de la SBP
# jeràrquica (no supervisada) va ser 'hc_b10'.
balance_interes_nom <- "hc_b10"

# ------------------------------------------------------------------
# PAS 2: Assegurar que la matriu SBP existeix
# ------------------------------------------------------------------
# Aquest objecte 'SBP_hc' s'hauria d'haver creat en el teu script anterior
if (!exists("SBP_hc")) {
  stop("Error: L'objecte 'SBP_hc' no s'ha trobat. 
       Assegura't d'haver executat el codi de clustering jeràrquic primer.")
}

# ------------------------------------------------------------------
# PAS 3: Extreure els tàxons d'aquest balanç
# ------------------------------------------------------------------
# Extraiem la columna (el balanç) que ens interessa
balanc_hc10_vector <- SBP_hc[, balance_interes_nom]

# Mirem quins tàxons estan a cada grup
taxa_grup_positiu <- names(balanc_hc10_vector[balanc_hc10_vector == 1])
taxa_grup_negatiu <- names(balanc_hc10_vector[balanc_hc10_vector == -1])
taxa_grup_zero <- names(balanc_hc10_vector[balanc_hc10_vector == 0])

# ------------------------------------------------------------------
# PAS 4: Mostrar els resultats
# ------------------------------------------------------------------
cat("\n======================================================\n")
cat("TRADUCCIÓ TÈCNICA DEL BALANÇ:", balance_interes_nom, "\n")
cat("======================================================\n")
cat("Aquest balanç (p-value =", 
    round(results_hc_all[balance_interes_nom, "p_value"], 5), 
    ")\nrepresenta el log-ràtio entre dos grups:\n\n")

cat("--- GRUP +1 (Numerador) ---\n")
if (length(taxa_grup_positiu) > 0) {
  print(taxa_grup_positiu)
} else {
  cat("[Cap tàxon en aquest grup]\n")
}

cat("\n--- GRUP -1 (Denominador) ---\n")
if (length(taxa_grup_negatiu) > 0) {
  print(taxa_grup_negatiu)
} else {
  cat("[Cap tàxon en aquest grup]\n")
}

cat("\n--- Tàxons no inclosos en aquest balanç ---\n")
if (length(taxa_grup_zero) > 0) {
  print(taxa_grup_zero)
} else {
  cat("[Tots els tàxons estan inclosos en el balanç]\n")
}
cat("\n======================================================\n")

