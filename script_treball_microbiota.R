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
# CREEM LA SBP
#--------------------------------------------------

#--------------------------------------------------
# PAS 1: Transformació CLR i Clustering Jeràrquic
#--------------------------------------------------
n_taxa <- ncol(microbiota_prop)
taxa_names <- colnames(microbiota_prop)

# Imputació de zeros (CZM)
microbiota_prop_nz <- tryCatch({
  cmultRepl(microbiota_prop, method = "CZM")
}, error = function(e){
  microbiota_prop
})

# Transformació CLR
microbiota_clr_cols <- clr(as.matrix(microbiota_prop_nz))

# Càlcul de distàncies i clustering
taxa_dist <- dist(t(microbiota_clr_cols), method = "euclidean")
hc_taxa <- hclust(taxa_dist, method = "ward.D2")

##################################################################
# GRÀFIC 1: DENDROGRAMA 
##################################################################
plot(hc_taxa, main = "GRÀFIC 1: Dendrograma de Tàxons (Basat en CLR)")


#--------------------------------------------------
# PAS 2: Construcció de la SBP i Coordenades OLR
#--------------------------------------------------

# Funció per extreure fulles (tàxons) d'un node
merges <- hc_taxa$merge
get_leaves <- function(node_index) {
  if (node_index < 0) return(-node_index)
  left <- merges[node_index, 1]
  right <- merges[node_index, 2]
  c(get_leaves(left), get_leaves(right))
}

# Creació de la Matriu SBP (Sequential Binary Partition)
SBP_hc <- matrix(0, nrow = n_taxa, ncol = n_taxa - 1,
                 dimnames = list(taxa_names, paste0("hc_b", 1:(n_taxa - 1))))

for (i in 1:(n_taxa - 1)) {
  left_node <- merges[i, 1]
  right_node <- merges[i, 2]
  left_leaves <- taxa_names[get_leaves(left_node)]
  right_leaves <- taxa_names[get_leaves(right_node)]
  
  SBP_hc[left_leaves, i] <- 1
  SBP_hc[right_leaves, i] <- -1
}

# Càlcul de les coordenades OLR (Orthogonal Log-Ratio)
basis_hc <- sbp_basis(SBP_hc)
microbiota_olr_hc <- coordinates(as.matrix(microbiota_prop_nz), basis = basis_hc)
microbiota_olr_hc_df <- as.data.frame(microbiota_olr_hc)
colnames(microbiota_olr_hc_df) <- paste0("hc_b", 1:ncol(microbiota_olr_hc_df))
microbiota_olr_hc_df$obesity <- data$obesity


#--------------------------------------------------
# PAS 3: Test Estadístic de TOTS els balanços
#--------------------------------------------------
# Fem un t-test per cada balanç (cada columna OLR) vs obesitat
results_hc_all <- data.frame(
  Balance = colnames(microbiota_olr_hc_df)[-ncol(microbiota_olr_hc_df)],
  p_value = sapply(colnames(microbiota_olr_hc_df)[-ncol(microbiota_olr_hc_df)], function(b){
    t.test(microbiota_olr_hc_df[[b]] ~ microbiota_olr_hc_df$obesity)$p.value
  })
)
results_hc_all$SBP <- "Jeràrquica"


##################################################################
# GRÀFIC 2: RÀNQUING DE BALANÇOS 

##################################################################
results_hc_all$log10_p <- -log10(results_hc_all$p_value)
results_hc_all$Significant <- ifelse(results_hc_all$p_value < 0.05, "Significatiu (p<0.05)", "No Significatiu")

p_ranking <- ggplot(results_hc_all, aes(x = reorder(Balance, log10_p), y = log10_p, color = Significant)) +
  geom_segment(aes(xend = reorder(Balance, log10_p), yend = 0), linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  coord_flip() +
  scale_color_manual(values = c("Significatiu (p<0.05)" = "#E63946", "No Significatiu" = "grey")) +
  labs(title = "GRÀFIC 2: Significació dels Balanços Estructurals ",
       x = "Balanç (Ordenat per p-value)",
       y = expression(-log[10](p-value))) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_ranking)


#--------------------------------------------------
# PAS 4: Selecció i Traducció del Balanç Guanyador
#--------------------------------------------------
# Basant-nos en el Gràfic 2, seleccionem el balanç més significatiu.
# (En aquest exemple, suposem que és 'hc_b10')
balance_interes_nom <- "hc_b10"

# Extraiem els tàxons que componen aquest balanç
balanc_hc10_vector <- SBP_hc[, balance_interes_nom]
taxa_grup_positiu <- names(balanc_hc10_vector[balanc_hc10_vector == 1])
taxa_grup_negatiu <- names(balanc_hc10_vector[balanc_hc10_vector == -1])

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

##################################################################
# GRÀFIC 3: BOXPLOT DEL BALANÇ CLAU 
##################################################################
p_boxplot_clau <- ggplot(microbiota_olr_hc_df, 
                         aes(x = obesity, y = .data[[balance_interes_nom]], fill = obesity)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
  labs(title = paste("GRÀFIC 3: Resum Estadístic (Balanç", balance_interes_nom, ")"),
       subtitle = paste("p-value =", round(results_hc_all[balance_interes_nom, "p_value"], 5)),
       x = "Estat d'Obesitat",
       y = paste("Coordenada OLR (", balance_interes_nom, ")")) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("No obesitat" = "#457B9D", "Obesitat" = "#E63946"))
print(p_boxplot_clau)


##################################################################
# GRÀFIC 4: GRÀFIC REPRESENTATIU 
##################################################################
gm_mean <- function(x) {
  if (any(x <= 0)) return(NA)
  exp(mean(log(x)))
}
gm_positiu <- apply(microbiota_prop_nz[, taxa_grup_positiu, drop = FALSE], 1, gm_mean)
gm_negatiu <- apply(microbiota_prop_nz[, taxa_grup_negatiu, drop = FALSE], 1, gm_mean)

df_plot_representatiu <- data.frame(
  log_gm_positiu = log(gm_positiu),
  log_gm_negatiu = log(gm_negatiu),
  obesity = data$obesity
)

p_representatiu <- ggplot(df_plot_representatiu, 
                          aes(x = log_gm_negatiu, y = log_gm_positiu, 
                              color = obesity, shape = obesity)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = paste("GRÀFIC 4: El Balanç", balance_interes_nom),
       x = "Log( Mitjana Geomètrica Grup -1 [Denominador] )",
       y = "Log( Mitjana Geomètrica Grup +1 [Numerador] )") +
  scale_color_manual(values = c("No obesitat" = "#457B9D", "Obesitat" = "#E63946")) +
  scale_shape_manual(values = c("No obesitat" = 16, "Obesitat" = 17)) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_representatiu)

##################################################################
# GRÀFIC 5  BOXPLOTS DE SUPORT 
##################################################################
taxa_interes <- c(taxa_grup_positiu, taxa_grup_negatiu)
data_suport <- cbind(microbiota_prop_nz[, taxa_interes, drop = FALSE], obesity = data$obesity)

data_suport_long <- data_suport %>%
  pivot_longer(cols = -obesity, names_to = "Taxon", values_to = "Abundancia") %>%
  mutate(Grup_Balanç = ifelse(Taxon %in% taxa_grup_positiu, 
                              "Grup +1 (Numerador)", 
                              "Grup -1 (Denominador)"))

p_suport <- ggplot(data_suport_long, aes(x = obesity, y = Abundancia, fill = obesity)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, height = 0) +
  facet_wrap(~ Grup_Balanç + Taxon, scales = "free_y", ncol = 4) +
  labs(title = "GRÀFIC 6 : Evidència (Abundància Tàxons Individuals)",
       x = "Estat d'Obesitat",
       y = "Abundància Relativa ") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("No obesitat" = "#457B9D", "Obesitat" = "#E63946"))
print(p_suport)
