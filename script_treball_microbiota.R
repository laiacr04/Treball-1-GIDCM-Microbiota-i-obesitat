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

# 1) Definició SBP seqüencial equilibrada
n_taxa <- ncol(microbiota_prop)
taxa_names <- colnames(microbiota_prop)
SBP_seq <- matrix(0, nrow = n_taxa, ncol = n_taxa - 1,
                  dimnames = list(taxa_names, paste0("seq_b", 1:(n_taxa-1))))
for (i in 1:(n_taxa - 1)) {
  SBP_seq[1:i, i] <- 1
  SBP_seq[(i+1):n_taxa, i] <- -1
}

# Mostrem les  columnes de la SBP per inspecció
print(SBP_seq[, 1:min(12, ncol(SBP_seq))])

# 2) Base i coordenades OLR (sobre les proporcions ja calculades)
basis_seq <- tryCatch(sbp_basis(SBP_seq), warning = function(w) { message("Avis sbp_basis: ", w$message); sbp_basis(SBP_seq) })
microbiota_olr_seq <- coordinates(as.matrix(microbiota_prop), basis = basis_seq)

# 3) Data.frame 
microbiota_olr_seq_df <- as.data.frame(microbiota_olr_seq)
colnames(microbiota_olr_seq_df) <- paste0("seq_b", 1:ncol(microbiota_olr_seq_df))
microbiota_olr_seq_df$obesity <- data$obesity

# 4) Resum
results_seq_all <- data.frame(
  Balance = colnames(microbiota_olr_seq_df)[-ncol(microbiota_olr_seq_df)],
  p_value = sapply(colnames(microbiota_olr_seq_df)[-ncol(microbiota_olr_seq_df)], function(b){
    t.test(microbiota_olr_seq_df[[b]] ~ microbiota_olr_seq_df$obesity)$p.value
  })
)
results_seq_all
results_seq_all$SBP <- "Seqüencial"


# 5) Gràfic boxplot 
p_seq_b1 <- ggplot(microbiota_olr_seq_df, aes(x = obesity, y = seq_b1, fill = obesity)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "SBP Seqüencial — Balanç 1 vs Obesitat", y = "OLR (seq_b1)", x = "Obesitat") +
  theme_minimal() +
  theme(legend.position = "none")

print(p_seq_b1)
#--------------------------------------------------
# SBP jerarquica 
#--------------------------------------------------

# Transformació CLR i clustering
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

#--------------------------------------------------
# SBP hibrida
#--------------------------------------------------


# --- 1) Transformació CLR per evitar zeros ---
microbiota_prop_nz <- tryCatch({
  cmultRepl(microbiota_prop, method = "CZM")
}, error = function(e) microbiota_prop)

microbiota_clr <- clr(as.matrix(microbiota_prop_nz))

# --- 2) Diferències mitjanes per grup ---
group <- data$obesity
mean_diff <- apply(microbiota_clr, 2, function(x) mean(x[group=="Obesitat"]) - mean(x[group=="No obesitat"]))

# --- 3) Selecció tàxons per balanç +1 i -1 ---
# Exemple: escollim els 5 tàxons més positius i 5 més negatius
top_positive <- names(sort(mean_diff, decreasing = TRUE)[1:5])
top_negative <- names(sort(mean_diff, decreasing = FALSE)[1:5])

taxa_names <- colnames(microbiota_prop_nz)
SBP_opt <- matrix(0, nrow = length(taxa_names), ncol = 1, dimnames = list(taxa_names, "Opt_b1"))

SBP_opt[top_positive, 1] <- 1
SBP_opt[top_negative, 1] <- -1

# --- 4) Base OLR i coordenades ---
basis_opt <- sbp_basis(SBP_opt)
microbiota_olr_opt <- coordinates(as.matrix(microbiota_prop_nz), basis = basis_opt)
microbiota_olr_opt_df <- data.frame(OLR = microbiota_olr_opt[,1], obesity = group)

# --- 5) Resum i T-test ---
cat("\nTàxons +1:\n"); print(top_positive)
cat("\nTàxons -1:\n"); print(top_negative)

tt_opt <- t.test(OLR ~ obesity, data = microbiota_olr_opt_df)
print(tt_opt)

# --- 6) Boxplot per visualització ---
ggplot(microbiota_olr_opt_df, aes(x = obesity, y = OLR, fill = obesity)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Balanç híbrid optimitzat — separació per obesitat",
       x = "Obesitat", y = "Coordenada OLR") +
  theme_minimal() +
  theme(legend.position = "none")
opt_hybrid_results <- data.frame(
  Balance = colnames(SBP_opt),  # "Opt_b1"
  p_value = t.test(microbiota_olr_opt_df$OLR ~ microbiota_olr_opt_df$obesity)$p.value,
  SBP = "Híbrid optimitzat"
)


# -----------------------------
# 1) Combinar tots els resultats
# -----------------------------
all_results <- rbind(results_seq_all, results_hc_all, opt_hybrid_results)
all_results <- all_results[order(all_results$p_value), ]

cat("\n===== Resum comparatiu de tots els balanços =====\n")
print(head(all_results, 10))
all_results
# -----------------------------
# 2) Gràfic: distribució p-values per SBP
# -----------------------------
ggplot(all_results, aes(x = SBP, y = -log10(p_value), fill = SBP)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Comparació de la força estadística dels balanços",
       x = "Tipus de SBP", y = expression(-log[10](p))) +
  scale_fill_manual(values = c("#5DADE2", "#58D68D", "#F5B041"))

# -----------------------------
# 3) Millor balanç global
# -----------------------------
best_balance <- all_results[which.min(all_results$p_value), ]
cat("\n===== Millor balanç global =====\n")
print(best_balance)
cat("\nLa millor SBP segons el contrast d'obesitat és:",
    best_balance$SBP, "amb el balanç", best_balance$Balance,
    " (p-value =", round(best_balance$p_value, 5), ")\n")

# -----------------------------
# 4) Balanços significatius (p < 0.05)
# -----------------------------
sig_balances <- subset(all_results, p_value < 0.05)
if (nrow(sig_balances) > 0) {
  ggplot(sig_balances, aes(x = Balance, y = -log10(p_value), fill = SBP)) +
    geom_col(position = "dodge", alpha = 0.8) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Balanços significatius (p < 0.05)",
         x = "Balanç", y = expression(-log[10](p))) +
    scale_fill_manual(values = c("#5DADE2", "#58D68D", "#F5B041"))
} else {
  cat("\nCap balanç significatiu (p < 0.05) trobat.\n")
}

# -----------------------------
# 5) Resum estadístic per tipus de SBP
# -----------------------------
summary_SBP <- aggregate(p_value ~ SBP, data = all_results, 
                         FUN = function(x) c(min = min(x), mean = mean(x)))
summary_SBP <- do.call(data.frame, summary_SBP)

cat("\n===== Estadístics resum per tipus de SBP =====\n")
print(summary_SBP)

ggplot(summary_SBP, aes(x = SBP, y = p_value.min, fill = SBP)) +
  geom_col(alpha = 0.8) +
  theme_minimal() +
  labs(title = "Comparació del p-value mínim per tipus de SBP",
       x = "Tipus de SBP", y = "p-value mínim (menor = millor)") +
  scale_fill_manual(values = c("#5DADE2", "#58D68D", "#F5B041"))

# -----------------------------
# 6) Detalls del millor balanç
# -----------------------------
b_name <- best_balance$Balance
balance_vec <- SBP_hc[, b_name]

taxa_positive <- names(balance_vec[balance_vec == 1])
taxa_negative <- names(balance_vec[balance_vec == -1])

cat("\n========== Informació del balanç", b_name, "==========\n")
cat("\nTàxons (+1):\n"); print(taxa_positive)
cat("\nTàxons (−1):\n"); print(taxa_negative)
cat("\nNombre de tàxons per grup:\n")
cat("+1:", length(taxa_positive), " |  −1:", length(taxa_negative), "\n")

# -----------------------------
# 7) Boxplot del millor balanç
# -----------------------------
ggplot(microbiota_olr_hc_df, aes(x = obesity, y = .data[[b_name]], fill = obesity)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = paste("Balanç", b_name, "— separació per obesitat"),
       x = "Obesitat", y = "Coordenada OLR") +
  theme_minimal() +
  theme(legend.position = "none")



