

####- 1-- Packages -------------------------------####
#-##################################################-#
pacman::p_load(
  GEOquery, tidyverse, ggrepel, limma, oligo, DT, pheatmap, tidyplots,  affy, oligoClasses, testit,
  hgu133plus2.db, org.Hs.eg.db, reshape2, survminer, arrayQualityMetrics, testit, biomaRt, AgiMicroRna, 
  Biobase, marray, pd.bovgene.1.0.st, readxl, VennDiagram, RColorBrewer,ggVennDiagram, ggplot2, fgsea,msigdbr)

#-##################################################-#
####- 2-- Set up your working directory ----------####
#-##################################################-#

setwd("C:/Users/ATL2024_03/OneDrive - Kagoshima University (1)/Maestría/BLV PROYECT/New attempt")

#setwd("C:/Users/user/OneDrive - Kagoshima University (1)/Maestría/BLV PROYECT/New attempt")




#-##################################################-#
####- 3-- Setting the Array list -----------------####
#-##################################################-#

raw <- list.files(pattern = "\\.CEL$", ignore.case = TRUE)

#-##################################################-#
####- 4-- Reading the list -----------------------####
#-##################################################-#

raw_data <- read.celfiles(raw)

#-##################################################-#
####- 5-- Verifyng ChipType and Information ------####
#-##################################################-#

#annotation(raw_data)
#featureNames(raw_data)
#probeNames(raw_data)
#class(raw_data)
#featureData(raw_data)
#pData(raw_data)
#protocolData(raw_data)



#-##################################################-#
####- 6-- QC Before Normalization ----------------####
#-##################################################-#

## If The array do not support NUSE and RLE in the arrayQualityMetrics QC, we must add them

arrayQualityMetrics(expressionset = raw_data,
                    outdir = "QC_before",
                    force= T,
                    do.logtransform = T) # Due to no normalization


#-##################################################-#
####- 7-- Normalization by RMA -------------------####
#-##################################################-#

normalized_data <- oligo::rma(raw_data, background=T,normalize=T)

#-##################################################-#
####- 8-- QC After Normalization -----------------####
#-##################################################-#

## If The array do not support NUSE and RLE in the arrayQualityMetrics QC, we must add them

arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "QC_after",
                    force= T) # No need logtransform, already normalized


#-##################################################-#
####- 9-- Filter Low Variance Transcripts---------####
#-##################################################-#
processed_data <- as.data.frame(exprs(normalized_data))

# exprs es tu matriz genes x muestras
var_genes <- apply(processed_data, 1, var) # varianza por gen
threshold <- quantile(var_genes, 0.10)     # percentil 10%
reprocessed_data <- processed_data[var_genes > threshold, ]

#-##################################################-#
####- 10-- Saving Expression Data -----------------####
#-##################################################-#

write.table(reprocessed_data,file="expression_probes_data.txt", quote=FALSE, sep="\t")


#-###################################################-#
####- 11-- Obtain Annotation Data -----------------####
#-###################################################-#

annotation_data <- getGEO("GPL16500")
Table(annotation_data)
annot <- Table(annotation_data)[,c("ID", "probeset_id","total_probes", "gene_assignment")]





#-###################################################-#
####- 12-- Merging (Annotation + Expression) Data -####
#-###################################################-#

reprocessed_data <- reprocessed_data|>
  as.data.frame ()|>
  tibble::rownames_to_column(var ="ID")


expression_annotation <- merge(annot, reprocessed_data, by ="ID")
str(expression_annotation)


#-###################################################-#
####- 13-- Modifying Gene Names   -----------------####
#-###################################################-#


expression_annotation <- expression_annotation %>%
  as.data.frame() %>%
  mutate(
    gene_assignment = case_when(
      is.na(gene_assignment) ~ NA_character_,                        # si es NA, queda NA
      str_detect(gene_assignment, " // ") ~ str_split(gene_assignment, " // ", simplify = TRUE)[,2], # si tiene slashes, toma el segundo
      TRUE ~ gene_assignment                                          # si no tiene slashes, deja igual
    ),
    # Eliminar si está vacío o no tiene letras
    gene_assignment = ifelse(
      is.na(gene_assignment) | gene_assignment == "" | !grepl("[A-Za-z]", gene_assignment),
      NA,
      gene_assignment
    )
  )%>%
  filter(!is.na(gene_assignment)) 




#-############################################################-#
####- 14-- Handle Multiple Probes Mapping to a Single Gene -####
#-############################################################-#

# Several strategies exist

#1. Retain probe with highest expression or variance
#2. Average or summarize  probe signals    ### selected ###
#4. Remove duplicate probes to maintain one row per gene

#Collapsing all based on gene_assignment

# Converting all to matrix to use it with avereps 
expression_matrix <- as.matrix(expression_annotation[,-(1:4)])  # Quitar columnas de anotación
exprs_averaged <- limma::avereps(expression_matrix, ID = expression_annotation$gene_assignment)

# Ver resultado
dim(expression_annotation)
dim(exprs_averaged)


#-###################################################-#
####- 15-- Saving (Annotation + Expression) Data --####
#-###################################################-#


write.table(exprs_averaged,file="expression_annotation_data.txt", quote=FALSE, sep="\t")

#-###################################################-#
####- 15-- DEG ------------------------------------####
#-###################################################-#

# Renaming the samples
colnames(exprs_averaged) <- sub("\\(.*", "", colnames(exprs_averaged))



# Opening sample data
BLV_sample <- read_excel("BLV_samples.xlsx")


# Create design matrix based on final_subtype

design <- model.matrix(~ 0 + factor(BLV_sample$Status))
design
rownames(design) <- BLV_sample$Patients
colnames(design) <- levels(factor(BLV_sample$Status))
head(design, 20)


# Fit linear model

fit <- lmFit(exprs_averaged, design)
# Make contrasts (adjust based on your comparisons)
levels(factor(BLV_sample$Status))

contrast.matrix <- makeContrasts(
  EBL_AC_LP = EBL - AC_LP,
  EBL_HD = EBL - HD,
  AC_LP_HD = AC_LP- HD,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
fit2

# Get significant probes (FDR < 0.01)
top_probes <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "B", p.value = 0.05) #mayores a 1 valor absoluto
write.table(top_probes,file="BLV_EBL_LN_DEG.txt", quote=FALSE, sep="\t")


summary(decideTests(fit2, adjust.method="BH", p.value=0.05, lfc=1))


#-###################################################-#
####- 16-- Venn Diagram ---------------------------####
#-###################################################-#

#top_probes<-read.table("TCGA_dichotomized_STAD_MSI.txt",
              # sep="\t", check.names = FALSE,header=TRUE)

# Getting individual coeficient  EBL vs. AC_LP / EBL vs. HD / AC vs. LP_HD

EBL_AC_LP <- topTable(fit2, coef = "EBL_AC_LP", number = Inf,
                      adjust.method = "BH", p.value = 0.05, lfc = 1, sort.by="B")  #1467
EBL_HD    <- topTable(fit2, coef = "EBL_HD", number = Inf,
                      adjust.method = "BH", p.value = 0.05, lfc = 1, sort.by="B")  #1891
AC_LP_HD  <- topTable(fit2, coef = "AC_LP_HD", number = Inf,
                      adjust.method = "BH", p.value = 0.05, lfc = 1, sort.by="B")  #0


# Getting gene list for Venn Diagram plotting 
genes_EBL_AC_LP <- rownames(EBL_AC_LP)
genes_EBL_HD    <- rownames(EBL_HD)
genes_AC_LP_HD   <- rownames(AC_LP_HD)

#### Specific EBL genes
#### - Between EBL_vs_AC_LP and EBL_vs_HD
#### - No in AC_LP_vs_HD (Exclude initial infection genes)

# List of items
x <- list(genes_EBL_AC_LP, genes_EBL_HD,genes_AC_LP_HD)

# Plotting Venn diagram with custom category names
ggVennDiagram(x,category.names = c("         EBL vs. nonEBL",
                                  "EBL vs. HD",
                                  "NonEBL vs HD"))+ 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")


# Plotting Upregulated and Downregulated genes 


# UP
EBL_up_AC_LP<-EBL_AC_LP[EBL_AC_LP$adj.P.Val < 0.05, ]
EBL_up_AC_LP<-EBL_AC_LP[EBL_AC_LP$logFC > 1, ]
# Down
EBL_down_AC_LP<-EBL_AC_LP[EBL_AC_LP$adj.P.Val < 0.05, ]
EBL_down_AC_LP<-EBL_AC_LP[EBL_AC_LP$logFC < 1, ]




# UP
EBL_up_HD<-EBL_HD[EBL_HD$adj.P.Val < 0.05, ]
EBL_up_HD<-EBL_HD[EBL_HD$logFC > 1, ]
# Down
EBL_down_HD<-EBL_HD[EBL_HD$adj.P.Val < 0.05, ]
EBL_down_HD<-EBL_HD[EBL_HD$logFC < 1, ]



# Getting gene list for Venn Diagram plotting 

# UP
genes_up_EBL_AC_LP <- rownames(EBL_up_AC_LP)
genes_up_EBL_HD    <- rownames(EBL_up_HD)

# Down

genes_down_EBL_AC_LP <- rownames(EBL_down_AC_LP)
genes_down_EBL_HD    <- rownames(EBL_down_HD)



# List of items
z <- liz <- list(genes_up_EBL_AC_LP, genes_up_EBL_HD)
w <- list(genes_down_EBL_AC_LP, genes_down_EBL_HD)

#### Plotting Venn diagram with custom category names
ggVennDiagram(w, category.names = c("A",
                                   "B" ))+ 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")







#-###################################################-#
####- 17-- HeatMap -###################################################-#
####- 18-- GSEA -----------------------------------####
#-###################################################-#


EBL_AC_LP <- topTable(fit2, coef = "EBL_AC_LP", number = Inf,
                      adjust.method = "BH", p.value = 0.05, sort.by="B")  #4307
EBL_HD    <- topTable(fit2, coef = "EBL_HD", number = Inf,
                      adjust.method = "BH", p.value = 0.05, sort.by="B")  #5425


### EBL vs AC_LP (nonEBL)

# Create the vector based on ranking (use t o logFC) DEG-pvalue is not important, however we selected significant genes
ranked_genes <- EBL_HD$t  # Here we are using t
names(ranked_genes) <- rownames(EBL_HD)

### EBL vs HD 

# Ordering 
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

View(ranked_genes)
# Obtener los gene sets inmunológicos (C5)
msig_c5 <- msigdbr(species = "Bos taurus", category = "C5")

#msig_h <- msigdbr(species = "Bos taurus", category = "H")    # Hallmark


# Preparar la lista de pathways para fgsea
c5_pathways <- split(msig_c5$gene_symbol, msig_c5$gs_name)
#h_pathways <- split(msig_h$gene_symbol, msig_h$gs_name)
View(c5_pathways)
#

# Ejecutar GSEA
set.seed(123)  # Para reproducibilidad
fgsea_c5 <- fgsea(pathways = c5_pathways,
                  stats = ranked_genes,
                  minSize = 15,
                  maxSize = 500,
                  nperm = 10000)

#fgsea_h <- fgsea(pathways = h_pathways,
 #                stats = ranked_genes,
  #               minSize = 15,
   #              maxSize = 500,
    #             nperm = 10000)
print("all")
View(fgsea_c5)

# Ordenar por NES y valor p ajustado
fgsea_c5 <- fgsea_c5[order(fgsea_c5$padj, -abs(fgsea_c5$NES)),]

# Filtrar resultados significativos (FDR < 0.05)
#sig_5_1 <- fgsea_c5[fgsea_c5$padj < 0.05,]
#sig_5_1 <- fgsea_c5[fgsea_c5$padj < 0.05,]

sig_5_2 <- fgsea_c5[fgsea_c5$padj < 0.05,]

# Merging GSEA venn Diagram, we are going to merge just common EBL vs -
# nonEBL and EBL vs HD GSEA pathways


genes_sig_5_1 <- rownames(sig_5_1)  #EBL vs nonEBL
genes_sig_5_2 <- rownames(sig_5_2)  #EBL vs HD

# Converting to a list, ready for plotting
x <- list(genes_sig_5_2,genes_sig_5_1)

# Plotting Venn diagram with custom category names
ggVennDiagram(x,category.names = c("A",
                                   "B"
                                   ))+ 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")


# Actual Merging 
sig_5 <- merge(genes_sig_5_1,genes_sig_5_2, by="pathway")


# Separar pathways up y down
top_up <- sig_5 %>% 
  filter(NES > 0) %>% 
  arrange(desc(NES)) %>% 
  head(10)

top_down <- sig_5 %>% 
  filter(NES < 0) %>% 
  arrange(NES) %>%  # Ordenar de más negativo a menos
  head(10)

# Combinar ambos grupos
top_pathways <- bind_rows(top_up, top_down) %>%
  mutate(
    Direction = ifelse(NES > 0, "Upregulated", "Downregulated"),
    log10_padj = -log10(padj)  # Para tamaño/transparencia
  )


library(ggplot2)


#Modificando GOBP

top_pathways <- top_pathways %>%
  mutate(pathway = str_replace_all(pathway, "^GOBP_", "GOBP: "),
         pathway = str_replace_all(pathway, "^HP_", "HP: "),
         pathway = str_replace_all(pathway, "^GOCC_", "GOCC: "),
         pathway = str_replace_all(pathway, "^GOMF_", "GOMF: "),
         # pathway = str_replace_all(pathway, "^GOMF_", "GO:MF ")
         pathway = str_replace_all(pathway, "_", " "),
         pathway = str_replace(pathway, "^(.*?: )(.*)$", function(x) {
           prefix <- str_match(x, "^(.*?: )")[,1]
           body   <- str_match(x, "^.*?: (.*)$")[,2]
           paste0(prefix, str_to_sentence(tolower(body)))
         })
  )



## PLOT 
ggplot(top_pathways, aes(x = NES, y = reorder(pathway, NES))) +
  geom_col(
    aes(fill = Direction),
    width = 0.8,
    color = "gray40",
    linewidth  = 0.3
  ) +
  scale_fill_manual(
    values = c("Upregulated" = "#FC8D62", "Downregulated" = "#8DA0CB"),
    name = "FDR < 0.05"  # Cambiar título de la leyenda
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    x = "Normalized Enrichment Score",
    y = NULL,
    title = "Top 20 common Pathways"
  ) +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial"),
    axis.text.y = element_text(size = 11, face = "bold", color = "gray20"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "right",
    plot.title = element_text(size=14, face= "bold")
  )

#920 #660






