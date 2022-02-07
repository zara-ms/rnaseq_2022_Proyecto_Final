## Paquetes a utilizar
# Obtención de datos del proyecto
library("recount3")

# Normalización de datos
library("edgeR")

# Gráficas
library("ggplot2")

# Análisis diferencial
library("limma")

# Creación de un heatmap
library("pheatmap")

# Revisar los proyectos con datos de humano con recount3
human_projects <- available_projects()

# Utilizar un proyecto de interés
# Proyecto sobre el canal iónico PIEZO1 y su relación con el glioma
proj_info <- subset(
  human_projects,
  project == "SRP140640" & project_type == "data_sources"
)

# Obtener infomarción del proyecto creando un objeto de tipo
# RangedSummarizedExperiment (RSE)
rse_gene_SRP140640 <- create_rse(proj_info)

# Obtener los valores del número de cuentas normales
assay(rse_gene_SRP140640, "counts") <- compute_read_counts(rse_gene_SRP140640)

# Procesar y expandir la información así como la provee sra
rse_gene_SRP140640 <- expand_sra_attributes(rse_gene_SRP140640)

# Visualizar la información obtenida tras la expansión
colData(rse_gene_SRP140640)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP140640)))
]

table(rse_gene_SRP140640$sra_attribute.cell_line)
table(rse_gene_SRP140640$sra_attribute.shRNA)
table(rse_gene_SRP140640$sra_attribute.source_name)
table(rse_gene_SRP140640$sra_attribute.cell_type)

# Crear variables categóricas para el modelo estadístico
rse_gene_SRP140640$sra_attribute.cell_line <-
  factor(rse_gene_SRP140640$sra_attribute.cell_line)
rse_gene_SRP140640$sra_attribute.shRNA <-
  factor(rse_gene_SRP140640$sra_attribute.shRNA)
rse_gene_SRP140640$sra_attribute.source_name <-
  factor(rse_gene_SRP140640$sra_attribute.source_name)

# Visualizar información
summary(as.data.frame(colData(rse_gene_SRP140640)[
  ,
  grepl("^sra_attribute.[cell_line|shRNA|source_name]",
        colnames(colData(rse_gene_SRP140640)))]))

# Calcular la proporción de lecturas asignadas a genes
rse_gene_SRP140640$assigned_gene_prop <-
  rse_gene_SRP140640$recount_qc.gene_fc_count_all.assigned /
  rse_gene_SRP140640$recount_qc.gene_fc_count_all.total

## Exploración de datos para decidir si se deben eliminar muestras
summary(rse_gene_SRP140640$assigned_gene_prop)
hist(rse_gene_SRP140640$assigned_gene_prop)
table(rse_gene_SRP140640$assigned_gene_prop < 0.5)

# Observar si existe diferencia entre lineas celulares
with(colData(rse_gene_SRP140640), tapply(assigned_gene_prop,
                                         sra_attribute.cell_line, summary))

# Observar si existe diferencia entre diferentes shRNAs usados
with(colData(rse_gene_SRP140640), tapply(assigned_gene_prop,
                                         sra_attribute.shRNA, summary))

# Calcular niveles medios de expresión
gene_means <- rowMeans(assay(rse_gene_SRP140640, "counts"))
summary(gene_means)

## Normalización de datos para el análisis diferencial
dge <- DGEList(
  counts = assay(rse_gene_SRP140640, "counts"),
  genes = rowData(rse_gene_SRP140640)
)
dge <- calcNormFactors(dge)

# Visualización de datos
ggplot(as.data.frame(colData(rse_gene_SRP140640)),
       aes(y = assigned_gene_prop, x = sra_attribute.cell_line)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Cell Line")

ggplot(as.data.frame(colData(rse_gene_SRP140640)),
       aes(y = assigned_gene_prop, x = sra_attribute.shRNA)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("shRNA")

# Modelo estadístico
mod <- model.matrix(~ sra_attribute.shRNA + sra_attribute.cell_line
                    + assigned_gene_prop,
                    data = colData(rse_gene_SRP140640)
)
colnames(mod)

# Análisis de expresión diferencial
vGene <- voom(dge, mod, plot = TRUE)

eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_results,
  coef = 2,  # Valores para el shRNA
  number = nrow(rse_gene_SRP140640),
  sort.by = "none"
)
dim(de_results)

# Genes diferencialmente expresados entre diferentes shRNAs
table(de_results$adj.P.Val < 0.05)

# Visualizar resultados estadísticos
plotMA(eb_results, coef = 2)
volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)

# Crear un heatmap con los genes más diferencialmente expresados
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]

dataf <- as.data.frame(colData(rse_gene_SRP140640)[, c("sra_attribute.shRNA",
                                                       "sra_attribute.cell_line"
                                                       )])

colnames(dataf) <- c("shRNA", "CellLine")
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = dataf
)
