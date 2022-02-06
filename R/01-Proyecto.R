# Paquete a utilizar
library("recount3")

# Revisar los proyectos con datos de human con recount3
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

