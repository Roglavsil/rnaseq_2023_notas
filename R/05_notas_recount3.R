## Load recount3 R package
library("recount3")
available_projects(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")


## Revisemos todos los proyectos con datos de humano en recount3
human_projects <- available_projects()

## Vemos el url actual.
getOption("recount3_url", "http://duffel.rail.bio/recount3")

# Cambiar URL.
options(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")

# Comprobar.
getOption("recount3_url", "http://duffel.rail.bio/recount3")

## Encuentra tu proyecto de interés. Aquí usaremos
## SRP009615 de ejemplo
proj_info <- subset(
  human_projects,
  project == "SRP009615" & project_type == "data_sources"
)
## Crea un objetio de tipo RangedSummarizedExperiment (RSE)
## con la información a nivel de genes
rse_gene_SRP009615 <- create_rse(proj_info)

## Convirtamos las cuentas por nucleotido a cuentas por lectura
## usando compute_read_counts().
## Para otras transformaciones como RPKM y TPM, revisa transform_counts().
assay(rse_gene_SRP009615, "counts") <- compute_read_counts(rse_gene_SRP009615)

## Para este estudio en específico, hagamos más fácil de usar la
## información del experimento
rse_gene_SRP009615 <- expand_sra_attributes(rse_gene_SRP009615)
colData(rse_gene_SRP009615)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP009615)))
]

