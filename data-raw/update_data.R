# -----------------------------------------------------------------------------
# As of 2022-11-20
# + by @luciorq and @jdblischak

# Original files:
# + data/dds_all_march2021.RDS ✅
# + data/fit_all_march2021.RDS ✅
# + data/fit_pdxOnly_march2021.RDS ✅
# + data/fit_primaryOnly_march2021.RDS ✅
# + data/fusions_2019.RDS ✅
# + data/gsva.RDS ✅
# + data/ssgsea.RDS ✅

# -----------------------------------------------------------------------------
# DDS object
# Remove DLBCL sample data from main objects
x <- readr::read_rds("data/dds_all_march2021.RDS")
sample_annot <- colData(x)
samples_to_keep <- colnames(
  x[, !stringr::str_detect(sample_annot$type, "DLBCL")]
)

y <- x[, !stringr::str_detect(sample_annot$type, "DLBCL")]


sample_annot <- colData(y)
sample_cols <- colnames(sample_annot)

# Remove empty column
sample_annot <- sample_annot[, !sample_cols %in% "X9"]
sample_cols <- colnames(sample_annot)
colnames(sample_annot)
unique(sample_annot$patient)
sum(sort(table(sample_annot$patient)))

# Replace "patient" column with "sample"
levels(sample_annot$type)
for (i in seq_along(sample_cols)) {
  if (isTRUE(is.factor(sample_annot[[i]]))) {
    sample_annot[[i]] <- forcats::fct_drop(sample_annot[[i]])
  }
}
y@colData <- sample_annot
y |>
  readr::write_rds("data/dds_all_nov2022.rds")
rm(x, y, i)

# -----------------------------------------------------------------------------
# LM Fit object ALL
x <- readr::read_rds("data/fit_all_march2021.RDS")

#' Remove DBLCL samples from fit object
remove_dlbcl_fit <- function(x) {
  for (i in seq_along(x)) {
    if (any(stringr::str_detect(colnames(x[[i]]), "DLBCL"))) {
      x[[i]] <- x[[i]][, !stringr::str_detect(colnames(x[[i]]), "DLBCL")]
    }
    if (any(stringr::str_detect(rownames(x[[i]]), "DLBCL"))) {
      x[[i]] <- x[[i]][!stringr::str_detect(rownames(x[[i]]), "DLBCL"), ]
    }
  }
  x[[10]] <- x[[10]][rownames(x[[10]]) %in% samples_to_keep, ]
  return(x)
}
remove_dlbcl_fit(x) |>
  readr::write_rds("data/fit_all_nov2022.rds")
rm(x)
# -----------------------------------------------------------------------------
# LM Fit object PDX
x <- readr::read_rds("data/fit_pdxOnly_march2021.RDS")
remove_dlbcl_fit(x) |>
  readr::write_rds("data/fit_pdx_only_nov2022.rds")
rm(x)
# -----------------------------------------------------------------------------
# LM Fit object PRIMARY
x <- readr::read_rds("data/fit_primaryOnly_march2021.RDS")
remove_dlbcl_fit(x) |>
  readr::write_rds("data/fit_primary_only_nov2022.rds")
rm(x)
# -----------------------------------------------------------------------------
# Fusions object
# TODO: @luciorq Change Patient column
# + and sample names
x <- readr::read_rds("data/fusions_2019.RDS")
names(x)
str(x$list)
x |>
  readr::write_rds("data/fusions_nov2022.rds")
rm(x)
# -----------------------------------------------------------------------------
# GSVA object
#' Remove DBLCL samples from GSVA object
remove_dlbcl_gsva <- function(x) {
  for (i in seq_along(x)) {
    x[[i]] <- x[[i]][, colnames(x[[i]]) %in% samples_to_keep]
  }
  return(x)
}
x <- readr::read_rds("data/gsva.RDS")
x |>
  remove_dlbcl_gsva() |>
  readr::write_rds("data/gsva_nov2022.rds")
rm(x)
# -----------------------------------------------------------------------------
# SSGSEA object
x <- readr::read_rds("data/ssgsea.RDS")
x |>
  remove_dlbcl_gsva() |>
  readr::write_rds("data/ssgsea_nov2022.rds")
rm(x)
# -----------------------------------------------------------------------------
# MSigDB GMTs
# -----------------------------------------------------------------------------
# + "data/msigdb/*"
# TODO: @luciorq Update msigdb, maybe use the R package to download the gmt
fs::dir_ls("data/msigdb/")
