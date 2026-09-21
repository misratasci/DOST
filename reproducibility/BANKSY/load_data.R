
load_DLPFC_sample <- function(slice.id, dir.input) {
  filename <- paste0(slice.id, "_filtered_feature_bc_matrix.h5")
  data.dir <- file.path(dir.input, slice.id)
  sp_data <- Seurat::Load10X_Spatial(data.dir, filename = filename, filter.matrix = FALSE)
  #add the annotations
  df_meta <- read.table(file.path(data.dir, 'gt', 'tissue_positions_list_GTs.txt'),
                        sep = ",", row.names = 1)
  common_cells <- colnames(sp_data[["Spatial"]]) %in% rownames(df_meta)
  sp_data <- sp_data[, common_cells]
  layer.data <- data.frame()
  layers <- c('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'WM')
  for (l in layers) {
    filename <- paste0(slice.id, "_", l, "_barcodes.txt")
    filename <- file.path(data.dir, 'gt', 'layered', filename)
    if (!file.exists(filename)) next
    data.temp <- read.table(filename)
    data.temp <- data.frame(barcode = data.temp[,1], layer = l, row.names = data.temp[,1])
    layer.data <- rbind(layer.data, data.temp)
  }
  sp_data <- SeuratObject::AddMetaData(sp_data,
                                       metadata = df_meta['V3'],
                                       col.name = 'row')
  sp_data <- SeuratObject::AddMetaData(sp_data,
                                       metadata = df_meta['V4'],
                                       col.name = 'col')
  sp_data <- SeuratObject::AddMetaData(sp_data,
                                       metadata = layer.data['layer'],
                                       col.name = 'layers')
  return (sp_data)
}

load_BC_sample <- function(dir.input, section_no = 1) {
  if (section_no != 1 & section_no != 2) {
    stop("section_no must be 1 or 2")
  }
  filename <- paste0("section", section_no, "_filtered_feature_bc_matrix.h5")
  sp_data <- Seurat::Load10X_Spatial(file.path(dir.input, paste0("section", section_no)),
                                     filename = filename, filter.matrix = FALSE)
  if (section_no == 1) {
    df_meta <- read.table(file.path(dir.input, paste0("section", section_no),
                                    "gt", "gold_metadata.tsv"),
                          sep = "\t", row.names = 1, header = TRUE)
    sp_data <- SeuratObject::AddMetaData(sp_data,
                                         metadata = df_meta[1],
                                         col.name = 'annot_type')
    sp_data <- SeuratObject::AddMetaData(sp_data,
                                         metadata = df_meta[2],
                                         col.name = 'fine_annot_type')
  }
  return(sp_data)
}

load_mMAMP_sample <- function(dir.input, section = "MA") {
  if (!(section %in% c("MA", "MP"))) {
    stop("section must be 'MA' or 'MP'")
  }
  dir.input <- file.path(dir.input, section)
  filename <- paste0(section, "_filtered_feature_bc_matrix.h5")
  sp_data <- Seurat::Load10X_Spatial(dir.input, filename = filename, filter.matrix = FALSE)
  df_meta <- read.table(file.path(dir.input, "metadata.tsv"),
                        sep = "\t", row.names = 1, header = TRUE)
  for (col in colnames(df_meta)) {
    sp_data <- SeuratObject::AddMetaData(sp_data,
                                         metadata = df_meta[col],
                                         col.name = col)
  }
  return(sp_data)
}
