suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(harmony)
})

# Full post-QC condition-level integration workflow for the CRC single-cell atlas.
# This uses the full study-condition objects, not the 20k subsampled atlas.
#
# Input preprocessing is preserved:
# - GEO objects: previously LogNormalized; do not NormalizeData() again.
# - CELLxGENE_Chen_2021: supplied processed/scaled matrix is stored in RNA counts;
#   do not NormalizeData() or ScaleData() that supplied matrix.
# - CELLxGENE_Moorman2025: already normalized; do not NormalizeData() again.
#
# CCA uses Seurat anchor-based integration.
# Harmony learns PCA from non-Chen normalized data. When Chen is present,
# Chen cells are projected into the non-Chen PCA space using the PCA loadings.

atlas_root <- Sys.getenv(
  "CRC_ATLAS_ROOT",
  unset = "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas"
)

config <- list(
  source_dir = file.path(
    atlas_root,
    "full_seurat_QC_normalized_by_study_condition"
  ),
  hvg_dir = file.path(
    atlas_root,
    "full_atlas_HVG_2500"
  ),
  output_root = file.path(
    atlas_root,
    "Integration_by_condition"
  ),
  conditions = c(
    "Adenoma_polyp",
    "CRC",
    "Normal",
    "TAN"
  ),
  npcs = 50L,
  integration_dims = 1:30,
  downstream_dims = 1:30,
  resolutions = c(0.3, 0.5),
  run_cca = TRUE,
  run_harmony = TRUE,
  save_intermediate = TRUE,
  overwrite = FALSE
)

studies_by_condition <- list(
  Adenoma_polyp = c(
    "CELLxGENE_Chen_2021",
    "GSE201348",
    "GSE261388"
  ),
  CRC = c(
    "CELLxGENE_Chen_2021",
    "CELLxGENE_Moorman2025",
    "GSE132465",
    "GSE178341",
    "GSE200997",
    "GSE201348",
    "GSE221575",
    "GSE225857",
    "GSE231559",
    "GSE261388",
    "GSE294300",
    "GSE302903",
    "GSE330797"
  ),
  Normal = c(
    "CELLxGENE_Chen_2021",
    "GSE201348",
    "GSE221575",
    "GSE231559",
    "GSE261388"
  ),
  TAN = c(
    "CELLxGENE_Moorman2025",
    "GSE132465",
    "GSE178341",
    "GSE200997",
    "GSE225857",
    "GSE231559",
    "GSE294300"
  )
)

expected_cells <- c(
  Adenoma_polyp = 306396L,
  CRC = 702029L,
  Normal = 151684L,
  TAN = 287894L
)

chen_study <- "CELLxGENE_Chen_2021"
moorman_study <- "CELLxGENE_Moorman2025"

assert <- function(ok, message) {
  if (!isTRUE(ok)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

safe_token <- function(x) {
  gsub(
    "[^A-Za-z0-9._-]+",
    "_",
    as.character(x)
  )
}

source_file <- function(study, condition) {
  file.path(
    config$source_dir,
    paste0(
      study,
      "_",
      condition,
      "_QC_normalized_seurat.rds"
    )
  )
}

hvg_file <- function(condition) {
  file.path(
    config$hvg_dir,
    paste0(
      condition,
      "_2500_HVGs.csv"
    )
  )
}

condition_output_dir <- function(condition) {
  file.path(
    config$output_root,
    condition
  )
}

cca_output_file <- function(condition) {
  file.path(
    condition_output_dir(condition),
    paste0(
      condition,
      "_CCA_integrated.rds"
    )
  )
}

harmony_output_file <- function(condition) {
  file.path(
    condition_output_dir(condition),
    paste0(
      condition,
      "_Harmony_integrated.rds"
    )
  )
}

load_hvgs <- function(condition) {
  path <- hvg_file(condition)

  assert(
    file.exists(path),
    paste0(
      "HVG file not found: ",
      path
    )
  )

  x <- fread(path)

  assert(
    "gene" %in% names(x),
    paste0(
      condition,
      ": HVG file must contain a column named 'gene'"
    )
  )

  genes <- unique(
    as.character(
      x$gene
    )
  )

  genes <- genes[
    !is.na(genes) &
      trimws(genes) != ""
  ]

  assert(
    length(genes) > 0,
    paste0(
      condition,
      ": HVG file contains no usable genes"
    )
  )

  genes
}

read_source <- function(study, condition) {
  path <- source_file(
    study,
    condition
  )

  assert(
    file.exists(path),
    paste0(
      "Source object not found: ",
      path
    )
  )

  obj <- readRDS(path)

  assert(
    inherits(obj, "Seurat"),
    paste0(
      study,
      " / ",
      condition,
      ": source file is not a Seurat object"
    )
  )

  assert(
    "RNA" %in% names(obj@assays),
    paste0(
      study,
      " / ",
      condition,
      ": RNA assay is missing"
    )
  )

  assert(
    identical(
      rownames(obj@meta.data),
      colnames(obj)
    ),
    paste0(
      study,
      " / ",
      condition,
      ": metadata and cell order do not match"
    )
  )

  obj
}

selected_expression_layers <- function(obj, study) {
  layers <- SeuratObject::Layers(
    obj[["RNA"]]
  )

  if (identical(study, chen_study)) {
    selected <- grep(
      "^counts($|\\.)",
      layers,
      value = TRUE
    )

    assert(
      length(selected) > 0,
      paste0(
        study,
        ": expected Chen processed/scaled matrix in an RNA counts layer"
      )
    )

    return(selected)
  }

  selected <- grep(
    "^data($|\\.)",
    layers,
    value = TRUE
  )

  if (
    length(selected) == 0 &&
    identical(study, moorman_study)
  ) {
    selected <- grep(
      "^counts($|\\.)",
      layers,
      value = TRUE
    )
  }

  assert(
    length(selected) > 0,
    paste0(
      study,
      ": no normalized RNA expression layer found"
    )
  )

  selected
}

available_expression_features <- function(obj, study) {
  layers <- selected_expression_layers(
    obj,
    study
  )

  feature_sets <- vector(
    "list",
    length(layers)
  )

  for (i in seq_along(layers)) {
    mat <- SeuratObject::LayerData(
      object = obj,
      assay = "RNA",
      layer = layers[i]
    )

    feature_sets[[i]] <- rownames(mat)

    rm(mat)
    invisible(gc())
  }

  Reduce(
    intersect,
    feature_sets
  )
}

extract_expression <- function(
  obj,
  study,
  features
) {
  layers <- selected_expression_layers(
    obj,
    study
  )

  matrices <- vector(
    "list",
    length(layers)
  )

  seen_cells <- character()

  for (i in seq_along(layers)) {
    layer_name <- layers[i]

    mat <- SeuratObject::LayerData(
      object = obj,
      assay = "RNA",
      layer = layer_name
    )

    mat <- as(
      mat,
      "CsparseMatrix"
    )

    missing_features <- setdiff(
      features,
      rownames(mat)
    )

    assert(
      length(missing_features) == 0,
      paste0(
        study,
        ": ",
        length(missing_features),
        " integration features are missing from layer ",
        layer_name
      )
    )

    mat <- mat[
      features,
      ,
      drop = FALSE
    ]

    duplicate_cells <- intersect(
      seen_cells,
      colnames(mat)
    )

    assert(
      length(duplicate_cells) == 0,
      paste0(
        study,
        ": duplicate cells occur across selected expression layers"
      )
    )

    seen_cells <- c(
      seen_cells,
      colnames(mat)
    )

    matrices[[i]] <- mat
  }

  expression <- do.call(
    cbind,
    matrices
  )

  assert(
    !anyDuplicated(colnames(expression)),
    paste0(
      study,
      ": duplicated cell IDs after joining expression layers"
    )
  )

  expected_cells <- colnames(obj)

  assert(
    setequal(
      colnames(expression),
      expected_cells
    ),
    paste0(
      study,
      ": selected expression layers do not represent exactly the source-object cells"
    )
  )

  expression[
    features,
    expected_cells,
    drop = FALSE
  ]
}

shared_condition_hvgs <- function(
  condition,
  studies
) {
  candidate_hvgs <- load_hvgs(
    condition
  )

  feature_sets <- vector(
    "list",
    length(studies)
  )

  for (i in seq_along(studies)) {
    study <- studies[i]

    message(
      "Checking integration features: ",
      condition,
      " / ",
      study
    )

    obj <- read_source(
      study,
      condition
    )

    feature_sets[[i]] <- available_expression_features(
      obj,
      study
    )

    rm(obj)
    invisible(gc())
  }

  common_features <- Reduce(
    intersect,
    feature_sets
  )

  features <- candidate_hvgs[
    candidate_hvgs %in%
      common_features
  ]

  assert(
    length(features) >= 500,
    paste0(
      condition,
      ": only ",
      length(features),
      " condition HVGs are shared across all contributing studies"
    )
  )

  message(
    condition,
    ": ",
    length(features),
    " shared integration HVGs"
  )

  features
}

make_minimal_expression_object <- function(
  study,
  condition,
  features
) {
  source <- read_source(
    study,
    condition
  )

  expression <- extract_expression(
    source,
    study,
    features
  )

  cells <- colnames(
    expression
  )

  metadata <- source@meta.data[
    cells,
    ,
    drop = FALSE
  ]

  metadata$integration_study <- study
  metadata$integration_condition <- condition

  zero_counts <- Matrix::sparseMatrix(
    i = integer(0),
    j = integer(0),
    dims = dim(expression),
    dimnames = dimnames(expression)
  )

  obj <- Seurat::CreateSeuratObject(
    counts = zero_counts,
    assay = "RNA",
    project = paste0(
      safe_token(condition),
      "_",
      safe_token(study)
    ),
    meta.data = metadata
  )

  obj <- SeuratObject::SetAssayData(
    object = obj,
    assay = "RNA",
    layer = "data",
    new.data = expression
  )

  SeuratObject::VariableFeatures(
    obj[["RNA"]]
  ) <- features

  rm(
    source,
    expression,
    zero_counts,
    metadata
  )

  invisible(gc())

  obj
}

combine_metadata <- function(metadata_list) {
  pieces <- lapply(
    metadata_list,
    function(md) {
      as.data.table(
        md,
        keep.rownames = "cell"
      )
    }
  )

  combined <- rbindlist(
    pieces,
    use.names = TRUE,
    fill = TRUE
  )

  assert(
    !anyDuplicated(combined$cell),
    "Duplicate cell IDs detected while combining metadata"
  )

  combined <- as.data.frame(
    combined,
    stringsAsFactors = FALSE
  )

  rownames(combined) <- combined$cell
  combined$cell <- NULL

  combined
}

build_combined_normalized_object <- function(
  studies,
  condition,
  features
) {
  expression_list <- vector(
    "list",
    length(studies)
  )

  metadata_list <- vector(
    "list",
    length(studies)
  )

  seen_cells <- character()

  for (i in seq_along(studies)) {
    study <- studies[i]

    message(
      "Loading normalized expression: ",
      condition,
      " / ",
      study
    )

    source <- read_source(
      study,
      condition
    )

    expression <- extract_expression(
      source,
      study,
      features
    )

    cells <- colnames(
      expression
    )

    duplicate_cells <- intersect(
      seen_cells,
      cells
    )

    assert(
      length(duplicate_cells) == 0,
      paste0(
        condition,
        ": duplicate cell IDs across studies"
      )
    )

    seen_cells <- c(
      seen_cells,
      cells
    )

    metadata <- source@meta.data[
      cells,
      ,
      drop = FALSE
    ]

    metadata$integration_study <- study
    metadata$integration_condition <- condition

    expression_list[[i]] <- expression
    metadata_list[[i]] <- metadata

    rm(
      source,
      expression,
      metadata
    )

    invisible(gc())
  }

  expression <- do.call(
    cbind,
    expression_list
  )

  metadata <- combine_metadata(
    metadata_list
  )

  metadata <- metadata[
    colnames(expression),
    ,
    drop = FALSE
  ]

  assert(
    identical(
      rownames(metadata),
      colnames(expression)
    ),
    paste0(
      condition,
      ": combined metadata and expression order do not match"
    )
  )

  zero_counts <- Matrix::sparseMatrix(
    i = integer(0),
    j = integer(0),
    dims = dim(expression),
    dimnames = dimnames(expression)
  )

  obj <- Seurat::CreateSeuratObject(
    counts = zero_counts,
    assay = "RNA",
    project = paste0(
      safe_token(condition),
      "_combined"
    ),
    meta.data = metadata
  )

  obj <- SeuratObject::SetAssayData(
    object = obj,
    assay = "RNA",
    layer = "data",
    new.data = expression
  )

  SeuratObject::VariableFeatures(
    obj[["RNA"]]
  ) <- features

  rm(
    expression_list,
    metadata_list,
    expression,
    metadata,
    zero_counts
  )

  invisible(gc())

  obj
}

run_cca <- function(
  condition,
  studies,
  features
) {
  out_dir <- condition_output_dir(
    condition
  )

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  output_file <- cca_output_file(
    condition
  )

  if (
    file.exists(output_file) &&
    !config$overwrite
  ) {
    message(
      "CCA output exists; skipping: ",
      output_file
    )

    return(
      invisible(output_file)
    )
  }

  message(
    "\nCCA integration: ",
    condition
  )

  object_list <- vector(
    "list",
    length(studies)
  )

  names(object_list) <- studies

  for (i in seq_along(studies)) {
    study <- studies[i]

    message(
      "Preparing CCA input: ",
      study
    )

    object_list[[i]] <- make_minimal_expression_object(
      study,
      condition,
      features
    )
  }

  all_cells <- unlist(
    lapply(
      object_list,
      colnames
    ),
    use.names = FALSE
  )

  assert(
    !anyDuplicated(all_cells),
    paste0(
      condition,
      ": duplicate cell IDs across CCA input objects"
    )
  )

  anchors <- Seurat::FindIntegrationAnchors(
    object.list = object_list,
    anchor.features = features,
    reduction = "cca",
    dims = config$integration_dims,
    normalization.method = "LogNormalize",
    verbose = TRUE
  )

  if (config$save_intermediate) {
    saveRDS(
      anchors,
      file.path(
        out_dir,
        paste0(
          condition,
          "_CCA_anchors.rds"
        )
      ),
      compress = FALSE
    )
  }

  integrated <- Seurat::IntegrateData(
    anchorset = anchors,
    new.assay.name = "integrated",
    features.to.integrate = features,
    dims = config$integration_dims,
    normalization.method = "LogNormalize",
    verbose = TRUE
  )

  rm(
    anchors,
    object_list
  )

  invisible(gc())

  SeuratObject::DefaultAssay(
    integrated
  ) <- "integrated"

  integrated <- Seurat::ScaleData(
    integrated,
    features = features,
    verbose = FALSE
  )

  integrated <- Seurat::RunPCA(
    integrated,
    assay = "integrated",
    features = features,
    npcs = config$npcs,
    reduction.name = "pca.cca",
    reduction.key = "CCAPC_",
    verbose = TRUE
  )

  integrated <- Seurat::RunUMAP(
    integrated,
    reduction = "pca.cca",
    dims = config$downstream_dims,
    reduction.name = "umap.cca",
    reduction.key = "CCAUMAP_",
    verbose = TRUE
  )

  integrated <- Seurat::FindNeighbors(
    integrated,
    reduction = "pca.cca",
    dims = config$downstream_dims,
    graph.name = c(
      "cca_nn",
      "cca_snn"
    ),
    verbose = TRUE
  )

  integrated <- Seurat::FindClusters(
    integrated,
    graph.name = "cca_snn",
    resolution = config$resolutions,
    verbose = TRUE
  )

  assert(
    ncol(integrated) ==
      unname(expected_cells[condition]),
    paste0(
      condition,
      ": CCA contains ",
      ncol(integrated),
      " cells; expected ",
      unname(expected_cells[condition])
    )
  )

  saveRDS(
    integrated,
    output_file,
    compress = FALSE
  )

  message(
    "CCA saved: ",
    output_file
  )

  rm(integrated)
  invisible(gc())

  invisible(output_file)
}

build_harmony_pca <- function(
  condition,
  studies,
  features
) {
  has_chen <- chen_study %in%
    studies

  non_chen_studies <- setdiff(
    studies,
    chen_study
  )

  assert(
    length(non_chen_studies) > 0,
    paste0(
      condition,
      ": no non-Chen study is available for PCA learning"
    )
  )

  non_chen <- build_combined_normalized_object(
    non_chen_studies,
    condition,
    features
  )

  non_chen <- Seurat::ScaleData(
    non_chen,
    features = features,
    verbose = FALSE
  )

  non_chen <- Seurat::RunPCA(
    non_chen,
    assay = "RNA",
    features = features,
    npcs = config$npcs,
    reduction.name = "pca",
    reduction.key = "PC_",
    verbose = TRUE
  )

  pca_non_chen <- SeuratObject::Embeddings(
    non_chen,
    reduction = "pca"
  )

  loadings <- SeuratObject::Loadings(
    non_chen,
    reduction = "pca"
  )

  loadings <- loadings[
    features,
    seq_len(config$npcs),
    drop = FALSE
  ]

  metadata_list <- list(
    non_chen@meta.data
  )

  if (has_chen) {
    message(
      "Projecting Chen cells into non-Chen PCA space..."
    )

    chen <- read_source(
      chen_study,
      condition
    )

    chen_scaled <- extract_expression(
      chen,
      chen_study,
      features
    )

    chen_metadata <- chen@meta.data[
      colnames(chen_scaled),
      ,
      drop = FALSE
    ]

    chen_metadata$integration_study <-
      chen_study

    chen_metadata$integration_condition <-
      condition

    chen_pca <- Matrix::t(
      chen_scaled
    ) %*%
      loadings

    chen_pca <- as.matrix(
      chen_pca
    )

    colnames(chen_pca) <- colnames(
      pca_non_chen
    )

    pca_all <- rbind(
      pca_non_chen,
      chen_pca
    )

    metadata_list[[2]] <- chen_metadata

    rm(
      chen,
      chen_scaled,
      chen_metadata,
      chen_pca
    )
  } else {
    pca_all <- pca_non_chen
  }

  metadata <- combine_metadata(
    metadata_list
  )

  metadata <- metadata[
    rownames(pca_all),
    ,
    drop = FALSE
  ]

  assert(
    identical(
      rownames(metadata),
      rownames(pca_all)
    ),
    paste0(
      condition,
      ": Harmony PCA and metadata order do not match"
    )
  )

  rm(
    non_chen,
    pca_non_chen
  )

  invisible(gc())

  list(
    embeddings = pca_all,
    loadings = loadings,
    metadata = metadata
  )
}

run_harmony <- function(
  condition,
  studies,
  features
) {
  out_dir <- condition_output_dir(
    condition
  )

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  output_file <- harmony_output_file(
    condition
  )

  if (
    file.exists(output_file) &&
    !config$overwrite
  ) {
    message(
      "Harmony output exists; skipping: ",
      output_file
    )

    return(
      invisible(output_file)
    )
  }

  message(
    "\nHarmony integration: ",
    condition
  )

  pca_bundle <- build_harmony_pca(
    condition,
    studies,
    features
  )

  pca_embeddings <- getElement(
    pca_bundle,
    "embeddings"
  )

  pca_loadings <- getElement(
    pca_bundle,
    "loadings"
  )

  metadata <- getElement(
    pca_bundle,
    "metadata"
  )

  cells <- rownames(
    pca_embeddings
  )

  dummy_counts <- Matrix::sparseMatrix(
    i = integer(0),
    j = integer(0),
    dims = c(
      1,
      length(cells)
    ),
    dimnames = list(
      "integration_placeholder",
      cells
    )
  )

  obj <- Seurat::CreateSeuratObject(
    counts = dummy_counts,
    assay = "RNA",
    project = paste0(
      safe_token(condition),
      "_Harmony"
    ),
    meta.data = metadata
  )

  obj[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = pca_embeddings,
    loadings = pca_loadings,
    key = "PC_",
    assay = "RNA"
  )

  rm(
    pca_bundle,
    pca_embeddings,
    pca_loadings,
    metadata,
    dummy_counts
  )

  invisible(gc())

  obj <- harmony::RunHarmony(
    object = obj,
    group.by.vars = "integration_study",
    reduction.use = "pca",
    dims.use = config$integration_dims,
    reduction.save = "harmony",
    verbose = TRUE
  )

  obj <- Seurat::RunUMAP(
    obj,
    reduction = "harmony",
    dims = config$downstream_dims,
    reduction.name = "umap.harmony",
    reduction.key = "HARMONYUMAP_",
    verbose = TRUE
  )

  obj <- Seurat::FindNeighbors(
    obj,
    reduction = "harmony",
    dims = config$downstream_dims,
    graph.name = c(
      "harmony_nn",
      "harmony_snn"
    ),
    verbose = TRUE
  )

  obj <- Seurat::FindClusters(
    obj,
    graph.name = "harmony_snn",
    resolution = config$resolutions,
    verbose = TRUE
  )

  assert(
    ncol(obj) ==
      unname(expected_cells[condition]),
    paste0(
      condition,
      ": Harmony contains ",
      ncol(obj),
      " cells; expected ",
      unname(expected_cells[condition])
    )
  )

  saveRDS(
    obj,
    output_file,
    compress = FALSE
  )

  message(
    "Harmony saved: ",
    output_file
  )

  rm(obj)
  invisible(gc())

  invisible(output_file)
}

validate_inputs <- function() {
  assert(
    dir.exists(config$source_dir),
    paste0(
      "Source directory not found: ",
      config$source_dir
    )
  )

  assert(
    dir.exists(config$hvg_dir),
    paste0(
      "HVG directory not found: ",
      config$hvg_dir
    )
  )

  for (condition in config$conditions) {
    studies <- getElement(
      studies_by_condition,
      condition
    )

    for (study in studies) {
      path <- source_file(
        study,
        condition
      )

      assert(
        file.exists(path),
        paste0(
          "Missing source file: ",
          path
        )
      )
    }

    assert(
      file.exists(
        hvg_file(condition)
      ),
      paste0(
        "Missing HVG file for ",
        condition
      )
    )
  }

  invisible(TRUE)
}

run_condition <- function(condition) {
  studies <- getElement(
    studies_by_condition,
    condition
  )

  assert(
    length(studies) >= 2,
    paste0(
      condition,
      ": cross-study integration requires at least two studies"
    )
  )

  message(
    "\nCondition: ",
    condition
  )

  message(
    "Studies: ",
    paste(
      studies,
      collapse = ", "
    )
  )

  features <- shared_condition_hvgs(
    condition,
    studies
  )

  out_dir <- condition_output_dir(
    condition
  )

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  fwrite(
    data.table(
      rank = seq_along(features),
      gene = features
    ),
    file.path(
      out_dir,
      paste0(
        condition,
        "_shared_integration_HVGs.csv"
      )
    )
  )

  if (config$run_cca) {
    run_cca(
      condition,
      studies,
      features
    )
  }

  if (config$run_harmony) {
    run_harmony(
      condition,
      studies,
      features
    )
  }

  invisible(gc())
}

validate_inputs()

dir.create(
  config$output_root,
  recursive = TRUE,
  showWarnings = FALSE
)

manifest <- rbindlist(
  lapply(
    config$conditions,
    function(condition) {
      studies <- getElement(
        studies_by_condition,
        condition
      )

      data.table(
        condition = condition,
        study = studies,
        source_file = vapply(
          studies,
          function(study) {
            source_file(
              study,
              condition
            )
          },
          character(1)
        )
      )
    }
  )
)

fwrite(
  manifest,
  file.path(
    config$output_root,
    "full_atlas_integration_manifest.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    config$output_root,
    "full_atlas_integration_sessionInfo.txt"
  )
)

for (condition in config$conditions) {
  run_condition(
    condition
  )
}

message(
  "\nFull condition-level CCA/Harmony integration completed."
)

message(
  paste0(
    "Serrated_polyp is intentionally excluded from cross-study integration ",
    "because only one study contributes to that condition."
  )
)
