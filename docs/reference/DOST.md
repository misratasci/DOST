# DOST: Distance-preserving Optimization for Spatial Transcriptomics

Runs DOST to identify spatial domains in spatial transcriptomics data.

## Usage

``` r
DOST(
  X,
  coords,
  R,
  selected_genes = "HVG",
  nGenes = 3000,
  neighborhood_threshold = 1,
  embedding_dim = 20,
  lambda = 0.025,
  lr = 16,
  max_iterations = 20,
  refinement = TRUE,
  refine_k = 30
)
```

## Arguments

- X:

  A sparse matrix or matrix of gene expression counts (Genes x Spots).

- coords:

  A data frame or matrix of spatial coordinates (Spots x 2). Rows of
  `coords` must match the columns of `X`.

- R:

  Integer. The number of clusters (domains) to identify.

- selected_genes:

  Character. Gene selection method: `"HVG"` (default), `"SVG"`, or
  `"all"`.

- nGenes:

  Integer. Number of genes to use if selection is enabled. Default is
  3000.

- neighborhood_threshold:

  Spatial adjacency threshold multiplier. Default is 1.

- embedding_dim:

  Dimension of the initial embedding. Default is 20.

- lambda:

  Regularization weight for spatial coherence. Default is 0.03.

- lr:

  Optimization learning rate. Default is 16.

- max_iterations:

  Max optimization steps. Default is 20.

- refinement:

  Logical. Whether to smooth the final labels using spatial neighbors.
  Recommended for tissues with laminar structures (e.g., DLPFC, mPFC).
  Default is `TRUE`.

- refine_k:

  Number of neighbors for label refinement. Default is 30.

## Value

A list containing:

- `Z`: The low-dimensional embedding matrix.

- `labels`: The final cluster labels.

- `losses`: The loss function history.

## Details

**Recommended Settings by Technology:**

The default parameters are optimized for **10x Visium** data. For other
technologies, we recommend the following adjustments based on our
benchmarks:

- **10x Visium:**

  - `embedding_dim = 20`

  - `neighborhood_threshold = 1`

  - `lambda = 0.03`

  - `nGenes = 3000` (HVG or SVG)

- **MERFISH:**

  - `embedding_dim = 10`

  - `neighborhood_threshold = 3` (Due to higher spatial resolution)

  - `lambda = 0.2`

  - `selected_genes = "all"` (If gene count is low)

- **STARmap:**

  - `embedding_dim = 10`

  - `neighborhood_threshold = 1`

  - `lambda = 0.15`

  - `selected_genes = "all"`

**Refinement Note:** The optional refinement step (`refinement = TRUE`)
effectively smooths boundaries and is highly recommended for tissues
with known laminar or layered structures (e.g., Cortex). For tissues
with discrete, punctate domains, you may wish to set this to `FALSE`.
