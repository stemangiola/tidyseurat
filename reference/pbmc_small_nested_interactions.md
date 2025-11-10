# Intercellular ligand-receptor interactions for 38 ligands from a single cell RNA-seq cluster.

A dataset containing ligand-receptor interactions within a sample. There
are 38 ligands from a single cell cluster versus 35 receptors in 6 other
clusters.

## Usage

``` r
data(pbmc_small_nested_interactions)
```

## Format

A \`tibble\` containing 100 rows and 9 columns. Cells are a subsample of
the PBMC dataset of 2,700 single cells. Cell interactions were
identified with \`SingleCellSignalR\`.

- sample:

  sample identifier

- ligand:

  cluster and ligand identifier

- receptor:

  cluster and receptor identifier

- ligand.name:

  ligand name

- receptor.name:

  receptor name

- origin:

  cluster containing ligand

- destination:

  cluster containing receptor

- interaction.type:

  type of interation, paracrine or autocrine

- LRscore:

  interaction score

## Source

<https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html>

## Value

\`tibble\`
