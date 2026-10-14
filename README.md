# S1: neutral binary trait with population structure

## Overview

This workflow simulates population structure and a neutral, genetically
influenced binary trait. Here, *neutral* means that the trait and its causal
QTLs do not affect fitness or reproductive success. It does not mean that the
trait is independent of genotype.

The workflow has three stages:

1. `S1_neutral_structure.slim` simulates the population history, QTLs,
   continuous liability, and binary phenotype.
2. `export_neutral_snps.py` simplifies the tree sequences to the selected
   individuals, overlays neutral SNPs, inserts QTLs at their genomic
   positions, and writes the genotype matrix.
3. `make_vae_input.py` aligns the genotype matrix, binary phenotype, and
   population labels by pedigree ID and checks the resulting arrays.

The current setup is a reduced two-chromosome test model rather than the final
eight-linkage-group Laupala model.

## Files

| File | Purpose |
|---|---|
| `S1_neutral_structure.slim` | Forward simulation and phenotype generation |
| `export_neutral_snps.py` | Tree-sequence simplification, recapitation, neutral SNP overlay, and genotype export |
| `make_vae_input.py` | Alignment and validation of `X`, `y`, and population labels |

## Current simulation settings

| Parameter | Value |
|---|---:|
| Random seed | 1024 |
| Simulation length | 5,000 generations |
| Populations | 10 serially diverged populations |
| Old population size | 4,000 |
| Middle population size | 3,000 |
| Young population size | 2,000 |
| Migration rate between neighbours | 0.1 |
| Chromosomes | 2 autosomes |
| Chromosome length | 2,500,000 positions |
| Recombination rate | `1e-8` per position per generation |
| QTLs | 6 per chromosome; 12 total |
| QTL allelic effects | `-0.05` or `+0.05` |
| QTL fitness effect | Neutral; fitness multiplier 1.0 |
| Target liability-scale heritability | 0.5 |
| Binary threshold | 0.0 |
| Sampled individuals | 300 per population; 3,000 total |
| Final variant columns | 25,000 per chromosome; 50,000 total |

The serial population history is:

```text
p1 -> p2 -> p3 -> p4 -> p5 -> p6 -> p7 -> p8 -> p9 -> p10
```

Populations `p1-p3`, `p4-p6`, and `p7-p10` represent the old, middle,
and young groups, respectively. Migration occurs only between neighbouring
populations.

The present migration rate of `0.1` and recombination rate of `1e-8` are the
values used for the current test. They should be reviewed before final
biological replicates. With a chromosome length of 2,500,000, `1e-8`
corresponds to approximately 0.025 Morgan, or 2.5 cM, per chromosome.

## Trait generation

For each sampled individual, SLiM calculates the additive genetic value from
the 12 QTLs and adds an independent environmental deviation:

```text
liability = genetic_value + environmental_deviation
```

The environmental variance is chosen to target liability-scale heritability
`h2 = 0.5`. The binary phenotype is then defined using the fixed threshold
zero:

```text
liability <= 0.0  -> phenotype_binary = 0
liability >  0.0  -> phenotype_binary = 1
```

Both liability and binary phenotype are retained. The binary phenotype is the
trait supplied to the current VAE workflow; liability is retained to audit the
simulation and realised heritability.

## Requirements

- SLiM 5
- Python 3
- NumPy
- pandas
- tskit
- PySLiM
- msprime

One possible environment setup is:

```bash
conda create -n slim-vae python=3.12 -y
conda activate slim-vae
pip install numpy pandas tskit pyslim msprime
```

## Step 1: run SLiM

From the directory containing the scripts:

```bash
slim S1_neutral_structure.slim
```

The default SLiM output directory is:

```text
results/S1_neutral_structure_10pop/
```

Important outputs include:

| Output | Description |
|---|---|
| `final_trees/` | One tree sequence per chromosome |
| `individuals.tsv` | Sample IDs, populations, genetic values, environmental deviations, liabilities, and binary phenotypes |
| `qtl_plus_allele_dosage.tsv` | QTL dosage for every sampled individual |
| `qtl_map.tsv` | QTL chromosome, position, effect, and neutrality status |
| `trait_summary.tsv` | Target and realised heritability, threshold, and overall class counts |
| `binary_population_summary.tsv` | Class 0 and class 1 counts within each population |
| `population_map.tsv` | Population split history |
| `population_summary.tsv` | Population sizes and genetic-value summaries through time |

Parameters may be overridden at runtime. For example:

```bash
slim -d SEED=1025 -d SAMPLE_PER_POP=300 S1_neutral_structure.slim
```

Use a different `OUT_DIR` for each replicate to prevent overwriting previous
results.

## Step 2: generate the genotype matrix

Run the post-processing program after SLiM finishes:

```bash
python export_neutral_snps.py \
    results/S1_neutral_structure_10pop/final_trees \
    --output-dir results/S1_neutral_structure_10pop/neutral_snps \
    --sample-individuals-tsv \
    results/S1_neutral_structure_10pop/individuals.tsv \
    --snps-per-chromosome 25000 \
    --recombination-rate 1e-8 \
    --qtl-dosage-tsv \
    results/S1_neutral_structure_10pop/qtl_plus_allele_dosage.tsv \
    --qtl-map-tsv \
    results/S1_neutral_structure_10pop/qtl_map.tsv \
    --write-combined-dosage
```

`--sample-individuals-tsv` ensures that Python uses exactly the individuals
selected by SLiM. Python does not draw a second, different sample.

### Composition of the genotype matrix

`--snps-per-chromosome 25000` specifies the final number of variant columns
per chromosome, including QTLs:

```text
Chromosome 1: 24,994 neutral SNPs + 6 QTLs = 25,000 columns
Chromosome 2: 24,994 neutral SNPs + 6 QTLs = 25,000 columns

Total: 49,988 neutral SNPs + 12 QTLs = 50,000 columns
```

The QTLs are inserted among the neutral SNPs according to chromosome and
position. They are not appended as the final 12 columns.

Post-processing outputs include:

| Output | Description |
|---|---|
| `combined_neutral_and_qtl_dosage.npy` | `uint8` genotype matrix with shape `(individuals, 50000)` and values 0, 1, or 2 |
| `combined_feature_map.tsv` | Exact matrix column, chromosome, position, feature type, and feature ID |
| `sampled_individuals.tsv` | Exact row order of the genotype matrix |
| `causal_qtl_dosage.npy` | Twelve QTL dosage columns before positional insertion |
| `causal_qtl_ids.tsv` | QTL dosage-column identifiers |
| `chromosome_1_neutral.trees` | Simplified chromosome 1 tree sequence with neutral sites |
| `chromosome_2_neutral.trees` | Simplified chromosome 2 tree sequence with neutral sites |
| `manifest.tsv` | Number of neutral sites and individuals per chromosome |

An older matrix with shape `(3000, 50012)` was produced by appending 12 QTLs
after 50,000 neutral SNPs. That layout is obsolete. The corrected matrix must
have shape `(3000, 50000)`.

## Step 3: prepare VAE inputs

Run:

```bash
python make_vae_input.py
```

The script aligns phenotype rows to genotype rows by matching:

```text
sampled_individuals.tsv: pedigree_id
individuals.tsv:         individual_id
```

The expected arrays are:

```text
X genotype:  (3000, 50000)
y phenotype: (3000,)
population:  (3000,)
```

No transpose is required. Rows are individuals and columns are genomic
features. For row `i`:

```text
X[i, :]       = genotype of individual i
y[i]          = binary phenotype of individual i
population[i] = population of individual i
```

The joint genotype-trait VAE receives both genotype and binary phenotype. If
the model accepts separate arguments, batches should have shapes:

```text
genotype_batch: (batch_size, 50000)
trait_batch:    (batch_size, 1)
```

If the model explicitly concatenates genotype and phenotype before encoding,
the joint input has shape `(batch_size, 50001)`.

## Locating the causal QTLs

QTL matrix columns must be read from `combined_feature_map.tsv`, because QTLs
are inserted at their genomic positions:

```python
import numpy as np
import pandas as pd

X = np.load(
    "results/S1_neutral_structure_10pop/neutral_snps/"
    "combined_neutral_and_qtl_dosage.npy",
    mmap_mode="r",
)

feature_map = pd.read_csv(
    "results/S1_neutral_structure_10pop/neutral_snps/"
    "combined_feature_map.tsv",
    sep="\t",
)

qtl_features = feature_map[
    feature_map["feature_type"] == "causal_qtl"
]
qtl_columns = qtl_features["column"].to_numpy()
X_qtl = X[:, qtl_columns]

print(qtl_features)
print(X_qtl.shape)  # (3000, 12)
```

Use this mapping to compare downstream feature importance or SHAP values with
the known causal variants.

## Current single-run result

The current run used 3,000 sampled individuals and took approximately four
minutes for simulation and post-processing. The realised trait statistics
were:

```text
target h2:          0.5
realised h2:        0.495203
binary threshold:   0.0
class 0:            1,224 individuals
class 1:            1,776 individuals
class 1 proportion: 0.592
```

The class 1 proportions across populations ranged from 0.550 to 0.643 in this
run. This is one run from one seed. No replicate simulations have yet been
completed, so these values do not quantify between-replicate variation.

## Reproducibility checks

Before VAE training, verify:

```python
assert X.shape == (3000, 50000)
assert y.shape == (3000,)
assert population.shape == (3000,)
assert np.isin(X, [0, 1, 2]).all()
assert np.isin(y, [0, 1]).all()
assert np.unique(population).size == 10
```

For each replicate, retain the random seed, SLiM settings, trait summary,
population-specific class proportions, QTL map, feature map, and exact sample
manifest. Train, validation, and test splits must use the same row indices for
`X`, `y`, and population labels.
