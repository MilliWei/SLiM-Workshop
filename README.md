# S1: neutral binary trait with population structure

## Overview

This workflow simulates population structure and a neutral, genetically
influenced binary trait. Here, *neutral* means that the trait and its causal
QTLs do not affect fitness or reproductive success.
The workflow has three stages:

1. `S1_neutral_structure.slim` simulates the population history, QTLs,
   continuous liability, and binary phenotype.
2. `export_neutral_snps.py` simplifies the tree sequences to the selected
   individuals, overlays neutral SNPs, inserts QTLs at their genomic
   positions, and writes the genotype matrix.
3. `make_vae_input.py` aligns the genotype matrix, binary phenotype, and
   population labels by pedigree ID and checks the resulting arrays.

## Files

| File | Purpose |
|---|---|
| `S1_neutral_structure.slim` | Forward simulation and phenotype generation |
| `export_neutral_snps.py` | Tree-sequence simplification, recapitation, neutral SNP overlay, and genotype export |
| `make_vae_input.py` | Alignment and validation of `x`, `y`, and population labels |

## Current simulation settings

| Parameter | Value |
|---|---:|
| Random seed | 1024 |
| Simulation length | 5,000 generations |
| Populations | 10 serially diverged populations |
| Old population size | 4,000 |
| Middle population size | 3,000 |
| Young population size | 2,000 |
| Migration rate between neighbours | 0.01 |
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

The present migration rate of `0.01` and recombination rate of `1e-8` are the
values used for the current test. With a chromosome length of 2,500,000, `1e-8`
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

Both liability and binary phenotype are retained.

## Requirements

- SLiM 5
- Python 3
- NumPy
- pandas
- tskit
- PySLiM
- msprime

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
position.

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
