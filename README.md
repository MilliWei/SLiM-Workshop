#  SLiM population-structure and quantitative-trait simulations

## Overview

This workflow simulates population structure and a neutral, genetically
influenced binary trait. Here, *neutral* means that the trait and its causal
QTLs do not affect fitness or reproductive success.
The workflow generates:

- a genotype matrix `X` with individuals in rows and variants in columns;
- a continuous liability and globally balanced binary phenotype `y`;
- population labels;
- known causal-QTL positions and effects; and
- diagnostics for population structure, heritability and class balance.

## Current scenarios

| Scenario | Trait selection | Population-specific optimum | Main purpose |
|---|---|---|---|
| S1 | None; QTLs are neutral for fitness | None | Neutral-trait benchmark under population structure |
| S2 Design A | Gaussian stabilising selection | `theta_i = 0` | Global-adaptation benchmark without an imposed environmental gradient |
| S2 Design B | Gaussian stabilising selection | `theta_i = (i - 1) * delta` | Local-adaptation benchmark with demographic confounding |

The current S1 script is `S1_neutral_structure_10pop_percent_qtl.slim`. The
S2 script is `S2_trait_no_imposed_confounding.slim`. The S2 filename is retained
for compatibility, although the script now supports both Designs A and B.

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
| QTL allelic effects | `-0.05` or `+0.05` or `-4` or `+4` |
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

| Populations | Group | Population size |
|---|---|---:|
| p1--p3 | Old | 4,000 |
| p4--p6 | Middle | 3,000 |
| p7--p10 | Young | 2,000 |

The present migration rate of `0.01` and recombination rate of `1e-8` are the
values used for the current test. With a chromosome length of 2,500,000, `1e-8`
corresponds to approximately 0.025 Morgan, or 2.5 cM, per chromosome.

## S1: neutral trait with population structure

S1 uses 12 fixed QTLs with allelic effects `-0.05` and `+0.05`. The QTLs
contribute to the recorded phenotype but have no effect on survival or
reproduction. Consequently, the population history evolves neutrally with
respect to the trait.

At the final generation, an environmental deviation is added to the additive
genetic value:

<p align="center"><i>P</i><sub>j</sub> = <i>G</i><sub>j</sub> + <i>E</i><sub>j</sub></p>

Environmental variance targets liability-scale `h2=0.5`. The sampled
liabilities are ranked globally, with the lower half assigned to class 0 and
the upper half assigned to class 1. Population labels are not used to create
the binary phenotype.

S1 is used to test whether population structure alone produces spurious SNP
signals when the simulated trait does not affect fitness.

## S2: quantitative trait under stabilising selection

S2 retains the same demographic history and 12 known causal-QTL positions.
Each QTL has fixed allelic effects `-4` and `+4`. With 12 diploid QTLs, the
theoretical additive genetic-value range is `[-96, 96]`.

For individual $j$, the expressed phenotype is:

| Quantity | Definition |
|---|---|
| Additive genetic value | <i>G</i><sub>j</sub> = sum of the 12 diploid QTL effects |
| Environmental variance | <i>V</i><sub>E</sub> = <i>V</i><sub>A, within</sub> × (1 − <i>h</i><sup>2</sup>) / <i>h</i><sup>2</sup> |
| Environmental deviation | <i>E</i><sub>j</sub> ∼ Normal(0, <i>V</i><sub>E</sub>) |
| Expressed phenotype | <i>P</i><sub>j</sub> = <i>G</i><sub>j</sub> + <i>E</i><sub>j</sub> |

The target within-population heritability is `h2=0.5`. Environmental variance
is calculated from pooled within-population additive genetic variance, so
adaptive differences among population means do not inflate `V_E`.

Fitness follows standard Gaussian stabilising selection:

<p align="center">
<i>w</i><sub>ij</sub> = exp[−(<i>P</i><sub>ij</sub> − θ<sub>i</sub>)<sup>2</sup> / (2<i>V</i><sub>S</sub>)]
</p>

`SELECTION_VARIANCE` represents <i>V</i><sub>S</sub>, not its square root. Smaller values
produce stronger selection around the optimum.

### Design A: common optimum

All populations experience the same optimum:

```text
theta_i = 0
V_S in {10, 22, 50}
```

Fifty replicates are planned for each value of `V_S`, giving 150 Design A
replicates.

### Design B: linear environmental gradient

The optimum increases by a fixed amount between adjacent populations:

<p align="center">
θ<sub>i</sub> = (i − 1)δ, &nbsp;&nbsp;&nbsp;
δ ∈ {4.6, 10}, &nbsp;&nbsp;&nbsp;
<i>V</i><sub>S</sub> ∈ {10, 22, 50}
</p>

Therefore, p<sub>1</sub> begins at θ<sub>1</sub> = 0:

| δ | Population optima from p<sub>1</sub> to p<sub>10</sub> |
|---:|---|
| 4.6 | 0, 4.6, 9.2, 13.8, 18.4, 23.0, 27.6, 32.2, 36.8, 41.4 |
| 10 | 0, 10, 20, 30, 40, 50, 60, 70, 80, 90 |

The ±4 QTL effects make the largest optimum of 90 genetically reachable.
Each of the six δ × <i>V</i><sub>S</sub> combinations has 50 planned replicates,
giving 300 Design B replicates.

| Design | δ | <i>V</i><sub>S</sub> | Replicates |
|---|---:|---:|---:|
| A | 0 | 10, 22, 50 | 3 × 50 = 150 |
| B | 4.6, 10 | 10, 22, 50 | 2 × 3 × 50 = 300 |

The complete S2 design therefore contains 450 planned replicates.

## Binary phenotype

The continuous liability is used during evolution and retained in the output.
At final sampling, all 3,000 individuals are ranked together by liability:

```text
lower 1,500 liabilities -> phenotype_binary = 0
upper 1,500 liabilities -> phenotype_binary = 1
```

This guarantees global 50/50 balance. It does not impose balance separately
within each population. Population-specific class proportions are recorded in
`binary_population_summary.tsv` and are an important confounding diagnostic.


## Repository files

| File | Purpose |
|---|---|
| `S1_neutral_structure.slim` | Neutral-trait benchmark |
| `S2_trait_no_imposed_confounding.slim` | Designs A and B under Gaussian selection |
| `export_neutral_snps.py` | Tree-sequence simplification, recapitation, neutral-SNP overlay and genotype export |

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
