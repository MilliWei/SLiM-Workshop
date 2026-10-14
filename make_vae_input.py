from pathlib import Path

import numpy as np
import pandas as pd

X = np.load("S1_neutral_structure_10pop/neutral_snps/combined_neutral_and_qtl_dosage.npy",
    mmap_mode="r",
)

samples = pd.read_csv("S1_neutral_structure_10pop/neutral_snps/sampled_individuals.tsv",
    sep="\t",
).sort_values("row")

phenotypes = pd.read_csv("S1_neutral_structure_10pop/individuals.tsv",
    sep="\t",
).set_index("individual_id")

# 按genotype矩阵的行顺序重新排列phenotype
aligned = phenotypes.loc[samples["pedigree_id"].to_numpy()]

y = aligned["phenotype_binary"].to_numpy(
    dtype=np.int64
)
population = aligned["population"].to_numpy()
pedigree_id = samples["pedigree_id"].to_numpy()


assert X.shape[0] == len(samples)
assert len(y) == X.shape[0]
assert np.array_equal(
    aligned.index.to_numpy(),
    pedigree_id,
)

print("genotype:", X.shape)
print("phenotype:", y.shape)
print("population:", population.shape)

assert X.shape[0] == len(samples)
assert len(y) == X.shape[0]
assert np.array_equal(
    aligned.index.to_numpy(),
    pedigree_id,
)

print("genotype:", X.shape)
print("phenotype:", y.shape)
print("population:", population.shape)

assert X.shape[0] == y.shape[0] == population.shape[0]
assert X.shape[1] == 50_000
assert np.isin(X, [0, 1, 2]).all()
assert np.isfinite(y).all()
assert np.unique(population).size == 10

import numpy as np

print("X shape:", X.shape)
print("y shape:", y.shape)

print("\ntop 30 individual × top 30 SNP：")
print(X[:30, :30])

print("\ntop 30 individual phenotype：")
print(y[:30])

summary = pd.crosstab(
    phenotypes["population"],
    phenotypes["phenotype_binary"],
    margins=False,
)

summary = summary.reindex(
    columns=[0, 1],
    fill_value=0,
)

proportions = summary.div(
    summary.sum(axis=1),
    axis=0,
)

print("Counts:")
print(summary)

print("\nProportions:")
print(proportions)

print("\nPercentages:")
print(proportions * 100)