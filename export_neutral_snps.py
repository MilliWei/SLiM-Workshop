#!/usr/bin/env python3
"""Overlay an exact number of neutral SNPs on a SLiM 5 trees archive.
python export_neutral_snps.py \
    S1_neutral_structure_10pop/final_trees \
    --output-dir S1_neutral_structure_10pop/neutral_snps \
    --sample-individuals-tsv S1_neutral_structure_10pop/individuals.tsv \
    --snps-per-chromosome 25000 \
    --recombination-rate 1e-8 \
    --qtl-dosage-tsv \
    S1_neutral_structure_10pop/qtl_plus_allele_dosage.tsv \
    --write-combined-dosage
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Optional

import msprime
import numpy as np
import pyslim
import tskit


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("trees_archive", type=Path)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--snps-per-chromosome",
        type=int,
        default=25_000,
        help=(
            "Final variant columns per chromosome. When QTL dosage and map "
            "are supplied, this total includes the QTLs."
        ),
    )
    parser.add_argument(
        "--individuals-per-population",
        type=int,
        default=100,
        help="0 keeps every final individual; default 100 controls matrix size.",
    )
    parser.add_argument(
        "--sex",
        choices=("all", "M", "F"),
        default="all",
        help="Restrict sampling to one sex; S0/S1 use M for male-only analysis.",
    )
    parser.add_argument(
        "--sample-individuals-tsv",
        type=Path,
        default=None,
        help=(
            "Use the individual_id and population columns in a SLiM-written "
            "individuals.tsv instead of drawing another random sample."
        ),
    )
    parser.add_argument("--ancestral-ne", type=float, default=4000.0)
    parser.add_argument(
        "--x-chromosome-index",
        type=int,
        default=0,
        help=(
            "One-based chromosome index for an X chromosome; 0 means none. "
            "Its recapitation Ne is set to 0.75 times --ancestral-ne."
        ),
    )
    parser.add_argument("--recombination-rate", type=float, default=2e-6)
    parser.add_argument("--seed", type=int, default=1024)
    parser.add_argument(
        "--write-dosage",
        action="store_true",
        help="Write individual x SNP uint8 .npy matrices (0/1/2).",
    )
    parser.add_argument(
        "--qtl-dosage-tsv",
        type=Path,
        default=None,
        help="SLiM-written causal-QTL dosage table to align by pedigree ID.",
    )
    parser.add_argument(
        "--qtl-map-tsv",
        type=Path,
        default=None,
        help=(
            "SLiM-written QTL map. Required with --write-combined-dosage so "
            "QTLs can be inserted at their genomic positions."
        ),
    )
    parser.add_argument(
        "--write-combined-dosage",
        action="store_true",
        help=(
            "Write one individual x (neutral SNP + causal QTL) uint8 matrix. "
            "Requires --qtl-dosage-tsv."
        ),
    )
    return parser.parse_args()


def population_label(ts: tskit.TreeSequence, individual_id: int) -> int:
    metadata = ts.individual(individual_id).metadata
    if isinstance(metadata, dict) and "subpopulation" in metadata:
        return int(metadata["subpopulation"])
    nodes = ts.individual(individual_id).nodes
    if len(nodes) == 0:
        raise ValueError(f"Individual {individual_id} has no nodes")
    pop_id = int(ts.node(nodes[0]).population)
    pop_metadata = ts.population(pop_id).metadata
    if isinstance(pop_metadata, dict) and "slim_id" in pop_metadata:
        return int(pop_metadata["slim_id"])
    return pop_id


def pedigree_id(ts: tskit.TreeSequence, individual_id: int) -> int:
    metadata = ts.individual(individual_id).metadata
    if isinstance(metadata, dict):
        for key in ("pedigree_id", "pedigreeID"):
            if key in metadata:
                return int(metadata[key])
    return int(individual_id)


def choose_individuals(
    ts: tskit.TreeSequence, per_population: int, seed: int, sex: str
) -> tuple[list[int], list[int]]:
    rng = np.random.default_rng(seed)
    alive = np.asarray(pyslim.individuals_alive_at(ts, 0), dtype=int)
    grouped: dict[int, list[int]] = {}
    requested_sex = {
        "M": pyslim.INDIVIDUAL_TYPE_MALE,
        "F": pyslim.INDIVIDUAL_TYPE_FEMALE,
    }.get(sex)
    for individual_id in alive:
        metadata = ts.individual(int(individual_id)).metadata
        if requested_sex is not None and metadata.get("sex") != requested_sex:
            continue
        pop = population_label(ts, int(individual_id))
        grouped.setdefault(pop, []).append(int(individual_id))

    chosen: list[int] = []
    labels: list[int] = []
    for pop in sorted(grouped):
        candidates = np.asarray(grouped[pop], dtype=int)
        if per_population > 0:
            if len(candidates) < per_population:
                raise ValueError(
                    f"Population {pop} has {len(candidates)} final individuals, "
                    f"fewer than requested {per_population}."
                )
            candidates = rng.choice(candidates, per_population, replace=False)
        candidates = np.sort(candidates)
        chosen.extend(candidates.tolist())
        labels.extend([pop] * len(candidates))
    return chosen, labels


def parse_population(value: str) -> int:
    """Accept either SLiM's numeric population ID or a label such as p10."""
    value = value.strip()
    if value.lower().startswith("p"):
        value = value[1:]
    return int(value)


def choose_individuals_from_tsv(
    ts: tskit.TreeSequence, path: Path
) -> tuple[list[int], list[int]]:
    """Match a SLiM sample manifest to alive tree-sequence individuals."""
    alive = np.asarray(pyslim.individuals_alive_at(ts, 0), dtype=int)
    by_pedigree: dict[int, int] = {}
    for individual_id in alive:
        pid = pedigree_id(ts, int(individual_id))
        if pid in by_pedigree:
            raise ValueError(f"Duplicate alive pedigree ID {pid} in tree sequence")
        by_pedigree[pid] = int(individual_id)

    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"individual_id", "population"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"{path} must contain individual_id and population columns"
            )
        requested = [
            (int(row["individual_id"]), parse_population(row["population"]))
            for row in reader
        ]

    requested_ids = [pid for pid, _ in requested]
    if len(requested_ids) != len(set(requested_ids)):
        raise ValueError(f"Duplicate individual_id in {path}")

    chosen: list[int] = []
    labels: list[int] = []
    for pid, expected_population in requested:
        if pid not in by_pedigree:
            raise ValueError(f"Pedigree ID {pid} from {path} is not alive in the trees")
        individual_id = by_pedigree[pid]
        observed_population = population_label(ts, individual_id)
        if observed_population != expected_population:
            raise ValueError(
                f"Population mismatch for pedigree ID {pid}: "
                f"{expected_population} != {observed_population}"
            )
        chosen.append(individual_id)
        labels.append(observed_population)
    if len(chosen) == 0:
        raise ValueError(f"No sampled individuals found in {path}")
    return chosen, labels


def sample_nodes(
    ts: tskit.TreeSequence, individuals: list[int]
) -> tuple[list[int], list[int]]:
    nodes: list[int] = []
    ploidies: list[int] = []
    for individual_id in individuals:
        ind_nodes = [
            int(node_id)
            for node_id in ts.individual(individual_id).nodes
            if not pyslim.node_is_vacant(ts, ts.node(int(node_id)))
        ]
        if len(ind_nodes) not in (1, 2):
            raise ValueError(
                f"Expected one or two sampled nodes for individual {individual_id}; "
                f"found {len(ind_nodes)}."
            )
        nodes.extend(ind_nodes)
        ploidies.append(len(ind_nodes))
    return nodes, ploidies


def mutation_opportunity(ts: tskit.TreeSequence) -> float:
    return sum(tree.total_branch_length * tree.span for tree in ts.trees())


def overlay_exact_sites(
    ts: tskit.TreeSequence, target: int, seed: int
) -> tskit.TreeSequence:
    opportunity = mutation_opportunity(ts)
    if opportunity <= 0:
        raise ValueError("Tree sequence has no branch length for mutation overlay")

    rate = 1.05 * target / opportunity
    for attempt in range(20):
        mutated = msprime.sim_mutations(
            ts,
            rate=rate,
            random_seed=seed + attempt,
            model=msprime.BinaryMutationModel(),
            discrete_genome=False,
            keep=False,
        )
        if mutated.num_sites >= target:
            break
        rate *= max(1.10, 1.05 * target / max(mutated.num_sites, 1))
    else:
        raise RuntimeError("Could not generate the requested number of neutral SNPs")

    if mutated.num_sites > target:
        rng = np.random.default_rng(seed + 10_000)
        keep = np.sort(rng.choice(mutated.num_sites, target, replace=False))
        remove = np.setdiff1d(np.arange(mutated.num_sites), keep, assume_unique=True)
        mutated = mutated.delete_sites(remove)
    if mutated.num_sites != target:
        raise AssertionError((mutated.num_sites, target))
    return mutated


def fill_dosage(
    ts: tskit.TreeSequence,
    dosage: np.ndarray,
    column_offset: int,
    ploidies: list[int],
) -> None:
    samples = ts.samples()
    if len(samples) != sum(ploidies):
        raise ValueError("Simplified tree sequence has an unexpected sample count")
    for site_index, variant in enumerate(ts.variants(samples=samples)):
        calls = np.asarray(variant.genotypes, dtype=np.int16)
        if np.any(calls < 0):
            raise ValueError(f"Missing genotype at site {site_index}")
        offset = 0
        for row, ploidy in enumerate(ploidies):
            dosage[row, column_offset + site_index] = calls[
                offset : offset + ploidy
            ].sum(dtype=np.uint8)
            offset += ploidy


def fill_dosage_columns(
    ts: tskit.TreeSequence,
    dosage: np.ndarray,
    columns: list[int],
    ploidies: list[int],
) -> None:
    """Write site dosages into arbitrary output columns."""
    if len(columns) != ts.num_sites:
        raise ValueError("One output column is required for every tree-sequence site")
    samples = ts.samples()
    if len(samples) != sum(ploidies):
        raise ValueError("Simplified tree sequence has an unexpected sample count")
    for site_index, variant in enumerate(ts.variants(samples=samples)):
        calls = np.asarray(variant.genotypes, dtype=np.int16)
        if np.any(calls < 0):
            raise ValueError(f"Missing genotype at site {site_index}")
        offset = 0
        column = columns[site_index]
        for row, ploidy in enumerate(ploidies):
            dosage[row, column] = calls[offset : offset + ploidy].sum(
                dtype=np.uint8
            )
            offset += ploidy


def write_dosage(
    ts: tskit.TreeSequence, path: Path, ploidies: list[int]
) -> None:
    dosage = np.lib.format.open_memmap(
        path,
        mode="w+",
        dtype=np.uint8,
        shape=(len(ploidies), ts.num_sites),
    )
    fill_dosage(ts, dosage, 0, ploidies)
    dosage.flush()


def load_qtl_dosage(
    path: Path, pedigree_ids: list[int], populations: list[int]
) -> tuple[list[str], np.ndarray]:
    with path.open(newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        if header[:2] != ["individual_id", "population"]:
            raise ValueError(f"Unexpected QTL dosage header in {path}")
        qtl_ids = header[2:]
        rows: dict[int, tuple[int, list[int]]] = {}
        for fields in reader:
            pedigree_id = int(fields[0])
            if pedigree_id in rows:
                raise ValueError(f"Duplicate pedigree ID {pedigree_id} in {path}")
            rows[pedigree_id] = (int(fields[1]), [int(value) for value in fields[2:]])

    dosage = np.empty((len(pedigree_ids), len(qtl_ids)), dtype=np.uint8)
    for row_index, (pedigree_id, expected_population) in enumerate(
        zip(pedigree_ids, populations)
    ):
        if pedigree_id not in rows:
            raise ValueError(f"Pedigree ID {pedigree_id} is missing from {path}")
        population, values = rows[pedigree_id]
        if population != expected_population:
            raise ValueError(
                f"Population mismatch for pedigree ID {pedigree_id}: "
                f"{population} != {expected_population}"
            )
        if len(values) != len(qtl_ids):
            raise ValueError(f"Wrong QTL count for pedigree ID {pedigree_id}")
        dosage[row_index] = values
    return qtl_ids, dosage


def load_qtl_map(
    path: Path, qtl_ids: list[str], chromosome_count: int
) -> list[dict[str, object]]:
    """Load and validate QTL chromosome positions in dosage-column order."""
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"qtl", "chromosome", "position"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain qtl, chromosome, and position")
        by_id: dict[str, dict[str, object]] = {}
        for row in reader:
            qtl_id = row["qtl"]
            if qtl_id in by_id:
                raise ValueError(f"Duplicate QTL ID {qtl_id} in {path}")
            chromosome = int(row["chromosome"])
            if chromosome < 1 or chromosome > chromosome_count:
                raise ValueError(
                    f"QTL {qtl_id} chromosome {chromosome} is outside 1.."
                    f"{chromosome_count}"
                )
            by_id[qtl_id] = {
                "qtl": qtl_id,
                "chromosome": chromosome,
                "position": float(row["position"]),
            }

    missing = [qtl_id for qtl_id in qtl_ids if qtl_id not in by_id]
    extra = [qtl_id for qtl_id in by_id if qtl_id not in set(qtl_ids)]
    if missing or extra:
        raise ValueError(
            f"QTL map/dosage mismatch; missing={missing}, extra={extra}"
        )

    records = []
    for qtl_index, qtl_id in enumerate(qtl_ids):
        record = dict(by_id[qtl_id])
        record["qtl_index"] = qtl_index
        records.append(record)
    return records


def main() -> None:
    args = parse_args()
    archive = args.trees_archive.resolve()
    files = sorted(archive.glob("chromosome_*.trees"))
    if len(files) == 0:
        raise ValueError(f"No chromosome_*.trees files found in {archive}")
    if args.x_chromosome_index < 0 or args.x_chromosome_index > len(files):
        raise ValueError("--x-chromosome-index is outside the chromosome range")
    if args.write_combined_dosage and (
        args.qtl_dosage_tsv is None or args.qtl_map_tsv is None
    ):
        raise ValueError(
            "--write-combined-dosage requires both --qtl-dosage-tsv and "
            "--qtl-map-tsv"
        )

    output = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else archive.parent / "neutral_snps"
    )
    output.mkdir(parents=True, exist_ok=True)

    first = tskit.load(str(files[0]))
    if args.sample_individuals_tsv is not None:
        individuals, populations = choose_individuals_from_tsv(
            first, args.sample_individuals_tsv.resolve()
        )
    else:
        individuals, populations = choose_individuals(
            first, args.individuals_per_population, args.seed, args.sex
        )
    pedigree_ids = [pedigree_id(first, individual_id) for individual_id in individuals]

    qtl_ids: list[str] = []
    qtl_dosage: Optional[np.ndarray] = None
    if args.qtl_dosage_tsv is not None:
        qtl_ids, qtl_dosage = load_qtl_dosage(
            args.qtl_dosage_tsv.resolve(), pedigree_ids, populations
        )
        np.save(output / "causal_qtl_dosage.npy", qtl_dosage)
        with (output / "causal_qtl_ids.tsv").open("w", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t")
            writer.writerow(["qtl_dosage_column", "qtl"])
            writer.writerows(enumerate(qtl_ids))

    qtl_records: list[dict[str, object]] = []
    if args.qtl_map_tsv is not None:
        if qtl_dosage is None:
            raise ValueError("--qtl-map-tsv also requires --qtl-dosage-tsv")
        qtl_records = load_qtl_map(
            args.qtl_map_tsv.resolve(), qtl_ids, len(files)
        )

    combined_dosage = None
    if args.write_combined_dosage:
        assert qtl_dosage is not None
        combined_dosage = np.lib.format.open_memmap(
            output / "combined_neutral_and_qtl_dosage.npy",
            mode="w+",
            dtype=np.uint8,
            shape=(
                len(individuals),
                len(files) * args.snps_per_chromosome,
            ),
        )

    with (output / "sampled_individuals.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["row", "tree_individual_id", "pedigree_id", "population"])
        for row, (ind, pedigree, pop) in enumerate(
            zip(individuals, pedigree_ids, populations)
        ):
            writer.writerow([row, ind, pedigree, pop])

    manifest_rows = []
    feature_map_rows: list[tuple[object, ...]] = []
    for chromosome_index, tree_path in enumerate(files, start=1):
        ts = tskit.load(str(tree_path))
        nodes, ploidies = sample_nodes(ts, individuals)
        sampled = ts.simplify(
            samples=nodes,
            filter_individuals=True,
            keep_input_roots=True,
        )
        recapitation_ne = args.ancestral_ne
        if chromosome_index == args.x_chromosome_index:
            recapitation_ne *= 0.75
        recapitated = pyslim.recapitate(
            sampled,
            ancestral_Ne=recapitation_ne,
            recombination_rate=args.recombination_rate,
            random_seed=args.seed + chromosome_index * 100,
        )
        chromosome_qtls = [
            record
            for record in qtl_records
            if record["chromosome"] == chromosome_index
        ]
        neutral_target = args.snps_per_chromosome
        if combined_dosage is not None:
            neutral_target -= len(chromosome_qtls)
            if neutral_target < 0:
                raise ValueError(
                    f"Chromosome {chromosome_index} has more QTLs than the "
                    "requested final SNP count"
                )

        neutral = overlay_exact_sites(
            recapitated,
            neutral_target,
            args.seed + chromosome_index * 1000,
        )
        neutral_path = output / f"chromosome_{chromosome_index}_neutral.trees"
        neutral.dump(neutral_path)
        if args.write_dosage:
            write_dosage(
                neutral,
                output / f"chromosome_{chromosome_index}_dosage.npy",
                ploidies,
            )
        if combined_dosage is not None:
            assert qtl_dosage is not None
            chromosome_start = (
                chromosome_index - 1
            ) * args.snps_per_chromosome
            features: list[tuple[float, str, int, str]] = []
            for neutral_index, site in enumerate(neutral.sites()):
                features.append(
                    (
                        float(site.position),
                        "neutral_snp",
                        neutral_index,
                        f"chr{chromosome_index}_neutral_{neutral_index + 1}",
                    )
                )
            for record in chromosome_qtls:
                features.append(
                    (
                        float(record["position"]),
                        "causal_qtl",
                        int(record["qtl_index"]),
                        str(record["qtl"]),
                    )
                )
            features.sort(key=lambda feature: (feature[0], feature[1]))
            if len(features) != args.snps_per_chromosome:
                raise AssertionError(
                    (chromosome_index, len(features), args.snps_per_chromosome)
                )

            neutral_columns = [-1] * neutral.num_sites
            for local_column, (position, feature_type, source_index, feature_id) in enumerate(
                features
            ):
                combined_column = chromosome_start + local_column
                feature_map_rows.append(
                    (
                        combined_column,
                        chromosome_index,
                        position,
                        feature_type,
                        feature_id,
                    )
                )
                if feature_type == "neutral_snp":
                    neutral_columns[source_index] = combined_column
                else:
                    combined_dosage[:, combined_column] = qtl_dosage[
                        :, source_index
                    ]
            if any(column < 0 for column in neutral_columns):
                raise AssertionError("A neutral SNP was not assigned an output column")
            fill_dosage_columns(
                neutral,
                combined_dosage,
                neutral_columns,
                ploidies,
            )
        manifest_rows.append(
            (chromosome_index, neutral.num_sites, len(individuals), neutral_path.name)
        )

    if combined_dosage is not None:
        combined_dosage.flush()
        with (output / "combined_feature_blocks.tsv").open("w", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t")
            writer.writerow(["start_column", "end_column_exclusive", "feature_type", "source"])
            for chromosome_index in range(1, len(files) + 1):
                start = (chromosome_index - 1) * args.snps_per_chromosome
                writer.writerow(
                    [
                        start,
                        start + args.snps_per_chromosome,
                        "genomic_order",
                        f"chromosome_{chromosome_index}",
                    ]
                )
        with (output / "combined_feature_map.tsv").open("w", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t")
            writer.writerow(
                ["column", "chromosome", "position", "feature_type", "feature_id"]
            )
            writer.writerows(feature_map_rows)

    with (output / "manifest.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["chromosome", "neutral_snps", "individuals", "tree_file"])
        writer.writerows(manifest_rows)


if __name__ == "__main__":
    main()
