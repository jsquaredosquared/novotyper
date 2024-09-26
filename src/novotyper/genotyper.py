import re
from pathlib import Path

import pandas as pd
import pyranges as pr
import polars as pl
import bioframe as bf
import numpy as np
import altair as alt


alt.data_transformers.enable("vegafusion")


def read_alts_fasta_descriptions(alts_fasta: Path) -> pl.DataFrame:
    with open(alts_fasta, "r") as file:
        records = re.findall(r">([\d\w_]+) LN:(\d+) .+ rg:([\d\w:-]+)", file.read())

        fields: dict[str, list] = {
            "alt_contig": [],
            "ref_loc": [],
            "alt_contig_len": [],
        }

        for record in records:
            fields["alt_contig"].append(record[0])
            fields["alt_contig_len"].append(record[1])
            fields["ref_loc"].append(record[2])

        data = pl.DataFrame(fields)

    return data


def extract_ref_locs_of_alts(data: pl.DataFrame) -> pr.PyRanges:
    ref_locs_bf = bf.from_any(data.to_series(1).to_list())
    ref_locs_bf["alt_contig"] = pd.Series(data.to_series(0).to_list())
    ref_locs_of_alts = pr.PyRanges(
        ref_locs_bf.rename(
            columns={"chrom": "Chromosome", "start": "Start", "end": "End"}
        )
    )

    return ref_locs_of_alts


def extract_sv_locs(data: pl.DataFrame) -> pr.PyRanges:
    sv_locs = pr.PyRanges(
        data.select("alt_contig", "alt_contig_len")
        .with_columns(
            pl.col("alt_contig")
            .str.extract_groups(r"(\d+|\w+)_(\d+)_(\d+)_\d")
            .struct.rename_fields(["Chromosome", "Start", "End"])
            .alias("alt_contig_locs"),
            pl.col("alt_contig_len").str.to_integer(),
        )
        .unnest("alt_contig_locs")
        .with_columns(pl.col(["Start", "End"]).str.to_integer())
        .to_pandas()
    )

    return sv_locs


def join_ref_and_sv_locs(
    ref_locs_of_alts: pr.PyRanges, sv_locs: pr.PyRanges
) -> pd.DataFrame:
    locs = sv_locs.df.join(ref_locs_of_alts.df, rsuffix="_ref").drop(
        columns=["alt_contig_ref", "Chromosome_ref"]
    )

    return locs


def define_alt_test_locs(locs: pd.DataFrame) -> pd.DataFrame:
    test_locs = locs.copy()

    # Get size of context to determine where on alt to start testing
    test_locs["Start_alt_test"] = test_locs["Start"] - test_locs["Start_ref"]

    # Determine where to end testing
    test_locs["End_alt_test"] = np.where(
        test_locs["alt_contig_len"] == 1801,
        test_locs["Start"] + 1 - test_locs["Start_ref"],
        test_locs["Start_alt_test"] + test_locs["alt_contig_len"] - 1801,
    )

    return test_locs


def extract_info_from_vcf(vcf_file: Path, sample: str) -> pd.DataFrame:
    vcf = (
        pl.read_csv(
            vcf_file,
            separator="\t",
            comment_prefix="#",
            has_header=False,
            new_columns=[
                "CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
                sample,
            ],
            schema_overrides={"CHROM": pl.String},
        )
        .with_columns(
            pl.col("INFO")
            .str.extract_groups(r"SVLEN=([-\d]+);")
            .cast(pl.Int32)
            .alias("SVLEN")
            .struct.rename_fields(["SVLEN"]),
            pl.col("INFO")
            .str.extract_groups(r"SVTYPE=(\w+);")
            .alias("SVTYPE")
            .struct.rename_fields(["SVTYPE"]),
            pl.col(sample).str.extract(r"^([\d.\/\|]*):").alias("GT"),
        )
        .unnest(["SVLEN", "SVTYPE"])
        .select(["FILTER", "GT", "SVLEN", "SVTYPE"])
        .to_pandas()
    )

    return vcf


def add_vcf_info_to_test_locs(locs: pd.DataFrame, vcf: pd.DataFrame) -> pd.DataFrame:
    locs_with_info = pd.concat([locs, vcf], axis=1)
    locs_with_info.insert(1, "alt_start", [0] * locs_with_info.shape[0])
    return locs_with_info


# TODO: Change this to accept a specified range of filters and genotypes.
def get_pass_variants(locs: pd.DataFrame) -> pd.DataFrame:
    return locs.query('FILTER == "PASS" or FILTER == "."')


def get_ref_test_locs(test_locs: pd.DataFrame) -> pd.DataFrame:
    return test_locs.copy()[["alt_contig", "Chromosome", "Start", "End"]]


def get_alt_test_locs(test_locs: pd.DataFrame) -> pd.DataFrame:
    alt_test_locs = test_locs.copy()
    alt_test_locs["Chromosome"] = alt_test_locs["alt_contig"]
    alt_test_locs = alt_test_locs.copy()[
        ["alt_contig", "Chromosome", "Start_alt_test", "End_alt_test"]
    ].rename(columns={"Start_alt_test": "Start", "End_alt_test": "End"})

    return alt_test_locs


def read_mapq_bedgraph(bedgraph_file: Path) -> pr.PyRanges:
    bedgraph: pd.DataFrame = pd.read_csv(
        bedgraph_file,
        sep=" ",
        header=0,
        skiprows=1,
        names=["Chromosome", "Start", "End", "MAPQ"],
        dtype={"Chromosome": str},
    )

    return pr.PyRanges(bedgraph)


def calculate_mapq(test_locs: pd.DataFrame, mapq_bedgraph: pd.DataFrame):
    return (
        pr.PyRanges(test_locs)
        .extend(100)
        .join(mapq_bedgraph)
        .df.groupby("alt_contig")
        .agg({"MAPQ": "mean"})
    )


def join_and_calculate_mapq_ratio(
    locs_with_info: pd.DataFrame, ref_mapq: pd.DataFrame, alt_mapq: pd.DataFrame
):
    results = (
        locs_with_info.set_index("alt_contig")
        .join(ref_mapq, rsuffix="_ref")
        .join(alt_mapq, rsuffix="_alt")
        .reset_index()
        .set_index(["Chromosome", "Start", "End"])
    )

    results["ratio_alt_to_ref"] = results["MAPQ_alt"] / results["MAPQ"]

    results["sqrt_ratio"] = np.sqrt(results["ratio_alt_to_ref"])
    results["cbrt_ratio"] = np.cbrt(results["ratio_alt_to_ref"])

    return results


def predict_genotype(
    results: pd.DataFrame, lower_het_bound: float, upper_het_bound: float
) -> pd.DataFrame:
    predictions = results.copy()

    predictions["prediction"] = np.where(
        results["ratio_alt_to_ref"].isna(),
        "./.",
        np.where(
            results["ratio_alt_to_ref"] >= upper_het_bound,
            "1/1",
            np.where(results["ratio_alt_to_ref"] < lower_het_bound, "0/0", "0/1"),
        ),
    )

    predictions["GT_unphased"] = predictions["GT"].map(
        {"0|0": "0/0", "0|1": "0/1", "1|0": "0/1", "1|1": "1/1"}
    )
    
    predictions["GT_concordance"] = predictions["prediction"] == predictions["GT_unphased"]

    return predictions


def calculate_performance(predictions: pd.DataFrame, out_file: str) -> None:
    recall = 100 * sum(predictions["GT_concordance"]) / predictions.shape[0]

    correct_predictions_by_svtype = predictions.groupby("SVTYPE").agg(
        {"GT_concordance": ["sum", "count"]}
    )
    correct_predictions_by_svtype["pct"] = (
        100
        * correct_predictions_by_svtype[("GT_concordance", "sum")]
        / correct_predictions_by_svtype[("GT_concordance", "count")]
    )

    contingency_table = pd.crosstab(
        predictions["GT_unphased"], predictions["prediction"], normalize="index"
    )

    with open(out_file, "w") as file:
        print(
            f"recall: {recall}",
            contingency_table.to_markdown(),
            correct_predictions_by_svtype.to_markdown(),
            sep="\n",
            file=file,
        )


def plot_mapq_distribution(results: pd.DataFrame, out_dir: str) -> None:
    inf = np.inf

    overview_chart = (
        alt.Chart(
            results.query("not ratio_alt_to_ref.isna() and ratio_alt_to_ref != @inf")
        )
        .mark_line()
        .encode(
            alt.X("ratio_alt_to_ref:Q").bin(maxbins=100),
            y="count()",
            color=alt.Color("GT"),
        )
        .interactive()
    )

    sqrt_chart = (
        alt.Chart(
            results.query("not ratio_alt_to_ref.isna() and ratio_alt_to_ref != @inf")
        )
        .mark_line()
        .encode(
            alt.X("sqrt_ratio:Q").bin(maxbins=250), y="count()", color=alt.Color("GT")
        )
        .interactive()
    )

    cbrt_chart = (
        alt.Chart(
            results.query("not ratio_alt_to_ref.isna() and ratio_alt_to_ref != @inf")
        )
        .mark_line()
        .encode(
            alt.X("cbrt_ratio:Q").bin(maxbins=250), y="count()", color=alt.Color("GT")
        )
        .interactive()
    )

    overview_chart.save(f"{out_dir}/overview_chart.svg")
    sqrt_chart.save(f"{out_dir}/sqrt_chart.svg")
    cbrt_chart.save(f"{out_dir}/cbrt_chart.svg")


def say_hello():
    print("Hello from genotyper!")
