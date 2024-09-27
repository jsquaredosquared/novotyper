import marimo

__generated_with = "0.8.18"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def __(mo):
    mo.md(r"""# Report""")
    return


@app.cell(hide_code=True)
def __(mo):
    mo.md(
        r"""
        ## Background

        - Investigating the feasibility of an **alignment-based SV genotyper**.
        - The functionality is based on Novocraft's HLA algorithm.
        - A list of known structural variants is used to construct alternate contigs.
        - A new reference consisting of the **reference genome + SV alt contigs** is used for alignment.
        - If an SV from the list is present, more reads from this region are expected to align to the corresponding SV alt contig with a better mapping quality compared to the corresponding reference location.
        - Thus, **MAPQ** along with other information may be used to select the combination of alleles which best account for the observed reads.
        """
    )
    return


@app.cell(hide_code=True)
def __(mo):
    mo.md(
        r"""
        ## Assumptions

        - Let $R_{\text{MAPQ}} = \frac{\text{MAPQ}_{\text{alt}}}{\text{MAPQ}_{\text{ref}}}$.
        - If the sample is **homozygous alternate**, we expect $\text{MAPQ}_\text{alt} \gg \text{MAPQ}_{\text{ref}}$, and therefore the MAPQ ratio should be very large ($R_{\text{MAPQ}} \gg 1$).
        - If the sample is **homozygous reference**, we expect $\text{MAPQ}_\text{alt} \ll \text{MAPQ}_{\text{ref}}$, and therefore the MAPQ ratio should be very small ($R_{\text{MAPQ}} \ll 1$).
        - If the sample is **heterozygous**, we expect $\text{MAPQ}_\text{alt} \approx \text{MAPQ}_{\text{ref}}$ and therefore the MAPQ ratio should be close to 1 ($R_{\text{MAPQ}} \approx 1$).
        """
    )
    return


@app.cell(hide_code=True)
def __(mo):
    mo.md("""## Observations (in HG002)""")
    return


@app.cell(hide_code=True)
def __(mo):
    mo.md(r"""### MAPQ ratio distribution""")
    return


@app.cell
def __(overview_charts):
    overview_charts
    return


@app.cell(hide_code=True)
def __(mo):
    mo.md(
        r"""
        - Transformations (e.g., $\sqrt{R_{\text{MAPQ}}}$ or $\sqrt[3]{R_{\text{MAPQ}}}$) can make the distribution easier to visualize at a glance.
        - How will the distribution of true negatives (homozygous reference) affect genotyping?
            - If it does not overlap significantly with the peak, 0/0 should be resolvable and performance characteristics may not be drastically affected.
        """
    )
    return


@app.cell
def __(root_charts):
    root_charts
    return


@app.cell
def __(mo):
    mo.md(
        r"""- Zooming in on distirbution of $R_{\text{MAPQ}}$ to select reasonable cutoffs."""
    )
    return


@app.cell
def __(initial_charts):
    initial_charts
    return


@app.cell(hide_code=True)
def __(mo):
    mo.md(
        r"""
        - The $R_{\text{MAPQ}}$ distribution has a mode at around 1, as expected.
        - However, there is some overlap where SVs have an unexpected $R_{\text{MAPQ}}$ for the reported genotype.
        - Using the following cutoff values yields a recall of 78% (genotype called correctly):

        $$
        \text{NaN} = ./. \\
        [0,0.2) = 0/0 \\
        [0.2,2.8) = 0/1 \\
        [2.8,\infty) = 1/1
        $$
        """
    )
    return


@app.cell
def __(contingency_table, mo):
    mo.md(contingency_table.to_markdown())
    return


@app.cell
def __(correct_predictions_by_sizecat):
    correct_predictions_by_sizecat
    return


@app.cell
def __(correct_predictions_by_svtype):
    correct_predictions_by_svtype
    return


@app.cell(hide_code=True)
def __(mo):
    mo.md(
        r"""
        ## Things to consider

        ### Overlaps $R_{\text{MAPQ}}$ distributions

        - Ideally the distributions for 0/0, 0/1, and 1/1 should be separable to distinguish genotypes.
        - Overlaps mean **some incorrectly called genotypes may be inevitable**.
            - Also, MAPQ involves probability.
        - Is there a systematic cause that can be addressed?
            - MAPQ ratio close to 1 (implying 0/1) but actually 1/1.
            - MAPQ ratio very low (implying 0/0) but actually 0/1 or 1/1.
            - MAPQ ratio very high (implying 1/1) but actually 0/0 or 0/1.
        - Most of them tend to be in what look like repetitive or low-complexity regions.
        - [ ] TODO Is there a reason for the "spikes" in the MAPQ distributions near 0?

        ### Generalizability of cutoffs

        - Can the cutoffs be applied to other samples?
            - The distributions could be affected by things such as sequencing parameters (coverage, etc.).
        - To investigate cutoffs, the HGSVC2 benchmarks can be used to see if the distributions are similar, or to combine all samples and see the overall distribution.

        ### Precision

        - To assess precision, we need to include true negative SVs (homozygous reference).
            - Use the HGSVC2 benchmark VCF, which contains 0/0 variants.
            - Use simulated reads, so you are completely certain about the alleles.
            - Use the calls from an SV caller and see whether precision can be improved without negatively affecting recall.

        ### Other considerations

        - Split multiallelic VCF files prior to alt contig generation (so that the only possible genotypes are 0/0, 0/1, 1/1, or ./.).
        - Should there be a minimum MAPQ?
        - Should repeat and low-complexity regions be masked?
        - How will it perform on other SV types?
            - Duplications and inversions shouldn't be too difficult to implement.
            - Translocations?
        """
    )
    return


@app.cell(hide_code=True)
def __(mo):
    mo.md(
        r"""
        ## Comparison to other genotypers

        - The novoSV method has strong similarities to graph-based SV genotypers.
            - e.g., NPSV, SV2, Paragraph, vg, SVTyper, SVJedi-graph, GraphTyper, BayesTyper.
            - These use nodes of a graph instead of alt contigs.
            - They are usually much faster.
            - Some also struggle with nearby and overlapping SVs.
        - The novoSV method can take inspiration from some of their ideas:
            - Taking into consideration the uniqueness of read mapping and alignment counts for genotyping.
        - What would be the advantages of novoSV?
            - Can it be used to improve the precision of other SV callers by a significant amount?
     
        ## Limitations

        - This method is limited to trying to detect and genotype known SVs.
            - It cannot detect _de novo_ SVs.
            - A library of known SVs of interest (e.g., those associated with a phenotype or disease, or those previously called by another SV caller) will be useful for this method.
        - NovoAlign takes too long (especially slowing down development).
        """
    )
    return


@app.cell
def __(mo):
    mo.md(r"""# Setup""")
    return


@app.cell
def __():
    import marimo as mo

    return (mo,)


@app.cell
def __():
    import pandas as pd
    import pyranges as pr
    import polars as pl
    import bioframe as bf
    import numpy as np
    import altair as alt

    return alt, bf, np, pd, pl, pr


@app.cell
def __(alt):
    alt.data_transformers.enable("vegafusion")
    return


@app.cell(hide_code=True)
def __(mo):
    mo.md(r"""## Load data""")
    return


@app.cell(hide_code=True)
def __(pl):
    alt_to_ref = pl.read_csv(
        "../previous_outputs/HG002/alt_to_ref.txt",
        separator=" ",
        has_header=False,
        new_columns=["alt_contig", "ref_loc", "alt_len"],
    ).with_columns(
        pl.col("alt_contig").str.replace(">", ""),
        pl.col("ref_loc").str.replace("rg:", ""),
    )

    alt_to_ref.head()
    return (alt_to_ref,)


@app.cell
def __(mo):
    mo.md("""### Extract reference regions""")
    return


@app.cell
def __(alt_to_ref, bf, pd, pr):
    ref_locs = bf.from_any(alt_to_ref.to_series(1).to_list())
    ref_locs["alt_contig"] = pd.Series(alt_to_ref.to_series(0).to_list())
    ref_locs_pr = pr.PyRanges(
        ref_locs.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    )

    ref_locs_pr
    return ref_locs, ref_locs_pr


@app.cell(hide_code=True)
def __(mo):
    mo.md("""### Extract SV regions""")
    return


@app.cell
def __(alt_to_ref, pl, pr):
    sv_locs_pr = pr.PyRanges(
        alt_to_ref.select("alt_contig", "alt_len")
        .with_columns(
            pl.col("alt_contig")
            .str.extract_groups(r"(\d+|\w+)_(\d+)_(\d+)_\d")
            .struct.rename_fields(["Chromosome", "Start", "End"])
            .alias("alt_contig_locs"),
            pl.col("alt_len").str.extract(r"LN:(\d+)").str.to_integer(),
        )
        .unnest("alt_contig_locs")
        .with_columns(pl.col(["Start", "End"]).str.to_integer())
        .to_pandas()
    )

    sv_locs_pr
    return (sv_locs_pr,)


@app.cell(hide_code=True)
def __(mo):
    mo.md(r"""### Join regions""")
    return


@app.cell
def __(ref_locs_pr, sv_locs_pr):
    locs = sv_locs_pr.df.join(ref_locs_pr.df, rsuffix="_ref").drop(
        columns=["alt_contig_ref", "Chromosome_ref"]
    )

    locs
    return (locs,)


@app.cell
def __(mo):
    mo.md(r"""### Get subregions over which MAPQ will be tested""")
    return


@app.cell
def __(locs, np):
    test_locs_pr = locs.copy()

    # Get size of context to determine where on alt to start testing
    test_locs_pr["Start_alt_test"] = test_locs_pr["Start"] - test_locs_pr["Start_ref"]

    # Determine where to end testing
    test_locs_pr["End_alt_test"] = np.where(
        test_locs_pr["alt_len"] == 1801,
        test_locs_pr["Start"] + 1 - test_locs_pr["Start_ref"],
        test_locs_pr["Start_alt_test"] + (test_locs_pr["alt_len"] - 1801),
    )

    test_locs_pr
    return (test_locs_pr,)


@app.cell
def __(mo):
    mo.md(r"""### Read in VCF to filter for only PASS variants""")
    return


@app.cell
def __(pl):
    vcf = (
        pl.read_csv(
            "../resources/sv-vcf-files/HG002_SVs_Tier1_v0.6.ALL.vcf",
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
                "HG002",
            ],
            schema_overrides={"CHROM": pl.String},
        )
        .with_columns(
            pl.col("INFO")
            .str.extract_groups(r"sizecat=([\w\d]+);")
            .alias("SIZECAT")
            .struct.rename_fields(["SIZECAT"]),
            pl.col("INFO")
            .str.extract_groups(r"SVLEN=([-\d]+);")
            .cast(pl.Int32)
            .alias("SVLEN")
            .struct.rename_fields(["SVLEN"]),
            pl.col("INFO")
            .str.extract_groups(r"SVTYPE=(\w+);")
            .alias("SVTYPE")
            .struct.rename_fields(["SVTYPE"]),
            pl.col("HG002").str.extract(r"^([\d|.|\/]*):").alias("GT"),
        )
        .unnest(["SIZECAT", "SVLEN", "SVTYPE"])
    )

    vcf
    return (vcf,)


@app.cell
def __(vcf):
    vcf_filter_and_gt = vcf.select(
        ["FILTER", "GT", "SVLEN", "SIZECAT", "SVTYPE"]
    ).to_pandas()

    vcf_filter_and_gt
    return (vcf_filter_and_gt,)


@app.cell
def __(pd, test_locs_pr, vcf_filter_and_gt):
    test_locs_with_gt = pd.concat([test_locs_pr, vcf_filter_and_gt], axis=1)
    test_locs_with_gt.insert(1, "alt_start", [0] * test_locs_with_gt.shape[0])

    test_locs_with_gt = test_locs_with_gt.query('FILTER == "PASS"')

    test_locs_with_gt
    return (test_locs_with_gt,)


@app.cell
def __(mo):
    mo.md(r"""### Separate ref and alt test locs""")
    return


@app.cell
def __(test_locs_with_gt):
    ref_test_locs = test_locs_with_gt[["alt_contig", "Chromosome", "Start", "End"]]
    ref_test_locs
    return (ref_test_locs,)


@app.cell
def __(test_locs_with_gt):
    alt_test_locs = test_locs_with_gt.copy()

    alt_test_locs["Chromosome"] = alt_test_locs["alt_contig"]

    alt_test_locs = alt_test_locs[
        ["alt_contig", "Chromosome", "Start_alt_test", "End_alt_test"]
    ].rename(columns={"Start_alt_test": "Start", "End_alt_test": "End"})

    alt_test_locs
    return (alt_test_locs,)


@app.cell
def __(mo):
    mo.md(r"""## Prepare MAPQ""")
    return


@app.cell(hide_code=True)
def __(pd, pr):
    mapq_bedgraph = pd.read_csv(
        "../previous_outputs/HG002/HG002.sv.sorted.bedgraph",
        sep=" ",
        header=0,
        skiprows=1,
        names=["Chromosome", "Start", "End", "MAPQ"],
        dtype={
            "Chromosome": str,
        },
    )

    mapq_bedgraph_pr = pr.PyRanges(mapq_bedgraph)

    mapq_bedgraph_pr
    return mapq_bedgraph, mapq_bedgraph_pr


@app.cell
def __(mo):
    mo.md(r"""### Calculate MAPQ over REF test locs""")
    return


@app.cell
def __(mapq_bedgraph_pr, pr, ref_test_locs):
    ref_mapq = (
        pr.PyRanges(ref_test_locs)
        .extend(100)
        .join(mapq_bedgraph_pr)
        .df.groupby("alt_contig")
        .agg({"MAPQ": "mean"})
    )

    ref_mapq
    return (ref_mapq,)


@app.cell
def __(mo):
    mo.md(r"""### Calculate MAPQ over ALT test locs""")
    return


@app.cell
def __(alt_test_locs, mapq_bedgraph_pr, pr):
    alt_mapq = (
        pr.PyRanges(alt_test_locs)
        .extend(100)
        .join(mapq_bedgraph_pr)
        .df.groupby("alt_contig")
        .agg({"MAPQ": "mean"})
    )

    alt_mapq
    return (alt_mapq,)


@app.cell
def __(mo):
    mo.md(r"""## Join all results and calculate $R_{\text{MAPQ}}$""")
    return


@app.cell(hide_code=True)
def __(alt_mapq, np, ref_mapq, test_locs_with_gt):
    results = (
        test_locs_with_gt.set_index("alt_contig")
        .join(ref_mapq, rsuffix="_ref")
        .join(alt_mapq, rsuffix="_alt")
        .reset_index()
        .set_index(["Chromosome", "Start", "End"])
    )

    results["ratio_alt_to_ref"] = results["MAPQ_alt"] / results["MAPQ"]

    results["sqrt_ratio"] = np.sqrt(results["ratio_alt_to_ref"])
    results["cbrt_ratio"] = np.cbrt(results["ratio_alt_to_ref"])

    results
    return (results,)


@app.cell
def __(mo):
    mo.md(r"""## Plot distributions of MAPQ ratio""")
    return


@app.cell
def __(np):
    inf = np.inf
    return (inf,)


@app.cell
def __(alt, results):
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
    return (overview_chart,)


@app.cell
def __(alt, results):
    overview_chart_more_detail = (
        alt.Chart(
            results.query("not ratio_alt_to_ref.isna() and ratio_alt_to_ref != @inf")
        )
        .mark_line()
        .encode(
            alt.X("ratio_alt_to_ref:Q").bin(maxbins=1000),
            y="count()",
            color=alt.Color("GT"),
        )
        .interactive()
    )
    return (overview_chart_more_detail,)


@app.cell
def __(alt, results):
    overview_chart_even_more_detail = (
        alt.Chart(
            results.query("not ratio_alt_to_ref.isna() and ratio_alt_to_ref != @inf")
        )
        .mark_line()
        .encode(
            alt.X("ratio_alt_to_ref:Q").bin(maxbins=5000),
            y="count()",
            color=alt.Color("GT"),
        )
        .interactive()
    )
    return (overview_chart_even_more_detail,)


@app.cell
def __(alt, results):
    ratio_lt_5 = (
        alt.Chart(
            results.query(
                "not ratio_alt_to_ref.isna() and ratio_alt_to_ref != @inf and ratio_alt_to_ref < 5"
            )
        )
        .mark_line()
        .encode(
            alt.X("ratio_alt_to_ref:Q").bin(maxbins=100),
            y="count()",
            color=alt.Color("GT"),
        )
        .interactive()
    )
    return (ratio_lt_5,)


@app.cell
def __(alt, results):
    ratio_ge_5 = (
        alt.Chart(
            results.query(
                "not ratio_alt_to_ref.isna() and ratio_alt_to_ref != @inf and ratio_alt_to_ref > 5"
            )
        )
        .mark_line()
        .encode(
            alt.X("ratio_alt_to_ref:Q").bin(maxbins=100),
            y="count()",
            color=alt.Color("GT"),
        )
        .interactive()
    )
    return (ratio_ge_5,)


@app.cell
def __(alt, results):
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
    return (sqrt_chart,)


@app.cell
def __(alt, results):
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
    return (cbrt_chart,)


@app.cell
def __(
    cbrt_chart,
    overview_chart,
    overview_chart_even_more_detail,
    overview_chart_more_detail,
    ratio_ge_5,
    ratio_lt_5,
    sqrt_chart,
):
    overview_charts = (
        overview_chart | overview_chart_more_detail | overview_chart_even_more_detail
    )

    initial_charts = ratio_lt_5 | ratio_ge_5

    root_charts = sqrt_chart | cbrt_chart
    return initial_charts, overview_charts, root_charts


@app.cell(hide_code=True)
def __(mo):
    mo.md(
        r"""
        ## Make predictions based on $R_{\text{MAPQ}}$

        $$
        \text{NaN} = ./. \\
        [0,0.2) = 0/0 \\
        [0.2,2.8) = 0/1 \\
        [2.8,\infty) = 1/1 \\
        $$
        """
    )
    return


@app.cell
def __(np, results):
    predictions = results.copy()

    predictions["prediction"] = np.where(
        results["ratio_alt_to_ref"].isna(),
        "./.",
        np.where(
            results["ratio_alt_to_ref"] >= 2.8,
            "1/1",
            np.where(results["ratio_alt_to_ref"] < 0.2, "0/0", "0/1"),
        ),
    )

    predictions["is_correct"] = predictions["prediction"] == predictions["GT"]

    predictions
    return (predictions,)


@app.cell
def __(predictions):
    correct_predictions_by_sizecat = predictions.groupby("SIZECAT").agg(
        {"is_correct": ["sum", "count"]}
    )

    correct_predictions_by_sizecat["pct"] = (
        100
        * correct_predictions_by_sizecat[("is_correct", "sum")]
        / correct_predictions_by_sizecat[("is_correct", "count")]
    )

    correct_predictions_by_sizecat
    return (correct_predictions_by_sizecat,)


@app.cell
def __(predictions):
    correct_predictions_by_svtype = predictions.groupby("SVTYPE").agg(
        {"is_correct": ["sum", "count"]}
    )

    correct_predictions_by_svtype["pct"] = (
        100
        * correct_predictions_by_svtype[("is_correct", "sum")]
        / correct_predictions_by_svtype[("is_correct", "count")]
    )

    correct_predictions_by_svtype
    return (correct_predictions_by_svtype,)


@app.cell
def __(pd, predictions):
    contingency_table = pd.crosstab(
        predictions["GT"], predictions["prediction"], normalize="index"
    )
    contingency_table
    return (contingency_table,)


if __name__ == "__main__":
    app.run()
