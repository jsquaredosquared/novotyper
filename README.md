# Novotyper

This repository contains the code used during the development of a prototype structural variant genotyper tentatively called Novotyper.

## Usage

1. Obtain a VCF of the structural variants you would like to genotype.

2. Run the Snakemake workflow found in `workflow/Snakefile`.
    - Modify the `config/config.yaml` file to point to the correct locations for each input file and executable.
    - The workflow should perform the following steps (but you could run them yourself):
        - `novoSV vcf2altscaffold`: Use the SVs from the VCF file and the reference `.nix` file to generate a FASTA file containing the alt scaffolds.
        - `novoindex`: Index the reference fasta + alt scaffolds fasta to create a new spiked `.nix` file.
        - `novoalign`: Produce the BAM file by aligning the reads using the spiked `.nix` file.
        - `novoSV zabedgraph`: Use the alt scaffolds FASTA and the BAM file to obtain the MAPQ as a bedgraph.
        - `novoutil bgzf` [OPTIONAL]: Combine reference fasta + alt scaffolds fasta into a new fasta file, bgzip it, and index it with `samtools faidx`. Do this if you would like to visualize the alignments in IGV.

3. Run `novotyper` to predict genotypes and perform some benchmarking.
    - This step requires your SV VCF, plus the alt scaffolds FASTA and MAPQ bedgraph generated in step 2 above.
    - Clone this repository, and run `pip install .` in the virtual environment in which you wish to install `novotyper`.
      - If successfully installed, typing `novotyper genotype --help` should bring up the list of parameters required to run the program.

## Report

- For a basic description of the reasoning behind each step and the findings made during prototyping, please see `notebooks/report.md`.
- If you want to know how `novotyper` works step by step, check out `notebooks/novotyper.ipynb`. Each cell corresponds to a function which you can run to see its output.
- The functions are defined in `src/novotyper/genotyper.py` and are run in the following order:

    ```coconut
    vcf_info = gt.extract_info_from_vcf(vcf_path, sample)
    mapq_bedgraph = gt.read_mapq_bedgraph(bedgraph_file)

    test_locs = (
        fasta_file
        |> gt.read_alts_fasta_descriptions
        |> lift(,)(gt.extract_ref_locs_of_alts, gt.extract_sv_locs)
        |*> gt.join_ref_and_sv_locs
        |> gt.add_vcf_info_to_test_locs$(?, vcf_info)
        |> gt.get_pass_variants
    )

    ratio_results = (
        test_locs
        |> lift(,)(gt.get_ref_test_locs, gt.get_alt_test_locs)
        |> map$(gt.calculate_mapq$(?, mapq_bedgraph))
        |*> gt.join_and_calculate_mapq_ratio$(test_locs)
        |> gt.predict_genotype$(?, 0.2, 2.8)
    )

    gt.calculate_performance(ratio_results, f"{out_dir}/performance.md")
    gt.plot_mapq_distribution(ratio_results, out_dir)
    ```
