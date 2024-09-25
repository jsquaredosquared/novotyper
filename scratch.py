from cyclopts import App
import genotyper

novotyper = App()


@novotyper.command
def genotype(
    alts_fasta,
    sv_vcf,
    mapq_bedgraph,
):
    """Predict the genotype based on the MAPQ ratio."""
    pass


if __name__ == "__main__":
    print("Hello from novotyper.py")
    genotyper.say_hello()
