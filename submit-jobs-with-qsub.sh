snakemake \
--sdm conda \
--executor cluster-generic \
--cluster-generic-submit-cmd "qsub -V -sync n -cwd -pe slave_pe {threads} -q slave.q -j y -N {rule}" \
--jobs unlimited \
--rerun-incomplete