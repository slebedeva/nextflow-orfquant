This is a nextflow wrapper for [ORFquant](https://github.com/lcalviell/ORFquant). This is a tool to call translated regions from ribosome profiling data.

### Prerequisites:

- install Docker or apptainer (no conda support)
- install [nextflow](https://www.nextflow.io/)

### How to run

To run ORFquant analysis, you first need the output of the RiboseQC pipeline. You can run the original pipeline or the [nextflow wrapper](https://github.com/slebedeva/nextflow-riboseqc).

Minimal input parameters:
- `input_dir` : directory with results of RiboseQC (minimal required file: `*_for_ORFquant`)
- `gtf` : unzipped gtf file of your annotation 
- `fasta` : unzipped genome sequence of your annotation

You can use either `local` profile to run on your computer or `slurm` if you have access to a cluster with SLURM scheduler.


For example, to run locally:

```bash
nextflow run slebedeva/nextflow-orfquant \
-profile local,docker -resume \
--input_dir "$INPUT_DIR" \
--gtf $GTF \
--fasta $FASTA 
```

If you have already generated ORFquant annotation , you can provide it with `--rannot` argument.