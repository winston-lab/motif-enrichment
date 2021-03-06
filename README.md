
# motif enrichment pipeline

## description

A pipeline for testing enrichment and depletion of a database of motifs in 'condition' regions versus control regions.

- search genome for motif occurrences with [FIMO](http://meme-suite.org/doc/fimo.html) 
- find motif occurrences overlapping condition and control regions
- test motif enrichment using Fisher's exact test

## requirements

### required software

- Unix-like operating system (tested on CentOS 7.2.1511)
- Git
- [conda](https://conda.io/docs/user-guide/install/index.html)

### required files

- a FASTA file of the genome

- motif databases in [MEME](http://meme-suite.org/doc/meme-format.html) format

- [BED6](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format annotation files of the regions to search for motifs

## instructions
**0**. Clone this repository.

```bash
git clone https://github.com/winston-lab/motif-enrichment.git
```

**1**. Create and activate the `snakemake_default` virtual environment for the pipeline using conda. The virtual environment creation can take a while. If you've already created the `snakemake_default` environment from another one of my pipelines, this is the same environment, so you can skip creating the environment and just activate it.

```bash
# navigate into the pipeline directory
cd motif-enrichment

# create the snakemake_default environment
conda env create -v -f envs/snakemake_default.yaml

# activate the environment
source activate snakemake_default

# to deactivate the environment
# source deactivate
```

**2**. Make a copy of the configuration file template `config_template.yaml` called `config.yaml`, and edit `config.yaml` to suit your needs.

```bash
# make a copy of the configuration template file
cp config_template.yaml config.yaml

# edit the configuration file
vim config.yaml    # or use your favorite editor
```

**3**. With the `snakemake_default` environment activated, do a dry run of the pipeline to see what files will be created.

```bash
snakemake -p --use-conda --dryrun
```

**4**. If running the pipeline on a local machine, you can run the pipeline using the above command, omitting the `--dryrun` flag. You can also use N cores by specifying the `--cores N` flag. The first time the pipeline is run, conda will create separate virtual environments for some of the jobs to operate in. Running the pipeline on a local machine can take a long time, especially for many samples, so it's recommended to use an HPC cluster if possible. On the HMS O2 cluster, which uses the SLURM job scheduler, entering `sbatch slurm_submit.sh` will submit the pipeline as a single job which spawns individual subjobs as necessary. This can be adapted to other job schedulers and clusters by adapting `slurm_submit.sh`, which submits the pipeline to the cluster, `slurm_status.sh`, which handles detection of job status on the cluster, and `cluster.yaml`, which specifies the resource requests for each type of job.

