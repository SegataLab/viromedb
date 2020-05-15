# viromeDB
Scripts for the ViromeDB project


## Virome Assembly Pipeline ##

Viromes with a ViromeQC enrichment score > 50 were assembled into contigs with this pipeline.

### assemble_sample.sh

- Reads cleaned to remove low quality reads with Trim-Galore
- Human hg19 removal (Bowtie2)
- split-and-sort.py is used to recover read-pairs and put them in the correct order
- Metagenomic assembly performed with (contigs longer than 500 bp are kept)
	- Spades (if the original reads are paired-end)
	- Megahit (if the original reads are single-end)
- Contigs are fed into Prokka with `--kingdom Viruses`

`launch_assembly_on_dataset.sh` Is the launcher script for the HPC cluster.

## Contigs Mapping ##

### blast_contig.sh / hmm_contig.sh

These scripts organize the mappings of each reconstructed contig against:

- NCBI80k (List of 80,000 reference genomes from NCBI)
- SGBs (The database of 154k genomes from __Pasolli et al.__, 2019)
- Viral Genomes (from RefSeq release 91)

And performs HMMscan mapping against:

- pFAM
- vFAM 2014
- vFAM 2019 (based on RefSeq Viral Proteins, release 91)

### do_bread.sh / bread.py / bread_unify.py

These scripts analyze the BLAST outputs against multiple databases and, for each contig, report the best match for each db. 

For each contig and each database (i.e. NCBI80k, SGBs and RefSeq_Viral_Genomes), the script reports:
- The best hit `label`, `perc. of identity`, `bitscore` and `alignment length` 
- The total breadth of coverage against that database

- The best hit is the one with the highest bitscore
- BLAST hits are filtered for identity (> 80%) and length (> 1000)
- The overall breadth of coverage is calculated by merging together all the hits that align on the contig. The script hence reports the overall breadth.


### prokka_contig.sh

Runs Prokka with `--kingdom Viruses` on the assembled contigs.

### 

## Contigs Filtering ##

analyze_contigs_merge_largethreshold2.py
analyze_contigs_merge_largethreshold.py

This step of the pipeline uses `analyze_contigs.py` and `analyze_contigs_merge.py` to process the contigs assembled from viromes and to structure the data available in different files for each contig. `analyze_contigs.py` can be parallelized on multiple FASTA files to speed up the process. `analyze_contigs_merge.py` unifies the data from different runs of `analyze_contigs.py`.

The resulting output is a (filtered) selection of all the assembled contigs, together with:
- The metadata for each sample
- Information on the best hit in each of the bacterial/viral databases (RefSeq, SGBs etc)
- a set of statistics on the prevalence of each contig across metagenomes and virome datasets (including the >9,000 metagenoms described in __Pasolli et al., 2019__, and the screened viromes described in __Zolfo et al__, 2019)


Then, `extract_contigs_from_vdb_report.py` reads the filtered contigs list and extracts the associated sequences in `FASTA` format.

## Contigs Clustering ##

WIP

## Utility scripts ##

### HPC Launchers 

Internal Scripts to launch jobs on the PBS HPC cluster of the University of Trento are located in the `hpc_launchers` folder:

Scripts to build the ViromeDB from Virome-assembled-contigs
```
vdb_build_launcher.sh
vdb_build_launch.sh
```

Launches Virome assembly on all samples of a folder

```
launch_assembly_on_dataset.sh
```


### sequenceEdit.py

edits FASTA and FASTQ files to modify Sequences IDs (e.g. to propagate metadata and information through the various step of the pipeline)

### sequenceSplit.py

creates a separate FASTA file for each entry of a bigger FASTA file 