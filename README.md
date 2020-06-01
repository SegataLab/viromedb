# viromeDB

The code organized in this repo recapitulates the computational steps for the assembly and curation of the sequences of ViromeDB.

The code provided here is not optimized for universal use and is released for information and reproducibility puroposes only. This code cannot be run "as is" and needs to be adapted to your storage and computational architecture to be used.

## ▶ Step 0: Virome Assembly ##

Viromes with a [ViromeQC](https://github.com/SegataLab/viromeqc) enrichment score > 50 are assembled into contigs and screened with this pipeline (`asemble_sample.sh`)

![](https://github.com/SegataLab/viromedb/blob/master/doc/img/vlp_viromes_1.jpg)

- Reads are cleaned to remove low quality reads with Trim-Galore
- Human hg19 removal (Bowtie2)
- split-and-sort.py is used to recover read-pairs and put them in the correct order
- Metagenomic assembly performed with (contigs longer than 500 bp are kept)
	- Spades (if the original reads are paired-end)
	- Megahit (if the original reads are single-end)
- Contigs are fed into Prokka with `--kingdom Viruses`

We use `launch_assembly_on_dataset.sh` to submit jobs to the HPC cluster.

## ▶ Step 1: Contigs Mapping ##

In this step, the assembled contigs are mapped against a list of interesting targets (e.g. known viruses, viral proteins...) and unwanted targets (e.g. bacterial genomes, MAGs)

### BLASTn against target databases

The `blast_contig.sh` script organize the mappings of each reconstructed contig against:

- NCBI80k (List of 80,000 reference genomes from NCBI)
- SGBs (The database of 154k genomes from __Pasolli et al.__, 2019 plus all the unbinned contigs not released in the original study)
- Viral Genomes (from RefSeq release 91)

### HMMscan against target protein models

Each reconstructed contig is mapped with `hmm_contig.sh` against the following models:

- [pFAM.A](https://academic.oup.com/nar/article/47/D1/D427/5144153)
- [vFAM 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4139300/)
- vFAM 2019 (custom built basing on RefSeq Viral Proteins, release 91)

### BLAST and HMM profiles merging

A set of three scripts (`do_bread.sh`, `bread.py`, `bread_unify.py`) analyzes the BLAST outputs against multiple databases and, for each contig, report the best match for each db. 

For each contig and each database (i.e. NCBI80k, SGBs and RefSeq_Viral_Genomes), the script reports:

- The best hit `label`, `perc. of identity`, `bitscore` and `alignment length` 
- The total `breadth of coverage`against that database
- The best hit is the one with the highest bitscore
- BLAST hits are filtered for identity (> 80%) and length (> 1000)
- The overall breadth of coverage is calculated by merging together all the hits that align on the contig. The script hence reports the overall breadth.


## ▶ Step 2: Contigs Filtering ##

This step of the pipeline uses `analyze_contigs.py` and `analyze_contigs_merge.py` to process the contigs assembled from viromes and to structure the data available in different files for each contig.

`analyze_contigs.py` can be parallelized on multiple FASTA files to speed up the process. `analyze_contigs_merge.py` unifies the data from different runs of `analyze_contigs.py`.

The resulting output is a (filtered) selection of all the assembled contigs, merged with:

- The metadata for each sample
- Information on the best hit in each of the bacterial/viral databases (RefSeq, SGBs etc)
- a set of statistics on the prevalence of each contig across metagenomes and virome datasets (including the >9,000 metagenoms described in __Pasolli et al., 2019__, and the screened viromes described in __Zolfo et al__, 2019)

Finally, `extract_contigs_from_vdb_report.py` reads the filtered contigs list and extracts the associated sequences in `FASTA` format.

## ▶ Step 3: Contigs Clustering ##

![](https://github.com/SegataLab/viromedb/blob/master/doc/img/vlp_viromes_3.jpg)

This steps clusters the filtered contis, then runs a BLAST search on assembled metagenomes and viromes to retrieve homologues in there. The clustering is then performed again to produce multiple sequence alignments.

1. Contigs from highly-enriched viromes are clustered with vsearch at 90% identity (i.e. **high enrinchment contigs**)
2. Contigs are then mapped against metagenomes (n~=9000) and viromes (n~=3100) to find homologous sequences, producing an extended contigs repertoire (i.e. **extended contigs**)
3. The extended contigs are then compared to the centroids of the initial clustering with mash. Contigs with a distance < 10% are kept.
4. A second clustering is performed internally within each cluster, by adding the extended contigs to each cluster. Contigs that still fall in the same cluster of the **high enrinchment contigs** are kept. (i.e. **the final clusters**)
5. Within each final cluster, only sequences with a length +/- 20% of the median length of thee cluster go through MSA and a phylogenetic tree is produced.

*Documentation: Work in Progress*

## Utility data and scripts ##

### HPC Launchers 

Internal Scripts to launch jobs on the PBS HPC cluster of the University of Trento are named `*_launcher`. These scripts gradually submit jobs to the HPC cluster and tack the execution of each job.

examples are:

```
vdb_contig_filtering/hpc_launchers/vdb_build_launcher.sh
vdb_contig_filtering/hpc_launchers/vdb_build_launch.sh
vdb_assembly/launch_assembly_on_dataset.sh
```

### sequenceEdit.py

Edits FASTA and FASTQ files to modify Sequences IDs (e.g. to propagate metadata and information through the various step of the pipeline)

### sequenceSplit.py

Creates a separate FASTA file for each entry of a bigger FASTA file 

### split_and_sort2.py 

Splits reads into `R1` and `R2` files, then checks intact read-pairs and puts singletons into `UP`. This ensures that all the quality-filtered reads are used to produce an assembly, and that unpaired reads are not discarded.
