# GraftM Gene Package Builder

This Snakemake workflow builds gene-specific [GraftM](https://github.com/geronimp/graftM) packages from seed protein sequences and seed taxonomy assignments. It is inteneded to extend and refine the annotation capabilities for genes that are typically poorly identified using typical HMMs and annotation tools. 

For each gene of interest, the workflow:

1. Searches seed sequences against UniRef90 using MMseqs2.
2. Extracts matching UniRef90 sequences.
3. Combines UniRef hits with the original seed sequences.
4. Creates a combined taxonomy file.
5. Builds a draft GraftM package.
6. Supports manual intervention when automatic tree rooting is not possible.
7. Can later be extended to incorporate manual tree annotation and final package testing.

The workflow can run either:

- locally, using a configurable number of cores; or
- on a SLURM cluster using the Snakemake SLURM executor.

## Attribution

This workflow was developed by Hannah Holland-Moritz and is based on an original SLURM/bash workflow written by **Sam Aroney** for building GraftM gene packages.

The original workflow logic for UniRef searching, sequence extraction, taxonomy construction, and GraftM package creation was written by Sam Aroney and subsequently modified and reorganized into this Snakemake workflow.

## Directory structure

A typical project should look like:

```text
graftm-package-builder/
├── Snakefile
├── config.yaml
├── README.md
├── inputs/
│   ├── genE_seeds.faa
│   ├── genE_seeds_tax.tsv
│   ├── adhA_seeds.faa
│   └── adhA_seeds_tax.tsv
├── envs/
│   ├── mmseqs.yaml
│   ├── mfqe.yaml
│   └── graftm.yaml
├── profiles/
│   ├── local/
│   │   └── config.yaml
│   └── slurm/
│       └── config.yaml
└── results/

## References

This Snakemake implementation is based on the original gene-package construction script written by Sam Aroney and subsequently modified and reorganized by Hannah Holland-Moritz.

If you find this workflow useful, please consider citing the references below:

The original version of this GraftM package-building workflow was developed for analyses used in:

Cronin, D. R., Holland-Moritz, H., Smith, D. A., Aroney, S. T. N., Hodgkins, S., Borton, M., Li, Y.-F., Healy, K., IsoGenie Field & Analytic Teams 2010–2017, EMERGE Institute Coordinators, Tfaily, M. M., Crill, P., McCalley, C. K., Wrighton, K., Varner, R. K., Tyson, G. W., Woodcroft, B., Bagby, S. C., Ernakovich, J., & Rich, V. I.  
**Stable states in an unstable landscape: microbial resistance at the front line of climate change.**  
bioRxiv (2025).  
https://doi.org/10.1101/2025.02.07.636677

### GraftM

GraftM is used to construct and apply phylogenetically informed gene packages.

Boyd, J. A., Woodcroft, B. J., & Tyson, G. W. (2018).  
**GraftM: a tool for scalable, phylogenetically informed classification of genes within metagenomes.**  
*Nucleic Acids Research*, 46(10), e59.  
https://doi.org/10.1093/nar/gky174

### MMseqs2

MMseqs2 is used to search seed sequences against the UniRef90 database.

Mirdita, M., Steinegger, M., & Söding, J. (2019).  
**MMseqs2 desktop and local web server app for fast, interactive sequence searches.**  
*Bioinformatics*, 35(16), 2856–2858.  
https://doi.org/10.1093/bioinformatics/bty1057

### mfqe

`mfqe` is used to extract selected protein sequences from the UniRef90 FASTA file.

The software is distributed through Bioconda and maintained at:

https://github.com/wwood/mfqe

### UniRef90

UniRef90 provides the protein sequence database searched by MMseqs2.

Suzek, B. E., Wang, Y., Huang, H., McGarvey, P. B., Wu, C. H., & UniProt Consortium. (2015).  
**UniRef clusters: a comprehensive and scalable alternative for improving sequence similarity searches.**  
*Bioinformatics*, 31(6), 926–932.  
https://doi.org/10.1093/bioinformatics/btu739

Users should also record the UniRef90 release or download date used for each analysis.