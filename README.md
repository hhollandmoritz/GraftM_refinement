# GraftM Gene Package Builder

This Snakemake workflow builds gene-specific [GraftM](https://github.com/geronimp/graftM) packages from seed protein sequences and seed taxonomy assignments. It is inteneded to extend and refine the annotation capabilities for genes that are typically poorly identified using typical HMMs and annotation tools. 

For each gene of interest, the workflow:

1. Searches seed sequences against UniRef90 using MMseqs2.
2. Extracts matching UniRef90 sequences.
3. Combines UniRef hits with the original seed sequences.
4. Creates a combined taxonomy file.
5. Builds a draft GraftM package.
6. Notifies user that manual intervention is needed when automatic tree rooting is not possible.
7. Can incorporate manually-rooted trees to create a final package.

## Installation and Setup

### Database
This workflow requires the Uniref90 database. If you do not already have it available, you will need to download it.

### Software
The pipeline comes packaged with instructions for creating the necessary environment `envs/graftm.yaml`. Before running the pipeline, run the installation on a node that has access to the internet: 

```bash
snakemake \
    --profile profiles/local \
    --conda-create-envs-only
```
If you want the conda environment installed somewhere other than `.snakemake/conda`, make sure to change the conda-prefix setting in each of the `profiles/*/config.yaml` files. 
## Usage

The workflow can be run either locally or on a SLURM cluster using the Snakemake profiles.

### 1\. Configure the workflow

Before running, edit `config.yaml` to specify:

* the UniRef90 MMseqs2 database path,
* the UniRef90 FASTA path,
* thread defaults,
* and the genes/packages to build.

For example:

```yaml
uniref_db: "/path/to/uniref90/uniref90"
uniref_fasta: "/path/to/uniref90/uniref90.fasta"

threads:
  mmseqs: 10
  graftm: 4

genes:
  genE:
    seeds: "inputs/genE_seeds.faa"
    taxonomy: "inputs/genE_seeds_tax.tsv"
```

Either the folder "inputs" will be searched for seed and taxonomy files, or each gene listed under `genes:` will be processed separately. When genes are declared in the config file, they take precedence over the search of the input directory. 

The file names are sensitive to misspellings. 

All seed fasta files should be named: `[gene]_seeds.faa`
All seed taxonomy files should be named: `[gene]_seeds_tax.tsv`

### 2\. Run workflow

To run on a local workstation:

```bash
snakemake \\
    --profile profiles/local
```

The maximum number of cores available to the workflow is controlled in:

```text
profiles/local/config.yaml
```

To run on a SLURM cluster:

```bash
snakemake \\
    --profile profiles/slurm
```

Cluster settings such as job limits, accounts, partitions, and optional resource overrides can be set in:

```text
profiles/slurm/config.yaml
```

### 3\. Manual rooting

If GraftM cannot automatically root a reference tree, the draft package step will fail and the corresponding log should be inspected:

```text
results/<gene>/GraftM\_draft.log
```

Inspect the generated alignment and tree, then create a manually rooted tree and save it as:

```text
results/<gene>/rooted.tree
```

In development: Rerunning Snakemake will allow the pipeline to continue from that point without repeating the UniRef search.

## Directory structure

```text
graftm_refinement/
├── Snakefile
├── config.yaml
├── README.md
├── inputs/
│   ├── genE_seeds.faa
│   ├── genE_seeds_tax.tsv
│   ├── adhA_seeds.faa
│   └── adhA_seeds_tax.tsv
├── envs/
│   └── graftm.yaml
├── profiles/
│   ├── local/
│   │   └── config.yaml
│   └── slurm/
│       └── config.yaml
└── results/
```

## References

This Snakemake implementation is based on the original gene-package construction script written by Sam Aroney and subsequently modified and reorganized by Hannah Holland-Moritz.

If you find this workflow useful, please consider citing the references below:

The original version of this GraftM package-building workflow was developed for analyses used in:

Cronin, D. R., Holland-Moritz, H., Smith, D. A., Aroney, S. T. N., Hodgkins, S., Borton, M., Li, Y.-F., Healy, K., IsoGenie Field & Analytic Teams 2010–2017, EMERGE Institute Coordinators, Tfaily, M. M., Crill, P., McCalley, C. K., Wrighton, K., Varner, R. K., Tyson, G. W., Woodcroft, B., Bagby, S. C., Ernakovich, J., & Rich, V. I.  **Stable states in an unstable landscape: microbial resistance at the front line of climate change.**  
bioRxiv (2025). https://doi.org/10.1101/2025.02.07.636677

### GraftM

GraftM is used to construct and apply phylogenetically informed gene packages.

Boyd, J. A., Woodcroft, B. J., & Tyson, G. W. (2018). **GraftM: a tool for scalable, phylogenetically informed classification of genes within metagenomes.** *Nucleic Acids Research*, 46(10), e59. https://doi.org/10.1093/nar/gky174

### MMseqs2

MMseqs2 is used to search seed sequences against the UniRef90 database.

Mirdita, M., Steinegger, M., & Söding, J. (2019). **MMseqs2 desktop and local web server app for fast, interactive sequence searches.**  
*Bioinformatics*, 35(16), 2856–2858. https://doi.org/10.1093/bioinformatics/bty1057

### mfqe

`mfqe` is used to extract selected protein sequences from the UniRef90 FASTA file.

The software is distributed through Bioconda and maintained at: https://github.com/wwood/mfqe

### UniRef90

UniRef90 provides the protein sequence database searched by MMseqs2.

Suzek, B. E., Wang, Y., Huang, H., McGarvey, P. B., Wu, C. H., & UniProt Consortium. (2015). **UniRef clusters: a comprehensive and scalable alternative for improving sequence similarity searches.** *Bioinformatics*, 31(6), 926–932.  
https://doi.org/10.1093/bioinformatics/btu739

Users should also record the UniRef90 release or download date used for each analysis.
