# Snakefile
#
# Build GraftM gene packages from seed sequences and UniRef90 searches.
#
# Expected config.yaml structure:
#
# uniref_db: "/path/to/uniref90"
# uniref_fasta: "/path/to/uniref90.fasta"
#
# threads:
#   mmseqs: 10
#
# genes:
#   genE:
#     seeds: "inputs/genE_seeds.faa"
#     taxonomy: "inputs/genE_seeds_tax.tsv"
#
#   adhA:
#     seeds: "inputs/adhA_seeds.faa"
#     taxonomy: "inputs/adhA_seeds_tax.tsv"

from pathlib import Path

########################
# Configuration
########################
configfile: "config.yaml"

INPUT_DIR = Path(config.get("input_dir", "inputs"))

UNIREF_DB = config["uniref_db"]
UNIREF_FASTA = config["uniref_fasta"]

# Resources for each tool.  These will be overridden by 
# anything in config.yaml.
MMSEQS_DB_CREATE_THREADS = (
    config.get("threads", {})
    .get("mmseqs_db_create", 20)
)
MMSEQS_THREADS = config.get("threads", {}).get("mmseqs", 10)
GRAFTM_THREADS = config.get("threads", {}).get("graftm", 4)

MMSEQS_DB_CREATE_MEM = (config.get("resources", {}).get("mmseqs_db_create", {}).get("mem_mb", 16000))
MMSEQS_MEM = config.get("resources", {}).get("mmseqs", {}).get("mem_mb", 16000)
GRAFTM_MEM = config.get("resources", {}).get("graftm", {}).get("mem_mb", 16000)


########################
# Find genes in input directory
########################
configured_genes = config.get("genes", [])

if configured_genes:
    # If individual genes are listed, that list takes precedence.
    GENES = configured_genes

else:
    # Otherwise find genes in the input directory.
    seed_files = sorted(INPUT_DIR.glob("*_seeds.faa"))

    GENES = [
        seed_file.name.removesuffix("_seeds.faa")
        for seed_file in seed_files
    ]

    if not GENES:
        raise ValueError(
            f"No genes were specified in config.yaml and no seed files "
            f"were found in {INPUT_DIR}.\n"
            f"Expected files matching:\n"
            f"  {INPUT_DIR}/<gene>_seeds.faa"
        )

########################
# Helper functions
########################

def get_seeds(wildcards):
    return str(INPUT_DIR / f"{wildcards.gene}_seeds.faa")


def get_seed_taxonomy(wildcards):
    return str(INPUT_DIR / f"{wildcards.gene}_seeds_tax.tsv")

def final_package(gene):
    """
    Use the rerooted package as the final target if a manually
    rooted tree exists; otherwise build the draft package.
    """
    rooted_tree = Path(f"results/{gene}/rooted.tree")

    if rooted_tree.exists():
        return f"results/{gene}/{gene}_rooted.gpkg"

    return f"results/{gene}/{gene}_draft.gpkg"

########################
# Main targets
########################

rule all:
    input:
        [final_package(gene) for gene in GENES]

##############################
# Build UniRef MMseqs database
##############################

rule build_uniref_db:
    input:
        fasta=UNIREF_FASTA
    output:
        db=UNIREF_DB,
        dbtype=UNIREF_DB + ".dbtype",
        index_complete=touch(UNIREF_DB + ".index_complete")
    params:
        tmp=UNIREF_DB + "_tmp"
    threads:
        MMSEQS_DB_CREATE_THREADS
    resources:
        mem_mb=MMSEQS_DB_CREATE_MEM,
    conda:
        "envs/graftm.yaml"
    log:
        "logs/mmseqs_build_uniref.log"
    shell:
        r"""
        mkdir -p $(dirname {output.db})
        mkdir -p {params.tmp}
        mkdir -p $(dirname {log})

        mmseqs createdb \
            {input.fasta} \
            {output.db} \
            --threads {threads} \
            >> {log} 2>&1

        mmseqs createindex \
            {output.db} \
            {params.tmp} \
            --threads {threads} \
            >> {log} 2>&1

        touch {output.index_complete}
        """

########################
# Search UniRef90
########################

rule search_uniref:
    input:
        seeds=get_seeds
        db=rules.build_uniref_db.output.db,
        dbtype=rules.build_uniref_db.output.dbtype,
        index_complete=rules.build_uniref_db.output.index_complete
    output:
        search="results/{gene}/uniref90_search.m8"
    params:
        db=UNIREF_FASTA,
        tmp="results/{gene}/tmp"
    threads:
        MMSEQS_THREADS
    resources: 
        mem_mb=MMSEQS_MEM
    conda:
        "envs/graftm.yaml"
    log:
        "results/{gene}/mmseqs.log"
    shell:
        r"""
        mkdir -p {params.tmp}

        mmseqs easy-search \
            --threads {threads} \
            {input.seeds} \
            {params.db} \
            {output.search} \
            {params.tmp} \
            &> {log}
        """


########################
# Extract UniRef IDs
########################

rule extract_uniref_ids:
    input:
        search="results/{gene}/uniref90_search.m8"
    output:
        ids="results/{gene}/uniref90_search_ids.txt"
    shell:
        r"""
        awk -F '\t' '{{print $2}}' {input.search} \
            | sort -u \
            > {output.ids}
        """


########################
# Extract UniRef sequences
########################

rule extract_uniref_sequences:
    input:
        ids="results/{gene}/uniref90_search_ids.txt"
    output:
        fasta="results/{gene}/uniref90_search.faa"
    params:
        uniref_fasta=UNIREF_FASTA
    conda:
        "envs/graftm.yaml"
    log:
        "results/{gene}/mfqe.log"
    shell:
        r"""
        mfqe \
            --sequence-name-lists {input.ids} \
            --input-fasta {params.uniref_fasta} \
            --output-fasta-files {output.fasta} \
            --output-uncompressed \
            &> {log}
        """


########################
# Combine seed + UniRef sequences
########################

rule combine_sequences:
    input:
        seeds=get_seeds,
        uniref="results/{gene}/uniref90_search.faa"
    output:
        fasta="results/{gene}/combined.faa"
    shell:
        r"""
        cat \
            {input.seeds} \
            {input.uniref} \
            > {output.fasta}
        """


########################
# Create taxonomy for UniRef hits
########################

rule create_uniref_taxonomy:
    input:
        ids="results/{gene}/uniref90_search_ids.txt"
    output:
        taxonomy="results/{gene}/uniref90_search_tax.txt"
    shell:
        r"""
        awk \
            -v gene="{wildcards.gene}" \
            'BEGIN {{OFS="\t"}} {{print $0, gene}}' \
            {input.ids} \
            > {output.taxonomy}
        """


########################
# Combine seed + UniRef taxonomy
########################

rule combine_taxonomy:
    input:
        seeds_tax=get_seed_taxonomy,
        uniref_tax="results/{gene}/uniref90_search_tax.txt"
    output:
        taxonomy="results/{gene}/combined_tax.txt"
    shell:
        r"""
        cat \
            {input.seeds_tax} \
            {input.uniref_tax} \
            > {output.taxonomy}
        """


########################
# Build draft GraftM package
########################

rule graftm_draft:
    input:
        sequences="results/{gene}/combined.faa",
        taxonomy="results/{gene}/combined_tax.txt"
    output:
        package=directory(
            "results/{gene}/{gene}_draft.gpkg"
        )
    log:
        "results/{gene}/GraftM_draft.log"
    conda:
        "envs/graftm.yaml"
    threads:
        GRAFTM_THREADS
    resources: 
        mem_mb=GRAFTM_MEM
    shell:
        r"""
        set +e

        cd results/{wildcards.gene}

        graftM create \
            --sequences combined.faa \
            --taxonomy combined_tax.txt \
            --output {wildcards.gene}_draft.gpkg \
            &> GraftM_draft.log

        status=$?

        if [ "$status" -ne 0 ]; then

            echo "" >&2
            echo "==================================================" >&2
            echo "GraftM package creation failed for {wildcards.gene}" >&2
            echo "==================================================" >&2
            echo "" >&2

            if grep -Eqi \
                'root|reroot|outgroup' \
                GraftM_draft.log
            then
                echo "Possible automatic-rooting failure." >&2
                echo "" >&2
                echo "Inspect:" >&2
                echo "    results/{wildcards.gene}/GraftM_draft.log" >&2
                echo "" >&2
                echo "If manual rooting is required, create:" >&2
                echo "" >&2
                echo "    results/{wildcards.gene}/rooted.tree" >&2
                echo "" >&2
                echo "and then rerun Snakemake." >&2
            else
                echo "The failure does not obviously appear to be" >&2
                echo "related to tree rooting." >&2
                echo "" >&2
                echo "Inspect:" >&2
                echo "    results/{wildcards.gene}/GraftM_draft.log" >&2
            fi

            exit "$status"
        fi
        """

########################
# Build GraftM package with rerooted tree
########################
rule graftm_rerooted:
    input:
        sequences="results/{gene}/combined.faa",
        rooted_tree="results/{gene}/rooted.tree",
        taxonomy_csv="results/{gene}/graftm_create_taxonomy.combined.csv",
        seqinfo_csv="results/{gene}/graftm_create_seqinfo.combined.csv",
        alignment="results/{gene}/graftm_create_alignment.combined.faa"
    output:
        package=directory(
            "results/{gene}/{gene}_rooted.gpkg"
        )
    threads:
        GRAFTM_THREADS
    resources:
        mem_mb=GRAFTM_MEM
    conda:
        "envs/graftm.yaml"
    log:
        "results/{gene}/GraftM_reroot.log"
    shell:
        r"""
        graftM create \
            --taxtastic_taxonomy {input.taxonomy_csv} \
            --taxtastic_seqinfo {input.seqinfo_csv} \
            --alignment {input.alignment} \
            --rerooted_tree {input.rooted_tree} \
            --sequences {input.sequences} \
            --output {output.package} \
            &> {log}
        """