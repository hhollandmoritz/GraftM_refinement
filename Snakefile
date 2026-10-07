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


########################
# Configuration
########################

GENES = list(config["genes"].keys())

UNIREF_DB = config["uniref_db"]
UNIREF_FASTA = config["uniref_fasta"]

MMSEQS_THREADS = config.get("threads", {}).get("mmseqs", 10)


########################
# Helper functions
########################

def get_seeds(wildcards):
    return config["genes"][wildcards.gene]["seeds"]


def get_seed_taxonomy(wildcards):
    return config["genes"][wildcards.gene]["taxonomy"]


########################
# Main targets
########################

rule all:
    input:
        expand(
            "results/{gene}/{gene}_draft.gpkg",
            gene=GENES
        )


########################
# Search UniRef90
########################

rule search_uniref:
    input:
        seeds=get_seeds
    output:
        search="results/{gene}/uniref90_search.m8"
    params:
        db=UNIREF_DB,
        tmp="results/{gene}/tmp"
    threads:
        MMSEQS_THREADS
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
    shell:
        r"""
        set +e

        graftM create \
            --sequences {input.sequences} \
            --taxonomy {input.taxonomy} \
            --output {output.package} \
            &> {log}

        status=$?

        if [ "$status" -ne 0 ]; then

            echo "" >&2
            echo "==================================================" >&2
            echo "GraftM package creation failed for {wildcards.gene}" >&2
            echo "==================================================" >&2
            echo "" >&2

            if grep -Eqi \
                'root|reroot|outgroup' \
                {log}
            then
                echo "Possible automatic-rooting failure." >&2
                echo "" >&2
                echo "Inspect:" >&2
                echo "    {log}" >&2
                echo "" >&2
                echo "If manual rooting is required, create:" >&2
                echo "" >&2
                echo "    results/{wildcards.gene}/manual_rooted.tree" >&2
                echo "" >&2
                echo "and then run the manual-rooting rule." >&2
            else
                echo "The failure does not obviously appear to be" >&2
                echo "related to tree rooting." >&2
                echo "" >&2
                echo "Inspect:" >&2
                echo "    {log}" >&2
            fi

            exit "$status"
        fi
        """