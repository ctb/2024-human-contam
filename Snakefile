KSIZE=51

DB = f"gtdb-rs214-k{KSIZE}.zip"
HG38 = "hg38-entire.sig.zip"

###

# load/discover FASTQ files
FILES, = glob_wildcards("ont-duplex/{name}_pass.fastq.gz")
print(f"found {len(FILES)} files.")

rule all:
    input:
        expand("sketches/{name}_pass.sig.zip", name=FILES),
        expand("output/{name}_pass.gather.csv", name=FILES),
        "output/all-merged.gather.csv",
        "sketches/all-merged.with-hg38.sig.zip"

rule sketch:
    input:
        expand("sketches/{name}_pass.sig.zip", name=FILES),
        expand("output/{name}_pass.gather.csv", name=FILES),

rule subtract:
    input:
        expand("sketches/{name}_pass.sub-hg38.sig.zip", name=FILES),
        "sketches/all-merged.sub-hg38.sig.zip"

rule gather:
    input:
        expand("sketches/{name}_pass.sig.zip", name=FILES),
        expand("output/{name}_pass.gather.csv", name=FILES),
        "output/all-merged.gather.csv",

rule sketch_sample:
    input:
        "ont-duplex/{name}_pass.fastq.gz",
    output:
        "sketches/{name}_pass.sig.zip",
    shell: """
        sourmash sketch dna -p k=21,k=31,k=51,abund {input} -o {output} \
            --name {wildcards.name:q}
    """

### worker rules

rule hg38_extract_k21:
    input:
        HG38
    output:
        f"hg38-entire.k{KSIZE}.sig.gz"
    shell:
        "sourmash sig cat {input} -o {output} -k {KSIZE}"

# subtract hg38, doing the necessary sig.gz/flatten/subtract/inflate/rename
# stuff. ugh.
rule subtract_hg38:
    input:
        orig="sketches/{name}_pass.sig.zip",
        hg38=f"hg38-entire.k{KSIZE}.sig.gz",
    output:
        sub="sketches/{name}_pass.sub-hg38.sig.zip",
        sig_gz=temporary("sketches/{name}_pass.sig.gz"),
        sub_gz=temporary("sketches/{name}_pass.sub-h38.sig.gz"),
    shell: """
        # convert .sig.zip to .sig.gz b/c subtract wants it, sigh.
        sourmash sig cat {input.orig} -o {output.sig_gz}

        # subtract hg38 from sketch, retaining abundances.
        sourmash sig subtract {output.sig_gz} {input.hg38} -k {KSIZE} --flatten \
           --abundances-from {output.sig_gz} -o {output.sub_gz} 

        # rename to something pleasing.
        sourmash sig rename {output.sub_gz} {wildcards.name:q} -o {output.sub}
    """

rule gather_sample:
    input:
        sub="sketches/{name}_pass.sub-hg38.sig.zip",
        db=DB,
    output:
        csv="output/{name}_pass.gather.csv",
        out="output/{name}_pass.gather.out",
        err="output/{name}_pass.gather.err",
    threads: 64
    shell: """
        /usr/bin/time -v sourmash scripts fastgather {input.sub} \
            -k {KSIZE} -c {threads} -t 10000 -s 10000 \
            {input.db} -o {output.csv} > {output.out} 2> {output.err}
    """

rule merge_samples:
    input:
        expand("sketches/{name}_pass.sig.zip", name=FILES)
    output:
        "sketches/all-merged.with-hg38.sig.zip"
    shell: """
        sourmash sig merge -k {KSIZE} {input} -o {output} --name "all-WGS-merged"
    """

rule merge_samples_sub:
    input:
        expand("sketches/{name}_pass.sub-hg38.sig.zip", name=FILES),
    output:
        "sketches/all-merged.sub-hg38.sig.zip"
    shell: """
        sourmash sig merge -k {KSIZE} {input} -o {output} --name "all-WGS-merged"
    """

rule gather_merged_samples_sub:
    input:
        query="sketches/all-merged.sub-hg38.sig.zip",
        db=DB,
    output:
        csv="output/all-merged.gather.csv",
        out="output/all-merged.gather.out",
        err="output/all-merged.gather.err",
    threads: 128
    shell: """
        /usr/bin/time -v sourmash scripts fastgather {input.query} \
            -k {KSIZE} -c {threads} -t 10000 -s 10000 \
            {input.db} -o {output.csv} > {output.out} 2> {output.err}
    """

rule prepare_tax_sqldb:
    input:
        csv='gtdb-rs214.lineages.csv.gz'
    output:
        sqldb='gtdb-rs214.lineages.sqldb'
    shell: """
        sourmash tax prepare -t {input} -o {output}
    """

rule summarize_tax_on_merged:
    input:
        csv="output/all-merged.gather.csv",
        sqldb='gtdb-rs214.lineages.sqldb',
    shell: """
        sourmash tax metagenome -t {input.sqldb} -g {input.csv} -F human
    """
