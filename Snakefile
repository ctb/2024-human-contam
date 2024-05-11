DB = "gtdb-rs214-k21.zip"

###

# load/discover FASTQ files
FILES, = glob_wildcards("ont-duplex/{name}_pass.fastq.gz")
print(f"found {len(FILES)} files.")

rule all:
    input:
        expand("sketches/{name}_pass.sig.zip", name=FILES),
        expand("output/{name}_pass.gather.csv", name=FILES),

rule sketch:
    input:
        expand("sketches/{name}_pass.sig.zip", name=FILES),
        expand("output/{name}_pass.gather.csv", name=FILES),

rule gather:
    input:
        expand("sketches/{name}_pass.sig.zip", name=FILES),
        expand("output/{name}_pass.gather.csv", name=FILES),

rule sketch_sample:
    input:
        "ont-duplex/{name}_pass.fastq.gz",
    output:
        "sketches/{name}_pass.sig.zip",
    shell: """
        sourmash sketch dna -p k=21,abund {input} -o {output} \
            --name {wildcards.name:q}
    """

# subtract hg38, doing the necessary sig.gz/flatten/subtract/inflate/rename
# stuff. ugh.
rule subtract_hg38:
    input:
        orig="sketches/{name}_pass.sig.zip",
        hg38="hg38-entire.k21.sig.gz",
    output:
        sub="sketches/{name}_pass.sub-hg38.sig.zip",
        sig_gz=temporary("sketches/{name}_pass.sig.gz"),
        flat_sub_gz=temporary("sketches/{name}_pass.sub-h38.flat.sig.gz"),
    shell: """
        # convert .sig.zip to .sig.gz b/c subtract wants it, sigh.
        sourmash sig cat {input.orig} -o {output.sig_gz}

        # subtract hg38 from sketch; remove abundances, because subtract wants
        # it, sigh.
        sourmash sig subtract {output.sig_gz} {input.hg38} -k 21 --flatten \
           -o {output.flat_sub_gz}

        # re-inflate subtracted sketch w/abundance from original sketch,
        # and rename to something pleasing.
        sourmash sig inflate {output.sig_gz} {output.flat_sub_gz} -o - |
           sourmash sig rename - {wildcards.name:q} -o {output.sub}
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
            -k 21 -t {threads} \
            {input.db} -o {output.csv} > {output.out} 2> {output.err}
    """
