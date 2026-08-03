"""
Call genes on unitigs, then use HMMER to scan for bacterial single-copy marker genes.
Use mmseqs2 to scan for PHROGs in unitigs.

Gene calling is a separate rule from the marker-gene search so the caller can be
swapped via --genecaller without touching anything downstream. Both callers emit
the same `{seqid}_{start}_{end}_{strand}` protein id convention, which
gene_utils.get_smg_unitigs relies on to recover the unitig name.
"""

PROTEINS_FILE = EDGES_FILE + ".frag.faa"


if GC == "pyrodigal-gv":

    rule call_genes:
        input:
            genome = EDGES_FILE,
        threads:
            config["resources"]["jobCPU"]
        resources:
            mem_mb = config["resources"]["jobMem"],
            mem = str(config["resources"]["jobMem"]) + "MB"
        output:
            faa = PROTEINS_FILE
        log:
            os.path.join(LOGSDIR, "gene_call_pyrodigal_gv.log")
        conda:
            os.path.join("..", "envs", "genecall.yaml")
        script:
            os.path.join("..", "scripts", "gene_caller.py")

else:

    rule call_genes:
        input:
            genome = EDGES_FILE,
        threads:
            config["resources"]["jobCPU"]
        resources:
            mem_mb = config["resources"]["jobMem"],
            mem = str(config["resources"]["jobMem"]) + "MB"
        output:
            faa = PROTEINS_FILE
        params:
            frag = EDGES_FILE + ".frag",
        log:
            out = os.path.join(LOGSDIR, "gene_call_fraggenescan_out.log"),
            err = os.path.join(LOGSDIR, "gene_call_fraggenescan_err.log"),
        conda:
            os.path.join("..", "envs", "smg.yaml")
        shell:
            """
                run_FragGeneScan.pl -genome={input.genome} -out={params.frag} -complete=0 -train=complete -thread={threads} 1>{log.out} 2>{log.err}
            """


rule scan_smg:
    input:
        faa = PROTEINS_FILE,
        hmm = os.path.join(DBPATH, "marker.hmm"),
    threads:
        config["resources"]["jobCPU"]
    resources:
        mem_mb = config["resources"]["jobMem"],
        mem = str(config["resources"]["jobMem"]) + "MB"
    output:
        hmmout = os.path.join(OUTDIR, "preprocess", "edges.fasta.hmmout")
    log:
        hmm_out=os.path.join(LOGSDIR, "smg_scan_hmm_out.log"),
        hmm_err=os.path.join(LOGSDIR, "smg_scan_hmm_err.log")
    conda:
        os.path.join("..", "envs", "smg.yaml")
    shell:
        """
            hmmsearch --domtblout {output.hmmout} --cut_tc --cpu {threads} {input.hmm} {input.faa} 1>{log.hmm_out} 2> {log.hmm_err}
        """


rule scan_phrogs:
    input:
        genome = EDGES_FILE,
        db = os.path.join(DBPATH,"phrogs_mmseqs_db","phrogs_profile_db")
    threads:
        config["resources"]["jobCPU"]
    resources:
        mem_mb = config["resources"]["jobMem"],
        mem = str(config["resources"]["jobMem"]) + "MB"
    output:
        os.path.join(OUTDIR, "preprocess", "phrogs_annotations.tsv")
    params:
        out_path = os.path.join(OUTDIR, "preprocess", "phrogs"),
        target_seq = os.path.join(OUTDIR, "preprocess", "phrogs", "target_seq"),
        results_mmseqs = os.path.join(OUTDIR, "preprocess", "phrogs", "results_mmseqs"),
        tmp = os.path.join(OUTDIR, "preprocess", "phrogs", "tmp"),
    log:
        os.path.join(LOGSDIR, "phrogs_scan.log")
    conda: 
        os.path.join("..", "envs", "mmseqs.yaml")
    shell:
        """
        mkdir -p {params.out_path}
        mmseqs createdb {input} {params.target_seq} > {log}
        mmseqs search {params.target_seq} {input.db} {params.results_mmseqs} {params.tmp} --threads {threads} -s 7 > {log}
        mmseqs createtsv {params.target_seq} {input.db} {params.results_mmseqs} {output} --threads {threads} --full-header > {log}
        rm -rf {params.out_path}
        """