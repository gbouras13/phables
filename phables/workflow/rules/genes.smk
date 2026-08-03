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


"""
Structural phage-gene detection: ProstT5 -> 3Di -> foldseek search against a
hallmark PHROG structure subDB, selected via --phage-detection prostt5-foldseek.

Gated behind `if PD == "prostt5-foldseek"`, unlike scan_phrogs above -- not just an
optimisation. hallmark_db/prostt5_checkpoint/etc default to empty in config.yaml
when unused, and Snakemake parses every rule's input/output at load time regardless
of whether it ends up in the target list, so an un-gated `hallmark_db = config["hallmark_db"]`
input fails with `SyntaxError: Input and output files have to be specified as
strings or lists of strings` on a plain `phables run` even when this path is never
selected. Only relying on the target list (as scan_phrogs does) isn't enough here
because these params, not just outputs, can be structurally invalid at parse time.

These consume PROTEINS_FILE, the same protein set the marker-gene search uses
(genes.smk top), unifying what used to be two independent detection paths (see
the fork audit, section 2.1: scan_phrogs previously ran mmseqs directly on the
nucleotide edges FASTA and never shared gene calls with the marker-gene path).
"""

if PD == "prostt5-foldseek":

    QUERY_3DI = os.path.join(OUTDIR, "preprocess", "hallmark", "proteins_3di.fasta")


    rule predict_3di:
        input:
            faa = PROTEINS_FILE,
        threads:
            config["resources"]["jobCPU"]
        resources:
            mem_mb = config["resources"]["jobMem"],
            mem = str(config["resources"]["jobMem"]) + "MB"
        output:
            threedi = QUERY_3DI
        params:
            checkpoint = config["prostt5_checkpoint"],
            model_name = config["prostt5_model"],
            model_dir = config["prostt5_model_dir"],
            half_precision = config["prostt5_half_precision"],
            cpu = config["prostt5_cpu"],
            max_residues = config["prostt5_max_residues"],
            max_seq_len = config["prostt5_max_seq_len"],
            max_batch = config["prostt5_max_batch"],
        log:
            os.path.join(LOGSDIR, "predict_3di.log")
        conda:
            os.path.join("..", "envs", "prostt5.yaml")
        script:
            os.path.join("..", "scripts", "predict_3di.py")


    rule build_hallmark_query_db:
        input:
            faa = PROTEINS_FILE,
            threedi = QUERY_3DI,
        output:
            db = os.path.join(OUTDIR, "preprocess", "hallmark", "query_db"),
        params:
            out_prefix = os.path.join(OUTDIR, "preprocess", "hallmark", "query_db"),
        log:
            os.path.join(LOGSDIR, "build_hallmark_query_db.log")
        conda:
            os.path.join("..", "envs", "foldseek.yaml")
        script:
            os.path.join("..", "scripts", "build_foldseek_query_db.py")


    rule scan_hallmark:
        input:
            db = os.path.join(OUTDIR, "preprocess", "hallmark", "query_db"),
            hallmark_db = config["hallmark_db"],
        threads:
            config["resources"]["jobCPU"]
        resources:
            mem_mb = config["resources"]["jobMem"],
            mem = str(config["resources"]["jobMem"]) + "MB"
        output:
            os.path.join(OUTDIR, "preprocess", "hallmark_hits.tsv")
        params:
            query_prefix = os.path.join(OUTDIR, "preprocess", "hallmark", "query_db"),
            result = os.path.join(OUTDIR, "preprocess", "hallmark", "result"),
            tmp = os.path.join(OUTDIR, "preprocess", "hallmark", "tmp"),
        log:
            os.path.join(LOGSDIR, "scan_hallmark.log")
        conda:
            os.path.join("..", "envs", "foldseek.yaml")
        shell:
            """
            foldseek search {params.query_prefix} {input.hallmark_db} {params.result} {params.tmp} \
                --threads {threads} -s 7 -e {config[hallmark_evalue]} > {log}
            foldseek convertalis {params.query_prefix} {input.hallmark_db} {params.result} {output} \
                --format-output query,target,evalue,bits,fident --threads {threads} >> {log}
            rm -rf {params.tmp}
            """