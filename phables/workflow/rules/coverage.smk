"""
Use raw_coverage to map to calculate coverage of unitigs.
Use combine_cov to combine the coverage values of multiple samples into one file.

The three rules below (coverm_map, coverm_bam2counts, coverm_combine) replace the
former single `koverage run coverm` call (PLAN.md §4.8, Vijini's recommendation:
"Drop the koverage wrapper — call CoverM directly"). They reproduce Koverage's own
real "coverm" mode chain exactly -- confirmed against its actual source
(github.com/beardymcjohnface/Koverage, workflow/rules/coverm.smk and its shipped
default config.yaml), not reimplemented from guesswork:
  - same minimap2 mapping command (-ax sr --secondary=no) and samtools
    sort/view -F 4/index chain (coverm_map_pe there),
  - same `coverm contig` --methods list and order: count, rpkm, tpm, mean,
    covered_fraction, variance (Koverage's own shipped default, which phables
    never overrode -- confirmed no `--profile`/coverm-params override exists
    anywhere in this repo),
  - same per-sample -> long-format Sample/Contig/<method...> TSV reshape
    (coverm_combine there).
Byte-for-byte the same `sample_coverm_coverage.tsv` shape as before, so
run_combine_cov below and coverage_utils.py's BAM globbing (both downstream
consumers) need no changes. Only the preprocess step above is affected --
postprocess.smk's own `koverage_genomes` rule uses Koverage's separate native
(non-CoverM) coverage engine, not this wrapper, and is out of scope here.

Note: the phables-side docstring on run_combine_cov used to say "Covered_bases"
for column 7 -- that's actually Covered_fraction (Koverage's real 6th --methods
value, confirmed above). Harmless either way: only column 6 (Mean) is consumed
by the awk below, per notes/phables_audit.md §3.2.

Also note: like Koverage's own rule, this always maps with -ax sr (short-read
preset) regardless of `config["longreads"]` -- a pre-existing gap inherited
unchanged from Koverage, not introduced by this swap. Long-read runs' unitig
coverage here is thus mapped with a short-read minimap2 preset; worth a
follow-up, out of scope for this change.
"""

rule koverage_tsv:
    """Generate TSV of samples and reads -- still needed by postprocess.smk's
    own (still Koverage-based) koverage_genomes rule, which reads this same
    file; not needed by coverm_map below, which uses SAMPLE_READS directly."""
    output:
        os.path.join(OUTDIR, "preprocess", "phables.samples.tsv")
    params:
        SAMPLE_READS
    run:
        from metasnek import fastq_finder
        fastq_finder.write_samples_tsv(params[0], output[0])


rule coverm_map:
    """Map each sample's reads to the unitig edges with minimap2 -> sorted,
    unmapped-filtered, indexed BAM. Matches Koverage's own coverm_map_pe rule
    exactly (see module docstring above)."""
    input:
        ref = EDGES_FILE,
        r1 = lambda wildcards: SAMPLE_READS[wildcards.sample]["R1"],
    params:
        # Koverage's own rule passes "" for single-end samples (R2 None) --
        # minimap2 then maps r1 alone. Preserved as-is.
        r2 = lambda wildcards: SAMPLE_READS[wildcards.sample]["R2"] or "",
    output:
        bam = os.path.join(OUTDIR, "preprocess", "temp", "{sample}.bam"),
        bai = os.path.join(OUTDIR, "preprocess", "temp", "{sample}.bam.bai"),
    threads:
        config["resources"]["jobCPU"]
    resources:
        mem_mb = config["resources"]["jobMem"]
    conda:
        os.path.join("..", "envs", "coverm.yaml")
    log:
        os.path.join(LOGSDIR, "coverm_map.{sample}.log")
    shell:
        """
        {{ minimap2 -t {threads} -ax sr --secondary=no {input.ref} {input.r1} {params.r2} \
            | samtools sort -T {wildcards.sample} -@ {threads} - \
            | samtools view -F 4 > {output.bam} ; \
        samtools index {output.bam} ; }} 2> {log}
        """


rule coverm_bam2counts:
    """Per-sample coverage stats. Koverage's own rule doesn't pass a thread
    count to coverm either (only Snakemake's own scheduling uses `threads:`
    here) -- preserved as-is rather than adding an unverified -t flag."""
    input:
        os.path.join(OUTDIR, "preprocess", "temp", "{sample}.bam")
    output:
        os.path.join(OUTDIR, "preprocess", "temp", "{sample}.cov")
    conda:
        os.path.join("..", "envs", "coverm.yaml")
    log:
        os.path.join(LOGSDIR, "coverm_bam2counts.{sample}.log")
    shell:
        """
        coverm contig -b {input} \
            -m count -m rpkm -m tpm -m mean -m covered_fraction -m variance \
            > {output} 2> {log}
        """


rule coverm_combine:
    """Reshape per-sample coverm TSVs (each header column is
    "<bam filename> <method>") into one long-format Sample/Contig/<method...>
    table -- same reshape as Koverage's own coverm_combine rule."""
    input:
        expand(os.path.join(OUTDIR, "preprocess", "temp", "{sample}.cov"), sample=SAMPLE_NAMES)
    output:
        os.path.join(OUTDIR, "preprocess", "results", "sample_coverm_coverage.tsv")
    run:
        with open(input[0]) as f:
            header = f.readline().rstrip("\n").split("\t")
        header = [" ".join(col.split()[1:]) for col in header]
        header[0] = "Contig"
        with open(output[0], "w") as out:
            out.write("Sample\t" + "\t".join(header) + "\n")
            for sample, cov_file in zip(SAMPLE_NAMES, input):
                with open(cov_file) as f:
                    next(f)  # skip this sample's own header line
                    for line in f:
                        out.write(f"{sample}\t{line}")


rule run_combine_cov:
    """Sample\tContig\tCount\tRPKM\tTPM\tMean\tCovered_fraction\tVariance\n"""
    input:
        os.path.join(OUTDIR, "preprocess", "results", "sample_coverm_coverage.tsv")
    output:
        os.path.join(OUTDIR, "preprocess", "coverage.tsv")
    shell:
        # NR>1 skips the header inside awk itself, rather than the old `sed -i '1d'
        # {input}` which mutated the input in place. That made the rule non-
        # idempotent: a retry after partial failure would see an already-header-
        # stripped input, silently drop its first real data row as if it were still
        # the header, and produce a wrong-but-plausible coverage.tsv rather than
        # erroring. This version never touches {input} at all.
        """
        awk -F '\t' 'NR>1 {{ sum[$2] += $6 }} END {{ for (key in sum) print key, sum[key] }}' {input} > {output}
        """
