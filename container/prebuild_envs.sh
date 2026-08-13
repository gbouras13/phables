#!/bin/bash
# Pre-builds EVERY per-rule conda env into the image, at Docker build time.
# Run from container/Dockerfile; not meant to be run on a host.
#
# Envs go into phables' own DEFAULT --conda-prefix (snake_base("workflow/conda"),
# i.e. inside the installed package) deliberately: a user inside the resulting
# container then runs plain `phables run ...` with no --conda-prefix flag and
# Snakemake resolves these exact envs. Snakemake names an env by hashing its
# file content together with the conda prefix path -- both identical at build
# time and run time here, since they're the same paths in the same image -- so
# the hashes match and nothing is rebuilt at runtime. That matters more than
# it sounds: a .sif is READ-ONLY when running, so an attempted rebuild is a
# hard failure, not just a slow path.
#
# --conda-create-envs-only builds a DAG's envs without running any of it. The
# DAG still has to RESOLVE, which requires the database files to EXIST -- but
# only to exist, since no job runs. Empty placeholder files are therefore
# enough. That was verified for real (dry-running every flag combination below
# against zero-byte placeholder DB files) before this script was written, and
# it's what keeps the multi-GB PHROGs/hallmark databases OUT of the image:
# mount the real ones at runtime via --databases, exactly as outside a
# container.
#
# One invocation per flag combination, because which envs a DAG needs depends
# on the flags -- gene caller, phage-detection mode and GPU backend each select
# different rules and env files. Together these cover every env under
# workflow/envs/ that any `phables run` (or `phables install`) can reach.

set -euxo pipefail

SRC="${1:-/opt/phables_src}"
DB=/tmp/placeholder_db

mkdir -p "$DB/phrogs_mmseqs_db" "$DB/hallmark_db"
touch "$DB/marker.hmm" \
      "$DB/phrog_annot_v4.tsv" \
      "$DB/phrogs_mmseqs_db/phrogs_profile_db" \
      "$DB/hallmark_db/hallmark_db" \
      "$DB/hallmark_db/hallmark_categories.tsv"

GFA="${SRC}/tests/data/ERR1301161/assembly_graph_after_simplification.gfa"
READS="${SRC}/tests/data/ERR1301161"
COMMON=(--input "$GFA" --reads "$READS" --databases "$DB" --threads 1)

# 1. default path -> coverm, genecall (FragGeneScan), smg (HMMER), mmseqs, phables
phables run "${COMMON[@]}" --output /tmp/envbuild1 --conda-create-envs-only

# 2. the other gene caller
phables run "${COMMON[@]}" --output /tmp/envbuild2 \
    --genecaller pyrodigal-gv --conda-create-envs-only

# 3. ProstT5 + foldseek detection on ROCm -- the reason this image exists
phables run "${COMMON[@]}" --output /tmp/envbuild3 \
    --phagedetection prostt5-foldseek --gpu-backend rocm --conda-create-envs-only

# 4. Same, CPU backend. Built because gpu_backend's own config default is `cpu`:
#    without this, forgetting `--gpu-backend rocm` inside the container would hit
#    a missing env on a read-only filesystem. (cuda is deliberately NOT built --
#    this is a ROCm image for Setonix, and the CUDA wheels would add several GB
#    that could never be used on this hardware.)
phables run "${COMMON[@]}" --output /tmp/envbuild4 \
    --phagedetection prostt5-foldseek --gpu-backend cpu --conda-create-envs-only

# 5. Optional phylogenetic tree -> phylotree env (MAFFT + cogent3/piqtree).
#    Note conda resolves cogent3 fine, unlike pip, where every published
#    release is a prerelease and a plain version range matches nothing.
phables run "${COMMON[@]}" --output /tmp/envbuild5 \
    --build-tree --conda-create-envs-only

# 6. install.smk's own env (curl), so `phables install` works inside here too
phables install --output /tmp/envbuild6 --databases "$DB" --conda-create-envs-only

rm -rf /tmp/envbuild1 /tmp/envbuild2 /tmp/envbuild3 /tmp/envbuild4 /tmp/envbuild5 /tmp/envbuild6 "$DB"
conda clean -a -y

echo "=== pre-built conda envs ==="
CONDA_PREFIX_DIR="$(python -c 'import phables, os; print(os.path.join(os.path.dirname(phables.__file__), "workflow", "conda"))')"
ls -1 "$CONDA_PREFIX_DIR"
