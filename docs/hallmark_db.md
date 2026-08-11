# Building the hallmark structure database

`--phagedetection prostt5-foldseek` needs two inputs that aren't downloaded by
`phables install`: a Foldseek structure subDB of phage hallmark proteins
(`--hallmark-db`) and its matching PHROG-category table
(`--hallmark-categories`). This page explains what they are, how
`build_hallmark_db.py` builds them, and why they're a separate, one-time admin
step rather than something `phables install` fetches automatically.

## What "hallmark" means here

`phables_utils/component_utils.py`'s `get_components` screens phage-like
components for structural/virion evidence — specifically the PHROG functional
categories **head and packaging**, **connector**, **tail**, and **lysis**.
Those four categories are what "hallmark" refers to throughout this feature;
everything else in the PHROG catalogue (metabolism, DNA/RNA processing,
unknown function, ...) is irrelevant to this specific structural-evidence
check and is dropped.

**Integration and excision** (integrases, recombinases) is kept as a
*separate* id list/subDB rather than folded into the hallmark set.
Integrase-family structures are a real lysogeny signal, but integrases and
recombinases are abundant on bacterial chromosomes and mobile genetic elements
too — mixing them into the structural-hallmark evidence would make that
signal noisier, not stronger.

## Source data: phold's structure database

The subDB is built by subsetting the **full phold structure database**
([`gbouras13/phold`](https://github.com/gbouras13/phold)) — not by predicting
structures from scratch. Two files from that download are needed:

- `all_phold_structures` (the Foldseek structure DB itself, plus its
  `_ss`/`_h` companion files and `.lookup` index — ~1.36M entries, several GB)
- `phold_annots.tsv` (the PHROG annotation table mapping each structure to a
  PHROG id, product name, and functional category)

`phold_annots.tsv` uses **CRLF line endings** — every field comparison in
`build_hallmark_db.py` strips a trailing `\r` explicitly. Skipping this
silently matches zero rows instead of erroring, which is exactly the failure
mode that was hit building the reference DB below — worth knowing if you ever
modify this script.

## Why a representative cap, and why a stride sample

Some PHROG families have far more structures in phold's DB than others.
Keeping every single one would bloat the subDB for no real benefit (Foldseek
search quality saturates well before "every known example"), so
`--max-per-phrog` (default 30) caps each PHROG family's contribution.

The cap is applied via a **deterministic stride sample** over the family's
sorted `.lookup` row indices — not the first N. Row order in phold's `.lookup`
file tends to cluster by source dataset, so taking the first N would bias
toward whatever dataset happened to be indexed first for that PHROG, rather
than spreading the cap across the diversity that actually exists in the
family. A stride (`ordered[int(i * len(ordered)/max_per_phrog)]` for
`i in range(max_per_phrog)`) picks evenly across the whole sorted range
instead.

## Running it

```bash
mamba env create -f phables/workflow/envs/foldseek.yaml -n foldseek
conda activate foldseek

python phables/workflow/scripts/build_hallmark_db.py \
    --phold-db-prefix /path/to/all_phold_structures \
    --annots /path/to/phold_annots.tsv \
    --out-dir hallmark_db/ \
    --max-per-phrog 30
```

This writes, under `hallmark_db/`:

| File | What it is |
|---|---|
| `hallmark_db*` (`hallmark_db`, `hallmark_db_ss`, `hallmark_db_h` + `.index`/`.dbtype`) | The Foldseek subDB — pass its prefix (`hallmark_db/hallmark_db`) to `--hallmark-db` |
| `hallmark_categories.tsv` | `phrog_id\tcategory` — pass to `--hallmark-categories` |
| `hallmark_ids.tsv` | The `.lookup` row indices `createsubdb` was given (intermediate; not needed at run time) |
| `integration_excision_db*`, `integration_excision_categories.tsv`, `integration_excision_ids.tsv` | The separate integration/excision channel discussed above — not currently wired into `--hallmark-db` (phables only consumes the hallmark channel today), kept in case a future check wants it |

Pass `--skip-createsubdb` first if you just want to sanity-check the category
counts before committing to the `createsubdb` step, which reads through the
full multi-GB structure DB and is the slow part.

## Reference build — real numbers, so you can sanity-check your own

Built once against a real, complete phold structure DB download (not a
subsample) as part of validating this feature end-to-end:

- **5,389** hallmark PHROGs (head and packaging + connector + tail + lysis)
- **436** integration-and-excision PHROGs
- **49,869** hallmark structures after capping (`--max-per-phrog 30`)
- **2,610** integration-and-excision structures after capping
- **137 MB** total subDB size — comfortably inside "low hundreds of
  thousands of structures, page-cacheable" territory
- Verified **queryable**, not just built: a self-search of the hallmark subDB
  returned biologically sensible cross-hits (`phrog_2` ↔ `phrog_5653`, both
  independently annotated "terminase large subunit / head and packaging" in
  `phold_annots.tsv`)

If your own build's PHROG-category counts differ substantially from the first
two numbers above, that's worth investigating before trusting the result —
those two counts depend only on `phold_annots.tsv`'s categories, not on which
structures happen to be in your particular phold DB snapshot, so they should
be stable across phold DB versions.

## Why this isn't part of `phables install`

`phables install` downloads small, purpose-built assets (`marker.hmm`, the
PHROGs MMseqs profile DB) directly from their own stable URLs. The hallmark
subDB's source — phold's full structure database — is a multi-GB, separately
versioned/distributed asset with its own release cadence, not something this
project should silently re-host or auto-fetch a pinned copy of. Building the
hallmark subDB is a deliberate, one-time admin step: get phold's structure DB
yourself (see phold's own docs for the current download instructions), then
run `build_hallmark_db.py` once and reuse the output across every
`--phagedetection prostt5-foldseek` run.
