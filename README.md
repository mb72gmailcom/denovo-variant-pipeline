# make-vcf-pipeline

Packaging of a 3-stage VCF workflow.

Currently implemented:
- Stage 1: classify and save ordered patient ID lists by class.
- Stage 2: process patient files in batches and save `dVars` / `dVars_inh` per batch.
- Stage 3: combine stage 2 denovo outputs, apply patient denovo-cap filtering, split/sort by chromosome, and write VCF-like outputs.

## Install (editable)

```bash
pip install -e .
```

Each stage writes a parameters JSON in its output directory:

| Stage | Path | Used by later stages for |
|-------|------|--------------------------|
| 1 | `<output-dir>/stage1/stage1_parameters.json` | Input dir, class map caps, suffixes per class |
| 2 | `<output-dir>/stage2/stage2_parameters.json` | Collect bucket, classes, suffixes, batch settings |
| 3 | `<output-dir>/stage3_vN/stage3_parameters.json` | Caps, SNV/filter settings for the run |

## Stage 1

`--output-dir` is the **pipeline root**; stage 1 writes under `<output-dir>/stage1/`.

Uniform cohort:

```bash
make-vcf stage1 --input-dir /data/in --output-dir /data/out
```

Classed cohort by patient ID prefix:

```bash
make-vcf stage1 --input-dir /data/in --output-dir /data/out --classes SSC,ABC
```

Outputs under `<output-dir>/stage1/`:

**Run root** (three files only):
- `stage1_patient_ids_by_class.json` (all patients by class)
- `stage1_patient_ids_by_class_filtered_<cap_min>_<cap_max>.json` (default caps 22000–75000; used by stage 2/3)
- `stage1_parameters.json` (run parameters: classes, caps, suffixes per class; read by stage 2/3 when CLI flags are omitted)

**Per class** (`class_<name>/`):
- `patient_ids.txt` — ordered patient IDs for that class
- When `--suffix` or `--class-suffix` is set (size filtering enabled):
  - `stats_<size|counts>.json` — per-patient file statistic
  - `stats_mtime.json` — per-patient modification time (Unix epoch seconds)
  - `missed.txt` — patient IDs with no file at the expected path
  - `small_<cap_min>.txt` — patient IDs with file size below `--cap-min`

Size filtering (`--cap-min`, `--cap-max`, defaults 22000 and 75000): when `--suffix` or `--class-suffix` is set, keep a patient only if that class's suffix file exists and its byte size is within `[cap_min, cap_max]`. Without a suffix, the filtered file matches the full class map.

Optional file statistics (when `--suffix` or `--class-suffix` is provided):

- `--stats size` (default): byte size of each file
- `--stats counts`: line count (supports `.gz`)

Paths checked: `input_dir/<patient_id>/<patient_id>.<suffix>` using the class default (`--suffix`) or per-class override (`--class-suffix`).

Only patients with existing files are included in the stats JSON files.

Example:

```bash
make-vcf stage1 \
  --input-dir /data/in \
  --output-dir /data/out \
  --classes SSC,SP \
  --class-suffix SSC=.final.vcf.gz \
  --class-suffix SP=.denovo.final.vcf.gz \
  --stats counts
```

Default suffix for all classes, with one override:

```bash
make-vcf stage1 \
  --input-dir /data/in \
  --output-dir /data/out \
  --classes SSC,SP \
  --suffix .final.vcf.gz \
  --class-suffix SP=.denovo.final.vcf.gz
```

Comma-separated pairs in one flag (same as stage 2):

```bash
--class-suffix SSC=.final.vcf.gz,SP=.denovo.final.vcf.gz
```

`stage1_parameters.json` records the full run configuration, including resolved suffixes per class:

```json
{
  "version": 1,
  "input_dir": "/data/in",
  "pipeline_root": "/data/out",
  "classes": ["SSC", "SP"],
  "stats": "size",
  "cap_min": 22000,
  "cap_max": 75000,
  "default_suffix": ".final.vcf.gz",
  "class_suffix_overrides": { "SP": ".denovo.final.vcf.gz" },
  "suffix_per_class": {
    "SSC": ".final.vcf.gz",
    "SP": ".denovo.final.vcf.gz",
    "unmatched": ".final.vcf.gz"
  },
  "suffix_filtering_enabled": true,
  "class_names": ["SSC", "SP", "unmatched"],
  "class_map": "stage1_patient_ids_by_class.json",
  "filtered_class_map": "stage1_patient_ids_by_class_filtered_22000_75000.json"
}
```

Stage 2 and stage 3 load this file automatically (when present) to pick the filtered class map caps; stage 2 also uses `suffix_per_class` and `class_names` from it.

## Stage 2

`--output-dir` is the **pipeline root** (writes under `<output-dir>/stage2/`).  
`--class-map-json` is **optional**; default is `<output-dir>/stage1/stage1_patient_ids_by_class_filtered_<cap_min>_<cap_max>.json`. When `stage1_parameters.json` exists, stage 2/3 use the caps recorded there for that default path.

Stage 2 reads **`input_dir`**, **suffixes**, and the default **class list** from `<output-dir>/stage1/stage1_parameters.json` (written by stage 1). You do not pass `--input-dir`, `--suffix`, or `--class-suffix` to stage 2.

Use `--classes` to process only a subset (comma-separated class names, e.g. `SSC,ABC`). If omitted, all classes listed in the stage 1 parameters file (present in the filtered class map) are processed.

```bash
make-vcf stage2 \
  --output-dir /data/out \
  --batch-size 1000
```

Process only selected classes:

```bash
make-vcf stage2 \
  --output-dir /data/out \
  --classes SSC,ABC
```

Expected file pattern (paths from stage 1 parameters: `input_dir` and suffix per class):
- `input_dir/patient_id/patient_id.<suffix>`

Batch outputs (depends on `--collect` and `--save-all`):
- `stage2_parameters.json` (run parameters: collect, save_all, suffixes per class, batch size, etc.)
- `class_<name>/batch_00001_dVars.pkl`
- `class_<name>/batch_00001_dVars_inh.pkl`
- `class_<name>/patient_dVars_counts.json` (variant records per patient for the class, written once after all batches)
- `class_<name>/patient_dVars_inh_counts.json` (when inherited variants are collected)

Histograms (default on; use `--stage2-nostats` to skip; primary `--collect` bucket only; AB files omitted if all bins are zero):

- `children_qt_hist.json`, `parents_qt_hist.json`, `children_dp_hist.json`, `parents_dp_hist.json` (all genotypes pooled)
- `children_00_ab_hist.json`, `parents_00_ab_hist.json`, `children_01_ab_hist.json`, `parents_01_ab_hist.json`, `children_11_ab_hist.json`, `parents_11_ab_hist.json` (AB = `ada/(adr+ada)`)

Optional patient-ID sidecar lists (only written when that batch or class has at least one ID):

- `--collect denovo` (default): `batch_*_patients_nonempty_dvars_inh.txt` per batch (if non-empty), and `patients_nonempty_dvars_inh.txt` for the class (if non-empty).
- `--collect inherited`: `batch_*_patients_nonempty_dvars.txt` per batch (if non-empty), and `patients_nonempty_dvars.txt` for the class (if non-empty).

Collection/saving defaults in stage 2:

- `--collect denovo` (default): collect/save `dVars` only. Use `--save-all` to also collect/save `dVars_inh`.
- `--collect inherited`: collect/save `dVars_inh` only. Use `--save-all` to also collect/save `dVars`.

Denovo classification (stage 2):

- Default: `if_denovo()` (child alt count > 0, both parents 0).
- `--use-ext-denovo`: use `if_denovo_ext()` instead (your extended rules on the same child/mother/father alt counts).

## Stage 3

Use class map from stage 1 and stage 2 batch outputs.  
`--output-dir` is the **pipeline root**; results go under a stage3 directory chosen by **`--collect`** (default: value from `stage2_parameters.json`) and optional **`--stage3-classes`**:

- **`collect=denovo`:** read `batch_*_dVars.pkl` per class → `<output-dir>/stage3_vN/`
- **`collect=inherited`:** read `batch_*_dVars_inh.pkl` per class → `<output-dir>/stage3_inherited_vN/`

Standalone `make-vcf stage3` reads **`collect`** and **classes** from `stage2/stage2_parameters.json` when those flags are omitted. Pass **`--collect`** or **`--stage3-classes`** to override. Class subset does not change the output directory name; each run gets the next `stage3_vN` (or `stage3_inherited_vN`).

Each stage3 run uses the next monotonic version `N` for that stem (e.g. `stage3_v0`, then `stage3_v1`; deleted intermediate versions are not reused). Variants are **not** merged across classes; each class is written to its own subdirectory with separate summary files.

`--class-map-json` defaults to `<output-dir>/stage1/stage1_patient_ids_by_class_filtered_<cap_min>_<cap_max>.json`.  
`--stage2-dir` defaults to `<output-dir>/stage2`.  
Use the same `--cap-min` / `--cap-max` as stage 1 (defaults 22000 / 75000).

Within each `class_<name>/` directory, per-chromosome `variants.vcf` files are written (headerless; variant set depends on `--snv`, `--autosomal`, and other stage3 filters). By default, only autosomal chromosomes (chr1–chr22) are included; use `--no-autosomal` to also include chrX, chrY, and chrM.

```bash
make-vcf stage3 \
  --output-dir /data/out \
  --denovo-cap 100 \
  --class-cap SSC=80 \
  --class-cap ABC=120 \
  --collect denovo
```

Inherited variants from stage 2 `dVars_inh` batches (example):

```bash
make-vcf stage3 \
  --output-dir /data/out \
  --denovo-cap 100 \
  --collect inherited
```

Process only some classes (each run still writes to the next `stage3_vN`; classes processed are recorded in `stage3_parameters.json`):

```bash
make-vcf stage3 \
  --output-dir /data/out \
  --denovo-cap 100 \
  --stage3-classes SSC,ABC
```

You can also provide caps as one comma-separated list:

```bash
--class-cap SSC=80,ABC=120
```

Stage3 filter order (recorded in each class `stage3_summary.json`):

1. **Per-patient variant cap** — `--denovo-cap` / `--class-cap` (excludes patients with too many variants)
2. **Chromosome allowlist** — chr1–chr22 by default (`--autosomal`); chr1–chr22 + chrX/chrY/chrM with `--no-autosomal`
3. **SNV-only** — `--snv` (default on; use `--no-snv` to keep indels/multi-base alleles)
4. **Genotype QC** — `--filter` (optional; `is_good` on mother, father, and child)

Defaults when `--filter` is set: `--filter-dp 20`, `--filter-qt 90`, `--filter-abHom0 0.05`, `--filter-abHom1 0.05`, `--filter-abHet 0.30`. Depth and quality pass when `dp >= dp_cap` and `qt >= qt_cap`.

```bash
make-vcf stage3 \
  --output-dir /data/out \
  --denovo-cap 100 \
  --filter
```

Override any threshold explicitly, e.g. `--filter-dp 15 --filter-abHet 0.25`.

Without `--filter`, stage3 behavior is unchanged (cap only).

Outputs (under the chosen stage3 run directory, one set per class):

```
stage3_vN/
  stage3_parameters.json  # run parameters: collect, autosomal, chromosomes, caps, snv/filter settings
  class_SSC/
    patient_variant_counts.json   # variants per patient (class-level; use --stage3-nostats to skip)
    variant_patient_counts.json   # patients per variant key (class-level; use --stage3-nostats to skip)
    chr1/
      variants.vcf
      variant_patients.json       # patient IDs per variant on chr1 (always written)
    chr2/
      variants.vcf
      variant_patients.json
    ...
    stage3_summary.json   # cap, snv_only, filter_caps for this class
    stage3_titv.json      # Ti/Tv per chromosome and overall for this class
  class_ABC/
    ...
```

Each `stage3_summary.json` includes `patient_variant_cap` (the cap applied to that class), `snv_only`, `filter_caps` when `--filter` is used, and a `filter_pipeline` array recording variant-key counts after each filter step (original → patient cap → SNV if enabled → genotype QC if enabled → any extra filters).

## Run all stages

Combined commands (`run12`, `run123`, `run.py`) pass stage 1 settings to later stages automatically:

- **Stage 2** uses suffixes and the default class list from stage 1 (in memory when stage 1 ran in the same command, otherwise from `stage1_parameters.json`).
- **Stage 3** uses **`collect`** and the class list from stage 2 (in memory when stage 2 ran in the same command, otherwise from `stage2_parameters.json`), unless `--collect` or `--stage3-classes` overrides.

Optional: restrict stage 2 with `--stage2-classes`; stage 3 will match unless `--stage3-classes` is set.

```bash
make-vcf run123 \
  --input-dir /data/in \
  --output-dir /data/out \
  --classes SSC,ABC \
  --suffix vcf.gz \
  --class-suffix ABC=vcf \
  --batch-size 1000 \
  --denovo-cap 100 \
  --class-cap SSC=80
```

(`--classes` applies to stage 1; `--suffix` / `--class-suffix` apply to stage 1 only. Stage 2 and stage 3 pick up classes and suffixes from the stage 1 parameters file automatically in the same command. Use `--stage2-classes` in `run12`/`run123`/`run.py` to limit stage 2 to a subset; stage 3 then follows stage 2's classes unless `--stage3-classes` overrides.)

Optional: restrict stage 3 to a subset of classes with `--stage3-classes`. Output still goes to the next `stage3_vN` under the pipeline root. Add to the `run123` command:

```bash
--stage3-classes SSC,ABC
```

## `run.py` helper

You can also run stages one by one with explicit parameters:

```bash
python run.py \
  --run-stage1 \
  --run-stage2 \
  --run-stage3 \
  --input-dir /data/in \
  --output-dir /data/out \
  --classes SSC,ABC \
  --suffix vcf.gz \
  --class-suffix ABC=vcf \
  --batch-size 1000 \
  --denovo-cap 100 \
  --class-cap SSC=80 \
  --stage3-classes SSC,ABC
```

Run only selected stages by providing only the needed flags:
- `--run-stage1` requires `--input-dir`
- `--run-stage2` and `--run-stage3` do not use `--input-dir` (read from `stage1_parameters.json` and stage2 outputs under `--output-dir`)

Stage3-only example:

```bash
python run.py --run-stage3 --output-dir /data/out --denovo-cap 100
```
