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

## Stage 1

`--output-dir` is the **pipeline root**; stage 1 writes under `<output-dir>/stage1/`.

Uniform cohort:

```bash
make-vcf stage1 --input-dir /data/in --output-dir /data/out
```

Classed cohort by prefixes:

```bash
make-vcf stage1 --input-dir /data/in --output-dir /data/out --prefixes SSC,ABC
```

Outputs (under `<output-dir>/stage1/`):
- `stage1_patient_ids_by_class.json` (all patients by class)
- `stage1_patient_ids_by_class_filtered_<cap_min>_<cap_max>.json` (default caps 22000–75000; used by stage 2/3)
- `stage1_<class>_patient_ids.txt`

Size filtering (`--cap-min`, `--cap-max`, defaults 22000 and 75000): when `--suffixes` is set, keep a patient only if **every** listed suffix file exists and its byte size is within `[cap_min, cap_max]`. Without `--suffixes`, the filtered file matches the full class map.

Optional file statistics (when `--suffixes` is provided):

- `--stats size` (default): byte size of each file
- `--stats counts`: line count (supports `.gz`)

Paths checked: `input_dir/<patient_id>/<patient_id>.<suffix>` for each suffix.

One JSON per class prefix and suffix for `--stats` (`stage1_stats_size_...` or `stage1_stats_counts_...`) and always for modification time (`stage1_stats_mtime_SSC_final_vcf_gz.json`). Mtime values are Unix epoch seconds (`st_mtime`). Only patients with existing files are included.

Missing files: patient IDs with no file at the expected path are listed in `{prefix}_{suffix_label}.missed.txt`, e.g. `SSC_final_vcf_gz.missed.txt` (always written; empty if none missing).

Files smaller than `--cap-min`: patient IDs with an existing file below the minimum size cap are listed in `{prefix}_{suffix_label}.small_{cap_min}.txt`, e.g. `SSC_final_vcf_gz.small_22000.txt` (always written; empty if none).

Example:

```bash
make-vcf stage1 \
  --input-dir /data/in \
  --output-dir /data/out \
  --prefixes SSC,SP \
  --suffixes .final.vcf.gz,.denovo.final.vcf.gz \
  --stats counts
```

## Stage 2

`--output-dir` is the **pipeline root** (writes under `<output-dir>/stage2/`).  
`--class-map-json` is **optional**; default is `<output-dir>/stage1/stage1_patient_ids_by_class_filtered_<cap_min>_<cap_max>.json` (same `--cap-min` / `--cap-max` as stage 1, defaults 22000 / 75000).

Use stage 1 class map, one default suffix for all classes:

```bash
make-vcf stage2 \
  --input-dir /data/in \
  --output-dir /data/out \
  --suffix vcf.gz \
  --batch-size 1000
```

Run only selected classes in stage 2:

```bash
make-vcf stage2 \
  --input-dir /data/in \
  --output-dir /data/out \
  --suffix vcf.gz \
  --classes-to-process SSC,ABC
```

Per-class suffix override:

```bash
make-vcf stage2 \
  --input-dir /data/in \
  --output-dir /data/out \
  --suffix vcf.gz \
  --class-suffix SSC=vcf.gz \
  --class-suffix ABC=vcf \
  --batch-size 1000
```

You can also provide a comma-separated list in one argument:

```bash
--class-suffix SSC=vcf.gz,ABC=vcf
```

Expected file pattern:
- `input_dir/patient_id/patient_id.<suffix>`

Batch outputs (depends on `--task` and save flags):
- `class_<name>/batch_00001_dVars.pkl`
- `class_<name>/batch_00001_dVars_inh.pkl`
- `class_<name>/patient_dVars_counts.json` (variant records per patient for the class, written once after all batches)
- `class_<name>/patient_dVars_inh_counts.json` (when inherited variants are collected)

Optional patient-ID sidecar lists (only written when that batch or class has at least one ID):

- `--task denovo` (default): `batch_*_patients_nonempty_dvars_inh.txt` per batch (if non-empty), and `patients_nonempty_dvars_inh.txt` for the class (if non-empty).
- `--task inherited`: `batch_*_patients_nonempty_dvars.txt` per batch (if non-empty), and `patients_nonempty_dvars.txt` for the class (if non-empty).

Collection/saving defaults in stage 2:

- `--task denovo`: collect/save `dVars` only. Use `--save-inh` to also collect/save `dVars_inh`.
- `--task inherited`: collect/save `dVars_inh` only. Use `--save-denovo` (or `--save_denovo`) to also collect/save `dVars`.

Denovo classification (stage 2):

- Default: `if_denovo()` (child alt count > 0, both parents 0).
- `--use-ext-denovo`: use `if_denovo_ext()` instead (your extended rules on the same child/mother/father alt counts).

## Stage 3

Use class map from stage 1 and stage 2 batch outputs.  
`--output-dir` is the **pipeline root**; results go under a stage3 directory chosen by `--task` and optional `--classes-to-process`:

- **`--task denovo` (default):** merge `batch_*_dVars.pkl` → `<output-dir>/stage3/` or `<output-dir>/stage3_<classes>/`
- **`--task inherited`:** merge `batch_*_dVars_inh.pkl` → `<output-dir>/stage3_inherited/` or `<output-dir>/stage3_inherited_<classes>/`

`--class-map-json` defaults to `<output-dir>/stage1/stage1_patient_ids_by_class_filtered_<cap_min>_<cap_max>.json`.  
`--stage2-dir` defaults to `<output-dir>/stage2`.  
Use the same `--cap-min` / `--cap-max` as stage 1 (defaults 22000 / 75000).

Within each directory, per-chromosome `variants_snv_nohead.vcf` files are written the same way as for denovo.

```bash
make-vcf stage3 \
  --output-dir /data/out \
  --denovo-cap 100 \
  --class-cap SSC=80 \
  --class-cap ABC=120 \
  --task denovo
```

Inherited variants from stage 2 `dVars_inh` batches (example):

```bash
make-vcf stage3 \
  --output-dir /data/out \
  --denovo-cap 100 \
  --task inherited
```

Process only some classes (e.g. `/data/out/stage3_SSC_ABC/` with `--task denovo`, or `/data/out/stage3_inherited_SSC_ABC/` with `--task inherited`):

```bash
make-vcf stage3 \
  --output-dir /data/out \
  --denovo-cap 100 \
  --classes-to-process SSC,ABC
```

You can also provide caps as one comma-separated list:

```bash
--class-cap SSC=80,ABC=120
```

Outputs (under the chosen stage3 directory):
- `chr1/variants_snv_nohead.vcf` (and `chr2` ... `chr22`, `chrX`)
- `stage3_summary.json`
- `stage3_titv.json` (transitions, transversions, and Ti/Tv ratio per chromosome and overall)

## Run all stages

```bash
make-vcf run123 \
  --input-dir /data/in \
  --output-dir /data/out \
  --prefixes SSC,ABC \
  --suffix vcf.gz \
  --class-suffix ABC=vcf \
  --batch-size 1000 \
  --denovo-cap 100 \
  --class-cap SSC=80
```

Optional: restrict stage 3 to a subset of classes (writes `/data/out/stage3_SSC_ABC/`). Add to the `run123` command:

```bash
--classes-to-process SSC,ABC
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
  --prefixes SSC,ABC \
  --suffix vcf.gz \
  --class-suffix ABC=vcf \
  --batch-size 1000 \
  --denovo-cap 100 \
  --class-cap SSC=80 \
  --classes-to-process SSC,ABC
```

Run only selected stages by providing only the needed flags:
- `--run-stage1`
- `--run-stage2`
- `--run-stage3`
