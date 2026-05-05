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

Uniform cohort:

```bash
make-vcf stage1 --input-dir /data/in --output-dir /data/out
```

Classed cohort by prefixes:

```bash
make-vcf stage1 --input-dir /data/in --output-dir /data/out --prefixes SSC,ABC
```

Outputs:
- `stage1_patient_ids_by_class.json`
- `stage1_<class>_patient_ids.txt`

## Stage 2

Use stage 1 class map, one default suffix for all classes:

```bash
make-vcf stage2 \
  --input-dir /data/in \
  --output-dir /data/out/stage2 \
  --class-map-json /data/out/stage1_patient_ids_by_class.json \
  --suffix vcf.gz \
  --batch-size 1000
```

Run only selected classes in stage 2:

```bash
make-vcf stage2 \
  --input-dir /data/in \
  --output-dir /data/out/stage2 \
  --class-map-json /data/out/stage1_patient_ids_by_class.json \
  --suffix vcf.gz \
  --classes-to-process SSC,ABC
```

Per-class suffix override:

```bash
make-vcf stage2 \
  --input-dir /data/in \
  --output-dir /data/out/stage2 \
  --class-map-json /data/out/stage1_patient_ids_by_class.json \
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

Batch outputs:
- `class_<name>/batch_00001_dVars.pkl`
- `class_<name>/batch_00001_dVars_inh.pkl`

Optional patient-ID sidecar lists (only written when that batch or class has at least one ID):

- `--task denovo` (default): `batch_*_patients_nonempty_dvars_inh.txt` per batch (if non-empty), and `patients_nonempty_dvars_inh.txt` for the class (if non-empty).
- `--task inherited`: `batch_*_patients_nonempty_dvars.txt` per batch (if non-empty), and `patients_nonempty_dvars.txt` for the class (if non-empty).

## Stage 3

Use class map from stage 1 and stage 2 batch outputs.  
`--output-dir` is the **pipeline root**; results go under a stage3 directory chosen by `--task` and optional `--classes-to-process`:

- **`--task denovo` (default):** merge `batch_*_dVars.pkl` → `<output-dir>/stage3/` or `<output-dir>/stage3_<classes>/`
- **`--task inherited`:** merge `batch_*_dVars_inh.pkl` → `<output-dir>/stage3_inherited/` or `<output-dir>/stage3_inherited_<classes>/`

Within each directory, per-chromosome `variants_snv_nohead.vcf` files are written the same way as for denovo.

```bash
make-vcf stage3 \
  --stage2-dir /data/out/stage2 \
  --output-dir /data/out \
  --class-map-json /data/out/stage1_patient_ids_by_class.json \
  --denovo-cap 100 \
  --class-cap SSC=80 \
  --class-cap ABC=120 \
  --task denovo
```

Inherited variants from stage 2 `dVars_inh` batches (example):

```bash
make-vcf stage3 \
  --stage2-dir /data/out/stage2 \
  --output-dir /data/out \
  --class-map-json /data/out/stage1_patient_ids_by_class.json \
  --denovo-cap 100 \
  --task inherited
```

Process only some classes (e.g. `/data/out/stage3_SSC_ABC/` with `--task denovo`, or `/data/out/stage3_inherited_SSC_ABC/` with `--task inherited`):

```bash
make-vcf stage3 \
  --stage2-dir /data/out/stage2 \
  --output-dir /data/out \
  --class-map-json /data/out/stage1_patient_ids_by_class.json \
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
