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

## Stage 3

Use class map from stage 1 and stage 2 batch outputs:

```bash
make-vcf stage3 \
  --stage2-dir /data/out/stage2 \
  --output-dir /data/out/stage3 \
  --class-map-json /data/out/stage1_patient_ids_by_class.json \
  --denovo-cap 100 \
  --class-cap SSC=80 \
  --class-cap ABC=120
```

You can also provide caps as one comma-separated list:

```bash
--class-cap SSC=80,ABC=120
```

Outputs:
- `<output_dir>/chr1/variants_snv_nohead.vcf` (and `chr2` ... `chr22`, `chrX`)
- `<output_dir>/stage3_summary.json`

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
  --class-cap SSC=80
```

Run only selected stages by providing only the needed flags:
- `--run-stage1`
- `--run-stage2`
- `--run-stage3`
