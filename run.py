from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List

# Allow running this script without installing the package.
ROOT_DIR = Path(__file__).resolve().parent
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from make_vcf_pipeline.stage1 import run_stage1
from make_vcf_pipeline.stage2 import parse_class_suffix_pairs, run_stage2
from make_vcf_pipeline.stage3 import parse_class_cap_pairs, run_stage3


def _split_csv(value: str) -> List[str]:
    if not value.strip():
        return []
    return [x.strip() for x in value.split(",") if x.strip()]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run stage1/stage2/stage3 one by one with explicit parameters."
    )

    # Which stages to execute, in fixed order (1 -> 2 -> 3)
    parser.add_argument("--run-stage1", action="store_true", help="Execute stage 1")
    parser.add_argument("--run-stage2", action="store_true", help="Execute stage 2")
    parser.add_argument("--run-stage3", action="store_true", help="Execute stage 3")

    # Shared paths
    parser.add_argument("--input-dir", type=Path, required=True, help="Input directory with patient subdirs")
    parser.add_argument("--output-dir", type=Path, required=True, help="Root output directory")

    # Stage 1 parameters
    parser.add_argument(
        "--prefixes",
        type=str,
        default="",
        help="Stage1: comma-separated patient class prefixes (empty means one class).",
    )

    # Stage 2 parameters
    parser.add_argument(
        "--suffix",
        type=str,
        default=None,
        help="Stage2: default file suffix (e.g. vcf.gz).",
    )
    parser.add_argument(
        "--class-suffix",
        action="append",
        default=[],
        help="Stage2: per-class suffix override, repeat as class=suffix.",
    )
    parser.add_argument("--batch-size", type=int, default=1000, help="Stage2: patient batch size")

    # Stage 3 parameters
    parser.add_argument(
        "--denovo-cap",
        type=int,
        default=None,
        help="Stage3: default cap of denovo variants per patient.",
    )
    parser.add_argument(
        "--class-cap",
        action="append",
        default=[],
        help="Stage3: per-class denovo cap override, repeat as class=cap.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if not (args.run_stage1 or args.run_stage2 or args.run_stage3):
        raise ValueError("Select at least one stage flag: --run-stage1/--run-stage2/--run-stage3")

    stage1_out_dir = args.output_dir
    stage1_json_path = stage1_out_dir / "stage1_patient_ids_by_class.json"
    stage2_out_dir = args.output_dir / "stage2"
    stage3_out_dir = args.output_dir / "stage3"

    class_map = None

    if args.run_stage1:
        prefixes = _split_csv(args.prefixes)
        stage1_result = run_stage1(args.input_dir, stage1_out_dir, prefixes)
        class_map = stage1_result.classes_to_patients
        print(f"[stage1] done -> {stage1_result.output_json}")

    if args.run_stage2:
        if class_map is None:
            if not stage1_json_path.is_file():
                raise FileNotFoundError(
                    f"Stage1 class map not found at '{stage1_json_path}'. "
                    "Run with --run-stage1 first or create this file."
                )
            class_map = json.loads(stage1_json_path.read_text(encoding="utf-8"))

        class_to_suffix = parse_class_suffix_pairs(args.class_suffix)
        stage2_outputs = run_stage2(
            input_dir=args.input_dir,
            output_dir=stage2_out_dir,
            class_map=class_map,
            default_suffix=args.suffix,
            class_to_suffix=class_to_suffix,
            batch_size=args.batch_size,
        )
        print(f"[stage2] done -> {len(stage2_outputs)} batch outputs in {stage2_out_dir}")

    if args.run_stage3:
        if class_map is None:
            if not stage1_json_path.is_file():
                raise FileNotFoundError(
                    f"Stage1 class map not found at '{stage1_json_path}'. "
                    "Run with --run-stage1 first or create this file."
                )
            class_map = json.loads(stage1_json_path.read_text(encoding="utf-8"))

        class_to_cap = parse_class_cap_pairs(args.class_cap)
        stage3_result = run_stage3(
            stage2_dir=stage2_out_dir,
            output_dir=stage3_out_dir,
            class_map=class_map,
            default_cap=args.denovo_cap,
            class_to_cap=class_to_cap,
        )
        print(f"[stage3] done -> {stage3_result.output_dir}")


if __name__ == "__main__":
    main()
