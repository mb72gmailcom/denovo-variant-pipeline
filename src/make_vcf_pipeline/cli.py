from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List

from .stage1 import run_stage1
from .stage2 import (
    load_stage1_class_map,
    parse_class_suffix_pairs,
    parse_classes_to_process as parse_stage2_classes_to_process,
    run_stage2,
)
from .stage3 import parse_class_cap_pairs, parse_classes_to_process, resolve_stage3_output, run_stage3


def _split_csv(value: str) -> List[str]:
    if not value.strip():
        return []
    return [x.strip() for x in value.split(",") if x.strip()]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="make-vcf",
        description="Pipeline runner for stage1, stage2, and stage3.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    p1 = subparsers.add_parser("stage1", help="Create ordered patient ID lists by class")
    p1.add_argument("--input-dir", type=Path, required=True)
    p1.add_argument("--output-dir", type=Path, required=True)
    p1.add_argument(
        "--prefixes",
        type=str,
        default="",
        help="Comma-separated class prefixes. Empty means uniform cohort.",
    )

    p2 = subparsers.add_parser("stage2", help="Process patient files in batches")
    p2.add_argument("--input-dir", type=Path, required=True)
    p2.add_argument("--output-dir", type=Path, required=True)
    p2.add_argument(
        "--class-map-json",
        type=Path,
        required=True,
        help="Path to stage1_patient_ids_by_class.json",
    )
    p2.add_argument(
        "--suffix",
        type=str,
        default=None,
        help="Default file suffix (e.g. vcf.gz). Can be overridden by --class-suffix.",
    )
    p2.add_argument(
        "--class-suffix",
        action="append",
        default=[],
        help="Per-class suffix override. Repeat format class=suffix.",
    )
    p2.add_argument("--batch-size", type=int, default=1000)
    p2.add_argument(
        "--classes-to-process",
        action="append",
        default=[],
        help="Optional: only these classes in stage2 (comma-separated or repeat). Omit for all classes.",
    )
    p2.add_argument(
        "--task",
        choices=("denovo", "inherited"),
        default="denovo",
        help="denovo: list patients with non-empty inherited (sidecar txt only if non-empty). "
        "inherited: list patients with non-empty denovo (sidecar txt only if non-empty).",
    )
    p2.add_argument(
        "--save-inh",
        action="store_true",
        help="Stage2: with --task denovo, also collect/save dVars_inh.",
    )
    p2.add_argument(
        "--save-denovo",
        "--save_denovo",
        action="store_true",
        help="Stage2: with --task inherited, also collect/save dVars.",
    )
    p2.add_argument(
        "--use-ext-denovo",
        action="store_true",
        help="Stage2: classify denovo with if_denovo_ext() instead of if_denovo().",
    )

    p12 = subparsers.add_parser("run12", help="Run stage1 then stage2")
    p12.add_argument("--input-dir", type=Path, required=True)
    p12.add_argument("--output-dir", type=Path, required=True)
    p12.add_argument("--prefixes", type=str, default="")
    p12.add_argument("--suffix", type=str, default=None)
    p12.add_argument("--class-suffix", action="append", default=[])
    p12.add_argument("--batch-size", type=int, default=1000)
    p12.add_argument(
        "--classes-to-process",
        action="append",
        default=[],
        help="Optional: only these classes in stage2 (comma-separated or repeat). Omit for all classes.",
    )
    p12.add_argument(
        "--task",
        choices=("denovo", "inherited"),
        default="denovo",
        help="Stage2 task (same meaning as make-vcf stage2 --task).",
    )
    p12.add_argument("--save-inh", action="store_true", help="Run12 stage2: also save inherited while task=denovo.")
    p12.add_argument(
        "--save-denovo",
        "--save_denovo",
        action="store_true",
        help="Run12 stage2: also save denovo while task=inherited.",
    )
    p12.add_argument(
        "--use-ext-denovo",
        action="store_true",
        help="Run12 stage2: use extended denovo definition (if_denovo_ext).",
    )

    p3 = subparsers.add_parser("stage3", help="Post-process stage2 outputs into chromosome files")
    p3.add_argument("--stage2-dir", type=Path, required=True)
    p3.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Pipeline output root; stage3 writes under stage3/ or stage3_<classes>/ (denovo), or stage3_inherited/... (inherited).",
    )
    p3.add_argument("--class-map-json", type=Path, required=True)
    p3.add_argument(
        "--denovo-cap",
        type=int,
        default=None,
        help="Default denovo cap per patient.",
    )
    p3.add_argument(
        "--class-cap",
        action="append",
        default=[],
        help="Per-class denovo cap override. Repeat format class=cap.",
    )
    p3.add_argument(
        "--classes-to-process",
        action="append",
        default=[],
        help="Optional: only these classes (comma-separated or repeat). Omit for all classes.",
    )
    p3.add_argument(
        "--task",
        choices=("denovo", "inherited"),
        default="denovo",
        help="denovo: merge batch dVars. inherited: merge batch dVars_inh; output under stage3_inherited/...",
    )

    p123 = subparsers.add_parser("run123", help="Run stage1, stage2, then stage3")
    p123.add_argument("--input-dir", type=Path, required=True)
    p123.add_argument("--output-dir", type=Path, required=True)
    p123.add_argument("--prefixes", type=str, default="")
    p123.add_argument("--suffix", type=str, default=None)
    p123.add_argument("--class-suffix", action="append", default=[])
    p123.add_argument("--batch-size", type=int, default=1000)
    p123.add_argument("--denovo-cap", type=int, default=None)
    p123.add_argument("--class-cap", action="append", default=[])
    p123.add_argument(
        "--classes-to-process",
        action="append",
        default=[],
        help="Optional: only these classes in stage2 and stage3 (comma-separated or repeat). Omit for all.",
    )
    p123.add_argument(
        "--task",
        choices=("denovo", "inherited"),
        default="denovo",
        help="Applies to stage2 (sidecar lists) and stage3 (which pickle set and output dir prefix).",
    )
    p123.add_argument("--save-inh", action="store_true", help="Run123 stage2: also save inherited while task=denovo.")
    p123.add_argument(
        "--save-denovo",
        "--save_denovo",
        action="store_true",
        help="Run123 stage2: also save denovo while task=inherited.",
    )
    p123.add_argument(
        "--use-ext-denovo",
        action="store_true",
        help="Run123 stage2: use extended denovo definition (if_denovo_ext).",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "stage1":
        prefixes = _split_csv(args.prefixes)
        result = run_stage1(args.input_dir, args.output_dir, prefixes)
        print(f"Stage1 complete: {result.output_json}")
        return

    if args.command == "stage2":
        class_map = load_stage1_class_map(args.class_map_json)
        class_to_suffix = parse_class_suffix_pairs(args.class_suffix)
        classes_to_process = parse_stage2_classes_to_process(args.classes_to_process)
        outputs = run_stage2(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            class_map=class_map,
            default_suffix=args.suffix,
            class_to_suffix=class_to_suffix,
            batch_size=args.batch_size,
            task=args.task,
            classes_to_process=classes_to_process,
            save_inh=args.save_inh,
            save_denovo=args.save_denovo,
            use_ext_denovo=args.use_ext_denovo,
        )
        print(json.dumps([o.__dict__ for o in outputs], indent=2, default=str))
        return

    if args.command == "run12":
        prefixes = _split_csv(args.prefixes)
        stage1_result = run_stage1(args.input_dir, args.output_dir, prefixes)
        class_to_suffix = parse_class_suffix_pairs(args.class_suffix)
        classes_to_process = parse_stage2_classes_to_process(args.classes_to_process)
        outputs = run_stage2(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            class_map=stage1_result.classes_to_patients,
            default_suffix=args.suffix,
            class_to_suffix=class_to_suffix,
            batch_size=args.batch_size,
            task=args.task,
            classes_to_process=classes_to_process,
            save_inh=args.save_inh,
            save_denovo=args.save_denovo,
            use_ext_denovo=args.use_ext_denovo,
        )
        print(f"Stage1 file: {stage1_result.output_json}")
        print(json.dumps([o.__dict__ for o in outputs], indent=2, default=str))
        return

    if args.command == "stage3":
        class_map = load_stage1_class_map(args.class_map_json)
        class_to_cap = parse_class_cap_pairs(args.class_cap)
        requested = parse_classes_to_process(args.classes_to_process)
        class_map_run, stage3_dir = resolve_stage3_output(
            args.output_dir, class_map, requested, task=args.task
        )
        result = run_stage3(
            stage2_dir=args.stage2_dir,
            output_dir=stage3_dir,
            class_map=class_map_run,
            default_cap=args.denovo_cap,
            class_to_cap=class_to_cap,
            task=args.task,
        )
        print(json.dumps(result.__dict__, indent=2, default=str))
        return

    if args.command == "run123":
        prefixes = _split_csv(args.prefixes)
        stage1_result = run_stage1(args.input_dir, args.output_dir, prefixes)
        class_to_suffix = parse_class_suffix_pairs(args.class_suffix)
        stage2_dir = args.output_dir / "stage2"

        stage2_outputs = run_stage2(
            input_dir=args.input_dir,
            output_dir=stage2_dir,
            class_map=stage1_result.classes_to_patients,
            default_suffix=args.suffix,
            class_to_suffix=class_to_suffix,
            batch_size=args.batch_size,
            task=args.task,
            classes_to_process=parse_stage2_classes_to_process(args.classes_to_process),
            save_inh=args.save_inh,
            save_denovo=args.save_denovo,
            use_ext_denovo=args.use_ext_denovo,
        )
        class_to_cap = parse_class_cap_pairs(args.class_cap)
        requested = parse_classes_to_process(args.classes_to_process)
        class_map_run, stage3_dir = resolve_stage3_output(
            args.output_dir, stage1_result.classes_to_patients, requested, task=args.task
        )
        stage3_result = run_stage3(
            stage2_dir=stage2_dir,
            output_dir=stage3_dir,
            class_map=class_map_run,
            default_cap=args.denovo_cap,
            class_to_cap=class_to_cap,
            task=args.task,
        )
        print(f"Stage1 file: {stage1_result.output_json}")
        print(json.dumps([o.__dict__ for o in stage2_outputs], indent=2, default=str))
        print(json.dumps(stage3_result.__dict__, indent=2, default=str))
        return


if __name__ == "__main__":
    main()
