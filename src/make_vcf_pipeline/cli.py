from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List

from .stage1 import (
    DEFAULT_CAP_MAX,
    DEFAULT_CAP_MIN,
    infer_pipeline_root,
    load_class_map,
    parse_suffixes,
    resolve_class_map_json,
    resolve_stage2_input_dir,
    resolve_stage2_output_dir,
    run_stage1,
)
from .stage2 import (
    parse_class_suffix_pairs,
    parse_classes_to_process as parse_stage2_classes_to_process,
    run_stage2,
)
from .stage3 import parse_class_cap_pairs, parse_classes_to_process, resolve_stage3_output, run_stage3


def _split_csv(value: str) -> List[str]:
    if not value.strip():
        return []
    return [x.strip() for x in value.split(",") if x.strip()]


def _add_cap_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--cap-min",
        type=int,
        default=DEFAULT_CAP_MIN,
        help="Stage1 size filter lower bound (bytes); used for filtered class map filename and filtering.",
    )
    parser.add_argument(
        "--cap-max",
        type=int,
        default=DEFAULT_CAP_MAX,
        help="Stage1 size filter upper bound (bytes); used for filtered class map filename and filtering.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="make-vcf",
        description="Pipeline runner for stage1, stage2, and stage3.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    p1 = subparsers.add_parser("stage1", help="Create ordered patient ID lists by class")
    p1.add_argument("--input-dir", type=Path, required=True)
    p1.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Pipeline output root; stage1 writes under <output-dir>/stage1/.",
    )
    p1.add_argument(
        "--prefixes",
        type=str,
        default="",
        help="Comma-separated class prefixes. Empty means uniform cohort.",
    )
    p1.add_argument(
        "--stats",
        choices=("size", "counts"),
        default="size",
        help="Per-patient file statistic when --suffixes is set (default: size).",
    )
    p1.add_argument(
        "--suffixes",
        type=str,
        default="",
        help="Comma-separated file suffixes for stage1 stats, e.g. .vcf.gz,.final.vcf.gz",
    )
    _add_cap_args(p1)

    p2 = subparsers.add_parser("stage2", help="Process patient files in batches")
    p2.add_argument("--input-dir", type=Path, required=True)
    p2.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Pipeline output root or stage2 dir; writes under <root>/stage2/ when root is given.",
    )
    p2.add_argument(
        "--class-map-json",
        type=Path,
        default=None,
        help="Optional; default <pipeline-root>/stage1/stage1_patient_ids_by_class_filtered_<cap_min>_<cap_max>.json",
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
    _add_cap_args(p2)

    p12 = subparsers.add_parser("run12", help="Run stage1 then stage2")
    p12.add_argument("--input-dir", type=Path, required=True)
    p12.add_argument("--output-dir", type=Path, required=True)
    p12.add_argument("--prefixes", type=str, default="")
    p12.add_argument("--stats", choices=("size", "counts"), default="size")
    p12.add_argument("--suffixes", type=str, default="")
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
    _add_cap_args(p12)

    p3 = subparsers.add_parser("stage3", help="Post-process stage2 outputs into chromosome files")
    p3.add_argument(
        "--stage2-dir",
        type=Path,
        default=None,
        help="Optional; default <output-dir>/stage2",
    )
    p3.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Pipeline output root; stage3 writes under stage3/ or stage3_<classes>/ (denovo), or stage3_inherited/... (inherited).",
    )
    p3.add_argument("--class-map-json", type=Path, default=None)
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
    _add_cap_args(p3)

    p123 = subparsers.add_parser("run123", help="Run stage1, stage2, then stage3")
    p123.add_argument("--input-dir", type=Path, required=True)
    p123.add_argument("--output-dir", type=Path, required=True)
    p123.add_argument("--prefixes", type=str, default="")
    p123.add_argument("--stats", choices=("size", "counts"), default="size")
    p123.add_argument("--suffixes", type=str, default="")
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
    _add_cap_args(p123)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "stage1":
        prefixes = _split_csv(args.prefixes)
        result = run_stage1(
            args.input_dir,
            args.output_dir,
            prefixes,
            stats=args.stats,
            suffixes=parse_suffixes(args.suffixes),
            cap_min=args.cap_min,
            cap_max=args.cap_max,
        )
        print(f"Stage1 complete: {result.output_dir}")
        print(f"Stage1 class map: {result.output_json}")
        print(f"Stage1 filtered class map: {result.filtered_output_json}")
        if result.stats_json_paths:
            print(json.dumps({k: str(v) for k, v in result.stats_json_paths.items()}, indent=2))
        if result.mtime_json_paths:
            print(json.dumps({k: str(v) for k, v in result.mtime_json_paths.items()}, indent=2))
        if result.missed_txt_paths:
            print(json.dumps({k: str(v) for k, v in result.missed_txt_paths.items()}, indent=2))
        if result.small_txt_paths:
            print(json.dumps({k: str(v) for k, v in result.small_txt_paths.items()}, indent=2))
        return

    if args.command == "stage2":
        pipeline_root = infer_pipeline_root(args.output_dir)
        class_map_path = resolve_class_map_json(
            args.class_map_json, pipeline_root, cap_min=args.cap_min, cap_max=args.cap_max
        )
        class_map = load_class_map(explicit=class_map_path)
        class_to_suffix = parse_class_suffix_pairs(args.class_suffix)
        classes_to_process = parse_stage2_classes_to_process(args.classes_to_process)
        stage2_out_dir = resolve_stage2_output_dir(args.output_dir)
        outputs = run_stage2(
            input_dir=args.input_dir,
            output_dir=stage2_out_dir,
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
        stage1_result = run_stage1(
            args.input_dir,
            args.output_dir,
            prefixes,
            stats=args.stats,
            suffixes=parse_suffixes(args.suffixes),
            cap_min=args.cap_min,
            cap_max=args.cap_max,
        )
        class_to_suffix = parse_class_suffix_pairs(args.class_suffix)
        classes_to_process = parse_stage2_classes_to_process(args.classes_to_process)
        stage2_out_dir = resolve_stage2_output_dir(args.output_dir)
        outputs = run_stage2(
            input_dir=args.input_dir,
            output_dir=stage2_out_dir,
            class_map=stage1_result.filtered_classes_to_patients,
            default_suffix=args.suffix,
            class_to_suffix=class_to_suffix,
            batch_size=args.batch_size,
            task=args.task,
            classes_to_process=classes_to_process,
            save_inh=args.save_inh,
            save_denovo=args.save_denovo,
            use_ext_denovo=args.use_ext_denovo,
        )
        print(f"Stage1 file: {stage1_result.filtered_output_json}")
        print(json.dumps([o.__dict__ for o in outputs], indent=2, default=str))
        return

    if args.command == "stage3":
        pipeline_root = args.output_dir
        class_map_path = resolve_class_map_json(
            args.class_map_json, pipeline_root, cap_min=args.cap_min, cap_max=args.cap_max
        )
        class_map = load_class_map(explicit=class_map_path)
        class_to_cap = parse_class_cap_pairs(args.class_cap)
        requested = parse_classes_to_process(args.classes_to_process)
        class_map_run, stage3_dir = resolve_stage3_output(
            args.output_dir, class_map, requested, task=args.task
        )
        stage2_dir = resolve_stage2_input_dir(args.stage2_dir, pipeline_root)
        result = run_stage3(
            stage2_dir=stage2_dir,
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
        stage1_result = run_stage1(
            args.input_dir,
            args.output_dir,
            prefixes,
            stats=args.stats,
            suffixes=parse_suffixes(args.suffixes),
            cap_min=args.cap_min,
            cap_max=args.cap_max,
        )
        class_to_suffix = parse_class_suffix_pairs(args.class_suffix)
        stage2_dir = resolve_stage2_output_dir(args.output_dir)
        filtered_map = stage1_result.filtered_classes_to_patients

        stage2_outputs = run_stage2(
            input_dir=args.input_dir,
            output_dir=stage2_dir,
            class_map=filtered_map,
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
            args.output_dir, filtered_map, requested, task=args.task
        )
        stage3_result = run_stage3(
            stage2_dir=stage2_dir,
            output_dir=stage3_dir,
            class_map=class_map_run,
            default_cap=args.denovo_cap,
            class_to_cap=class_to_cap,
            task=args.task,
        )
        print(f"Stage1 file: {stage1_result.filtered_output_json}")
        print(json.dumps([o.__dict__ for o in stage2_outputs], indent=2, default=str))
        print(json.dumps(stage3_result.__dict__, indent=2, default=str))
        return


if __name__ == "__main__":
    main()
