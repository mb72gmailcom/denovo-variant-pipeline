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

from make_vcf_pipeline.stage1 import (
    DEFAULT_CAP_MAX,
    DEFAULT_CAP_MIN,
    load_class_map,
    load_stage1_parameters_from_result,
    print_stage1_summary,
    resolve_stage2_input_dir,
    resolve_stage2_output_dir,
    resolve_stage2_run_config,
    run_stage1,
    stage1_parameters_path,
)
from make_vcf_pipeline.stage2 import (
    parse_class_suffix_pairs,
    print_stage2_summary,
    run_stage2,
)
from make_vcf_pipeline.stage3 import (
    DEFAULT_FILTER_AB_HET,
    DEFAULT_FILTER_AB_HOM0,
    DEFAULT_FILTER_AB_HOM1,
    DEFAULT_FILTER_DP,
    DEFAULT_FILTER_QT,
    filter_caps_from_args,
    parse_class_cap_pairs,
    parse_classes_to_process,
    print_stage3_summary,
    resolve_stage3_classes,
    resolve_stage3_output,
    run_stage3,
)


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
        "--classes",
        type=str,
        default="",
        help="Stage1: comma-separated class name prefixes (empty means one class).",
    )
    parser.add_argument(
        "--stats",
        choices=("size", "counts"),
        default="size",
        help="Stage1: per-patient file statistic when --suffix or --class-suffix is set (default: size).",
    )
    parser.add_argument(
        "--cap-min",
        type=int,
        default=DEFAULT_CAP_MIN,
        help="Stage1 size filter lower bound (bytes); stage2/3 read filtered class map with these caps.",
    )
    parser.add_argument(
        "--cap-max",
        type=int,
        default=DEFAULT_CAP_MAX,
        help="Stage1 size filter upper bound (bytes); stage2/3 read filtered class map with these caps.",
    )
    parser.add_argument(
        "--suffix",
        type=str,
        default=None,
        help="Stage1: default file suffix (e.g. vcf.gz). Overridden by --class-suffix per class.",
    )
    parser.add_argument(
        "--class-suffix",
        action="append",
        default=[],
        help="Stage1: per-class file suffix. Repeat format class=suffix.",
    )

    # Stage 2 parameters
    parser.add_argument(
        "--stage2-classes",
        type=str,
        default="",
        help="Stage2: optional subset of classes (comma-separated). Default: all from stage1 parameters.",
    )
    parser.add_argument("--batch-size", type=int, default=1000, help="Stage2: patient batch size")
    parser.add_argument(
        "--task",
        choices=("denovo", "inherited"),
        default="denovo",
        help="Stage2 sidecar lists + stage3 source (dVars vs dVars_inh) and output dir (stage3 vs stage3_inherited).",
    )
    parser.add_argument(
        "--save-inh",
        action="store_true",
        help="Stage2: with --task denovo, also collect/save dVars_inh.",
    )
    parser.add_argument(
        "--save-denovo",
        "--save_denovo",
        action="store_true",
        help="Stage2: with --task inherited, also collect/save dVars.",
    )
    parser.add_argument(
        "--use-ext-denovo",
        action="store_true",
        help="Stage2: classify denovo with if_denovo_ext() instead of if_denovo().",
    )
    parser.add_argument(
        "--stage2-nohist",
        action="store_true",
        help="Stage2: skip qt/dp/ab histogram JSON files (computed by default per class).",
    )

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
    parser.add_argument(
        "--classes-to-process",
        action="append",
        default=[],
        help="Stage3: only these classes (comma-separated or repeat). Omit to process all classes.",
    )
    parser.add_argument(
        "--snv",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Stage3: drop non-SNV variant keys (default: enabled). Use --no-snv to disable.",
    )
    parser.add_argument(
        "--filter",
        action="store_true",
        help="Stage3: after per-patient cap, apply is_good QC on mother, father, and child.",
    )
    parser.add_argument("--filter-dp", type=int, default=DEFAULT_FILTER_DP, metavar="DP_CAP")
    parser.add_argument("--filter-qt", type=int, default=DEFAULT_FILTER_QT, metavar="QT_CAP")
    parser.add_argument(
        "--filter-abHom0", dest="filter_ab_hom0", type=float, default=DEFAULT_FILTER_AB_HOM0, metavar="AB_CAP"
    )
    parser.add_argument(
        "--filter-abHom1", dest="filter_ab_hom1", type=float, default=DEFAULT_FILTER_AB_HOM1, metavar="AB_CAP"
    )
    parser.add_argument(
        "--filter-abHet", dest="filter_ab_het", type=float, default=DEFAULT_FILTER_AB_HET, metavar="AB_CAP"
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if not (args.run_stage1 or args.run_stage2 or args.run_stage3):
        raise ValueError("Select at least one stage flag: --run-stage1/--run-stage2/--run-stage3")

    stage2_out_dir = resolve_stage2_output_dir(args.output_dir)

    class_map = None
    stage1_result = None
    stage2_classes: List[str] | None = None

    if args.run_stage1:
        classes = _split_csv(args.classes)
        class_to_suffix = parse_class_suffix_pairs(args.class_suffix)
        stage1_result = run_stage1(
            args.input_dir,
            args.output_dir,
            classes,
            stats=args.stats,
            default_suffix=args.suffix,
            class_to_suffix=class_to_suffix,
            cap_min=args.cap_min,
            cap_max=args.cap_max,
        )
        class_map = stage1_result.filtered_classes_to_patients
        print(f"[stage1] done -> {stage1_result.output_dir}")
        print(f"[stage1] filtered class map -> {stage1_result.filtered_output_json}")
        print(f"[stage1] parameters -> {stage1_result.parameters_json}")
        if stage1_result.stats_json_paths:
            print(f"[stage1] stats -> {len(stage1_result.stats_json_paths)} files")
        if stage1_result.mtime_json_paths:
            print(f"[stage1] mtime -> {len(stage1_result.mtime_json_paths)} files")
        if stage1_result.missed_txt_paths:
            print(f"[stage1] missed -> {len(stage1_result.missed_txt_paths)} files")
        if stage1_result.small_txt_paths:
            print(f"[stage1] small -> {len(stage1_result.small_txt_paths)} files")
        print_stage1_summary(stage1_result)

    if args.run_stage2:
        if class_map is None:
            class_map = load_class_map(pipeline_root=args.output_dir, cap_min=args.cap_min, cap_max=args.cap_max)

        stage1_params = (
            load_stage1_parameters_from_result(stage1_result) if stage1_result is not None else None
        )
        classes_override = _split_csv(args.stage2_classes) if args.stage2_classes.strip() else None
        stage2_config = resolve_stage2_run_config(
            args.output_dir,
            class_map,
            classes_override,
            stage1_params=stage1_params,
        )
        stage2_classes = stage2_config.classes
        stage2_result = run_stage2(
            input_dir=args.input_dir,
            output_dir=stage2_out_dir,
            class_map=class_map,
            suffix_per_class=stage2_config.suffix_per_class,
            classes=stage2_config.classes,
            batch_size=args.batch_size,
            task=args.task,
            save_inh=args.save_inh,
            save_denovo=args.save_denovo,
            use_ext_denovo=args.use_ext_denovo,
            write_hist=not args.stage2_nohist,
            stage1_parameters_path=(
                stage1_result.parameters_json if stage1_result is not None else stage1_parameters_path(args.output_dir)
            ),
        )
        print(f"[stage2] done -> {len(stage2_result.outputs)} batch outputs in {stage2_out_dir}")
        print(f"[stage2] parameters -> {stage2_result.parameters_json}")
        print_stage2_summary(stage2_result, task=args.task)

    if args.run_stage3:
        if class_map is None:
            class_map = load_class_map(pipeline_root=args.output_dir, cap_min=args.cap_min, cap_max=args.cap_max)

        class_names = resolve_stage3_classes(
            class_map,
            classes_override=parse_classes_to_process(args.classes_to_process) or None,
            stage2_classes=stage2_classes,
            pipeline_root=args.output_dir,
        )
        class_map_run, stage3_out_dir = resolve_stage3_output(
            args.output_dir, class_map, class_names, task=args.task
        )
        stage3_result = run_stage3(
            stage2_dir=resolve_stage2_input_dir(None, args.output_dir),
            output_dir=stage3_out_dir,
            class_map=class_map_run,
            default_cap=args.denovo_cap,
            class_to_cap=parse_class_cap_pairs(args.class_cap),
            filter_caps=filter_caps_from_args(
                filter_enabled=args.filter,
                filter_dp=args.filter_dp,
                filter_qt=args.filter_qt,
                filter_ab_hom0=args.filter_ab_hom0,
                filter_ab_hom1=args.filter_ab_hom1,
                filter_ab_het=args.filter_ab_het,
            ),
            snv_only=args.snv,
            task=args.task,
        )
        print(f"[stage3] done -> {stage3_result.output_dir}")
        print(f"[stage3] parameters -> {stage3_result.parameters_json}")
        print_stage3_summary(stage3_result, task=args.task)


if __name__ == "__main__":
    main()
