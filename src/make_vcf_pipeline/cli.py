from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List

from .stage1 import (
    DEFAULT_CAP_MAX,
    DEFAULT_CAP_MIN,
    Stage2RunConfig,
    infer_pipeline_root,
    load_class_map,
    load_stage1_parameters_from_result,
    print_stage1_summary,
    resolve_class_map_json,
    resolve_stage2_input_dir,
    resolve_stage2_output_dir,
    resolve_stage2_run_config,
    run_stage1,
    stage1_parameters_path,
)
from .stage2 import (
    parse_class_suffix_pairs,
    print_stage2_summary,
    run_stage2,
)
from .stage3 import (
    filter_caps_from_args,
    parse_class_cap_pairs,
    parse_classes_to_process,
    print_stage3_summary,
    resolve_stage3_classes,
    resolve_stage3_output,
    run_stage3,
    DEFAULT_FILTER_AB_HET,
    DEFAULT_FILTER_AB_HOM0,
    DEFAULT_FILTER_AB_HOM1,
    DEFAULT_FILTER_DP,
    DEFAULT_FILTER_QT,
)


def _split_csv(value: str) -> List[str]:
    if not value.strip():
        return []
    return [x.strip() for x in value.split(",") if x.strip()]


def _add_stage2_hist_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--stage2-nohist",
        action="store_true",
        help="Stage2: skip qt/dp/ab histogram JSON files (computed by default per class_<name>/).",
    )


def _add_stage3_snv_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--snv",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Stage3: drop non-SNV variant keys before output (default: enabled). Use --no-snv to disable.",
    )


def _add_stage3_filter_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--filter",
        action="store_true",
        help="Stage3: after per-patient cap, apply is_good QC on mother, father, and child.",
    )
    parser.add_argument(
        "--filter-dp",
        type=int,
        default=DEFAULT_FILTER_DP,
        metavar="DP_CAP",
        help=f"Stage3 with --filter: minimum depth (keep if dp >= dp_cap; default {DEFAULT_FILTER_DP}).",
    )
    parser.add_argument(
        "--filter-qt",
        type=int,
        default=DEFAULT_FILTER_QT,
        metavar="QT_CAP",
        help=f"Stage3 with --filter: minimum genotype quality (keep if qt >= qt_cap; default {DEFAULT_FILTER_QT}).",
    )
    parser.add_argument(
        "--filter-abHom0",
        dest="filter_ab_hom0",
        type=float,
        default=DEFAULT_FILTER_AB_HOM0,
        metavar="AB_CAP",
        help=f"Stage3 with --filter: 0/0 alt fraction ab=ada/(adr+ada) must be < this (default {DEFAULT_FILTER_AB_HOM0}).",
    )
    parser.add_argument(
        "--filter-abHom1",
        dest="filter_ab_hom1",
        type=float,
        default=DEFAULT_FILTER_AB_HOM1,
        metavar="AB_CAP",
        help=f"Stage3 with --filter: 1/1 must have ab > 1 - ab_capHom1 (default {DEFAULT_FILTER_AB_HOM1}).",
    )
    parser.add_argument(
        "--filter-abHet",
        dest="filter_ab_het",
        type=float,
        default=DEFAULT_FILTER_AB_HET,
        metavar="AB_CAP",
        help=f"Stage3 with --filter: 0/1 must have ab_capHet < ab < 1 - ab_capHet (default {DEFAULT_FILTER_AB_HET}).",
    )


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
        "--classes",
        type=str,
        default="",
        help="Comma-separated class name prefixes for patient ID matching. Empty means uniform cohort.",
    )
    p1.add_argument(
        "--stats",
        choices=("size", "counts"),
        default="size",
        help="Per-patient file statistic when --suffix or --class-suffix is set (default: size).",
    )
    p1.add_argument(
        "--suffix",
        type=str,
        default=None,
        help="Default file suffix for size filtering/stats (e.g. vcf.gz). Overridden by --class-suffix.",
    )
    p1.add_argument(
        "--class-suffix",
        action="append",
        default=[],
        help="Per-class file suffix for size filtering/stats. Repeat format class=suffix.",
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
        "--classes",
        type=str,
        default="",
        help="Optional subset of classes to process (comma-separated). "
        "Default: all classes from stage1_parameters.json.",
    )
    p2.add_argument("--batch-size", type=int, default=1000)
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
    _add_stage2_hist_args(p2)

    p12 = subparsers.add_parser("run12", help="Run stage1 then stage2")
    p12.add_argument("--input-dir", type=Path, required=True)
    p12.add_argument("--output-dir", type=Path, required=True)
    p12.add_argument("--classes", type=str, default="")
    p12.add_argument("--stats", choices=("size", "counts"), default="size")
    p12.add_argument("--suffix", type=str, default=None)
    p12.add_argument("--class-suffix", action="append", default=[])
    p12.add_argument("--batch-size", type=int, default=1000)
    p12.add_argument(
        "--stage2-classes",
        type=str,
        default="",
        help="Stage2: optional subset of classes (comma-separated). Default: all from stage1 parameters.",
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
    _add_stage2_hist_args(p12)

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
        help="Pipeline output root; stage3 writes under stage3_vN/ or stage3_<classes>_vN/ (denovo), or stage3_inherited_vN/...",
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
        help="denovo: merge batch dVars. inherited: merge batch dVars_inh; output under stage3_inherited_vN/...",
    )
    _add_cap_args(p3)
    _add_stage3_snv_args(p3)
    _add_stage3_filter_args(p3)

    p123 = subparsers.add_parser("run123", help="Run stage1, stage2, then stage3")
    p123.add_argument("--input-dir", type=Path, required=True)
    p123.add_argument("--output-dir", type=Path, required=True)
    p123.add_argument("--classes", type=str, default="")
    p123.add_argument("--stats", choices=("size", "counts"), default="size")
    p123.add_argument("--suffix", type=str, default=None)
    p123.add_argument("--class-suffix", action="append", default=[])
    p123.add_argument("--batch-size", type=int, default=1000)
    p123.add_argument(
        "--stage2-classes",
        type=str,
        default="",
        help="Stage2: optional subset of classes (comma-separated). Default: all from stage1 parameters.",
    )
    p123.add_argument("--denovo-cap", type=int, default=None)
    p123.add_argument("--class-cap", action="append", default=[])
    p123.add_argument(
        "--classes-to-process",
        action="append",
        default=[],
        help="Optional: only these classes in stage3 (comma-separated or repeat). Omit for all classes.",
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
    _add_stage2_hist_args(p123)
    _add_stage3_snv_args(p123)
    _add_stage3_filter_args(p123)
    return parser


def _stage3_filter_caps(args: argparse.Namespace):
    return filter_caps_from_args(
        filter_enabled=args.filter,
        filter_dp=args.filter_dp,
        filter_qt=args.filter_qt,
        filter_ab_hom0=args.filter_ab_hom0,
        filter_ab_hom1=args.filter_ab_hom1,
        filter_ab_het=args.filter_ab_het,
    )


def _stage2_classes_override(args: argparse.Namespace, *, combined: bool) -> List[str] | None:
    csv = args.stage2_classes if combined else args.classes
    override = _split_csv(csv)
    return override if override else None


def _resolve_stage2_run(
    pipeline_root: Path,
    class_map: dict[str, list[str]],
    args: argparse.Namespace,
    *,
    combined: bool = False,
    stage1_params=None,
    stage1_params_path: Path | None = None,
) -> tuple[Stage2RunConfig, Path]:
    config = resolve_stage2_run_config(
        pipeline_root,
        class_map,
        _stage2_classes_override(args, combined=combined),
        stage1_params=stage1_params,
    )
    params_path = stage1_params_path or stage1_parameters_path(pipeline_root)
    return config, params_path


def _run_stage3_from_context(
    args: argparse.Namespace,
    pipeline_root: Path,
    class_map: dict[str, list[str]],
    *,
    stage2_classes: List[str] | None = None,
    stage2_dir: Path | None = None,
):
    class_names = resolve_stage3_classes(
        class_map,
        classes_override=parse_classes_to_process(args.classes_to_process) or None,
        stage2_classes=stage2_classes,
        pipeline_root=pipeline_root,
    )
    class_map_run, stage3_dir = resolve_stage3_output(
        pipeline_root, class_map, class_names, task=args.task
    )
    result = run_stage3(
        stage2_dir=stage2_dir or resolve_stage2_input_dir(None, pipeline_root),
        output_dir=stage3_dir,
        class_map=class_map_run,
        default_cap=args.denovo_cap,
        class_to_cap=parse_class_cap_pairs(args.class_cap),
        filter_caps=_stage3_filter_caps(args),
        snv_only=args.snv,
        task=args.task,
    )
    return result


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "stage1":
        classes = _split_csv(args.classes)
        class_to_suffix = parse_class_suffix_pairs(args.class_suffix)
        result = run_stage1(
            args.input_dir,
            args.output_dir,
            classes,
            stats=args.stats,
            default_suffix=args.suffix,
            class_to_suffix=class_to_suffix,
            cap_min=args.cap_min,
            cap_max=args.cap_max,
        )
        print(f"Stage1 complete: {result.output_dir}")
        print(f"Stage1 class map: {result.output_json}")
        print(f"Stage1 filtered class map: {result.filtered_output_json}")
        print(f"Stage1 parameters: {result.parameters_json}")
        if result.stats_json_paths:
            print(json.dumps({k: str(v) for k, v in result.stats_json_paths.items()}, indent=2))
        if result.mtime_json_paths:
            print(json.dumps({k: str(v) for k, v in result.mtime_json_paths.items()}, indent=2))
        if result.missed_txt_paths:
            print(json.dumps({k: str(v) for k, v in result.missed_txt_paths.items()}, indent=2))
        if result.small_txt_paths:
            print(json.dumps({k: str(v) for k, v in result.small_txt_paths.items()}, indent=2))
        print_stage1_summary(result)
        return

    if args.command == "stage2":
        pipeline_root = infer_pipeline_root(args.output_dir)
        class_map_path = resolve_class_map_json(
            args.class_map_json, pipeline_root, cap_min=args.cap_min, cap_max=args.cap_max
        )
        class_map = load_class_map(explicit=class_map_path)
        stage2_config, stage1_params_path = _resolve_stage2_run(
            pipeline_root, class_map, args, combined=False
        )
        stage2_out_dir = resolve_stage2_output_dir(args.output_dir)
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
            stage1_parameters_path=stage1_params_path,
        )
        print(json.dumps([o.__dict__ for o in stage2_result.outputs], indent=2, default=str))
        print(f"Stage2 parameters: {stage2_result.parameters_json}")
        print_stage2_summary(stage2_result, task=args.task)
        return

    if args.command == "run12":
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
        stage1_params = load_stage1_parameters_from_result(stage1_result)
        stage2_config, stage1_params_path = _resolve_stage2_run(
            args.output_dir,
            stage1_result.filtered_classes_to_patients,
            args,
            combined=True,
            stage1_params=stage1_params,
            stage1_params_path=stage1_result.parameters_json,
        )
        stage2_out_dir = resolve_stage2_output_dir(args.output_dir)
        stage2_result = run_stage2(
            input_dir=args.input_dir,
            output_dir=stage2_out_dir,
            class_map=stage1_result.filtered_classes_to_patients,
            suffix_per_class=stage2_config.suffix_per_class,
            classes=stage2_config.classes,
            batch_size=args.batch_size,
            task=args.task,
            save_inh=args.save_inh,
            save_denovo=args.save_denovo,
            use_ext_denovo=args.use_ext_denovo,
            write_hist=not args.stage2_nohist,
            stage1_parameters_path=stage1_params_path,
        )
        print(f"Stage1 file: {stage1_result.filtered_output_json}")
        print_stage1_summary(stage1_result)
        print(json.dumps([o.__dict__ for o in stage2_result.outputs], indent=2, default=str))
        print(f"Stage2 parameters: {stage2_result.parameters_json}")
        print_stage2_summary(stage2_result, task=args.task)
        return

    if args.command == "stage3":
        pipeline_root = args.output_dir
        class_map_path = resolve_class_map_json(
            args.class_map_json, pipeline_root, cap_min=args.cap_min, cap_max=args.cap_max
        )
        class_map = load_class_map(explicit=class_map_path)
        result = _run_stage3_from_context(
            args, pipeline_root, class_map, stage2_dir=resolve_stage2_input_dir(args.stage2_dir, pipeline_root)
        )
        print(json.dumps(result.__dict__, indent=2, default=str))
        print(f"Stage3 parameters: {result.parameters_json}")
        print_stage3_summary(result, task=args.task)
        return

    if args.command == "run123":
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
        stage2_dir = resolve_stage2_output_dir(args.output_dir)
        filtered_map = stage1_result.filtered_classes_to_patients

        stage1_params = load_stage1_parameters_from_result(stage1_result)
        stage2_config, stage1_params_path = _resolve_stage2_run(
            args.output_dir,
            filtered_map,
            args,
            combined=True,
            stage1_params=stage1_params,
            stage1_params_path=stage1_result.parameters_json,
        )
        stage2_result = run_stage2(
            input_dir=args.input_dir,
            output_dir=stage2_dir,
            class_map=filtered_map,
            suffix_per_class=stage2_config.suffix_per_class,
            classes=stage2_config.classes,
            batch_size=args.batch_size,
            task=args.task,
            save_inh=args.save_inh,
            save_denovo=args.save_denovo,
            use_ext_denovo=args.use_ext_denovo,
            write_hist=not args.stage2_nohist,
            stage1_parameters_path=stage1_params_path,
        )
        stage3_result = _run_stage3_from_context(
            args,
            args.output_dir,
            filtered_map,
            stage2_classes=stage2_config.classes,
            stage2_dir=stage2_dir,
        )
        print(f"Stage1 file: {stage1_result.filtered_output_json}")
        print_stage1_summary(stage1_result)
        print(json.dumps([o.__dict__ for o in stage2_result.outputs], indent=2, default=str))
        print(f"Stage2 parameters: {stage2_result.parameters_json}")
        print_stage2_summary(stage2_result, task=args.task)
        print(json.dumps(stage3_result.__dict__, indent=2, default=str))
        print(f"Stage3 parameters: {stage3_result.parameters_json}")
        print_stage3_summary(stage3_result, task=args.task)
        return


if __name__ == "__main__":
    main()
