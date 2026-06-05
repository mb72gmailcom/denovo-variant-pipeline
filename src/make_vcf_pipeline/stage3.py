from __future__ import annotations

import json
import pickle
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List

from .stage2 import merge


VarsMap = Dict[str, List[str]]
FilterFn = Callable[[VarsMap, int], VarsMap]
ExtraFilter = FilterFn | tuple[str, FilterFn]

FILTER_STEP_LABELS: Dict[str, str] = {
    "original": "original variants",
    "patient_variant_cap": "after patient variant cap",
    "snv_only": "after SNV filter",
    "genotype_qc": "after genotype QC (depth/quality/AB)",
}

TRANSITION_PAIRS = frozenset({("A", "G"), ("G", "A"), ("T", "C"), ("C", "T")})
STAGE3_CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

DEFAULT_FILTER_DP = 20
DEFAULT_FILTER_QT = 90
DEFAULT_FILTER_AB_HOM0 = 0.05
DEFAULT_FILTER_AB_HOM1 = 0.05
DEFAULT_FILTER_AB_HET = 0.30


def parse_classes_to_process(entries: List[str]) -> List[str]:
    """Flatten repeatable/comma-separated class names; preserve first-seen order, dedupe."""
    out: List[str] = []
    seen: set[str] = set()
    for raw in entries:
        for chunk in raw.split(","):
            name = chunk.strip()
            if name and name not in seen:
                seen.add(name)
                out.append(name)
    return out


def resolve_stage3_stem(
    class_map: Dict[str, List[str]],
    classes_to_process: List[str],
    task: str = "denovo",
) -> str:
    """Base directory name before ``_vN`` (e.g. ``stage3``, ``stage3_SSC_ABC``)."""
    if task not in ("denovo", "inherited"):
        raise ValueError("task must be 'denovo' or 'inherited'")

    base = "stage3" if task == "denovo" else "stage3_inherited"

    if not classes_to_process:
        return base

    missing = [c for c in classes_to_process if c not in class_map]
    if missing:
        raise ValueError(
            f"Unknown --classes-to-process: {missing}. Known classes: {sorted(class_map.keys())}"
        )

    selected_keys = set(classes_to_process)
    if selected_keys == set(class_map.keys()):
        return base

    tag = "_".join(sorted(selected_keys))
    return f"{base}_{tag}"


def allocate_monotonic_versioned_dir(output_root: Path, stem: str) -> tuple[Path, int]:
    """
    Pick ``{stem}_v{N}`` where N is one greater than the highest existing version (no gap fill).
    """
    pattern = re.compile(rf"^{re.escape(stem)}_v(\d+)$")
    max_version = -1
    if output_root.is_dir():
        for child in output_root.iterdir():
            if not child.is_dir():
                continue
            match = pattern.match(child.name)
            if match:
                max_version = max(max_version, int(match.group(1)))

    version = max_version + 1
    return output_root / f"{stem}_v{version}", version


def parse_versioned_dir_name(dir_name: str) -> tuple[str, int]:
    stem, sep, version_text = dir_name.rpartition("_v")
    if sep and version_text.isdigit():
        return stem, int(version_text)
    return dir_name, 0


def resolve_stage3_output(
    output_root: Path,
    class_map: Dict[str, List[str]],
    classes_to_process: List[str],
    task: str = "denovo",
) -> tuple[Dict[str, List[str]], Path]:
    """
    Choose classes and a versioned output directory under output_root.

    Stems: ``stage3``, ``stage3_inherited``, ``stage3_<classes>``, ``stage3_inherited_<classes>``.
    Each run writes to ``{stem}_vN`` with N = max(existing N) + 1 (monotonic; gaps are not reused).
    """
    if not classes_to_process:
        class_map_run = dict(class_map)
    else:
        missing = [c for c in classes_to_process if c not in class_map]
        if missing:
            raise ValueError(
                f"Unknown --classes-to-process: {missing}. Known classes: {sorted(class_map.keys())}"
            )
        class_map_run = {k: class_map[k] for k in classes_to_process}

    stem = resolve_stage3_stem(class_map, classes_to_process, task=task)
    output_dir, _version = allocate_monotonic_versioned_dir(output_root, stem)
    return class_map_run, output_dir


def parse_class_cap_pairs(pairs: List[str]) -> Dict[str, int]:
    mapping: Dict[str, int] = {}
    for raw in pairs:
        chunks = [chunk.strip() for chunk in raw.split(",") if chunk.strip()]
        for pair in chunks:
            if "=" not in pair:
                raise ValueError(f"Invalid --class-cap '{pair}'. Expected format class=cap")
            class_name, cap = pair.split("=", 1)
            mapping[class_name.strip()] = int(cap.strip())
    return mapping


def _load_vars_file(path: Path) -> VarsMap:
    if path.suffix == ".json":
        return json.loads(path.read_text(encoding="utf-8"))
    with path.open("rb") as f:
        return pickle.load(f)


def load_combined_dvars_for_class(stage2_dir: Path, class_name: str) -> VarsMap:
    class_dir = stage2_dir / f"class_{class_name}"
    combined: VarsMap = {}
    if not class_dir.is_dir():
        return combined

    batch_files = sorted(list(class_dir.glob("batch_*_dVars.pkl")) + list(class_dir.glob("batch_*_dVars.json")))
    for file_path in batch_files:
        dvars = _load_vars_file(file_path)
        merge(combined, dvars)
    return combined


def load_combined_dvars_inh_for_class(stage2_dir: Path, class_name: str) -> VarsMap:
    class_dir = stage2_dir / f"class_{class_name}"
    combined: VarsMap = {}
    if not class_dir.is_dir():
        return combined

    batch_files = sorted(
        list(class_dir.glob("batch_*_dVars_inh.pkl")) + list(class_dir.glob("batch_*_dVars_inh.json"))
    )
    for file_path in batch_files:
        dvars_inh = _load_vars_file(file_path)
        merge(combined, dvars_inh)
    return combined


def denovo_counts_by_patient(dvars: VarsMap) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for values in dvars.values():
        for v in values:
            patient_id = v.split("-")[0]
            counts[patient_id] = counts.get(patient_id, 0) + 1
    return counts


@dataclass(frozen=True)
class FilterCaps:
    dp_cap: int
    qt_cap: int
    ab_cap_hom0: float
    ab_cap_hom1: float
    ab_cap_het: float


def filter_caps_from_args(
    *,
    filter_enabled: bool,
    filter_dp: int | None,
    filter_qt: int | None,
    filter_ab_hom0: float | None,
    filter_ab_hom1: float | None,
    filter_ab_het: float | None,
) -> FilterCaps | None:
    if not filter_enabled:
        return None

    return FilterCaps(
        dp_cap=filter_dp if filter_dp is not None else DEFAULT_FILTER_DP,
        qt_cap=filter_qt if filter_qt is not None else DEFAULT_FILTER_QT,
        ab_cap_hom0=filter_ab_hom0 if filter_ab_hom0 is not None else DEFAULT_FILTER_AB_HOM0,
        ab_cap_hom1=filter_ab_hom1 if filter_ab_hom1 is not None else DEFAULT_FILTER_AB_HOM1,
        ab_cap_het=filter_ab_het if filter_ab_het is not None else DEFAULT_FILTER_AB_HET,
    )


def is_good(
    gt: str,
    gex: str,
    dp_cap: int,
    qt_cap: int,
    ab_cap_hom0: float,
    ab_cap_hom1: float,
    ab_cap_het: float,
) -> bool:
    parts = gex.split(":")
    if len(parts) != 4:
        return False
    try:
        qt, dp, adr, ada = (int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]))
    except ValueError:
        return False

    denom = adr + ada
    if denom <= 0:
        return False
    ab = ada / denom

    if qt < qt_cap or dp < dp_cap:
        return False

    if gt == "0/0" and ab < ab_cap_hom0:
        return True
    if gt == "1/1" and ab > 1 - ab_cap_hom1:
        return True
    if gt == "0/1" and ab > ab_cap_het and ab < 1 - ab_cap_het:
        return True
    return False


def _compact_gtex_record(gtex: str) -> str:
    parts = gtex.split("-", 4)
    if len(parts) >= 2:
        return f"{parts[0]}-{parts[1]}"
    return parts[0] if parts else gtex


def passes_triplet_filter(gtex: str, caps: FilterCaps) -> bool:
    parts = gtex.split("-", 4)
    if len(parts) != 5:
        return False
    _, gt, mex, fex, cex = parts
    gt_parts = gt.split(":", 2)
    if len(gt_parts) != 3:
        return False
    mt, ft, ch = gt_parts
    return (
        is_good(mt, mex, caps.dp_cap, caps.qt_cap, caps.ab_cap_hom0, caps.ab_cap_hom1, caps.ab_cap_het)
        and is_good(ft, fex, caps.dp_cap, caps.qt_cap, caps.ab_cap_hom0, caps.ab_cap_hom1, caps.ab_cap_het)
        and is_good(ch, cex, caps.dp_cap, caps.qt_cap, caps.ab_cap_hom0, caps.ab_cap_hom1, caps.ab_cap_het)
    )


def filter_by_patient_denovo_cap(
    dvars: VarsMap,
    denovo_cap: int,
    *,
    keep_full_gtex: bool = False,
) -> VarsMap:
    counts = denovo_counts_by_patient(dvars)
    filtered: VarsMap = {}
    for key, values in dvars.items():
        for v in values:
            vv = v.split("-", 4)
            patient_id = vv[0]
            if counts.get(patient_id, 0) < denovo_cap:
                vs = v if keep_full_gtex else _compact_gtex_record(v)
                filtered.setdefault(key, []).append(vs)
    return filtered


def filter_by_is_good(dvars: VarsMap, caps: FilterCaps) -> VarsMap:
    filtered: VarsMap = {}
    for key, values in dvars.items():
        for v in values:
            if passes_triplet_filter(v, caps):
                filtered.setdefault(key, []).append(_compact_gtex_record(v))
    return filtered


def parse_snv_key(key: str) -> tuple[str, str, str] | None:
    parts = key.split("_")
    if len(parts) < 4:
        return None
    chrom, ref, alt = parts[0], parts[2], parts[3]
    if chrom not in STAGE3_CHROMOSOMES:
        return None
    if len(ref) != 1 or len(alt) != 1:
        return None
    return chrom, ref.upper(), alt.upper()


def count_snv_variant_keys(dvars: VarsMap) -> int:
    """Count variant keys with single-base ref and alt (same rule as VCF / Ti-Tv)."""
    return sum(1 for key in dvars if parse_snv_key(key) is not None)


def count_variant_keys(dvars: VarsMap) -> int:
    return len(dvars)


def filter_snv_variants(dvars: VarsMap) -> VarsMap:
    return {key: values for key, values in dvars.items() if parse_snv_key(key) is not None}


@dataclass(frozen=True)
class FilterStepCount:
    step: str
    variants_remaining: int

    @property
    def label(self) -> str:
        return FILTER_STEP_LABELS.get(self.step, self.step)


@dataclass(frozen=True)
class FilterPipelineResult:
    dvars: VarsMap
    steps: List[FilterStepCount]


def _normalize_extra_filter(item: ExtraFilter, index: int) -> tuple[str, FilterFn]:
    if isinstance(item, tuple):
        name, fn = item
        if not name.strip():
            raise ValueError("Extra filter name must be non-empty")
        return name.strip(), fn
    return f"extra_filter_{index}", item


def apply_filters(
    dvars: VarsMap,
    denovo_cap: int,
    *,
    snv_only: bool = True,
    filter_caps: FilterCaps | None = None,
    extra_filters: Iterable[ExtraFilter] | None = None,
) -> FilterPipelineResult:
    steps: List[FilterStepCount] = []
    steps.append(FilterStepCount("original", count_variant_keys(dvars)))

    # 1) per-patient variant cap (stage1 file-size selection is upstream)
    current = filter_by_patient_denovo_cap(
        dvars, denovo_cap, keep_full_gtex=filter_caps is not None
    )
    steps.append(FilterStepCount("patient_variant_cap", count_variant_keys(current)))

    # 2) SNV-only variant keys
    if snv_only:
        current = filter_snv_variants(current)
        steps.append(FilterStepCount("snv_only", count_variant_keys(current)))

    # 3) genotype quality (is_good), optional
    if filter_caps is not None:
        current = filter_by_is_good(current, filter_caps)
        steps.append(FilterStepCount("genotype_qc", count_variant_keys(current)))

    if extra_filters:
        for index, item in enumerate(extra_filters):
            step_name, fn = _normalize_extra_filter(item, index)
            current = fn(current, denovo_cap)
            steps.append(FilterStepCount(step_name, count_variant_keys(current)))

    return FilterPipelineResult(dvars=current, steps=steps)


def split_keys_by_chromosome(dvars: VarsMap) -> Dict[str, List[str]]:
    dchr: Dict[str, List[str]] = {chrom: [] for chrom in STAGE3_CHROMOSOMES}

    for key in sorted(dvars.keys()):
        row = key.split("_")
        if len(row) < 2:
            continue
        chrom = row[0]
        if chrom in dchr:
            dchr[chrom].append(key)
    return dchr


def sort_chromosome_keys(dchr: Dict[str, List[str]]) -> Dict[str, List[str]]:
    sorted_map: Dict[str, List[str]] = {}
    for chrom, keys in dchr.items():
        positions = [int(k.split("_")[1]) for k in keys]
        sorted_keys = [x for _, x in sorted(zip(positions, keys))]
        sorted_map[chrom] = sorted_keys
    return sorted_map


def write_vcf(dvars: VarsMap, output_dir: Path, dchr_sorted: Dict[str, List[str]]) -> None:
    """Write headerless per-chromosome VCF from filtered variant keys (no #CHROM line)."""
    for chrom, keys in dchr_sorted.items():
        chrom_dir = output_dir / chrom
        chrom_dir.mkdir(parents=True, exist_ok=True)
        out_file = chrom_dir / "variants.vcf"

        prev_pos = None
        prev_row = None
        with out_file.open("w", encoding="utf-8") as f:
            for key in keys:
                parts = key.split("_", 3)
                if len(parts) < 4:
                    continue
                row_chrom, pos, ref, alt = parts

                if pos != prev_pos:
                    if prev_row is not None:
                        f.write("\t".join(prev_row) + "\n")
                    prev_row = [row_chrom, pos, "rsXXX", ref, alt]
                    prev_pos = pos
                else:
                    prev_row[4] += "," + alt

            if prev_row is not None:
                f.write("\t".join(prev_row) + "\n")


def is_transition(ref: str, alt: str) -> bool:
    return (ref.upper(), alt.upper()) in TRANSITION_PAIRS


def titv_counts_for_chromosome(keys: Iterable[str]) -> Dict[str, int | float | None]:
    transitions = 0
    transversions = 0
    for key in keys:
        parsed = parse_snv_key(key)
        if parsed is None:
            continue
        _, ref, alt = parsed
        if ref == alt:
            continue
        if is_transition(ref, alt):
            transitions += 1
        else:
            transversions += 1
    snv_sites = transitions + transversions
    titv_ratio = transitions / transversions if transversions else None
    return {
        "transitions": transitions,
        "transversions": transversions,
        "snv_sites": snv_sites,
        "titv_ratio": titv_ratio,
    }


def compute_titv_summary(dchr_sorted: Dict[str, List[str]]) -> Dict[str, object]:
    by_chromosome: Dict[str, Dict[str, int | float | None]] = {}
    total_transitions = 0
    total_transversions = 0

    for chrom in STAGE3_CHROMOSOMES:
        counts = titv_counts_for_chromosome(dchr_sorted.get(chrom, []))
        by_chromosome[chrom] = counts
        total_transitions += int(counts["transitions"])
        total_transversions += int(counts["transversions"])

    overall_snv_sites = total_transitions + total_transversions
    overall_titv = total_transitions / total_transversions if total_transversions else None

    return {
        "by_chromosome": by_chromosome,
        "overall": {
            "transitions": total_transitions,
            "transversions": total_transversions,
            "snv_sites": overall_snv_sites,
            "titv_ratio": overall_titv,
        },
    }


@dataclass(frozen=True)
class Stage3ClassSummary:
    class_name: str
    variants_left: int
    transitions: int
    transversions: int
    titv_ratio: float | None
    output_dir: Path
    filter_steps: List[FilterStepCount]


@dataclass(frozen=True)
class Stage3Result:
    classes_processed: int
    variants_kept: int
    output_dir: Path
    class_summaries: List[Stage3ClassSummary]


def print_stage3_summary(result: Stage3Result, *, task: str) -> None:
    print(f"[stage3 summary] task={task} output={result.output_dir}")
    for row in result.class_summaries:
        titv_display = "n/a" if row.titv_ratio is None else f"{row.titv_ratio:.6f}"
        print(
            f"  {row.class_name}: variants_left={row.variants_left} "
            f"Ti/Tv={titv_display} (Ti={row.transitions} Tv={row.transversions}) "
            f"-> {row.output_dir}"
        )
        for step in row.filter_steps:
            print(f"    {step.label}: {step.variants_remaining}")


def _filter_steps_to_json(steps: List[FilterStepCount]) -> List[Dict[str, object]]:
    return [
        {
            "step": step.step,
            "label": step.label,
            "variants_remaining": step.variants_remaining,
        }
        for step in steps
    ]


def _write_class_stage3_outputs(
    class_dir: Path,
    class_name: str,
    task: str,
    dvars_filtered: VarsMap,
    dchr_sorted: Dict[str, List[str]],
    *,
    run_dir_name: str,
    patient_variant_cap: int,
    snv_only: bool,
    filter_caps: FilterCaps | None,
    filter_steps: List[FilterStepCount],
) -> None:
    class_dir.mkdir(parents=True, exist_ok=True)
    write_vcf(dvars_filtered, class_dir, dchr_sorted)

    titv_summary = compute_titv_summary(dchr_sorted)
    (class_dir / "stage3_titv.json").write_text(
        json.dumps(titv_summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    summary: Dict[str, object] = {
        "task": task,
        "class": class_name,
        "output_dir": class_dir.name,
        "run_dir": run_dir_name,
        "variants_kept": len(dvars_filtered),
        "chromosome_counts": {chrom: len(keys) for chrom, keys in dchr_sorted.items()},
        "patient_variant_cap": patient_variant_cap,
        "snv_only": snv_only,
        "filter_enabled": filter_caps is not None,
        "filter_pipeline": _filter_steps_to_json(filter_steps),
    }
    if filter_caps is not None:
        summary["filter_caps"] = {
            "dp_cap": filter_caps.dp_cap,
            "qt_cap": filter_caps.qt_cap,
            "ab_cap_hom0": filter_caps.ab_cap_hom0,
            "ab_cap_hom1": filter_caps.ab_cap_hom1,
            "ab_cap_het": filter_caps.ab_cap_het,
        }

    (class_dir / "stage3_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )


def run_stage3(
    stage2_dir: Path,
    output_dir: Path,
    class_map: Dict[str, List[str]],
    default_cap: int | None,
    class_to_cap: Dict[str, int],
    filter_caps: FilterCaps | None = None,
    snv_only: bool = True,
    extra_filters: Iterable[ExtraFilter] | None = None,
    task: str = "denovo",
) -> Stage3Result:
    if task not in ("denovo", "inherited"):
        raise ValueError("task must be 'denovo' or 'inherited'")

    output_dir.mkdir(parents=True, exist_ok=True)

    load_combo = load_combined_dvars_for_class if task == "denovo" else load_combined_dvars_inh_for_class

    classes_processed = 0
    class_summaries: List[Stage3ClassSummary] = []
    variants_kept_total = 0

    for class_name in class_map.keys():
        cap = class_to_cap.get(class_name, default_cap)
        if cap is None:
            raise ValueError(f"No denovo cap provided for class '{class_name}'")

        dvars_combo = load_combo(stage2_dir, class_name)
        filter_result = apply_filters(
            dvars_combo,
            cap,
            snv_only=snv_only,
            filter_caps=filter_caps,
            extra_filters=extra_filters,
        )
        dvars_filtered = filter_result.dvars
        dchr_class_sorted = sort_chromosome_keys(split_keys_by_chromosome(dvars_filtered))
        titv_class = compute_titv_summary(dchr_class_sorted)
        overall = titv_class["overall"]
        variants_left = len(dvars_filtered)

        class_dir = output_dir / f"class_{class_name}"
        _write_class_stage3_outputs(
            class_dir,
            class_name,
            task,
            dvars_filtered,
            dchr_class_sorted,
            run_dir_name=output_dir.name,
            patient_variant_cap=cap,
            snv_only=snv_only,
            filter_caps=filter_caps,
            filter_steps=filter_result.steps,
        )

        variants_kept_total += variants_left
        class_summaries.append(
            Stage3ClassSummary(
                class_name=class_name,
                variants_left=variants_left,
                transitions=int(overall["transitions"]),
                transversions=int(overall["transversions"]),
                titv_ratio=overall["titv_ratio"],
                output_dir=class_dir,
                filter_steps=filter_result.steps,
            )
        )
        classes_processed += 1

    return Stage3Result(
        classes_processed=classes_processed,
        variants_kept=variants_kept_total,
        output_dir=output_dir,
        class_summaries=class_summaries,
    )
