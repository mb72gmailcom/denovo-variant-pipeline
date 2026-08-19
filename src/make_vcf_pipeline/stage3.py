from __future__ import annotations

import json
import pickle
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Collection, Dict, Iterable, List

from .stage2 import build_patient_counts_json, merge, try_load_stage2_parameters


VarsMap = Dict[str, List[str]]
FilterFn = Callable[[VarsMap, int], VarsMap]
ExtraFilter = FilterFn | tuple[str, FilterFn]

FILTER_STEP_LABELS: Dict[str, str] = {
    "original": "original variants",
    "patient_variant_cap": "after patient variant cap",
    "chromosome": "after chromosome allowlist",
    "snv_only": "after SNV filter",
    "genotype_qc": "after genotype QC (depth/quality/AB)",
}

TRANSITION_PAIRS = frozenset({("A", "G"), ("G", "A"), ("T", "C"), ("C", "T")})
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
AUTOSOMAL_CHROMOSOMES = [f"chr{i}" for i in range(1, 23)]


def resolve_stage3_chromosomes(*, autosomal: bool) -> List[str]:
    return list(AUTOSOMAL_CHROMOSOMES if autosomal else CHROMOSOMES)

DEFAULT_FILTER_DP = 20
DEFAULT_FILTER_QT = 90
DEFAULT_FILTER_AB_HOM0 = 0.05
DEFAULT_FILTER_AB_HOM1 = 0.05
DEFAULT_FILTER_AB_HET = 0.30


def parse_stage3_classes(value: str) -> List[str]:
    """Parse comma-separated stage3 class names; preserve first-seen order, dedupe."""
    out: List[str] = []
    seen: set[str] = set()
    for chunk in value.split(","):
        name = chunk.strip()
        if name and name not in seen:
            seen.add(name)
            out.append(name)
    return out


def resolve_stage3_classes(
    class_map: Dict[str, List[str]],
    *,
    classes_override: List[str] | None = None,
    stage2_classes: List[str] | None = None,
    pipeline_root: Path | None = None,
) -> List[str]:
    """
    Choose class names for a stage3 run.

    Priority: explicit --stage3-classes > stage2 classes from same session >
    stage2_parameters.json. Does not default to all classes in the class map.
    """
    if classes_override:
        missing = [c for c in classes_override if c not in class_map]
        if missing:
            raise ValueError(
                f"Unknown --stage3-classes: {missing}. Known classes: {sorted(class_map.keys())}"
            )
        return classes_override

    if stage2_classes:
        selected = [c for c in stage2_classes if c in class_map]
        if not selected:
            raise ValueError("No stage2 classes found in the class map for stage3.")
        return selected

    if pipeline_root is not None:
        stage2_params = try_load_stage2_parameters(pipeline_root)
        if stage2_params and stage2_params.classes:
            selected = [c for c in stage2_params.classes if c in class_map]
            if selected:
                return selected

    raise ValueError(
        "No stage3 classes specified. Pass --stage3-classes or run stage2 first "
        "(stage3 reads classes from stage2/stage2_parameters.json)."
    )


def resolve_stage3_collect(
    *,
    collect_override: str | None = None,
    stage2_collect: str | None = None,
    pipeline_root: Path | None = None,
) -> str:
    """
    Choose collect bucket for stage3.

    Priority: explicit --collect > stage2 collect from same session >
    stage2_parameters.json.
    """
    if collect_override is not None:
        if collect_override not in ("denovo", "inherited"):
            raise ValueError("collect must be 'denovo' or 'inherited'")
        return collect_override

    if stage2_collect is not None:
        if stage2_collect not in ("denovo", "inherited"):
            raise ValueError("collect must be 'denovo' or 'inherited'")
        return stage2_collect

    if pipeline_root is not None:
        stage2_params = try_load_stage2_parameters(pipeline_root)
        if stage2_params is not None:
            return stage2_params.collect

    raise ValueError(
        "No stage3 collect specified. Pass --collect or run stage2 first "
        "(stage3 reads collect from stage2/stage2_parameters.json)."
    )


def resolve_stage3_stem(collect: str = "denovo") -> str:
    """Base directory name before ``_vN`` (``stage3`` or ``stage3_inherited``)."""
    if collect not in ("denovo", "inherited"):
        raise ValueError("collect must be 'denovo' or 'inherited'")
    return "stage3" if collect == "denovo" else "stage3_inherited"


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
    collect: str = "denovo",
) -> tuple[Dict[str, List[str]], Path]:
    """
    Choose classes and a versioned output directory under output_root.

    Stems: ``stage3`` or ``stage3_inherited`` (class subset does not affect the name).
    Each run writes to ``{stem}_vN`` with N = max(existing N) + 1 (monotonic; gaps are not reused).
    """
    if not classes_to_process:
        class_map_run = dict(class_map)
    else:
        missing = [c for c in classes_to_process if c not in class_map]
        if missing:
            raise ValueError(
                f"Unknown --stage3-classes: {missing}. Known classes: {sorted(class_map.keys())}"
            )
        class_map_run = {k: class_map[k] for k in classes_to_process}

    stem = resolve_stage3_stem(collect=collect)
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


def variant_patient_counts(dvars: VarsMap) -> Dict[str, int]:
    return {key: len(values) for key, values in sorted(dvars.items())}


def patients_per_variant_for_keys(dvars: VarsMap, keys: List[str]) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for key in keys:
        values = dvars.get(key)
        if not values:
            continue
        patients = sorted({v.split("-", 1)[0] for v in values if v})
        if patients:
            out[key] = patients
    return out


def write_class_level_stage3_stats(
    class_dir: Path,
    *,
    patient_ids: List[str],
    dvars_filtered: VarsMap,
) -> None:
    (class_dir / "patient_variant_counts.json").write_text(
        json.dumps(
            build_patient_counts_json(patient_ids, denovo_counts_by_patient(dvars_filtered)),
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    (class_dir / "variant_patient_counts.json").write_text(
        json.dumps(variant_patient_counts(dvars_filtered), indent=2, sort_keys=True),
        encoding="utf-8",
    )


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


def parse_snv_key(key: str, *, allowed_chromosomes: Collection[str]) -> tuple[str, str, str] | None:
    parts = key.split("_")
    if len(parts) < 4:
        return None
    chrom, ref, alt = parts[0], parts[2], parts[3]
    if chrom not in allowed_chromosomes:
        return None
    if len(ref) != 1 or len(alt) != 1:
        return None
    return chrom, ref.upper(), alt.upper()


def count_snv_variant_keys(dvars: VarsMap, *, allowed_chromosomes: Collection[str]) -> int:
    """Count variant keys with single-base ref and alt (same rule as VCF / Ti-Tv)."""
    return sum(1 for key in dvars if parse_snv_key(key, allowed_chromosomes=allowed_chromosomes) is not None)


def count_variant_keys(dvars: VarsMap) -> int:
    return len(dvars)


def filter_by_chromosome(dvars: VarsMap, chromosomes: Collection[str]) -> VarsMap:
    allowed = set(chromosomes)
    return {
        key: values
        for key, values in dvars.items()
        if (key.split("_", 1)[0] if key else "") in allowed
    }


def filter_snv_variants(dvars: VarsMap, *, allowed_chromosomes: Collection[str]) -> VarsMap:
    return {
        key: values
        for key, values in dvars.items()
        if parse_snv_key(key, allowed_chromosomes=allowed_chromosomes) is not None
    }


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
    chromosomes: List[str],
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

    # 2) chromosome allowlist (autosomal-only or full CHROMOSOMES set)
    current = filter_by_chromosome(current, chromosomes)
    steps.append(FilterStepCount("chromosome", count_variant_keys(current)))

    # 3) SNV-only variant keys
    if snv_only:
        current = filter_snv_variants(current, allowed_chromosomes=chromosomes)
        steps.append(FilterStepCount("snv_only", count_variant_keys(current)))

    # 4) genotype quality (is_good), optional
    if filter_caps is not None:
        current = filter_by_is_good(current, filter_caps)
        steps.append(FilterStepCount("genotype_qc", count_variant_keys(current)))

    if extra_filters:
        for index, item in enumerate(extra_filters):
            step_name, fn = _normalize_extra_filter(item, index)
            current = fn(current, denovo_cap)
            steps.append(FilterStepCount(step_name, count_variant_keys(current)))

    return FilterPipelineResult(dvars=current, steps=steps)


def split_keys_by_chromosome(dvars: VarsMap, chromosomes: List[str]) -> Dict[str, List[str]]:
    dchr: Dict[str, List[str]] = {chrom: [] for chrom in chromosomes}

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
    """Write headerless per-chromosome VCF from filtered variant keys (no #CHROM line).

    Each variant key is one row. Keys that share a position are not merged.
    """
    for chrom, keys in dchr_sorted.items():
        chrom_dir = output_dir / chrom
        chrom_dir.mkdir(parents=True, exist_ok=True)
        out_file = chrom_dir / "variants.vcf"

        with out_file.open("w", encoding="utf-8") as f:
            for key in keys:
                parts = key.split("_", 3)
                if len(parts) < 4:
                    continue
                row_chrom, pos, ref, alt = parts
                f.write("\t".join([row_chrom, pos, "rsXXX", ref, alt]) + "\n")

        variant_patients = patients_per_variant_for_keys(dvars, keys)
        (chrom_dir / "variant_patients.json").write_text(
            json.dumps(variant_patients, indent=2, sort_keys=True),
            encoding="utf-8",
        )


def is_transition(ref: str, alt: str) -> bool:
    return (ref.upper(), alt.upper()) in TRANSITION_PAIRS


def titv_counts_for_chromosome(
    keys: Iterable[str],
    *,
    allowed_chromosomes: Collection[str],
) -> Dict[str, int | float | None]:
    transitions = 0
    transversions = 0
    for key in keys:
        parsed = parse_snv_key(key, allowed_chromosomes=allowed_chromosomes)
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


def compute_titv_summary(
    dchr_sorted: Dict[str, List[str]],
    chromosomes: List[str],
) -> Dict[str, object]:
    by_chromosome: Dict[str, Dict[str, int | float | None]] = {}
    total_transitions = 0
    total_transversions = 0
    allowed = set(chromosomes)

    for chrom in chromosomes:
        counts = titv_counts_for_chromosome(
            dchr_sorted.get(chrom, []),
            allowed_chromosomes=allowed,
        )
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
    parameters_json: Path
    class_summaries: List[Stage3ClassSummary]


STAGE3_PARAMETERS_FILENAME = "stage3_parameters.json"
STAGE3_PARAMETERS_VERSION = 1


def build_cap_per_class(
    class_map: Dict[str, List[str]],
    default_cap: int | None,
    class_to_cap: Dict[str, int],
) -> Dict[str, int]:
    resolved: Dict[str, int] = {}
    for class_name in class_map:
        cap = class_to_cap.get(class_name, default_cap)
        if cap is not None:
            resolved[class_name] = cap
    return dict(sorted(resolved.items()))


def build_stage3_parameters_payload(
    *,
    stage2_dir: Path,
    output_dir: Path,
    collect: str,
    default_cap: int | None,
    class_to_cap: Dict[str, int],
    cap_per_class: Dict[str, int],
    snv_only: bool,
    autosomal: bool,
    chromosomes: List[str],
    filter_caps: FilterCaps | None,
    class_names: List[str],
    classes_processed: int,
    variants_kept: int,
    write_stats: bool,
) -> Dict[str, object]:
    payload: Dict[str, object] = {
        "version": STAGE3_PARAMETERS_VERSION,
        "stage2_dir": str(stage2_dir.resolve()),
        "output_dir": output_dir.name,
        "output_path": str(output_dir.resolve()),
        "collect": collect,
        "default_cap": default_cap,
        "class_cap_overrides": dict(sorted(class_to_cap.items())),
        "cap_per_class": cap_per_class,
        "snv_only": snv_only,
        "autosomal": autosomal,
        "chromosomes": chromosomes,
        "write_stats": write_stats,
        "filter_enabled": filter_caps is not None,
        "class_names": class_names,
        "classes_processed": classes_processed,
        "variants_kept": variants_kept,
    }
    if filter_caps is not None:
        payload["filter_caps"] = {
            "dp_cap": filter_caps.dp_cap,
            "qt_cap": filter_caps.qt_cap,
            "ab_cap_hom0": filter_caps.ab_cap_hom0,
            "ab_cap_hom1": filter_caps.ab_cap_hom1,
            "ab_cap_het": filter_caps.ab_cap_het,
        }
    return payload


def print_stage3_summary(result: Stage3Result, *, collect: str) -> None:
    print(f"[stage3 summary] collect={collect} output={result.output_dir}")
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
    collect: str,
    dvars_filtered: VarsMap,
    dchr_sorted: Dict[str, List[str]],
    *,
    run_dir_name: str,
    patient_variant_cap: int,
    snv_only: bool,
    autosomal: bool,
    chromosomes: List[str],
    filter_caps: FilterCaps | None,
    filter_steps: List[FilterStepCount],
    patient_ids: List[str],
    write_stats: bool,
) -> None:
    class_dir.mkdir(parents=True, exist_ok=True)
    write_vcf(dvars_filtered, class_dir, dchr_sorted)

    if write_stats:
        write_class_level_stage3_stats(
            class_dir,
            patient_ids=patient_ids,
            dvars_filtered=dvars_filtered,
        )

    titv_summary = compute_titv_summary(dchr_sorted, chromosomes)
    (class_dir / "stage3_titv.json").write_text(
        json.dumps(titv_summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    summary: Dict[str, object] = {
        "collect": collect,
        "class": class_name,
        "output_dir": class_dir.name,
        "run_dir": run_dir_name,
        "variants_kept": len(dvars_filtered),
        "chromosome_counts": {chrom: len(keys) for chrom, keys in dchr_sorted.items()},
        "patient_variant_cap": patient_variant_cap,
        "snv_only": snv_only,
        "autosomal": autosomal,
        "chromosomes": chromosomes,
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
    autosomal: bool = True,
    extra_filters: Iterable[ExtraFilter] | None = None,
    collect: str = "denovo",
    write_stats: bool = True,
) -> Stage3Result:
    if collect not in ("denovo", "inherited"):
        raise ValueError("collect must be 'denovo' or 'inherited'")

    output_dir.mkdir(parents=True, exist_ok=True)
    chromosomes = resolve_stage3_chromosomes(autosomal=autosomal)

    load_combo = load_combined_dvars_for_class if collect == "denovo" else load_combined_dvars_inh_for_class

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
            chromosomes=chromosomes,
            snv_only=snv_only,
            filter_caps=filter_caps,
            extra_filters=extra_filters,
        )
        dvars_filtered = filter_result.dvars
        dchr_class_sorted = sort_chromosome_keys(split_keys_by_chromosome(dvars_filtered, chromosomes))
        titv_class = compute_titv_summary(dchr_class_sorted, chromosomes)
        overall = titv_class["overall"]
        variants_left = len(dvars_filtered)

        class_dir = output_dir / f"class_{class_name}"
        _write_class_stage3_outputs(
            class_dir,
            class_name,
            collect,
            dvars_filtered,
            dchr_class_sorted,
            run_dir_name=output_dir.name,
            patient_variant_cap=cap,
            snv_only=snv_only,
            autosomal=autosomal,
            chromosomes=chromosomes,
            filter_caps=filter_caps,
            filter_steps=filter_result.steps,
            patient_ids=class_map[class_name],
            write_stats=write_stats,
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

    cap_per_class = build_cap_per_class(class_map, default_cap, class_to_cap)
    parameters_json = output_dir / STAGE3_PARAMETERS_FILENAME
    parameters_json.write_text(
        json.dumps(
            build_stage3_parameters_payload(
                stage2_dir=stage2_dir,
                output_dir=output_dir,
                collect=collect,
                default_cap=default_cap,
                class_to_cap=class_to_cap,
                cap_per_class=cap_per_class,
                snv_only=snv_only,
                autosomal=autosomal,
                chromosomes=chromosomes,
                filter_caps=filter_caps,
                class_names=sorted(class_map.keys()),
                classes_processed=classes_processed,
                variants_kept=variants_kept_total,
                write_stats=write_stats,
            ),
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    return Stage3Result(
        classes_processed=classes_processed,
        variants_kept=variants_kept_total,
        output_dir=output_dir,
        parameters_json=parameters_json,
        class_summaries=class_summaries,
    )
