from __future__ import annotations

import gzip
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List


@dataclass(frozen=True)
class Stage1ClassSummary:
    class_name: str
    missed_files: int
    small_files: int
    huge_files: int
    selected_patients: int


@dataclass(frozen=True)
class Stage1Result:
    classes_to_patients: Dict[str, List[str]]
    filtered_classes_to_patients: Dict[str, List[str]]
    class_summaries: List[Stage1ClassSummary]
    output_dir: Path
    output_json: Path
    filtered_output_json: Path
    parameters_json: Path
    stats_json_paths: Dict[str, Path]
    mtime_json_paths: Dict[str, Path]
    missed_txt_paths: Dict[str, Path]
    small_txt_paths: Dict[str, Path]


@dataclass(frozen=True)
class Stage1Parameters:
    version: int
    input_dir: Path
    pipeline_root: Path
    classes: List[str]
    stats: str
    cap_min: int
    cap_max: int
    default_suffix: str | None
    class_suffix_overrides: Dict[str, str]
    suffix_per_class: Dict[str, str]
    suffix_filtering_enabled: bool
    class_names: List[str]
    class_map: str
    filtered_class_map: str

    @classmethod
    def from_dict(cls, data: Dict[str, object]) -> Stage1Parameters:
        version = int(data.get("version", 0))
        if version != STAGE1_PARAMETERS_VERSION:
            raise ValueError(
                f"Unsupported stage1 parameters version {version}; expected {STAGE1_PARAMETERS_VERSION}"
            )
        class_names_raw = data.get("class_names")
        if class_names_raw is not None:
            classes_raw = data.get("classes", [])
        elif "prefixes" in data:
            classes_raw = data["prefixes"]
            class_names_raw = data.get("classes", [])
        else:
            classes_raw = data.get("classes", [])
            class_names_raw = data.get("classes", [])
        return cls(
            version=version,
            input_dir=Path(str(data["input_dir"])),
            pipeline_root=Path(str(data["pipeline_root"])),
            classes=[str(p) for p in classes_raw],
            stats=str(data.get("stats", "size")),
            cap_min=int(data["cap_min"]),
            cap_max=int(data["cap_max"]),
            default_suffix=data.get("default_suffix") if data.get("default_suffix") is not None else None,
            class_suffix_overrides={
                str(k): str(v) for k, v in dict(data.get("class_suffix_overrides", {})).items()
            },
            suffix_per_class={str(k): str(v) for k, v in dict(data.get("suffix_per_class", {})).items()},
            suffix_filtering_enabled=bool(data.get("suffix_filtering_enabled", False)),
            class_names=[str(c) for c in class_names_raw],
            class_map=str(data["class_map"]),
            filtered_class_map=str(data["filtered_class_map"]),
        )


def print_stage1_summary(result: Stage1Result) -> None:
    print("[stage1 summary]")
    for row in result.class_summaries:
        print(
            f"  {row.class_name}: missed_files={row.missed_files} "
            f"small_files={row.small_files} huge_files={row.huge_files} "
            f"selected_patients={row.selected_patients}"
        )


def classify_patient_ids(input_dir: Path, classes: List[str]) -> Dict[str, List[str]]:
    """Collect and sort patient IDs by class prefix."""
    patient_ids = sorted([p.name for p in input_dir.iterdir() if p.is_dir()])

    if not classes:
        return {"all": patient_ids}

    classes_to_patients: Dict[str, List[str]] = {class_prefix: [] for class_prefix in classes}
    classes_to_patients["unmatched"] = []

    for patient_id in patient_ids:
        matched = False
        for class_prefix in classes:
            if patient_id.startswith(class_prefix):
                classes_to_patients[class_prefix].append(patient_id)
                matched = True
                break
        if not matched:
            classes_to_patients["unmatched"].append(patient_id)

    return classes_to_patients


STAGE1_SUBDIR = "stage1"
STAGE2_SUBDIR = "stage2"
CLASS_MAP_FILENAME = "stage1_patient_ids_by_class.json"
PARAMETERS_FILENAME = "stage1_parameters.json"
DEFAULT_CAP_MIN = 22000
DEFAULT_CAP_MAX = 75000
STAGE1_PARAMETERS_VERSION = 1


def filtered_class_map_filename(cap_min: int, cap_max: int) -> str:
    return f"stage1_patient_ids_by_class_filtered_{cap_min}_{cap_max}.json"


def stage1_class_dir(stage1_dir: Path, class_name: str) -> Path:
    return stage1_dir / f"class_{class_name}"


def stage1_class_map_path(pipeline_root: Path) -> Path:
    return pipeline_root / STAGE1_SUBDIR / CLASS_MAP_FILENAME


def stage1_filtered_class_map_path(pipeline_root: Path, cap_min: int, cap_max: int) -> Path:
    return pipeline_root / STAGE1_SUBDIR / filtered_class_map_filename(cap_min, cap_max)


def stage1_parameters_path(pipeline_root: Path) -> Path:
    return pipeline_root / STAGE1_SUBDIR / PARAMETERS_FILENAME


def resolve_class_suffix(
    class_name: str,
    default_suffix: str | None,
    class_to_suffix: Dict[str, str],
) -> str | None:
    return class_to_suffix.get(class_name, default_suffix)


def suffix_filtering_enabled(
    default_suffix: str | None,
    class_to_suffix: Dict[str, str],
) -> bool:
    return default_suffix is not None or bool(class_to_suffix)


def build_suffix_per_class(
    classes: Iterable[str],
    default_suffix: str | None,
    class_to_suffix: Dict[str, str],
) -> Dict[str, str]:
    resolved: Dict[str, str] = {}
    for class_name in classes:
        suffix = resolve_class_suffix(class_name, default_suffix, class_to_suffix)
        if suffix is not None:
            resolved[class_name] = suffix
    return dict(sorted(resolved.items()))


def build_stage1_parameters_payload(
    *,
    input_dir: Path,
    pipeline_root: Path,
    classes: List[str],
    stats: str,
    cap_min: int,
    cap_max: int,
    default_suffix: str | None,
    class_to_suffix: Dict[str, str],
    class_names: List[str],
) -> Dict[str, object]:
    suffix_per_class = build_suffix_per_class(class_names, default_suffix, class_to_suffix)
    return {
        "version": STAGE1_PARAMETERS_VERSION,
        "input_dir": str(input_dir.resolve()),
        "pipeline_root": str(pipeline_root.resolve()),
        "classes": classes,
        "stats": stats,
        "cap_min": cap_min,
        "cap_max": cap_max,
        "default_suffix": default_suffix,
        "class_suffix_overrides": dict(sorted(class_to_suffix.items())),
        "suffix_per_class": suffix_per_class,
        "suffix_filtering_enabled": suffix_filtering_enabled(default_suffix, class_to_suffix),
        "class_names": class_names,
        "class_map": CLASS_MAP_FILENAME,
        "filtered_class_map": filtered_class_map_filename(cap_min, cap_max),
    }


def try_load_stage1_parameters(pipeline_root: Path) -> Stage1Parameters | None:
    path = stage1_parameters_path(pipeline_root)
    if not path.is_file():
        return None
    return Stage1Parameters.from_dict(json.loads(path.read_text(encoding="utf-8")))


def load_stage1_parameters(
    *,
    explicit: Path | None = None,
    pipeline_root: Path | None = None,
) -> Stage1Parameters:
    if explicit is not None:
        path = explicit
    elif pipeline_root is not None:
        path = stage1_parameters_path(pipeline_root)
    else:
        raise ValueError("Provide explicit parameters path or pipeline_root")
    if not path.is_file():
        raise FileNotFoundError(
            f"Stage1 parameters not found at '{path}'. Run stage1 or pass --stage1-parameters-json explicitly."
        )
    return Stage1Parameters.from_dict(json.loads(path.read_text(encoding="utf-8")))


def resolve_file_suffix_for_class(
    class_name: str,
    *,
    default_suffix: str | None = None,
    class_to_suffix: Dict[str, str] | None = None,
    suffix_per_class: Dict[str, str] | None = None,
) -> str | None:
    """CLI overrides first, then stage1 suffix_per_class fallback."""
    class_to_suffix = class_to_suffix or {}
    if class_name in class_to_suffix:
        return class_to_suffix[class_name]
    if default_suffix is not None:
        return default_suffix
    if suffix_per_class and class_name in suffix_per_class:
        return suffix_per_class[class_name]
    return None


@dataclass(frozen=True)
class Stage2RunConfig:
    input_dir: Path
    suffix_per_class: Dict[str, str]
    classes: List[str]


def resolve_stage2_run_config(
    pipeline_root: Path,
    class_map: Dict[str, List[str]],
    classes_override: List[str] | None = None,
    *,
    stage1_params: Stage1Parameters | None = None,
) -> Stage2RunConfig:
    """
    Load input_dir, suffixes, and class list for stage2 from stage1 parameters.

    When classes_override is empty/None, use class_names from stage1 parameters
    restricted to keys present in class_map.

    Pass stage1_params when stage1 ran in the same session (avoids re-read from disk).
    """
    params = stage1_params
    if params is None:
        params = try_load_stage1_parameters(pipeline_root)
    if params is None:
        path = stage1_parameters_path(pipeline_root)
        raise FileNotFoundError(
            f"Stage1 parameters not found at '{path}'. Run stage1 before stage2."
        )

    return _resolve_stage2_run_config_from_params(params, class_map, classes_override)


def _resolve_stage2_run_config_from_params(
    params: Stage1Parameters,
    class_map: Dict[str, List[str]],
    classes_override: List[str] | None,
) -> Stage2RunConfig:
    suffix_per_class = dict(params.suffix_per_class)
    if not suffix_per_class:
        raise ValueError(
            "stage1_parameters.json has no suffix_per_class. "
            "Run stage1 with --suffix and/or --class-suffix."
        )

    if classes_override:
        missing_map = [c for c in classes_override if c not in class_map]
        if missing_map:
            raise ValueError(
                f"Unknown --classes: {missing_map}. Known classes: {sorted(class_map.keys())}"
            )
        missing_suffix = [c for c in classes_override if c not in suffix_per_class]
        if missing_suffix:
            raise ValueError(
                f"No suffix in stage1 parameters for class(es): {missing_suffix}"
            )
        return Stage2RunConfig(
            input_dir=params.input_dir,
            suffix_per_class=suffix_per_class,
            classes=classes_override,
        )

    default_classes = [c for c in params.class_names if c in class_map]
    if not default_classes:
        default_classes = sorted(class_map.keys())

    missing_suffix = [c for c in default_classes if c not in suffix_per_class]
    if missing_suffix:
        raise ValueError(
            f"No suffix in stage1 parameters for class(es): {missing_suffix}"
        )
    return Stage2RunConfig(
        input_dir=params.input_dir,
        suffix_per_class=suffix_per_class,
        classes=default_classes,
    )


def load_stage1_parameters_from_result(result: Stage1Result) -> Stage1Parameters:
    return Stage1Parameters.from_dict(
        json.loads(result.parameters_json.read_text(encoding="utf-8"))
    )


def subset_class_map(class_map: Dict[str, List[str]], classes: List[str]) -> Dict[str, List[str]]:
    missing = [c for c in classes if c not in class_map]
    if missing:
        raise ValueError(
            f"Unknown class(es): {missing}. Known classes: {sorted(class_map.keys())}"
        )
    return {class_name: class_map[class_name] for class_name in classes}


def effective_stage1_caps(
    pipeline_root: Path,
    cap_min: int,
    cap_max: int,
) -> tuple[int, int]:
    """Use caps recorded in stage1_parameters.json when that file exists."""
    params = try_load_stage1_parameters(pipeline_root)
    if params is None:
        return cap_min, cap_max
    return params.cap_min, params.cap_max


def infer_pipeline_root(path: Path) -> Path:
    """Return pipeline root from a stage output directory or the root itself."""
    if path.name in (STAGE1_SUBDIR, STAGE2_SUBDIR) or path.name.startswith("stage3"):
        return path.parent
    return path


def resolve_class_map_json(
    explicit: Path | None,
    pipeline_root: Path,
    cap_min: int = DEFAULT_CAP_MIN,
    cap_max: int = DEFAULT_CAP_MAX,
) -> Path:
    cap_min, cap_max = effective_stage1_caps(pipeline_root, cap_min, cap_max)
    path = (
        explicit
        if explicit is not None
        else stage1_filtered_class_map_path(pipeline_root, cap_min, cap_max)
    )
    if not path.is_file():
        raise FileNotFoundError(
            f"Stage1 class map not found at '{path}'. "
            "Run stage1 or pass --class-map-json explicitly."
        )
    return path


def load_class_map(
    explicit: Path | None = None,
    pipeline_root: Path | None = None,
    cap_min: int = DEFAULT_CAP_MIN,
    cap_max: int = DEFAULT_CAP_MAX,
) -> Dict[str, List[str]]:
    if explicit is not None:
        path = explicit
    elif pipeline_root is not None:
        cap_min, cap_max = effective_stage1_caps(pipeline_root, cap_min, cap_max)
        path = stage1_filtered_class_map_path(pipeline_root, cap_min, cap_max)
    else:
        raise ValueError("Provide explicit class map path or pipeline_root")
    if not path.is_file():
        raise FileNotFoundError(
            f"Stage1 class map not found at '{path}'. "
            "Run stage1 or pass --class-map-json explicitly."
        )
    return json.loads(path.read_text(encoding="utf-8"))


def resolve_stage2_output_dir(output_dir: Path) -> Path:
    if output_dir.name == STAGE2_SUBDIR:
        return output_dir
    return output_dir / STAGE2_SUBDIR


def resolve_stage2_input_dir(explicit: Path | None, pipeline_root: Path) -> Path:
    if explicit is not None:
        return explicit
    return pipeline_root / STAGE2_SUBDIR


def patient_filename(patient_id: str, suffix: str) -> str:
    """Build <patient_id>.<suffix>, or append suffix as-is when it starts with '.' or '_'."""
    suffix = suffix.strip()
    if not suffix:
        raise ValueError("suffix must be non-empty")
    if suffix[0] in "._":
        return f"{patient_id}{suffix}"
    return f"{patient_id}.{suffix}"


def build_patient_file_path(input_dir: Path, patient_id: str, suffix: str) -> Path:
    return input_dir / patient_id / patient_filename(patient_id, suffix)


def file_stat_value(file_path: Path, stats: str) -> int:
    if stats == "size":
        return file_path.stat().st_size
    if stats == "counts":
        open_fn = gzip.open if file_path.suffix == ".gz" else open
        with open_fn(file_path, "rt", encoding="utf-8", errors="replace") as f:
            return sum(1 for _ in f)
    raise ValueError("stats must be 'size' or 'counts'")


def process_patients_for_suffix(
    input_dir: Path,
    patient_ids: List[str],
    suffix: str,
    stats: str,
    cap_min: int,
    cap_max: int,
) -> tuple[Dict[str, int], Dict[str, int], List[str], List[str], List[str], Dict[str, int]]:
    """One pass: stats, mtime, missed, small, huge, and byte sizes for existing files."""
    stats_dict: Dict[str, int] = {}
    mtime_dict: Dict[str, int] = {}
    sizes: Dict[str, int] = {}
    missed: List[str] = []
    small: List[str] = []
    huge: List[str] = []

    for patient_id in patient_ids:
        file_path = build_patient_file_path(input_dir, patient_id, suffix)
        if not file_path.is_file():
            missed.append(patient_id)
            continue

        st = file_path.stat()
        size = st.st_size
        sizes[patient_id] = size
        mtime_dict[patient_id] = int(st.st_mtime)

        if stats == "size":
            stats_dict[patient_id] = size
        else:
            stats_dict[patient_id] = file_stat_value(file_path, stats)

        if size < cap_min:
            small.append(patient_id)
        elif size > cap_max:
            huge.append(patient_id)

    return stats_dict, mtime_dict, missed, small, huge, sizes


def filter_class_by_sizes(
    patient_ids: List[str],
    sizes: Dict[str, int],
    cap_min: int,
    cap_max: int,
) -> List[str]:
    return [
        patient_id
        for patient_id in patient_ids
        if patient_id in sizes and cap_min <= sizes[patient_id] <= cap_max
    ]


def run_stage1(
    input_dir: Path,
    output_dir: Path,
    classes: List[str],
    stats: str = "size",
    default_suffix: str | None = None,
    class_to_suffix: Dict[str, str] | None = None,
    cap_min: int = DEFAULT_CAP_MIN,
    cap_max: int = DEFAULT_CAP_MAX,
) -> Stage1Result:
    if stats not in ("size", "counts"):
        raise ValueError("stats must be 'size' or 'counts'")

    class_to_suffix = class_to_suffix or {}

    stage1_dir = output_dir / "stage1"
    stage1_dir.mkdir(parents=True, exist_ok=True)

    classes_to_patients = classify_patient_ids(input_dir=input_dir, classes=classes)
    output_json = stage1_dir / "stage1_patient_ids_by_class.json"

    with output_json.open("w", encoding="utf-8") as f:
        json.dump(classes_to_patients, f, indent=2, sort_keys=True)

    for class_name, patient_ids in classes_to_patients.items():
        class_dir = stage1_class_dir(stage1_dir, class_name)
        class_dir.mkdir(parents=True, exist_ok=True)
        (class_dir / "patient_ids.txt").write_text(
            "\n".join(patient_ids) + ("\n" if patient_ids else ""),
            encoding="utf-8",
        )

    stats_json_paths: Dict[str, Path] = {}
    mtime_json_paths: Dict[str, Path] = {}
    missed_txt_paths: Dict[str, Path] = {}
    small_txt_paths: Dict[str, Path] = {}
    filtered_classes: Dict[str, List[str]] = {}
    class_summaries: List[Stage1ClassSummary] = []

    if suffix_filtering_enabled(default_suffix, class_to_suffix):
        for class_name, patient_ids in classes_to_patients.items():
            suffix = resolve_class_suffix(class_name, default_suffix, class_to_suffix)
            if not suffix:
                raise ValueError(f"No suffix provided for class '{class_name}'")

            class_dir = stage1_class_dir(stage1_dir, class_name)
            class_dir.mkdir(parents=True, exist_ok=True)
            stats_dict, mtime_dict, missed, small, huge, sizes = process_patients_for_suffix(
                input_dir, patient_ids, suffix, stats, cap_min, cap_max
            )

            stats_path = class_dir / f"stats_{stats}.json"
            stats_path.write_text(json.dumps(stats_dict, indent=2, sort_keys=True), encoding="utf-8")
            stats_json_paths[class_name] = stats_path

            mtime_path = class_dir / "stats_mtime.json"
            mtime_path.write_text(json.dumps(mtime_dict, indent=2, sort_keys=True), encoding="utf-8")
            mtime_json_paths[class_name] = mtime_path

            missed_path = class_dir / "missed.txt"
            missed_path.write_text(
                "\n".join(sorted(missed)) + ("\n" if missed else ""),
                encoding="utf-8",
            )
            missed_txt_paths[class_name] = missed_path

            small_path = class_dir / f"small_{cap_min}.txt"
            small_path.write_text(
                "\n".join(sorted(small)) + ("\n" if small else ""),
                encoding="utf-8",
            )
            small_txt_paths[class_name] = small_path

            filtered_classes[class_name] = filter_class_by_sizes(
                patient_ids, sizes, cap_min, cap_max
            )
            class_summaries.append(
                Stage1ClassSummary(
                    class_name=class_name,
                    missed_files=len(missed),
                    small_files=len(small),
                    huge_files=len(huge),
                    selected_patients=len(filtered_classes[class_name]),
                )
            )
    else:
        filtered_classes = {
            class_name: list(patient_ids) for class_name, patient_ids in classes_to_patients.items()
        }
        for class_name, patient_ids in classes_to_patients.items():
            class_summaries.append(
                Stage1ClassSummary(
                    class_name=class_name,
                    missed_files=0,
                    small_files=0,
                    huge_files=0,
                    selected_patients=len(patient_ids),
                )
            )

    filtered_output_json = stage1_dir / filtered_class_map_filename(cap_min, cap_max)
    with filtered_output_json.open("w", encoding="utf-8") as f:
        json.dump(filtered_classes, f, indent=2, sort_keys=True)

    parameters_json = stage1_dir / PARAMETERS_FILENAME
    parameters_payload = build_stage1_parameters_payload(
        input_dir=input_dir,
        pipeline_root=output_dir,
        classes=classes,
        stats=stats,
        cap_min=cap_min,
        cap_max=cap_max,
        default_suffix=default_suffix,
        class_to_suffix=class_to_suffix,
        class_names=sorted(classes_to_patients.keys()),
    )
    parameters_json.write_text(
        json.dumps(parameters_payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    return Stage1Result(
        classes_to_patients=classes_to_patients,
        filtered_classes_to_patients=filtered_classes,
        class_summaries=class_summaries,
        output_dir=stage1_dir,
        output_json=output_json,
        filtered_output_json=filtered_output_json,
        parameters_json=parameters_json,
        stats_json_paths=stats_json_paths,
        mtime_json_paths=mtime_json_paths,
        missed_txt_paths=missed_txt_paths,
        small_txt_paths=small_txt_paths,
    )
