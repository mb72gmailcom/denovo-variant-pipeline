from __future__ import annotations

import gzip
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List


@dataclass(frozen=True)
class Stage1Result:
    classes_to_patients: Dict[str, List[str]]
    filtered_classes_to_patients: Dict[str, List[str]]
    output_dir: Path
    output_json: Path
    filtered_output_json: Path
    stats_json_paths: Dict[str, Path]
    mtime_json_paths: Dict[str, Path]
    missed_txt_paths: Dict[str, Path]
    small_txt_paths: Dict[str, Path]


def classify_patient_ids(input_dir: Path, prefixes: List[str]) -> Dict[str, List[str]]:
    """Collect and sort patient IDs by class prefix."""
    patient_ids = sorted([p.name for p in input_dir.iterdir() if p.is_dir()])

    if not prefixes:
        return {"all": patient_ids}

    classes_to_patients: Dict[str, List[str]] = {prefix: [] for prefix in prefixes}
    classes_to_patients["unmatched"] = []

    for patient_id in patient_ids:
        matched = False
        for prefix in prefixes:
            if patient_id.startswith(prefix):
                classes_to_patients[prefix].append(patient_id)
                matched = True
                break
        if not matched:
            classes_to_patients["unmatched"].append(patient_id)

    return classes_to_patients


STAGE1_SUBDIR = "stage1"
STAGE2_SUBDIR = "stage2"
CLASS_MAP_FILENAME = "stage1_patient_ids_by_class.json"
DEFAULT_CAP_MIN = 22000
DEFAULT_CAP_MAX = 75000


def filtered_class_map_filename(cap_min: int, cap_max: int) -> str:
    return f"stage1_patient_ids_by_class_filtered_{cap_min}_{cap_max}.json"


def stage1_class_map_path(pipeline_root: Path) -> Path:
    return pipeline_root / STAGE1_SUBDIR / CLASS_MAP_FILENAME


def stage1_filtered_class_map_path(pipeline_root: Path, cap_min: int, cap_max: int) -> Path:
    return pipeline_root / STAGE1_SUBDIR / filtered_class_map_filename(cap_min, cap_max)


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


def parse_suffixes(value: str) -> List[str]:
    if not value.strip():
        return []
    return [s.strip() for s in value.split(",") if s.strip()]


def suffix_label(suffix: str) -> str:
    return suffix.lstrip(".").replace(".", "_").replace("/", "_")


def build_patient_file_path(input_dir: Path, patient_id: str, suffix: str) -> Path:
    return input_dir / patient_id / f"{patient_id}.{suffix.lstrip('.')}"


def file_stat_value(file_path: Path, stats: str) -> int:
    if stats == "size":
        return file_path.stat().st_size
    if stats == "counts":
        open_fn = gzip.open if file_path.suffix == ".gz" else open
        with open_fn(file_path, "rt", encoding="utf-8", errors="replace") as f:
            return sum(1 for _ in f)
    raise ValueError("stats must be 'size' or 'counts'")


def compute_file_metadata_by_patient(
    input_dir: Path,
    patient_ids: List[str],
    suffix: str,
    stats: str,
) -> tuple[Dict[str, int], Dict[str, int], List[str]]:
    stats_dict: Dict[str, int] = {}
    mtime_dict: Dict[str, int] = {}
    missed: List[str] = []
    for patient_id in patient_ids:
        file_path = build_patient_file_path(input_dir, patient_id, suffix)
        if not file_path.is_file():
            missed.append(patient_id)
            continue
        st = file_path.stat()
        mtime_dict[patient_id] = int(st.st_mtime)
        if stats == "size":
            stats_dict[patient_id] = st.st_size
        else:
            stats_dict[patient_id] = file_stat_value(file_path, stats)
    return stats_dict, mtime_dict, missed


def compute_stats_by_patient(
    input_dir: Path,
    patient_ids: List[str],
    suffix: str,
    stats: str,
) -> tuple[Dict[str, int], List[str]]:
    stats_dict, _, missed = compute_file_metadata_by_patient(input_dir, patient_ids, suffix, stats)
    return stats_dict, missed


def patient_file_size_in_range(
    input_dir: Path,
    patient_id: str,
    suffix: str,
    cap_min: int,
    cap_max: int,
) -> bool:
    file_path = build_patient_file_path(input_dir, patient_id, suffix)
    if not file_path.is_file():
        return False
    size = file_path.stat().st_size
    return cap_min <= size <= cap_max


def patients_with_file_smaller_than_cap(
    input_dir: Path,
    patient_ids: List[str],
    suffix: str,
    cap_min: int,
) -> List[str]:
    small: List[str] = []
    for patient_id in patient_ids:
        file_path = build_patient_file_path(input_dir, patient_id, suffix)
        if not file_path.is_file():
            continue
        if file_path.stat().st_size < cap_min:
            small.append(patient_id)
    return sorted(small)


def filter_classes_by_file_size(
    input_dir: Path,
    classes_to_patients: Dict[str, List[str]],
    suffixes: List[str],
    cap_min: int,
    cap_max: int,
) -> Dict[str, List[str]]:
    if not suffixes:
        return {class_name: list(patient_ids) for class_name, patient_ids in classes_to_patients.items()}

    filtered: Dict[str, List[str]] = {}
    for class_name, patient_ids in classes_to_patients.items():
        kept = [
            patient_id
            for patient_id in patient_ids
            if all(
                patient_file_size_in_range(input_dir, patient_id, suffix, cap_min, cap_max)
                for suffix in suffixes
            )
        ]
        filtered[class_name] = kept
    return filtered


def run_stage1(
    input_dir: Path,
    output_dir: Path,
    prefixes: List[str],
    stats: str = "size",
    suffixes: List[str] | None = None,
    cap_min: int = DEFAULT_CAP_MIN,
    cap_max: int = DEFAULT_CAP_MAX,
) -> Stage1Result:
    if stats not in ("size", "counts"):
        raise ValueError("stats must be 'size' or 'counts'")

    stage1_dir = output_dir / "stage1"
    stage1_dir.mkdir(parents=True, exist_ok=True)

    classes_to_patients = classify_patient_ids(input_dir=input_dir, prefixes=prefixes)
    output_json = stage1_dir / "stage1_patient_ids_by_class.json"

    with output_json.open("w", encoding="utf-8") as f:
        json.dump(classes_to_patients, f, indent=2, sort_keys=True)

    for class_name, patient_ids in classes_to_patients.items():
        txt_path = stage1_dir / f"stage1_{class_name}_patient_ids.txt"
        txt_path.write_text("\n".join(patient_ids) + ("\n" if patient_ids else ""), encoding="utf-8")

    stats_json_paths: Dict[str, Path] = {}
    mtime_json_paths: Dict[str, Path] = {}
    missed_txt_paths: Dict[str, Path] = {}
    small_txt_paths: Dict[str, Path] = {}
    if suffixes:
        for class_name, patient_ids in classes_to_patients.items():
            for suffix in suffixes:
                label = f"{class_name}_{suffix_label(suffix)}"
                stats_dict, mtime_dict, missed = compute_file_metadata_by_patient(
                    input_dir, patient_ids, suffix, stats
                )

                stats_path = stage1_dir / f"stage1_stats_{stats}_{label}.json"
                stats_path.write_text(json.dumps(stats_dict, indent=2, sort_keys=True), encoding="utf-8")
                stats_json_paths[label] = stats_path

                mtime_path = stage1_dir / f"stage1_stats_mtime_{label}.json"
                mtime_path.write_text(json.dumps(mtime_dict, indent=2, sort_keys=True), encoding="utf-8")
                mtime_json_paths[label] = mtime_path

                missed_path = stage1_dir / f"{label}.missed.txt"
                missed_path.write_text("\n".join(sorted(missed)) + ("\n" if missed else ""), encoding="utf-8")
                missed_txt_paths[label] = missed_path

                small_patients = patients_with_file_smaller_than_cap(
                    input_dir, patient_ids, suffix, cap_min
                )
                small_path = stage1_dir / f"{label}.small_{cap_min}.txt"
                small_path.write_text(
                    "\n".join(small_patients) + ("\n" if small_patients else ""),
                    encoding="utf-8",
                )
                small_txt_paths[label] = small_path

    filtered_classes = filter_classes_by_file_size(
        input_dir, classes_to_patients, suffixes or [], cap_min, cap_max
    )
    filtered_output_json = stage1_dir / filtered_class_map_filename(cap_min, cap_max)
    with filtered_output_json.open("w", encoding="utf-8") as f:
        json.dump(filtered_classes, f, indent=2, sort_keys=True)

    return Stage1Result(
        classes_to_patients=classes_to_patients,
        filtered_classes_to_patients=filtered_classes,
        output_dir=stage1_dir,
        output_json=output_json,
        filtered_output_json=filtered_output_json,
        stats_json_paths=stats_json_paths,
        mtime_json_paths=mtime_json_paths,
        missed_txt_paths=missed_txt_paths,
        small_txt_paths=small_txt_paths,
    )
