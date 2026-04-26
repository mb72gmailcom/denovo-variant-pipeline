from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List


@dataclass(frozen=True)
class Stage1Result:
    classes_to_patients: Dict[str, List[str]]
    output_json: Path


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


def run_stage1(input_dir: Path, output_dir: Path, prefixes: List[str]) -> Stage1Result:
    output_dir.mkdir(parents=True, exist_ok=True)

    classes_to_patients = classify_patient_ids(input_dir=input_dir, prefixes=prefixes)
    output_json = output_dir / "stage1_patient_ids_by_class.json"

    with output_json.open("w", encoding="utf-8") as f:
        json.dump(classes_to_patients, f, indent=2, sort_keys=True)

    for class_name, patient_ids in classes_to_patients.items():
        txt_path = output_dir / f"stage1_{class_name}_patient_ids.txt"
        txt_path.write_text("\n".join(patient_ids) + ("\n" if patient_ids else ""), encoding="utf-8")

    return Stage1Result(classes_to_patients=classes_to_patients, output_json=output_json)
