from __future__ import annotations

import gzip
import json
import pickle
from dataclasses import dataclass
from itertools import islice
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Tuple


def get_details(pp: str) -> Tuple[str, int, int, int, int]:
    dpp = pp.split(":")
    if len(dpp) >= 5 and dpp[1] != "." and "." not in dpp[2] and dpp[4] != ".":
        gt = dpp[0]
        dp = int(dpp[1])
        adr = int(dpp[2].split(",")[0])
        ada = int(dpp[2].split(",")[1])
        gq = int(dpp[4])
        return gt, gq, dp, adr, ada
    return "x/x", -1, -1, -1, -1


def get_child_index(line: str, hid: str) -> int:
    dd = line.split("\t")
    if hid == dd[-1]:
        return 1
    if hid == dd[-2]:
        return 2
    if hid == dd[-3]:
        return 3
    return -1


def get_triplet(line: str, ind: int) -> Tuple[str, str, str, str]:
    dd = line.strip().split("\t")
    kk = dd[2]
    if ind == 1:
        return dd[-3], dd[-2], dd[-1], kk
    if ind == 2:
        return dd[-3], dd[-1], dd[-2], kk
    if ind == 3:
        return dd[-2], dd[-1], dd[-3], kk
    return ".", ".", ".", ","


def update(dvars: Dict[str, List[str]], kk: str, gt: str) -> None:
    if kk in dvars:
        dvars[kk].append(gt)
    else:
        dvars[kk] = [gt]


def merge(target: Dict[str, List[str]], source: Dict[str, List[str]]) -> None:
    for kk, values in source.items():
        if kk in target:
            target[kk].extend(values)
        else:
            target[kk] = list(values)


def if_denovo(chl: int, mot: int, fat: int) -> bool:
    return chl > 0 and mot == 0 and fat == 0


def get_vars(file_path: Path, patient_id: str) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    dvars_inh: Dict[str, List[str]] = {}
    dvars: Dict[str, List[str]] = {}
    open_fn = gzip.open if file_path.suffix == ".gz" else open

    with open_fn(file_path, "rt", encoding="utf-8", errors="replace") as f:
        ind = -1
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHR"):
                ind = get_child_index(line.strip(), patient_id)
                continue
            if ind < 0:
                continue

            mother, father, child, kk = get_triplet(line, ind)
            if mother == "." or father == "." or child == "." or ";" in kk:
                continue

            ch, ch_gq, ch_dp, ch_adr, ch_ada = get_details(child)
            mt, mt_gq, mt_dp, mt_adr, mt_ada = get_details(mother)
            ft, ft_gq, ft_dp, ft_adr, ft_ada = get_details(father)
            if ch_dp == -1 or mt_dp == -1 or ft_dp == -1:
                continue

            gt = f"{mt}:{ft}:{ch}"
            mex = ":".join([str(mt_gq), str(mt_dp), str(mt_adr), str(mt_ada)])
            fex = ":".join([str(ft_gq), str(ft_dp), str(ft_adr), str(ft_ada)])
            cex = ":".join([str(ch_gq), str(ch_dp), str(ch_adr), str(ch_ada)])
            gtex = "-".join([patient_id, gt, mex, fex, cex])

            try:
                chl = sum(int(a) > 0 for a in ch.split("/"))
                mot = sum(int(a) > 0 for a in mt.split("/"))
                fat = sum(int(a) > 0 for a in ft.split("/"))
                if if_denovo(chl, mot, fat):
                    update(dvars, kk, gtex)
                else:
                    update(dvars_inh, kk, gtex)
            except Exception:
                continue
    return dvars, dvars_inh


def chunked(items: Iterable[str], size: int) -> Iterator[List[str]]:
    iterator = iter(items)
    while True:
        chunk = list(islice(iterator, size))
        if not chunk:
            return
        yield chunk


def load_stage1_class_map(path: Path) -> Dict[str, List[str]]:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_class_suffix_pairs(pairs: List[str]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for raw in pairs:
        chunks = [chunk.strip() for chunk in raw.split(",") if chunk.strip()]
        for pair in chunks:
            if "=" not in pair:
                raise ValueError(f"Invalid --class-suffix '{pair}'. Expected format class=suffix")
            class_name, suffix = pair.split("=", 1)
            mapping[class_name.strip()] = suffix.strip()
    return mapping


def parse_classes_to_process(entries: List[str]) -> List[str]:
    """Flatten repeatable/comma-separated class names; preserve order, dedupe."""
    out: List[str] = []
    seen: set[str] = set()
    for raw in entries:
        for chunk in raw.split(","):
            name = chunk.strip()
            if name and name not in seen:
                seen.add(name)
                out.append(name)
    return out


@dataclass(frozen=True)
class BatchOutput:
    class_name: str
    batch_index: int
    patient_count: int
    dvars_path: Path
    dvars_inh_path: Path


def build_file_path(input_dir: Path, patient_id: str, suffix: str) -> Path:
    return input_dir / patient_id / f"{patient_id}.{suffix.lstrip('.')}"


def run_stage2(
    input_dir: Path,
    output_dir: Path,
    class_map: Dict[str, List[str]],
    default_suffix: str | None,
    class_to_suffix: Dict[str, str],
    batch_size: int = 1000,
    task: str = "denovo",
    classes_to_process: List[str] | None = None,
) -> List[BatchOutput]:
    if task not in ("denovo", "inherited"):
        raise ValueError("task must be 'denovo' or 'inherited'")

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: List[BatchOutput] = []

    if classes_to_process:
        missing = [c for c in classes_to_process if c not in class_map]
        if missing:
            raise ValueError(
                f"Unknown --classes-to-process: {missing}. Known classes: {sorted(class_map.keys())}"
            )
        class_items = [(k, class_map[k]) for k in classes_to_process]
    else:
        class_items = list(class_map.items())

    for class_name, patient_ids in class_items:
        suffix = class_to_suffix.get(class_name, default_suffix)
        if not suffix:
            raise ValueError(f"No suffix provided for class '{class_name}'")

        class_dir = output_dir / f"class_{class_name}"
        class_dir.mkdir(parents=True, exist_ok=True)

        # Sidecar lists: denovo task -> patients with non-empty inherited; inherited task -> patients with non-empty denovo.
        patients_nonempty_sidecar: set[str] = set()

        for batch_index, batch_patient_ids in enumerate(chunked(patient_ids, batch_size), start=1):
            dVars: Dict[str, List[str]] = {}
            dVars_inh: Dict[str, List[str]] = {}
            batch_patients_nonempty_inh: set[str] = set()
            batch_patients_nonempty_denovo: set[str] = set()

            for patient_id in batch_patient_ids:
                file_path = build_file_path(input_dir, patient_id, suffix)
                if not file_path.is_file():
                    continue
                dvars, dvars_inh = get_vars(file_path, patient_id)
                merge(dVars, dvars)
                merge(dVars_inh, dvars_inh)
                if dvars_inh:
                    batch_patients_nonempty_inh.add(patient_id)
                if dvars:
                    batch_patients_nonempty_denovo.add(patient_id)

            if task == "denovo":
                batch_sidecar = batch_patients_nonempty_inh
                batch_list_name = f"batch_{batch_index:05d}_patients_nonempty_dvars_inh.txt"
                class_list_name = "patients_nonempty_dvars_inh.txt"
            else:
                batch_sidecar = batch_patients_nonempty_denovo
                batch_list_name = f"batch_{batch_index:05d}_patients_nonempty_dvars.txt"
                class_list_name = "patients_nonempty_dvars.txt"

            patients_nonempty_sidecar.update(batch_sidecar)

            dvars_path = class_dir / f"batch_{batch_index:05d}_dVars.pkl"
            dvars_inh_path = class_dir / f"batch_{batch_index:05d}_dVars_inh.pkl"

            with dvars_path.open("wb") as f:
                pickle.dump(dVars, f)
            with dvars_inh_path.open("wb") as f:
                pickle.dump(dVars_inh, f)

            if batch_sidecar:
                (class_dir / batch_list_name).write_text(
                    "\n".join(sorted(batch_sidecar)) + "\n",
                    encoding="utf-8",
                )

            outputs.append(
                BatchOutput(
                    class_name=class_name,
                    batch_index=batch_index,
                    patient_count=len(batch_patient_ids),
                    dvars_path=dvars_path,
                    dvars_inh_path=dvars_inh_path,
                )
            )

        if patients_nonempty_sidecar:
            (class_dir / class_list_name).write_text(
                "\n".join(sorted(patients_nonempty_sidecar)) + "\n",
                encoding="utf-8",
            )

    return outputs
