from __future__ import annotations

import csv
import gzip
import json
import pickle
import statistics
from dataclasses import dataclass, field
from itertools import islice
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Tuple

from .stage1 import load_class_map

HIST_BIN_COUNT = 100
AB_GT_KEYS = ("00", "01", "11")
GT_TO_AB_KEY = {"0/0": "00", "0/1": "01", "1/1": "11"}
AB_KEY_TO_GT = {"00": "0/0", "01": "0/1", "11": "1/1"}


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


@dataclass(frozen=True)
class FamilyColumns:
    child: str
    mother: str
    father: str
    map: str | None = None


@dataclass(frozen=True)
class Trio:
    child: str
    mother: str
    father: str
    child_sample: str
    mother_sample: str
    father_sample: str


def _load_family_columns_data(spec: str | Path) -> tuple[Dict[str, object], Path | None]:
    text = str(spec).strip()
    if not text:
        raise ValueError("Family column spec is empty")
    if text.startswith("{"):
        try:
            data = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid inline family column JSON: {exc}") from exc
        if not isinstance(data, dict):
            raise ValueError("Family column JSON must be an object")
        return data, None

    path = Path(text)
    if not path.is_file():
        raise FileNotFoundError(
            f"Family column JSON not found at '{path}'. "
            "Pass a JSON file path or an inline object such as "
            '\'{"child": "spid", "mother": "mother_id", "father": "father_id"}\'.'
        )
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid family column JSON in '{path}': {exc}") from exc
    if not isinstance(data, dict):
        raise ValueError("Family column JSON must be an object")
    return data, path


def load_family_columns(spec: str | Path, *, qmap: bool) -> tuple[FamilyColumns, Path | None]:
    data, source_path = _load_family_columns_data(spec)
    missing = [key for key in ("child", "mother", "father") if not str(data.get(key, "")).strip()]
    if missing:
        raise ValueError(f"Family column JSON missing keys: {missing}")
    map_col = str(data["map"]).strip() if data.get("map") is not None else ""
    if qmap and not map_col:
        raise ValueError("Family column JSON must include 'map' when --qmap is set")
    return (
        FamilyColumns(
            child=str(data["child"]).strip(),
            mother=str(data["mother"]).strip(),
            father=str(data["father"]).strip(),
            map=map_col if qmap else None,
        ),
        source_path,
    )


def _family_delimiter(header_line: str) -> str:
    return "\t" if "\t" in header_line else ","


def _family_id(value: object) -> str:
    text = str(value).strip() if value is not None else ""
    if text in {"", ".", "NA", "na", "NaN", "nan", "None", "none"}:
        return ""
    return text


def load_trios(family_file: Path, columns: FamilyColumns, *, qmap: bool) -> Dict[str, Trio]:
    if not family_file.is_file():
        raise FileNotFoundError(f"Family file not found at '{family_file}'")

    with family_file.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        first = handle.readline()
        if not first:
            raise ValueError(f"Family file is empty: '{family_file}'")
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter=_family_delimiter(first))
        if reader.fieldnames is None:
            raise ValueError(f"Family file has no header: '{family_file}'")
        needed = [columns.child, columns.mother, columns.father]
        if qmap and columns.map:
            needed.append(columns.map)
        missing_cols = [name for name in needed if name not in reader.fieldnames]
        if missing_cols:
            raise ValueError(
                f"Family file missing column(s) {missing_cols}. Found: {list(reader.fieldnames)}"
            )
        rows = list(reader)

    person_to_sample: Dict[str, str] = {}
    if qmap and columns.map:
        for row in rows:
            person = _family_id(row.get(columns.child))
            sample = _family_id(row.get(columns.map))
            if not person or not sample:
                continue
            previous = person_to_sample.get(person)
            if previous is not None and previous != sample:
                raise ValueError(
                    f"Conflicting sample IDs for '{person}': {previous!r} vs {sample!r}"
                )
            person_to_sample[person] = sample

    trios: Dict[str, Trio] = {}
    for row in rows:
        child = _family_id(row.get(columns.child))
        mother = _family_id(row.get(columns.mother))
        father = _family_id(row.get(columns.father))
        if not child or not mother or not father:
            continue
        if qmap:
            child_sample = person_to_sample.get(child, "")
            mother_sample = person_to_sample.get(mother, "")
            father_sample = person_to_sample.get(father, "")
            if not child_sample or not mother_sample or not father_sample:
                continue
        else:
            child_sample, mother_sample, father_sample = child, mother, father

        trio = Trio(
            child=child,
            mother=mother,
            father=father,
            child_sample=child_sample,
            mother_sample=mother_sample,
            father_sample=father_sample,
        )
        existing = trios.get(child)
        if existing is not None and existing != trio:
            raise ValueError(
                f"Conflicting trio rows for child '{child}': "
                f"({existing.mother}, {existing.father}) vs ({mother}, {father})"
            )
        trios[child] = trio
    return trios


def trio_column_indices(header_fields: List[str], trio: Trio) -> Tuple[int, int, int] | None:
    try:
        return (
            header_fields.index(trio.mother_sample),
            header_fields.index(trio.father_sample),
            header_fields.index(trio.child_sample),
        )
    except ValueError:
        return None


def get_triplet_at(
    line: str,
    mother_i: int,
    father_i: int,
    child_i: int,
) -> Tuple[str, str, str, str]:
    dd = line.strip().split("\t")
    kk = dd[2] if len(dd) > 2 else ","
    if max(mother_i, father_i, child_i) >= len(dd):
        return ".", ".", ".", ","
    return dd[mother_i], dd[father_i], dd[child_i], kk


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


def _empty_hist_bins() -> List[int]:
    return [0] * HIST_BIN_COUNT


def _bin_integer_metric(value: int, counts: List[int]) -> None:
    if value < 1:
        return
    if value >= 100:
        counts[99] += 1
    else:
        counts[value - 1] += 1


def _bin_ab_ratio(adr: int, ada: int, counts: List[int]) -> None:
    denom = adr + ada
    if denom <= 0:
        return
    ab = ada / denom
    if ab >= 0.99:
        counts[99] += 1
    elif ab >= 0:
        counts[int(ab / 0.01)] += 1


@dataclass
class Stage2HistCollector:
    """Per-class histograms for variants in the active stage2 collect bucket."""

    children_qt: List[int] = field(default_factory=_empty_hist_bins)
    parents_qt: List[int] = field(default_factory=_empty_hist_bins)
    children_dp: List[int] = field(default_factory=_empty_hist_bins)
    parents_dp: List[int] = field(default_factory=_empty_hist_bins)
    children_ab: Dict[str, List[int]] = field(
        default_factory=lambda: {key: _empty_hist_bins() for key in AB_GT_KEYS}
    )
    parents_ab: Dict[str, List[int]] = field(
        default_factory=lambda: {key: _empty_hist_bins() for key in AB_GT_KEYS}
    )

    def record_parent(self, gt: str, gq: int, dp: int, adr: int, ada: int) -> None:
        _bin_integer_metric(gq, self.parents_qt)
        _bin_integer_metric(dp, self.parents_dp)
        ab_key = GT_TO_AB_KEY.get(gt)
        if ab_key is not None:
            _bin_ab_ratio(adr, ada, self.parents_ab[ab_key])

    def record_child(self, gt: str, gq: int, dp: int, adr: int, ada: int) -> None:
        _bin_integer_metric(gq, self.children_qt)
        _bin_integer_metric(dp, self.children_dp)
        ab_key = GT_TO_AB_KEY.get(gt)
        if ab_key is not None:
            _bin_ab_ratio(adr, ada, self.children_ab[ab_key])

    def write_json_files(self, class_dir: Path, collect: str) -> List[Path]:
        written: List[Path] = []
        specs = [
            ("children_qt_hist.json", "children", "qt", None, self.children_qt),
            ("parents_qt_hist.json", "parents", "qt", None, self.parents_qt),
            ("children_dp_hist.json", "children", "dp", None, self.children_dp),
            ("parents_dp_hist.json", "parents", "dp", None, self.parents_dp),
        ]
        for filename, role, metric, genotype, counts in specs:
            if sum(counts) == 0:
                continue
            path = class_dir / filename
            path.write_text(
                json.dumps(
                    _hist_payload(collect, role, metric, genotype, counts, integer=True),
                    indent=2,
                    sort_keys=True,
                ),
                encoding="utf-8",
            )
            written.append(path)

        for ab_suffix in AB_GT_KEYS:
            for role, ab_map in (("children", self.children_ab), ("parents", self.parents_ab)):
                counts = ab_map[ab_suffix]
                if sum(counts) == 0:
                    continue
                filename = f"{role}_{ab_suffix}_ab_hist.json"
                path = class_dir / filename
                path.write_text(
                    json.dumps(
                        _hist_payload(
                            collect,
                            role,
                            "ab",
                            ab_suffix,
                            counts,
                            integer=False,
                        ),
                        indent=2,
                        sort_keys=True,
                    ),
                    encoding="utf-8",
                )
                written.append(path)
        return written


def _hist_payload(
    collect: str,
    role: str,
    metric: str,
    genotype: str | None,
    counts: List[int],
    *,
    integer: bool,
) -> Dict[str, object]:
    payload: Dict[str, object] = {
        "schema": "stage2_hist_v1",
        "collect": collect,
        "role": role,
        "metric": metric,
        "counts": counts,
    }
    if genotype is not None:
        payload["genotype"] = AB_KEY_TO_GT[genotype]
    if integer:
        payload["bin_spec"] = {
            "type": "integer",
            "buckets": HIST_BIN_COUNT,
            "labels": "1..99 map to index 0..98; index 99 is >=100",
        }
    else:
        payload["bin_spec"] = {
            "type": "ab_ratio",
            "buckets": HIST_BIN_COUNT,
            "width": 0.01,
            "last_bin": "[0.99, 1.0]",
            "other_bins": "[k*0.01, (k+1)*0.01) for k=0..98",
            "ab_formula": "ada/(adr+ada)",
        }
    return payload


def if_denovo_ext(child: int, mother: int, father: int) -> bool:
    if child == 0 and (mother == 2 or father == 2):
        return True
    if child == 1 and (mother + father == 4 or mother + father == 0):
        return True
    if child == 2 and (mother == 0 or father == 0):
        return True
    return False


def get_vars(
    file_path: Path,
    patient_id: str,
    trio: Trio,
    collect_denovo: bool = True,
    collect_inherited: bool = True,
    use_ext_denovo: bool = False,
    hist: Stage2HistCollector | None = None,
    record_hist_denovo: bool = False,
    record_hist_inherited: bool = False,
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]], bool]:
    dvars_inh: Dict[str, List[str]] = {}
    dvars: Dict[str, List[str]] = {}
    open_fn = gzip.open if file_path.suffix == ".gz" else open

    with open_fn(file_path, "rt", encoding="utf-8", errors="replace") as f:
        mother_i = father_i = child_i = -1
        header_ok = False
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHR"):
                indices = trio_column_indices(line.strip().split("\t"), trio)
                if indices is not None:
                    mother_i, father_i, child_i = indices
                    header_ok = True
                continue
            if not header_ok:
                continue

            mother, father, child, kk = get_triplet_at(line, mother_i, father_i, child_i)
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
                is_denovo = if_denovo_ext(chl, mot, fat) if use_ext_denovo else if_denovo(chl, mot, fat)
                if is_denovo:
                    if collect_denovo:
                        update(dvars, kk, gtex)
                        if hist is not None and record_hist_denovo:
                            hist.record_parent(mt, mt_gq, mt_dp, mt_adr, mt_ada)
                            hist.record_parent(ft, ft_gq, ft_dp, ft_adr, ft_ada)
                            hist.record_child(ch, ch_gq, ch_dp, ch_adr, ch_ada)
                else:
                    if collect_inherited:
                        update(dvars_inh, kk, gtex)
                        if hist is not None and record_hist_inherited:
                            hist.record_parent(mt, mt_gq, mt_dp, mt_adr, mt_ada)
                            hist.record_parent(ft, ft_gq, ft_dp, ft_adr, ft_ada)
                            hist.record_child(ch, ch_gq, ch_dp, ch_adr, ch_ada)
            except Exception:
                continue
    return dvars, dvars_inh, header_ok


def chunked(items: Iterable[str], size: int) -> Iterator[List[str]]:
    iterator = iter(items)
    while True:
        chunk = list(islice(iterator, size))
        if not chunk:
            return
        yield chunk


def load_stage1_class_map(path: Path) -> Dict[str, List[str]]:
    return load_class_map(explicit=path)


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


def count_patient_variant_records(dvars: Dict[str, List[str]]) -> int:
    return sum(len(values) for values in dvars.values())


def build_patient_counts_json(
    patient_ids: List[str],
    counts: Dict[str, int],
) -> Dict[str, int]:
    return {patient_id: counts.get(patient_id, 0) for patient_id in patient_ids}


@dataclass(frozen=True)
class BatchOutput:
    class_name: str
    batch_index: int
    patient_count: int
    dvars_path: Path | None
    dvars_inh_path: Path | None


@dataclass(frozen=True)
class Stage2ClassSummary:
    class_name: str
    mean_variants_per_patient: float
    stdev_variants_per_patient: float
    missing_trio: int
    missing_vcf_samples: int


@dataclass(frozen=True)
class Stage2Result:
    output_dir: Path
    parameters_json: Path
    outputs: List[BatchOutput]
    class_summaries: List[Stage2ClassSummary]


STAGE2_PARAMETERS_FILENAME = "stage2_parameters.json"
STAGE2_PARAMETERS_VERSION = 1


@dataclass(frozen=True)
class Stage2Parameters:
    classes: List[str]
    suffix_per_class: Dict[str, str]
    collect: str
    input_dir: Path | None

    @classmethod
    def from_dict(cls, data: Dict[str, object]) -> Stage2Parameters:
        raw_collect = data.get("collect") or data.get("task") or "denovo"
        collect = str(raw_collect)
        if collect not in ("denovo", "inherited"):
            raise ValueError(f"Invalid collect in stage2 parameters: {collect!r}")
        return cls(
            classes=[str(c) for c in data.get("classes", [])],
            suffix_per_class={
                str(k): str(v) for k, v in dict(data.get("suffix_per_class", {})).items()
            },
            collect=collect,
            input_dir=Path(str(data["input_dir"])) if data.get("input_dir") else None,
        )


def stage2_parameters_path(pipeline_root: Path) -> Path:
    return pipeline_root / "stage2" / STAGE2_PARAMETERS_FILENAME


def try_load_stage2_parameters(pipeline_root: Path) -> Stage2Parameters | None:
    path = stage2_parameters_path(pipeline_root)
    if not path.is_file():
        return None
    return Stage2Parameters.from_dict(json.loads(path.read_text(encoding="utf-8")))


def resolve_collect_flags(collect: str, save_all: bool) -> tuple[bool, bool]:
    if collect not in ("denovo", "inherited"):
        raise ValueError("collect must be 'denovo' or 'inherited'")
    collect_denovo = collect == "denovo" or save_all
    collect_inherited = collect == "inherited" or save_all
    return collect_denovo, collect_inherited


def build_stage2_parameters_payload(
    *,
    input_dir: Path,
    output_dir: Path,
    stage1_parameters_path: Path,
    collect: str,
    save_all: bool,
    batch_size: int,
    classes: List[str],
    use_ext_denovo: bool,
    write_hist: bool,
    collect_denovo: bool,
    collect_inherited: bool,
    suffix_per_class: Dict[str, str],
    family_file: Path,
    family_columns_path: Path | None,
    family_columns: FamilyColumns,
    qmap: bool,
    trio_count: int,
) -> Dict[str, object]:
    return {
        "version": STAGE2_PARAMETERS_VERSION,
        "input_dir": str(input_dir.resolve()),
        "output_dir": str(output_dir.resolve()),
        "stage1_parameters": stage1_parameters_path.name,
        "stage1_parameters_path": str(stage1_parameters_path.resolve()),
        "collect": collect,
        "save_all": save_all,
        "batch_size": batch_size,
        "classes": classes,
        "use_ext_denovo": use_ext_denovo,
        "write_hist": write_hist,
        "collect_denovo": collect_denovo,
        "collect_inherited": collect_inherited,
        "suffix_per_class": dict(sorted(suffix_per_class.items())),
        "family_file": str(family_file.resolve()),
        "family_columns_path": str(family_columns_path.resolve()) if family_columns_path else None,
        "family_columns": {
            "child": family_columns.child,
            "mother": family_columns.mother,
            "father": family_columns.father,
            "map": family_columns.map,
        },
        "qmap": qmap,
        "trio_count": trio_count,
    }


def print_stage2_summary(result: Stage2Result, *, collect: str) -> None:
    print(f"[stage2 summary] collect={collect} output={result.output_dir}")
    for row in result.class_summaries:
        print(
            f"  {row.class_name}: mean_variants_per_patient={row.mean_variants_per_patient:.4f} "
            f"stdev={row.stdev_variants_per_patient:.4f} "
            f"missing_trio={row.missing_trio} missing_vcf_samples={row.missing_vcf_samples}"
        )


def _variants_per_patient_stats(counts: List[int]) -> tuple[float, float]:
    if not counts:
        return 0.0, 0.0
    mean = statistics.mean(counts)
    stdev = statistics.stdev(counts) if len(counts) > 1 else 0.0
    return mean, stdev


def build_file_path(input_dir: Path, patient_id: str, suffix: str) -> Path:
    return input_dir / patient_id / f"{patient_id}.{suffix.lstrip('.')}"


def run_stage2(
    input_dir: Path,
    output_dir: Path,
    class_map: Dict[str, List[str]],
    suffix_per_class: Dict[str, str],
    classes: List[str],
    family_file: Path,
    family_columns_spec: str,
    batch_size: int = 1000,
    collect: str = "denovo",
    save_all: bool = False,
    use_ext_denovo: bool = False,
    write_hist: bool = True,
    qmap: bool = False,
    stage1_parameters_path: Path | None = None,
) -> Stage2Result:
    collect_denovo, collect_inherited = resolve_collect_flags(collect, save_all)
    family_columns, family_columns_path = load_family_columns(family_columns_spec, qmap=qmap)
    trios = load_trios(family_file, family_columns, qmap=qmap)

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: List[BatchOutput] = []
    class_summaries: List[Stage2ClassSummary] = []

    missing = [c for c in classes if c not in class_map]
    if missing:
        raise ValueError(
            f"Unknown class(es) for stage2: {missing}. Known classes: {sorted(class_map.keys())}"
        )
    class_items = [(class_name, class_map[class_name]) for class_name in classes]

    for class_name, patient_ids in class_items:
        suffix = suffix_per_class.get(class_name)
        if not suffix:
            raise ValueError(f"No suffix provided for class '{class_name}'")

        class_dir = output_dir / f"class_{class_name}"
        class_dir.mkdir(parents=True, exist_ok=True)

        # Sidecar lists: collect=denovo -> patients with non-empty inherited; collect=inherited -> non-empty denovo.
        patients_nonempty_sidecar: set[str] = set()
        class_dvars_counts: Dict[str, int] = {}
        class_dvars_inh_counts: Dict[str, int] = {}
        hist_collector = Stage2HistCollector() if write_hist else None
        record_hist_denovo = write_hist and collect == "denovo"
        record_hist_inherited = write_hist and collect == "inherited"
        missing_trio: List[str] = []
        missing_vcf_samples: List[str] = []

        for batch_index, batch_patient_ids in enumerate(chunked(patient_ids, batch_size), start=1):
            dVars: Dict[str, List[str]] = {}
            dVars_inh: Dict[str, List[str]] = {}
            batch_patients_nonempty_inh: set[str] = set()
            batch_patients_nonempty_denovo: set[str] = set()

            for patient_id in batch_patient_ids:
                trio = trios.get(patient_id)
                if trio is None:
                    missing_trio.append(patient_id)
                    continue
                file_path = build_file_path(input_dir, patient_id, suffix)
                if not file_path.is_file():
                    continue
                dvars, dvars_inh, header_ok = get_vars(
                    file_path,
                    patient_id,
                    trio,
                    collect_denovo=collect_denovo,
                    collect_inherited=collect_inherited,
                    use_ext_denovo=use_ext_denovo,
                    hist=hist_collector,
                    record_hist_denovo=record_hist_denovo,
                    record_hist_inherited=record_hist_inherited,
                )
                if not header_ok:
                    missing_vcf_samples.append(patient_id)
                    continue
                if collect_denovo:
                    merge(dVars, dvars)
                    class_dvars_counts[patient_id] = (
                        class_dvars_counts.get(patient_id, 0) + count_patient_variant_records(dvars)
                    )
                if collect_inherited:
                    merge(dVars_inh, dvars_inh)
                    class_dvars_inh_counts[patient_id] = (
                        class_dvars_inh_counts.get(patient_id, 0) + count_patient_variant_records(dvars_inh)
                    )
                if collect_inherited and dvars_inh:
                    batch_patients_nonempty_inh.add(patient_id)
                if collect_denovo and dvars:
                    batch_patients_nonempty_denovo.add(patient_id)

            if collect == "denovo":
                batch_sidecar = batch_patients_nonempty_inh
                batch_list_name = f"batch_{batch_index:05d}_patients_nonempty_dvars_inh.txt"
                class_list_name = "patients_nonempty_dvars_inh.txt"
            else:
                batch_sidecar = batch_patients_nonempty_denovo
                batch_list_name = f"batch_{batch_index:05d}_patients_nonempty_dvars.txt"
                class_list_name = "patients_nonempty_dvars.txt"

            patients_nonempty_sidecar.update(batch_sidecar)

            dvars_path: Path | None = None
            dvars_inh_path: Path | None = None

            if collect_denovo:
                dvars_path = class_dir / f"batch_{batch_index:05d}_dVars.pkl"
                with dvars_path.open("wb") as f:
                    pickle.dump(dVars, f)
            if collect_inherited:
                dvars_inh_path = class_dir / f"batch_{batch_index:05d}_dVars_inh.pkl"
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

        if collect_denovo:
            dvars_counts_path = class_dir / "patient_dVars_counts.json"
            dvars_counts_path.write_text(
                json.dumps(
                    build_patient_counts_json(patient_ids, class_dvars_counts),
                    indent=2,
                    sort_keys=True,
                ),
                encoding="utf-8",
            )
        if collect_inherited:
            dvars_inh_counts_path = class_dir / "patient_dVars_inh_counts.json"
            dvars_inh_counts_path.write_text(
                json.dumps(
                    build_patient_counts_json(patient_ids, class_dvars_inh_counts),
                    indent=2,
                    sort_keys=True,
                ),
                encoding="utf-8",
            )

        if hist_collector is not None:
            hist_collector.write_json_files(class_dir, collect=collect)

        if missing_trio:
            (class_dir / "patients_missing_trio.txt").write_text(
                "\n".join(sorted(missing_trio)) + "\n",
                encoding="utf-8",
            )
        if missing_vcf_samples:
            (class_dir / "patients_missing_vcf_samples.txt").write_text(
                "\n".join(sorted(missing_vcf_samples)) + "\n",
                encoding="utf-8",
            )

        count_map = class_dvars_counts if collect == "denovo" else class_dvars_inh_counts
        per_patient = [count_map.get(patient_id, 0) for patient_id in patient_ids]
        mean_v, stdev_v = _variants_per_patient_stats(per_patient)
        class_summaries.append(
            Stage2ClassSummary(
                class_name=class_name,
                mean_variants_per_patient=mean_v,
                stdev_variants_per_patient=stdev_v,
                missing_trio=len(missing_trio),
                missing_vcf_samples=len(missing_vcf_samples),
            )
        )

    if stage1_parameters_path is None:
        raise ValueError("stage1_parameters_path is required for stage2 parameter recording")

    parameters_json = output_dir / STAGE2_PARAMETERS_FILENAME
    parameters_json.write_text(
        json.dumps(
            build_stage2_parameters_payload(
                input_dir=input_dir,
                output_dir=output_dir,
                stage1_parameters_path=stage1_parameters_path,
                collect=collect,
                save_all=save_all,
                batch_size=batch_size,
                classes=classes,
                use_ext_denovo=use_ext_denovo,
                write_hist=write_hist,
                collect_denovo=collect_denovo,
                collect_inherited=collect_inherited,
                suffix_per_class={name: suffix_per_class[name] for name in classes},
                family_file=family_file,
                family_columns_path=family_columns_path,
                family_columns=family_columns,
                qmap=qmap,
                trio_count=len(trios),
            ),
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    return Stage2Result(
        output_dir=output_dir,
        parameters_json=parameters_json,
        outputs=outputs,
        class_summaries=class_summaries,
    )
