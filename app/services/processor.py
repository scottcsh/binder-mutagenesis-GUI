import csv
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

from Bio.PDB import PDBIO, PDBParser, Select, is_aa


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PISA_SCRIPT_PATH = PROJECT_ROOT / "scripts" / "pisa_batch_dir_to_csv.py"
SUPPORTED_AA = set("ACDEFGHIKLMNPQRSTVWY")
AA1_TO_AA3 = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}
AA3_TO_AA1 = {value: key for key, value in AA1_TO_AA3.items()}


class MultiChainProteinSelect(Select):
    def __init__(self, chain_ids: Sequence[str]) -> None:
        super().__init__()
        self.chain_ids = set(chain_ids)

    def accept_chain(self, chain) -> bool:
        return chain.id in self.chain_ids

    def accept_residue(self, residue) -> bool:
        return residue.id[0] == " " and is_aa(residue, standard=False)


def process_mutagenesis(
    input_pdb: str,
    output_dir: str,
    target_chain: str,
    design_chain: str,
    residue_range: str,
    replicates: int,
    settings: Dict[str, str],
    log_callback: Optional[Callable[[str], None]] = None,
    input_result_callback: Optional[Callable[[Dict[str, object]], None]] = None,
    result_callback: Optional[Callable[[Dict[str, object]], None]] = None,
    mutagenesis_result_callback: Optional[Callable[[List[Dict[str, object]]], None]] = None,
    progress_callback: Optional[Callable[[int], None]] = None,
    should_cancel_callback: Optional[Callable[[], bool]] = None,
    register_process_callback: Optional[Callable[[subprocess.Popen], None]] = None,
    unregister_process_callback: Optional[Callable[[subprocess.Popen], None]] = None,
) -> List[Dict[str, object]]:
    input_path = Path(input_pdb).expanduser().resolve()
    out_dir = Path(output_dir).expanduser().resolve()
    target_chain = target_chain.strip()
    design_chain = design_chain.strip()
    replicates = int(replicates)

    if not input_path.exists() or input_path.suffix.lower() != ".pdb":
        raise ValueError("Input PDB must be an existing .pdb file.")
    if not target_chain:
        raise ValueError("Target chain is required.")
    if not design_chain:
        raise ValueError("Design chain is required.")
    if target_chain == design_chain:
        raise ValueError("Target chain and design chain must be different.")
    if replicates <= 0:
        raise ValueError("ProteinMPNN replicates must be a positive integer.")

    out_dir.mkdir(parents=True, exist_ok=True)
    run_dir = out_dir / f"{input_path.stem}_mutagenesis"
    run_dir.mkdir(parents=True, exist_ok=True)
    variants_dir = run_dir / "variants"
    variants_dir.mkdir(parents=True, exist_ok=True)

    _emit(progress_callback, 0)
    _log(log_callback, f"[INFO] Input PDB: {input_path}")
    _log(log_callback, f"[INFO] Run directory: {run_dir}")

    residues_to_screen = parse_residue_range(residue_range)
    structure_data = load_structure_data(input_path, target_chain, design_chain)
    screen_residues = select_existing_design_residues(
        design_residue_numbers=structure_data["design_residues"],
        requested_residue_numbers=residues_to_screen,
    )
    if not screen_residues:
        raise ValueError("No design-chain residues were found in the requested range.")

    _log(log_callback, "[INFO] Screen residues: " + ", ".join(str(r) for r in screen_residues))
    _emit(progress_callback, 3)

    native_sequence = structure_data["design_sequence"]
    input_result = {
        "input pdb": input_path.name,
        "residue range": format_residue_list(screen_residues),
        "sequence": sequence_for_residues(native_sequence, structure_data["design_residues"], screen_residues),
        "EvoEF2 binding": "",
        "Interface area": "",
        "Interface solvation E": "",
        "p-value": "",
    }
    _emit_input_result(input_result_callback, input_result)

    _log(log_callback, "[INFO] Running EvoEF2 ComputeBinding for input PDB")
    input_result["EvoEF2 binding"] = run_evoef2_compute_binding(
        variant_path=input_path,
        settings=settings,
        log_callback=log_callback,
        should_cancel_callback=should_cancel_callback,
        register_process_callback=register_process_callback,
        unregister_process_callback=unregister_process_callback,
    )
    _emit_input_result(input_result_callback, input_result)
    _emit(progress_callback, 8)

    _log(log_callback, "[INFO] Running PISA for input PDB")
    input_pisa_scores = run_pisa_batch(
        variant_paths=[input_path],
        run_dir=run_dir,
        target_chain=target_chain,
        design_chain=design_chain,
        settings=settings,
        log_callback=log_callback,
        should_cancel_callback=should_cancel_callback,
        register_process_callback=register_process_callback,
        unregister_process_callback=unregister_process_callback,
        batch_name="input_pisa",
        progress_callback=progress_callback,
        progress_start=8,
        progress_end=15,
    )
    input_pisa = input_pisa_scores.get(input_path.name, {})
    input_result["Interface area"] = input_pisa.get("Interface area", "")
    input_result["Interface solvation E"] = input_pisa.get("Interface solvation E", "")
    input_result["p-value"] = input_pisa.get("p-value", "")
    _emit_input_result(input_result_callback, input_result)
    _emit(progress_callback, 15)

    all_rows: List[Dict[str, object]] = []
    _raise_if_cancelled(should_cancel_callback)
    mpnn_input_dir = run_dir / "input"
    mpnn_input_dir.mkdir(parents=True, exist_ok=True)
    mpnn_input_pdb = mpnn_input_dir / input_path.name
    shutil.copy2(input_path, mpnn_input_pdb)
    _emit(progress_callback, 18)

    fixed_positions = build_fixed_positions_for_mutable_residues(
        design_residue_numbers=structure_data["design_residues"],
        mutable_residues=screen_residues,
    )
    _log(log_callback, "[INFO] Running ProteinMPNN for requested residue range")
    _emit(progress_callback, 20)
    mpnn_output_dir = run_proteinmpnn(
        input_dir=mpnn_input_dir,
        output_dir=run_dir / "proteinmpnn",
        design_chain=design_chain,
        replicates=replicates,
        fixed_positions=fixed_positions,
        settings=settings,
        log_callback=log_callback,
        should_cancel_callback=should_cancel_callback,
        register_process_callback=register_process_callback,
        unregister_process_callback=unregister_process_callback,
    )

    designed_sequences = read_mpnn_sequences(
        mpnn_output_dir,
        input_path.stem,
        expected_design_length=len(structure_data["design_residues"]),
    )
    if not designed_sequences:
        raise RuntimeError("ProteinMPNN did not produce sequences for the requested residue range.")
    _emit(progress_callback, 35)

    selected_sequences = designed_sequences[:replicates]
    for replicate_index, sequence in enumerate(selected_sequences, start=1):
        _raise_if_cancelled(should_cancel_callback)
        residue_to_new_aa = {
            residue_number: sequence_aa_at_residue(sequence, structure_data["design_residues"], residue_number)
            for residue_number in screen_residues
        }
        designed_range_sequence = sequence_for_residues(sequence, structure_data["design_residues"], screen_residues)
        variant_name = f"{input_path.stem}_mpnn_rep{replicate_index}.pdb"
        variant_path = variants_dir / variant_name
        write_multi_point_variant_pdb(input_path, variant_path, design_chain, residue_to_new_aa)
        row = {
            "input pdb": input_path.name,
            "variant pdb": variant_path.name,
            "variant pdb path": str(variant_path),
            "sequence": designed_range_sequence,
            "EvoEF2 binding": "",
            "Interface area": "",
            "Interface solvation E": "",
            "p-value": "",
            "status": "variant ready",
        }
        all_rows.append(row)
        _emit_result(result_callback, row)
        _emit(progress_callback, 35 + int((replicate_index / max(len(selected_sequences), 1)) * 10))

    _emit(progress_callback, 45)

    _log(log_callback, f"[INFO] Running EvoEF2 ComputeBinding for {len(all_rows)} variant(s)")
    for index, row in enumerate(all_rows, start=1):
        _raise_if_cancelled(should_cancel_callback)
        variant_path = Path(str(row["variant pdb path"]))
        score = run_evoef2_compute_binding(
            variant_path=variant_path,
            settings=settings,
            log_callback=log_callback,
            should_cancel_callback=should_cancel_callback,
            register_process_callback=register_process_callback,
            unregister_process_callback=unregister_process_callback,
        )
        row["EvoEF2 binding"] = score if score is not None else ""
        row["status"] = "EvoEF2 done"
        _emit_result(result_callback, row)
        _emit(progress_callback, 45 + int((index / max(len(all_rows), 1)) * 25))

    _log(log_callback, f"[INFO] Running PISA for {len(all_rows)} variant(s)")
    pisa_scores = run_pisa_batch(
        variant_paths=[Path(str(row["variant pdb path"])) for row in all_rows],
        run_dir=run_dir,
        target_chain=target_chain,
        design_chain=design_chain,
        settings=settings,
        log_callback=log_callback,
        should_cancel_callback=should_cancel_callback,
        register_process_callback=register_process_callback,
        unregister_process_callback=unregister_process_callback,
        batch_name="variant_pisa",
        progress_callback=progress_callback,
        progress_start=70,
        progress_end=95,
    )
    for row in all_rows:
        result = pisa_scores.get(str(row["variant pdb"]), {})
        row["Interface area"] = result.get("Interface area", "")
        row["Interface solvation E"] = result.get("Interface solvation E", "")
        row["p-value"] = result.get("p-value", "")
        row["status"] = "completed"
        _emit_result(result_callback, row)

    mutagenesis_rows = summarize_improved_mutations(
        rows=all_rows,
        input_result=input_result,
        native_sequence=native_sequence,
        residue_numbers=structure_data["design_residues"],
        screen_residues=screen_residues,
        design_chain=design_chain,
    )
    _emit_mutagenesis_result(mutagenesis_result_callback, mutagenesis_rows)
    _emit(progress_callback, 98)

    output_csv = run_dir / "mutagenesis_results.csv"
    write_results_csv(output_csv, all_rows)
    mutagenesis_csv = run_dir / "mutation_counts.csv"
    write_mutation_counts_csv(mutagenesis_csv, mutagenesis_rows)
    _log(log_callback, f"[INFO] Wrote results: {output_csv}")
    _log(log_callback, f"[INFO] Wrote mutation counts: {mutagenesis_csv}")
    _emit(progress_callback, 100)
    return all_rows


def parse_residue_range(text: str) -> List[int]:
    values = set()
    for token in re.split(r"[,\s]+", text.strip()):
        if not token:
            continue
        match = re.fullmatch(r"(-?\d+)\s*-\s*(-?\d+)", token)
        if match:
            start = int(match.group(1))
            end = int(match.group(2))
            step = 1 if end >= start else -1
            values.update(range(start, end + step, step))
            continue
        if re.fullmatch(r"-?\d+", token):
            values.add(int(token))
            continue
        raise ValueError(f"Invalid residue range token: {token}")
    if not values:
        raise ValueError("Residue range is required.")
    return sorted(values)


def load_structure_data(input_path: Path, target_chain: str, design_chain: str) -> Dict[str, object]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("input", str(input_path))
    model = next(structure.get_models())
    if target_chain not in model:
        raise ValueError(f"Target chain was not found: {target_chain}")
    if design_chain not in model:
        raise ValueError(f"Design chain was not found: {design_chain}")

    design_residues = []
    design_sequence = []
    for residue in model[design_chain].get_residues():
        if residue.id[0] != " " or not is_aa(residue, standard=False):
            continue
        resseq = int(residue.id[1])
        aa1 = AA3_TO_AA1.get(str(residue.resname).upper(), "X")
        design_residues.append(resseq)
        design_sequence.append(aa1)

    if not design_residues:
        raise ValueError(f"No amino-acid residues were found in design chain {design_chain}.")

    return {
        "structure": structure,
        "model": model,
        "design_residues": design_residues,
        "design_sequence": "".join(design_sequence),
    }


def select_existing_design_residues(
    design_residue_numbers: Sequence[int],
    requested_residue_numbers: Sequence[int],
) -> List[int]:
    requested = set(int(number) for number in requested_residue_numbers)
    return [int(number) for number in design_residue_numbers if int(number) in requested]


def build_fixed_positions_for_mutable_residues(
    design_residue_numbers: Sequence[int],
    mutable_residues: Sequence[int],
) -> str:
    mutable = set(int(number) for number in mutable_residues)
    fixed = [str(number) for number in design_residue_numbers if int(number) not in mutable]
    return " ".join(fixed)


def format_residue_list(residue_numbers: Sequence[int]) -> str:
    return ",".join(str(number) for number in residue_numbers)


def sequence_for_residues(sequence: str, residue_numbers: Sequence[int], screen_residues: Sequence[int]) -> str:
    index_by_residue = {int(residue_number): index for index, residue_number in enumerate(residue_numbers)}
    selected = []
    for residue_number in screen_residues:
        index = index_by_residue[int(residue_number)]
        selected.append(sequence[index])
    return "".join(selected)


def summarize_improved_mutations(
    rows: Sequence[Dict[str, object]],
    input_result: Dict[str, object],
    native_sequence: str,
    residue_numbers: Sequence[int],
    screen_residues: Sequence[int],
    design_chain: str,
) -> List[Dict[str, object]]:
    counts: Dict[str, Dict[str, object]] = {}
    input_evoef2 = _optional_float(input_result.get("EvoEF2 binding", ""))
    input_area = _optional_float(input_result.get("Interface area", ""))
    input_solvation = _optional_float(input_result.get("Interface solvation E", ""))

    for row in rows:
        evoef2_improved = is_lower_better_improved(row.get("EvoEF2 binding", ""), input_evoef2)
        area_improved = is_higher_better_improved(row.get("Interface area", ""), input_area)
        solvation_improved = is_lower_better_improved(row.get("Interface solvation E", ""), input_solvation)
        if not (evoef2_improved or area_improved or solvation_improved):
            continue

        designed_sequence = str(row.get("sequence", ""))
        for index, residue_number in enumerate(screen_residues):
            if index >= len(designed_sequence):
                continue
            old_aa = sequence_aa_at_residue(native_sequence, residue_numbers, residue_number)
            new_aa = designed_sequence[index].upper()
            if old_aa == new_aa or new_aa not in SUPPORTED_AA:
                continue
            mutation = f"{old_aa}{residue_number}{new_aa}"
            entry = counts.setdefault(
                mutation,
                {
                    "mutation": mutation,
                    "count": 0,
                    "EvoEF2 improved": 0,
                    "Interface area improved": 0,
                    "Solvation E improved": 0,
                },
            )
            entry["count"] = int(entry["count"]) + 1
            if evoef2_improved:
                entry["EvoEF2 improved"] = int(entry["EvoEF2 improved"]) + 1
            if area_improved:
                entry["Interface area improved"] = int(entry["Interface area improved"]) + 1
            if solvation_improved:
                entry["Solvation E improved"] = int(entry["Solvation E improved"]) + 1

    return sorted(
        counts.values(),
        key=lambda item: (
            -int(item.get("count", 0)),
            str(item.get("mutation", "")),
        ),
    )


def is_lower_better_improved(value: object, baseline: object) -> bool:
    candidate = _optional_float(value)
    if not isinstance(candidate, float) or not isinstance(baseline, float):
        return False
    return candidate < baseline


def is_higher_better_improved(value: object, baseline: object) -> bool:
    candidate = _optional_float(value)
    if not isinstance(candidate, float) or not isinstance(baseline, float):
        return False
    return candidate > baseline


def run_proteinmpnn(
    input_dir: Path,
    output_dir: Path,
    design_chain: str,
    replicates: int,
    fixed_positions: str,
    settings: Dict[str, str],
    log_callback: Optional[Callable[[str], None]],
    should_cancel_callback: Optional[Callable[[], bool]],
    register_process_callback: Optional[Callable[[subprocess.Popen], None]],
    unregister_process_callback: Optional[Callable[[subprocess.Popen], None]],
) -> Path:
    proteinmpnn_root = Path(settings.get("proteinmpnn_dir", "")).expanduser()
    if not proteinmpnn_root.exists():
        raise ValueError("ProteinMPNN directory is not configured or does not exist.")

    output_dir.mkdir(parents=True, exist_ok=True)
    parse_script = proteinmpnn_root / "helper_scripts" / "parse_multiple_chains.py"
    assign_script = proteinmpnn_root / "helper_scripts" / "assign_fixed_chains.py"
    fixed_script = proteinmpnn_root / "helper_scripts" / "make_fixed_positions_dict.py"
    mpnn_script = proteinmpnn_root / "protein_mpnn_run.py"
    for script in [parse_script, assign_script, fixed_script, mpnn_script]:
        if not script.exists():
            raise FileNotFoundError(f"ProteinMPNN script was not found: {script}")

    run_script = output_dir / "run_proteinmpnn.sh"
    env_setup = settings.get("proteinmpnn_env_setup", "").strip()
    lines = ["#!/bin/bash", "set -euo pipefail", ""]
    if env_setup:
        lines.extend(env_setup.splitlines())
        lines.append("")
    lines.extend(
        [
            f'folder_with_pdbs="{input_dir}"',
            f'output_dir="{output_dir}"',
            'parsed_jsonl="$output_dir/parsed_pdbs.jsonl"',
            'assigned_jsonl="$output_dir/assigned_pdbs.jsonl"',
            'fixed_jsonl="$output_dir/fixed_pdbs.jsonl"',
            f'python -u "{parse_script}" --input_path="$folder_with_pdbs" --output_path="$parsed_jsonl"',
            f'python -u "{assign_script}" --input_path="$parsed_jsonl" --output_path="$assigned_jsonl" --chain_list "{design_chain}"',
            f'python -u "{fixed_script}" --input_path="$parsed_jsonl" --output_path="$fixed_jsonl" --chain_list "{design_chain}" --position_list "{fixed_positions}"',
            f'python -u "{mpnn_script}" --jsonl_path "$parsed_jsonl" --chain_id_jsonl "$assigned_jsonl" --fixed_positions_jsonl "$fixed_jsonl" --out_folder "$output_dir" --num_seq_per_target {replicates} --sampling_temp 0.4 --seed {settings.get("proteinmpnn_seed", "37") or "37"} --batch_size {settings.get("proteinmpnn_batch_size", "1") or "1"}',
            'echo "ProteinMPNN run completed."',
        ]
    )
    run_script.write_text("\n".join(lines) + "\n", encoding="utf-8")
    run_script.chmod(0o755)

    command = ["bash", str(run_script)]
    _run_streaming_process(
        command=command,
        cwd=output_dir,
        log_callback=log_callback,
        should_cancel_callback=should_cancel_callback,
        register_process_callback=register_process_callback,
        unregister_process_callback=unregister_process_callback,
    )
    return output_dir


def read_mpnn_sequences(mpnn_output_dir: Path, pdb_stem: str, expected_design_length: int) -> List[str]:
    seq_dir = mpnn_output_dir / "seqs"
    candidates = [seq_dir / f"{pdb_stem}.fa", seq_dir / f"{pdb_stem}.fasta"]
    fasta_path = next((path for path in candidates if path.exists()), None)
    if fasta_path is None:
        fasta_files = sorted(seq_dir.glob("*.fa")) + sorted(seq_dir.glob("*.fasta"))
        fasta_path = fasta_files[0] if fasta_files else None
    if fasta_path is None:
        return []
    records = []
    header = None
    seq_lines = []
    for line in fasta_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_lines).strip()))
            header = line
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq_lines).strip()))
    sequences = []
    for _, sequence in records[1:]:
        sequence = sequence.strip()
        if not sequence:
            continue
        if "/" in sequence:
            segments = [segment.strip() for segment in sequence.split("/") if segment.strip()]
            matching = [segment for segment in segments if len(segment) == expected_design_length]
            sequence = matching[-1] if matching else max(segments, key=len)
        sequences.append(sequence)
    return sequences


def sequence_aa_at_residue(sequence: str, residue_numbers: Sequence[int], residue_number: int) -> str:
    index = list(residue_numbers).index(int(residue_number))
    if index >= len(sequence):
        raise ValueError(f"ProteinMPNN sequence is too short for residue {residue_number}.")
    aa = sequence[index].upper()
    if aa not in SUPPORTED_AA:
        raise ValueError(f"Unsupported amino acid from ProteinMPNN at residue {residue_number}: {aa}")
    return aa


def write_multi_point_variant_pdb(
    input_path: Path,
    output_path: Path,
    chain_id: str,
    residue_to_new_aa: Dict[int, str],
) -> None:
    residue_to_resname = {int(residue_number): AA1_TO_AA3[new_aa] for residue_number, new_aa in residue_to_new_aa.items()}
    changed = False
    with input_path.open("r", encoding="utf-8", errors="replace") as source, output_path.open("w", encoding="utf-8", newline="") as target:
        for line in source:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 26:
                line_chain = line[21].strip()
                try:
                    line_resseq = int(line[22:26])
                except ValueError:
                    line_resseq = None
                if line_chain == chain_id and line_resseq in residue_to_resname:
                    line = line[:17] + residue_to_resname[line_resseq].rjust(3) + line[20:]
                    changed = True
            target.write(line)
    if not changed:
        raise ValueError(f"No requested residues were found in chain {chain_id} in the PDB text.")


def prepare_selected_chains_pdb(source_path: Path, output_pdb_path: Path, chain_ids: Sequence[str]) -> None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(source_path.stem, str(source_path))
    model = next(structure.get_models())
    normalized_chain_ids = []
    for chain_id in chain_ids:
        if chain_id and chain_id not in normalized_chain_ids:
            normalized_chain_ids.append(chain_id)
    if not normalized_chain_ids:
        raise ValueError(f"No chains were selected for {source_path.name}")

    available_chain_ids = {chain.id for chain in model}
    missing = [chain_id for chain_id in normalized_chain_ids if chain_id not in available_chain_ids]
    if missing:
        raise ValueError(
            "Selected chain(s) {0} were not found in {1}. Available chains: {2}".format(
                ", ".join(missing),
                source_path.name,
                ", ".join(sorted(available_chain_ids)) or "(none)",
            )
        )

    for chain_id in normalized_chain_ids:
        chain = model[chain_id]
        has_protein_residue = any(residue.id[0] == " " and is_aa(residue, standard=False) for residue in chain)
        if not has_protein_residue:
            raise ValueError(f"Chain '{chain_id}' in {source_path.name} has no protein residues.")

    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb_path), MultiChainProteinSelect(normalized_chain_ids))


def run_evoef2_compute_binding(
    variant_path: Path,
    settings: Dict[str, str],
    log_callback: Optional[Callable[[str], None]],
    should_cancel_callback: Optional[Callable[[], bool]],
    register_process_callback: Optional[Callable[[subprocess.Popen], None]],
    unregister_process_callback: Optional[Callable[[subprocess.Popen], None]],
) -> Optional[str]:
    executable = settings.get("evoef2_executable", "").strip() or "EvoEF2"
    command = [executable, "--command=ComputeBinding", f"--pdb={variant_path}"]
    output = _run_capture_process(
        command=command,
        cwd=variant_path.parent,
        env_updates={"OMP_NUM_THREADS": settings.get("evoef2_threads", "1") or "1"},
        log_callback=log_callback,
        should_cancel_callback=should_cancel_callback,
        register_process_callback=register_process_callback,
        unregister_process_callback=unregister_process_callback,
    )
    return parse_evoef2_total_score(output)


def parse_evoef2_total_score(text: str) -> Optional[str]:
    match = re.search(r"Total\s*=\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)", text)
    if not match:
        return None
    return match.group(1)


def run_pisa_batch(
    variant_paths: Sequence[Path],
    run_dir: Path,
    target_chain: str,
    design_chain: str,
    settings: Dict[str, str],
    log_callback: Optional[Callable[[str], None]],
    should_cancel_callback: Optional[Callable[[], bool]],
    register_process_callback: Optional[Callable[[subprocess.Popen], None]],
    unregister_process_callback: Optional[Callable[[subprocess.Popen], None]],
    batch_name: str = "pisa",
    progress_callback: Optional[Callable[[int], None]] = None,
    progress_start: int = 0,
    progress_end: int = 100,
) -> Dict[str, Dict[str, object]]:
    if not PISA_SCRIPT_PATH.exists():
        raise FileNotFoundError(f"PISA script was not found: {PISA_SCRIPT_PATH}")
    pisa_run_dir = run_dir / batch_name
    pisa_input_dir = pisa_run_dir / "pisa_input"
    pisa_output_dir = pisa_run_dir / "pisa_output"
    pisa_input_dir.mkdir(parents=True, exist_ok=True)
    pisa_output_dir.mkdir(parents=True, exist_ok=True)
    for existing in pisa_input_dir.glob("*.pdb"):
        existing.unlink()
    total_files = max(len(variant_paths), 1)
    for variant_path in variant_paths:
        prepared_pdb = pisa_input_dir / variant_path.name
        prepare_selected_chains_pdb(variant_path, prepared_pdb, [target_chain, design_chain])
        _log(log_callback, f"[INFO] Prepared PISA input PDB: {prepared_pdb}")

    command = [
        sys.executable,
        str(PISA_SCRIPT_PATH),
        "--input",
        str(pisa_input_dir),
        "--output-dir",
        str(pisa_output_dir),
        "--concurrency",
        "5",
    ]
    _log(log_callback, "[INFO] Running PISA batch for {0} file(s): {1}".format(len(variant_paths), " ".join(command)))
    completed_files = 0

    def on_pisa_line(line: str) -> None:
        nonlocal completed_files
        if line.startswith("[DONE]") or line.startswith("[ERROR] Failed for "):
            completed_files += 1
            progress_value = progress_start + int((completed_files / total_files) * (progress_end - progress_start))
            _emit(progress_callback, min(progress_end, progress_value))

    _run_streaming_process(
        command=command,
        cwd=pisa_run_dir,
        log_callback=log_callback,
        should_cancel_callback=should_cancel_callback,
        register_process_callback=register_process_callback,
        unregister_process_callback=unregister_process_callback,
        line_callback=on_pisa_line,
    )
    _emit(progress_callback, progress_end)
    return parse_pisa_scores(pisa_output_dir / "PISA_interface.csv")


def parse_pisa_scores(csv_path: Path) -> Dict[str, Dict[str, object]]:
    if not csv_path.exists():
        raise FileNotFoundError(f"PISA output CSV was not found: {csv_path}")
    scores: Dict[str, Dict[str, object]] = {}
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            file_name = Path(str(row.get("pdb_file", "")).strip()).name
            if not file_name:
                continue
            candidate = {
                "Interface area": _optional_float(row.get("int_area", "")),
                "Interface solvation E": _optional_float(row.get("int_solv_energy", "")),
                "p-value": _optional_float(row.get("pvalue", "")),
            }
            current = scores.get(file_name)
            if current is None:
                scores[file_name] = candidate
                continue
            candidate_energy = candidate.get("Interface solvation E", "")
            current_energy = current.get("Interface solvation E", "")
            if isinstance(candidate_energy, float) and (
                not isinstance(current_energy, float) or candidate_energy < current_energy
            ):
                scores[file_name] = candidate
    return scores


def write_results_csv(output_csv: Path, rows: Sequence[Dict[str, object]]) -> None:
    columns = [
        "input pdb",
        "variant pdb",
        "sequence",
        "EvoEF2 binding",
        "Interface area",
        "Interface solvation E",
        "p-value",
        "status",
        "variant pdb path",
    ]
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_mutation_counts_csv(output_csv: Path, rows: Sequence[Dict[str, object]]) -> None:
    columns = [
        "mutation",
        "count",
        "EvoEF2 improved",
        "Interface area improved",
        "Solvation E improved",
    ]
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def _run_streaming_process(
    command: Sequence[str],
    cwd: Path,
    log_callback: Optional[Callable[[str], None]],
    should_cancel_callback: Optional[Callable[[], bool]],
    register_process_callback: Optional[Callable[[subprocess.Popen], None]],
    unregister_process_callback: Optional[Callable[[subprocess.Popen], None]],
    line_callback: Optional[Callable[[str], None]] = None,
) -> None:
    _log(log_callback, "[INFO] Running command: " + " ".join(str(part) for part in command))
    process = subprocess.Popen(
        [str(part) for part in command],
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        start_new_session=True,
    )
    _register_process(register_process_callback, process)
    try:
        if process.stdout is not None:
            for raw_line in process.stdout:
                if _is_cancelled(should_cancel_callback):
                    _terminate_process(process)
                    break
                line = raw_line.rstrip()
                if line:
                    _log(log_callback, line)
                    if line_callback is not None:
                        line_callback(line)
        return_code = process.wait()
        if return_code != 0:
            if _is_cancelled(should_cancel_callback):
                raise RuntimeError("Cancelled")
            raise RuntimeError(f"Command failed with exit code {return_code}: {' '.join(command)}")
    finally:
        _unregister_process(unregister_process_callback, process)


def _run_capture_process(
    command: Sequence[str],
    cwd: Path,
    env_updates: Dict[str, str],
    log_callback: Optional[Callable[[str], None]],
    should_cancel_callback: Optional[Callable[[], bool]],
    register_process_callback: Optional[Callable[[subprocess.Popen], None]],
    unregister_process_callback: Optional[Callable[[subprocess.Popen], None]],
) -> str:
    _log(log_callback, "[INFO] Running command: " + " ".join(str(part) for part in command))
    env = dict(os.environ)
    env.update(env_updates)
    process = subprocess.Popen(
        [str(part) for part in command],
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        env=env,
        start_new_session=True,
    )
    _register_process(register_process_callback, process)
    output_lines = []
    try:
        if process.stdout is not None:
            for raw_line in process.stdout:
                if _is_cancelled(should_cancel_callback):
                    _terminate_process(process)
                    break
                line = raw_line.rstrip()
                output_lines.append(line)
                if line:
                    _log(log_callback, line)
        return_code = process.wait()
        if return_code != 0:
            if _is_cancelled(should_cancel_callback):
                raise RuntimeError("Cancelled")
            raise RuntimeError(f"Command failed with exit code {return_code}: {' '.join(command)}")
        return "\n".join(output_lines)
    finally:
        _unregister_process(unregister_process_callback, process)


def _optional_float(value: object) -> object:
    text = str(value).strip()
    if not text:
        return ""
    try:
        return float(text)
    except ValueError:
        return ""


def _emit(callback: Optional[Callable[[int], None]], value: int) -> None:
    if callback is not None:
        callback(value)


def _emit_result(callback: Optional[Callable[[Dict[str, object]], None]], row: Dict[str, object]) -> None:
    if callback is not None:
        callback(dict(row))


def _emit_input_result(callback: Optional[Callable[[Dict[str, object]], None]], row: Dict[str, object]) -> None:
    if callback is not None:
        callback(dict(row))


def _emit_mutagenesis_result(
    callback: Optional[Callable[[List[Dict[str, object]]], None]],
    rows: List[Dict[str, object]],
) -> None:
    if callback is not None:
        callback([dict(row) for row in rows])


def _log(callback: Optional[Callable[[str], None]], message: str) -> None:
    if callback is not None:
        callback(message)


def _is_cancelled(callback: Optional[Callable[[], bool]]) -> bool:
    return bool(callback is not None and callback())


def _raise_if_cancelled(callback: Optional[Callable[[], bool]]) -> None:
    if _is_cancelled(callback):
        raise RuntimeError("Cancelled")


def _register_process(callback: Optional[Callable[[subprocess.Popen], None]], process: subprocess.Popen) -> None:
    if callback is not None:
        callback(process)


def _unregister_process(callback: Optional[Callable[[subprocess.Popen], None]], process: subprocess.Popen) -> None:
    if callback is not None:
        callback(process)


def _terminate_process(process: subprocess.Popen) -> None:
    if process.poll() is not None:
        return
    try:
        process.terminate()
    except Exception:
        pass
