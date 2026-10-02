#!/usr/bin/env python3
"""Independently verify the passive-carrier nine-concept certificates."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path


CAMPAIGN = Path(__file__).resolve().parent
DEFAULT_MANIFEST = CAMPAIGN / "damp_nine_passive_carrier_cnf_manifest.json"
DEFAULT_OUTPUT = CAMPAIGN / "damp_nine_passive_carrier_proof_verification.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def dimacs_header(path: Path) -> tuple[int, int]:
    with path.open(encoding="ascii") as stream:
        for line in stream:
            if line.startswith("p cnf "):
                _, _, variables, clauses = line.split()
                return int(variables), int(clauses)
    raise ValueError(f"no DIMACS header in {path}")


def run_json_script(path: Path) -> dict[str, object]:
    completed = subprocess.run(
        [sys.executable, str(path)],
        cwd=CAMPAIGN,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"local certificate failed for {path.name}:\n"
            f"{completed.stdout}{completed.stderr}"
        )
    return json.loads(completed.stdout)


def parse_checker_output(output: str) -> dict[str, int | float | str]:
    normalized = output.replace("\r", "\n")
    original = re.search(r"(\d+) of (\d+) clauses in core", normalized)
    lemmas = re.search(
        r"(\d+) of (\d+) lemmas in core using (\d+) resolution steps",
        normalized,
    )
    rat = re.search(
        r"(\d+) RAT lemmas in core; (\d+) redundant literals", normalized
    )
    elapsed = re.search(r"verification time: ([0-9.]+) seconds", normalized)
    if "s VERIFIED" not in normalized or not all((original, lemmas, rat, elapsed)):
        raise AssertionError(f"unexpected drat-trim output:\n{normalized}")
    assert original is not None and lemmas is not None
    assert rat is not None and elapsed is not None
    return {
        "status": "VERIFIED",
        "core_original_clauses": int(original.group(1)),
        "original_clauses": int(original.group(2)),
        "core_lemmas": int(lemmas.group(1)),
        "proof_lemmas": int(lemmas.group(2)),
        "resolution_steps": int(lemmas.group(3)),
        "core_rat_lemmas": int(rat.group(1)),
        "core_redundant_literals": int(rat.group(2)),
        "verification_seconds": float(elapsed.group(1)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--drat-trim", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    manifest = json.loads(arguments.manifest.read_text())

    for source in manifest["sources"]:
        path = CAMPAIGN / source["file"]
        if sha256(path) != source["sha256"]:
            raise AssertionError(f"source hash changed: {path.name}")

    local_records: dict[str, object] = {}
    for name, certificate in manifest["local_certificates"].items():
        path = CAMPAIGN / certificate["script"]
        if sha256(path) != certificate["script_sha256"]:
            raise AssertionError(f"local-certificate source changed: {path.name}")
        record = run_json_script(path)
        if record.get("status") != "VERIFIED":
            raise AssertionError(f"local certificate did not verify: {path.name}")
        if record.get("mathematical_digest") != certificate["digest"]:
            raise AssertionError(f"local-certificate digest changed: {path.name}")
        local_records[name] = {
            "script": path.name,
            "script_sha256": sha256(path),
            "mathematical_digest": record["mathematical_digest"],
            "status": record["status"],
        }

    checker = arguments.drat_trim.resolve()
    profile_records = []
    for profile in manifest["profiles"]:
        cnf = CAMPAIGN / profile["cnf"]
        proof = CAMPAIGN / profile["core_proof"]
        generator_record_path = CAMPAIGN / profile["generator_record"]
        generator_record = json.loads(generator_record_path.read_text())
        variables, clauses = dimacs_header(cnf)

        expected = {
            "profile": profile["profile"],
            "status": "unsat",
            "dimacs_sha256": profile["cnf_sha256"],
            "dimacs_clause_count": profile["clauses"],
            "smt2_sha256": profile["smt2_sha256"],
        }
        for key, value in expected.items():
            if generator_record.get(key) != value:
                raise AssertionError(
                    f"generator record changed for {profile['profile']}: {key}"
                )
        if (variables, clauses) != (profile["variables"], profile["clauses"]):
            raise AssertionError(f"DIMACS header changed for {profile['profile']}")
        if sha256(cnf) != profile["cnf_sha256"]:
            raise AssertionError(f"CNF hash changed for {profile['profile']}")
        if proof.stat().st_size != profile["core_proof_bytes"]:
            raise AssertionError(f"proof size changed for {profile['profile']}")
        if sha256(proof) != profile["core_proof_sha256"]:
            raise AssertionError(f"proof hash changed for {profile['profile']}")

        completed = subprocess.run(
            [str(checker), str(cnf), str(proof)],
            cwd=CAMPAIGN,
            check=False,
            capture_output=True,
            text=True,
        )
        transcript = completed.stdout + completed.stderr
        if completed.returncode != 0:
            raise RuntimeError(
                f"DRAT verification failed for {profile['profile']}:\n{transcript}"
            )
        verification = parse_checker_output(transcript)
        if verification["original_clauses"] != clauses:
            raise AssertionError("checker clause count disagrees with DIMACS header")
        profile_records.append(
            {
                "profile": profile["profile"],
                "variables": variables,
                "clauses": clauses,
                "cnf": cnf.name,
                "cnf_sha256": sha256(cnf),
                "core_proof": proof.name,
                "core_proof_bytes": proof.stat().st_size,
                "core_proof_sha256": sha256(proof),
                "generator_record": generator_record_path.name,
                "verification": verification,
            }
        )
        print(f"profile={profile['profile']}: VERIFIED")

    mathematics = {
        "schema": "damp-nine-passive-carrier-proof-verification-v1",
        "claim": manifest["claim"],
        "scope": manifest["scope"],
        "local_certificates": local_records,
        "profiles": profile_records,
        "all_verified": True,
        "consequence": (
            "No deletion-minimal trapped ample interval has nine concepts; "
            "therefore the relative corner threshold is at least ten."
        ),
    }
    stable_mathematics = json.loads(json.dumps(mathematics))
    for profile in stable_mathematics["profiles"]:
        profile["verification"].pop("verification_seconds", None)
    canonical = json.dumps(
        stable_mathematics, sort_keys=True, separators=(",", ":")
    )
    output = {
        "mathematics": mathematics,
        "checker": {
            "name": "drat-trim",
            "path_used": str(checker),
            "binary_sha256": sha256(checker),
        },
        "metadata": {
            "generator": Path(__file__).name,
            "manifest": arguments.manifest.name,
            "manifest_sha256": sha256(arguments.manifest),
            "mathematical_sha256": hashlib.sha256(canonical.encode()).hexdigest(),
        },
    }
    arguments.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"mathematical_sha256={output['metadata']['mathematical_sha256']}")
    print(f"wrote={arguments.output}")


if __name__ == "__main__":
    main()
