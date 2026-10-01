#!/usr/bin/env python3
"""Validate schemas, records, references, and Partial Cubes semantics."""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any


EXPECTED_TABLES = {
    "classes",
    "characterizations",
    "characterization_equivalences",
    "relations",
    "operations",
    "operation_results",
    "invariants",
    "invariant_values",
    "obstruction_families",
    "obstruction_bases",
    "examples",
    "results",
}
SUPPORTED_TYPES = {"string", "integer", "boolean", "text", "latex", "enum", "array", "reference"}
SHORT_NAME = re.compile(r"^[a-z][a-z0-9]*(?:_[a-z0-9]+)*$")
REFERENCE = re.compile(r"^#([a-z][a-z0-9_]*)/([a-z0-9_]+)$")

REFERENCE_FIELDS = {
    ("classes", "primary_characterization_id"): "characterizations",
    ("characterizations", "class_id"): "classes",
    ("characterizations", "invariant_ids"): "invariants",
    ("characterizations", "obstruction_basis_id"): "obstruction_bases",
    ("characterization_equivalences", "class_id"): "classes",
    ("characterization_equivalences", "characterization_ids"): "characterizations",
    ("relations", "subject_id"): "classes",
    ("relations", "object_id"): "classes",
    ("relations", "witness_ids"): "examples",
    ("operation_results", "class_id"): "classes",
    ("operation_results", "operation_id"): "operations",
    ("operation_results", "target_class_id"): "classes",
    ("operation_results", "witness_ids"): "examples",
    ("invariants", "domain_class_id"): "classes",
    ("invariant_values", "invariant_id"): "invariants",
    ("invariant_values", "class_id"): "classes",
    ("invariant_values", "example_id"): "examples",
    ("obstruction_families", "ambient_class_id"): "classes",
    ("obstruction_families", "catalogued_class_id"): "classes",
    ("obstruction_bases", "target_class_id"): "classes",
    ("obstruction_bases", "ambient_class_id"): "classes",
    ("obstruction_bases", "minor_operation_id"): "operations",
    ("obstruction_bases", "obstruction_family_ids"): "obstruction_families",
    ("examples", "obstruction_family_ids"): "obstruction_families",
    ("results", "class_ids"): "classes",
    ("results", "invariant_ids"): "invariants",
    ("results", "operation_ids"): "operations",
    ("results", "witness_ids"): "examples",
}

CLAIM_TABLES = {
    "characterization_equivalences",
    "relations",
    "operation_results",
    "invariant_values",
    "obstruction_bases",
    "results",
}


class Validator:
    def __init__(self, data_dir: Path):
        self.data_dir = data_dir
        self.schemas: dict[str, dict[str, Any]] = {}
        self.records: dict[str, list[dict[str, Any]]] = defaultdict(list)
        self.files: dict[tuple[str, int], Path] = {}
        self.by_id: dict[str, dict[int, dict[str, Any]]] = defaultdict(dict)
        self.by_short: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
        self.errors: list[str] = []

    def error(self, message: str) -> None:
        self.errors.append(message)

    def read_json(self, path: Path) -> Any:
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            self.error(f"{path}: cannot parse JSON: {exc}")
            return None

    def load(self) -> None:
        actual = {
            path.name
            for path in self.data_dir.iterdir()
            if path.is_dir() and (path / "schema.json").exists()
        }
        missing = EXPECTED_TABLES - actual
        extra = actual - EXPECTED_TABLES
        if missing:
            self.error(f"missing schema tables: {', '.join(sorted(missing))}")
        if extra:
            self.error(f"unexpected schema tables: {', '.join(sorted(extra))}")

        for table in sorted(actual):
            schema_path = self.data_dir / table / "schema.json"
            schema = self.read_json(schema_path)
            if not isinstance(schema, dict):
                continue
            self.schemas[table] = schema
            self.validate_schema(table, schema, schema_path)
            for path in sorted((self.data_dir / table).glob("*.json")):
                if path.name == "schema.json":
                    continue
                record = self.read_json(path)
                if isinstance(record, dict):
                    self.load_record(table, record, path)
        self.validate_reverse_reference_columns()

    def validate_reverse_reference_columns(self) -> None:
        for source_table, schema in self.schemas.items():
            for column in schema.get("columns", []):
                if column.get("type") != "reference":
                    continue
                target_table = column.get("table")
                target_column = column.get("column")
                target_schema = self.schemas.get(target_table, {})
                target_names = {item.get("name") for item in target_schema.get("columns", [])}
                if target_column not in target_names:
                    self.error(
                        f"{self.data_dir / source_table / 'schema.json'}: "
                        f"reverse reference {column.get('name')} targets missing "
                        f"{target_table}.{target_column}"
                    )

    def validate_schema(self, table: str, schema: dict[str, Any], path: Path) -> None:
        if schema.get("table_name") != table:
            self.error(f"{path}: table_name must be {table!r}")
        columns = schema.get("columns")
        if not isinstance(columns, list):
            self.error(f"{path}: columns must be an array")
            return
        names = []
        for index, column in enumerate(columns):
            if not isinstance(column, dict):
                self.error(f"{path}: column {index} must be an object")
                continue
            name = column.get("name")
            if not isinstance(name, str):
                self.error(f"{path}: column {index} has no string name")
                continue
            names.append(name)
            if column.get("type") not in SUPPORTED_TYPES:
                self.error(f"{path}: {name} has unsupported type {column.get('type')!r}")
            if column.get("type") == "enum":
                choices = column.get("enum")
                if not isinstance(choices, list) or not choices:
                    self.error(f"{path}: enum column {name} has no choices")
                elif len({item.get("value") for item in choices if isinstance(item, dict)}) != len(choices):
                    self.error(f"{path}: enum column {name} has duplicate or malformed choices")
            if column.get("type") == "reference":
                if column.get("table") not in EXPECTED_TABLES or not isinstance(column.get("column"), str):
                    self.error(f"{path}: reverse reference {name} has invalid table/column")
        if len(names) != len(set(names)):
            self.error(f"{path}: duplicate column names")
        for required in ("id", "short_name", "name"):
            if required not in names:
                self.error(f"{path}: missing core column {required}")

    def load_record(self, table: str, record: dict[str, Any], path: Path) -> None:
        identifier = record.get("id")
        short_name = record.get("short_name")
        if not isinstance(identifier, int) or isinstance(identifier, bool):
            self.error(f"{path}: id must be an integer")
            return
        if identifier <= 0:
            self.error(f"{path}: id must be positive")
        if identifier in self.by_id[table]:
            self.error(f"{path}: duplicate {table} id {identifier}")
        if not isinstance(short_name, str) or not SHORT_NAME.fullmatch(short_name):
            self.error(f"{path}: short_name must be lower snake case")
            return
        if short_name in self.by_short[table]:
            self.error(f"{path}: duplicate {table} short_name {short_name!r}")
        expected_name = f"{identifier:03d}_{short_name}.json"
        if path.name != expected_name:
            self.error(f"{path}: filename must be {expected_name}")
        self.records[table].append(record)
        self.by_id[table][identifier] = record
        self.by_short[table][short_name] = record
        self.files[(table, identifier)] = path

    def validate_records(self) -> None:
        for table, records in self.records.items():
            schema = self.schemas.get(table)
            if not schema:
                continue
            columns = {column["name"]: column for column in schema["columns"]}
            stored_columns = {
                name: column for name, column in columns.items() if column.get("type") != "reference"
            }
            for record in records:
                path = self.files[(table, record["id"])]
                unknown = set(record) - set(stored_columns)
                if unknown:
                    self.error(f"{path}: unknown or computed fields: {', '.join(sorted(unknown))}")
                for name, column in stored_columns.items():
                    value = record.get(name)
                    if column.get("required") and (value is None or value == "" or value == []):
                        self.error(f"{path}: required field {name} is missing or empty")
                        continue
                    if value is None:
                        continue
                    if not self.value_has_type(value, column["type"]):
                        self.error(f"{path}: field {name} must have type {column['type']}")
                    if column["type"] == "enum" and isinstance(value, str):
                        choices = {
                            item["value"]
                            for item in column["enum"]
                            if isinstance(item, dict) and isinstance(item.get("value"), str)
                        }
                        if value not in choices:
                            self.error(f"{path}: field {name} has invalid enum value {value!r}")
                self.validate_references(table, record, path)
                self.validate_provenance(table, record, path)

    @staticmethod
    def value_has_type(value: Any, kind: str) -> bool:
        if kind == "integer":
            return isinstance(value, int) and not isinstance(value, bool)
        if kind == "boolean":
            return isinstance(value, bool)
        if kind == "array":
            return isinstance(value, list)
        return isinstance(value, str)

    def resolve(self, reference: Any, expected_table: str | None = None) -> tuple[str | None, dict[str, Any] | None]:
        if not isinstance(reference, str):
            return None, None
        match = REFERENCE.fullmatch(reference)
        if not match:
            return None, None
        table, key = match.groups()
        if expected_table and table != expected_table:
            return table, None
        if key.isdigit():
            return table, self.by_id[table].get(int(key))
        return table, self.by_short[table].get(key)

    def validate_references(self, table: str, record: dict[str, Any], path: Path) -> None:
        for (source_table, field), target_table in REFERENCE_FIELDS.items():
            if source_table != table or field not in record or record[field] is None:
                continue
            values = record[field] if isinstance(record[field], list) else [record[field]]
            for value in values:
                parsed_table, target = self.resolve(value, target_table)
                if parsed_table is None:
                    self.error(f"{path}: {field} contains malformed reference {value!r}")
                elif parsed_table != target_table:
                    self.error(f"{path}: {field} must reference #{target_table}/...")
                elif target is None:
                    self.error(f"{path}: {field} references missing record {value!r}")

    def validate_provenance(self, table: str, record: dict[str, Any], path: Path) -> None:
        if record.get("status") != "established":
            return
        needs_provenance = table in CLAIM_TABLES
        if needs_provenance and not any(record.get(field) for field in ("proof", "references", "proof_source")):
            self.error(f"{path}: established claim needs proof, references, or proof_source")

    def validate_semantics(self) -> None:
        self.validate_class_characterizations()
        self.validate_equivalences()
        self.validate_relations()
        self.validate_invariant_values()
        self.validate_obstructions()
        self.validate_results()
        self.validate_closure_summaries()
        self.validate_citations()

    def canonical_ref(self, table: str, record: dict[str, Any]) -> str:
        return f'#{table}/{record["short_name"]}'

    def validate_class_characterizations(self) -> None:
        established_graph: dict[str, set[str]] = defaultdict(set)
        for equivalence in self.records["characterization_equivalences"]:
            if equivalence.get("status") != "established":
                continue
            resolved = [self.resolve(ref, "characterizations")[1] for ref in equivalence.get("characterization_ids", [])]
            refs = [self.canonical_ref("characterizations", item) for item in resolved if item]
            for first in refs:
                established_graph[first].update(second for second in refs if second != first)

        for class_record in self.records["classes"]:
            path = self.files[("classes", class_record["id"])]
            _, primary = self.resolve(class_record.get("primary_characterization_id"), "characterizations")
            if not primary:
                continue
            if self.resolve(primary.get("class_id"), "classes")[1] is not class_record:
                self.error(f"{path}: primary characterization points back to another class")
            primary_ref = self.canonical_ref("characterizations", primary)
            reachable = {primary_ref}
            pending = [primary_ref]
            while pending:
                current = pending.pop()
                for neighbor in established_graph[current] - reachable:
                    reachable.add(neighbor)
                    pending.append(neighbor)
            for characterization in self.records["characterizations"]:
                if characterization.get("status") != "established":
                    continue
                if self.resolve(characterization.get("class_id"), "classes")[1] is not class_record:
                    continue
                reference = self.canonical_ref("characterizations", characterization)
                if reference != primary_ref and reference not in reachable:
                    char_path = self.files[("characterizations", characterization["id"])]
                    self.error(f"{char_path}: established alternative is not connected to the primary characterization")

        for characterization in self.records["characterizations"]:
            path = self.files[("characterizations", characterization["id"])]
            kind = characterization.get("characterization_type")
            if kind == "invariant_condition" and not characterization.get("invariant_ids"):
                self.error(f"{path}: invariant_condition needs invariant_ids")
            if kind == "forbidden_minor" and not characterization.get("obstruction_basis_id"):
                self.error(f"{path}: forbidden_minor needs obstruction_basis_id")

    def validate_equivalences(self) -> None:
        for record in self.records["characterization_equivalences"]:
            path = self.files[("characterization_equivalences", record["id"])]
            references = record.get("characterization_ids", [])
            if len(references) < 2 or len(set(references)) < 2:
                self.error(f"{path}: characterization_ids needs at least two distinct references")
            _, declared_class = self.resolve(record.get("class_id"), "classes")
            for reference in references:
                _, characterization = self.resolve(reference, "characterizations")
                if characterization and self.resolve(characterization.get("class_id"), "classes")[1] is not declared_class:
                    self.error(f"{path}: every characterization must belong to the declared class")

    def validate_relations(self) -> None:
        adjacency: dict[str, set[str]] = defaultdict(set)
        relation_paths = {}
        for record in self.records["relations"]:
            path = self.files[("relations", record["id"])]
            _, subject = self.resolve(record.get("subject_id"), "classes")
            _, object_ = self.resolve(record.get("object_id"), "classes")
            if subject and subject is object_:
                self.error(f"{path}: a class relation cannot have identical endpoints")
            if (
                subject and object_
                and record.get("status") == "established"
                and record.get("relation_type") in {"inclusion", "strict_inclusion"}
            ):
                source = self.canonical_ref("classes", subject)
                target = self.canonical_ref("classes", object_)
                adjacency[source].add(target)
                relation_paths[(source, target)] = path

        state: dict[str, int] = {}
        def visit(node: str) -> None:
            state[node] = 1
            for successor in adjacency[node]:
                if state.get(successor) == 1:
                    self.error(f"{relation_paths[(node, successor)]}: established inclusions contain a directed cycle")
                elif state.get(successor, 0) == 0:
                    visit(successor)
            state[node] = 2
        for node in list(adjacency):
            if state.get(node, 0) == 0:
                visit(node)

    def validate_invariant_values(self) -> None:
        for record in self.records["invariant_values"]:
            path = self.files[("invariant_values", record["id"])]
            subjects = int(bool(record.get("class_id"))) + int(bool(record.get("example_id")))
            if subjects != 1:
                self.error(f"{path}: exactly one of class_id and example_id must be set")
            if record.get("example_id") and record.get("quantifier") != "individual":
                self.error(f"{path}: example values must use quantifier individual")

    def validate_obstructions(self) -> None:
        for record in self.records["obstruction_bases"]:
            path = self.files[("obstruction_bases", record["id"])]
            if len(set(record.get("obstruction_family_ids", []))) != len(record.get("obstruction_family_ids", [])):
                self.error(f"{path}: obstruction_family_ids contains duplicates")
            _, operation = self.resolve(record.get("minor_operation_id"), "operations")
            if operation and operation.get("category") != "minor":
                self.error(f"{path}: minor_operation_id must name an operation in category minor")

        for characterization in self.records["characterizations"]:
            if characterization.get("characterization_type") != "forbidden_minor":
                continue
            path = self.files[("characterizations", characterization["id"])]
            _, basis = self.resolve(characterization.get("obstruction_basis_id"), "obstruction_bases")
            _, target = self.resolve(characterization.get("class_id"), "classes")
            if basis and self.resolve(basis.get("target_class_id"), "classes")[1] is not target:
                self.error(f"{path}: obstruction basis target does not match characterization class")

    def validate_results(self) -> None:
        for record in self.records["results"]:
            path = self.files[("results", record["id"])]
            if not any(record.get(field) for field in ("class_ids", "invariant_ids", "operation_ids")):
                self.error(f"{path}: result needs at least one class, invariant, or operation subject")

    def validate_closure_summaries(self) -> None:
        operation_to_field = {"projection": "p_closed", "conditioning": "c_closed", "pc_minor": "pc_closed"}
        field_to_status = {
            "p_closed": "p_closure_status",
            "c_closed": "c_closure_status",
            "pc_closed": "pc_closure_status",
        }
        authoritative: dict[tuple[int, str], set[bool]] = defaultdict(set)
        for result in self.records["operation_results"]:
            if result.get("status") != "established" or result.get("result_type") not in {"closed_under", "not_closed_under"}:
                continue
            _, class_record = self.resolve(result.get("class_id"), "classes")
            _, operation = self.resolve(result.get("operation_id"), "operations")
            if not class_record or not operation or operation.get("short_name") not in operation_to_field:
                continue
            field = operation_to_field[operation["short_name"]]
            authoritative[(class_record["id"], field)].add(result["result_type"] == "closed_under")

        for class_record in self.records["classes"]:
            path = self.files[("classes", class_record["id"])]
            for field in operation_to_field.values():
                facts = authoritative[(class_record["id"], field)]
                status = class_record.get(field_to_status[field])
                if len(facts) > 1:
                    self.error(f"{path}: contradictory established operation results for {field}")
                elif field in class_record and class_record[field] is not None:
                    if not facts:
                        self.error(f"{path}: cached {field} has no established operation result")
                    elif class_record[field] not in facts:
                        self.error(f"{path}: cached {field} disagrees with its operation result")
                    expected_status = "closed" if class_record[field] else "not_closed"
                    if status != expected_status:
                        self.error(f"{path}: {field_to_status[field]} disagrees with cached {field}")
                elif facts:
                    self.error(f"{path}: established operation result is not cached in {field}")
                elif status in {"closed", "not_closed"}:
                    self.error(f"{path}: established status {status!r} has no operation result or Boolean cache")

            p = class_record.get("p_closed")
            c = class_record.get("c_closed")
            pc = class_record.get("pc_closed")
            if pc is True and (p is False or c is False):
                self.error(f"{path}: pc_closed=true is incompatible with a failed p/c closure")
            if p is True and c is True and pc is not True:
                self.error(f"{path}: p_closed=true and c_closed=true imply pc_closed=true")

    def validate_citations(self) -> None:
        bib_path = self.data_dir / "latex" / "references.bib"
        try:
            bibliography = bib_path.read_text(encoding="utf-8")
        except OSError as exc:
            self.error(f"{bib_path}: cannot read bibliography: {exc}")
            return
        keys = set(re.findall(r"@\w+\s*\{\s*([^,\s]+)", bibliography))
        citation_patterns = (
            re.compile(r"\\cite\w*\{([^}]+)\}"),
            re.compile(r"\[@([A-Za-z0-9_:.+-]+)\]"),
        )
        for table, records in self.records.items():
            for record in records:
                for value in record.values():
                    if not isinstance(value, str):
                        continue
                    cited = set()
                    for match in citation_patterns[0].finditer(value):
                        cited.update(item.strip() for item in match.group(1).split(","))
                    cited.update(match.group(1) for match in citation_patterns[1].finditer(value))
                    missing = cited - keys
                    if missing:
                        path = self.files[(table, record["id"])]
                        self.error(f"{path}: missing bibliography keys: {', '.join(sorted(missing))}")

    def run(self) -> int:
        self.load()
        self.validate_records()
        self.validate_semantics()
        if self.errors:
            for error in self.errors:
                print(f"ERROR: {error}", file=sys.stderr)
            print(f"Validation failed with {len(self.errors)} error(s).", file=sys.stderr)
            return 1
        count = sum(len(records) for records in self.records.values())
        print(f"Validated {len(self.schemas)} schemas and {count} records.")
        return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    args = parser.parse_args()
    return Validator(args.data_dir.resolve()).run()


if __name__ == "__main__":
    raise SystemExit(main())
