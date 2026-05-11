"""Tests for ``src.ssu_bundler``.

Pure transformation tests on small fake per-MAG FASTA files plus a fake
routing log. No subprocess, no network. The bundler is exercised in three
modes: full SILVA bundle, partial bundle (one per-MAG FASTA empty), and
upsert-into-existing routing-log behaviour.
"""

from __future__ import annotations

from pathlib import Path

from src.ssu_bundler import (
    RoutingDecision,
    bundle_silva_sina,
    read_routing_log,
    upsert_routing_log,
    write_routing_log,
)
from src.ssu_extractor import (
    ROUTING_LOG_COLUMNS,
    SSU_DECISION_GTDBTK,
    SSU_DECISION_SILVA,
)


def _fake_silva_decision(
    accession: str,
    *,
    ssu_fasta_path: Path,
    n_16S: int = 1,
) -> RoutingDecision:
    """Build a fake SILVA-decision routing row for tests."""
    return RoutingDecision(
        accession=accession,
        source="ENA",
        fasta_path=Path("/tmp/fake.fasta.gz"),
        n_contigs_scanned=10,
        n_16S_predicted=n_16S,
        decision=SSU_DECISION_SILVA,
        output_path=ssu_fasta_path,
        barrnap_version="barrnap 0.9",
        notes="",
    )


def _fake_gtdbtk_decision(accession: str, fasta_path: Path) -> RoutingDecision:
    """Build a fake GTDB-Tk-decision routing row for tests."""
    return RoutingDecision(
        accession=accession,
        source="ENA",
        fasta_path=fasta_path,
        n_contigs_scanned=10,
        n_16S_predicted=0,
        decision=SSU_DECISION_GTDBTK,
        output_path=fasta_path,
        barrnap_version="barrnap 0.9",
        notes="",
    )


def _write_per_mag_ssu_fasta(path: Path, *, accession: str, n_records: int = 1) -> None:
    """Write a small per-MAG 16S FASTA with the requested record count."""
    body_lines = []
    for i in range(n_records):
        body_lines.append(f">{accession}|ctg{i}|1-20|+|len=20")
        body_lines.append("ACGT" * 5)
    path.write_text("\n".join(body_lines) + ("\n" if body_lines else ""), encoding="utf-8")


def test_write_routing_log_uses_canonical_columns(tmp_path: Path):
    """The header row matches :data:`ROUTING_LOG_COLUMNS` exactly."""
    log = tmp_path / "routing_log.tsv"
    n = write_routing_log(
        [
            _fake_silva_decision("MAG_A", ssu_fasta_path=tmp_path / "MAG_A.16S.fasta"),
            _fake_gtdbtk_decision("MAG_B", tmp_path / "MAG_B.fasta"),
        ],
        log,
    )
    assert n == 2
    text = log.read_text(encoding="utf-8")
    header = text.splitlines()[0].split("\t")
    assert tuple(header) == ROUTING_LOG_COLUMNS


def test_upsert_routing_log_replaces_in_place(tmp_path: Path):
    """Re-running upsert with the same accession overwrites that row only."""
    log = tmp_path / "routing_log.tsv"
    write_routing_log(
        [
            _fake_silva_decision("MAG_A", ssu_fasta_path=tmp_path / "MAG_A.16S.fasta"),
            _fake_gtdbtk_decision("MAG_B", tmp_path / "MAG_B.fasta"),
        ],
        log,
    )
    # Replace MAG_A's row with a GTDB-Tk decision.
    upsert_routing_log(
        log,
        [_fake_gtdbtk_decision("MAG_A", tmp_path / "MAG_A.fasta")],
    )
    rows = read_routing_log(log)
    by_acc = {r["accession"]: r for r in rows}
    assert by_acc["MAG_A"]["decision"] == SSU_DECISION_GTDBTK
    assert by_acc["MAG_B"]["decision"] == SSU_DECISION_GTDBTK


def test_bundle_silva_sina_concatenates_only_silva_decisions(tmp_path: Path):
    """The bundle pulls in every SILVA-routed MAG and ignores GTDB-Tk rows."""
    ssu_dir = tmp_path / "ssu"
    ssu_dir.mkdir()
    fasta_a = ssu_dir / "MAG_A.16S.fasta"
    fasta_b = ssu_dir / "MAG_B.16S.fasta"
    _write_per_mag_ssu_fasta(fasta_a, accession="MAG_A", n_records=2)
    _write_per_mag_ssu_fasta(fasta_b, accession="MAG_B", n_records=1)
    log = ssu_dir / "routing_log.tsv"
    write_routing_log(
        [
            _fake_silva_decision("MAG_A", ssu_fasta_path=fasta_a, n_16S=2),
            _fake_silva_decision("MAG_B", ssu_fasta_path=fasta_b, n_16S=1),
            _fake_gtdbtk_decision("MAG_C", ssu_dir / "MAG_C.fasta"),
        ],
        log,
    )

    bundle = ssu_dir / "silva_sina_input.fasta"
    result = bundle_silva_sina(log, ssu_dir, bundle)

    assert result.n_mags_included == 2
    assert result.n_records_total == 3
    assert result.skipped_accessions == []
    body = bundle.read_text(encoding="utf-8")
    assert body.count(">") == 3
    assert "MAG_A" in body
    assert "MAG_B" in body
    assert "MAG_C" not in body


def test_bundle_silva_sina_skips_empty_per_mag_fastas(tmp_path: Path):
    """SILVA-routed MAGs whose 16S FASTA is missing or empty are reported."""
    ssu_dir = tmp_path / "ssu"
    ssu_dir.mkdir()
    fasta_a = ssu_dir / "MAG_A.16S.fasta"
    fasta_b = ssu_dir / "MAG_B.16S.fasta"
    _write_per_mag_ssu_fasta(fasta_a, accession="MAG_A", n_records=1)
    fasta_b.write_text("", encoding="utf-8")
    log = ssu_dir / "routing_log.tsv"
    write_routing_log(
        [
            _fake_silva_decision("MAG_A", ssu_fasta_path=fasta_a),
            _fake_silva_decision("MAG_B", ssu_fasta_path=fasta_b),
        ],
        log,
    )
    bundle = ssu_dir / "silva_sina_input.fasta"
    result = bundle_silva_sina(log, ssu_dir, bundle)
    assert result.n_mags_included == 1
    assert result.skipped_accessions == ["MAG_B"]


def test_read_routing_log_returns_empty_when_missing(tmp_path: Path):
    """A missing log returns an empty list rather than raising."""
    assert read_routing_log(tmp_path / "nonexistent.tsv") == []
