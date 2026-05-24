#!/usr/bin/env python3
"""
Metagenome reads → SQLite database
====================================
Converts the simulated metagenome BLAST nucleotide database into a
queryable SQLite database, joining sequence data with the read-level
metadata from read_labels.tsv.

Why SQLite:
  - Universally renderable (DB Browser for SQLite, pandas, R DBI, etc.)
  - Supports full SQL queries on read metadata (condition, genome, is_pathogen)
  - Compresses well for Zenodo archival
  - No BLAST installation required to inspect data

Extraction method:
  blastdbcmd -outfmt "%t\t%s" recovers original FASTA headers via the
  sequence title field (%t), which preserves the InSilicoSeq read IDs
  (format: <contig>_<readN>_<pairN>/<direction>). These IDs match
  read_labels.tsv exactly — confirmed by auditing 10 records before build.

Output:
  data/results/blast_screen/metagenome_reads.db

Tables:
  reads    — read_id TEXT PK, sequence TEXT, seq_len INTEGER
  metadata — read_id TEXT PK, condition TEXT, genome TEXT, is_pathogen INTEGER
  summary  — aggregate stats per (condition, genome)

Run inside the fp_pipeline conda environment:
    conda run -n fp_pipeline python scripts/analysis/build_metagenome_sqlite.py
"""

import logging
import sqlite3
import subprocess
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

REPO_ROOT   = Path(__file__).parents[2]
BLAST_DIR   = REPO_ROOT / "data" / "results" / "blast_screen"
META_DIR    = REPO_ROOT / "data" / "metagenome"
DB_PATH     = BLAST_DIR / "metagenome_reads.db"
BLAST_DB    = BLAST_DIR / "metagenome_db"
LABELS_TSV  = META_DIR / "read_labels.tsv"

BIN_DIR     = Path(sys.executable).parent
BLASTDBCMD  = BIN_DIR / "blastdbcmd"
if not BLASTDBCMD.exists():
    BLASTDBCMD = Path("blastdbcmd")  # fall back to PATH


def audit_id_fidelity(n_check: int = 20) -> None:
    """
    Verify that blastdbcmd %t output IDs match read_labels.tsv IDs.
    Raises RuntimeError if any mismatch is detected.
    """
    log.info("Auditing record ID fidelity (first %d records)…", n_check)

    # Dump all titles then slice — blastdbcmd has no -max_seqs flag
    dump = subprocess.run(
        [str(BLASTDBCMD), "-db", str(BLAST_DB), "-entry", "all", "-outfmt", "%t"],
        capture_output=True, text=True, check=True,
    )
    blast_ids = [ln.strip() for ln in dump.stdout.strip().splitlines()
                 if ln.strip()][:n_check]

    label_ids = []
    with open(LABELS_TSV) as fh:
        next(fh)  # skip header
        for line in fh:
            label_ids.append(line.split("\t")[0].strip())
            if len(label_ids) == n_check:
                break

    mismatches = [(b, l) for b, l in zip(blast_ids, label_ids) if b != l]
    if mismatches:
        raise RuntimeError(
            f"Record ID mismatch between BLAST DB and read_labels.tsv "
            f"({len(mismatches)}/{n_check} records differ). "
            f"First mismatch: BLAST={mismatches[0][0]!r}  TSV={mismatches[0][1]!r}"
        )
    log.info("ID fidelity confirmed — %d/%d records matched exactly.", n_check, n_check)


def create_schema(conn: sqlite3.Connection) -> None:
    conn.executescript("""
        DROP TABLE IF EXISTS reads;
        DROP TABLE IF EXISTS metadata;
        DROP TABLE IF EXISTS summary;

        CREATE TABLE reads (
            read_id  TEXT PRIMARY KEY,
            sequence TEXT NOT NULL,
            seq_len  INTEGER NOT NULL
        );

        CREATE TABLE metadata (
            read_id    TEXT PRIMARY KEY,
            condition  TEXT,   -- 'low' | 'mid' | 'high' (O157:H7 spike level)
            genome     TEXT,   -- source genome name (e.g. O157H7_Sakai)
            is_pathogen INTEGER  -- 1 = pathogenic origin, 0 = commensal
        );

        CREATE TABLE summary (
            condition   TEXT,
            genome      TEXT,
            is_pathogen INTEGER,
            n_reads     INTEGER,
            mean_len    REAL,
            PRIMARY KEY (condition, genome)
        );
    """)
    conn.commit()
    log.info("Schema created: tables reads, metadata, summary")


def load_reads(conn: sqlite3.Connection) -> int:
    """Stream sequences from BLAST DB via blastdbcmd and insert into reads table."""
    log.info("Dumping sequences from BLAST DB (this may take a minute)…")

    proc = subprocess.Popen(
        [str(BLASTDBCMD), "-db", str(BLAST_DB), "-entry", "all",
         "-outfmt", "%t\t%s"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )

    conn.execute("BEGIN")
    n = 0
    chunk = []
    CHUNK_SIZE = 50_000

    for line in proc.stdout:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        read_id, seq = parts[0].strip(), parts[1].strip()
        chunk.append((read_id, seq, len(seq)))
        n += 1
        if len(chunk) >= CHUNK_SIZE:
            conn.executemany("INSERT OR IGNORE INTO reads VALUES (?,?,?)", chunk)
            chunk.clear()
            log.info("  Inserted %d reads so far…", n)

    if chunk:
        conn.executemany("INSERT OR IGNORE INTO reads VALUES (?,?,?)", chunk)

    proc.wait()
    if proc.returncode != 0:
        err = proc.stderr.read()
        raise RuntimeError(f"blastdbcmd failed: {err}")

    conn.execute("COMMIT")
    log.info("Reads loaded: %d total", n)
    return n


def load_metadata(conn: sqlite3.Connection) -> None:
    """Stream read_labels.tsv and insert into metadata table."""
    log.info("Loading read metadata from %s…", LABELS_TSV)
    conn.execute("BEGIN")
    chunk = []
    CHUNK_SIZE = 50_000
    n = 0

    with open(LABELS_TSV) as fh:
        next(fh)  # skip header: read_id condition genome is_pathogen
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            read_id, condition, genome, is_pathogen = parts[:4]
            chunk.append((read_id.strip(), condition.strip(),
                          genome.strip(), int(is_pathogen)))
            n += 1
            if len(chunk) >= CHUNK_SIZE:
                conn.executemany("INSERT OR IGNORE INTO metadata VALUES (?,?,?,?)", chunk)
                chunk.clear()

    if chunk:
        conn.executemany("INSERT OR IGNORE INTO metadata VALUES (?,?,?,?)", chunk)

    conn.execute("COMMIT")
    log.info("Metadata loaded: %d records", n)


def build_summary(conn: sqlite3.Connection) -> None:
    """Materialise per-(condition, genome) aggregate stats into the summary table."""
    log.info("Building summary table…")
    conn.executescript("""
        INSERT INTO summary (condition, genome, is_pathogen, n_reads, mean_len)
        SELECT
            m.condition,
            m.genome,
            m.is_pathogen,
            COUNT(*)            AS n_reads,
            AVG(r.seq_len)      AS mean_len
        FROM metadata m
        JOIN reads r ON m.read_id = r.read_id
        GROUP BY m.condition, m.genome;
    """)
    conn.commit()

    rows = conn.execute(
        "SELECT condition, genome, is_pathogen, n_reads, ROUND(mean_len,1) "
        "FROM summary ORDER BY condition, is_pathogen DESC"
    ).fetchall()
    log.info("Summary table (%d rows):", len(rows))
    for row in rows:
        log.info("  condition=%-6s  genome=%-20s  pathogen=%d  n=%7d  mean_len=%.1f",
                 *row)


def build_indexes(conn: sqlite3.Connection) -> None:
    conn.executescript("""
        CREATE INDEX IF NOT EXISTS idx_reads_len     ON reads    (seq_len);
        CREATE INDEX IF NOT EXISTS idx_meta_genome   ON metadata (genome);
        CREATE INDEX IF NOT EXISTS idx_meta_pathogen ON metadata (is_pathogen);
        CREATE INDEX IF NOT EXISTS idx_meta_cond     ON metadata (condition);
    """)
    conn.commit()
    log.info("Indexes built.")


def main() -> None:
    log.info("=== Metagenome BLAST DB → SQLite  |  fp_pipeline env ===")
    log.info("Source BLAST DB : %s", BLAST_DB)
    log.info("Source labels   : %s", LABELS_TSV)
    log.info("Output SQLite   : %s", DB_PATH)

    if not BLAST_DB.with_suffix(".ndb").exists():
        raise FileNotFoundError(f"BLAST DB not found: {BLAST_DB}.*")
    if not LABELS_TSV.exists():
        raise FileNotFoundError(f"read_labels.tsv not found: {LABELS_TSV}")

    # Safety check: confirm IDs match before writing a byte to the DB
    audit_id_fidelity(n_check=20)

    if DB_PATH.exists():
        log.warning("Removing existing SQLite DB: %s", DB_PATH)
        DB_PATH.unlink()

    conn = sqlite3.connect(str(DB_PATH))
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")

    try:
        create_schema(conn)
        n_reads = load_reads(conn)
        load_metadata(conn)
        build_summary(conn)
        build_indexes(conn)

        # Final validation
        n_joined = conn.execute(
            "SELECT COUNT(*) FROM reads r JOIN metadata m ON r.read_id = m.read_id"
        ).fetchone()[0]
        log.info("Join validation: %d / %d reads have metadata (%.1f%%)",
                 n_joined, n_reads, 100 * n_joined / max(n_reads, 1))

    finally:
        conn.close()

    size_mb = DB_PATH.stat().st_size / 1e6
    log.info("SQLite DB written: %.1f MB → %s", size_mb, DB_PATH)
    log.info("=== Build complete — DB is renderable in DB Browser for SQLite, pandas, R DBI ===")


if __name__ == "__main__":
    main()
