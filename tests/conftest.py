"""Shared pytest fixtures for BioResilient AI tests.

Fixture dataset: 3 species (human, naked_mole_rat, ground_squirrel), 5 genes.
Uses an in-memory SQLite database — no PostgreSQL required for unit tests.
"""

import uuid
from pathlib import Path

import pytest
from sqlalchemy import create_engine, event
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

# ---------------------------------------------------------------------------
# Patch ARRAY → JSON before any model imports so SQLite can render the schema.
# This must happen at module level (before db.models is imported anywhere).
# ---------------------------------------------------------------------------
from sqlalchemy.dialects.postgresql import JSON as PG_JSON
import sqlalchemy.dialects.postgresql as _pg_dialect

if not getattr(_pg_dialect, "_array_patched_for_sqlite", False):
    class _ArrayAsJSON(PG_JSON):
        def __init__(self, *args, **kwargs):
            super().__init__()

    _pg_dialect.ARRAY = _ArrayAsJSON
    _pg_dialect._array_patched_for_sqlite = True

# ---------------------------------------------------------------------------
# SQLite in-memory engine (PostgreSQL arrays shimmed as JSON above)
# ---------------------------------------------------------------------------

_TEST_DB_URL = "sqlite://"


@pytest.fixture(scope="session")
def engine():
    """Session-scoped in-memory SQLite engine."""
    from db.models import Base

    eng = create_engine(
        _TEST_DB_URL,
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    Base.metadata.create_all(eng)
    yield eng
    Base.metadata.drop_all(eng)


@pytest.fixture
def session(engine):
    """Function-scoped test session with automatic rollback."""
    connection = engine.connect()
    transaction = connection.begin()
    Session = sessionmaker(bind=connection)
    sess = Session()

    yield sess

    sess.close()
    transaction.rollback()
    connection.close()


# ---------------------------------------------------------------------------
# Fixture data
# ---------------------------------------------------------------------------

FIXTURE_SPECIES = [
    {
        "id": "human",
        "taxid": 9606,
        "name": "Homo sapiens",
        "phenotype": ["baseline"],
        "genome_assembly": "GCF_000001405.40",
        "lineage_group": "Primates",
    },
    {
        "id": "naked_mole_rat",
        "taxid": 10181,
        "name": "Heterocephalus glaber",
        "phenotype": ["cancer_resistance", "longevity"],
        "genome_assembly": "GCF_000247695.1",
        "lineage_group": "Rodents",
    },
    {
        "id": "ground_squirrel",
        "taxid": 43179,
        "name": "Ictidomys tridecemlineatus",
        "phenotype": ["hibernation"],
        "genome_assembly": "GCF_000236235.1",
        "lineage_group": "Rodents",
    },
]

FIXTURE_GENES = [
    {"human_gene_id": "GENE001", "gene_symbol": "TP53",  "human_protein": "P04637"},
    {"human_gene_id": "GENE002", "gene_symbol": "BRCA1", "human_protein": "P38398"},
    {"human_gene_id": "GENE003", "gene_symbol": "PARP1", "human_protein": "P09874"},
    {"human_gene_id": "GENE004", "gene_symbol": "FOXO3", "human_protein": "O43524"},
    {"human_gene_id": "GENE005", "gene_symbol": "SIRT1", "human_protein": "Q96EB6"},
]

FIXTURE_SEQS = {
    "human":          "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDP",
    "naked_mole_rat": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSRAMDDLMLSPDDIEQWFTEDP",
    "ground_squirrel":"MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPAQAMDDLMLSPDDIEQWFTEDP",
}


@pytest.fixture
def populated_db(session):
    """Insert fixture species, genes and orthologs into the test session."""
    from db.models import Gene, Ortholog, Species

    # Insert species
    for entry in FIXTURE_SPECIES:
        s = Species(
            id=entry["id"],
            taxid=entry["taxid"],
            scientific_name=entry["name"],
            phenotypes=entry["phenotype"],
            genome_assembly=entry.get("genome_assembly"),
            lineage_group=entry.get("lineage_group"),
        )
        session.add(s)

    # Insert genes + orthologs
    for gene_entry in FIXTURE_GENES:
        gene_id = str(uuid.uuid4())
        gene = Gene(
            id=gene_id,
            human_gene_id=gene_entry["human_gene_id"],
            gene_symbol=gene_entry["gene_symbol"],
            human_protein=gene_entry["human_protein"],
        )
        session.add(gene)

        for species_id, seq in FIXTURE_SEQS.items():
            ortholog = Ortholog(
                gene_id=gene_id,
                species_id=species_id,
                protein_id=f"{species_id}|{gene_entry['human_protein']}",
                protein_seq=seq,
                orthofinder_og=f"OG_{gene_entry['human_gene_id']}",
                sequence_identity_pct=99.0 if species_id == "human" else 85.0,
            )
            session.add(ortholog)

    session.flush()
    return session
