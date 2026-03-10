"""Seed the database with species from config/species_registry.json.

Run after `alembic upgrade head`:
    python db/seed.py
"""

import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from db.models import Species
from db.session import get_session

logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")
log = logging.getLogger(__name__)

REGISTRY_PATH = Path(__file__).resolve().parents[1] / "config" / "species_registry.json"


def seed_species() -> None:
    with open(REGISTRY_PATH) as f:
        species_list = json.load(f)

    with get_session() as session:
        existing_ids = {s.id for s in session.query(Species).all()}
        added = 0
        for entry in species_list:
            if entry["id"] in existing_ids:
                log.info("  skip (already exists): %s", entry["id"])
                continue
            species = Species(
                id=entry["id"],
                taxid=entry["taxid"],
                scientific_name=entry["name"],
                phenotypes=entry.get("phenotype", []),
                genome_assembly=entry.get("genome_assembly"),
                lineage_group=entry.get("lineage_group"),
                geo_search_terms=entry.get("geo_search_terms", []),
                is_control=entry.get("is_control", False),
            )
            session.add(species)
            added += 1
            log.info("  added: %s (%s)", entry["id"], entry["name"])

    log.info("Seeded %d species (%d already present).", added, len(existing_ids))


if __name__ == "__main__":
    seed_species()
