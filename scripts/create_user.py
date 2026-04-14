#!/usr/bin/env python3
"""Create or update a BioResilient platform user.

Usage:
    python scripts/create_user.py --email admin@example.com --name "Admin User" --password changeme
    python scripts/create_user.py --email researcher@lab.org --name "Researcher" --role viewer
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from passlib.context import CryptContext

_pwd = CryptContext(schemes=["bcrypt"], deprecated="auto")


def create_user(email: str, name: str, password: str, role: str = "viewer") -> None:
    from db.session import get_session
    from db.models import User

    email = email.lower().strip()

    with get_session() as session:
        existing = session.query(User).filter_by(email=email).first()
        if existing:
            existing.name = name
            existing.hashed_password = _pwd.hash(password)
            existing.role = role
            existing.is_active = True
            action = "updated"
        else:
            session.add(User(
                email=email,
                name=name,
                hashed_password=_pwd.hash(password),
                role=role,
                is_active=True,
            ))
            action = "created"

    print(f"[ok] User {email} ({role}) {action}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create or update a BioResilient user")
    parser.add_argument("--email", required=True, help="User email address")
    parser.add_argument("--name", required=True, help="Display name")
    parser.add_argument("--password", required=True, help="Password (will be hashed)")
    parser.add_argument(
        "--role", default="viewer",
        choices=["admin", "operator", "viewer"],
        help="Role: admin (full access), operator (can run pipeline), viewer (read-only)",
    )
    args = parser.parse_args()
    create_user(args.email, args.name, args.password, args.role)
