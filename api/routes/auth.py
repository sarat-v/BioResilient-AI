"""Authentication routes for BioResilient AI.

POST /auth/login   — email + password → JWT access token
GET  /auth/me      — validate token, return current user info
"""

import os
from datetime import datetime, timedelta

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

router = APIRouter()

_ALGORITHM = "HS256"
_TOKEN_EXPIRE_DAYS = 7


def _secret() -> str:
    key = os.environ.get("JWT_SECRET_KEY", "")
    if not key:
        raise RuntimeError(
            "JWT_SECRET_KEY env var not set. "
            "Generate one with: python -c \"import secrets; print(secrets.token_hex(32))\""
        )
    return key


def _make_token(payload: dict) -> str:
    from jose import jwt
    data = {**payload, "exp": datetime.utcnow() + timedelta(days=_TOKEN_EXPIRE_DAYS)}
    return jwt.encode(data, _secret(), algorithm=_ALGORITHM)


def _decode_token(token: str) -> dict:
    from jose import jwt, JWTError
    try:
        return jwt.decode(token, _secret(), algorithms=[_ALGORITHM])
    except JWTError as exc:
        raise HTTPException(status_code=401, detail="Invalid or expired token") from exc


# ---------------------------------------------------------------------------
# Request / Response models
# ---------------------------------------------------------------------------

class LoginRequest(BaseModel):
    email: str
    password: str


class UserOut(BaseModel):
    id: str
    email: str
    name: str
    role: str


class TokenResponse(BaseModel):
    access_token: str
    token_type: str = "bearer"
    user: UserOut


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

@router.post("/login", response_model=TokenResponse)
def login(req: LoginRequest):
    from passlib.context import CryptContext
    from db.session import get_session
    from db.models import User

    pwd = CryptContext(schemes=["bcrypt"], deprecated="auto")

    with get_session() as session:
        user = session.query(User).filter_by(email=req.email.lower().strip()).first()
        if not user or not pwd.verify(req.password, user.hashed_password):
            raise HTTPException(status_code=401, detail="Invalid email or password")
        if not user.is_active:
            raise HTTPException(status_code=403, detail="Account is disabled")

        token = _make_token({
            "sub": user.id,
            "email": user.email,
            "name": user.name,
            "role": user.role,
        })
        return TokenResponse(
            access_token=token,
            user=UserOut(id=user.id, email=user.email, name=user.name, role=user.role),
        )


@router.get("/me", response_model=UserOut)
def me(request: Request):
    """Validate a stored token and return the current user's info.

    Called by the frontend on startup to confirm the stored token is still valid.
    The /auth prefix is exempt from JWTMiddleware, so we decode here directly.
    """
    auth = request.headers.get("Authorization", "")
    if not auth.startswith("Bearer "):
        raise HTTPException(status_code=401, detail="Token required")
    payload = _decode_token(auth[len("Bearer "):])
    return UserOut(
        id=payload["sub"],
        email=payload["email"],
        name=payload["name"],
        role=payload["role"],
    )
