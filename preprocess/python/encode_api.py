"""Small helper for the ENCODE REST API (https://www.encodeproject.org/help/rest-api/)."""

from __future__ import annotations

import time

import requests

BASE = "https://www.encodeproject.org"
HEADERS = {"accept": "application/json"}
_session = requests.Session()
_session.headers.update(HEADERS)


def _request(url: str, params=None, retries: int = 5, timeout: int = 300) -> dict:
    for attempt in range(retries):
        try:
            r = _session.get(url, params=params, timeout=timeout)
            if r.status_code == 404 and "/search/" in url:
                return {"@graph": []}  # ENCODE returns 404 for an empty search
            r.raise_for_status()
            return r.json()
        except (requests.RequestException, ValueError):
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def search(params: dict, released_only: bool = True) -> list[dict]:
    """Run a /search/ query and return every hit."""
    q = {"format": "json", "limit": "all", **params}
    if released_only:
        q.setdefault("status", "released")
    return _request(f"{BASE}/search/", params=q).get("@graph", [])


def search_by_dataset(accessions, params: dict, batch: int = 60) -> list[dict]:
    """Search File objects belonging to many experiments, batching the dataset filter."""
    accessions = sorted(set(accessions))
    out: list[dict] = []
    for i in range(0, len(accessions), batch):
        chunk = accessions[i : i + batch]
        out += search({**params, "type": "File", "dataset": [f"/experiments/{a}/" for a in chunk]})
    return out


def get(path: str) -> dict:
    """Fetch one object, e.g. get('/files/ENCFF590SSX/')."""
    return _request(f"{BASE}{path}", params={"format": "json"})


def audit_summary(obj: dict, levels=("ERROR", "NOT_COMPLIANT")) -> str:
    """Compact text summary of the serious audit flags on an object."""
    a = obj.get("audit", {}) or {}
    parts = sorted({f"{lvl}:{i.get('category')}" for lvl in levels for i in a.get(lvl, [])})
    return ";".join(parts)
