"""Central configuration and NCBI credentials for Helixora Bioinformatics Studio."""
from __future__ import annotations
import os
from dataclasses import dataclass
from typing import Any

DEFAULT_CONTACT_EMAIL = ""
APP_NAME = "Helixora Bioinformatics Studio"
DEVELOPER_NAME = "Mahdi Mirzaei"
DEVELOPER_EMAIL = DEFAULT_CONTACT_EMAIL
DEVELOPER_PHONE = ""


def _secret(key: str, default: str = "") -> str:
    try:
        import streamlit as st
        value = st.secrets.get(key, default)
        return str(value) if value is not None else default
    except Exception:
        return os.getenv(key, default)

@dataclass(frozen=True)
class AppSettings:
    ncbi_email: str = DEFAULT_CONTACT_EMAIL
    ncbi_api_key: str = ""
    request_timeout: int = 30
    app_name: str = APP_NAME

    @classmethod
    def load(cls) -> "AppSettings":
        return cls(
            ncbi_email=_secret("NCBI_EMAIL", DEFAULT_CONTACT_EMAIL),
            ncbi_api_key=_secret("NCBI_API_KEY", ""),
            request_timeout=int(_secret("REQUEST_TIMEOUT", "30")),
        )


def configure_entrez(email: str | None = None, api_key: str | None = None) -> None:
    """Apply one shared Entrez configuration used by all NCBI modules."""
    from Bio import Entrez
    settings = AppSettings.load()
    resolved_email = (email or settings.ncbi_email).strip()
    resolved_key = settings.ncbi_api_key if api_key is None else api_key.strip()
    if not resolved_email or "@" not in resolved_email:
        raise ValueError("A valid NCBI contact email is required. Configure it in Settings.")
    Entrez.email = resolved_email
    Entrez.tool = "BioinformaticsDashboardPro"
    if resolved_key:
        Entrez.api_key = resolved_key
    elif hasattr(Entrez, "api_key"):
        Entrez.api_key = None
