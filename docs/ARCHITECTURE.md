# Architecture Overview

- `app.py`: Streamlit entry point, navigation, shared session state and global presentation layer.
- `config/settings.py`: one source of truth for NCBI contact configuration and project metadata.
- `core/data_manager.py`: session-scoped transcript registry shared by transcript retrieval, MEF2C and primer workflows.
- `modules/ncbi_client.py`: common Entrez retrieval functions; delegates credential setup to central configuration.
- `modules/mef2c_analysis.py`: live MEF2C transcript catalogue, coding/protein mapping and conservative transcript-specificity feasibility assessment.
- `modules/primer_design.py`, `modules/rtqpcr_primer.py`, `modules/ncbi_primer_screen.py`: sequence alignment, Primer3 candidate design and NCBI BLAST-based screening.
- `modules/gckr_characterization.py`: gene/protein annotation integrations.
- `modules/mutation_analyzer.py`: coding SNV consequence calculations.

## Shared data lifecycle

Biopython `SeqRecord` objects are registered by accession in Streamlit session state. Modules can retain and reuse those records without downloading the same sequence again. The registry is session-scoped and is cleared when the Streamlit session ends.

## Credentials

NCBI email/API-key configuration is centralized in `config/settings.py`. The sidebar provides a session-level editor; environment variables or `.streamlit/secrets.toml` provide persistent deployment configuration. No API key is included in the source archive.
