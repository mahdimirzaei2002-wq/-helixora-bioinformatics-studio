# Helixora deployment (VPS + Docker + Caddy)

This bundle serves the bilingual public website at `/` and the Streamlit research studio at `/app/` on the same domain.

## Requirements

- A Linux VPS (Ubuntu 22.04/24.04 recommended) with at least 2 GB RAM; 4 GB is preferable for scientific Python dependencies.
- A domain name with DNS `A` (and, if used, `AAAA`) records pointing to the VPS public IP.
- Inbound TCP ports 80 and 443 open in the VPS firewall/security group.
- Docker Engine and the Docker Compose plugin.

## Deploy

1. Copy this project folder to the VPS (for example with `scp` or Git).
2. Install Docker Engine and Docker Compose on the VPS.
3. Copy `.env.example` to `.env` and set `DOMAIN=your-domain.example`. Optionally add `NCBI_EMAIL` and `NCBI_API_KEY`. Do not publish `.env` or commit it to Git.
4. Start the services from this directory:

   ```bash
   docker compose up -d --build
   ```

5. Open `https://your-domain.example`. Caddy obtains and renews TLS certificates automatically when the domain DNS and ports are configured correctly. The app is available at `https://your-domain.example/app/`.
6. Check logs if needed:

   ```bash
   docker compose logs -f app
   docker compose logs -f web
   ```

## Local smoke test

Set `DOMAIN=localhost` in `.env`, then run `docker compose up --build`. Open `http://localhost` and `http://localhost/app/`.

## Operations and security

- Back up the `.env` file securely and keep it outside public repositories.
- Configure a firewall to expose only SSH, HTTP and HTTPS. The Streamlit port is internal to the Docker network and is not published directly.
- Use SSH keys, disable password-based root login where practical, and keep the OS and Docker images updated.
- The app currently uses Streamlit session state for its shared active-record workflow; it is not a multi-user database or persistent project store.
- External NCBI/UniProt/GTEx/EBI/BLAST services require outbound internet and may enforce rate limits.
- This is an academic research/education tool, not a clinical decision-support system.

## What is not included

This repository contains deployment configuration, not a purchased VPS, registered domain, DNS configuration, or cloud account. Those resources must be supplied and authorized by the project owner.
