# Qathet: Myrtle Rock Bird and Bat Detections

Live webpage displaying automated 48-hour and Year-to-Date (YTD) acoustic bird and bat detection collages from two BirdWeather stations in Myrtle Rocks (Powell River, BC).

## Live Website
- **GitHub Pages URL:** `https://github.io`

## Architecture & Files
- **`update_collage.R`**: Fetches detections from two BirdWeather station tokens (`BIRDWEATHER_TOKEN` and `BIRDWEATHER_TOKEN_2`), retrieves woodblock bird illustrations and eBird/iNaturalist taxonomy links, packs birds into a central flock and bats into a bottom-right corner box, and exports 10 PNG collages plus 10 matching JSON hotspot maps (`48h` and `ytd` across `All`, `60%`, `70%`, `80%`, and `90%` confidence).
- **`index.html`**: Renders the warm cream webpage with top toggle pills for **Confidence (`All | 60% | 70% | 80% | 90%`)** and **Timeframe (`Last 48 Hours | Year to Date`)**, plus transparent clickable links over every bird (opening **eBird**) and bat (opening **iNaturalist**).
- **`.github/workflows/update.yml`**: Runs `update_collage.R` inside `rocker/tidyverse:4.4.1` and commits updated PNGs and JSON maps to `main`.

## Automation & Secrets
- **GitHub Repository Secrets** (*Settings > Secrets and variables > Actions*):
  - `BIRDWEATHER_TOKEN` (Station 1 API token)
  - `BIRDWEATHER_TOKEN_2` (Station 2 API token)
- **External Scheduler (`cron-job.org`)**:
  - Sends a `POST` request every 15 minutes (`:00`, `:15`, `:30`, `:45`) to `https://github.com` with body `{"ref":"main"}` using a fine-grained GitHub Personal Access Token (`Actions: Read and write`).

## Safe Testing Workflow
- Create a temporary branch (e.g., `test-feature`) from `main` before testing major layout or script changes, verify the rendered PNGs in GitHub Actions, and apply verified code changes to `main`.
