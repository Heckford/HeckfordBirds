# qathet: Myrtle Rock Bird and Bat Detections

Live webpage displaying automated 48-hour and Year-to-Date (YTD) acoustic bird and bat detection collages from two BirdWeather stations in Myrtle Rocks (qathet Regional District, BC).

## Live Website
- **GitHub Pages URL:** `https://github.io`

## Architecture & Files
- **`update_collage.R`** (organized into 5 annotated sections):
  - **Section 1**: Fetches full species summaries (`All` detections) and uncapped timestamped detections (`>= 60%` confidence) from two BirdWeather station tokens (`BIRDWEATHER_TOKEN` and `BIRDWEATHER_TOKEN_2`).
  - **Section 2**: Downloads the official eBird taxonomy table (`speciesCode` and taxonomic `order`), caches woodblock bird illustrations, and defines fallback bird/bat vector silhouettes.
  - **Section 3**: Maps every species into deterministic taxonomic groups (**Songbirds**, **Raptors**, **Waterfowl & Seabirds**, **Shorebirds & Gulls**, **Woodpeckers & Other Birds**, and **Bats**) with zero dropped species, and packs them into ranked horizontal field-guide rows.
  - **Section 4 (4A & 4B)**: Renders framed sepia panels with two-line stacked labels (`Common Name` + formatted `(count)`) and exports 10 PNG collages plus 10 matching JSON clickable hotspot maps (`48h` and `ytd` across `All`, `60%`, `70%`, `80%`, and `90%` confidence).
  - **Section 5**: Sets `America/Vancouver` timeframe cutoffs and executes the full 10-collage rendering loop.
- **`index.html`**: Renders the warm cream webpage with top toggle pills for **Confidence (`All | 60% | 70% | 80% | 90%`)** and **Timeframe (`Last 48 Hours | Year to Date`)**, plus transparent clickable links over every bird (**eBird**) and bat (**iNaturalist**).
- **`.github/workflows/update.yml`**: Runs `update_collage.R` inside `rocker/tidyverse:4.4.1` and commits updated PNGs (`HeckfordBirds_*.png`) and JSON maps (`species_map_*.json`) to `main`.

## Automation & Secrets
- **GitHub Repository Secrets** (*Settings > Secrets and variables > Actions*):
  - `BIRDWEATHER_TOKEN` (Station 1 API token)
  - `BIRDWEATHER_TOKEN_2` (Station 2 API token)
- **External Scheduler (`cron-job.org`)**:
  - Sends a `POST` request every 15 minutes (`:00`, `:15`, `:30`, `:45`) to `https://github.com` with body `{"ref":"main"}` using a fine-grained GitHub Personal Access Token (`Actions: Read and write`).

## Safe Testing Workflow
- Create a temporary branch (e.g., `test-feature`) from `main` before testing major layout or script changes, verify the rendered PNGs in GitHub Actions, apply verified code changes to `main`, and delete the temporary branch once finished.
