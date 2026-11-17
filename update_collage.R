.libPaths(c(.libPaths(), "/usr/lib/R/site-library"))
library(httr)
library(jsonlite)
library(dplyr)
library(ggplot2)
library(png)
library(grid)

station_token <- Sys.getenv("BIRDWEATHER_TOKEN")
stopifnot("Missing BIRDWEATHER_TOKEN secret" = nchar(station_token) > 0)

# 1. Fetch ALL 24-hour species counts directly from BirdWeather's /species endpoint
species_url <- paste0(
  "https://", paste(c("app", "birdweather", "com"), collapse = "."), "/",
  paste(c("api", "v1", "stations", station_token, "species"), collapse = "/")
)

res <- GET(
  url = species_url,
  query = list(period = "day", limit = 200),
  add_headers("Accept" = "application/json", "User-Agent" = "Mozilla/5.0")
)
stop_for_status(res)
parsed <- fromJSON(content(res, as = "text", encoding = "UTF-8"), flatten = TRUE)

# Handle list structure returned by /species endpoint
sp_df <- if (!is.null(parsed$species)) parsed$species else as.data.frame(parsed)
message("Columns returned by /species: ", paste(names(sp_df), collapse = ", "))

common_col <- grep("commonName|common_name", names(sp_df), value = TRUE)
sci_col    <- grep("scientificName|scientific_name", names(sp_df), value = TRUE)
count_col  <- grep("detections_count|detections.total|count|total", names(sp_df), value = TRUE)

flock_df <- sp_df %>%
  transmute(
    common_name = as.character(.data[[common_col]]),
    sci_name    = as.character(.data[[sci_col]]),
    detections  = as.numeric(.data[[count_col]])
  ) %>%
  filter(!is.na(common_name), !is.na(sci_name), detections > 0) %>%
  mutate(clean_sci = sub("^(\\S+\\s+\\S+).*", "\\1", sci_name)) %>%
  group_by(common_name, clean_sci) %>%
  summarise(detections = sum(detections), .groups = "drop") %>%
  arrange(desc(detections)) %>%
  mutate(
    pose = ifelse(row_number() %% 2 == 1, 1, 2),
    area_score = detections^0.45
  )

message("24-hour unique species fetched: ", nrow(flock_df))

# 2. Fetch full-color woodblock PNGs (keeps ALL species even if a cutout is missing)
get_bird_art <- function(sci_name, common_name, pose = 1) {
  cutout_url <- paste0(
    "https://", paste(c("bird", "onethreenine", "net"), collapse = "."),
    "/avian/api/cutout.php"
  )
  res <- try(GET(
    url = cutout_url,
    query = list(sci = sci_name, com = common_name, pose = pose),
    add_headers("User-Agent" = "Mozilla/5.0"),
    timeout(10)
  ), silent = TRUE)
  
  if (!inherits(res, "try-error") && status_code(res) == 200) {
    img <- try(readPNG(content(res, as = "raw")), silent = TRUE)
    if (!inherits(img, "try-error")) return(img)
  }
  return(NULL)
}

rasters <- mapply(get_bird_art, flock_df$clean_sci, flock_df$common_name, flock_df$pose, SIMPLIFY = FALSE)
has_art <- !sapply(rasters, is.null)
message("Woodblock art found for ", sum(has_art), " of ", nrow(flock_df), " species")

# 3. True Non-Overlapping Elliptical Spiral Packing (accounts for both bird + text label)
pack_spiral <- function(areas) {
  r <- sqrt(areas / pi)
  n <- length(r)
  x <- numeric(n)
  y <- numeric(n)
  if (n > 1) {
    angles <- seq(0, 40 * pi, length.out = 5000)
    for (i in 2:n) {
      step <- r[i] * 0.02
      rad <- 0
      for (ang in angles) {
        cx <- rad * cos(ang) * 1.45
        cy <- rad * sin(ang)
        # Enforce strict spacing: horizontal gap >= 1.25x sum of radii, vertical >= 1.35x sum of radii
        dx <- abs(cx - x[1:(i - 1)]) / ((r[i] + r[1:(i - 1)]) * 1.25)
        dy <- abs(cy - y[1:(i - 1)]) / ((r[i] + r[1:(i - 1)]) * 1.35)
        if (all(sqrt(dx^2 + dy^2) >= 1.0)) {
          x[i] <- cx
          y[i] <- cy
          break
        }
        rad <- rad + step
      }
    }
  }
  tibble(x_wide = x, y = y, radius = r)
}

packing <- pack_spiral(flock_df$area_score)
flock_df <- bind_cols(flock_df, packing) %>% mutate(has_art = has_art)

# 4. Render warm cream woodblock collage
bc_time <- format(as.POSIXct(Sys.time(), tz = "America/Vancouver"), "%b %d, %I:%M %p %Z")
total_calls <- sum(flock_df$detections)

p <- ggplot(flock_df, aes(x = x_wide, y = y)) +
  geom_blank() +
  coord_equal(expand = TRUE) +
  scale_x_continuous(expand = expansion(mult = 0.18)) +
  scale_y_continuous(expand = expansion(mult = 0.18)) +
  labs(
    title = "Bird Detections",
    subtitle = paste0("Last 24 Hours: ", nrow(flock_df), " Species (", total_calls, " Detections) • Updated ", bc_time)
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "#f6f0e2", color = NA),
    panel.background = element_rect(fill = "#f6f0e2", color = NA),
    plot.title = element_text(
      hjust = 0.5, family = "serif", size = 18, face = "bold",
      color = "#2b241d", margin = margin(t = 15, b = 3)
    ),
    plot.subtitle = element_text(
      hjust = 0.5, family = "serif", size = 11, color = "#6e6051", margin = margin(b = 8)
    ),
    plot.margin = margin(15, 15, 15, 15)
  )

# Draw woodblock art (sized inside 80% of padded cell so birds & labels never touch)
for (i in seq_len(nrow(flock_df))) {
  if (flock_df$has_art[i]) {
    img <- rasters[[i]]
    h_px <- nrow(img)
    w_px <- ncol(img)
    scale_factor <- (flock_df$radius[i] * 1.45) / sqrt(h_px^2 + w_px^2)
    half_w <- (w_px * scale_factor) / 2
    half_h <- (h_px * scale_factor) / 2
    
    p <- p + annotation_raster(
      img,
      xmin = flock_df$x_wide[i] - half_w,
      xmax = flock_df$x_wide[i] + half_w,
      ymin = flock_df$y[i] - half_h + flock_df$radius[i] * 0.12,
      ymax = flock_df$y[i] + half_h + flock_df$radius[i] * 0.12
    )
  } else {
    p <- p + annotate(
      "point",
      x = flock_df$x_wide[i], y = flock_df$y[i] + flock_df$radius[i] * 0.12,
      size = flock_df$radius[i] * 4.0, shape = 21,
      fill = "#e3d6be", color = "#8c7a65", stroke = 0.8
    )
  }
}

# Place every species label cleanly in the reserved bottom space of its cell
p <- p + geom_text(
  data = flock_df,
  aes(x = x_wide, y = y - radius * 0.85, label = paste0(common_name, " (", detections, ")")),
  size = 2.6, family = "serif", color = "#3b3228"
)

ggsave("avian_visitors.png", plot = p, width = 12, height = 8, dpi = 200, bg = "#f6f0e2")
