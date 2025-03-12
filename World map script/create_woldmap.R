library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(dplyr)
library(readr)
library(tidyr)
library(rworldmap)

file_path <- "./prism_cleaned.csv"
data_points <- read_csv(file_path)

data_points <- data_points %>%
  mutate(country = case_when(
    country == "United States" ~ "United States of America",
    TRUE ~ country
  ))

data_points <- data_points %>%
  group_by(country, diag_defn) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = diag_defn, values_from = count, values_fill = list(count = 0)) 

colnames(data_points) <- tolower(colnames(data_points))  
data_points <- data_points %>%
  rename(gold_count = gold, lln_count = lln)  

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  select(iso_a3, name, geometry)

data_points <- data_points %>%
  left_join(world, by = c("country" = "name"))

data_points <- data_points %>% filter(!is.na(geometry))

generate_random_points <- function(geometry, n, buffer_dist) {
  if (n > 0 && !is.null(geometry)) {
    buffered_geometry <- st_buffer(geometry, buffer_dist)  
    sampled_points <- st_sample(buffered_geometry, n) 
    return(st_coordinates(st_sample(geometry, n))) 
  } else {
    return(matrix(NA, nrow = 0, ncol = 2)) 
  }
}

gold_points <- do.call(rbind, mapply(generate_random_points, data_points$geometry, data_points$gold_count, MoreArgs = list(buffer_dist = -0.3), SIMPLIFY = FALSE))
lln_points <- do.call(rbind, mapply(generate_random_points, data_points$geometry, data_points$lln_count, MoreArgs = list(buffer_dist = -0.3), SIMPLIFY = FALSE))

#gold_points <- do.call(rbind, mapply(generate_random_points, data_points$geometry, data_points$gold_count, SIMPLIFY = FALSE))
#lln_points <- do.call(rbind, mapply(generate_random_points, data_points$geometry, data_points$lln_count, SIMPLIFY = FALSE))

gold_df <- data.frame(longitude = gold_points[, 1], latitude = gold_points[, 2], type = "GOLD")
lln_df <- data.frame(longitude = lln_points[, 1], latitude = lln_points[, 2], type = "LLN")

plot_data <- bind_rows(gold_df, lln_df)


world <- ne_countries(scale = "large", returnclass = "sf")
world <- world %>% filter(continent != "Antarctica")


who_regions <- data.frame(
  iso_a3 = c(
    # African Region (AFRO) - 47 countries
    "DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD", "COM",
    "COD", "COG", "CIV", "GNQ", "ERI", "SWZ", "ETH", "GAB", "GMB", "GHA", "GIN",
    "GNB", "KEN", "LSO", "LBR", "MDG", "MWI", "MLI", "MRT", "MUS", "MOZ", "NAM",
    "NER", "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "ZAF", "SSD", "TGO", "UGA",
    "TZA", "ZMB", "ZWE",
    
    # Region of the Americas (PAHO) - 35 countries
    "ATG", "ARG", "BHS", "BRB", "BLZ", "BOL", "BRA", "CAN", "CHL", "COL", "CRI",
    "CUB", "DMA", "DOM", "ECU", "SLV", "GRD", "GTM", "GUY", "HTI", "HND", "JAM",
    "MEX", "NIC", "PAN", "PRY", "PER", "KNA", "LCA", "VCT", "SUR", "TTO", "USA",
    "URY", "VEN",
    
    # South-East Asia Region (SEARO) - 11 countries
    "BGD", "BTN", "IND", "IDN", "MDV", "MMR", "NPL", "LKA", "THA", "TLS", "VNM",
    
    # European Region (EURO) - 53 countries
    "ALB", "AND", "ARM", "AUT", "AZE", "BLR", "BEL", "BIH", "BGR", "HRV", "CYP",
    "CZE", "DNK", "EST", "FIN", "FRA", "GEO", "DEU", "GRC", "HUN", "ISL", "IRL",
    "ISR", "ITA", "KAZ", "KGZ", "LVA", "LIE", "LTU", "LUX", "MLT", "MCO", "MNE",
    "NLD", "MKD", "NOR", "POL", "POR", "MDA", "ROU", "RUS", "SMR", "SRB", "SVK",
    "SVN", "ESP", "SWE", "CHE", "TJK", "TKM", "UKR", "GBR", "UZB",
    
    # Eastern Mediterranean Region (EMRO) - 21 countries
    "AFG", "BHR", "DJI", "EGY", "IRN", "IRQ", "JOR", "KWT", "LBN", "LBY", "MAR",
    "OMN", "PAK", "QAT", "SAU", "SOM", "SDN", "SYR", "TUN", "ARE", "YEM",
    
    # Western Pacific Region (WPRO) - 27 countries
    "AUS", "BRN", "KHM", "CHN", "FJI", "PYF", "GUM", "JPN", "KIR", "PRK", "KOR",
    "LAO", "MYS", "MHL", "FSM", "MNG", "NRU", "NCL", "NZL", "PLW", "PNG", "PHL",
    "WSM", "SGP", "SLB", "TWN", "TON"
  ),
  
  who_region = c(
    rep("African region (AFRO)", 47),
    rep("Region of the Americas (PAHO)", 35),
    rep("South-East Asia region (SEARO)", 11),
    rep("European region (EURO)", 53),
    rep("Eastern Mediterranean region (EMRO)", 21),
    rep("Western Pacific region (WPRO)", 27)
  )
)
world <- left_join(world, who_regions, by = "iso_a3")
who_colors <- c(
  "African region (AFRO)" = "gold",
  "Region of the Americas (PAHO)" = "lightcoral",
  "South-East Asia region (SEARO)" = "lightgreen",
  "European region (EURO)" = "lightpink",
  "Eastern Mediterranean region (EMRO)" = "darkorange2",
  "Western Pacific region (WPRO)" = "lightblue",
  "Other" = "lightyellow"
)
world$who_region[is.na(world$who_region)] <- "Other"

p <- ggplot(data = world) +
  geom_sf(aes(fill = who_region), color = "white", size = 0.1) + 
  scale_fill_manual(
    values = who_colors, 
    name = "WHO Regional Classification"  
  ) +  
  coord_sf(ylim = c(-55, 90)) + 
  ggtitle("Location of studies reporting PRISm and/or RSP prevalence using GOLD and LLN definitions") +
  
  
  geom_point(data = plot_data %>% filter(type == "GOLD"), aes(x = longitude, y = latitude), 
             shape = 21, fill = "darkred", color = "black", alpha = 0.8, size = 2.0) +
  
  geom_point(data = plot_data %>% filter(type == "LLN"), aes(x = longitude, y = latitude), 
             shape = 24, fill = "darkgreen", color = "black", alpha = 0.8, size = 1.8) +
  
  
  annotate("rect", xmin = -170, xmax = -140, ymin = -61, ymax = -40, fill = "white", color = "black") +

  annotate("text", x = -169, y = -43, label = "Data Type", hjust = 0, size = 3.2, fontface = "bold") +
  annotate("text", x = -160, y = -50, label = "LLN", hjust = 0, size = 3.2, fontface = "bold") +
  annotate("text", x = -160, y = -57, label = "GOLD", hjust = 0, size = 3.2, fontface = "bold") +
  
  geom_point(aes(x = -164, y = -57), shape = 21, fill = "darkred", color = "black", size = 3) +
  geom_point(aes(x = -164, y = -51), shape = 24, fill = "darkgreen", color = "black", size = 2) +
  
  
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    
    plot.title = element_text(
      hjust = 0.5,       
      size = 10,         
      face = "bold",     
      color = "black" 
    ),
    
    legend.position = c(0.95, -0.14),  
    legend.justification = c(1, 0),  
    legend.box = "vertical",
    legend.background = element_rect(fill = "white", color = "black"), 
    legend.key.size = unit(0.4, "cm"),  
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 8)  
  ) +
  guides(fill = guide_legend(title = "WHO Region", ncol = 2, byrow = TRUE)) 




print(p)