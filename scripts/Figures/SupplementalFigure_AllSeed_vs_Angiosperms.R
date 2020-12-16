library(dplyr)
library(ggplot2)
library(sf)
library(raster)
library(pals)

#' Script to generate Suplemental Figure, which Maps PD, RPD, and PE of all 
#' seed plants and displays them next to angiosperms
#' 
#' Note the degree of similarity


# First make the all seed plant plots
###### PD AND RPD PLOTS ################
# na shapefile
na_maps <- maps::map(regions=c("usa", "mexico", "canada"), plot = FALSE, fill = TRUE) %>%
  st_as_sf()
na_but <- st_transform(na_maps, 
                       crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")

# PD 
pd <- raster("data/AllSeedPlants/spres_PD_P.tif")
crs(pd) <- crs(na_but)
pd_df <- as.data.frame(pd, xy = TRUE)

## clip out the non US, MEXICO, and CA areas
pd_sf <- st_as_sf(pd_df,
                  coords = c("x", "y"),
                  crs = st_crs(na_but))

pd_coords <- st_join(pd_sf, na_but, join=st_is_within_distance, dist = 1000) %>% 
  st_coordinates() %>% 
  as.data.frame()
pd_clip <- st_join(pd_sf, na_but, join=st_is_within_distance, dist = 1000) %>% 
  st_drop_geometry()
pd_clip <- dplyr::bind_cols(pd_clip, pd_coords)


pd_plot <- ggplot() +
  geom_tile(pd_clip, mapping = aes(x = X, y = Y, fill = spres_PD_P)) +
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) + 
  scale_fill_viridis_c(trans = "log", na.value = "transparent", 
                       option = "plasma") +
  labs(fill = "PD") + 
  theme_void()

## RPD Plot time
rpd <- raster("data/AllSeedPlants/spres_PHYLO_RPD2.tif")
crs(rpd) <- crs(na_but)
rpd_df <- as.data.frame(rpd, xy = TRUE)

## clip out the non US, MEXICO, and CA areas
rpd_sf <- st_as_sf(rpd_df,
                   coords = c("x", "y"),
                   crs = st_crs(na_but))

rpd_coords <- st_join(rpd_sf, na_but, join=st_is_within_distance, dist = 1000) %>% 
  st_coordinates() %>% 
  as.data.frame()
rpd_clip <- st_join(rpd_sf, na_but, join=st_is_within_distance, dist = 1000) %>% 
  st_drop_geometry()
rpd_clip <- dplyr::bind_cols(rpd_clip, rpd_coords)


rpd_plot <- ggplot() +
  geom_tile(rpd_clip, mapping = aes(x = X, y = Y, fill = spres_PHYLO_RPD2)) +
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) + 
  scale_fill_gradientn(trans = "log", na.value = "transparent",
                       colours = parula(25), guide = "colourbar") +
  labs(fill = "RPD") + 
  theme_void()

################ PE PLOTS ###############

# na shapefile
na_maps <- maps::map(regions=c("usa", "mexico", "canada"), plot = FALSE, fill = TRUE) %>%
  st_as_sf()
na_but <- st_transform(na_maps, 
                       crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")

# PE
pe <- raster("data/AllSeedPlants/spres_PE_WE_P.tif")
crs(pe) <- crs(na_but)
pe_df <- as.data.frame(pe, xy = TRUE)

## clip out the non US, MEXICO, and CA areas
pe_sf <- st_as_sf(pe_df,
                  coords = c("x", "y"),
                  crs = st_crs(na_but))

pe_coords <- st_join(pe_sf, na_but, join=st_is_within_distance, dist = 1000) %>% 
  st_coordinates() %>% 
  as.data.frame()
pe_clip <- st_join(pe_sf, na_but, join=st_is_within_distance, dist = 1000) %>% 
  st_drop_geometry()
pe_clip <- dplyr::bind_cols(pe_clip, pe_coords)

pe_plot <- ggplot() +
  geom_tile(pe_clip, mapping = aes(x = X, y = Y, fill = spres_PE_WE_P)) +
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) + 
  scale_fill_viridis_c(trans = "log", na.value = "transparent", 
                       option = "cividis") +
  labs(fill = "PE") +
  theme_void()


## Combine plots into 3 panel, vertical plot

cp_seeds <- cowplot::plot_grid(pd_plot, rpd_plot, pe_plot,
                               ncol = 1, labels = c("A", "B",
                                                    "C"),label_y = 0.85,
                               hjust = -3, vjust = 0.5)


## Now make the angiosperm plots

# PD 
pd_as <- raster("data/Angiosperms/exports/results_PD_P.tiff")
crs(pd_as) <- crs(na_but)
pd_as_df <- as.data.frame(pd_as, xy = TRUE)

## clip out the non US, MEXICO, and CA areas
pd_as_sf <- st_as_sf(pd_as_df,
                     coords = c("x", "y"),
                     crs = st_crs(na_but))

pd_as_coords <- st_join(pd_as_sf, na_but, join=st_is_within_distance, dist = 1000) %>% 
  st_coordinates() %>% 
  as.data.frame()
pd_as_clip <- st_join(pd_as_sf, na_but, join=st_is_within_distance, dist = 1000) %>% 
  st_drop_geometry()
pd_as_clip <- dplyr::bind_cols(pd_as_clip, pd_as_coords)


pd_as_plot <- ggplot() +
  geom_tile(pd_as_clip, mapping = aes(x = X, y = Y, fill = results_PD_P)) +
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) + 
  scale_fill_viridis_c(trans = "log", na.value = "transparent", 
                       option = "plasma") +
  labs(fill = "PD") + 
  theme_void()

## RPD Plot time
rpd_as <- raster("data/Angiosperms/exports/results_PHYLO_RPD2.tiff")
crs(rpd_as) <- crs(na_but)
rpd_as_df <- as.data.frame(rpd_as, xy = TRUE)

## clip out the non US, MEXICO, and CA areas
rpd_as_sf <- st_as_sf(rpd_as_df,
                      coords = c("x", "y"),
                      crs = st_crs(na_but))

rpd_as_coords <- st_join(rpd_as_sf, na_but, join=st_is_within_distance, dist = 1000) %>% 
  st_coordinates() %>% 
  as.data.frame()
rpd_as_clip <- st_join(rpd_as_sf, na_but, join=st_is_within_distance, dist = 1000) %>% 
  st_drop_geometry()
rpd_as_clip <- dplyr::bind_cols(rpd_as_clip, rpd_as_coords)


rpd_as_plot <- ggplot() +
  geom_tile(rpd_as_clip, mapping = aes(x = X, y = Y, fill = results_PHYLO_RPD2)) +
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) + 
  scale_fill_gradientn(trans = "log", na.value = "transparent",
                       colours = parula(25), guide = "colourbar", 
                       breaks = c(0.5,1,2, 3.3)) +
  labs(fill = "RPD") + 
  theme_void()

################ PE PLOTS ###############

# na shapefile
na_maps <- maps::map(regions=c("usa", "mexico", "canada"), plot = FALSE, fill = TRUE) %>%
  st_as_sf()
na_but <- st_transform(na_maps, 
                       crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")

# pe_as
pe_as <- raster("data/Angiosperms/exports/results_PE_WE_P.tiff")
crs(pe_as) <- crs(na_but)
pe_as_df <- as.data.frame(pe_as, xy = TRUE)

## clip out the non US, MEXICO, and CA areas
pe_as_sf <- st_as_sf(pe_as_df,
                     coords = c("x", "y"),
                     crs = st_crs(na_but))

pe_as_coords <- st_join(pe_as_sf, na_but, join=st_is_within_distance, dist = 1000) %>% 
  st_coordinates() %>% 
  as.data.frame()
pe_as_clip <- st_join(pe_as_sf, na_but, join=st_is_within_distance, dist = 1000) %>% 
  st_drop_geometry()
pe_as_clip <- dplyr::bind_cols(pe_as_clip, pe_as_coords)

pe_as_plot <- ggplot() +
  geom_tile(pe_as_clip, mapping = aes(x = X, y = Y, fill = results_PE_WE_P)) +
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) + 
  scale_fill_viridis_c(trans = "log", na.value = "transparent", 
                       option = "cividis") +
  labs(fill = "PE") +
  theme_void()


## Combine plots into 3 panel, vertical plot

cp_angiosperms <- cowplot::plot_grid(pd_as_plot, rpd_as_plot, pe_as_plot,
                                     ncol = 1, labels = c("D", "E",
                                                          "F"),label_y = 0.85,
                                     hjust = -3, vjust = 0.5)

## Combine two three panel plots
seeds_vs_angiosperms <- cowplot::plot_grid(cp_seeds, cp_angiosperms,
                                           labels = c("Seed Plants",
                                                      "Angiosperms"))

ggsave(filename = "figure_outputs/SupplementalFigure_seed_vs_angiosperms.png", plot = seeds_vs_angiosperms,
       dpi = 300, device = "png", width = 15, height = 12)
