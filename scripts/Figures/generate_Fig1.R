library(dplyr)
library(ggplot2)
library(sf)
library(raster)

#' Script to generate Figure 1, which Maps 
#' a) taxic richness
#' b) PD
#' c) RPD
#' d) PE

## Figure of species richness #########
fn <- readr::read_csv("data/Biodiverse_Inputs/cleaned_global_fishnet-may28-20.csv")

fn_rich <- fn %>% 
  group_by(X, Y) %>% 
  summarise(richness = n())

# na shapefile
na_maps <- maps::map(regions=c("usa", "mexico", "canada"), plot = FALSE, fill = TRUE) %>%
  st_as_sf()
na_but <- st_transform(na_maps, 
                       crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
na_but_buf <- st_buffer(na_but, dist = 10000)

ggplot() +
  geom_tile(fn_rich, mapping = aes(x = X, y = Y, fill = richness)) +
  geom_sf(na_but, mapping = aes(), fill = NA) + 
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) + 
  scale_fill_viridis_c(trans = "log")

## clip out the non US, MEXICO, and CA areas
fn_rich_sf <- st_as_sf(fn_rich,
                       coords = c("X", "Y"),
                       crs = st_crs(na_but))

fn_coords <- st_intersection(fn_rich_sf, na_but_buf) %>% 
  st_coordinates() %>% 
  as.data.frame()
fn_rich_clip <- st_intersection(fn_rich_sf, na_but_buf) %>% 
  st_drop_geometry()
fn_rich_clip <- dplyr::bind_cols(fn_rich_clip, fn_coords)

spp_rich <- ggplot() +
  geom_tile(fn_rich_clip, mapping = aes(x = X, y = Y, fill = richness), alpha = 0.95) +
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) + 
  scale_fill_distiller(palette = 13, na.value = "transparent", direction = 1,
                       values = scales::rescale((1:10)^4, c(0,1))) +
  labs(fill = "Richness") +
  theme_void()

###### PD AND RPD PLOTS ################
# na shapefile
na_maps <- maps::map(regions=c("usa", "mexico", "canada"), plot = FALSE, fill = TRUE) %>%
  st_as_sf()
na_but <- st_transform(na_maps, 
                       crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")

# PD 
pd <- raster("data/Biodiverse_Outputs/PD_P.tif")
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
  geom_tile(pd_clip, mapping = aes(x = X, y = Y, fill = PD_P)) +
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) + 
  scale_fill_distiller(palette = 4, na.value = "transparent", direction = 1,
                       values = scales::rescale((1:10)^4, c(0,1))) +
  labs(fill = "PD") + 
  theme_void()

## RPD Plot time
rpd <- raster("data/Biodiverse_Outputs/RPD2.tif")
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
  geom_tile(rpd_clip, mapping = aes(x = X, y = Y, fill = RPD2)) +
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) + 
  scale_fill_distiller(palette = 3, na.value = "transparent", direction = 1, 
                       values = scales::rescale((1:10)^4, c(0,1))) +
  labs(fill = "RPD") + 
  theme_void()

################ PE PLOTS ###############

# na shapefile
na_maps <- maps::map(regions=c("usa", "mexico", "canada"), plot = FALSE, fill = TRUE) %>%
  st_as_sf()
na_but <- st_transform(na_maps, 
                       crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")

# PE
pe <- raster("data/Biodiverse_Outputs/PE_WE_P.tif")
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
  geom_tile(pe_clip, mapping = aes(x = X, y = Y, fill = PE_WE_P)) +
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) + 
  scale_fill_distiller(palette = 2, na.value = "transparent", direction = 1, 
                       values = scales::rescale((1:10)^4, c(0,1))) +
  labs(fill = "PE") +
  theme_void()


## Combine plots into 4 panel plot

cp <- cowplot::plot_grid(spp_rich, pd_plot, rpd_plot, pe_plot,
                         ncol = 2, labels = c("A", "B",
                                              "C", "D"),label_y = 0.85,
                         hjust = -3, vjust = 0.5)

ggsave(filename = "figure_outputs/Figure1.png", plot = cp,
       dpi = 300, device = "png", width = 12, height = 12)
