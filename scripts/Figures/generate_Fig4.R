library(raster)
library(sf)
library(dplyr)
library(ggplot2)
library(scales)
library(gridExtra)

#' script to generate Figure 4 which displays RPD randomizations results for
#' a) angiosperms
#' b) butterflies

# na shapefile
na_maps <- maps::map(regions=c("usa", "mexico", "canada"), plot = FALSE, fill = TRUE) %>%
  st_as_sf()
na_but <- st_transform(na_maps, 
                       crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")

# READ in plant RPD data and butterfly RPD data
plant_rpd <- raster("data/Angiosperms/exports/results_PHYLO_RPD2.tiff")
crs(plant_rpd) <- crs(na_but)
but_rpd <- raster("data/Biodiverse_Outputs/RPD2.tif")
crs(but_rpd) <- crs(na_but)
but_rpd_df <- as.data.frame(but_rpd, xy = TRUE)
plant_rpd_df <- as.data.frame(plant_rpd, xy = TRUE)

plant_rand_rpd <- raster("data/Angiosperms/exports/p_rank_PHYLO_RPD2.tiff")
crs(plant_rand_rpd) <- crs(na_but)
but_rand_rpd <- raster("data/Biodiverse_Outputs/randomization_p_rank_RPD2.tif")
crs(but_rand_rpd) <- crs(na_but)
but_rand_rpd_df <- as.data.frame(but_rand_rpd, xy = TRUE)
plant_rand_rpd_df <- as.data.frame(plant_rand_rpd, xy = TRUE)

plant_rand_rpd2 <- left_join(plant_rand_rpd_df, plant_rpd_df)
but_rand_rpd2 <- dplyr::left_join(but_rand_rpd_df, but_rpd_df)

plant_rand_rpdc <- plant_rand_rpd2 %>% 
  mutate(Significance = case_when(p_rank_PHYLO_RPD2 <= 0.05 ~ "Low",
                                  p_rank_PHYLO_RPD2 > 0.05 & p_rank_PHYLO_RPD2 < 0.95 ~ "Non",
                                  p_rank_PHYLO_RPD2 >= 0.95 ~ "High")) %>% 
  filter(!is.na(results_PHYLO_RPD2))

plant_rand_rpdc2 <- plant_rand_rpdc %>% 
  mutate(Significance = ifelse(is.na(Significance),yes = "Non", no = Significance))


plant_rpd_plot <- ggplot() +
  geom_tile(plant_rand_rpdc2, mapping = aes(x = x, y = y, fill = factor(Significance, levels = c("High", "Non", "Low")))) + 
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  scale_fill_manual(values = c( "turquoise4", "grey93","violetred2")) +
  labs(fill = "Significance") +
  ggtitle("A") + 
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) +
  theme_void() + 
  theme(plot.title = element_text(size = 14, hjust = 0.25)) +
  theme(legend.position = c(1.1,0.5))

but_rand_rpdc <- but_rand_rpd2 %>% 
  mutate(Significance = case_when(randomization_p_rank_RPD2 <= 0.05 ~ "Low",
                                  randomization_p_rank_RPD2 > 0.05 & randomization_p_rank_RPD2 < 0.95 ~ "Non",
                                  randomization_p_rank_RPD2 >= 0.95 ~ "High")) %>% 
  filter(!is.na(RPD2)) 


but_rand_rpdc2 <- but_rand_rpdc %>% 
  mutate(Significance = ifelse(is.na(Significance),yes = "Non", no = Significance))

#plant plot
plant_rpd_plot2 <- ggplot() +
  geom_tile(plant_rand_rpdc2, mapping = aes(x = x, y = y, fill = factor(Significance, levels = c("High", "Non", "Low")))) + 
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  scale_fill_manual(values = c("dodgerblue2", "grey93", "red1")) +
  labs(fill = "Significance") +
  ggtitle("A") + 
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) +
  theme_void() + 
  theme(plot.title = element_text(size = 14, hjust = 0.25)) +
  theme(legend.position = c(1.1,0.5))

# butterfly plot
but_rpd_plot2 <- ggplot() +
  geom_tile(but_rand_rpdc2, mapping = aes(x = x, y = y, fill = factor(Significance, levels = c("High", "Non", "Low")))) + 
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  scale_fill_manual(values = c("dodgerblue2", "grey93", "red1")) +
  labs(fill = "Significance") +
  ggtitle("B") + 
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) +
  theme_void() + 
  theme(plot.title = element_text(size = 14, hjust = 0.25)) +
  theme(legend.position = c(1.1,0.5))

rpds2 <- cowplot::plot_grid(plant_rpd_plot2, but_rpd_plot2, ncol = 2)

ggsave(plot = rpds2, filename = "figure_outputs/Figure4.png",
       dpi = 450, device = "png", width = 10, height = 8)
