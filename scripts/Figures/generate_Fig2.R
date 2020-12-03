library(raster)
library(sf)
library(dplyr)
library(ggplot2)
library(scales)
library(gridExtra)

#' script to generate Figure 2 which displays PD randomizations results for
#' a) angiosperms
#' b) butterflies

# na shapefile
na_maps <- maps::map(regions=c("usa", "mexico", "canada"), plot = FALSE, fill = TRUE) %>%
  st_as_sf()
na_but <- st_transform(na_maps, 
                       crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")

# READ in plant pd data and butterfly pd data
plant_pd <- raster("data/Angiosperms/exports/results_PD_P.tiff")
crs(plant_pd) <- crs(na_but)
but_pd <- raster("data/Biodiverse_Outputs/PD_P.tif")
crs(but_pd) <- crs(na_but)
but_pd_df <- as.data.frame(but_pd, xy = TRUE)
plant_pd_df <- as.data.frame(plant_pd, xy = TRUE)

plant_rand_pd <- raster("data/Angiosperms/exports/p_rank_PD_P.tiff")
crs(plant_rand_pd) <- crs(na_but)
but_rand_pd <- raster("data/Biodiverse_Outputs/randomization_p_rank_PD_P.tif")
crs(but_rand_pd) <- crs(na_but)
but_rand_pd_df <- as.data.frame(but_rand_pd, xy = TRUE)
plant_rand_pd_df <- as.data.frame(plant_rand_pd, xy = TRUE)

plant_rand_pd2 <- left_join(plant_rand_pd_df, plant_pd_df)
but_rand_pd2 <- dplyr::left_join(but_rand_pd_df, but_pd_df)

plant_rand_pdc <- plant_rand_pd2 %>% 
  mutate(Significance = case_when(p_rank_PD_P <= 0.05 ~ "Low",
                                  p_rank_PD_P > 0.05 & p_rank_PD_P < 0.95 ~ "Non",
                                  p_rank_PD_P >= 0.95 ~ "High")) %>% 
  filter(!is.na(results_PD_P))

plant_rand_pdc2 <- plant_rand_pdc %>% 
  mutate(Significance = ifelse(is.na(Significance),yes = "Non", no = Significance))

but_rand_pdc <- but_rand_pd2 %>% 
  mutate(Significance = case_when(randomization_p_rank_PD_P <= 0.05 ~ "Low",
                                  randomization_p_rank_PD_P > 0.05 & randomization_p_rank_PD_P < 0.95 ~ "Non",
                                  randomization_p_rank_PD_P >= 0.95 ~ "High")) %>% 
  filter(!is.na(PD_P)) 


but_rand_pdc2 <- but_rand_pdc %>% 
  mutate(Significance = ifelse(is.na(Significance),yes = "Non", no = Significance))

###### Plot plant PD randomizations
plant_pd_plot2 <- ggplot() +
  geom_tile(plant_rand_pdc2, mapping = aes(x = x, y = y, fill = factor(Significance, levels = c("High", "Non", "Low")))) + 
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  scale_fill_manual(values = c("dodgerblue2", "grey93", "red1")) +
  labs(fill = "Significance") +
  ggtitle("A") + 
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) +
  theme_void() + 
  theme(plot.title = element_text(size = 14, hjust = 0.25)) +
  theme(legend.position = c(1.1,0.5))

###### Plot bfly PD randomizations
but_pd_plot2 <- ggplot() +
  geom_tile(but_rand_pdc2, mapping = aes(x = x, y = y, fill = factor(Significance, levels = c("High", "Non", "Low")))) + 
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  scale_fill_manual(values = c("dodgerblue2", "grey93", "red1")) +
  labs(fill = "Significance") +
  ggtitle("B") + 
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) +
  theme_void() + 
  theme(plot.title = element_text(size = 14, hjust = 0.25)) +
  theme(legend.position = c(1.1,0.5))

pds2 <- cowplot::plot_grid(plant_pd_plot2, but_pd_plot2, ncol = 2)

ggsave(plot = pds2, filename = "Figure_Outputs/Figure2.png",
       dpi = 450, device = "png", width = 10, height = 8)
