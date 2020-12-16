library(sf)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)

#' script to generate Figure 5 which displays CANAPE results for
#' a) angiosperms
#' b) butterflies

# na shapefile
na_maps <- maps::map(regions=c("usa", "mexico", "canada"), plot = FALSE, fill = TRUE) %>%
  st_as_sf()
na_but <- st_transform(na_maps, 
                       crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")

# Butterfly Canape results
can <- readr::read_csv("data/Biodiverse_Outputs/CANAPE_results.csv") %>% 
  rename(x = 1, y = 2, Significance = 3) %>% 
  mutate(Significance = ifelse(test = Significance == "Palaeo", yes = "Paleo",
                               no = Significance))

can$Significance <- factor(can$Significance, levels = c("Paleo", "Mixed", "Neo", "Not Sig"))

but_canape_plot <- ggplot() + 
  geom_tile(can, mapping = aes(x = x, y = y, fill = Significance)) + 
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  scale_fill_manual(values = c("royalblue1", "#CB7FFF", "red", "grey93")) +
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) +
  ggtitle("B") + 
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) +
  theme_void() + 
  theme(plot.title = element_text(size = 14, hjust = 0.25)) +
  theme(legend.position = c(1.1,0.5))



# CANAPE PLANT RESULTS
xx <- raster("data/Angiosperms/exports/p_rank_PE_WE_P.tiff")
xx_df <- as.data.frame(xx, xy = TRUE)

yy <- raster("data/Angiosperms/exports/p_rank_PHYLO_RPE_NULL2.tiff")
yy_df <- as.data.frame(yy, xy = TRUE)

zz <- raster("data/Angiosperms/exports/p_rank_PHYLO_RPE2.tiff")
zz_df <- as.data.frame(zz, xy = TRUE)

significance_fun <- function(x, y, z){
  #  simplify the logic below
  ifelse(test = is.na(x), yes = 0, no = x)
  ifelse(test = is.na(y), yes = 0, no = x)
  ifelse(test = is.na(z), yes = 0.5, no = x)
  
  ifelse(test = x <= 0.95 & y <= 0.95, yes = "Not Sig",
         no = 
           ifelse(z < 0.025, yes = "Neo",
                  no = 
                    ifelse(z > 0.975, yes = "Paleo",
                           no = "Mixed")
           ))
}

Significance <- significance_fun(x = xx_df$p_rank_PE_WE_P, y = yy_df$p_rank_PHYLO_RPE_NULL2, z = zz_df$p_rank_PHYLO_RPE2)

df2 <- left_join(xx_df, yy_df) 
df2 <- left_join(df2, zz_df)

df3 <- df2 %>% 
  mutate(PE_WE_P = replace_na(p_rank_PE_WE_P, 0),
         PHYLO_RPE_NULL2 = replace_na(p_rank_PHYLO_RPE_NULL2, 0),
         PHYLO_RPE2 = replace_na(p_rank_PHYLO_RPE2, 0.5)) 

Significance <- significance_fun(x = df3$PE_WE_P, y = df3$PHYLO_RPE_NULL2,
                                 z = df3$PHYLO_RPE2)

canape_csv <- cbind(df3, Significance) %>% 
  dplyr::select(x, y, Significance)

pd_plants <- raster("data/Angiosperms/exports/results_PD.tiff")
pd_plants_df <- as.data.frame(xy = TRUE, x = pd_plants)

canape_csv <- left_join(canape_csv, pd_plants_df) %>% 
  filter(!is.na(results_PD))

canape_csv$Significance <- factor(canape_csv$Significance, levels = c("Paleo", "Mixed", "Neo", "Not Sig"))

plant_canape_plot <- ggplot() + 
  geom_tile(canape_csv, mapping = aes(x = x, y = y, fill = Significance)) + 
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  scale_fill_manual(values = c("royalblue1", "#CB7FFF", "red", "grey93")) +
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) +
  ggtitle("A") + 
  theme_void() + 
  theme(plot.title = element_text(size = 14, hjust = 0.25)) +
  theme(legend.position = c(1.1,0.5))

ces <- cowplot::plot_grid(plant_canape_plot, but_canape_plot, ncol = 2, 
                          label_y = 0.85, 
                          label_x = 0.25)
ces


ggsave(plot = ces, filename = "figure_outputs/Figure5.png",
       dpi = 450, device = "png", width = 10, height = 8)
