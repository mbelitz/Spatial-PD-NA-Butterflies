library(raster)
library(tidyverse)
library(sf)

#' script to generate Figure 5 which displays spatial residuals of linear models
#' that use plant metrics to predict butterfly metrics. The three metrics are:
#' a) PD
#' b) RPE
#' c) PE


#crop extent
na <- rnaturalearth::ne_countries(continent = "North America", returnclass = "sp")
na_but <- spTransform(na, "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")

#projection
m_proj <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

# read in butterfly response variables
pd <- raster("data/Biodiverse_Outputs/PD_P.tif")
crs(pd) <- m_proj
rpd <- raster("data/Biodiverse_Outputs/RPD2.tif")
crs(rpd) <- m_proj
pe <- raster("data/Biodiverse_Outputs/PE_WE_P.tif")
crs(pe) <- m_proj
rpe <- raster("data/Biodiverse_Outputs/CANAPE_PHYLO_RPE2.tif")
crs(rpe) <- m_proj

# read in plant response variables
pd_plants <- raster("data/Angiosperms/exports/results_PD_P.tiff")
crs(pd_plants) <- m_proj
rpd_plants <- raster("data/Angiosperms/exports/results_PHYLO_RPD2.tiff")
crs(rpd_plants) <- m_proj
pe_plants <- raster("data/Angiosperms/exports/results_PE_WE_P.tiff")
crs(pe_plants) <- m_proj
rpe_plants <- raster("data/Angiosperms/exports/results_PHYLO_RPE2.tiff")
crs(rpe_plants) <- m_proj

## CORRELATION TIME ##

# PD
pd_plants2 <- projectRaster(from = pd_plants, to = pd)

pd_stack <- stack(pd_plants2, pd)
plot(pd_stack)

cor(values(pd_stack)[,1],
    values(pd_stack)[,2],
    use = "na.or.complete",
    method = "pearson") ############ 0.581

# RPD
rpd_plants2 <- projectRaster(from = rpd_plants, to = rpd)

rpd_stack <- stack(rpd_plants2, rpd)
plot(rpd_stack)

cor(values(rpd_stack)[,1],
    values(rpd_stack)[,2],
    use = "na.or.complete",
    method = "pearson") ############ 0.0988


# PE
pe_plants2 <- projectRaster(from = pe_plants, to = pe)

pe_stack <- stack(pe_plants2, pe)
plot(pe_stack)

cor(values(pe_stack)[,1],
    values(pe_stack)[,2],
    use = "na.or.complete",
    method = "pearson") ############ 0.30

# RPE
rpe_plants2 <- projectRaster(from = rpe_plants, to = rpe)

rpe_stack <- stack(rpe_plants2, rpe)
plot(rpe_stack)

cor(values(rpe_stack)[,1],
    values(rpe_stack)[,2],
    use = "na.or.complete",
    method = "pearson") ############ - 0.002


## Regression models ##
# na shapefile
na_maps <- maps::map(regions=c("usa", "mexico", "canada"), plot = FALSE, fill = TRUE) %>%
  st_as_sf()
na_but <- st_transform(na_maps, 
                       crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

## PD
pd_stack_df <- as.data.frame(pd_stack, xy = TRUE) %>% 
  na.omit()
# model is predicting butterfly PD based on plant PD
pd_reg <- lm(PD_P ~ results_PD_P, data = pd_stack_df)
summary(pd_reg) # results_PD_P 1.757559   0.054961   31.98   <2e-16 ***
pd_reg_res <- residuals(pd_reg)

pd_stack_df2 <- bind_cols(pd_stack_df, pd_reg_res) %>% 
  rename(residuals = 5)

pd_res_plot <- ggplot() +
  geom_tile(pd_stack_df2, mapping = aes(x = x, y = y, fill = residuals)) +
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  scale_fill_gradient2(high = "#593d9cff", mid = "gray90", low  = "#f68f46ff") +
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) + 
  labs(fill = "Residuals") +
  theme_void()

## RPD
rpd_stack_df <- as.data.frame(rpd_stack, xy = TRUE) %>% 
  na.omit()

rpd_stack_df2 <- rpd_stack_df %>% 
  mutate(plant_RPD_P = results_PHYLO_RPD2 / max(results_PHYLO_RPD2)) %>% 
  mutate(but_RPD_P = RPD2 / max (RPD2)) 

rpd_reg <- lm(but_RPD_P ~ plant_RPD_P, data = rpd_stack_df2)
summary(rpd_reg) #plant_RPD_P  0.12168    0.02738   4.443 9.35e-06 ***
rpd_reg_res <- residuals(rpd_reg)

rpd_stack_df3 <- bind_cols(rpd_stack_df2, rpd_reg_res) %>% 
  rename(residuals = 7)

hist(rpd_reg_res)
min(rpd_reg_res)

rpd_res_plot <- ggplot() +
  geom_tile(rpd_stack_df3, mapping = aes(x = x, y = y, fill = residuals)) +
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  scale_fill_gradient2(high = "#593d9cff", mid = "gray90", low  = "#f68f46ff") +
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) + 
  labs(fill = "Residuals") +
  theme_void()

## PE
pe_stack_df <- as.data.frame(pe_stack, xy = TRUE) %>% 
  na.omit()

pe_reg <- lm(PE_WE_P ~ results_PE_WE_P, data = pe_stack_df)
summary(pe_reg) # results_PE_WE_P 2.352e-01  1.617e-02   14.55   <2e-16 ***
pe_reg_res <- residuals(pe_reg)

pe_stack_df2 <- bind_cols(pe_stack_df, pe_reg_res) %>% 
  rename(residuals = 5)

hist(pe_reg_res)
min(pe_reg_res)

pe_res_plot <- ggplot() +
  geom_tile(pe_stack_df2, mapping = aes(x = x, y = y, fill = residuals)) +
  geom_sf(na_but, mapping = aes(), fill = NA, color = "grey25") + 
  scale_fill_gradient2(high = "#593d9cff", mid = "gray90", low  = "#f68f46ff") +
  coord_sf(xlim = c(-4850000, 3050000), ylim = c(-2950000, 4550000)) + 
  labs(fill = "Residuals") +
  theme_void()

#combine plots
reg_cps <- cowplot::plot_grid(pd_res_plot, rpd_res_plot, pe_res_plot,
                              ncol = 1, labels = c("A", "B", "C"),
                              label_x = 0.1)

ggsave(plot = reg_cps, filename = "Figure_Outputs/Figure5.png",
       width = 6, height = 8)