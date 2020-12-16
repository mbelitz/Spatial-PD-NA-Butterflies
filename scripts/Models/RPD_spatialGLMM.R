library(raster)
library(MuMIn)
library(dplyr)
library(car)
library(sf)
library(ggplot2)
library(geoR)
library(spdep)
library(DHARMa)
library(spaMM)

#' Script to run spatial model with spatially correlated random effects
#' 
#' First we check for spatial autocorrelation using Moran's I 
#' 
#' Then we run Spatial Model
#' 
#### WARNING ----- Spatial Models take awhile to fit ----- ~ 1 hour 
#### for this script using a 16 ram 4 core machine


#projection
m_proj <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

# read in butterfly response variables
pd <- raster("data/Biodiverse_Outputs/PD_P.tif")
crs(pd) <- m_proj
rpd <- raster("data/Biodiverse_Outputs/RPD2.tif")
crs(rpd) <- m_proj


pd_rand <- raster("data/Biodiverse_Outputs/randomization_p_rank_PD_P.tif")
crs(pd_rand) <- m_proj
rpd_rand <- raster("data/Biodiverse_Outputs/randomization_p_rank_RPD2.tif")
crs(rpd_rand) <- m_proj

## Read in model covariates 
temp <- raster("data/Model_Covariates/temp.tif")
prec <- raster("data/Model_Covariates/prec.tif")
tempSeas <- raster("data/Model_Covariates/tempSeas.tif")
precSeas <- raster("data/Model_Covariates/precSeas.tif")
tempStab <- raster("data/Model_Covariates/tempStab.tif")
precStab <- raster("data/Model_Covariates/precStab.tif")
elev <- raster("data/Model_Covariates/elev.tif")

# stack rasters
ras_stack <- stack(pd, rpd, pd_rand, rpd_rand,
                   temp, prec, tempSeas, precSeas, 
                   tempStab, precStab,
                   elev)

plot(ras_stack)
## Standardize the data to have mean of 0 and SD of 1
ras_stack_df <- as.data.frame(ras_stack, xy = TRUE)

ras_stack_df <- ras_stack_df %>% 
  dplyr::rename(PD = 3, RPD = 4, randPD = 5, randRPD = 6) %>% 
  mutate(temp = scale (temp),
         prec = scale(prec),
         precStab = scale(precStab),
         tempStab = scale(tempStab),
         tempSeas = scale(tempSeas),
         precSeas =scale(precSeas),
         elev = scale(elev))

## Start with PD Linear Models

## RPD
rpd_df <- ras_stack_df %>% 
  dplyr::select(-c(randPD, randRPD, PD)) 

rpd_df2 <- na.omit(rpd_df)

rpd_top <- lm(RPD ~ elev + prec + precSeas + 
                precStab + tempSeas + tempStab,
              data = rpd_df2, na.action = "na.fail")


summary(rpd_top)


## Let's look at variogram of pd_top model to examine for spatial autocorrelation
jitter_dist = 50000 # 12.5 km

rpd_top_resids <- residuals(rpd_top)
rpd_top_stand_resids <- rstandard(rpd_top)

rpd_df_resids <- bind_cols(rpd_df2, rpd_top_resids) %>% 
  bind_cols(rpd_top_stand_resids) %>% 
  rename(residuals = 11, stand_residuals = 12)

rpd_gdf <- as.geodata(rpd_df_resids, coords.col = c("x", "y"), data.col = "residuals")
rpd_gdf <-  geoR::jitterDupCoords(rpd_gdf, max = jitter_dist)
rpd_vario <- variog(rpd_gdf)
plot(rpd_vario)

ggplot(rpd_df_resids) + 
  geom_tile(mapping = aes(x = x, y = y, fill = residuals)) + 
  scale_fill_viridis_c()

# test for autocorrelation
testSpatialAutocorrelation(rpd_top_resids, x = rpd_df2$x, y = rpd_df2$y, plot = F)

# def showing autocorrelation. 
# We will make a spatial autoregressive model
# fit the spatial model

## Warning this step can take ~ 1 hour to fit

rpd_spamm <- fitme(RPD ~ elev + prec + precSeas + 
                     precStab + tempSeas + tempStab +
                     Matern(1 | x + y), 
                   data = rpd_df2)

summary(rpd_spamm)

## Model Output

## formula: RPD ~ elev + prec + precSeas + precStab + tempSeas + tempStab + 
##   Matern(1 | x + y)
## ML: Estimation of corrPars, lambda and phi by ML.
## Estimation of fixed effects by ML.
## Estimation of lambda and phi by 'outer' ML, maximizing p_v.
## Family: gaussian ( link = identity ) 
## ------------ Fixed effects (beta) ------------
##   Estimate Cond. SE t-value
## (Intercept)  0.916784 0.070155 13.0680
## elev        -0.002439 0.001175 -2.0761
## prec         0.004600 0.001573  2.9240
## precSeas    -0.002061 0.002073 -0.9941
## precStab    -0.001289 0.002847 -0.4527
## tempSeas     0.002998 0.003611  0.8303
## tempStab    -0.010336 0.005048 -2.0474
## --------------- Random effects ---------------
##   Family: gaussian ( link = identity ) 
## --- Correlation parameters:
##   1.nu        1.rho 
## 5.634083e-01 4.055155e-07 
## --- Variance parameters ('lambda'):
##   lambda = var(u) for u ~ Gaussian; 
## x + y  :  0.01703  
## # of obs: 2012; # of groups: x + y, 2012 
## ------------- Residual variance  -------------
##   phi estimate was 1.41465e-06 
## ------------- Likelihood values  -------------
##   logLik
## p_v(h) (marginal L): 4727.728

rpd_dd <- dist(rpd_df2[,c("x","y")])
rpd_mm <- MaternCorr(rpd_dd, nu = rpd_spamm$ranFix$corrPars$`1`$nu, 
                     rho = rpd_spamm$ranFix$corrPars$`1`$rho)
plot(as.numeric(rpd_dd), as.numeric(rpd_mm), 
     xlab = "Distance between pairs of location [in m]", 
     ylab = "Estimated correlation")

rpd_sims2 <- simulateResiduals(rpd_spamm)
plot(rpd_sims2)

testSpatialAutocorrelation(rpd_sims2, x = rpd_df2$x, rpd_df2$y, plot = T)
testSpatialAutocorrelation(rpd_top_resids, x = rpd_df2$x, y = rpd_df2$y, plot = F)


## check resids of new model
spam_rpd_top_resids <- residuals(rpd_spamm)

spam_rpd_df_resids <- bind_cols(rpd_df2,spam_rpd_top_resids) %>% 
  rename(residuals = 11)

spam_rpd_gdf <- as.geodata(spam_rpd_df_resids, coords.col = c("x", "y"), data.col = "residuals")
spam_rpd_gdf <-  geoR::jitterDupCoords(spam_rpd_gdf, max = jitter_dist)
spam_rpd_vario <- variog(spam_rpd_gdf)
plot(spam_rpd_vario)
plot(rpd_vario)

ggplot(rpd_df_resids) + 
  geom_tile(mapping = aes(x = x, y = y, fill = residuals)) + 
  scale_fill_viridis_c()

ggplot(spam_rpd_df_resids) + 
  geom_tile(mapping = aes(x = x, y = y, fill = residuals)) + 
  scale_fill_viridis_c()

mdf <- data.frame(Model = "RPD_SPAMM", Fixed_Effects = rpd_spamm$fixef)