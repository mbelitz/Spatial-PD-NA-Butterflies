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
#### WARNING ----- Spatial logistic regression Models take a long time to fit ----- 
####~ 6 hours for this script using a 16 ram 4 core machine

#projection
m_proj <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

# read in butterfly response variables
pd <- raster("data/Biodiverse_Outputs/PD_P.tif")
crs(pd) <- m_proj
pd <- raster("data/Biodiverse_Outputs/pd2.tif")
crs(pd) <- m_proj


pd_rand <- raster("data/Biodiverse_Outputs/randomization_p_rank_PD_P.tif")
crs(pd_rand) <- m_proj
pd_rand <- raster("data/Biodiverse_Outputs/randomization_p_rank_pd2.tif")
crs(pd_rand) <- m_proj

## Read in model covariates 
temp <- raster("data/Model_Covariates/temp.tif")
prec <- raster("data/Model_Covariates/prec.tif")
tempSeas <- raster("data/Model_Covariates/tempSeas.tif")
precSeas <- raster("data/Model_Covariates/precSeas.tif")
tempStab <- raster("data/Model_Covariates/tempStab.tif")
precStab <- raster("data/Model_Covariates/precStab.tif")
elev <- raster("data/Model_Covariates/elev.tif")

# stack rasters
ras_stack <- stack(pd, pd, pd_rand, pd_rand,
                   temp, prec, tempSeas, precSeas, 
                   tempStab, precStab,
                   elev)

plot(ras_stack)
## Standardize the data to have mean of 0 and SD of 1
ras_stack_df <- as.data.frame(ras_stack, xy = TRUE)

ras_stack_df <- ras_stack_df %>% 
  dplyr::rename(PD = 3, pd = 4, randPD = 5, randpd = 6) %>% 
  mutate(temp = scale (temp),
         prec = scale(prec),
         precStab = scale(precStab),
         tempStab = scale(tempStab),
         tempSeas = scale(tempSeas),
         precSeas =scale(precSeas),
         elev = scale(elev))

## Start with PD Linear Models
## PD Rand
pd_rand_df <- ras_stack_df %>% 
  dplyr::select(-c(pd, PD, randpd))

pd_rand_df2 <- filter(pd_rand_df, !is.na(randPD))
pd_rand_df2 <- na.omit(pd_rand_df2) %>% 
  mutate(Significance =  ifelse(test = randPD <= 0.05, yes = 0, no = 1))

pd_rand_top <- glm(Significance ~ elev + prec + temp + tempStab, 
                   data = pd_rand_df2, 
                   family = binomial(link = "logit"),
                   na.action = "na.fail")

summary(pd_rand_top)

## Let's look at variogram of pd_top model to examine for spatial autocorrelation
jitter_dist = 50000 # 12.5 km


pd_rand_top_resids <- residuals(pd_rand_top)
pd_rand_top_stand_resids <- rstandard(pd_rand_top)

pd_rand_df_resids <- bind_cols(pd_rand_df2, pd_rand_top_resids) %>% 
  bind_cols(pd_rand_top_stand_resids) %>% 
  rename(residuals = 12, stand_residuals = 13)
pd_rand_gdf <- as.geodata(pd_rand_df_resids, coords.col = c("x", "y"), data.col = "residuals")
pd_rand_gdf <-  geoR::jitterDupCoords(pd_rand_gdf, max = jitter_dist)
pd_rand_vario <- variog(pd_rand_gdf)
plot(pd_rand_vario)

ggplot(pd_rand_df_resids) + 
  geom_tile(mapping = aes(x = x, y = y, fill = residuals)) + 
  scale_fill_viridis_c()


# Now test for Moran's I, Null hypothesis is that there is no spatial dependency
# between the residuals in the LM
# test for autocorrelation
testSpatialAutocorrelation(pd_rand_top_resids, x = pd_rand_df2$x, y = pd_rand_df2$y, plot = F)

# def showing autocorrelation. 
# We will make a spatial autoregressive model
# fit the spatial model

## Warning this takes a very long time 
## Ran overnight

pd_rand_spamm <- fitme(Significance ~ elev  + prec + temp + tempStab + 
                         Matern(1 | x + y), 
                       data = pd_rand_df2,
                       family = binomial(link = "logit"))

## Model Output
#' 
#' formula: Significance ~ elev + prec + temp + tempStab + Matern(1 | x + 
# y)
# Estimation of corrPars and lambda by Laplace ML approximation (p_v).
# Estimation of fixed effects by Laplace ML approximation (p_v).
# Estimation of lambda by 'outer' ML, maximizing p_v.
# Family: binomial ( link = logit ) 
# ------------ Fixed effects (beta) ------------
#               Estimate Cond. SE t-value
# (Intercept)  -51.869    42.06 -1.2332
# elev           1.552    11.79  0.1316
# prec           2.981    18.82  0.1584
# temp          53.185    30.30  1.7553
# tempStab      -7.414    21.64 -0.3426
# --------------- Random effects ---------------
#   Family: gaussian ( link = identity ) 
# --- Correlation parameters:
#   1.nu        1.rho 
# 1.086866e+00 2.264443e-06 
# --- Variance parameters ('lambda'):
#   lambda = var(u) for u ~ Gaussian; 
# x + y  :  2439  
# # of obs: 1336; # of groups: x + y, 1336 
# ------------- Likelihood values  -------------
#   logLik
# p_v(h) (marginal L): -14.70234


pd_rand_dd <- dist(pd_rand_df2[,c("x","y")])
pd_rand_mm <- MaternCorr(pd_rand_dd, nu = pd_rand_spamm$ranFix$corrPars$`1`$nu, 
                         rho = pd_rand_spamm$ranFix$corrPars$`1`$rho)
plot(as.numeric(pd_rand_dd), as.numeric(pd_rand_mm), 
     xlab = "Distance between pairs of location [in m]", 
     ylab = "Estimated correlation")

pd_rand_sims2 <- simulateResiduals(pd_rand_spamm)
plot(pd_rand_sims2)

testSpatialAutocorrelation(pd_rand_sims2, x = pd_rand_df2$x, pd_rand_df2$y, plot = T)

## check resids of new model
spam_pd_rand_top_resids <- residuals(pd_rand_spamm)

spam_pd_rand_df_resids <- bind_cols(pd_rand_df2,spam_pd_rand_top_resids) %>% 
  rename(residuals = 12)

spam_pd_rand_gdf <- as.geodata(spam_pd_rand_df_resids, coords.col = c("x", "y"), data.col = "residuals")
spam_pd_rand_gdf <-  geoR::jitterDupCoords(spam_pd_rand_gdf, max = jitter_dist)
spam_pd_rand_vario <- variog(spam_pd_rand_gdf)
plot(spam_pd_rand_vario)

ggplot(pd_df_resids) + 
  geom_tile(mapping = aes(x = x, y = y, fill = residuals)) + 
  scale_fill_viridis_c()

ggplot(spam_pd_rand_df_resids) + 
  geom_tile(mapping = aes(x = x, y = y, fill = residuals)) + 
  scale_fill_viridis_c()

mdf <- data.frame(Model = "pd_rand_SPAMM", Fixed_Effects = pd_rand_spamm$fixef)
