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
## RPD RAND
## RPD Rand
rpd_rand_df <- ras_stack_df %>% 
  dplyr::select(-c(RPD, PD, randPD))

rpd_rand_df2 <- filter(rpd_rand_df, !is.na(randRPD))
rpd_rand_df2 <- na.omit(rpd_rand_df2) %>% 
  mutate(Significance =  ifelse(test = randRPD <= 0.05, yes = 0, no = 1))

rpd_rand_top <- glm(Significance ~ elev  + prec + temp + tempStab,
                    data = rpd_rand_df2, na.action = "na.fail",
                    family = binomial(link = "logit"))

summary(rpd_rand_top)


## Let's look at variogram of pd_top model to examine for spatial autocorrelation
jitter_dist = 50000 # 12.5 km

rpd_rand_top_resids <- residuals(rpd_rand_top)
rpd_rand_top_stand_resids <- rstandard(rpd_rand_top)

rpd_rand_df_resids <- bind_cols(rpd_rand_df2, rpd_rand_top_resids) %>% 
  bind_cols(rpd_rand_top_stand_resids) %>% 
  rename(residuals = 12, stand_residuals = 13)

rpd_rand_gdf <- as.geodata(rpd_rand_df_resids, coords.col = c("x", "y"), data.col = "residuals")
rpd_rand_gdf <-  geoR::jitterDupCoords(rpd_rand_gdf, max = jitter_dist)
rpd_rand_vario <- variog(rpd_rand_gdf)
plot(rpd_rand_vario)

ggplot(rpd_rand_df_resids) + 
  geom_tile(mapping = aes(x = x, y = y, fill = residuals)) + 
  scale_fill_viridis_c()

# Now test for Moran's I, Null hypothesis is that there is no spatial dependency
# test for autocorrelation
testSpatialAutocorrelation(rpd_rand_top_resids, x = rpd_rand_df2$x, y = rpd_rand_df2$y, plot = F)

# def showing autocorrelation. 
# We will make a spatial autoregressive model
# fit the spatial model

## Warning this takes a very long time to fit -- I ran it overnight

rpd_rand_spamm <- fitme(Significance ~ elev  + prec + temp + tempStab + 
                          Matern(1 | x + y), 
                        data = rpd_rand_df2,
                        family = binomial(link = "logit"))

summary(rpd_rand_spamm)

## Model Output

# Summary
#formula: Significance ~ elev + prec + temp + tempStab + Matern(1 | x + 
#                                                                 y)
#Estimation of corrPars and lambda by Laplace ML approximation (p_v).
#Estimation of fixed effects by Laplace ML approximation (p_v).
#Estimation of lambda by 'outer' ML, maximizing p_v.
#Family: binomial ( link = logit ) 
#------------ Fixed effects (beta) ------------
#  Estimate Cond. SE t-value
#(Intercept)  -10.375    84.27 -0.1231
#elev           7.601    29.54  0.2573
#prec           9.956    28.89  0.3447
#temp          57.056   133.24  0.4282
#tempStab     -11.297    72.40 -0.1560
#--------------- Random effects ---------------
#  Family: gaussian ( link = identity ) 
#--- Correlation parameters:
#  1.nu        1.rho 
#1.604900e+01 4.774708e-06 
#--- Variance parameters ('lambda'):
#  lambda = var(u) for u ~ Gaussian; 
#x + y  :  14700  
## of obs: 613; # of groups: x + y, 613 
#------------- Likelihood values  -------------
#  logLik
#p_v(h) (marginal L): -10.82063

rpd_rand_dd <- dist(rpd_rand_df2[,c("x","y")])
rpd_rand_mm <- MaternCorr(rpd_rand_dd, nu = rpd_rand_spamm$ranFix$corrPars$`1`$nu, 
                          rho = rpd_rand_spamm$ranFix$corrPars$`1`$rho)
plot(as.numeric(rpd_rand_dd), as.numeric(rpd_rand_mm), 
     xlab = "Distance between pairs of location [in m]", 
     ylab = "Estimated correlation")

rpd_rand_sims2 <- simulateResiduals(rpd_rand_spamm)
plot(rpd_rand_sims2)

testSpatialAutocorrelation(rpd_rand_sims2, x = rpd_rand_df2$x, rpd_rand_df2$y, plot = T)
testSpatialAutocorrelation(rpd_rand_top_resids, x = rpd_rand_df2$x, y = rpd_rand_df2$y, plot = F)


## check resids of new model
spam_rpd_rand_top_resids <- residuals(rpd_rand_spamm)

spam_rpd_rand_df_resids <- bind_cols(rpd_rand_df2,spam_rpd_rand_top_resids) %>% 
  rename(residuals = 12)

spam_rpd_rand_gdf <- as.geodata(spam_rpd_rand_df_resids, coords.col = c("x", "y"), data.col = "residuals")
spam_rpd_rand_gdf <-  geoR::jitterDupCoords(spam_rpd_rand_gdf, max = jitter_dist)
spam_rpd_rand_vario <- variog(spam_rpd_rand_gdf)
plot(spam_rpd_rand_vario)

ggplot(rpd_rand_df_resids) + 
  geom_tile(mapping = aes(x = x, y = y, fill = residuals)) + 
  scale_fill_viridis_c()

ggplot(spam_rpd_rand_df_resids) + 
  geom_tile(mapping = aes(x = x, y = y, fill = residuals)) + 
  scale_fill_viridis_c()

mdf <- data.frame(Model = "RPD_rand_SPAMM", Fixed_Effects = rpd_rand_spamm$fixef)
