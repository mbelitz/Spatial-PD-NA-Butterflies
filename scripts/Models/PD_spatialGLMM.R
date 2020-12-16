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
pd_df <- ras_stack_df %>% 
  dplyr::select(-c(randPD, randRPD, RPD)) 

pd_df2 <- na.omit(pd_df)

pd_top <- lm(PD ~ elev + prec + precSeas +
               temp + tempStab, 
             data = pd_df2, na.action = "na.fail")

summary(pd_top)

## Let's look at variogram of pd_top model to examine for spatial autocorrelation
jitter_dist = 50000 # 12.5 km

pd_top_resids <- residuals(pd_top)
pd_top_stand_resids <- rstandard(pd_top)

pd_df_resids <- bind_cols(pd_df2, pd_top_resids) %>% 
  bind_cols(pd_top_stand_resids) %>% 
  rename(residuals = 11, stand_residuals = 12)

pd_gdf <- as.geodata(pd_df_resids, coords.col = c("x", "y"), data.col = "residuals")
pd_gdf <-  geoR::jitterDupCoords(pd_gdf, max = jitter_dist)
pd_vario <- variog(pd_gdf)
plot(pd_vario)

ggplot(pd_df_resids) + 
  geom_tile(mapping = aes(x = x, y = y, fill = residuals)) + 
  scale_fill_viridis_c()

# Now test for Moran's I, Null hypothesis is that there is no spatial dependency
# between the residuals in the LM
pd_sims <- simulateResiduals(fittedModel = pd_top)
# test for autocorrelation
testSpatialAutocorrelation(pd_sims, x = pd_df2$x, y = pd_df2$y, plot = T)

# def showing autocorrelation. 
# We will make a spatial autoregressive model
# fit the spatial model

## Warning: This line takes a long time ~ 1 hour
pd_spamm <- fitme(PD ~ elev + prec + precSeas +
                    temp + tempStab + 
                    Matern(1 | x + y), 
                  data = pd_df2) # started 1243 # finished 1346

# summary 
summary(pd_spamm)

# Here is the model output
#
#formula: PD ~ elev + prec + precSeas + temp + tempStab + Matern(1 | x + 
#                                                                  y)
#ML: Estimation of corrPars, lambda and phi by ML.
#Estimation of fixed effects by ML.
#Estimation of lambda and phi by 'outer' ML, maximizing p_v.
#Family: gaussian ( link = identity ) 
#------------ Fixed effects (beta) ------------
#             Estimate Cond. SE t-value
#(Intercept)  0.102037 0.027325  3.7342
#elev         0.003235 0.001163  2.7811
#prec         0.010862 0.001190  9.1263
#precSeas     0.003106 0.001585  1.9599
#temp         0.014368 0.004631  3.1025
#tempStab    -0.003615 0.004737 -0.7632
#--------------- Random effects ---------------
#  Family: gaussian ( link = identity ) 
#--- Correlation parameters:
#  1.nu        1.rho 
#8.897215e-01 1.566233e-06 
#--- Variance parameters ('lambda'):
#  lambda = var(u) for u ~ Gaussian; 
#x + y  :  0.009016  
## of obs: 2012; # of groups: x + y, 2012 
#------------- Residual variance  -------------
#  phi estimate was 1.00003e-06 
#------------- Likelihood values  -------------
#  logLik
#p_v(h) (marginal L): 5251.543


pd_dd <- dist(pd_df2[,c("x","y")])
pd_mm <- MaternCorr(pd_dd, nu = pd_spamm$ranFix$corrPars$`1`$nu, 
                    rho = pd_spamm$ranFix$corrPars$`1`$rho)
plot(as.numeric(pd_dd), as.numeric(pd_mm), 
     xlab = "Distance between pairs of location [in m]", 
     ylab = "Estimated correlation")

pd_sims2 <- simulateResiduals(pd_spamm)
plot(pd_sims2)

testSpatialAutocorrelation(pd_sims2, x = pd_df2$x, pd_df2$y, plot = T)
testSpatialAutocorrelation(pd_sims, x = pd_df2$x, pd_df2$y, plot = T)



## check resids of new model
spam_pd_top_resids <- residuals(pd_spamm)

spam_pd_df_resids <- bind_cols(pd_df2,spam_pd_top_resids) %>% 
  rename(residuals = 11)

spam_pd_gdf <- as.geodata(spam_pd_df_resids, coords.col = c("x", "y"), data.col = "residuals")
spam_pd_gdf <-  geoR::jitterDupCoords(spam_pd_gdf, max = jitter_dist)
spam_pd_vario <- variog(spam_pd_gdf)
plot(spam_pd_vario)
plot(pd_vario)


ggplot(pd_df_resids) + 
  geom_tile(mapping = aes(x = x, y = y, fill = residuals)) + 
  scale_fill_viridis_c()

ggplot(spam_pd_df_resids) + 
  geom_tile(mapping = aes(x = x, y = y, fill = residuals)) + 
  scale_fill_viridis_c()

mdf <- data.frame(Model = "PD_SPAMM", Fixed_Effects = pd_spamm$fixef)
