library(raster)
library(MuMIn)
library(dplyr)
library(knitr)
library(car)
library(sf)
library(gt)

#' Script to run models that examine what the environmental drivers are for
#' 1) PD, 2)RPD, 3)PD-Rand, and 4)RPD-Rand

## KEY NOTES OF MODELS ##
#' Temp and Temp Seasonality are correlaed. Only one can be in top model

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
ras_stack_df <- data.frame(values(ras_stack))

ras_stack_df <- ras_stack_df %>% 
  dplyr::rename(PD = 1, RPD = 2, randPD = 3, randRPD = 4) %>% 
  mutate(temp = scale (temp),
         prec = scale(prec),
         precStab = scale(precStab),
         tempStab = scale(tempStab),
         tempSeas = scale(tempSeas),
         precSeas =scale(precSeas),
         elev = scale(elev))

## Start with PD
pd_df <- ras_stack_df %>% 
  dplyr::select(-c(randPD, randRPD, RPD)) 

pd_df2 <- na.omit(pd_df)

pd_global <- lm(PD ~ temp + prec + tempSeas +  precSeas +
                  tempStab + precStab +
                  elev, 
                data = pd_df2, na.action = "na.fail")

pd_dd <- dredge(pd_global)

pd_top <- lm(PD ~ elev + prec + precSeas +
               temp + tempStab, 
             data = pd_df2, na.action = "na.fail")

vif(pd_top)

pd_second <- lm(PD ~ prec + precSeas + precStab +   
                  temp + tempStab,
                data = pd_df2, na.action = "na.fail")

vif(pd_second)


AIC(pd_top, pd_second)
MuMIn::Weights(AIC(pd_top, pd_second))
car::vif(pd_top)
summary(pd_top)

## RPD
rpd_df <- ras_stack_df %>% 
  dplyr::select(-c(randPD, randRPD, PD)) 

rpd_df2 <- na.omit(rpd_df)

rpd_global <- lm(RPD ~ temp + prec + tempSeas +  precSeas +
                   tempStab + precStab +
                   elev, 
                 data = rpd_df2, na.action = "na.fail")

rpd_dd <- dredge(rpd_global)

rpd_top <- lm(RPD ~ elev + prec + precSeas + 
                precStab + temp + tempSeas, 
              data = rpd_df2, na.action = "na.fail")

car::vif(rpd_top) ## temp and tempSeas are correlated so need to remove one of those

rpd_new_top <- lm(RPD ~ elev + prec + precSeas + 
                    precStab + tempSeas + tempStab,
                  data = rpd_df2, na.action = "na.fail")

car::vif(rpd_new_top)

rpd_second <- lm(RPD ~ elev + prec + precSeas + 
                   precStab + tempSeas, 
                 data = rpd_df2, na.action = "na.fail")
car::vif(rpd_second)

AIC(rpd_new_top, rpd_second)
MuMIn::Weights(AIC(rpd_new_top, rpd_second))
summary(rpd_new_top)

## PD Rand
pd_rand_df <- ras_stack_df %>% 
  dplyr::select(-c(RPD, PD, randRPD))

pd_rand_df2 <- filter(pd_rand_df, !is.na(randPD))
pd_rand_df2 <- na.omit(pd_rand_df2) %>% 
  mutate(Significance =  ifelse(test = randPD <= 0.05, yes = 0, no = 1))

pd_rand_global <- lm(randPD ~ temp + prec + tempSeas +  precSeas +
                       tempStab + precStab +
                       elev, 
                     data = pd_rand_df2, na.action = "na.fail")

pd_rand_dd <- dredge(pd_rand_global)

car::vif(pd_rand_global) ## temp and tempSeas are correlated

pd_rand_top <- lm(randPD ~ elev + prec + precSeas +
                    precStab + temp + tempStab, 
                  data = pd_rand_df2, na.action = "na.fail")

car::vif(pd_rand_top)

pd_rand_second <- lm(randPD ~ elev + prec + precSeas +
                       temp + tempStab, 
                     data = pd_rand_df2, na.action = "na.fail")

car::vif(pd_rand_second)

AIC(pd_rand_top, pd_rand_second)
MuMIn::Weights(AIC(pd_rand_top, pd_rand_second))
summary(pd_rand_top)


## RPD Rand
rpd_rand_df <- ras_stack_df %>% 
  dplyr::select(-c(RPD, PD, randPD))

rpd_rand_df2 <- filter(rpd_rand_df, !is.na(randRPD))
rpd_rand_df2 <- na.omit(rpd_rand_df2) %>% 
  mutate(Significance =  ifelse(test = randRPD <= 0.05, yes = 0, no = 1))

rpd_rand_global <- lm(randRPD ~ temp + prec + tempSeas +  precSeas +
                        tempStab + precStab +
                        elev, 
                      data = rpd_rand_df2, na.action = "na.fail")

rpd_rand_dd <- dredge(rpd_rand_global)

car::vif(rpd_rand_global) ## temp and tempSeas are correlated

rpd_rand_top <- lm(randRPD ~ precSeas + precStab + temp,
                   data = rpd_rand_df2, na.action = "na.fail")

car::vif(rpd_rand_top)

rpd_rand_second <- lm(randRPD ~ elev + precSeas + temp,
                      data = rpd_rand_df2, na.action = "na.fail")

car::vif(rpd_rand_second)

AIC(rpd_rand_top, rpd_rand_second)
MuMIn::Weights(AIC(rpd_rand_top, rpd_rand_second))
summary(rpd_rand_top)

## Prepare table in a pretty way
Covariates <- c("Temp", "Prec", "Temp Seas", "Prec Seas", "Temp Stab", 
                "Precip Stab", "Elev", "Delta", "Weight")

## PD top model PD ~ temp + prec + NA + precSeas + Temp Stab + NA + elev

PD <- c(pd_top$coefficients['temp'], pd_top$coefficients['prec'],
        NA, pd_top$coefficients['precSeas'], pd_top$coefficients["tempStab"], NA,
        pd_top$coefficients["elev"], -29.75, 1)

RPD <- c(NA, rpd_new_top$coefficients['prec'], rpd_new_top$coefficients["tempSeas"],
         rpd_new_top$coefficients['precSeas'], rpd_new_top$coefficients["tempStab"],
         rpd_new_top$coefficients["precStab"], rpd_new_top$coefficients["elev"], -7.92, 0.981)

PD_Sig <- c(pd_rand_top$coefficients["temp"], pd_rand_top$coefficients["prec"],
            NA, pd_rand_top$coefficients["precSeas"], pd_rand_top$coefficients["tempStab"],
            pd_rand_top$coefficients["precStab"], pd_rand_top$coefficients["elev"], -12.50, 0.998)

RPD_Sig <- c(rpd_rand_top$coefficients["temp"], NA, NA, rpd_rand_top$coefficients["precSeas"],
             NA, rpd_rand_top$coefficients["precStab"], NA, -7.81, 0.98)

table_df <- data.frame(Covariates, PD, RPD, PD_Sig, RPD_Sig)

new_tab <- table_df %>% tidyr::pivot_longer(-Covariates, names_to = "Model", values_to = "count") %>%            
  tidyr::pivot_wider(names_from = "Covariates", values_from = "count") 


new_tab2 <- new_tab %>% 
  dplyr::rename("Model" = 1, "Temp" = 2, "Prec" = 4, 
                "Temp Seas" = 4, "Prec Seas" = 5, 
                "Temp Stab" = 6, "Prec Stab" = 7, "Elev" = 8, 
                "Delta" = 9, "Weight" = 10)

new_tab2


gt_tab <- new_tab2 %>% gt() %>% 
  fmt_number(columns = c("Temp", "Prec", "Temp Seas", "Prec Seas", "Temp Stab",
                         "Prec Stab", "Elev"), decimals = 3) %>% 
  fmt_number(columns = c("Weight", "Delta"), decimals = 2) %>% 
  cols_align(align = "center")

gt_tab
