#====Load Packages==========
library(pspearman)
library(dplyr)
library(tidyverse)
library(magrittr)
library(lme4)
library(FactoMineR)
library(factoextra)


#====Color Palettes=====
#palette for ComType c("#82DFB1", "#727CF1", "#A38EC9")
#palette for hilo c("#DF6957", "#3366CC")

#====Load Data==========
setwd("~/Documents/Projects/EcuadorDataProcessing/") #Set working directory
water_qual_perc <- read.csv("EnglishWorkingFile.csv", header = T, na.strings='999') #Read in English translated survey data

#====Create Quality Indicator======
water_qual_perc$TotalStreamQual <- (water_qual_perc$StreamColor + water_qual_perc$StreamOdor + water_qual_perc$StreamQuality)/3
hilototalstream <- vector()
for (i in 1:length(water_qual_perc$TotalStreamQual))
{
  if (water_qual_perc$TotalStreamQual[i] < 2.24) #2.24 is the mean of TotalStreamQual
  {
    hilototalstream[i] = 0
  }
  if (water_qual_perc$TotalStreamQual[i] >= 2.24)
  {
    hilototalstream[i] = 1
  }
}
water_qual_perc$hilo <- hilototalstream  #Assign this indicator to the dataframe
#====='Relevel' Usage Amount==========
water_qual_perc$UsageAmount <- water_qual_perc$DrinkCook + water_qual_perc$Irrigation + water_qual_perc$CattleTroughs + water_qual_perc$PersonalHygene + water_qual_perc$WashClothes + water_qual_perc$Recreation + water_qual_perc$RiverSports + water_qual_perc$PureAir + water_qual_perc$NaturalRelaxxation

#Due to uneven sample sizes, create three bins of usage:
for (i in 1:dim(water_qual_perc)[1]){
  if (water_qual_perc$UsageAmount[i] == 0){
    water_qual_perc$UsageBins[i] <- 0
  }
  if (water_qual_perc$UsageAmount[i] == 1 | water_qual_perc$UsageAmount[i] ==2 | water_qual_perc$UsageAmount[i] == 3){
    water_qual_perc$UsageBins[i] <- 1
  }
  if (water_qual_perc$UsageAmount[i] ==4 | water_qual_perc$UsageAmount[i] ==5 | water_qual_perc$UsageAmount[i] ==6 | water_qual_perc$UsageAmount[i] ==7 | water_qual_perc$UsageAmount[i] ==8| water_qual_perc$UsageAmount[i] ==9)
  { #Found a mistake
    water_qual_perc$UsageBins[i] <- 2
  }
}

#====Clean Data Columns==========
colnames(water_qual_perc)[1] <- 'ComCode' #Because the spanish accents in some areas give encoding errors

water_qual_perc[,c("ComCode", "Name", "Date", "ComName", "Other", "OtherProgram", "OtherSpec", "birdsamphibs", "X",
                   "X.1", "X.2", "Initiatives.of.Interest", "KindofChange", "Climate", "ChangeinClimate", "Quality.of.Water", 
                   "Profession", "Ninguno")] <- list(NULL)

unleveled_factors <- c("ComNumber", "Gender", "Residency", "PotableWater", "Piped.Water", "WaterQuebadra", "Sewer", 
                       "Telephone", "Internet", "PavedStreets", "GarbageCollection", "Internet", "Electricity", "DrinkCook",
                       "Irrigation", "CattleTroughs", "PersonalHygene", "WashClothes", "Recreation", "RiverSports", "PureAir",
                       "NaturalRelaxxation")
water_qual_perc %<>% mutate_at(unleveled_factors, factor)

factors_12345 <- c("Education", "WaterQuality")
water_qual_perc %<>% mutate_at(factors_12345, ~factor(., ordered = TRUE, levels = c(1,2,3,4,5)))

factors_0321 <- c("ProAgriChem", "ProTrash", "ProWaterContam", "ProAirContan", "Mosquitos", "Traffic", "Dust", "Noise", "Watewater",
                  "WaterServices", "AgriculturalUse", "LiveStock", "Debris", "Mining", "Deforestation", "ApplyLaws", "Fines", "Incentives",
                  "EducationPrograms")
water_qual_perc %<>% mutate_at(factors_0321, ~factor(., ordered =T, levels = c(0, 3, 2, 1)))

factors_123 <- c("ImpWaterQuality", "ImpElectricity", "ImpRoads", "ImpTrash", "ImpPublicSanitation", "ImpPublicSchool", "ImpAirQuality",
                 "ImpEnvironment", "ImpPublicTransport")
water_qual_perc %<>% mutate_at(factors_123, ~factor(., ordered = T, levels = c(1,2,3)))

#other special cases
water_qual_perc$ServiceLevel <- factor(water_qual_perc$ServiceLevel, ordered = TRUE, levels=c(1, 2))
water_qual_perc$Basic.Services <- factor(water_qual_perc$Basic.Services, ordered = TRUE, levels=c(3, 2, 1))
water_qual_perc$Age <- factor(water_qual_perc$Age, ordered = TRUE, levels = c(1, 2, 3, 4))
water_qual_perc$UsageAmount <- factor(water_qual_perc$UsageAmount, ordered = TRUE, levels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9))
water_qual_perc$UsageBins <- factor(water_qual_perc$UsageBins, ordered = T, levels=c(0, 1, 2))
water_qual_perc$hilo <- factor(water_qual_perc$hilo, levels = c(0, 1), labels=c("low", "high"))
water_qual_perc <- na.omit(water_qual_perc)

#=========Check for Collinearity ============
#Use Spearman's Rho to check for correlation

rho_mat <- matrix(NA, nrow = dim(water_qual_perc)[2], ncol = dim(water_qual_perc)[2]) #Create an empty nxn matrix
for (i in 1:length(water_qual_perc)){
  for (j in 1:length(water_qual_perc)){
    foo <- spearman.test(water_qual_perc[,i], water_qual_perc[,j]) #create a dummy variable
    if (foo$p.value <= 0.05) #If the correlation is significant...
    {
      rho_mat[i, j] <- foo$estimate #Assign it to the rho_matirx
    }
  }
}

rho_df <- as.data.frame(rho_mat) 
colnames(rho_df) <- colnames(water_qual_perc) #Assign names to the nxn dataframe
rownames(rho_df) <- colnames(water_qual_perc)


#====Determine Important Variables======
#First, create an MCA (multiple correspondence analysis)
mca_water <- water_qual_perc
#delete other columns that would affect the MCA
mca_water[,c("TotalStreamQual", "StreamColor", "StreamTurbidity", "StreamQuality", "StreamOdor", "PrimaryEnviroConcern",
                   "Concerns", "UsageAmount")] <- list(NULL)
res.mca <- MCA(mca_water, quanti.sup = c(3,4), quali.sup = c(58), graph = F)
eig.val <- get_eigenvalue(res.mca)
var <- get_mca_var(res.mca)

hilo_MCA <- fviz_mca_ind(res.mca, #explore the effects of different variables on the MCA
                         label = "none", # hide individual labels
                         habillage = c("hilo"), #insert variable of interest here
                         addEllipses = TRUE, ellipse.type = "confidence",
                         ggtheme = theme_minimal())
hilo_MCA + scale_color_manual(name="WQP", labels = c("High", "Low"),
                              values= c("#DF6957", "#3366CC")) #manually scale the colors 
#Get and sort contributions
important_vars <- res.mca$var$contrib
sum_impVars <- sort(rowSums(important_vars), decreasing = T)

quality_vars <- res.mca$var$cos2
sum_qualVars <- sort(rowSums(quality_vars), decreasing = T)

#Second, create MFA (multiple factor analysis)
mfa_water <- mca_water #make a new df, get rid of all aggregate counts
mfa_water[,c("Basic.Services", "ServiceLevel", "Piped.Water")] <- list(NULL)

#Groupings from survey
#ComNumber and ComType together? 2
#UTMX UTMY 2
#PotableWater, Piped.Water, WaterQuebadra, Sewer, 
#Telephone, Internet, PavedStreets, GarbageCollection, Electricity 9
#ImpWaterQuality, ImpElectricity, ImpRoads, ImpTrash, ImpPublicSanitation, ImpPublicTransport,
#ImpPublicSchool, ImpAirQuiality, ImpEnvironment 9
#ProAgriChem, ProTrash,ProWaterContam, ProAirContan, Mosquitos, Traffic, Dust, Noise, Watewater 9
#WaterQuality 1
#DrinkCook, Irrigation, CattleTroughs, Personal Hygene, WashClotehs, Recreation, RiverSports, PureAir, NaturalRelaxxation 9
#WaterServices, AgriculturalUse, LiveStock, Debris, Mining, Deforestation 6
#ApplyLaws, Fines, Incentives, EducationPrograms 4
#Gender, Age, Education, Residency 4
#UsageBins 1
#hilo 1
surveyQ <- c(2, 2, 8, 9, 9, 1, 9, 6, 4, 4, 1, 1)
group_names <- c("Community", "Geography", "Services", "Improvements", "Concerns", "DrinkingWQ", "Uses", "Contaminants", "Fixes", "Demography", "Usage Amount", "hilo")
group_types <- c("n", "s", "n", "n", "n", "n", "n", "n", "n", "n", "n", "n")

res.mfa <- MFA(mfa_water, 
               group = surveyQ, 
               type = group_types,
               name.group = group_names,
               num.group.sup = c(12),
               graph = FALSE)
fviz_screeplot(res.mfa)
fviz_mfa_var(res.mfa, "group")
fviz_cos2(res.mfa, choice = "quali.var", axes = 1)

fviz_mfa_ind(res.mfa, 
             habillage = "ComType", # color by groups 
             #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, ellipse.type = "confidence", 
             label ="none",
             repel = TRUE, # Avoid text overlapping
             axes = c(1,2)
) 

fviz_mfa_ind(res.mfa, partial = "all")
fviz_contrib(res.mfa, choice = "quali.var", axes = c(1,2), top = 30,
             palette = "jco")
fviz_mfa_quali_biplot(res.mfa, axes=c(1,2), select.ind = list(contrib = 5), 
                      select.var = list(contrib = 20),
                      repel = T,
                      ggtheme = theme_minimal())


#====Create Models======
mfa_model <- glmer(hilo~ComType+UsageBins+PersonalHygene+Irrigation+WashClothes+CattleTroughs+PureAir+PotableWater+Sewer+PavedStreets+
                      Debris+Deforestation+ImpRoads+Residency+
                      (1|ComNumber),
                    data=water_qual_perc, family = binomial(link="logit"),
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

mca_model <- glmer(hilo~UsageBins+PersonalHygene+DrinkCook+WashClothes+Sewer+
                     (1|ComNumber),
                   data=water_qual_perc, family = binomial(link="logit"),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))



