library(deSolve) 
library(ggplot2) 
library(tidyr)
library(ggsci)
library(wesanderson)
library(scales)
library(dplyr)
library(readr)
library(ggplot2)
library(mgcv)
library(deSolve)
library(tidyverse)
library(stats)
library(itsadug)
library(car)
library(broom)
library(viridis)
library(RColorBrewer)
library(pomp)
library(stringr)
library(glmmTMB)
library(cowplot)
library(bbmle)
library(ggeffects) # for model predictions


setwd("~/Desktop/Rscripts/Data")

TidalPhytoplankton <- read_csv("LivingResourcesReportedHUC8.csv",locale=locale(encoding="latin1"))
unique(TidalPhytoplankton$ScientificName)
#remove organisms named "NA"
TidalPhytoplankton <- dplyr::filter(TidalPhytoplankton,  !is.na(ScientificName))

TidalPhytoplanktonAbundances <- TidalPhytoplankton %>% group_by(Basin,ScientificName) %>% summarise(IndivPerBasin=sum(ReportingValue*15))

TidalRelComposition <- TidalPhytoplanktonAbundances %>%
  group_by(Basin) %>%
  mutate(TotalInBasin = sum(IndivPerBasin),
         RelAbundance = IndivPerBasin / TotalInBasin,
         RelAbunPercent = RelAbundance*100)

ggplot(TidalRelComposition, aes(x = Basin, y = RelAbundance, fill = ScientificName)) +
  geom_bar(stat = "identity") +
  ylab("Relative Abundance") +
  xlab("Basin") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(fill = "Taxon") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) 


threshold <- 0.05  # 5% threshold for rarity

TidalCompositionThreshold <- TidalRelComposition %>%
  mutate(ScientificName_Modified = ifelse(RelAbundance < threshold, "Other", ScientificName))

# Step 3: Plot with "Other" combined
ggplot(TidalCompositionThreshold, aes(x = Basin, y = RelAbundance, fill = ScientificName_Modified)) +
  geom_bar(stat = "identity") +
  ylab("Relative Abundance") +
  xlab("Basin") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


