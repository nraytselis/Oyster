library(deSolve) 
library(ggplot2) 
library(tidyr)
library(ggsci)
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
library(sf) #gis
library(ggspatial) #gis
library(vegan) #beta diversity

setwd("~/Desktop/Rscripts/Data")

TidalPhytoplankton <- read_csv("LivingResourcesReportedHUC8.csv",locale=locale(encoding="latin1"))
unique(TidalPhytoplankton$ScientificName)

#remove organisms named "NA"
TidalPhytoplankton <- dplyr::filter(TidalPhytoplankton,  !is.na(ScientificName))

#extract year
TidalPhytoplankton$SampleDate <- as.Date(TidalPhytoplankton$SampleDate , 
                                                   format = "%m/%d/%y")
TidalPhytoplankton <- separate(TidalPhytoplankton, SampleDate, c('year', 'month', 'day'), sep = "-",remove = FALSE)
TidalPhytoplankton[order(TidalPhytoplankton$SampleDate ),]

#calculate # of individuals per basin per year 
TidalPhytoplanktonAbundances <- TidalPhytoplankton %>% group_by(Basin,ScientificName,year) %>% 
  summarise(IndivPerBasin=sum(ReportingValue*15))


#summary table with relative abundances 
TidalRelComposition <- TidalPhytoplanktonAbundances %>%
  group_by(Basin,year) %>%
  mutate(TotalInBasin = sum(IndivPerBasin),
         RelAbundance = IndivPerBasin / TotalInBasin,
         RelAbunPercent = RelAbundance*100)

threshold <- 0.025  # 5% threshold for rarity

TidalCompositionThreshold <- TidalRelComposition %>%
  mutate(ScientificName_Modified = ifelse(RelAbundance < threshold, "Other", ScientificName))

species <- unique(TidalCompositionThreshold$ScientificName_Modified)
#species
#ggplot color palette
base_colors <- c(brewer.pal(12, "Set3"),
                 brewer.pal(8, "Dark2"),
                 brewer.pal(8, "Accent"),
                 brewer.pal(9, "Paired"))

nb.cols <- length(unique(TidalCompositionThreshold$ScientificName_Modified))

mycolors <- colorRampPalette(base_colors)(nb.cols)

# Plot with "Other" combined
ggplot(TidalCompositionThreshold, aes(x = Basin, y = RelAbundance, fill = ScientificName_Modified)) +
  geom_bar(stat = "identity") +
  ylab("Relative Abundance") +
  xlab("Basin") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) + scale_fill_manual(values = setNames(mycolors[1:length(species)], species)) + facet_wrap(~ year) +
  ggtitle("Relative abundances of phytoplankton across basins of the Chesapeake Bay (> 2.5% preveleance)") 


#genus level differences
TidalRelComposition$genus = str_extract(TidalRelComposition$ScientificName, "^[^ ]+")
unique(TidalRelComposition$genus)

#calculate the individuals per genus per basin per year
TidalRelCompositionGenus = TidalRelComposition %>% group_by(Basin,year,genus) %>% 
  summarise(IndivPerGenus = sum(IndivPerBasin)) 

#calculate the relative abundances of each genus
TidalRelCompositionGenus = TidalRelCompositionGenus %>%
  group_by(Basin, year) %>%  
  mutate(RelAbunGenus = IndivPerGenus / sum(IndivPerGenus)) 

#impose 2.5% threshold
threshold <- 0.025  # 5% threshold for rarity

TidalRelCompositionGenus <- TidalRelCompositionGenus %>%
  mutate(genus_updated = ifelse(RelAbunGenus < threshold, "Other", genus))

#update color palette 
genus <- unique(TidalRelCompositionGenus$genus_updated)

nb.cols.genus <- length(unique(TidalRelCompositionGenus$genus_updated))

mycolorsgenus <- colorRampPalette(base_colors)(nb.cols.genus)

# Plot with "Other" combined at genus level
ggplot(TidalRelCompositionGenus, aes(x = Basin, y = RelAbunGenus, fill = genus_updated)) +
  geom_bar(stat = "identity") +
  ylab("Relative Abundance Genus") +
  xlab("Basin") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right")+ facet_wrap(~ year) + scale_fill_manual(values = setNames(mycolorsgenus[1:length(genus)], genus)) +
  ggtitle("Relative abundances of phytoplankton genus across basins of the Chesapeake Bay (> 2.5% preveleance)") 






#diversity metrics
# Group by year, then pivot
TidalRelComposition_by_year <- TidalRelComposition %>%
  group_by(year) %>%
  group_split() %>% #splits the dataframe into separate data frames for each year in a list
  purrr::map(~ .x %>%
               group_by(ScientificName, Basin) %>%
               summarise(TotalIndividuals = sum(IndivPerBasin), .groups = "drop") %>%
               pivot_wider(names_from = Basin,
                           values_from = TotalIndividuals,
                           values_fill = 0) %>%
               column_to_rownames("ScientificName") %>%
               as.matrix()
  )


TidalRelComposition_by_year <- TidalPhytoplankton %>%
  mutate(IndivsPerSample = ReportingValue * 15) %>%  # Adjust for volume: individuals per sample
  group_by(year) %>%
  group_split() %>% #splits the dataframe into separate data frames for each year in a list
  purrr::map(~ .x %>% #applies the following to each year (element in list)
               group_by(ScientificName, Basin) %>%
               summarise(TotalIndivs = sum(IndivsPerSample), .groups = "drop") %>%  # Sum across all samples in each basin
               pivot_wider(
                 names_from = Basin,
                 values_from = TotalIndivs,  # ✅ Correct column name
                 values_fill = list(TotalIndivs = 0)
               ) %>%
               column_to_rownames("ScientificName") %>%
               as.matrix()
  )

#add years to the data
names(TidalRelComposition_by_year) <- TidalPhytoplankton %>%
  pull(year) %>%
  unique() %>%
  sort()

#calculate bray curtis dissimilatiry matrix 
bray_by_year <- purrr::map(TidalRelComposition_by_year, ~ vegdist(.x, method = "bray"))

bray_summary <- map_df(bray_by_year, ~ {
  data.frame(
    mean_dissimilarity = mean(.x),
    median_dissimilarity = median(.x),
    max_dissimilarity = max(.x),
    min_dissimilarity = min(.x)
  )
}, .id = "year")

print(bray_summary) #average dissimilarity between all basin-pairs within each year 



# Now, TidalRelComposition_by_year is a list of matrices, one per year.


#TidalRelComposition_by_year[[1]] 2020 dataframe

vegdist(TidalRelComposition_by_year,method="bray")






#import in shape file
#https://mdl.library.utoronto.ca/technology/tutorials/introduction-gis-using-r 
BayShapeFile = read_sf("Chesapeake_Bay_Watershed-shp") #https://data.virginia.gov/gl/dataset/chesapeake-bay-watershed/resource/3194945f-a972-4f1f-ba86-eb45a290ec7d
ggplot() + geom_sf(data=BayShapeFile)





