library(dplyr)

#load data
df1 <-read.csv("C:/Users/Sandra/OneDrive/Documentos/GitHub/MSc.-Communal-roosting-behav/Chapters/Brain size data and sources/AllDataSourcesCombined.csv")
df2 <- Bird_data <- read.csv("C:/Users/Sandra/OneDrive - UBC/PhD proposal/Chapter 2/Updated database/Chapter_2_PhD_NEWdata_updated.csv")

#Delete duplicates and merge data

#Delete 1
df1_clean <- df1 %>%
  group_by(Species) %>%
  # drop rows that are in a duplicated‐Species group AND have no Mass
  filter(!(n() > 1 & is.na(body.mass))) %>%
  ungroup()

#Delete 2
df2_clean <- df2 %>%
  group_by(Species) %>%
  # drop rows that are in a duplicated‐Species group AND have no Mass
  filter(!(n() > 1 & is.na(Raw_CRB))) %>%
  ungroup()

#Merge
df2_updated <- df2_clean %>%
  left_join(
    df1_clean %>% select(Species, brain.mass.g),
    by = "Species"
  )


#Save the data for brain data modeling
write.csv(df2_updated, file = "C:/Users/Sandra/OneDrive - UBC/PhD proposal/Chapter 2/Brain size data and sources/Bird_data.csv") #save model object for future use

