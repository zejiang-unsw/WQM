## code to prepare `DATASET` dataset goes here
## ACCESS NWP model----
library(dplyr)
library(data.table)
load(paste0("MWR2022/data-raw/dat.tmp.Rdata"))

#summary(dat.tmp)
station.no <- colnames(dat.tmp)[-c(1:3)]
## sampled station
Ind.samp <- c(34,75)
station.no[Ind.samp]

#------------------------------------------------------------------------------#
###sample station
dat.tmp.n <- dat.tmp %>% subset(Lead==1)
NWP.rain.samp <- dat.tmp.n[, c(1:3, Ind.samp+3)] %>% data.table() %>% melt(id=1:3, variable="Station", value="Rain") %>%
  dcast(Lead+Date+Station~Model, value.var = 'Rain') %>% mutate_if(is.factor,as.character) #%>%

data.h.df <- NWP.rain.samp %>% na.omit() %>% arrange(Date) %>%
  group_by(Station) %>%  group_split() %>% lapply(FUN=as.data.frame)

usethis::use_data(data.h.df, overwrite = TRUE)
###all stations
NWP.rain.h <- dat.tmp %>% subset(Lead==1) %>% data.table() %>% melt(id=1:3, variable = "Station", value = "Rain") %>%
  dcast(Lead+Date+Station~Model, value.var = 'Rain') %>% mutate_if(is.factor,as.character)

usethis::use_data(NWP.rain.h, overwrite = TRUE)
