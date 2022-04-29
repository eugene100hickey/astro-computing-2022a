library(tidyverse)
library(glue)
library(RCurl)
library(gtools)
# library(caret)
# library(GGally)

# SQL that downloads some info on the chosen target from SDSS.
# ObjID from SDSS specifies the target
# set search parameters
N <- 10
sub_N <- 10
delta <- 0.3
# delta_gr_mag <- runif(n = 1, min = -2, max = +2)
# delta_ri_mag <- runif(n = 1, min = -2, max = +2)
bands_min <- 15
bands_max <- 20


get_spectrum <- function(object, wavelength_lower_limit = 5500, wavelength_upper_limit = 7000){
  plate <- object$plate
  mjd <- object$mjd
  fiber <- object$fiberid
  url_spect <- glue("http://dr12.sdss.org/csvSpectrum?plateid={plate}", 
                    "&mjd={mjd}&fiber={fiber}&reduction2d=v5_7_0")
  spectrum <- read_csv(file = url_spect)
  spectrum %>% 
    filter(between(Wavelength, wavelength_lower_limit, wavelength_upper_limit)) %>% 
    select(Wavelength, BestFit)
}

master_target_SqlQuery <- glue("SELECT top {N} p.ra, p.dec, ",
                               "p.u, p.g, p.r, p.i, p.z, p.objid, ", 
                               "s.specobjid, s.class, s.subclass, s.survey, ", 
                               "s.plate, s.mjd, s.fiberid ", 
                               "FROM photoObj AS p ", 
                               "JOIN SpecObj AS s ON s.bestobjid = p.objid ",
                               "WHERE p.g BETWEEN {bands_min} AND {bands_max} ",
                               "AND p.r BETWEEN {bands_min} AND {bands_max} ", 
                               "AND p.i BETWEEN {bands_min} AND {bands_max} ", 
                               "AND s.class = 'STAR' ",
                               "AND p.clean = 1 AND (p.calibStatus_r & 1) != 0",
                               "AND s.survey != 'eboss'",
                               "AND p.psfmagerr_u < 0.05",
                               "AND p.psfmagerr_g < 0.05",
                               "AND p.psfmagerr_r < 0.05",
                               "AND p.psfmagerr_i < 0.05",
                               "AND p.psfmagerr_z < 0.05" )

# downloads target data
# dataframe target has necessary info
master_target_SqlQuery <- str_squish(master_target_SqlQuery)
urlBase <- "http://skyserver.sdss.org/dr15/SkyserverWS/SearchTools/SqlSearch?"
X <- getForm(urlBase, cmd = master_target_SqlQuery, format = "csv")
master_targets <- read.table(text = X, header = TRUE, sep = ",", dec = ".", comment.char = "#") %>% 
  mutate(objid = as.character(objid),
         specobjid = as.character(specobjid))


# Fri Apr 24 18:55:25 2020 ------------------------------
# index <- createDataPartition(master_targets$survey, 
#                              p = sub_N/N, list = F)
# master_targets <- master_targets[index,]
# Fri Apr 24 18:55:36 2020 ------------------------------

get_correlation <- function(index = 8){
  progress_bar$tick()$print()
  delta_gr_mag <- master_targets$g[index] - master_targets$r[index]
  delta_ri_mag <- master_targets$r[index] - master_targets$i[index]
  plate <- master_targets$plate[index]
  
  match_SqlQuery <- glue("SELECT top 100 p.ra, p.dec, ",
                         "p.u, p.g, p.r, p.i, p.z, p.objid, ", 
                         "s.specobjid, s.class, s.subclass, s.survey, ", 
                         "s.plate, s.mjd, s.fiberid ", 
                         "FROM photoObj AS p ", 
                         "JOIN SpecObj AS s ON s.bestobjid = p.objid ",
                         "WHERE p.g BETWEEN {bands_min} AND {bands_max} ",
                         "AND p.r BETWEEN {bands_min} AND {bands_max} ", 
                         "AND p.i BETWEEN {bands_min} AND {bands_max} ", 
                         "AND s.class = 'STAR' ",
                         "AND s.scienceprimary = 1",
                         "AND s.zWarning = 0",
                         "AND p.clean = 1 AND (p.calibStatus_r & 1) != 0",
                         "AND p.psfmagerr_u < 0.05",
                         "AND p.psfmagerr_g < 0.05",
                         "AND p.psfmagerr_r < 0.05",
                         "AND p.psfmagerr_i < 0.05",
                         "AND p.psfmagerr_z < 0.05",
                         "AND (p.g-p.r) BETWEEN {delta_gr_mag} AND {delta_gr_mag+delta} ",
                         "AND (p.r-p.i) BETWEEN {delta_ri_mag} AND {delta_ri_mag+delta} ",
                         "AND s.survey != 'eboss' ", 
                         "AND s.plate != {plate}" )
  match_SqlQuery <- str_squish(match_SqlQuery)
  X <- getForm(urlBase, cmd = match_SqlQuery, format = "csv")
  match <- read.table(text = X, header = TRUE, sep = ",", dec = ".", comment.char = "#")
  match <- match[sample(1:100, size = 1),]
  spect1 <- get_spectrum(master_targets[index,])
  spect2 <- get_spectrum(match)
  cor_value <- cor(spect1$BestFit, spect2$BestFit)
  bind_cols(star_1=match, star_2=master_targets[index,], cor = cor_value)
}

# Wed Apr 22 18:41:31 2020 ------------------------------
progress_bar <- progress_estimated(nrow(master_targets))

get_correlation_safely <- safely(get_correlation)

z <- map(1:nrow(master_targets), ~get_correlation_safely(.x)) %>% 
  map_df("result") %>% 
  compact()

z <- z %>% 
  mutate(cor_logit = logit(cor, min = -1, max = 1) %>% 
           scale()) %>% 
  janitor::clean_names() %>% 
  mutate(specobjid_9 = as.character(specobjid_9),
         specobjid_24 = as.character(specobjid_24))