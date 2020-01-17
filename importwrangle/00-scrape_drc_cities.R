#..............................................................
# Scrape All DRC Cities from Wikipedia
#..............................................................
library(rvest)
library(tidyverse)


#..............................................................
# scrape
#..............................................................
url <- "https://en.wikipedia.org/wiki/List_of_cities_in_the_Democratic_Republic_of_the_Congo"
cities <- url %>%
  read_html() %>%
  html_nodes(xpath='//*[@id="mw-content-text"]/div/table[2]') %>%
  html_table()
cities <- cities[[1]]

#..............................................................
# Tidy up
#..............................................................
# coords
cities <- cities %>%
  magrittr::set_colnames(c("city", "population", "coords", "province", "formname")) %>%
  dplyr::mutate(
    latnum = stringr::str_split_fixed(coords, "/", n = 3)[,3],
    latnum = stringr::str_split_fixed(latnum, ";", n = 2)[,1],

    longnum = stringr::str_split_fixed(coords, "/", n = 3)[,3],
    longnum = stringr::str_split_fixed(longnum, ";", n = 2)[,2],
    longnum = gsub(" ", "", longnum),
    longnum = stringr::str_extract(longnum, "^-?[0-9]+\\d*(\\.\\d+)?"),

    longnum = as.numeric(longnum),
    latnum = as.numeric(latnum),

    population = gsub(",", "", population),
    population = as.numeric(population)
    ) %>%
  dplyr::select(c("city", "population", "province", "longnum", "latnum"))



# drop any with missing data
cities <- na.omit(cities)

#..............................................................
# save out
#..............................................................
readr::write_csv(cities,
                 path = "~/Documents/GitHub/Space_the_Final_Miptier/data/map_bases/DRC_city_coordinates.csv")





