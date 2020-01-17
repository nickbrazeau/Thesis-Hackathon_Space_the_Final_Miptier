distdatgood <- readxl::read_excel("~/Documents/GitHub/Space_the_Final_Miptier/ignore/distdatgood.xlsx")
nick <- get_distance_geno_likelihood(name = "nick",
                             distdata = distdatgood,
                             scalar = 1,
                             imputclst = 99999,
                             global_fii = 99999)


distdatbad <- readxl::read_excel("~/Documents/GitHub/Space_the_Final_Miptier/ignore/distdatbad.xlsx")
nick2 <- get_distance_geno_likelihood(name = "nick",
                             distdata = distdatbad,
                             scalar = 1,
                             imputclst = 99999,
                             global_fii = 99999)

nick <- unlist(nick$pij)
nick <- nick[!is.na(nick)]


nick2 <- unlist(nick2$pij)
nick2 <- nick2[!is.na(nick2)]


