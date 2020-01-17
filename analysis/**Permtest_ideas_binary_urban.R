
#..............................................................
# I have 43 btwn cluster pairs
# --> potentially 86 draws from urban rural from all clst
# and then i have my real porportion
#..............................................................
# 221 rural
# 130 urban clsts
set.seed(48)
# 0 is rural, 1 is urban
reps <- 1e4
ret <- sapply(1:reps, function(x){
  return(mean(sample(x = c(0, 1), size = 86, prob = c(221,130), replace = T)))
  })
summary(ret)
hist(ret)
obs <- 14/(14+34)
obs

quantile(ret, 0.025)

#..............................................................
# contigency tables
# but not these are dependent on cells within
#..............................................................
reps <- 1000
stat <- lapply(1:reps, function(x) return(list()))
for (i in 1:reps) {
  c1 <- sample(x = c("R", "U"), size = 43, prob = c(221,130), replace = T)
  c2 <- sample(x = c("R", "U"), size = 43, prob = c(221,130), replace = T)

  ret <- data.frame(c1 = c1, c2 = c2)
  ret <- ifelse(ret$c1 == "U" & ret$c2 == "U", "UU",
                ifelse(ret$c1 == "R" & ret$c2 == "R", "RR",
                       ifelse(ret$c1 == "R" & ret$c2 == "U" |
                                ret$c2 == "R" & ret$c1 == "U", "RU", NA)))
  ret <- factor(ret, levels = c("RR", "RU", "UU"))
  stat[[i]] <- table(ret)

}



chisq.test(table(ibD.meiotic.urb$urbnclss))
