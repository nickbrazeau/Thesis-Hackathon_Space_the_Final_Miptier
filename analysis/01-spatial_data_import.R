#.............
# Accessibility from Weiss
#.............
# see explanation in the DHS GC manual
# NOTE, this is from 2015
summary(dt$travel_times_2015)
hist(dt$travel_times_2015)
dt <- dt %>%
  dplyr::mutate(travel_times_2015_scale = scale(log(travel_times_2015 + tol), center = T, scale = T))

summary(dt$travel_times_2015_scale)
hist(dt$travel_times_2015_scale) # many, many 0s
hist(dt$travel_times_2015_scale) # standardization doesn't look as good, but should capture urban v. rural well
summary(dt$travel_times_2015_scale); sd(dt$travel_times_2015_scale)

