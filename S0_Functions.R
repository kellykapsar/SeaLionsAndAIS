




# A function equivalent to += in Python 
`%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))


# This function takes a date and converts into a season and year
# Season is calculted based on the meteorological definition (e.g., equinox/solstice)
# Winter dates in December are attributed to the following year
# (i.e., December 21-31, 2019 is classified as part of Winter 2020)
getSeasonYear <- function(input.date){
  numeric.date <- 100*month(input.date)+day(input.date)
  yrs <- year(input.date)
  yearsnew <- unlist(lapply(1:length(yrs), function(x) ifelse(numeric.date[x] > 1220, yrs[x] %+=% 1, yrs[x])))
  ## input Seasons upper limits in the form MMDD in the "break =" option:
  cuts <- base::cut(numeric.date, breaks = c(0,319,0620,0921,1220,1231)) 
  # rename the resulting groups (could've been done within cut(...levels=) if "Winter" wasn't double
  levels(cuts) <- c("Winter","Spring","Summer","Fall","Winter")
  seasonyear <- paste0(cuts, "-", yearsnew)
  return(seasonyear)
}