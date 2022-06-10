# Check to see if sf package is installed and, if not, install it 
list.of.packages <- c("sf", "ggplot2","googledrive")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load packages 
library(sf)
library(ggplot2)
library(googledrive)

# See above link for info on how to fill out

drive_auth()

drive_download("https://drive.google.com/file/d/1LmNGwAdt1y_Pwb49hd6k99CDvVBe5250/view?usp=sharing", type = "csv")
test <- read.csv("test.csv")
