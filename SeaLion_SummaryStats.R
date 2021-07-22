# Tag duration values from Progress Report: 
# https://drive.google.com/drive/u/0/folders/1O4Pmj-SGIdaV3tSbhV7BIahFfnsSEWLi

tags <- c(235, 161, 209, 248, 158,184, 173, 164, 205, 202, 262)  
mean(tags)
sd(tags)
sd(tags)/sqrt(length(tags)) # Standard error

# I'm reporting standard deviation because I'm more interested in the variability around the 
# sample mean than I am in estimating how accurately the sample mean reflects the true population mean. 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1255808/


