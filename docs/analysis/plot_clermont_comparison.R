# set working directory to clermontpcr
setwd("~/clermontpcr/docs/analysis/")

metadata <- read.csv("data/metadata_clermont13.csv", stringsAsFactors = F)

expdata <- read.csv("2018-04-20-ecor_results.txt", sep="\t", col.names = c("long_acc", "exp_group"), stringsAsFactors = F)
expdata$acc <- gsub( "GCF", "GCA", gsub("(.*\\.1)_(.*)","\\1", expdata$long_acc))

keyfile <- read.csv("data/accession_name.txt", sep="\t", col.names = c("acc", "ecor"), stringsAsFactors = F)

dat <- merge(keyfile, expdata, by="acc",)
dat <- merge(dat, metadata, by.x="ecor", by.y="strain")

table(dat$exp_group == dat$quardiplex_pcr_base)
