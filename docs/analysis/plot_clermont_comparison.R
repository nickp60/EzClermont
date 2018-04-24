# set working directory to clermontpcr
setwd("~/GitHub/clermontpcr/docs/analysis/")

setwd("~/GitHub/Nolan_et_al_2018/")

metadata <- read.csv("data/metadata_clermont13.csv", stringsAsFactors = F)

expdata <- read.csv("2018-04-20-combined-redo.txt", sep="\t", col.names = c("long_acc", "exp_group"), stringsAsFactors = F, header=F)
expdata$acc <- gsub("(GCA_.*?\\.\\d)_(.*)","\\1", expdata$long_acc)

keyfile <- read.csv("data/accession_name.txt", sep="\t", col.names = c("acc", "ecor"), stringsAsFactors = F, header=F)
keyfile$acc[!keyfile$acc %in% expdata$acc]
dat <- merge(keyfile, expdata, by="acc")
dat <- merge(dat, metadata, by.x="ecor", by.y="strain")
dat$match <- dat$exp_group == dat$quardiplex_pcr_base
table(dat$exp_group == dat$quardiplex_pcr_base)
View(dat[,c("acc", "ecor", "quardiplex_pcr_base", "exp_group", "match")])
