# set working directory to clermontpcr
setwd("~/GitHub/clermontpcr/docs/analysis/")

#setwd("~/GitHub/Nolan_et_al_2018/")

metadata <- read.csv("validate/metadata_clermont15.csv", stringsAsFactors = F)

expdata <- read.csv("validate/2018-04-30.txt", sep="\t", col.names = c("long_acc", "exp_group"), stringsAsFactors = F, header=F)
expdata$acc <- gsub("(GCA_.*?\\.\\d)_(.*)","\\1", expdata$long_acc)

keyfile <- read.csv("validate/accession_name.txt", sep="\t", col.names = c("acc", "ecor"), stringsAsFactors = F, header=F)
keyfile$acc[!keyfile$acc %in% expdata$acc]
dat <- merge(keyfile, expdata, by="acc")
dat <- merge(dat, metadata, by.x="ecor", by.y="strain")


dat$raw_match <- dat$exp_group == dat$quardiplex_pcr_base

# fix mismatch of particular clade type
dat$match <- ifelse(grepl("Clade", dat$quardiplex_pcr_base) & 
                          grepl("cryptic", dat$exp_group) & !dat$raw_match,
                        T,
                        dat$raw_match)



table(dat$raw_match)
table(dat$match)

#View(dat[,c("acc", "ecor", "quardiplex_pcr_base", "exp_group", "match")])

write.table(file = "~/Desktop/matches.txt", dat[,c("acc", "ecor", "quardiplex_pcr_base", "exp_group", "match")])

## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
## biocLite("BiocUpgrade") ## you may need this
#biocLite("ggtree")
library(ggtree)



tree <- read.newick("./validate/kSNP_output_k19/tree_AlleleCounts.parsimony.tre")

dat$label <- substr(dat$acc, 1, 13)

labdf <- merge(data.frame(label=tree$tip.label), dat, by="label")
#labdf$quardiplex_pcr_base_base <- ifelse(labdf$quardiplex_pcr_basemerge(data.frame(label=tree$tip.label), dat, by="label")

# add labels for where ezclermont went wrong
labdf$label2 <- ifelse(labdf$match,
                       labdf$ecor,
                       paste0("(", labdf$exp_group, ") ", labdf$ecor))
# combine all the clades into one
labdf$quardiplex_pcr_base_tidy <- ifelse(grepl("Clade", labdf$quardiplex_pcr_base),
                       "Cryptic", labdf$quardiplex_pcr_base)
labdf$matchsize <- ifelse(labdf$match,  1, 5)


colors = c(A="grey", B1="darkgreen", B2="purple",   C="orange", "Clade I"="yellow", "Cryptic"="yellow", "Clade II"="yellow", "Clade IV"="yellow", "Clade V"="yellow", D="red", E="green", F="steelblue")

(tp<- ggtree(tree, layout="circular", branch.length = 'none') %<+% labdf + 
    geom_tiplab2(aes(label=label2), size = 3, offset=1) + #, color=tmatch)) +
    # geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
    geom_hilight(node=138, fill=colors['F'], alpha=.6) + # F
    geom_hilight(node=145, fill=colors['D'], alpha=.6) + # D
    geom_hilight(node=150, fill=colors['Clade II'], alpha=.6) + # clades
    geom_hilight(node=109, fill=colors['E'], alpha=.6) + # E
    geom_hilight(node=154, fill=colors['B1'], alpha=.6) + # B1
    geom_hilight(node=170, fill=colors['C'], alpha=.6) + # C
    geom_hilight(node=116, fill=colors['B2'], alpha=.6) + # B2
    geom_tippoint(aes(size=matchsize , shape=match, color=quardiplex_pcr_base_tidy), alpha=0.9) + 
    scale_color_manual(values = colors) + 
    scale_size_area(guide=F)  +
    scale_shape_discrete(guide=F) +
    labs(title="", color="Phylogroup" ) +
    # 
    # labs(title="Parsimony cladogram of strains from Clermont, et al 2015", color="Phylogroup", subtitle="Tree generated with kSNP3 (k=19). Enlarged circular tips show where EzClermont \ndiffered from reported phylogroup (EzClermont type show in brackets)." ) +
    theme(legend.position="right") 
  
)
png(units = "in", filename = "cladogram.png", width=9, height = 9,res = 300)
print(tp)
dev.off()
pdf(file = "cladogram.pdf", width=9, height = 9)
print(tp)
dev.off()
