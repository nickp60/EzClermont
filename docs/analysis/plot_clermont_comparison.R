# set working directory to clermontpcr
#  setwd("~/GitHub/ezclermont//docs/analysis/")


metadata <- read.csv("validate/validation_metadata.csv", stringsAsFactors = F)
metadata$simple <-  gsub("(.*?) (.*)","\\1", metadata$reported_phylogroup)


expdata <- read.csv("validate/2020-04-06-ezclermont.txt", 
                    sep="\t", col.names = c("long_acc", "ezClermont_phylogroup"), 
                    stringsAsFactors = F, header=F)
expdata$accession <- gsub("(GCA_.*?\\.\\d)_(.*)","\\1", expdata$long_acc)
together <- merge(metadata, expdata, by="accession", all.y = T)


together$raw_match <- together$simple == together$ezClermont_phylogroup

# fix mismatch of particular clade type
together$match <- ifelse(grepl("Clade", together$reported_phylogroup) & 
                          grepl("cryptic", together$ezClermont_phylogroup) & !together$raw_match,
                        T,
                        together$raw_match)



table(together$raw_match)
table(together$match)

#View(dat[,c("acc", "ecor", "quardiplex_pcr_base", "exp_group", "match")])

write.table(file = file.path(".", "matched_results.txt"), together[,c("accession", "strain", "reported_phylogroup", "ezClermont_phylogroup", "match")])

## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
## biocLite("BiocUpgrade") ## you may need this
#biocLite("ggtree")
library(ggplot2)
library(ggtree)



tree <- read.tree("./validate/kSNP_output_k19/tree_AlleleCounts.parsimony.tre")

together$label <- substr(together$acc, 1, 13)

labdf <- merge(data.frame(label=tree$tip.label), together, by="label")
#labdf$quardiplex_pcr_base_base <- ifelse(labdf$quardiplex_pcr_basemerge(data.frame(label=tree$tip.label), dat, by="label")

# add labels for where ezclermont went wrong
labdf$label2 <- ifelse(labdf$match,
                       labdf$strain,
                       paste0("(", labdf$ezClermont_phylogroup, ") ", labdf$strain))
# combine all the clades into one
labdf$quardiplex_pcr_base_tidy <- ifelse(grepl("Clade", labdf$simple),
                       "Cryptic", labdf$simple)
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
    guides(colour = guide_legend(override.aes = list(size=3))) +
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
