# set working directory to clermontpcr
#  setwd("~/GitHub/ezclermont//docs/analysis/")
library(tidyverse)
library(patchwork)

metadata <- read.csv("validate/validation_metadata.csv", stringsAsFactors = F)
metadata$simple <-  gsub("(.*?) (.*)","\\1", metadata$reported_phylogroup)
metadata <- metadata %>% select(-dataset, -reported_phylogroup) %>% filter(accession != "" )

expdata <- read.csv("validate/2020-04-09-ezclermont.txt", 
                    sep="\t", col.names = c("long_acc", "ezClermont_phylogroup"), 
                    stringsAsFactors = F, header=F)
expdata$accession <- gsub(".scaffold.fasta", "", 
                          gsub("(GCA_.*?\\.\\d)_(.*)","\\1", expdata$long_acc))
together <- left_join(metadata, expdata %>% select(-long_acc), by="accession")


together$raw_match <- together$simple == together$ezClermont_phylogroup

# fix mismatch of particular clade type
together$match <- ifelse(grepl("Clade", together$simple) & 
                          grepl("cryptic", together$ezClermont_phylogroup) & !together$raw_match,
                        T,
                        together$raw_match)

table(together$raw_match)
table(together$match)
together$raw_match <- NULL

expdata_ct <- read.csv("validate/2020-04-09-CT-results.txt", 
                    sep="\t", col.names = c("accession", "alleles", "profile", "al", "ct_phylotype", "mashpath"), 
                    stringsAsFactors = F, header=F) %>%
  select(accession, ct_phylotype) %>%
  transform(accession= gsub(".scaffold.fasta", "", gsub("(.*?\\..)_(.*)", "\\1", accession)))
expdata_ct$ct_phylotype <- gsub("(.*?) (.*)","\\1", expdata_ct$ct_phylotype)


together <- left_join(together, expdata_ct, by="accession")

write.table(file = file.path(".", "matched_results.txt"), together[,c("accession", "strain", "reported_phylogroup", "ezClermont_phylogroup", "match")])

## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
## biocLite("BiocUpgrade") ## you may need this
#biocLite("ggtree")
library(ggplot2)
library(ggtree)
library(ape)



#tree2 <- read.tree("./validate/kSNP_output_k19/tree_AlleleCounts.parsimony.tre")
tree <- read.tree("./validate/alignment/parsnp.clean.phy_phyml_tree.txt")

# we had to triim the ma,e becuse of the Phylip format
tree$tip.label <- ifelse(
  !startsWith(tree$tip.label, "0"),
  paste0("ESC_", tree$tip.label),
  paste0("GCA_", tree$tip.label))
together$label <- substr(together$acc, 1, 13)

labdf <- left_join(data.frame(label=tree$tip.label, stringsAsFactors = FALSE), together, by="label")
#labdf$quardiplex_pcr_base_base <- ifelse(labdf$quardiplex_pcr_basemerge(data.frame(label=tree$tip.label), dat, by="label")

# add labels for where ezclermont went wrong
labdf$label2 <- ifelse(labdf$match,
                       labdf$strain,
                       paste0("(", labdf$ezClermont_phylogroup, ") ", labdf$strain))
# combine all the clades into one
labdf$quardiplex_pcr_base_tidy <- ifelse(grepl("lade", labdf$simple),
                       "Cryptic", labdf$simple)
labdf$ct_phylotype <- ifelse(grepl("lade", labdf$ct_phylotype),
                                         "Cryptic", ifelse(grepl("Unknown", labdf$ct_phylotype), 
                                                           "U", labdf$ct_phylotype))
labdf$ezClermont_phylogroup <- ifelse(grepl("cryptic", labdf$ezClermont_phylogroup), "Cryptic", labdf$ezClermont_phylogroup)
labdf$matchsize <- ifelse(labdf$match,  1, 5)



true_phy <- read.table(stringsAsFactors = F, header=T,sep=",", text="
strain,Phylogeny,Note
APEC01,A,
ECOR04,A,
ECOR12,A,
ECOR40,F,
ECOR44,D,ArpA1_r mutationn position 17
ECOR46,F/G,contaminated assembly
ECOR51,B2,
ECOR58,B1,
ECOR70,C,
ECOR71,C,
ECOR43,E,
ECOR23,B2,
ECOR07,B1,
ECOR49,D,
ECOR72,B1,contaminaed assembly
SMS-3-5,F,reported as phylogroup F in Vangchhia et al 2016
")

results_table <- labdf %>%
  select(-label, -label2, -match, -simple, -matchsize) %>%
    rename(
      "Strain"=strain,
      "Accession"=accession,
      "Reported"=quardiplex_pcr_base_tidy,
      "ClermonTyping"=ct_phylotype,
      "EzClermont"=ezClermont_phylogroup) %>%
    select(Accession,Strain, Reported, ClermonTyping, EzClermont) 

dput(
  left_join(results_table, 
            true_phy, by=c("Strain"="strain")) %>% 
       filter (Reported != EzClermont | 
                 Reported != ClermonTyping | 
                 EzClermont != ClermonTyping) %>%
  select(Strain, Accession, Reported, Phylogeny, ClermonTyping, EzClermont, Note) 
)


colors = c(A="grey", B1="darkgreen", B2="purple",   
           C="orange", "Clade I"="yellow", "Cryptic"="yellow",
           "Clade II"="yellow", "Clade IV"="yellow", 
           "Clade V"="yellow", D="red", E="green", 
           F="steelblue", G="black", U="pink")

# 
# for ( i in 126:200){
#   readline(prompt=paste("Press [enter] to show reroot at", i))
#   print(ggtree(root.phylo(tree, node = i)) + geom_tiplab())
# }
# ggtree(tree, layout="circular", branch.length = 'none' ) %<+% labdf + 
#   geom_tiplab2(aes(label=label2, color=quardiplex_pcr_base_tidy), size = 3, offset=5) + #, color=tmatch)) +
#   scale_color_manual(values = colors) + 
#   geom_text2(aes(subset=!isTip, label=node))
tp <- gheatmap( ggtree(ape::root.phylo(tree, node = 142), layout="circular", branch.length = 'none' ) %<+% labdf + 
    #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
    geom_hilight(node=145, fill=colors['A'], alpha=.6) + # A
    geom_hilight(node=179, fill=colors['B1'], alpha=.6) + # B1
    geom_hilight(node=218, fill=colors['B2'], alpha=.6) + # B2
    geom_hilight(node=199, fill=colors['C'], alpha=.6) + # C
    geom_hilight(node=212, fill=colors['D'], alpha=.6) + # D
    #geom_hilight(node=186, fill=colors['Clade II'], alpha=.6) + # clades are made the root
    geom_hilight(node=202, fill=colors['E'], alpha=.6) + # E
    geom_hilight(node=135, fill=colors['F'], alpha=.6) + # F
    geom_hilight(node=132, fill=colors['G'], alpha=.6) + # G 
    # geom_tippoint(aes(x=x, size=matchsize , shape=match, color=quardiplex_pcr_base_tidy), alpha=1) + 
    # geom_tippoint(aes(x=x+3, size=matchsize , shape=match, color=ct_phylotype), alpha=0.7) + 
    # geom_tippoint(aes(x=x+6, size=matchsize , shape=match,  color=ezClermont_phylogroup), alpha=0.4) + 
    geom_tiplab2(aes(label=label2), size = 3, offset=5) + #, color=tmatch)) +
    scale_color_manual(values = colors) + 
    scale_size_area(guide=F)  +
    scale_shape_discrete(guide=F) +
    labs(title="", color="Phylogroup" ) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="right")     ,
    
    data= results_table %>% 
      transform(Accession=substring(first = 1, last = 13, text = Accession)) %>% 
      select(-Strain) %>%
      column_to_rownames("Accession") ,
    color="black", 
    colnames_position="top", 
    offset=.004,
    colnames_angle=90, colnames_offset_y = 1, 
    width=.2,
    hjust=.5, font.size=2
) +  scale_fill_manual(values=colors) + labs(fill="Phylogroup")


ggsave(tp, filename = "cladogram.png", width=10, height = 10, dpi = 300)

# use this tree to figure out where to root the tree after droping cryptic clades for clarity
ggtree(tree) %<+% labdf + 
  geom_tiplab(aes(label=label2), size = 2, offset=0, align=T) + #, color=tmatch)) +
  geom_tippoint(aes(shape=match, color=quardiplex_pcr_base_tidy), alpha=0.9)+
  scale_color_manual(values = colors) + 
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 
plot(  tree)
plot(   ape::drop.tip(tree, cryptics))
nodelabels()

cryptics <- c("GCA_000190955", "GCA_001660175", "GCA_002110245", "GCA_002109985")
(recttp <- ggtree(ape::root.phylo( 
  ape::drop.tip(tree, cryptics),
  node=137)) %<+% labdf + 
       geom_tiplab(aes(label=label2), size = 2, offset=.0, align=T) + #, color=tmatch)) +
    #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
     geom_tippoint(aes(shape=match, color=quardiplex_pcr_base_tidy), alpha=0.9) + 
    scale_color_manual(values = colors) + 
    scale_size_area(guide=F)  +
    scale_shape_discrete(guide=F) +
    labs(title="", color="Phylogroup" ) +
  geom_treescale(y=-2, fontsize=3, linesize=1.5, offset=1) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
    # 
    # labs(title="Parsimony cladogram of strains from Clermont, et al 2015", color="Phylogroup", subtitle="Tree generated with kSNP3 (k=19). Enlarged circular tips show where EzClermont \ndiffered from reported phylogroup (EzClermont type show in brackets)." ) +
    theme(legend.position="right") 
  
)
ggsave(filename = "tree.png", width=9, height = 9,dpi = 300,

gheatmap(recttp, results_table %>% 
           transform(Accession=substring(first = 1, last = 13, text = Accession)) %>% 
           select(-Strain) %>%
           column_to_rownames("Accession") ,
           color="black", 
         colnames_position="top", 
         offset=.004,
         colnames_angle=0, colnames_offset_y = 1, 
         width=.2,
         hjust=.5, font.size=2
         ) +  scale_fill_manual(values=colors, guide="none")

)

#####################################################################################
altlabdf <- metadata %>% rename("PRJNA230969"="accession","PRJNA321606"="old_accession") %>% 
  gather(key = "acc", value = "lab", -strain, -simple) %>% filter(lab != "") %>% 
  mutate(label2=paste0(strain," (", acc, ")")) %>% filter(startsWith(strain, "ECOR")) %>%  select(lab, strain, label2, acc) 
combined_tree <- read.tree("./validate/alignment_combined/parsnp.tree")
combined_tree$tip.label <- substr(combined_tree$tip.label, 1, 13)
altlabdf$lab <- substr(altlabdf$lab, 1, 13)

combined_tree <- keep.tip(combined_tree, altlabdf$lab[altlabdf$lab %in%combined_tree$tip.label])

       
(bioprojA <- ggtree(combined_tree) %<+% altlabdf + 
  geom_tiplab(aes(label=" "), size = 2, offset=.001, align=T) +
  geom_tippoint(aes(shape=acc, color=!strain %in% paste0("ECOR", c( 7, 20, 21, 23, 71, 72, 32, 33, 43, 37, 46, 39)))) +
  scale_color_manual(values = c("orange", "forestgreen")) + 
  labs(title="", color="Confirmed Identity", shape="BioProject") +
  guides(colour = guide_legend(override.aes = list(size=3)),
         shape = guide_legend(override.aes = list(size=3))) +
  theme(legend.position="bottom")   + coord_cartesian(clip = 'off') +
  geom_treescale(y=-3, fontsize=3, linesize=1.5, offset=.5) 
)
(bioprojB <-   ggtree(combined_tree, branch.length = "none") %<+% altlabdf + 
  geom_tiplab(aes(label=strain), size = 2, offset=-1, hjust = 1) +
  geom_tippoint(aes(shape=acc, color=!strain %in% paste0("ECOR", c( 7, 20, 21, 23, 71, 72, 32, 33, 43, 37, 46, 39)))) +
  scale_color_manual(values = c("orange", "forestgreen")) + 
  labs(title="", color="Confirmed Identity", shape="BioProject") +
  guides(colour = guide_legend(override.aes = list(size=3)),
         shape = guide_legend(override.aes = list(size=3))) +
  geom_treescale(y=-3, x=-15, fontsize=3, linesize=1.5, offset=.5) +
  coord_cartesian(clip = 'off') +
theme(legend.position="bottom")  + scale_x_reverse()
)

layout <- "
AAAAAABBBBB
AAAAAABBBBB
AAAAAABBBBB
AAAAAABBBBB
AAAAAABBBBB
AAAAAABBBBB
AAAAAABBBBB
CCCCCCCCCCC
"
combined_plot <- bioprojA + bioprojB  +guide_area() +plot_layout(guides = 'collect', design = layout) + plot_annotation(tag_levels = 'A') 

ggsave(combined_plot, filename = "~/GitHub/ezclermont/docs/analysis/bioproject_comparison.png",
       width=9, height = 12,dpi = 300,
      )
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
  # geom_hilight(node=146, fill=colors['A'], alpha=.6) + # A
  # geom_hilight(node=182, fill=colors['B1'], alpha=.6) + # B1
  # geom_hilight(node=224, fill=colors['B2'], alpha=.6) + # B2
  # geom_hilight(node=204, fill=colors['C'], alpha=.6) + # C
  # geom_hilight(node=211, fill=colors['D'], alpha=.6) + # D
  #geom_hilight(node=186, fill=colors['Clade II'], alpha=.6) + # clades are made the root
  # geom_hilight(node=205, fill=colors['E'], alpha=.6) + # E
  # geom_hilight(node=216, fill=colors['F'], alpha=.6) + # F
  # geom_hilight(node=137, fill=colors['G'], alpha=.6) + # G 
  # geom_tippoint(aes(x=x, size=matchsize , shape=match, color=quardiplex_pcr_base_tidy), alpha=1) + 
  # geom_tippoint(aes(x=x+3, size=matchsize , shape=match, color=ct_phylotype), alpha=0.7) + 
  # geom_tippoint(aes(x=x+6, size=matchsize , shape=match,  color=ezClermont_phylogroup), alpha=0.4) + 

#####################################################################################





##  old
(tp<- ggtree(tree2, layout="circular", branch.length = 'none') %<+% labdf + 
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

expdata_ct_default <- read.csv("validate/default.tab", 
                       sep="\t", col.names = c("accession", "alleles", "profile", "al", "ct_phylotype", "mashpath"), 
                       stringsAsFactors = F, header=F) %>%
  select(accession, ct_phylotype) %>%
  transform(accession= gsub(".scaffold.fasta", "", gsub("(.*?\\..)_(.*)", "\\1", accession)))

expdata_ct_500 <- read.csv("validate/500bp.tab", 
                               sep="\t", col.names = c("accession", "alleles", "profile", "al", "ct_phylotype_500", "mashpath"), 
                               stringsAsFactors = F, header=F) %>%
  select(accession, ct_phylotype_500) %>%
  transform(accession= gsub(".scaffold.fasta", "", gsub("(.*?\\..)_(.*)", "\\1", accession)))
merge(expdata_ct_default, expdata_ct_500)
