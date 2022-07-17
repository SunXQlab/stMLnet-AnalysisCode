#############
## library ##
#############

library(OmnipathR)
library(dplyr)
library(readr)

library(ggvenn)
library(ggpubr)
library(ggsci)

rm(list = ls())
gc()

setwd("./stMLnet/prior_knowledge/")

###########
## color ##
###########

scales::show_col(pal_locuszoom(palette = "default", alpha = 0.8)(7))
mycolors_locus <- pal_locuszoom(palette = "default", alpha = 0.8)(7)

mycolor_key <- mycolors_locus[1:4]
scales::show_col(mycolor_key)

##########
## load ##
##########
## load Ligand-Receptor Interaction Network ####

## stMLnet
load("./output/scMLnet.human.LigRec.rda")
lr_Network_stMLnet <- hm.LigRec
lr_Network_stMLnet_Unique <- lr_Network_stMLnet %>%
  dplyr::distinct(source, target) %>% 
  dplyr::mutate(key = paste(source, target, sep = '_'))

lig_stMLnet <- lr_Network_stMLnet_Unique %>% dplyr::select(source) %>% unlist() %>% unique()
rec_stMLnet <- lr_Network_stMLnet_Unique %>% dplyr::select(target) %>% unlist() %>% unique()

## NicheNet
load("../other_method/NicheNet/data/lr_network.rda")
lr_Network_Nichenet <- lr_network
lr_Network_Nichenet_Unique <- lr_Network_Nichenet %>% 
  dplyr::distinct(from, to) %>% 
  dplyr::rename(source=from, target=to) %>% 
  dplyr::mutate(key = paste(source, target, sep = '_'))

lig_Nichenet <- lr_Network_Nichenet_Unique %>% dplyr::select(source) %>% unlist() %>% unique()
rec_Nichenet <- lr_Network_Nichenet_Unique %>% dplyr::select(target) %>% unlist() %>% unique()

## Omnipath
lr_Network_Omnipath <- readRDS(file = "../other_method/Omnipath/cleaned_data/LigRec_from_interactions.rds")
lr_Network_Omnipath_Unique <- lr_Network_Omnipath %>%
  dplyr::filter(weight > 0) %>% 
  dplyr::mutate(key = paste(source, target, sep = '_'))

lig_Omnipath <- lr_Network_Omnipath_Unique %>% dplyr::select(source) %>% unlist() %>% unique()
rec_Omnipath <- lr_Network_Omnipath_Unique %>% dplyr::select(target) %>% unlist() %>% unique()

## load TF-TG Interaction Network ####

## stMLnet
load("./output/scMLnet.human.TFTG.rda")
tftg_Network_stMLnet <- hm.TFTG
tftg_Network_stMLnet_Unique <- tftg_Network_stMLnet %>%
  dplyr::distinct(source, target) %>% 
  dplyr::mutate(key = paste(source, target, sep = '_'))

tf_stMLnet <- tftg_Network_stMLnet_Unique %>% dplyr::select(source) %>% unlist() %>% unique()
tg_stMLnet <- tftg_Network_stMLnet_Unique %>% dplyr::select(target) %>% unlist() %>% unique()

## NicheNet

load("../other_method/NicheNet/data/gr_network.rda")
tftg_Network_Nichenet <- gr_network
tftg_Network_Nichenet_Unique <- tftg_Network_Nichenet %>% 
  dplyr::distinct(from, to) %>% 
  dplyr::rename(source=from, target=to) %>% 
  dplyr::mutate(key = paste(source, target, sep = '_'))

tf_Nichenet <- tftg_Network_Nichenet_Unique %>% dplyr::select(source) %>% unlist() %>% unique()
tg_Nichenet <- tftg_Network_Nichenet_Unique %>% dplyr::select(target) %>% unlist() %>% unique()

## Omnipath

tftg_Network_Omnipath <- readRDS(file = "../other_method/Omnipath/cleaned_data/TFTG_from_interactions.rds")
tftg_Network_Omnipath_Unique <- tftg_Network_Omnipath %>%
  dplyr::filter(weight>0) %>%
  dplyr::mutate(key = paste(source, target, sep = '_'))

tf_Omnipath <- tftg_Network_Omnipath_Unique %>% dplyr::select(source) %>% unlist() %>% unique()
tg_Omnipath <- tftg_Network_Omnipath_Unique %>% dplyr::select(target) %>% unlist() %>% unique()


## load Signaling Interaction Network ####

## stMLnet
load("./output/scMLnet.human.pw.rda")
signal_Network_stMLnet <- hm.pw
signal_Network_stMLnet_Unique <- signal_Network_stMLnet %>%
  dplyr::distinct(source, target) %>% 
  dplyr::mutate(key = paste(source, target, sep = '_'))

signal_stMLnet <- unique(c(signal_Network_stMLnet_Unique$source,signal_Network_stMLnet_Unique$target))

## NicheNet
load("../other_method/NicheNet/nichenet.human.ppi.rda")
signal_Network_Nichenet <- hm.ppi
signal_Network_Nichenet_Unique <- signal_Network_Nichenet %>% 
  dplyr::distinct(source, target) %>% 
  dplyr::mutate(key = paste(source, target, sep = '_'))

signal_Nichenet <- unique(c(signal_Network_Nichenet_Unique$source,signal_Network_Nichenet_Unique$target))

## Omnipath
signal_Network_Omnipath <- readRDS(file = "../other_method/Omnipath/cleaned_data/RecTF_from_interactions.rds")
signal_Network_Omnipath_Unique <- signal_Network_Omnipath %>%
  dplyr::filter(weight>0) %>%
  dplyr::mutate(key = paste(source, target, sep = '_'))

signal_Omnipath <- unique(c(signal_Network_Omnipath_Unique$source,signal_Network_Omnipath_Unique$target))

###########################################
## veen for pair from different dataases ##
###########################################

## LigRec

listinput <- list(Nichenet = lr_Network_Nichenet_Unique$key,
                  Omnipath = lr_Network_Omnipath_Unique$key,
                  stMLnet = lr_Network_stMLnet_Unique$key)

p <- venn.diagram(listinput, resolution = 1500, imagetype = "tiff",
                  main.cex = 2, main.pos = c(0.5, 1.2), cat.cex = 1.4, cex = 1.4,
                  fill=c("light blue", "light green", "red"), alpha=c(0.5,0.5,0.5),
                  main="Intersection of LigRec pairs \ncollected from different databases",
                  filename = NULL)
pdf("./figure/veen_LigRec_from_databases.pdf",width = 7,height = 7)
grid.draw(p)
dev.off()

## TFTG

listinput <- list(Nichenet = tftg_Network_Nichenet_Unique$key,
                  Omnipath = tftg_Network_Omnipath_Unique$key,
                  stMLnet = tftg_Network_stMLnet_Unique$key)

p <- venn.diagram(listinput, resolution = 1500, imagetype = "tiff",
                  main.cex = 2, main.pos = c(0.5, 1.2), cat.cex = 1.4, cex = 1.4,
                  fill=c("light blue", "light green", "red"), alpha=c(0.5,0.5,0.5),
                  main="Intersection of TFTG pairs \ncollected from different databases",
                  filename = NULL)
pdf("./figure/veen_TFTG_from_databases.pdf",width = 7,height = 7)
grid.draw(p)
dev.off()

## Signaling

listinput <- list(Nichenet = signal_Network_Nichenet_Unique$key,
                  Omnipath = signal_Network_Omnipath_Unique$key,
                  stMLnet = signal_Network_stMLnet_Unique$key)

p <- venn.diagram(listinput, resolution = 1500, imagetype = "tiff",
                  main.cex = 2, main.pos = c(0.5, 1.2), cat.cex = 1.4, cex = 1.4,
                  fill=c("light blue", "light green", "red"), alpha=c(0.5,0.5,0.5),
                  main="Intersection of Signaling pairs \ncollected from different databases",
                  filename = NULL)
pdf("./figure/veen_Signaling_from_databases.pdf",width = 7,height = 7)
grid.draw(p)
dev.off()
