###### BTD ######
.libPaths(new = "/mnt/nfs/simo/rpack_4_4/")
library(scBubbletree)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Seurat)

btd_input <- get(load("data/GEX/processed/btd_input.RData"))
A <- btd_input$A

btd <- get_bubbletree_graph(x = A,
                            r = 0.05,
                            algorithm = "original",
                            n_start = 20,
                            iter_max = 100,
                            knn_k = 30,
                            cores = 20,
                            B = 1000,
                            N_eff = 200,
                            round_digits = 1,
                            show_simple_count = FALSE)

btd$tree

save(btd, file = "data/GEX/processed/btd.RData", compress = TRUE)




###### heatmaps ######
.libPaths(new = "/mnt/nfs/simo/rpack_4_4/")
library(scBubbletree)
library(ggplot2)
library(patchwork)
library(ClustIRR)

# BTD/GEX data
btd_input <- get(load("data/GEX/processed/btd_input.RData"))
d <- get(load("data/GEX/processed/d_final.RData"))
btd <- get(load("data/GEX/processed/btd.RData"))
e <- btd_input$e
symbols <- btd_input$symbols
m <- btd_input$m

# GCD data
gcd <- get(load(file = "res/gcd.RData"))
ns <- gcd$node_summary
ns$sample <- gsub(pattern = "EBV", replacement = 'E', x = ns$sample)
ns$sample <- gsub(pattern = "MART1", replacement = 'M', x = ns$sample)

# betas
ps <- get(load(file = "res/ps.RData"))
beta <- ps$beta
beta <- beta[order(beta$mean, decreasing = T),]

# find EBV and MART1 specific clones
ns <- get_ag_species_hits(node_summary = ns, db = "vdjdb", ag_species = "EBV")$node_summary
coms_e <- unique(ns$community[ns$EBV_CDR3a!=0|ns$EBV_CDR3b!=0])
coms_e <- beta$community[beta$sample == "EBV" & beta$community %in% coms_e][1:5]

ns <- get_ag_gene_hits(node_summary = ns, db = "vdjdb", ag_gene = c("MLANA"))$node_summary
coms_m <- unique(ns$community[ns$MLANA_CDR3a!=0|ns$MLANA_CDR3b!=0])
coms_m <- beta$community[beta$sample == "MART1" & beta$community %in% coms_m][1:5]

ns$E <- ifelse(test = ns$EBV_CDR3a!=0|ns$EBV_CDR3b!=0, yes = '+', no = '-')
ns$M <- ifelse(test = ns$MLANA_CDR3a!=0|ns$MLANA_CDR3b!=0, yes = '+', no = '-')

x <- data.frame(community = c(coms_e, coms_m), community_id = c(paste0("e",1:5), paste0("m",1:5)))
rm(coms_e, coms_m)

ns_sel <- ns[ns$community %in% x$community, ]
ns_sel$cells <- ns_sel$clone_size
ns_sel$clones <- 1
ns_sel <- merge(x = ns_sel, y = x, by = "community", all.x = TRUE)
ns_sel$community_id <- factor(x = ns_sel$community_id, levels = x$community_id)


g_s_v <- get_cat_tiles(btd = btd,
                       f = m$sample,
                       integrate_vertical = TRUE,
                       rotate_x_axis_labels = TRUE,
                       round_digits = 1,
                       tile_text_size = 2.25)

g_s_h <- get_cat_tiles(btd = btd,
                       f = m$sample,
                       integrate_vertical = FALSE,
                       rotate_x_axis_labels = TRUE,
                       round_digits = 1,
                       tile_text_size = 2.25)

# add cell cycle markers
g_phase <- get_cat_tiles(btd = btd, 
                     f = d@meta.data$Phase,
                     integrate_vertical = F,
                     round_digits = 1, 
                     tile_text_size = 2.25,
                     rotate_x_axis_labels = T)

# add T-cell activation
e_tca <- e[,colnames(e) %in% symbols$t_act_szabo]
e_tca <- data.frame(e_tca)
e_tca <- rowSums(e_tca)
e_tca <- data.frame(e = e_tca, sample = d@meta.data$sample, clonotype_id = d@meta.data$clonotype_id)
e_tca$i <- 1:nrow(e_tca)

g_tca <- get_num_violins(btd = btd, fs = e_tca$e, 
                         x_axis_name = 'GEX', 
                         rotate_x_axis_labels = F)

ns_sel <- aggregate(cells~clonotype_id+sample+community_id+E+M, data = ns_sel, FUN = sum)
e_tca <- merge(x = e_tca, 
               y = ns_sel[,c("clonotype_id",  "sample", "community_id", "E", "M", "cells")], 
               by = c("clonotype_id", "sample"), all.x = TRUE)
rm(ns_sel)



g_top <- (btd$tree|g_s_v$plot|g_s_h$plot|g_tca$plot|g_phase$plot)+
  patchwork::plot_layout(widths = c(1.4, 1.6, 1.6, 1.15, 1.1))



ggsave(plot = g_top,
       filename = "manuscript/Fig_2.pdf",
       device = "pdf",
       width = 9,
       height = 2.4)

rm(btd, btd_input, d, ds, e, e_as, e_tca, e_tcc, g_dac, g_general, g_phase, g_s_h, 
   g_s_v, g_tca, g_tcc, A, bs, gcd, m, ns, symbols, top_dacs, mart1_dacs, ebv_dacs,
   mart_dacs)
