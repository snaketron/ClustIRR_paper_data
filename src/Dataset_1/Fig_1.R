#### Fig 1A-C ####
library(ggplot2)
ggplot2::theme_set(new = theme_bw(base_size = 10))
library(ggforce)
library(ggrepel)
library(patchwork)
library(dplyr)
library(ClustIRR)

ps <- get(load(file = "ps.RData"))
bs <- ps$beta
bs$sample[bs$sample=="MART1"] <- "M"
bs$sample[bs$sample=="EBV"] <- "E"
bs <- bs %>% group_by(sample) %>% mutate(rank = rank(-mean, ties.method = "first")) %>% ungroup()

gcd <- get(load(file = "gcd.RData"))
ns <- gcd$node_summary
ns$sample[ns$sample=="MART1"] <- "M"
ns$sample[ns$sample=="EBV"] <- "E"

ns <- merge(x = ns, y = bs[,c("sample", "community", "mean", "rank")], 
            by = c("sample", "community"), all.x = TRUE)

# meta
m <- ns[regexpr(pattern = "EBV", text = ns$Ag_species)!=-1 & ns$sample %in% c("E"),]
m <- m[duplicated(m[, c("rank", "sample")])==F,]
m <- m[order(m$rank, decreasing = F),]
m <- m[1:5,]
m$cid <- paste0("e", 1:nrow(m))

v <- get_beta_violin_ag(beta = bs, node_summary = ns, ag = "EBV", ag_species = TRUE, db = "vdjdb")
Fig_1B <- v+
  ggtitle(label = '', subtitle = '')+
  scale_color_manual(name = "EBV-specific", values = c("orange", "steelblue"))+
  scale_radius(name = "EBV-specific\ncells", breaks = scales::pretty_breaks(4), range = c(0.6, 4))+
  guides(size = guide_legend(order = 2), colour = guide_legend(order = 1))+
  scale_y_continuous(limits = c(-2, 6.5))+
  xlab(label = "TCR repertoire")+
  ggrepel::geom_text_repel(data = m[m$sample == "E",],
                           aes(x = sample, y = mean, label = cid),
                           min.segment.length = 0, size = 2.5)
  
# meta
m <- ns[regexpr(pattern = "MLANA", text = ns$Ag_gene)!=-1 & ns$sample %in% c("M"),]
m <- m[duplicated(m[, c("rank", "sample")])==F,]
m <- m[order(m$rank, decreasing = F),]
m <- m[1:5,]
m$cid <- paste0("m", 1:nrow(m))

v <- get_beta_violin_ag(beta = bs, node_summary = ns, ag = c("MLANA"), 
                     ag_species = FALSE, db = "vdjdb")
Fig_1C <- v+
  ggtitle(label = '', subtitle = '')+
  scale_color_manual(name = "MART1-specific", values = c("orange", "purple"))+
  scale_radius(name = "MART1-specific\ncells", breaks = scales::pretty_breaks(4), range = c(0.6, 4))+
  guides(size = guide_legend(order = 2), colour = guide_legend(order = 1))+
  scale_y_continuous(limits = c(-2, 6.5))+
  xlab(label = "TCR repertoire")+
  ggrepel::geom_text_repel(data = m[m$sample == "M",],
                           aes(x = sample, y = mean, label = cid),
                           min.segment.length = 0, size = 2.5)


get_cs <- function(com) {
  ns <- colnames(com)
  ns[ns=="MART1"] <- "M"
  ns[ns=="EBV"] <- "E"
  
  cos_sim <- matrix(data = 0, nrow = ncol(com), ncol = ncol(com))
  for(i in 1:ncol(com)) {
    for(j in 1:ncol(com)) {
      v1 <- com[,i]
      v2 <- com[,j]
      cos_sim[i,j] <- sum(v1 * v2) / (sqrt(sum(v1^2)) * sqrt(sum(v2^2)))
    }
  }
  colnames(cos_sim) <- ns
  rownames(cos_sim) <- ns
  cos_sim <- reshape2::melt(cos_sim)
  colnames(cos_sim) <- c("i", "j", "CS")
  
  g <- ggplot(data = cos_sim)+
    geom_tile(aes(x = i, y = j, fill = CS), col = "white")+
    geom_text(aes(x = i, y = j, label = round(CS, digits = 2)), 
              size = 2.75, col = "black")+
    scale_fill_distiller(palette = "Spectral", 
                         limits = c(-0.01, 1.01),
                         breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                         labels = c(0, 0.25, 0.5, 0.75, 1.0))+
    theme(legend.position = "right")+
    guides(fill = guide_colourbar(barheight = 5, barwidth = 0.5))+
    xlab(label = '')+
    ylab(label = '')
  
  return(list(g = g, cs = cos_sim))
}

Fig_1A <- get_cs(com = gcd$community_occupancy_matrix)$g


#### Fig 1D ####
library(ggplot2)
ggplot2::theme_set(new = theme_bw(base_size = 10))
library(ggforce)
library(ggrepel)
library(patchwork)
library(dplyr)
library(ClustIRR)

gcd <- get(load(file = "gcd.RData"))
com <- gcd$community_occupancy_matrix
colnames(com)[colnames(com) == "EBV"] <- "E"
colnames(com)[colnames(com) == "MART1"] <- "M"

ps <- get(load(file = "ps.RData"))
bs <- ps$beta
bs$sample[bs$sample=="MART1"] <- "M"
bs$sample[bs$sample=="EBV"] <- "E"

get_beta_cprob_ag <- function(node_summary,
                              ag,
                              ag_species = TRUE,
                              db = "vdjdb") {
  
  getc <- function(x, h) {
    q <- h[h$sample == x, ]
    z <- cumsum(q$clone_size)
    p <- z/max(z)
    p[is.infinite(p)|is.nan(p)] <- 0
    return(data.frame(p = p, b = q$mean, sample = x))
  }
  
  if(ag_species) {
    h <- get_ag_species_hits(node_summary = node_summary, db = db, ag_species = ag)
  } else {
    h <- get_ag_gene_hits(node_summary = node_summary, db = db, ag_gene = ag)
  }
  
  h$node_summary$spec <- ifelse(test = apply(X = h$node_summary[,h$new_columns], 
                                             MARGIN=1, FUN=sum)==0, yes = FALSE, no = TRUE)
  h <- h$node_summary
  
  ha <- aggregate(clone_size~community+sample+mean+spec, data = h, FUN = sum)
  ha$clone_size[ha$spec == FALSE] <- 0
  ha <- ha[order(ha$mean, decreasing = TRUE),]
  ha <- do.call(rbind, lapply(X = unique(ha$sample), h = ha, FUN = getc))
  colnames(ha) <- c("p_ag", "b", "sample")
  
  hb <- aggregate(clone_size~community+sample+mean, data = h, FUN = sum)
  hb <- hb[order(hb$mean, decreasing = TRUE),]
  hb <- do.call(rbind, lapply(X = unique(hb$sample), h = hb, FUN = getc))
  colnames(hb) <- c("p_b", "b", "sample")
  
  hab <- merge(x = ha, y = hb, by = c("b", "sample"))
  hab <- hab[order(hab$b, decreasing = TRUE),]
  hab$Ag <- ag
  
  return(hab)
}

v <- rbind(get_beta_cprob_ag(node_summary = ns, ag = "MLANA", ag_species = F, db = "vdjdb"),
           get_beta_cprob_ag(node_summary = ns, ag = "EBV", ag_species = T, db = "vdjdb"))
v <- v[order(v$p_ag, decreasing = TRUE),]
v$Ag[v$Ag=="MLANA"] <- "MART1"

Fig_1D <- ggplot(data = v)+
  facet_wrap(~sample)+
  geom_line(aes(x = b, y = p_b, col = Ag),
            linetype = "dashed", size = .5, alpha = 0.7, col = "black")+
  geom_line(aes(x = b, y = p_ag, col = Ag), size = .5, alpha = 0.7)+
  xlab(label = expression(beta))+
  ylab(label = "Cumulative probability")+
  scale_color_manual(name = "Antigen",values = c("EBV"="steelblue", "MART1"="purple"))+
  theme(legend.position = "top", 
        strip.text.x = element_text(margin = margin(0.02,0,0.02,0, "cm")))


Fig_1 <- (Fig_1A|Fig_1B)/(Fig_1D|Fig_1C)

ggsave(filename = "Fig_1.pdf",
       plot = Fig_1,
       device = "pdf",
       width = 13,
       height = 3.5)

