library(ClustIRR)
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/home/simo/miniconda3/bin/", sep=.Platform$path.sep))

ds <- get(load("Dataset_3"))
ds$sample <- ds$rep
ds$dataset <- paste0("N", ds$N, "_D", ds$depth)
d <- ds[ds$dataset == "N5_D10000", ]
s_olga <- data.frame(CDR3a = d$CDR3a, CDR3b = d$CDR3b, clone_size = d$clone_size,
                     TRAV = d$Va, TRBV = d$Vb, TRAJ = d$Ja, TRBJ = d$Jb, sample = d$sample,
                     dataset = ifelse(d$species == "human", yes = "human_G", no = "mouse_G"),
                     pgena = d$pgena, pgenb = d$pgenb)
s_olga$sample <- gsub(pattern = "human_", replacement = "H", x = s_olga$sample)
s_olga$sample <- gsub(pattern = "mouse_", replacement = "M", x = s_olga$sample)


# 1.
c <- clustirr(s = s_olga[,c("sample", "clone_size", "CDR3a", "CDR3b")],
              meta = s_olga, control = list(blast_gmi = 0.8,
                                            blast_cores = 50))
save(c, file = "c.RData", compress = TRUE)

gcd <- detect_communities(graph = c$graph,
                          algorithm = "leiden",
                          metric = "average",
                          resolution = 1,
                          iterations = 1000,
                          chains = c("CDR3a", "CDR3b"))
save(gcd, file = "gcd.RData", compress = TRUE)


# prep for DCO
gcd <- get(load("gcd.RData"))
com <- gcd$community_occupancy_matrix
ns <- gcd$node_summary
ns <- ns[,c("sample", "dataset")]
ns <- ns[duplicated(ns)==F, ]
ns$dataset_num <- as.numeric(as.factor(ns$dataset))
meta <- c()
for(i in 1:ncol(com)) {
    meta <- rbind(meta, data.frame(dataset = ns$dataset[ns$sample == colnames(com)[i]],
                                   dataset_num = ns$dataset_num[ns$sample == colnames(com)[i]],
                                   sample = colnames(com)[i]))
}
rm(i, ns, gcd)
gc();gc();gc();gc();gc();gc()

# 3.
dco <- dco(community_occupancy_matrix = com,
           mcmc_control = list(mcmc_cores = 4,
                               mcmc_chains = 4,
                               mcmc_warmup = 750,
                               mcmc_iter = 1750,
                               max_treedepth = 10,
                               adapt_delta = 0.90),
           compute_delta=TRUE,
           groups = meta$dataset_num)
save(dco, file = "dco.RData", compress = TRUE)

# 4.
ps <- dco$posterior_summary
save(ps, file = "ps.RData", compress = TRUE)

