library(ClustIRR)
# make sure blast is installed: rBLAST depends on it
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/home/simo/miniconda3/bin/", sep=.Platform$path.sep))

# data
cs <- get(load(file = "Dataset_1.RData"))

# select samples
cs$sample <- gsub(pattern = "Control_T_", replacement = 'C', x = cs$sample)
cs <- cs[cs$sample %in% c("C1", "C2", "C3", "EBV", "MART1"), ]
colnames(cs)[3:4] <- c("CDR3a", "CDR3b")

# clustirr analysis
c <- ClustIRR::clustirr(s = cs[, c("CDR3a", "CDR3b", "clone_size", "sample")],
                        meta = cs,
                        control = list(blast_gmi = 0.8,
                                       blast_cores = 40))
save(c, file = "c.RData", compress = TRUE)

gcd <- ClustIRR::detect_communities(c$graph,
                                    algorithm = "leiden", 
                                    metric = "average",
                                    resolution = 1,
                                    iterations = 1000,
                                    chains = c("CDR3a", "CDR3b"))
save(gcd, file = "gcd.RData", compress = TRUE)

dco <- ClustIRR::dco(community_occupancy_matrix = gcd$community_occupancy_matrix,
                   mcmc_control = list(mcmc_chains = 4,
                                       mcmc_cores = 4,
                                       mcmc_warmup = 750,
                                       mcmc_iter = 1750,
                                       adapt_delta = 0.9,
                                       max_treedepth = 10))
save(dco, file = "dco.RData", compress = TRUE)

ps <- dco$posterior_summary
save(ps, file = "ps.RData", compress = T)