library(ClustIRR)
library(parallel)
library(dplyr)

Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/home/simo/miniconda3/bin/", sep=.Platform$path.sep))

# These script contains 3 steps:
# 1. clustirr analysis
# 2. dco
# 3. parsing outputs

# step 1
d <- get(load(file = "d.RData"))
for(pt in unique(d$patient)) {
    dir.create(pt)

    if(file.exists(paste0(pt, "/gcd_light.RData"))==FALSE) {

        s <- d[d$patient == pt, c("CDR3b", "sample", "clone_size")]
        m <- d[d$patient == pt, ]

        # main clustirr analysis
        c <- clustirr(s = s, meta = m, control = list(blast_gmi = 0.8, blast_cores = 40, trim_flank_aa = 3))
        save(c, file = paste0(pt, "/c.RData"))
        cat("DONE!\n")

        # detect_communities
        gcd <- detect_communities(graph = c$graph, algorithm = "leiden", metric = "average",
                                  resolution = 1, iterations = 1000, chains = "CDR3b")
        save(gcd, file = paste0(pt, "/gcd.RData"), compress = TRUE)
        cat("DONE!\n")

        # graphs
        graph <- gcd$graph
        save(graph, file = paste0(pt, "/graph.RData"), compress = TRUE)
        cat("DONE!\n")

        gcd_light <- gcd
        gcd_light$graph <- NULL
        save(gcd_light, file = paste0(pt, "/gcd_light.RData"), compress = TRUE)
    }
}

# step 2
d <- get(load(file = "Dataset_2.RData"))
pts <- unique(d$patient)
rm(d);gc();gc()
dco <- mclapply(mc.cores = length(pts), X = pts, FUN = function(x) {
    if(file.exists(paste0(x, "/dco.RData"))==FALSE) {
        gcd_light <- get(load(file = paste0(x, "/gcd_light.RData")))
        com <- gcd_light$community_occupancy_matrix
        rm(gcd_light);gc();gc();gc();
        dco <- dco(community_occupancy_matrix = com,
                   mcmc_control = list(mcmc_chains = 4,
                                       mcmc_cores = 1,
                                       mcmc_warmup = 750,
                                       mcmc_iter = 1750,
                                       adapt_delta = 0.9,
                                       max_treedepth = 10),
                   compute_delta = FALSE,
                   groups = NA)
        
        save(dco, file = paste0(x, "/dco.RData"), compress=T)
    }
})


# step 3

parse_out_beta <- function(x) {
    cat(x, "\n")
    y <- gsub(pattern = "dco.RData", replacement = "beta.RData", x = x)
    dco <- get(load(x))
    beta <- dco$posterior_summary$beta
    save(beta, file = y, compress = T)
}

fs <- list.files(pattern = "dco\\.RData", recursive = T, full.names = T)
lapply(X = fs, FUN = parse_out_beta)


parse_out_node_summary <- function(x) {
    cat(x, "\n")
    y <- gsub(pattern = "gcd_light\\.RData", replacement = "ns.RData", x = x)
    gcd_light <- get(load(x))
    ns <- gcd_light$node_summary
    save(ns, file = y, compress = T)
}

fs <- list.files(path = "results/", pattern = "gcd_light\\.RData", recursive = T, full.names = T)
lapply(X = fs, FUN = parse_out_node_summary)



parse_out_delta <- function(x) {
    cat(x, "\n")
    fo <- gsub(pattern = "dco\\.RData", replacement = "delta\\.RData", x = x)
    dco <- get(load(x))

    cs <- colnames(dco$community_occupancy_matrix)
    j_22 <- which(regexpr(pattern = "Day22", text = cs)!=-1)
    j_0 <- which(regexpr(pattern = "Day0", text = cs)!=-1)


    beta <- rstan::extract(dco$fit, par = "beta")
    beta <- beta$beta[,c(j_0,j_22),]

    d <- vapply(X = 1:dim(beta)[3], d = beta,
                FUN.VALUE = numeric(4),
                FUN = function(x, d) {
                    y <- numeric(4)
                    z <- d[,2,x]-d[,1,x]
                    zhdi <- get_hdi(z, hdi_level = 0.95)
                    y[1] <- mean(z)
                    y[2] <- zhdi[1]
                    y[3] <- zhdi[2]
                    y[4] <- get_pmax(x = z)
                    return(y)
                })
    d <- t(d)
    delta <- data.frame(d)
    colnames(delta) <- c("mean", "L95", "H95", "pmax")


    save(delta, file = fo, compress = T)
}

fs <- list.files(path = "", pattern = "dco\\.RData", recursive = T, full.names = T)
lapply(X = fs, FUN = parse_out_delta)
