# R-script to generate CDR3 sequences with OLGA (function olga-generate_sequences)
# OLGA has to be pre-installed
# Only one of the generated datasets: N=5 and D=10,000 is used for the downstream 
# analysis
get_irr <- function(N, D, species, seed, path) {
    d <- c()
    
    for(i in 1:N) {
        fs <- paste0(path, species, "_TRA_", i)
        cmd <- paste("olga-generate_sequences", 
                      "--seed", seed+100+i, 
                      paste0("--", species, "TRA"),
                      "--num_seqs", D, 
                      "--seq_type amino_acid",
                      "--outfile", fs,
                     sep =  " ")
        system(cmd)
        a <- read.csv(file = fs, sep = "\t", header = F)
        colnames(a) <- c("CDR3a", "Va", "Ja")
        a$seeda <- seed+100+i
        
        # quantify pgen
        fsp <- paste0(path, "pgen_", species, "_TRA_", i)
        cmd <- paste("olga-compute_pgen", 
                     "--infile", fs,
                     paste0("--", species, "TRA"),
                     "--outfile", fsp,
                     "--seq_in", 0,
                     "--v_mask_index", 1,
                     "--j_mask_index", 2,
                     sep = " ")
        system(cmd)
        a_pgen <- read.csv(file = fsp, sep = "\t", header = F)
        a$pgena <- a_pgen[,2]
        file.remove(fs, fsp)
        
        fs <- paste0(path, species, "_TRB_", i)
        cmd <- paste("olga-generate_sequences", 
                      "--seed", seed+500+i, 
                      paste0("--", species, "TRB"),
                      "--num_seqs", D, 
                      "--seq_type amino_acid",
                      "--outfile", fs,
                     sep = " ")
        system(cmd)
        b <- read.csv(file = fs, sep = "\t", header = F)
        colnames(b) <- c("CDR3b", "Vb", "Jb")
        b$seedb <- seed+500+i
        
        # quantify pgen
        fsp <- paste0(path, "pgen_", species, "_TRB_", i)
        cmd <- paste("olga-compute_pgen", 
                     "--infile", fs,
                     paste0("--", species, "TRB"),
                     "--outfile", fsp,
                     "--seq_in", 0,
                     "--v_mask_index", 1,
                     "--j_mask_index", 2,
                     sep = " ")
        system(cmd)
        b_pgen <- read.csv(file = fsp, sep = "\t", header = F)
        b$pgenb <- b_pgen[,2]
        file.remove(fs, fsp)
        
        
        # cbind
        ab <- cbind(a,b)
        ab$species <- species
        ab$subject <- i
        ab$rep <- paste0(species, "_", i)
        d <- rbind(d, ab)
    }
    d$depth <- D
    d$N <- N
    return(d)
}

dir.create(path = "data")
Ns <- c(1, 3, 5, 10)
Ds <- c(500, 1000, 5000, 10000)

ds <- c()
for(n in Ns) {
    for(d in Ds) {
        hs <- get_irr(N = n, D = d, species = "human", seed = sample(x = 1:10^7, size = 1), path = "data/")
        mm <- get_irr(N = n, D = d, species = "mouse", seed = sample(x = 1:10^7, size = 1), path = "data/")
        ds <- rbind(ds, rbind(hs, mm))
    }
}
rm(hs, mm, d, n)

ds$clone_size <- 1
ds$id <- paste0(ds$N, '|', ds$depth, '|', ds$rep)

# no duplicates at all
d <- aggregate(clone_size~Va+Ja+CDR3a+Vb+Jb+CDR3b+id, data = ds, FUN = sum)
save(d, file = "Dataset3.RData", compress = TRUE)
