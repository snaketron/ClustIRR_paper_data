# read individual repertoires, remove TCRs with unresolved V/J genes, keep only productive TCRs
d <- lapply(X = list.files("Dataset_2_raw/", pattern = "\\.tsv", full.names = TRUE), 
       FUN = function(x) {
    o <- read.csv(file = x,sep = "\t")
    o <- o[o$sequenceStatus=="In",]
    o <- o[o$vGeneName != "unresolved" ,]
    o <- o[o$jGeneName != "unresolved" ,]
    o$file <- x
    return(o)
})

# format the reperoires
d <- do.call(rbind, d)
d <- d[,c("aminoAcid", "count..templates.reads.", "file", "vGeneName", "jGeneName")]
colnames(d) <- c("CDR3b", "clone_size", "sample", "V", "J")
d$sample <- gsub(pattern = "data\\/\\/|\\.tsv", replacement = '', x = d$sample)
d <- d[!d$sample %in% c("Pt4_CR_PEP15_T-CELLS", "Pt4_CR_PEP16_T-CELLS"),]
d$sample[which(regexpr(pattern = "TumorA|TumorB|Tumor", text = d$sample)!=-1)] <- 
    paste0(d$sample[which(regexpr(pattern = "TumorA|TumorB|Tumor", text = d$sample)!=-1)], "_DayH")

# no V or J -> remove
d <- d[-which(d$V==""|d$J==""),]

m <- do.call(rbind, strsplit(x = d$sample, split = "_"))
d$patient <- m[,1]
d$outcome <- m[,2]
d$tissue <- m[,3]
d$timepoint <- m[,4]
rm(m)

# identical V+J+CDR3 -> clonotype
d <- aggregate(clone_size~., data = d, FUN = sum)

p <- d[d$tissue == "PBMC",]
t <- d[d$tissue != "PBMC",]

p$key <- paste0(p$patient, p$CDR3b, p$V, p$J)
t$key <- paste0(t$patient, t$CDR3b, t$V, t$J)

p$tumor <- ifelse(test = p$key %in% t$key, yes = "T+", no = "T-")
p$key <- NULL

d <- p

save(d, file = "Dataset_2.RData", compress = TRUE)
    