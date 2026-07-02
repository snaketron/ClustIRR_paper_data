library(scBubbletree)
library(data.table)

#### Process GEX ####
meta <- read.csv(file = "data/cell_metadata.csv", sep = ",")
meta$row <- 1:nrow(meta)
i <- which(meta$sample %in% c("Control_T_1", "Control_T_2", "Control_T_3", "EBV", "MART1"))
meta <- meta[i, ]
meta$row_new <- 1:nrow(meta)

# split count_matrix.mtx -l 1000000000 -> this will create 3 files: xaa, xab, xac
x1 <- fread(input = "data/GEX/xaa", sep = " ", header = FALSE, skip = 3,
            colClasses = c("integer", "integer", "integer"))
y <- x1$V1
y <- which(y %in% i)
x1 <- x1[y, ]
rm(y)
gc()

x2 <- fread(input = "data/GEX/xab", sep = " ", header = FALSE, skip = 0,
            colClasses = c("integer", "integer", "integer"))
y <- x2$V1
y <- which(y %in% i)
x2 <- x2[y, ]
rm(y)
gc()

x3 <- fread(input = "data/GEX/xac", sep = " ", header = FALSE, skip = 0,
            colClasses = c("integer", "integer", "integer"))
y <- x3$V1
y <- which(y %in% i)
x3 <- x3[y, ]
rm(y)
gc()

x <- rbind(x1, x2, x3)
rm(x1, x2, x3)
gc()



# write new .mtx
write.table(x = x, file = "data/GEX/count_matrix_sub.mtx", sep = " ", 
            row.names = F, col.names = F, quote = F)



# we need new cell IDs
u <- merge(x = x, y = meta[,c("row", "row_new")], by.x = "V1", by.y = "row", all.x = TRUE)
u$V1 <- u$row_new
u$row_new <- NULL
write.table(x = u, file = "data/GEX/count_matrix_sub2.mtx", sep = " ", 
            row.names = F, col.names = F, quote = F)
write.table(x = meta, file = "data/GEX/cell_metadata2.csv", sep = ",", 
            row.names = F, col.names = T, quote = F)




#### now create Seurat object #####
library(scBubbletree)
library(data.table)
library(Seurat)
library(Matrix)


features <- "data/GEX/all_genes.csv"

m <- ReadParseBio(data.dir = "data/GEX/processed/")

# Check to see if empty gene names are present, add name if so.
table(rownames(m) == "")
rownames(m)[rownames(m) == ""] <- "unknown"

# Read in cell meta data
cell_meta <- read.csv("data/GEX/processed/cell_metadata.csv", row.names = 1)

# Create object
d <- CreateSeuratObject(counts = m, 
                        min.features = 100, 
                        min.cells = 100,
                        names.field = 0, 
                        meta.data = cell_meta)

# Setting our initial cell class to a single type, this will changer after clustering. 
d@meta.data$orig.ident <- factor(rep("d", nrow(d@meta.data)))
Idents(d) <- d@meta.data$orig.ident
save(d, file = "data/GEX/processed/d_before_qc.RData", compress = TRUE)


d <- get(load("data/GEX/processed/d_before_qc.RData"))
d@meta.data$Barcode <- rownames(d@meta.data)
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "^MT-")
# plot <- VlnPlot(d, pt.size = 0.10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# get metadata to filter cells without TCR now
meta <- read.csv(file = "data/cell_metadata.csv", sep = ",")
br <- read.csv(file = "data/barcode_report.tsv", sep = "\t")
br <- br[which(br$isMultiplet == 0), ]
br <- br[which(is.na(br$clonotype_id)==FALSE),] # no TRAV TRBV
br <- br[-which(regexpr(pattern = "\\*", text = br$TRA_cdr3_aa)!=-1|regexpr(pattern = "\\*", text = br$TRB_cdr3_aa)!=-1),]
br <- merge(x = br, y = meta[,c("bc_wells", "sample")], by.x = "Barcode", by.y = "bc_wells", all.x = T)

# Perform the filtering
d <- subset(d, subset = nFeature_RNA < 5000 & nCount_RNA < 20000 & percent.mt < 15 
            & sample %in% c("Control_T_1", "Control_T_2", "Control_T_3", "EBV", "MART1")
            & Barcode %in% br$Barcode)

rm(meta)

d <- NormalizeData(d, normalization.method = "LogNormalize", scale.factor = 10000)

d <- FindVariableFeatures(d, selection.method = "vst", nfeatures = 3000)

d <- ScaleData(d, vars.to.regress = "percent.mt")

d <- RunPCA(d)

d <- CellCycleScoring(object = d, s.features = cc.genes$s.genes, 
                      g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

save(d, file = "data/GEX/processed/d_pca.RData", compress = TRUE)


#### update meta #####
d <- get(load("data/GEX/processed/d_pca.RData"))

meta <- read.csv(file = "data/cell_metadata.csv", sep = ",")
br <- read.csv(file = "data/barcode_report.tsv", sep = "\t")
br <- br[which(br$isMultiplet == 0), ]
br <- br[which(is.na(br$clonotype_id)==FALSE),] # no TRAV TRBV
br <- br[-which(regexpr(pattern = "\\*", text = br$TRA_cdr3_aa)!=-1|regexpr(pattern = "\\*", text = br$TRB_cdr3_aa)!=-1),]
br <- merge(x = br, y = meta[,c("bc_wells", "sample")], by.x = "Barcode", by.y = "bc_wells", all.x = T)
br <- br[,c("Barcode", "clonotype_id", "TRA_cdr3_aa", "TRB_cdr3_aa", 
           "TRA_V", "TRB_V", "TRA_J", "TRB_J")]

d@meta.data$i <- 1:nrow(d@meta.data)
m <- d@meta.data
m <- merge(x = m, y = br, by = "Barcode")
m <- m[order(m$i, decreasing = F),]
table(m$sample)

# merge community ID
m$sample <- gsub(pattern = "Control_T_", replacement = 'C', x = m$sample)
m$sample <- gsub(pattern = "EBV", replacement = 'E', x = m$sample)
m$sample <- gsub(pattern = "MART1", replacement = 'M', x = m$sample)

d@meta.data <- m
save(d, file = "data/GEX/processed/d_final.RData", compress = TRUE)
rm(m, b)



#### extract data for bubbletree analysis ####
library(Seurat)
d <- get(load("data/GEX/processed/d_final.RData"))

A <- d@reductions$pca@cell.embeddings[, 1:25]
m <- d@meta.data

markers <- read.csv(file = "data/GEX/processed/tcell_markers.csv", sep = ";")
tcell_symbols <- unlist(strsplit(x = markers$Symbol,  split = ', '))

tcell_activation <- read.csv(file = "data/GO/go_t_cell_activation.tsv", sep = "\t")
tcell_activation <- unique(tcell_activation$SYMBOL)

tcell_chemotaxis <- read.csv(file = "data/GO/go_t_cell_chemotaxis.tsv", sep = "\t")
tcell_chemotaxis <- unique(tcell_chemotaxis$SYMBOL)

oxphos <- read.csv(file = "data/GO/go_oxidative_phosphorylation.tsv", sep = "\t")
oxphos <- unique(oxphos$SYMBOL)

general_symbols <- c("IL7R", "CCR7", "S100A4", "CD14", "LYZ", "MS4A1", 
                     "CD8A", "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", 
                     "CST3", "PPBP", "IGHA2", "ZNF385D", "Ki-67")

# naive vs. memory
nm_symbols <- c("SELL", "CCR7", "IL7R", "TCF7", "CD27", "CD95")

genes_activated_t_cells <- c("CD69", "CD44", "SELL", "CD38","HLA-DRA", "CTLA4", "PDCD1", 
                             "GZMB", "PRF1", "IFNG", "IL2RA", "IL2RB", "MKI67")

tcr_symbols <- rownames(d@assays$RNA$data)[which(regexpr(
  pattern = "TRAV|TRAC|TRAJ|TRBV|TRBJ|TRBC", 
  text = rownames(d@assays$RNA$data))!=-1)]


# empirical evidence from
# https://www.nature.com/articles/s41467-019-12464-3#Sec27
require(readxl)
meta_1 <-readxl::read_xlsx("data/41467_2019_12464_MOESM9_ESM.xlsx", sheet = 1)
meta_1 <- meta_1[, c("gene", "cluster10", "cluster11")]
meta_1 <- meta_1[which(apply(X = meta_1[,-1], MARGIN = 1, FUN = function(x) {
  return(sum(x>0))
})>=1),]
meta_2 <-readxl::read_xlsx("data/41467_2019_12464_MOESM9_ESM.xlsx", sheet = 2)
meta_2 <- meta_2[, c("gene", "cluster8", "cluster10")]
meta_2 <- meta_2[which(apply(X = meta_2[,-1], MARGIN = 1, FUN = function(x) {
  return(sum(x>0))
})>=1),]
t_act_szabo <- unique(c(meta_1$gene, meta_2$gene))
rm(meta_1, meta_2)

# GO activated ab
t_ab_act <- read.csv("data/GO/QuickGO-annotations-1762437677769-20251106.tsv", sep = "\t")$SYMBOL



e <- t(as.matrix(d@assays$RNA$data[
  rownames(d@assays$RNA$data) %in% c(tcell_symbols,
                                     general_symbols, 
                                     genes_activated_t_cells,
                                     tcell_activation,
                                     tcell_chemotaxis,
                                     nm_symbols, 
                                     tcr_symbols,
                                     t_act_szabo,
                                     t_ab_act,
                                     oxphos), ]))

symbols <- list(tca = tcell_activation,
                tcc = tcell_chemotaxis,
                ts = tcell_symbols,
                gs = general_symbols,
                t_act = genes_activated_t_cells,
                t_act_szabo = t_act_szabo,
                t_ab_act = t_ab_act,
                oxphos = oxphos,
                ns = nm_symbols)

btd_input <- list(A = A, m = m, e = e, symbols = symbols)
save(btd_input, file = "data/GEX/processed/btd_input.RData", compress = TRUE)



#### bubbletree clustering resolution analysis ####
library(scBubbletree)

btd_input <- get(load("data/GEX/processed/btd_input.RData"))
A <- btd_input$A
m <- btd_input$m
e <- btd_input$e
symbols <- btd_input$symbols

br <- get_r(B_gap = 20,
             rs = 10^seq(from = -2, to = 1, by = 0.2),
             x = A,
             n_start = 10,
             iter_max = 50,
             algorithm = "original",
             knn_k = 30,
             cores = 30)

save(br, file = "data/GEX/processed/br.RData")
