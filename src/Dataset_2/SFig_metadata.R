fs <- list.files(path = "Dataset_2_raw/", pattern = "", full.names = T)
fs <- fs[-which(regexpr(pattern = "_NE_", text = fs)!=-1)]
fs <- fs[which(regexpr(pattern = "_PBMC_", text = fs)!=-1)]
pt <- gsub(pattern = "paper\\_RT\\/data\\/\\/|\\.tsv", replacement = '', x = fs)
pt <- do.call(rbind, strsplit(x = pt, split = "\\_"))[,c(1,2,4)]

meta <- data.frame(pt)
colnames(meta) <- c("patient", "response", "time")
meta$time <- factor(meta$time, levels = c("Day0", "Day22", "Day43", "Day64",
                                          "Day88", "Day172", "Day203"))
meta$tumor <- "-"
meta$tumor[meta$patient %in% c("Pt4", "Pt32", "Pt36", "Pt38")] <- "A"
meta$tumor[meta$patient %in% c("Pt4")] <- "A,B"
meta$response <- factor(x = meta$response, levels = c("CR", "PR", "SD", "PD"))
meta$patient_num <- gsub(pattern = "Pt", replacement = '', x = meta$patient)
meta$patient_num <- as.numeric(meta$patient_num)
sort(unique(meta$patient_num))
meta$patient <- factor(x = meta$patient, levels = paste0("Pt", sort(decreasing = T, unique(meta$patient_num))))

gmeta <- ggplot(data = meta)+
    geom_tile(aes(x = time, y = patient), col = "white", fill = "steelblue")+
    geom_tile(aes(x = "tumor", y = patient), fill = "white", col = "gray")+
    geom_text(data = meta[duplicated(meta$patient)==FALSE, ],
              aes(x = "tumor", y = patient, label = tumor), size = 2.5)+
    geom_tile(aes(x = "response", y = patient, fill = response), col = "white")+
    scale_fill_brewer(palette = "OrRd", type = "seq", direction = -1)+
    scale_x_discrete(name = "", expand = c(0,0))+
    theme_bw(base_size = 10)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

gmeta

ggsave(plot = gmeta,
       filename = "gmeta.pdf",
       device = "pdf",
       width = 3,
       height = 3.25)