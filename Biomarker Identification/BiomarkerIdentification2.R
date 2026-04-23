library(Seurat)
library(pheatmap)

# =====================================================
# 1. READ EXPRESSION FILE
# =====================================================
expr <- read.csv(
  "C:/Users/C11_LAB/Downloads/Merged_AML_log2RPM_BatchRemoved2.csv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expr <- expr[!is.na(expr[,1]) & expr[,1] != "", ]

gene_names <- expr[,1]
expr_mat <- as.matrix(expr[,-1])
expr_mat <- apply(expr_mat, 2, as.numeric)
rownames(expr_mat) <- gene_names

# Remove duplicates
expr_mat <- rowsum(expr_mat, group = rownames(expr_mat))

# =====================================================
# 2. CREATE SEURAT OBJECT
# =====================================================
aml <- CreateSeuratObject(counts = expr_mat)
aml <- SetAssayData(aml, layer = "data", new.data = expr_mat)

# =====================================================
# 3. Find variable features
# =====================================================
aml <- FindVariableFeatures(
  aml,
  selection.method = "vst",
  nfeatures = 2000
)

# =====================================================
# 4. Scale data
# =====================================================
aml <- ScaleData(aml, verbose = FALSE)

# =====================================================
# 5. PCA
# =====================================================
aml <- RunPCA(aml, npcs = 5, verbose = FALSE)
# =====================================================
# 6. Clustering
# =====================================================
aml <- FindNeighbors(aml, dims = 1:5)
aml <- FindClusters(aml, resolution = 1.3)

# =====================================================
# 7. UMAP
# =====================================================
aml <- RunUMAP(aml, dims = 1:3, n.neighbors = 5)


# =====================================================
# 3. READ RESPONSE METADATA
# =====================================================
meta <- read.csv(
  "C:/Users/C11_LAB/Downloads/response.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(meta)[1:2] <- c("Sample","Response_Group")
meta$Sample <- as.character(meta$Sample)

meta <- meta[meta$Sample %in% colnames(aml), ]
rownames(meta) <- meta$Sample

aml <- AddMetaData(aml, meta)

# =====================================================
# 4. KEEP ONLY CR & NR (REMOVE NA)
# =====================================================
aml_sub <- subset(
  aml,
  subset = Response_Group %in% c("CR","NR")
)

# Check
table(aml_sub@meta.data$Response_Group)

# =====================================================
# 5. EXPRESSION MATRIX
# =====================================================
expr_all <- GetAssayData(aml_sub, layer = "data")
inflam_genes <- c(
  "IGLL1","CD79B","CD79A","IGLC1","PAX5","IGHM","IGKC","HLA-DRA","CD81","PSMB6",
  "HLA-DPB1","IGLC2","PSMB9","SYK","PSMB8","PSMA6","IGHD","PSMA4","PSMA7","HLA-DRB1",
  "CD38","EZR","CD44","TYROBP","FTL","LYZ","ADGRE5","PAFAH1B2","MGST1","CTSA","CTSD",
  "RHOG","CFD","ITGB2","ATP6V0C","PRAM1","ATP8B4","ADAM17","TMEM30A","RNASET2","KMT2E",
  "MANBA","AZU1","SDCBP","GLIPR1","CD63","PTPRC","SLC2A3","MPO","CD36","CTSG","CXCL8",
  "UBC","ANXA1","NFKBIA","IL1RAP","SQSTM1","HIF1A","CEBPB","TFPI","PELI1","IRAK3",
  "PSME4","CCL5","MAP3K8","KLF2","TRIB1","ZFP36","PTGER4","PLCG2","VIM","JUNB","FOS",
  "RIOK3","CTNNB1","PIK3R1","RGCC","LRRFIP1","STAT3","EIF2AK2","DDIT3","DDX21","RUNX1",
  "HGF","HIPK1","CYLD","ELF1","ADAR","PARP14","SOCS3","IFI6","IFI44L","DDIT4","IFITM2",
  "ISG15","CALCOCO2","RBPJ","MX1","LYST","HOXA9","VSIR","ID2","ATP2B1","HOMER3","SRGN",
  "B2M","CST3","COTL1","FCER1G","FCN1","FGL2","CTSS","MNDA","CYBA","S100A11","ARPC5",
  "S100A9","SLC11A1","CLEC12A","GSTP1","CD14","S100A12","CPPED1","SERPINA1","GCA","CTSZ",
  "PYCARD","TNFRSF1B","S100A8","LILRB2","LAMTOR1","PECAM1","CAP1","NPC2","CTSH","C5AR1",
  "GDI2","TRAPPC1","DBNL","RETN","VAMP8","LTA4H","SELL","HSPA1A","DYNLL1","ALOX5",
  "IQGAP2","FGR","FPR1","CDA","LILRB3","CD74","AP2S1","HLA-DPA1","NCF2","CLEC4A",
  "HLA-DQB1","HLA-DMB","HLA-DMA","PSME1","PSME2","HLA-DQA1","HLA-E","HLA-DRB5","CAPZB",
  "HLA-F","NCF1","ACTB","ARPC3","FCGR3A","ARPC1B","FYB1","BRK1","ICAM3","HCK","FYN",
  "ARPC4","MEF2C","LILRB1","THEMIS2","FCER1A","CD48","POU2F2","CD86","TNFSF13B",
  "TNFSF13","IRF1","AIF1","CSF1R","LRRK2","SULF2","CASP1","PARK7","LILRA5","GPSM3",
  "CYP1B1","CNPY3","CD300LF","ARRB2","IRF2","SAMHD1","TSPO","CX3CR1","LST1","DUSP1",
  "TMEM176B","RIPOR2","TMEM176A","CLU","CDK6","IFI27","HMGA1","RNASE2","TRIM11",
  "CXCR4","TNFAIP3","OASL","IFNGR1","PMAIP1","IFIT3","THBD","MAPKAPK2","CXCL3","CXCL2",
  "ELANE","ADM","CCL3","LITAF","IL1RN","HMOX1","PRNP","CD84","BCL6","RABGEF1","CD83",
  "MYB","ATF4","ARHGEF2","HSPA5","PPP1R15B","PABPC1","P4HB","CHMP4B","MAFB","CCL7","CCL3L1"
)

# Keep genes present in dataset
inflam_genes <- intersect(inflam_genes, rownames(aml))

aml <- AddModuleScore(
  aml,
  features = list(inflam_genes),
  name = "Inflammation"
)

# =====================================================
# 1. Keep only NR and CR samples (remove NA)
# =====================================================
aml_sub <- subset(aml, subset = Response_Group %in% c("NR","CR"))

# =====================================================
# 2. Expression matrix
# =====================================================
expr_all <- GetAssayData(aml_sub, layer = "data")

# =====================================================
# 3. Only selected genes for heatmap
# =====================================================
genes_focus <- c(
  "SAMHD1","NCF2","CKAP4","POU2F2",
  "CTSH","KCTD12","LRRC25",
  "CD14","LILRA6"
)

genes_focus <- intersect(genes_focus, rownames(expr_all))

expr_focus <- expr_all[genes_focus,,drop=FALSE]

# Row-wise Z-score
expr_focus_z <- t(scale(t(expr_focus)))

# =====================================================
# 4. Order samples: NR first → then CR
# =====================================================
meta_ord <- aml_sub@meta.data[
  order(aml_sub@meta.data$Response_Group),
]

# Make sure NR is first
meta_ord$Response_Group <- factor(
  meta_ord$Response_Group,
  levels = c("NR","CR")
)

meta_ord <- meta_ord[order(meta_ord$Response_Group), ]

sample_order <- rownames(meta_ord)
expr_focus_z <- expr_focus_z[, sample_order]

# =====================================================
# 5. Annotation
# =====================================================
annotation_col <- data.frame(
  Response = meta_ord$Response_Group,
  Inflammation = meta_ord$Inflammation1,
  row.names = sample_order
)

ann_colors <- list(
  Response = c(CR="#27AE60", NR="#C0392B"),
  Inflammation = colorRampPalette(c("#ffffcc","maroon"))(100)
)

# =====================================================
# 6. Heatmap
# =====================================================
heat_breaks <- seq(-2, 2, length.out = 101)

pheatmap(
  expr_focus_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,   # Keep NR → CR order
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = heat_breaks,
  fontsize_row = 12,
  border_color = NA,
  main = "Selected Inflammatory Genes (NR → CR)"
)
