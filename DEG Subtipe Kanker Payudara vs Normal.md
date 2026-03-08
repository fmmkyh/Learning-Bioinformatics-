\#Modul: Analisis Ekspresi Gen Subtipe Payudara vs Normal 



\#Dataset: GSE45827 (Kanker Basal vs Normal)

\#Platform: Microarray (Affymetrix Human Genome U133 Plus 2.0 Array - GPL570)

\#Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG) antara jaringan kanker payudara dan jaringan normal.



\#PART A. PENGANTAR KONSEP

\#Analisis ekspresi gen bertujuan untuk membandingkan ekspresi gen

\#antara jaringan payudara normal dengan jaringan payudara dengan biopsi kanker (basal)



\#Pada modul ini kita menggunakan pendekatan statistik limma (Linear Models

\#for Microarray Data), yang merupakan standar emas untuk datamicroarray.



\#PART B. PERSIAPAN LINGKUNGAN KERJA (INSTALL \& LOAD PACKAGE)



\#Apa itu package?

\#Package adalah kumpulan fungsi siap pakai di R

\#Bioinformatika di R sangat bergantung pada package dari CRAN dan Bioconductor



\#1. Install BiocManager (manajer paket Bioconductor)

\#IF adalah struktur logika : “jika kondisi terpenuhi, lakukan aksi”

if (!require("BiocManager", quietly = TRUE)) {

&nbsp;	nstall.packages("BiocManager")

}



\# 2. Install paket Bioconductor (GEOquery \& limma)

\#GEOquery: mengambil data dari database GEO

\#limma: analisis statistik ekspresi gen

BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update =FALSE)



\#hgu133plus2.db: database anotasi khusus untuk GSE42568)

BiocManager::install("hgu133plus2.db", ask = FALSE, update =

FALSE)



\#3. Install paket CRAN untuk visualisasi dan manipulasi data

\#phetmap: heatmap ekspresi gen

\#ggplot2: grafik (volcano plot)

\#dplyr: manipulasi tabel data

install.packages(c("pheatmap", "ggplot2", "dplyr"))



\#umap: grafik (plot UMAP)

if (!requireNamespace("umap", quietly = TRUE)) {

install.packages("umap")

}



\#clusterProfiler : analisis enrichment

\#enrichplot → membuat enrichment plot

\#org.Hs.eg.db → anotasi gen manusia

BiocManager::install(c("clusterProfiler", "enrichplot", "org.Hs.eg.db"))



\#4. Memanggil library

\#library() digunakan agar fungsi di dalam package bisa digunakan

library(GEOquery)

library(limma)

library(pheatmap)

library(ggplot2)

library(dplyr)

library(hgu133plus2.db)

library(AnnotationDbi)

library(umap)

library(clusterProfiler)

library(enrichplot)

library(org.Hs.eg.db)



\#PART C. PENGAMBILAN DATA DARI GEO

\#GEO (Gene Expression Omnibus) adalah database publik milik NCBI

\#getGEO(): fungsi untuk mengunduh dataset berdasarkan ID GEO

\#GSEMatrix = TRUE -> data diambil dalam format ExpressionSet

\#AnnotGPL = TRUE -> anotasi gen (Gene Symbol) ikut diunduh



gset <- getGEO("GSE45827", GSEMatrix = TRUE, AnnotGPL = TRUE)\[\[1]]



\#ExpressionSet berisi:

\# - exprs() : matriks ekspresi gen

\# - pData() : metadata sampel

\# - fData() : metadata fitur (probe / gen)



\#PART D. PRE-PROCESSING DATA EKSPRESI

\# exprs(): mengambil matriks ekspresi gen

\# Baris = probe/gen

\# Kolom = sampel

ex <- exprs(gset)



\#Mengapa perlu log2 transformasi?

\#Data microarray mentah memiliki rentang nilai sangat besar.

\#Log2 digunakan untuk:

\#1. Menstabilkan varians

\#2. Mendekati asumsi model linear

\#3. Memudahkan interpretasi log fold change



\#quantile(): menghitung nilai kuantil (persentil)

\#as.numeric(): mengubah hasil quantile (yang berupa named vector)

\#menjadi vektor numerik biasa agar mudah dibandingkan

qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))



\#LogTransform adalah variabel logika (TRUE / FALSE)

\#Operator logika:

\#> : lebih besar dari

\#|| : OR (atau)

\#\&\& : AND (dan)

LogTransform <- (qx\[5] > 100) || (qx\[6] - qx\[1] > 50 \&\& qx\[2] > 0)



\#IF statement:

\#Jika LogTransform = TRUE, maka lakukan log2

if (LogTransform) {

\# Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA

&nbsp;	x\[ex <= 0] <- NA

&nbsp;	x <- log2(ex)

} 



\#Namun LogTransform pada data yang diambil bernilai FALSE, 

\#maka fungsi ini dilewati karena telah dilakukan log2 pada data



\#PART E. DEFINISI KELOMPOK SAMPEL



\#pData(): metadata sampel

\#source\_name\_ch1 berisi informasi kondisi biologis sampel

group\_info <- pData(gset)\[\["source\_name\_ch1"]]



\#make.names(): mengubah teks menjadi format valid untuk R

groups <- make.names(group\_info)



\#factor():

\#Mengubah data kategorik menjadi faktor

\#Faktor sangat penting untuk analisis statistik di R

gset$group <- factor(groups)



\#levels(): melihat kategori unik dalam faktor

nama\_grup <- levels(gset$group)

print(nama\_grup)



\#PART F. DESIGN MATRIX (KERANGKA STATISTIK)



\#model.matrix():

\#Membuat matriks desain untuk model linear

\#~0 berarti TANPA intercept (best practice limma)

design <- model.matrix(~0 + gset$group)



\#colnames(): memberi nama kolom agar mudah dibaca

colnames(design) <- levels(gset$group)



\#Menentukan perbandingan biologis

tumor\_basal <- nama\_grup\[1]

tumor\_kultur <- nama\_grup\[2] 

tumor\_HER2 <- nama\_grup\[3]

tumor\_luminalA <- nama\_grup\[4]

tumor\_luminalB <-nama\_grup\[5]

sel\_normal <- nama\_grup\[6]



sel\_tumor <- c(

&nbsp;	tumor\_basal,

&nbsp;	tumor\_kultur,

&nbsp;	tumor\_HER2,

&nbsp;	tumor\_luminalA,

&nbsp;	tumor\_luminalB

)



contrast\_formula <- paste(sel\_tumor, "-", sel\_normal)



print(paste("Kontras yang dianalisis:", contrast\_formula))



\#PART G. ANALISIS DIFFERENTIAL EXPRESSION (LIMMA)



\#lmFit():

\#Membangun model linear untuk setiap gen

fit <- lmFit(ex, design)



\#makeContrasts(): mendefinisikan perbandingan antar grup

contrast\_matrix <- makeContrasts(contrasts = contrast\_formula, levels = design)



\#contrasts.fit(): menerapkan kontras ke model

fit2 <- contrasts.fit(fit, contrast\_matrix)



\#eBayes():

\#Empirical Bayes untuk menstabilkan estimasi varians

fit2 <- eBayes(fit2)



\#topTable():

\#Mengambil hasil akhir DEG

\#adjust = "fdr" -> koreksi multiple testing

\#p.value = 0.01 -> gen sangat signifikan

topTableResults <- topTable(

&nbsp;	fit2,

&nbsp;	adjust = "fdr",

&nbsp;	sort.by = "B",

&nbsp;	number = Inf,

&nbsp;	p.value = 0.01

)



head(topTableResults)



\#PART H. ANOTASI NAMA GEN



\#Penting:

\#Pada data microarray Affymetrix, unit analisis awal adalah PROBE,

\#bukan gen. Oleh karena itu, anotasi ulang diperlukan menggunakan

\#database resmi Bioconductor.



\#Mengambil ID probe dari hasil DEG

probe\_ids <- rownames(topTableResults)



\#Mapping probe -> gene symbol \& gene name

gene\_annotation <- AnnotationDbi::select(

&nbsp;	hgu133plus2.db,

&nbsp;	keys = probe\_ids,

&nbsp;	columns = c("SYMBOL", "GENENAME"),

&nbsp;	keytype = "PROBEID"

)



\#Gabungkan dengan hasil limma

topTableResults$PROBEID <- rownames(topTableResults)



topTableResults <- merge(

&nbsp;	topTableResults,

&nbsp;	gene\_annotation,

&nbsp;	by = "PROBEID",

&nbsp;	all.x = TRUE

)



\#Cek hasil anotasi

head(topTableResults\[, c("PROBEID", "SYMBOL", "GENENAME")])



\#PART I.1 BOXPLOT DISTRIBUSI NILAI EKSPRESI



\#Boxplot digunakan untuk:

\#- Mengecek distribusi nilai ekspresi antar sampel

\#- Melihat apakah ada batch effect

\#- Mengevaluasi apakah normalisasi/log-transform sudah wajar



\#Set warna berdasarkan grup

group\_colors <- as.numeric(gset$group)



boxplot(

&nbsp;	ex,

&nbsp;	col = group\_colors,

&nbsp;	las = 2,

&nbsp;	outline = FALSE,

&nbsp;	main = "Boxplot Distribusi Nilai Ekspresi per Sampel",

&nbsp;	ylab = "Expression Value (log2)"

)



legend(

&nbsp;	"topright",

&nbsp;	legend = levels(gset$group),

&nbsp;	fill = unique(group\_colors),

&nbsp;	cex = 0.5

)



\#PART I.2 DISTRIBUSI NILAI EKSPRESI (DENSITY PLOT)



\#Density plot menunjukkan sebaran global nilai ekspresi gen

\#Digunakan untuk:

\#- Mengecek efek log-transform

\#- Membandingkan distribusi antar grup



\#Gabungkan ekspresi \& grup ke data frame

expr\_long <- data.frame(

&nbsp;	Expression = as.vector(ex),

&nbsp;	Group = rep(gset$group, each = nrow(ex))

)



ggplot(expr\_long, aes(x = Expression, color = Group)) +

&nbsp;	geom\_density(linewidth = 1) +

&nbsp;	theme\_minimal() +

&nbsp;	labs(

&nbsp;		title = "Distribusi Nilai Ekspresi Gen",

&nbsp;		x = "Expression Value (log2)",

&nbsp;		y = "Density"

)



\#PART I.3 UMAP (VISUALISASI DIMENSI RENDAH)



\#UMAP digunakan untuk:

\#- Mereduksi ribuan gen menjadi 2 dimensi

\#- Melihat pemisahan sampel secara global

\#- Alternatif PCA (lebih sensitif ke struktur lokal)



\#Transpose matriks ekspresi:

\#UMAP bekerja pada OBSERVATION = sampel

umap\_input <- t(ex)



\#Jalankan UMAP

umap\_result <- umap(umap\_input)



\#Simpan hasil ke data frame

umap\_df <- data.frame(

&nbsp;	UMAP1 = umap\_result$layout\[, 1],

&nbsp;	UMAP2 = umap\_result$layout\[, 2],

&nbsp;	Group = gset$group

)



\#Plot UMAP

ggplot(umap\_df, aes(x = UMAP1, y = UMAP2, color = Group)) +

&nbsp;	geom\_point(size = 3, alpha = 0.8) +

&nbsp;	theme\_minimal() +

&nbsp;	labs(

&nbsp;		title = "UMAP Plot Sampel Berdasarkan Ekspresi Gen",

&nbsp;		x = "UMAP 1",

&nbsp;		y = "UMAP 2"

)



\#PART J.1 VISUALISASI VOLCANO PLOT



\#Volcano plot menggabungkan:

\#- Log fold change (efek biologis)

\#- Signifikansi statistic



contrast\_cols <- colnames(topTableResults)\[2:6]



for (contrast in contrast\_cols) {

&nbsp;	volcano\_data <- data.frame(

&nbsp;		logFC = topTableResults\[\[contrast]],

&nbsp;		adj.P.Val = topTableResults$adj.P.Val,

&nbsp;		Gene = topTableResults$SYMBOL

&nbsp;	)



&nbsp;	#Klasifikasi status gen

&nbsp;	volcano\_data$status <- "NO"

&nbsp;	volcano\_data$status\[volcano\_data$logFC > 1 \& volcano\_data$adj.P.Val < 0.01] <- "UP"

&nbsp;	volcano\_data$status\[volcano\_data$logFC < -1 \& volcano\_data$adj.P.Val < 0.01] <- "DOWN"



&nbsp; 	p <- ggplot(volcano\_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) + geom\_point(alpha = 0.6) + scale\_color\_manual(values = c("DOWN"="blue","NO"="grey","UP"="red")) + geom\_vline(xintercept=c(-1,1), linetype="dashed") + geom\_hline(yintercept=-log10(0.01), linetype="dashed") + theme\_minimal() + ggtitle(paste("Volcano Plot:", contrast))

&nbsp; 	print(p)

}



\#PART J.2 VISUALISASI HEATMAP



\#Heatmap digunakan untuk melihat pola ekspresi gen

\#antar sampel berdasarkan gen-gen paling signifikan



\#Pilih 50 gen paling signifikan berdasarkan adj.P.Val

topTableResults <- topTableResults\[

&nbsp;	order(topTableResults$adj.P.Val),

]



top50 <- head(topTableResults, 50)



\#Ambil matriks ekspresi untuk gen terpilih

mat\_heatmap <- ex\[top50$PROBEID, ]



\#Gunakan Gene Symbol (fallback ke Probe ID)

gene\_label <- ifelse(

&nbsp;	is.na(top50$SYMBOL) | top50$SYMBOL == "",

&nbsp;	top50$PROBEID, # jika SYMBOL kosong → probe ID

&nbsp;	top50$SYMBOL # jika ada → gene symbol

)



rownames(mat\_heatmap) <- gene\_label



\#Pembersihan data (WAJIB agar tidak error hclust)

\#Hapus baris dengan NA

mat\_heatmap <- mat\_heatmap\[

&nbsp;	rowSums(is.na(mat\_heatmap)) == 0,

]



\#Hapus gen dengan varians nol

gene\_variance <- apply(mat\_heatmap, 1, var)

mat\_heatmap <- mat\_heatmap\[gene\_variance > 0, ]



\#Anotasi kolom (kelompok sampel)

annotation\_col <- data.frame(

&nbsp;	group = gset$group

)



rownames(annotation\_col) <- colnames(mat\_heatmap)



\#Visualisasi heatmap

pheatmap( 

&nbsp;	mat\_heatmap,

&nbsp;	scale = "row", # Z-score per gen

&nbsp;	annotation\_col = annotation\_col,

&nbsp;	show\_colnames = FALSE, # nama sampel dimatikan

&nbsp;	show\_rownames = TRUE,

&nbsp;	fontsize\_row = 7,

&nbsp;	clustering\_distance\_rows = "euclidean",

&nbsp;	clustering\_distance\_cols = "euclidean",

&nbsp;	clustering\_method = "complete",

&nbsp;	main = "Top 50 Differentially Expressed Genes"

)



\#PART K. MENYIMPAN HASIL DEG



\# write.csv(): menyimpan hasil analisis ke file CSV

write.csv(topTableResults, "Hasil\_GSE45827\_DEG.csv")



message("Analisis selesai. File hasil telah disimpan.")



\#PART L. ENRICHMENT GO



\#Memilih gen signifikan dari DEG

for (contrast in contrast\_cols) {

&nbsp;	deg\_genes <- topTableResults %>% filter(adj.P.Val < 0.01 \& abs(.data\[\[contrast]]) > 1)

}



\#Mengambil simbol

gene\_symbols <- unique(deg\_genes$SYMBOL)

gene\_symbols <- gene\_symbols\[!is.na(gene\_symbols)]



\#GO

ego <- enrichGO( 

&nbsp;	gene = gene\_symbols,

&nbsp;	OrgDb = org.Hs.eg.db,

&nbsp;     	keyType = "SYMBOL", 

      	ont = "BP",

&nbsp;  	pAdjustMethod = "BH",

&nbsp;   	pvalueCutoff = 0.01,

&nbsp;	qvalueCutoff = 0.01,

&nbsp;   	readable = TRUE

)



\#Visualisasi GO

dotplot(ego, showCategory = 15) +

&nbsp;   ggtitle("GO Biological Process Enrichment")



barplot(ego, showCategory = 15) +

&nbsp;   ggtitle("Top Enriched GO Terms")



cnetplot(ego, showCategory = 5)



\#Simpan hasil GO Enrichment

enrichment\_results <- as.data.frame(ego)



write.csv(

&nbsp;   enrichment\_results,

&nbsp;   "GO\_Enrichment\_GSE45827.csv",

&nbsp;   row.names = FALSE

)



\#PART M. KEGG PATHWAY ENRICHMENT

\#Konversi Gen Symbol -> Entrez ID

gene\_df <- bitr(

&nbsp;   gene\_symbols,

&nbsp;   fromType = "SYMBOL",

&nbsp;   toType = "ENTREZID",

&nbsp;   OrgDb = org.Hs.eg.db

)



entrez\_genes <- gene\_df$ENTREZID



\#Analisis

kegg\_enrich <- enrichKEGG(

&nbsp;   gene = entrez\_genes,

&nbsp;   organism = "hsa", #kode manusia

&nbsp;   pvalueCutoff = 0.01

)



\#Visualisasi KEGG

dotplot(kegg\_enrich, showCategory = 15) +

&nbsp;   ggtitle("KEGG Pathway Enrichment")



barplot(kegg\_enrich, showCategory = 15) +

&nbsp;   ggtitle("Top KEGG Pathways")



cnetplot(kegg\_enrich, showCategory = 5)



\#Simpan hasil KEGG pathway enrichment

kegg\_results <- as.data.frame(kegg\_enrich)



write.csv(

&nbsp;   kegg\_results,

&nbsp;   "KEGG\_Enrichment\_GSE45827.csv",

&nbsp;   row.names = FALSE

)



