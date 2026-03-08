#### Analisis Diferensial Ekspresi Gen pada Subtipe Kanker Payudara Menggunakan Dataset GSE45827



###### Pendahuluan 

Kanker payudara dapat dikelompokkan menjadi beberapa subtipe, seperti Basal (Triple Negative), HER2, Luminal A, dan Luminal B. Setiap subtipe tentunya memiliki karakter biologis yang berbeda-beda. Hal ini mendorong penggunaan teknologi microarray yang bisa menganalisis ekspresi ribuan gen sehingga gen-gen yang memiliki perubahan ekspresi signifikan pada penyakit kanker dapat dikenali.



Penelitian menggunakan datasetGSE45827 yang berisi data ekspresi gen dari keempat subtipe kanker payudara yang telah disebutkan, jaringan kultur sel, serta jaringan normal manusia. Data diperoleh melalui platform Affymetrix Human Genome U133 Plus 2.0 Array. Jenis eksperimen yang dilakuka untuk memperoleh data tersebut ialah *expression profilling by microarray* dengan sampel Basal (Triple Negative) sebanyak 41, HER2 sebanyak 30, Luminal A sebanyak 29, Luminal B sebanyak 30, jaringan normal sebanyak 11, dan cell line sebanyak 14 sampel.



###### Tujuan Penelitian

Mengidentifikasi Differentially Expressed Genes (DEG) antara jaringan kanker payudara dan jaringan normal.



###### Alur Analisis

Berikut ialah versi ringkas alur penelitian. Langkah-langkah selengkapnya terlampir.



Download GEO (GSE45827)

&nbsp;       ↓

Preprocessing 

&nbsp;       ↓

DEG analysis (limma)

&nbsp;       ↓

Volcano plot

&nbsp;       ↓

Heatmap

&nbsp;       ↓

GO enrichment

&nbsp;       ↓

KEGG enrichment

&nbsp;       ↓

Pathway enrichment map



Pada tahap pre-processing data, dilakukan:

1. Transformasi log dan normalisasi
2. Pemeriksaan distribusi gen menggunakan boxplot dan density plot
3. Klaster gen menggunakan UMAP



Adapun analisis diferensial ekspresi gen yang dilakukan melalui limma dalam bahasa R memiliki kriteria berupa |log2 Fold Change| > 1 dan adjusted p.value < 0.01



###### Hasil dan Interpretasi

1. Boxplot:
   - Distribusi nilai ekspresi antar sampel relatif seragam.
   - Normalisasi dan transformasi data dilakukan dengan baik.
   
2. Density Plot:
   - Distribusi nilai ekspresi sampel sebagian besar mirip.
   - Distribusi pada sel kanker relatif serupa, jika dibandingkan dengan jaringan normal dan cell line.
   
3. UMAP: 
   Sampel cell line dan jaringan normal berada relatif jauh dari kelompok tumor, sedangkan  sampel tumor cenderung membentuk klaster tersendiri. Subtipe tumor seperti Basal, HER2, Luminal A, dan Luminal B menunjukkan kecenderungan pengelompokan sesuai tipe masing-masing, meskipun terdapat beberapa overlap antar subtipe.
   
4. Volcano plot:
   Terdapat sejumlah gen yang mengalami perubahan ekspresi signifikan antara jaringan tumor dan jaringan normal. Artinya, gen-gen tersebut kemungkinan berperan dalamproses perkembangan kanker payudara.
   
5. Heatmap: 
   Sampel jaringan normal membentuk klaster yang berbeda dari sampel tumor. Cell line juga menunjukkan pola ekspresi gen yang berbeda dari jaringan tumor primer. Sementara itu, subtipe tumor menunjukkan pola ekspresi yang relatif mirip satu sama lain, meskipun masih terdapat beberapa variasi antar subtipe. 
   
6. GO Enrichment Analysis:
   Pada gen upregulated, beberapa fungsi biologis yang dominan antara lain: actin filament organization, cytoskeleton organization, dan cell migration. Sementara itu, pada  gen downregulated, beberapa fungsi biologis yang muncul antara lain: regulation of actin cytoskeleton dan cell adhesion processes.
   
7. KEGG Pathway Enrichment:
8. Hasil KEGG menunjukkan beberapa jalur biologis yang signifikan, di antaranya yaitu: Endocytosis pathway (upregulated) dan Human papillomavirus infection pathway (downregulated).



###### Kesimpulan

1. Data ekspresi gen telah ter-normalisasi dengan baik dan tidak menunjukkan adanya batch effect yang signifikan.
2. Analisis diferensial ekspresi berhasil mengidentifikasi sejumlah gen yang mengalami perubahan ekspresi signifikan antara jaringan tumor dan jaringan normal.
3. Visualisasi UMAP dan heatmap menunjukkan bahwa jaringan tumor memiliki pola ekspresi gen yang berbeda dari jaringan normal dan cell line.
4. Analisis GO enrichment menunjukkan bahwa gen yang mengalami perubahan ekspresi banyak terlibat dalam sitoskeleton dan sel.
5. Analisis KEGG pathway mengidentifikasi jalur biologis seperti endositosis yang berpotensi berperan dalam perkembangan kanker payudara.









###### Penelitian terkait:

&nbsp;	Gruosso, T., Mieulet, V., Cardon, M., Bourachot, B., Kieffer, Y., Devun, F., Dubois, T., Dutreix, M., Vincent-Salomon, A., Miller, K. M., \& Mechta-Grigoriou, F. (2016). Chronic oxidative stress promotes H2AX protein degradation and enhances chemosensitivity in breast cancer patients. EMBO molecular medicine, 8(5), 527–549. https://doi.org/10.15252/emmm.201505891



