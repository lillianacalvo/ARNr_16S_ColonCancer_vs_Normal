# Indicar directorio de trabajo
cd /home/bioinf-2/crc

# Activar ambiente bioinfo
conda activate bioinfo

# Acceso a SRA y descarga de datos del bioproyecto PRJNA284355:
esearch -db sra -query PRJNA284355 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR > srr_crc_ids.txt

# Creación de carpeta para resultados de FastQC en datos crudos
mkdir -p raw_fastqc

# Descarga de todos los archivos .sra, conversión a FASTQ y ejecución de FastQC:
for srr_id in $(cat srr_crc_ids.txt); do
  prefetch $srr_id # Descarga el archivo .sra
  fastq-dump --split-files --gzip $srr_id # Convierte a FASTQ

  # FastQC en los datos crudos (comprimidos en .gz)
  fastqc ${srr_id}_1.fastq.gz -o raw_fastqc # archivo forward
  fastqc ${srr_id}_2.fastq.gz -o raw_fastqc # archivo reverse
done

# MultiQC
mkdir -p /home/bioinf-2/crc/raw_fastqc/multiqc_report
multiqc /home/bioinf-2/crc/raw_fastqc/raw_fastqc -o /home/bioinf-2/crc/raw_fastqc/multiqc_report

# Filtrado y recorte de secuencias
# Directorio de entrada con los archivos fastq.gz
input_dir="/home/bioinf-2/crc/raw_data"

# Archivo de adaptadores
adapters="/home/bioinf-2/adapters.fa"

# Directorio de salida
mkdir /home/bioinf-2/crc/trimmed_data2
output_dir="/home/bioinf-2/crc/trimmed_data2"

# Iterar sobre los archivos forward (_1.fastq.gz)
for forward_file in "$input_dir"/*_1.fastq.gz; do
  # Obtener el nombre base del archivo (sin el sufijo _1.fastq
  base_name=$(basename "$forward_file" "_1.fastq.gz")
  # Construir el nombre del archivo reverse correspondiente (_
  reverse_file="$input_dir/${base_name}_2.fastq.gz"
  # Nombres de los archivos de salida en el directorio corresp
  trimmed_forward="$output_dir/trimmed_${base_name}_1.fastq.gz unpaired_forward="$output_dir/unpaired_${base_name}_1.fastq
  trimmed_reverse="$output_dir/trimmed_${base_name}_2.fastq.gz unpaired_reverse="$output_dir/unpaired_${base_name}_2.fastq
  # Ejecutar Trimmomatic para cada par de archivos
trimmomatic PE -threads 4 \
  "$forward_file" "$reverse_file" \
  "$trimmed_forward" "$unpaired_forward" \
  "$trimmed_reverse" "$unpaired_reverse" \
  ILLUMINACLIP:"$adapters":2:30:5 TRAILING:20 SLIDINGWINDO MINLEN:100
  # Imprimir un mensaje de estado
  echo "Trimming completado para: $base_name"
done

# Crear archivos de metadatos y manifiesto
# Manifiesto:
# Ruta a la carpeta con los datos
data_dir="/home/bioinf-2/crc/trimmed_data2"
manifest_file="manifiesto.tsv"
# Crear el archivo de manifiesto y escribir la cabecera
echo "#SampleID,Forward,Reverse" > $manifest_file
# Recorrer los archivos en la carpeta
for forward_file in "$data_dir"/*_1.fastq.gz; do
  if [[ -e "$forward_file" ]]; then
    # Extraer el SampleID
    sample_id=$(basename "$forward_file" _1.fastq.gz)
    reverse_file="$data_dir/${sample_id}_2.fastq.gz"
    # Verificar si el archivo reverse existe
    if [[ -e "$reverse_file" ]]; then
      # Escribir la línea en el manifiesto
      echo "$sample_id,$forward_file,$reverse_file" >> $ma
    else
      echo "Advertencia: No se encontró el archivo
      reverse para $sample_id"
    fi
  fi
done

#Se indica el directorio de trabajo
cd /home/bioinf-2/crc/qiime2/

#Se activa el ambiente QUIIME2
conda activate qiime2-2021.11

#Indicar el archivo de metadatos
qiime metadata tabulate --m-input-file metadatos.tsv --o-visualization tabulated-metadata.qzv
#Indicar a qiime el manifiesto
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifiesto.tsv --output-path paired-end-sequences --input-format PairedEndFastqManifestPhred33V2

#Análisis de calidad de las secuencias - denoising
qiime demux summarize --i-data paired-end-sequences.qza --o-visualization demux.qzv
qiime dada2 denoise-paired --i-demultiplexed-seqs paired-end-seq --p-trunc-len-f 250 --p-trim-left-f 0 --p-trunc-len-r 250 --p-trim-left-r 0 --o-representative-sequences rep-seqs.qza --o-table table.qza --o-denoising-stats stats.qza

#Se evalúa la calidad de los datos depurados
qiime metadata tabulate --m-input-file stats.qza --o-visualization stats-dada2.qzv

# Remoción de quimeras
# Identificación de quimeras de novo
qiime vsearch uchime-denovo --i-table table.qza --i-sequences rep-seqs.qza --output-dir uchime-dn-out
qiime metadata tabulate --m-input-file uchime-dn-out/stats.qza --o-visualization uchime-dn-out/stats.qzv

#Filtrar cuadro de abundancia "otu-table" (remoción de quimeras)
qiime feature-table filter-features --i-table table.qza --m-metadata-file uchime-dn-out/nonchimeras.qza --o-filtered-table uchime-dn-out/table-nonchimeric-wo-borderline.qza
qiime feature-table filter-seqs --i-data rep-seqs.qza --m-metadata-file uchime-dn-out/nonchimeras.qza --o-filtered-data uchime-dn-out/rep-seqs-nonchimeric-wo-borderline.qza
qiime feature-table summarize --i-table uchime-dn-out/table-nonchimeric-wo-borderline.qza --o-visualization uchime-dn-out/table-nonchimeric-wo-borderline.qzv


#Análisis filogenético
#El siguiente comando genera un árbol filogenético de los ASVs necesario para el análisis de rarefacción:
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences uchime-dn-out/rep-seqs-nonchimeric-wo-borderline.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza
# Visualización de los datos
qiime feature-table summarize \
--i-table uchime-dn-out/table-nonchimeric-wo-borderline.qza \
--o-visualization uchime-dn-out/table-nonchimeric-wo-borderline.qzv \
--m-sample-metadata-file metadatos.tsv
#El siguiente comando genera un archivo con la lista de ASVs con su secuencia y un link para buscar en BLAST:
qiime feature-table tabulate-seqs \
--i-data uchime-dn-out/rep-seqs-nonchimeric-wo-borderline.qza \
--o-visualization uchime-dn-out/rep-seqs.qzv

# Rarefacción
qiime diversity alpha-rarefaction \
 --i-table uchime-dn-out/table-nonchimeric-wo-borderline.qza \
 --i-phylogeny rooted-tree.qza \
 --p-max-depth 36694 \
 --m-metadata-file metadatos.tsv \
 --o-visualization alpha-rarefaction.qzv
 
qiime feature-table rarefy \
 --i-table uchime-dn-out/table-nonchimeric-wo-borderline.qza \
 --p-sampling-depth 5000 \
 --o-rarefied-table rarefied-table.qza

#El siguiente comando genera un archivo con estadísticas resumen de ASVs recortada:
qiime feature-table summarize --i-table rarefied-table.qza
--o-visualization rarefied-table.qzv --m-sample-metadata-fil

#Clasificación taxonómica
#Se descarga la base GSR: no tiene un clasificador pre entrenado disponible en qiime2.
wget https://manichanh.vhir.org/gsrdb/GSR-DB_full-16S.tar.gz
tar -xvzf GSR-DB_full-16S.tar.gz

# Entrenamiento naive-bayes del clasificador
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads GSR-DB_full-16S_filt_seqs.qza \
  --i-reference-taxonomy GSR-DB_full-16S_filt_taxa.qza \
  --o-classifier GSR-DB_full-16S_filt_classifier.qza
  
#Asignación de la taxonomía
qiime feature-classifier classify-sklearn \
  --i-classifier GSR-DB_full-16S_filt_classifier.qza \
  --i-reads uchime-dn-out/rep-seqs-nonchimeric-wo-borderline.qza \
  --o-classification taxonomy_nochimeras.qza
  
qiime metadata tabulate \
  --m-input-file taxonomy_nochimeras.qza \
  --o-visualization taxonomy_nochimeras.qzv
  
qiime taxa barplot
  --i-table rarefied-table.qza \
  --i-taxonomy taxonomy_nochimeras.qza \
  --m-metadata-file metadatos.tsv \
  --o-visualization taxa-bar-plots.qzv

qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted-tree.qza \
--i-table rarefied-table.qza \
--p-sampling-depth 5000 \
--m-metadata-file metadatos.tsv \
--output-dir core-metrics-results-3

#Diversidad alfa (variables categóricas)
qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results-3/faith_pd_vector.qza \
--m-metadata-file metadatos.tsv \
--o-visualization core-metrics-results-3/faith-pd-group-significance.qzv

#Diversidad Beta (variables categóricas)
qiime diversity beta-group-significance \
--i-distance-matrix core-metrics-results-3/unweighted_unifrac_distance_matrix.qza \
--m-metadata-file metadatos.tsv \
--m-metadata-column Description \
--o-visualization core-metrics-results-3/unweighted-unifrac-description-significance.qzv \
--p-pairwise

#Análisis de PCA unweighted
qiime emperor plot \
--i-pcoa core-metrics-results-3/unweighted_unifrac_pcoa_results.qza \
--m-metadata-file metadatos.tsv \
--o-visualization core-metrics-results-3/unweighted-unifrac-emperor.qzv

#Análisis de PCA weighted
qiime diversity pcoa \
--i-distance-matrix weighted_unifrac_distance_matrix.qza \
--o-pcoa weighted-unifrac-pcoa-results.qza

qiime emperor plot \
--i-pcoa weighted-unifrac-pcoa-results.qza \
--m-metadata-file metadatos.tsv \
--o-visualization weighted-unifrac-emperor-Description1.qzv
