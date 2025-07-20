# Dokumentacja pipeline Snakemake

## Plik konfiguracyjny: `config.yaml`
- Zawiera listę próbek (`samples`), numery odczytów (`reads`), ścieżki do danych surowych (`raw_dir`), wyników (`results_dir`) oraz ścieżkę do referencyjnego genomu (`ref_genome`).

---

## Zmienne globalne

- `SAMPLES` — lista próbek z `config.yaml`
- `READS` — lista numerów odczytów (np. `[1, 2]`)
- `RAW_DIR` — katalog z surowymi danymi FASTQ
- `RESULTS_DIR` — katalog na wyniki pipeline
- Podkatalogi wynikowe:  
  - `TRIMMED_DIR` — z przyciętymi odczytami  
  - `FASTQC_DIR` — wyniki kontroli jakości FastQC  
  - `MULTIQC_DIR` — wyniki MultiQC (zbiorcze raporty)  
  - `STAR_INDEX` — indeks genomu do alignera STAR  
  - `ALIGNMENTS_DIR` — wyniki mapowania (plik BAM)  
  - `COVERAGE_DIR` — profile pokrycia genomu (bedgraph)  
  - `VARIANTS_DIR` — wywołane warianty (plik VCF)  
- `REF_GENOME` — plik referencyjnego genomu FASTA

---

## Reguły (rules)

### all
**Opis:**  
Reguła zbiorcza, która wymusza wyprodukowanie finalnych plików dla wszystkich próbek i odczytów.  
**Wyjście:**  
- Raporty FastQC (zip) dla każdego odczytu  
- Przycięte pliki FASTQ (pary)  
- Posortowane pliki BAM i indeksy BAM  
- Profile pokrycia w formacie bedgraph  
- Pliki VCF z wariantami  
- Raport MultiQC

---

### prepare_dirs
**Opis:**  
Tworzy wszystkie katalogi wynikowe wymagane w pipeline, jeśli jeszcze nie istnieją.  
**Wyjście:**  
- Katalog `results_dir` i wszystkie podkatalogi

---

### fastqc
**Opis:**  
Analiza jakości surowych odczytów FASTQ za pomocą FastQC.  
**Wejście:**  
- Surowy plik FASTQ (`data/raw/{sample}_{read}.fastq`)  
**Wyjście:**  
- Raport HTML i zip FastQC  
**Kontener:**  
- `biocontainers/fastqc:v0.11.9_cv8`  
**Wątki:** 2

---

### multiqc
**Opis:**  
Konsoliduje raporty FastQC w jeden zbiorczy raport HTML.  
**Wejście:**  
- Wszystkie pliki zip FastQC  
**Wyjście:**  
- `multiqc_report.html`  
**Kontener:**  
- `quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0`

---

### trimmomatic
**Opis:**  
Przycinanie surowych odczytów (usuwanie adapterów i niskiej jakości fragmentów) z par FASTQ.  
**Wejście:**  
- Parowane surowe pliki FASTQ (`_1.fastq` i `_2.fastq`)  
**Wyjście:**  
- Przycięte pliki FASTQ (paired)  
- Tymczasowe pliki FASTQ (unpaired)  
**Kontener:**  
- `quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2`  
**Wątki:** 4  
**Parametry:**  
- Adaptery TruSeq3-PE  
- Min. długość 36, ustawienia jakościowe (LEADING, TRAILING, SLIDINGWINDOW)

---

### star_index
**Opis:**  
Budowa indeksu referencyjnego genomu dla alignera STAR.  
**Wejście:**  
- Plik FASTA referencyjnego genomu  
**Wyjście:**  
- Plik `SA` w katalogu indeksu STAR (dotykany jako marker)  
**Kontener:**  
- `quay.io/biocontainers/star:2.7.10a--h9ee0642_1`  
**Wątki:** 4

---

### star_align
**Opis:**  
Mapowanie przyciętych odczytów na genom referencyjny przy pomocy STAR.  
**Wejście:**  
- Przycięte pliki FASTQ (paired)  
- Indeks STAR  
**Wyjście:**  
- Posortowany plik BAM  
**Kontener:**  
- `quay.io/biocontainers/star:2.7.10a--h9ee0642_1`  
**Wątki:** 8

---

### index_bam
**Opis:**  
Indeksowanie pliku BAM (tworzenie pliku .bai) za pomocą samtools.  
**Wejście:**  
- BAM z mapowania STAR  
**Wyjście:**  
- Plik indeksu BAM `.bai`  
**Kontener:**  
- `quay.io/biocontainers/samtools:1.15.1--h1170115_0`

---

### coverage
**Opis:**  
Obliczanie pokrycia genomu na podstawie BAM za pomocą bedtools.  
**Wejście:**  
- BAM i jego indeks `.bai`  
**Wyjście:**  
- Plik bedgraph z profilem pokrycia  
**Kontener:**  
- `quay.io/biocontainers/bedtools:2.30.0--h468198e_3`

---

### variant_calling
**Opis:**  
Wykrywanie wariantów (SNP, indel) w oparciu o pliki BAM i genom referencyjny za pomocą freebayes.  
**Wejście:**  
- BAM i indeks BAM  
- Genom referencyjny  
**Wyjście:**  
- Plik VCF z wariantami  
**Kontener:**  
- `quay.io/biocontainers/freebayes:1.3.5--h2d02072_3`

---

## Podsumowanie

Pipeline realizuje kompleksową analizę danych sekwencjonowania:

- Kontrola jakości (FastQC + MultiQC)
- Przycinanie surowych odczytów (Trimmomatic)
- Indeksowanie genomu (STAR)
- Mapowanie (STAR)
- Indeksowanie BAM (samtools)
- Analiza pokrycia genomu (bedtools)
- Wywołanie wariantów (freebayes)
