## Checking FASTQ

For batch processing 
##### Step 1: Open a new script
```bash
nano count_reads.sh
```
##### Step 2: Paste the following code
This script counts reads in paired-end FASTQ files and saves results to a CSV.

```bash
#!/bin/bash
set -euo pipefail

INDIR="raw_data"
OUTDIR="csv_output"
mkdir -p "$OUTDIR"

OUTFILE="$OUTDIR/fastq_read_counts.csv"
echo "Sample,R1_reads,R2_reads" > "$OUTFILE"
echo "📊 Counting reads in FASTQ files from '$INDIR'..."

for R1 in "$INDIR"/*_1.fastq.gz "$INDIR"/*_R1.fastq.gz; do
    [[ -f "$R1" ]] || continue
    SAMPLE=$(basename "$R1" | sed -E 's/_R?1.*\.fastq\.gz//')
    R2=""
    for suffix in "_2.fastq.gz" "_R2.fastq.gz" "_R2_*.fastq.gz"; do
        [[ -f "$INDIR/${SAMPLE}${suffix}" ]] && R2="$INDIR/${SAMPLE}${suffix}" && break
    done
    R1_COUNT=$(( $(zcat "$R1" | wc -l) / 4 ))
    R2_COUNT=$([[ -n "$R2" ]] && echo $(( $(zcat "$R2" | wc -l) / 4 )) || echo "NA")
    echo "$SAMPLE,$R1_COUNT,$R2_COUNT" >> "$OUTFILE"
    echo "✅ $SAMPLE → R1: $R1_COUNT | R2: $R2_COUNT"
done

echo "🎉 All done! Read counts saved to '$OUTFILE'"
```

##### Step 3: Save and exit nano
Press Ctrl + O → Enter (to write the file)
Press Ctrl + X → Exit nano

##### Step 3: Make the script executable
```bash
chmod +x count_reads.sh
```
##### Step 4: Run the script
```bash
./count_reads.sh
```

### 3. Base composition

```bash
#!/bin/bash
for fq in raw_data/SRR28821350_1.fastq.gz raw_data/SRR28821350_2.fastq.gz; do
    [ -f "$fq" ] || continue
    echo "Counting bases in $fq..."
    zcat "$fq" | awk 'NR%4==2 { for(i=1;i<=length($0);i++) b[substr($0,i,1)]++ } 
    END { for(base in b) print base, b[base] }'
    echo "----------------------"
done
```

### 5.   Checking FASTQ Pairing 

We ensured all our FASTQ files are correctly paired before running any bioinformatics analysis.


##### Step 1: Create the script
```bash
nano check_fastq_pairs.sh
```
##### Step 2: Paste the following into `check_fastq_pairs.sh`
```bash
#!/bin/bash
set -euo pipefail

INDIR="raw_data"
[[ "$(basename "$PWD")" != "raw_data" ]] && cd "$INDIR" || { echo "❌ raw_data directory not found"; exit 1; }

echo "🔍 Checking FASTQ pairings in $PWD ..."

MISSING=false
PAIRED_COUNT=0
TOTAL_COUNT=0

for R1 in *_1.fastq.gz *_R1.fastq.gz *_R1_*.fastq.gz *_001.fastq.gz; do
    [[ -f "$R1" ]] || continue
    TOTAL_COUNT=$((TOTAL_COUNT+1))
    SAMPLE=${R1%_1.fastq.gz}; SAMPLE=${SAMPLE%_R1.fastq.gz}; SAMPLE=${SAMPLE%_R1_*.fastq.gz}; SAMPLE=${SAMPLE%_001.fastq.gz}; SAMPLE=${SAMPLE%_R1_001.fastq.gz}

    if [[ -f "${SAMPLE}_2.fastq.gz" || -f "${SAMPLE}_R2.fastq.gz" || -f "${SAMPLE}_R2_*.fastq.gz" || -f "${SAMPLE}_002.fastq.gz" ]]; then
        echo "✅ $SAMPLE — paired"
        PAIRED_COUNT=$((PAIRED_COUNT+1))
    else
        echo "❌ $SAMPLE — missing R2 file"
        MISSING=true
    fi
done

echo -e "\nTotal samples checked: $TOTAL_COUNT"
echo "Correctly paired samples: $PAIRED_COUNT"
$MISSING && echo "⚠ Some samples are missing pairs. Fix before running fastp." || echo "✅ All FASTQ files are correctly paired."
```

##### Step 4: Make the script executable
```bash
chmod +x check_fastq_pairs.sh
```
##### Step 5: Run the script
```bash
./check_fastq_pairs.sh
```
> **Tip:** Ensure all R1/R2 naming conventions in your directory match the patterns used in the script.  
> You can adjust the patterns (`*_1.fastq.gz`, `*_R1.fastq.gz`, etc.) if needed.


##### Step 1:  Open nano to create a new script
```bash
nano fastq_read_length_summary.sh
```
##### Step 2: Paste the following code into nano

```bash
#!/bin/bash
set -euo pipefail

FASTQ_DIR="raw_data"
OUTDIR="csv_output"
OUTPUT_CSV="${OUTDIR}/read_length_summary.csv"

mkdir -p "$OUTDIR"

echo "Sample,R1_min,R1_max,R1_avg,R2_min,R2_max,R2_avg" > "$OUTPUT_CSV"

for R1 in "$FASTQ_DIR"/*_1.fastq.gz "$FASTQ_DIR"/*_R1.fastq.gz; do
    [[ -f "$R1" ]] || continue
    SAMPLE=$(basename "$R1" | sed -E 's/_R?1.*\.fastq\.gz//')
    R2="$FASTQ_DIR/${SAMPLE}_2.fastq.gz"
    [[ -f "$R2" ]] || R2="$FASTQ_DIR/${SAMPLE}_R2.fastq.gz"

    if [[ -f "$R2" ]]; then
        echo "Processing sample $SAMPLE"

        calc_stats() {
            zcat "$1" | awk 'NR%4==2 {len=length($0); sum+=len; if(min==""){min=len}; if(len<min){min=len}; if(len>max){max=len}; count++} END{avg=sum/count; printf "%d,%d,%.2f", min, max, avg}'
        }

        STATS_R1=$(calc_stats "$R1")
        STATS_R2=$(calc_stats "$R2")

        echo "$SAMPLE,$STATS_R1,$STATS_R2" >> "$OUTPUT_CSV"
    else
        echo "⚠ Missing R2 for $SAMPLE, skipping."
    fi
done

echo "✅ Read length summary saved to $OUTPUT_CSV"

```
##### Step 3: Save and exit nano
Press Ctrl + O → Enter (to write the file)
Press Ctrl + X → Exit nano
##### Step 4: Make the script executable
```bash
chmod +x fastq_read_length_summary.sh
```
##### Step 5: Run the script
```bash
./fastq_read_length_summary.sh
```



# 3️⃣ FASTP – Quality Control and Trimming

### Steps to Run FASTP
##### Step 1: **Open nano to create the script**
```bash
nano run_fastp.sh
```
##### Step 2: Paste the following code into nano
```bash
#!/bin/bash
set -euo pipefail

INDIR="raw_data"
OUTDIR="fastp_results_min_50"
mkdir -p "$OUTDIR"

SAMPLES=()

for R1 in "$INDIR"/*_1.fastq.gz "$INDIR"/*_R1.fastq.gz "$INDIR"/*_001.fastq.gz "$INDIR"/*_R1_001.fastq.gz; do
    [[ -f "$R1" ]] || continue

    SAMPLE=$(basename "$R1")
    SAMPLE=${SAMPLE%%_1.fastq.gz}
    SAMPLE=${SAMPLE%%_R1.fastq.gz}
    SAMPLE=${SAMPLE%%_001.fastq.gz}
    SAMPLE=${SAMPLE%%_R1_001.fastq.gz}

    if   [[ -f "$INDIR/${SAMPLE}_2.fastq.gz" ]]; then R2="$INDIR/${SAMPLE}_2.fastq.gz"
    elif [[ -f "$INDIR/${SAMPLE}_R2.fastq.gz" ]]; then R2="$INDIR/${SAMPLE}_R2.fastq.gz"
    elif [[ -f "$INDIR/${SAMPLE}_002.fastq.gz" ]]; then R2="$INDIR/${SAMPLE}_002.fastq.gz"
    elif [[ -f "$INDIR/${SAMPLE}_R2_001.fastq.gz" ]]; then R2="$INDIR/${SAMPLE}_R2_001.fastq.gz"
    else
        echo "⚠ No R2 file found for $SAMPLE — skipping."
        continue
    fi

    if [[ -f "$OUTDIR/${SAMPLE}_1.trim.fastq.gz" && -f "$OUTDIR/${SAMPLE}_2.trim.fastq.gz" ]]; then
        echo "⏩ Skipping $SAMPLE (already processed)."
        continue
    fi

    SAMPLES+=("$SAMPLE,$R1,$R2")
done

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "❌ No paired FASTQ files found in $INDIR"
    exit 1
fi

THREADS=$(nproc)
FASTP_THREADS=$(( THREADS / 2 ))

run_fastp() {
    SAMPLE=$1
    R1=$2
    R2=$3
    echo "✅ Processing sample: $SAMPLE"
    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$OUTDIR/${SAMPLE}_1.trim.fastq.gz" \
        -O "$OUTDIR/${SAMPLE}_2.trim.fastq.gz" \
        -h "$OUTDIR/${SAMPLE}_fastp.html" \
        -j "$OUTDIR/${SAMPLE}_fastp.json" \
        --length_required 50 \
        --qualified_quality_phred 20 \
        --detect_adapter_for_pe \
        --thread $FASTP_THREADS \
        &> "$OUTDIR/${SAMPLE}_fastp.log"
}

export -f run_fastp
export OUTDIR FASTP_THREADS

printf "%s\n" "${SAMPLES[@]}" | parallel -j 3 --colsep ',' run_fastp {1} {2} {3}

echo "🎉 Completed fastp for $(ls "$OUTDIR"/*_fastp.json | wc -l) samples."
```

##### Step 3: Save & exit nano
Press CTRL+O, Enter (save)
Press CTRL+X (exit)
##### Step 4: Make the script executable
```bash
chmod +x run_fastp.sh
```
##### Step 5: Activate your conda env and run
```bash
conda activate fastp_env
./run_fastp.sh
```

Count R1 trimmed files
```bash
ls -lth fastp_results_min_50/*_1.trim.fastq.gz | wc -l
```
Count R2 trimmed files
```bash
ls -lth fastp_results_min_50/*_2.trim.fastq.gz | wc -l
```

View first 10 quality lines in trimmed FASTQ

```bash
 Show first 10 quality lines from R1
```bash
echo "🔹 First 10 quality lines from ET3_S55_1 (R1):"
zcat fastp_results_min_50/ET3_S55_1.trim.fastq.gz \
| sed -n '4~4p' \
| head -n 10 \
| awk '{print "✅ " $0}'
```
 Show first 10 quality lines from R2
```bash
echo "🔹 First 10 quality lines from ET3_S55_2 (R2):"
zcat fastp_results_min_50/ET3_S55_2.trim.fastq.gz \
| sed -n '4~4p' \
| head -n 10 \
| awk '{print "✅ " $0}'

```
<details>
  <summary>🔍 How it works</summary>

- `zcat` → Decompresses the trimmed FASTQ.  
- `sed -n '4~4p'` → Prints every 4th line starting from line 4 (the quality score line of each read).  
- `head -n 10` → Shows the first 10 quality lines for quick inspection.  

</details>

Count ASCII characters in quality lines:
```bash
    Count base composition in R1
```bash
zcat fastp_results_min_50/ET3_S55_1.trim.fastq.gz \
| sed -n '4~4p' \
| awk '{
    for(i=1;i<=length($0);i++){ q[substr($0,i,1)]++ }
} END {
    for (k in q) print k, q[k]
}' \
| awk '{print "✅ Base " $1 ": " $2 " occurrences"}'
```
   Count base composition in R2
```bash
zcat fastp_results_min_50/ET3_S55_2.trim.fastq.gz \
| sed -n '4~4p' \
| awk '{
    for(i=1;i<=length($0);i++){ q[substr($0,i,1)]++ }
} END {
    for (k in q) print k, q[k]
}' \
| awk '{print "✅ Base " $1 ": " $2 " occurrences"}'
```
<details>
  <summary>🔢 How it works</summary>

- `zcat` → Decompresses trimmed FASTQ.  
- `sed -n '4~4p'` → Selects every 4th line (quality score line).  
- `awk '{for(i=1;i<=length($0);i++){q[substr($0,i,1)]++}} END{for (k in q) print k,q[k]}'` → Counts occurrences of each ASCII character in the quality scores.  

This helps quickly identify if the quality encoding is correct (usually Phred+33 for Illumina) and whether trimming improved overall quality.

</details>

For batch processing trimmed FASTQ
##### Step 1: Open a new script
```bash
nano count_trimmed_reads.sh
```
##### Step 2: Paste the following code
This script counts reads in trimmed paired-end FASTQ files and saves results to a CSV.
```bash
#!/bin/bash
set -euo pipefail

INDIR="fastp_results_min_50"
OUTDIR="csv_output"
mkdir -p "$OUTDIR"

OUTFILE="$OUTDIR/trimmed_read_counts.csv"
echo "Sample,R1_reads,R2_reads" > "$OUTFILE"
echo "📊 Counting reads in trimmed FASTQ files from '$INDIR'..."

for R1 in "$INDIR"/*_1.trim.fastq.gz "$INDIR"/*_R1.trim.fastq.gz; do
    [[ -f "$R1" ]] || continue
    SAMPLE=$(basename "$R1" | sed -E 's/_R?1.*\.trim\.fastq\.gz//')
    R2=""
    for suffix in "_2.trim.fastq.gz" "_R2.trim.fastq.gz" "_R2_*.trim.fastq.gz"; do
        [[ -f "$INDIR/${SAMPLE}${suffix}" ]] && R2="$INDIR/${SAMPLE}${suffix}" && break
    done
    R1_COUNT=$(( $(zcat "$R1" | wc -l) / 4 ))
    R2_COUNT=$([[ -n "$R2" ]] && echo $(( $(zcat "$R2" | wc -l) / 4 )) || echo "NA")
    echo "$SAMPLE,$R1_COUNT,$R2_COUNT" >> "$OUTFILE"
    echo "✅ $SAMPLE → R1: $R1_COUNT | R2: $R2_COUNT"
done

echo "🎉 All done! Read counts saved to '$OUTFILE'"


```

##### Step 3: Save and exit nano
Press Ctrl + O → Enter (to write the file)
Press Ctrl + X → Exit nano

##### Step 4: Make the script executable
```bash
chmod +x count_trimmed_reads.sh
```
##### Step 5: Run the script
```bash
./count_trimmed_reads.sh
```

read length summary on trimmed FASTQ files
##### Step 1: Open nano to create a new script
```bash
nano trimmed_fastq_read_length_summary.sh
```
##### Step 2: Paste the following code into nano
```bash
#!/bin/bash
set -euo pipefail

FASTQ_DIR="fastp_results_min_50"
OUTDIR="csv_output"
OUTPUT_CSV="${OUTDIR}/trimmed_read_length_summary.csv"

mkdir -p "$OUTDIR"

echo "Sample,R1_min,R1_max,R1_avg,R2_min,R2_max,R2_avg" > "$OUTPUT_CSV"

for R1 in "$FASTQ_DIR"/*_1.trim.fastq.gz; do
    SAMPLE=$(basename "$R1" _1.trim.fastq.gz)
    R2="${FASTQ_DIR}/${SAMPLE}_2.trim.fastq.gz"

    if [[ -f "$R2" ]]; then
        echo "Processing sample $SAMPLE"

        calc_stats() {
            zcat "$1" | awk 'NR%4==2 {
                len=length($0)
                sum+=len
                if(min==""){min=len}
                if(len<min){min=len}
                if(len>max){max=len}
                count++
            } END {
                avg=sum/count
                printf "%d,%d,%.2f", min, max, avg
            }'
        }

        STATS_R1=$(calc_stats "$R1")
        STATS_R2=$(calc_stats "$R2")

        echo "$SAMPLE,$STATS_R1,$STATS_R2" >> "$OUTPUT_CSV"
    else
        echo "⚠ Missing R2 for $SAMPLE, skipping."
    fi
done

echo "✅ Trimmed read length summary saved to $OUTPUT_CSV"

```

##### Step 3: Save and exit nano
Press Ctrl + O → Enter
Press Ctrl + X → Exit
##### Step 4: Make the script executable
```
chmod +x trimmed_fastq_read_length_summary.sh
```
##### Step 5: Run the script
```bash
./trimmed_fastq_read_length_summary.sh
```

# 4️⃣ MultiQC

### Script: Run MultiQC
##### Step 1: **Open nano to create the script `run_multiqc.sh`
```bash
nano run_multiqc.sh
```
##### Step 2: Paste the following code into nano
```bash
#!/bin/bash
set -euo pipefail

INPUT_DIR="fastp_results_min_50"
OUTPUT_DIR="multiqc/fastp_multiqc"

mkdir -p "$OUTPUT_DIR"

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory '$INPUT_DIR' does not exist!"
    exit 1
fi

multiqc "$INPUT_DIR" -o "$OUTPUT_DIR"

echo "MultiQC report generated in '$OUTPUT_DIR'."

```

##### Step 3: Save & exit nano
Press CTRL+O, Enter (save)
Press CTRL+X (exit)
##### Step 4: Make the script executable
```bash
chmod +x run_multiqc.sh
```
##### Step 5: Activate your conda env and run
```bash
conda activate multiqc_env
./run_multiqc.sh
```

# 5️⃣ Snippy

##### Step 1: Create the script
```bash
nano run_snippy.sh
```
##### Step 2: Paste the following into `run_snippy.sh`
```bash
#!/bin/bash
set -euo pipefail

# -----------------------------
# Configuration
# -----------------------------
REF="Salmonella_LT2.fasta"        # Shortened reference name
FASTP_DIR="fastp_results_min_50"
OUTDIR="snippy_results"
THREADS=8
BWA_THREADS=30
JOBS=4

mkdir -p "$OUTDIR"

# -----------------------------
# Function to run Snippy on a sample
# -----------------------------
run_snippy_sample() {
    SAMPLE="$1"
    R1="${FASTP_DIR}/${SAMPLE}_1.trim.fastq.gz"
    R2="${FASTP_DIR}/${SAMPLE}_2.trim.fastq.gz"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "⚠ Missing R1/R2 for $SAMPLE"
        return
    fi

    # Skip if VCF already exists
    if [[ -f "${OUTDIR}/${SAMPLE}.vcf" ]]; then
        echo "ℹ ${SAMPLE} already processed. Skipping."
        return
    fi

    echo "🚀 Running Snippy on sample: $SAMPLE"
    TMP_DIR="${OUTDIR}/${SAMPLE}_tmp"
    mkdir -p "$TMP_DIR"

    snippy --cpus "$THREADS" \
           --outdir "$TMP_DIR" \
           --ref "$REF" \
           --R1 "$R1" \
           --R2 "$R2" \
           --force \
           --bwaopt "-T $BWA_THREADS"

    # Move outputs to main OUTDIR
    if [[ -f "$TMP_DIR/snps.vcf" ]]; then
        mv "$TMP_DIR/snps.vcf" "${OUTDIR}/${SAMPLE}.vcf"
    fi

    for f in "$TMP_DIR"/*; do
        base=$(basename "$f")
        case "$base" in
            *.consensus.fa) mv "$f" "${OUTDIR}/${SAMPLE}.consensus.fa" ;;
            *.bam) mv "$f" "${OUTDIR}/${SAMPLE}.bam" ;;
            *.bam.bai) mv "$f" "${OUTDIR}/${SAMPLE}.bam.bai" ;;
            *.tab) mv "$f" "${OUTDIR}/${SAMPLE}.snps.tab" ;;
        esac
    done

    rm -rf "$TMP_DIR"

    [[ -f "${OUTDIR}/${SAMPLE}.vcf" ]] && echo "✅ Full VCF generated for $SAMPLE" || echo "⚠ No VCF produced for $SAMPLE"
}

# -----------------------------
# Export variables and function
# -----------------------------
export -f run_snippy_sample
export REF FASTP_DIR OUTDIR THREADS BWA_THREADS

# -----------------------------
# Run Snippy in parallel
# -----------------------------
ls "${FASTP_DIR}"/*_1.trim.fastq.gz \
    | sed 's|.*/||; s/_1\.trim\.fastq\.gz//' \
    | parallel -j "$JOBS" run_snippy_sample {}

# -----------------------------
# Compare processed samples
# -----------------------------
ls "${FASTP_DIR}"/*_1.trim.fastq.gz \
    | sed 's|.*/||; s/_1\.trim\.fastq\.gz//' | sort > fastq_samples.txt

ls "${OUTDIR}"/*.vcf 2>/dev/null \
    | sed 's|.*/||; s/\.vcf//' | sort > snippy_samples.txt

echo "FASTQ pairs count: $(wc -l < fastq_samples.txt)"
echo "Snippy outputs count: $(wc -l < snippy_samples.txt)"

if diff fastq_samples.txt snippy_samples.txt >/dev/null; then
    echo "✅ All FASTQ pairs have corresponding Snippy results."
else
    echo "⚠ Missing samples detected:"
    diff fastq_samples.txt snippy_samples.txt || true
fi

rm -f fastq_samples.txt snippy_samples.txt

echo "🎯 All steps completed!"
echo "Snippy results are in: ${OUTDIR}/"


```
##### Step 3: Save and exit nano
Press Ctrl + O, then Enter (save)
Press Ctrl + X (exit)

##### Step 4: Make the script executable
```bash
chmod +x run_snippy.sh
```
##### Step 5: Activate environment and run
```bash
conda activate snippy_env
./run_snippy.sh
```

Check Snippy VCFs Before Variant Filtering
- `tb_variant_filter` relies on a correctly formatted VCF with the `#CHROM` line to parse variants.  
- VCFs missing this line will **fail** during filtering, causing errors like:  
  `Missing line starting with "#CHROM"`.  
- Running this check before variant filtering saves time and ensures downstream analysis runs smoothly.

This script ensures that all VCF files generated by **Snippy** contain the required `#CHROM` header line.
```bash
#!/bin/bash
OUTDIR="snippy_results"

echo "Checking Snippy VCFs in $OUTDIR ..."

for vcf in "$OUTDIR"/*.vcf; do
    SAMPLE=$(basename "$vcf")
    if grep -q "^#CHROM" "$vcf"; then
        echo "✅ $SAMPLE contains #CHROM line"
    else
        echo "⚠ $SAMPLE is missing #CHROM line"
    fi
done
```

# 6️⃣ Qualimap BAM QC


##### Step 1: Create the script
```bash
nano run_qualimap.sh
```
#####  Step 2: Paste the following into `run_qualimap.sh`
```bash
#!/bin/bash
set -euo pipefail

SNIPPY_DIR="snippy_results"
QUALIMAP_OUT="qualimap_reports"
mkdir -p "$QUALIMAP_OUT"

for bam in "$SNIPPY_DIR"/*.bam; do
    sample=$(basename "$bam" .bam)
    echo "Running Qualimap BAM QC for sample: $sample"

    outdir="${QUALIMAP_OUT}/${sample}"
    mkdir -p "$outdir"

    qualimap bamqc \
        -bam "$bam" \
        -outdir "$outdir" \
        -outformat pdf:html \
        --java-mem-size=4G
done
```


##### Step 3: Save and exit nano

Press Ctrl + O, then Enter (save)
Press Ctrl + X (exit)
##### Step 4: Make the script executable
```bash
chmod +x run_qualimap.sh
```
##### Step 5: Activate environment and install GNU Parallel into your `qualimap_env`:
```bash
conda activate qualimap_env
```
##### Step 6: run:
```bash
./run_qualimap.sh
```

# 7️⃣ MultiQC after Qualimap

### Run MultiQC on Qualimap outputs

##### Step 1: **Open nano to create the script `run_multiqc_qualimap.sh`
```bash
nano run_multiqc_qualimap.sh
```
##### Step 2: Paste the following code into nano
```bash
#!/bin/bash
set -euo pipefail

INPUT_DIR="qualimap_reports"
OUTPUT_DIR="multiqc/qualimap_multiqc"

mkdir -p "$OUTPUT_DIR"

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory '$INPUT_DIR' does not exist!"
    exit 1
fi

multiqc "$INPUT_DIR" -o "$OUTPUT_DIR"

echo "MultiQC report generated in '$OUTPUT_DIR'."

```


##### Step 3: Save & exit nano
Press CTRL+O, Enter (save)
Press CTRL+X (exit)
##### Step 4: Make the script executable
```bash
chmod +x run_multiqc_qualimap.sh
```
##### Step 5: Activate your conda env and run
```bash
conda activate multiqc_env
./run_multiqc_qualimap.sh
```


# 1️⃣0️⃣ BCFTools Consensus Generation

<details>
<summary>🧬 Generate Sample-Specific Consensus Sequences</summary>

After filtering VCFs with **tb_variant_filter**, we generate **consensus FASTA sequences** for each sample. These sequences represent the **full genome of each isolate**, including only **high-confidence variants** relative to the reference genome (*H37Rv*).

### Why consensus sequences matter
- Provide a **single representative genome** per sample for downstream analyses.  
- Used in **phylogenetic reconstruction**, outbreak investigation, and comparative genomics.  
- Incorporate **only reliable SNPs and indels**, minimizing noise from sequencing errors.  
- Standardize genome representations across multiple isolates, ensuring comparability.  

</details>
---

##### Step 1: Compress and index each filtered VCF
```bash
for vcf in tb_variant_filter_results/*.vcf; do
   
    if [ ! -s "$vcf" ]; then
        echo "Skipping empty VCF: $vcf"
        continue
    fi

    gz_file="${vcf}.gz"
    echo "Compressing $vcf → $gz_file"
    bgzip -c "$vcf" > "$gz_file"

    echo "Indexing $gz_file"
    bcftools index "$gz_file"
done
```
##### Step 2: Create a script to generate consensus sequences
```bash
nano generate_consensus_all.sh
```
##### Step 3: Paste this code
```bash
#!/bin/bash
set -euo pipefail

CURDIR=$(pwd)
VCFDIR="$CURDIR/tb_variant_filter_results"
OUTDIR="$CURDIR/consensus_sequences"
mkdir -p "$OUTDIR"

for vcf in "$VCFDIR"/*.vcf; do
    sample=$(basename "$vcf" .vcf)

    if [ $(grep -v '^#' "$vcf" | wc -l) -eq 0 ]; then
        echo "$sample VCF is empty. Skipping."
        continue
    fi

    gz_file="${vcf}.gz"
    bgzip -c "$vcf" > "$gz_file"
    bcftools index "$gz_file"
    bcftools consensus -f "$CURDIR/H37Rv.fasta" "$gz_file" | sed "1s/.*/>$sample/" > "$OUTDIR/${sample}.consensus.fasta"
done

```
<details>
<summary>🧬 VCF-to-Consensus Script Explanation</summary>

- `#!/bin/bash` → Run script with Bash.  
- `set -euo pipefail` → Exit on errors or undefined variables.  
- `CURDIR=$(pwd)` → Save current working directory.  
- `VCFDIR="$CURDIR/tb_variant_filter_results"` → Folder with filtered VCFs.  
- `OUTDIR="$CURDIR/consensus_sequences"` → Folder for consensus FASTA sequences.  
- `mkdir -p "$OUTDIR"` → Ensure output directory exists.  
- `for vcf in "$VCFDIR"/*.vcf; do ... done` → Loop through all filtered VCF files.  
- `sample=$(basename "$vcf" .vT2.fasta 
cf)` → Extract sample name.  
- `bgzip -c "$vcf" > "$vcf.gz"` → Compress VCF with bgzip.  
- `bcftools index "$vcf.gz"` → Index compressed VCF.  
- `bcftools consensus -f "$CURDIR/H37Rv.fasta" "$vcf.gz" | sed "1s/.*/>$sample/" > "$OUTDIR/${sample}.consensus.fasta"` → Generate consensus FASTA and replace header with sample name.  
- `echo "✅ $sample consensus generated"` → Confirmation per sample.  
- `echo "🎉 All consensus sequences saved in $OUTDIR"` → Final message.  

**⚠ Note:** Activate the `tb_consensus_env` before running this script.

</details>


##### Step 4: Save and exit nano
Press Ctrl + O → Enter (to write the file)
Press Ctrl + X → Exit nano

##### Step 5: Make the script executable
```bash
chmod +x generate_consensus_all.sh
```
##### Step 6: Run the script
```bash
conda activate tb_consensus_env
./generate_consensus_all.sh
```

# 1️⃣1️⃣ Check Consensus FASTA Lengths

After generating consensus sequences, it's important to **verify the genome length** for each sample.  
This ensures no sequences are truncated or incomplete due to missing coverage or filtering.

---

### 📏 Calculating Consensus Genome Lengths

We can check the length of each consensus FASTA sequence to ensure completeness and consistency.  
This helps verify that consensus sequences cover the full *M. tuberculosis* genome (~4.4 Mbp) and can reveal missing regions.

### Rename the FASTA files

```bash
#!/bin/bash
FASTA_DIR="consensus_sequences"

for f in "$FASTA_DIR"/*.filtered.consensus.fasta; do
    mv "$f" "${f/.filtered.consensus/}"
done

echo "✅ All consensus FASTA files have been renamed to .fasta."

```
###  Update headers inside the FASTA files
```bash
#!/bin/bash
FASTA_DIR="consensus_sequences"

for f in "$FASTA_DIR"/*.fasta; do
    sample=$(basename "$f" .fasta)
    awk -v s="$sample" '/^>/{print ">" s; next} {print}' "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"
    echo "✅ Updated header in: $(basename "$f")"
done
echo "🎉 All FASTA headers have been successfully updated."
```

### Using `grep` and `wc`
We remove the FASTA headers (`>` lines) and count the remaining nucleotides to get the total genome length:

```bash
for f in consensus_sequences/*.fasta; do
    sample=$(basename "$f")
    length=$(grep -v ">" "$f" | tr -d '\n' | wc -c)
    echo "$sample : $length bp"
done
```

to save the result in csv file 
```bash
#!/bin/bash
FASTA_DIR="consensus_sequences"
OUTDIR="csv_output"
mkdir -p "$OUTDIR"  
OUTPUT_CSV="${OUTDIR}/consensus_lengths.csv"

echo "Sample,Length_bp" > "$OUTPUT_CSV"

for f in "$FASTA_DIR"/*.fasta; do
    sample=$(basename "$f" .fasta)
    length=$(grep -v ">" "$f" | tr -d '\n' | wc -c)
    echo "$sample,$length" >> "$OUTPUT_CSV"
done

echo "✅ Consensus genome lengths saved to $OUTPUT_CSV"

```




# 1️⃣2️⃣ Multiple Sequence Alignment with MAFFT


##### Step 1: Merge all consensus FASTAs
   
```bash
cat consensus_sequences/*.fasta > consensus_sequences/all_consensus.fasta
```
##### Step 2: Paste the script into nano
```bash
#!/bin/bash
mkdir -p mafft_results
mafft --auto --parttree consensus_sequences/all_consensus.fasta > mafft_results/aligned_consensus.fasta
```

##### Step 3: Verify the alignment

A. Quickly inspect the top of the aligned FASTA:
```bash
head mafft_results/aligned_consensus.fast
```
B. Quick visual inspection with `less`
```bash
fold -w 100 mafft_results/aligned_consensus.fasta | less
```
C. Checking MAFFT alignment output to 
- Ensure the alignment file **exists** and contains all intended sequences.  
- Verify the **number of sequences** matches expectations.  
- Check **sequence lengths** (total, min, max, average) to confirm proper alignment.  
- Detect if sequences have **varying lengths**, which may indicate misalignment or excessive gaps.  
- A **good alignment** is essential for accurate phylogenetic tree inference with IQ-TREE or other downstream analyses.

```bash
#!/bin/bash
set -euo pipefail

FILE="mafft_results/aligned_consensus.fasta"

if [[ ! -f "$FILE" ]]; then
    echo "❌ $FILE not found!"
    exit 1
fi

SEQ_COUNT=$(grep -c ">" "$FILE")
LENGTHS=($(grep -v ">" "$FILE" | awk 'BEGIN{RS=">"} NR>1{print length($0)}'))
TOTAL_LENGTH=$(IFS=+; echo "$((${LENGTHS[*]}))")
MIN_LENGTH=$(printf "%s\n" "${LENGTHS[@]}" | sort -n | head -n1)
MAX_LENGTH=$(printf "%s\n" "${LENGTHS[@]}" | sort -n | tail -n1)
AVG_LENGTH=$(( TOTAL_LENGTH / SEQ_COUNT ))

echo "✅ $FILE exists"
echo "Sequences: $SEQ_COUNT"
echo "Total length (bp): $TOTAL_LENGTH"
echo "Min length: $MIN_LENGTH"
echo "Max length: $MAX_LENGTH"
echo "Avg length: $AVG_LENGTH"

if [[ $MIN_LENGTH -eq $MAX_LENGTH ]]; then
    echo "✅ All sequences have equal length, alignment looks good."
else
    echo "⚠️ Sequence lengths vary; check for gaps or misalignment."
fi

```
D. Using `seqkit` stats (recommended)
seqkit is a fast toolkit for FASTA/Q file summaries. It gives a detailed report of sequences in a file:
```bash
conda activate seqkit_env
seqkit stats mafft_results/aligned_consensus.fasta
```
E. Check for gaps / alignment columns

```bash
grep -v ">" mafft_results/aligned_consensus.fasta \
| awk '{print length($0)}' \
| sort -n \
| uniq -c \
| awk -v L=60 '{print ($2<L?"⚠️ ":"✅ ") $1 " sequences of length " $2 " bp"}'
```

F. Compute pairwise identity
```bash
awk '/^>/{if(seqlen){print seqlen}; seqlen=0; next} {seqlen+=length($0)} END{print seqlen}' mafft_results/aligned_consensus.fasta \
| sort -n \
| uniq -c \
| awk -v L=60 '{print ($2<L?"⚠️ ":"✅ ") $1 " sequences of length " $2 " bp"}'
```
G. Use AMAS (Alignment Manipulation and Summary)
AMAS is a Python tool to summarize alignments:
```bash
conda activate amas_env
AMAS.py summary -f fasta -d dna -i mafft_results/aligned_consensus.fasta
```
H. Use aliview or MEGA for GUI inspection
Load the FASTA alignment in AliView, MEGA, or Geneious.
Advantages:
  Can visually check gaps, conserved regions, and misaligned sequences.
  Highlight sequences that differ significantly.


# 1️⃣3️⃣ IQtree

 Steps 
##### Step 1: activate iqtree environment
```bash
conda activate iqtree_env
```
##### Step 2: Run the following script
```bash
mkdir -p iqtree_results

iqtree2 -s mafft_results/aligned_consensus.fasta \
        -m GTR+G \
        -bb 1000 \
        -nt 4 \
        -o SRR10828835 \
        -pre iqtree_results/aligned_consensus
```



