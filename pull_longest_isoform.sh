# Path to your input files
GFF_DIR="/core/projects/EBP/Wegrzyn/EVOME/thymallus_arcticus/03_pangenome/orthofinder_salmonid_family/protiens/longest_isoform/run_longest/tmp_grayling"
SINGULARITY_IMG="/isg/shared/databases/nfx_singularity_cache/depot.galaxyproject.org-singularity-agat-1.2.0--pl5321hdfd78af_0.img"

# Loop over every GFF file in the directory
for gff in "$GFF_DIR"/*.gff; do
    base=$(basename "$gff" .gff)
    fna="$GFF_DIR/${base}.fna"

    if [[ ! -f "$fna" ]]; then
        echo " No matching FASTA file found for $base — skipping."
        continue
    fi

    echo "🔹 Processing $base ..."

    # Step 1: Replace transcript → mRNA
    sed 's/transcript/mRNA/g' "$gff" > "${base}_mRNA.gff"

    # Step 2: Keep longest isoform
    singularity exec "$SINGULARITY_IMG" agat_sp_keep_longest_isoform.pl \
        -gff "${base}_mRNA.gff" -o "${base}_longest_isoform.gff"

    # Step 3: Extract protein sequences
    singularity exec "$SINGULARITY_IMG" agat_sp_extract_sequences.pl \
        -g "${base}_longest_isoform.gff" \
        -f "$fna" \
        -p -o "${base}_LI.pep"

done
