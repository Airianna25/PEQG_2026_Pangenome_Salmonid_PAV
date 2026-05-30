# PEQG 2026 Pangenome_Salmonid PAV Tutorial
This is a quick tutorial on finding gene presence absence variants (PAV) in family level analysis with alignment methods (BLAST and miniprot). 

## Introduction ##

This is one approach you can take to find putative gene PAVs. The first step is a quick BLAST/diamond run to understand sequence similarity, and then a protien to genome alignment (miniprot) to check if genes with low sequence similarity might exist in your other genomes. In my study, I am focused on finding gene PAVs in salmonids, particularly those that might exist in one genus (Thymallus) for which I have three assemblies, compared to other salmonids. Salmonids are a family of both anadromous and non-anadromous fishes.


If you want to do this on your own data you would need:

1. Protien files (.faa, .pep) for all of your indviduals - generated from annotation
2. Assemblies for each (.fasta)
3. A database created from the protiens (.dmnd) which I'll describe how to create

### Input data ###

I have uploaded all the protien files, and assemblies here for use. Here are the species (I have one assembly and annotation for each):
1. Arctic grayling (*Thymallus arcticus*)
2. European grayling (*Thymallus thymallus*)
3. Amur grayling (*Thymallus grubii*)
4. Northern pike (*Esox lucius*) - outgroup species
5. Cisco (*Coregonus artedi*)
6. Atlantic salmon (*Salmo salar*)
7. Arctic char (*Salvelinus alpinus*)
8. Brook trout (*Salvelinus fontinalus*)
9. Lake trout (*Salvelinus namaycush*)
10. Rainbow trout (*Oncorhynchus mykiss*)
11. Cutthrout trout (*Oncorhynchus clarkii*)
12. Coho salmon (*Oncorhynchus tshawytscha*)
13. Sockeye salmon (*Oncorhynchus nerka*)
14. Chum salmon (*Oncorhynchus keta*)
15. Pink salmon (*Onchorhynchus gorbsucha*)



### Required Software ###

- miniprot
- diamond (faster than BLAST)
- seqkit 

### Installation ###

miniprot
```
git clone https://github.com/lh3/miniprot
cd miniprot && make
```

diamond
```
wget http://github.com/bbuchfink/diamond/releases/download/v2.2.1/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
```
seqkit
```
conda install -c bioconda seqkit
```

### Step 1: Assess Sequence Similarity ###


First, we want to create a database of all the species we are going to compare to. In this case, I am BLASTing all *Thymallus* protiens against the protiens of all the other salmonids. This will give us an idea of which genes may or may not align, and provide insight into finding putative PAVs.
Download the ```all_non_grayling_salmonids.pep ``` file, which we will use to create the DIAMOND database.

Create the database:
```
#everything except graylings
diamond makedb \
  --in all_non_grayling_salmonids.pep \
  -d non_grayling_salmonids
```

The output of this is our database: ```non_grayling_salmonids.dmnd ```

Next, we will perform the sequence similarity search. ```grayling_all_query.pep``` is all of the ***Thymallus***, or grayling, protiens. It will align all of these to every other salmonid. 

```
#blast all grayling protiens against all other salmonids
diamond blastp \
    -q grayling_all_query.pep \
    -d non_grayling_salmonids.dmnd \
    -o all_grayling_against_salmonid_pangenome.txt \
    -e 1e-10 \
    --more-sensitive \
    --max-target-seqs 50 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue
```

![blast results](grayling_blast_results.png)

We now have our sequence similarity results. We will compare all the grayling genes that went into the search, with the BLAST results to _find genes that did not align whatsoever_. 

The steps for this include:
1. Pulling all grayling IDs (from ```grayling_all_query.pep ``` into a list (.txt)
2. Pull the first column of the blast (grayling IDs that did align).
3. Compare and output genes that didn't align

```
grep ">" grayling_all_query.pep > all_grayling_ids.txt
cut -f1 all_grayling_against_salmonid_pangenome.txt | sort | uniq > blast_query_ids.txt
comm -23 all_grayling_ids.txt blast_query_ids.txt > grayling_nohit_genes.txt
```
The last line will output a file that contains the genes that are in the original search set, but not in the BLAST results, meaning they didn't align. 
The output is just a list of genes, so I won't show it, but there are 377. These are our putative gene candidates. Next, we will pull these gene IDs from our origianal protein file to have a new, candidate gene protien set. 

```
seqkit grep -f grayling_nohit_genes.txt  grayling_all_query.pep > best_candidate_pav_genes.pep 
```

Now you have genes for miniprot - the next section!


### Step 2: Protien to Genome Alignmnet ###

We need to align our ```best_candidate_pav_genes.pep``` file to all of our salmonid genomes. Even if they didn't align, we still need to check with a protein to genome alignment. We will do this against all other salmonids, and the input set to compare grayling vs. grayling alignments. Specifically, I was interested to see if there were any genes present all graylings (Arctic, European, and Amur) but not in the other salmonids.

```
atlantic_salmon=/core/projects/EBP/Wegrzyn/EVOME/thymallus_arcticus/02_Analysis/compare_genomes/atlantic_salmon/atlantic_salmon_29chrs_renamed.fasta
euro_grayling=/core/projects/EBP/Wegrzyn/EVOME/thymallus_arcticus/02_Analysis/compare_genomes/european_grayling/european_grayling_chrs_renamed.fasta
arctic_grayling=/core/projects/EBP/Wegrzyn/EVOME/thymallus_arcticus/02_Analysis/repeats/repeatmasker_out/hap1_male_hifiasm_haphic_masked_renamed.fasta
amur_grayling=/core/projects/EBP/Wegrzyn/EVOME/thymallus_arcticus/02_Analysis/synteny/grayling_genus_synteny/genomes/amur_grayling_46chrs_renamed.fasta
arctic_char=/core/projects/EBP/Wegrzyn/EVOME/thymallus_arcticus/02_Analysis/compare_genomes/arctic_char/arctic_char_39_chrs_renamed.fasta
chum_salmon=/core/projects/EBP/Wegrzyn/EVOME/thymallus_arcticus/02_Analysis/compare_genomes/chum_salmon/chum_salmon_37chr_renamed.fasta
lake_trout=/core/projects/EBP/Wegrzyn/EVOME/thymallus_arcticus/02_Analysis/compare_genomes/lake_trout/lake_trout_42_chroms_renamed.fasta
pink_salmon=/core/projects/EBP/Wegrzyn/EVOME/thymallus_arcticus/02_Analysis/compare_genomes/pink_salmon/pink_salmon_27_chrs_renamed.fasta
rainbow_trout=/core/projects/EBP/Wegrzyn/EVOME/thymallus_arcticus/02_Analysis/compare_genomes/rainbow_trout/rainbow_trout_33chr_renamed.fasta
sockeye_salmon=/core/projects/EBP/Wegrzyn/EVOME/thymallus_arcticus/02_Analysis/compare_genomes/sockeye_salmon/sockeye_salmon_chroms_only_renamed.fasta
brook_trout=/core/projects/EBP/Wegrzyn/EVOME/thymallus_arcticus/02_Analysis/compare_genomes/brook_trout/brook_trout_42_chrs_renamed.fasta
cutthroat_trout=/core/projects/EBP/Wegrzyn/EVOME/thymallus_arcticus/02_Analysis/compare_genomes/cutthroat_trout/cutthroat_trout_33chrs_renamed.fasta
coregonus=/core/projects/EBP/Wegrzyn/EVOME/thymallus_arcticus/02_Analysis/compare_genomes/coregonus_artedi/coregonus_artedi_chrs_renamed.fasta
coho_salmon=/core/projects/EBP/Wegrzyn/EVOME/thymallus_arcticus/02_Analysis/compare_genomes/coho_salmon/coho_salmon_30_chrs_renamed.fasta
chinook_salmon=/core/projects/EBP/Wegrzyn/EVOME/thymallus_arcticus/02_Analysis/compare_genomes/chinook_salmon/chinook_34_chrs_renamed.fasta


miniprot -t8 --gff $atlantic_salmon best_candidate_pav_genes.pep > atlantic_grayling_pav_validation.gff
miniprot -t8 --gff $amur_grayling best_candidate_pav_genes.pep > amur_grayling_pav_validation.gff
miniprot -t8 --gff $arctic_char best_candidate_pav_genes.pep > arctic_char_grayling_pav_validation.gff
miniprot -t8 --gff $chum_salmon best_candidate_pav_genes.pep > chum_salmon_grayling_pav_validation.gff
miniprot -t8 --gff $lake_trout best_candidate_pav_genes.pep > lake_trout_grayling_pav_validation.gff
miniprot -t8 --gff $pink_salmon best_candidate_pav_genes.pep > pink_salmon_grayling_pav_validation.gff
miniprot -t8 --gff $rainbow_trout best_candidate_pav_genes.pep > rainbow_trout_grayling_pav_validation.gff
miniprot -t8 --gff $sockeye_salmon best_candidate_pav_genes.pep > sockeye_salmon_grayling_pav_validation.gff
miniprot -t8 --gff $brook_trout best_candidate_pav_genes.pep > brook_trout_grayling_pav_validation.gff
miniprot -t8 --gff $cutthroat_trout best_candidate_pav_genes.pep > cutthroat_trout_pav_validation.gff
miniprot -t8 --gff $coregonus best_candidate_pav_genes.pep > coregonus_grayling_pav_validation.gff
miniprot -t8 --gff $coho_salmon best_candidate_pav_genes.pep > coho_salmon_pav_validation.gff
miniprot -t8 --gff $chinook_salmon best_candidate_pav_genes.pep > chinook_salmon_pav_validation.gff

#grayling to grayling
miniprot -t8 --gff $arctic_grayling best_candidate_pav_genes.pep > arctic_grayling_pav_validation.gff
miniprot -t8 --gff $euro_grayling best_candidate_pav_genes.pep > european_grayling_pav_validation.gff

cat *gff > all_species_miniprot_hits.gff
 
```
Here is what it looks like:

![miniprot results](miniprot_results.png)

From left to right, the first column where the hit was identified, the location of the alignment, the percent identity, and the gene ID that aligned.

Now, we want to filter by the lowest percent identity alignments. Again, these are genes that didn't align during BLAST, and have low protien to genome alignments in the other salmonids, increasing the confidence that they could be PAVs.

Filter by the lowest alignments (by percent identity):

```
awk -F'\t' '$3=="mRNA" {
    match($9,/Identity=([0-9.]+)/,a)
    match($9,/Target=([^ ]+)/,b)
    print a[1] "\t" b[1]
}' all_species_miniprot_hits.gff  | sort -g -k1,1 | head -20 > lowest_alns_pav.gff

```

Results:

![lowest_alns](lowest_alns.png)

There are 20 candidates here. As a sanity check, we can grep these IDs and make sure the alignments are all low percent identity in all non-grayling salmonids. Since we have a GFF for each we can check. This one ```egapxtmp_007423-R1_arctic_grayling```  is low across the board after independently checking each species GFF.

Lets just test that ID for now:

```
grep "egapxtmp_034625-R1_arctic_grayling" *.gff > check_all_species_test_candidate.gff

```
[Checking results output] (check_all_alns.gff)

What we see, is that we have very high quality alignmnets to the arctic grayling (as expected!), european grayling, and the amur grayling. In both the rainbow trout and coho salmon, we see an alignment rate of ~45, but no alignments whatsoever to other species.

Now you have identified a potential PAV!





