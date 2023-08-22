**Phylotranscriptomic: phylogenetic approach from *De-novo* transcriptome dataset without reference genome dataset.**

For more details see the phylo transcriptomic papers below.

Kashimoto et al. **(2023)**, *Anemonefish are better taxonomists than humans:CURBIO19572 **Current biology** will appear online and with the issue on Monday, September 25.

[Kashimoto et al. **(2022)**, *Transcriptomes of Giant Sea Anemones from Okinawa as a Tool for Understanding Their Phylogeny and Symbiotic Relationships with Anemonefish:zs210111 **Zoological Science**](https://www.researchgate.net/publication/360318669_Transcriptomes_of_Giant_Sea_Anemones_from_Okinawa_as_a_Tool_for_Understanding_Their_Phylogeny_and_Symbiotic_Relationships_with_Anemonefish)
 
Installation\
1 Trimmomatic\
2 Trinity\
3 CD-HIT\
4 blast\
5 [get_gc_content.pl](https://github.com/jmeneghin/perl-for-reysenbach-lab)\
6 Transdecoder\
7 Orthofinder\
8 fasttree\
9 BUSCO\
10 iqtree2


**1 Trimmomatic: Quality trimming and adapter clipping**\

$java -jar /Version0.39binary/Trimmomatic-0.39/trimmomatic-0.39.jar PE R1_001.fastq R1_002.fastq output_forward_paired.fastq output_forward_unpaired.fastq output_reverse_paired.fastq output_reverse_unpaired.fastq ILLUMINNACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

**2 Trinity for *De-novo* assembely for array job**\

$module load python/3.7.3
/trinityrnaseq-v2.11.0/Trinity --seqType fq --left ${file}_R1_paired.fastq --right ${file}_R2_paired.fastq --max_memory 200G --min_kmer_cov 2 --min_contig_length 300 --bflyCPU 30 --output trinity_ID295_${file} --bflyHeapSpaceMax 1G --bflyGCThreads 2 --SS_lib_type RF --no_normalize_reads

You will get the final assembled reads **Trinity.fasta**.

**3 Gene(isoform) name change** You can delete TRINITY_DN change to tr or your preffered name and delete the information after the blank space.

$sed "s/>TRINITY_DN/>${file}_tr/" ${file}_Trinity.fasta > ${file}_Trinity.fasta.mod01
awk -F" " '{print $1}' ${file}_Trinity.fasta.mod01 > ${file}_Trinity.fasta.mod02

**4 CD-HIT**

$cd-hit-est -i ${file}_Trinity.fasta.mod02 -o ${file}.fasta.db98 -c 0.98 -n 11 -M 16000 -d 50 -T 8


**5 Concatanation - if you need to sort your target species from symbiotic organisms**

$for i in *_Trinity.fasta.mod02; do cat ${i} >> all_anemones_tr.fasta; done
	
**Blast in case you want to purify your target species transcriptome data from symbiotic organisms**\
In my paper, I separated symbiodinium and alga sequences from sea anemone dataset sue to we obtained the sample from tentacles.
	
	IN_DIR="/CDhit_210324"
	OUT_DIR="/blast_outfolder"
	BLAST_DIR="/ncbi-blast-2.10.1/bin/"
	
	# if the database in somewhere else
	export BLASTDB="/blast/dbs"
	
	 ${BLAST_DIR}/blastn -task blastn -query ${IN_DIR}/all_anemones_tr.fasta -db ${BLASTDB}/symb-A-B1-C_transcriptomes \
	 -num_descriptions 5 -num_alignments 5 -line_length 100 \
	 -out ${OUT_DIR}/all_anemones_tr_vs_symABCtr_blastn_1e-20 -evalue 1e-20 -outfmt 0 -num_threads 64

**GC content in case you want to purify your target species transcriptome data from symbiotic organisms**
#It takes bit time
nohup ./get_gc_content.pl ${file}.fasta.db98 ${file}.GC & 

#Check if you sorted correctly for the GC content
cat Trinity.fasta.mod02.GC | grep '>' | wc
cat Trinity.fasta | grep '>' | wc
# the number should be same between Trinity.fasta.mod02.GC file and Trinity.fasta file

# Concatanate all sample that you have analysised
for i in *_Trinity.fasta.mod02.GC; do cat ${i} >> all_anemones_tr.fasta.GC

# Get the GC content cloum that only you need to use for further steps
for i in all_anemones_tr.fasta.GC; do awk -F"\t" '{if (NR>1) {printf "%s\t%s\n",$1,$2}}' ${i} >> all_anemones_tr.GC; done

# Sort the result from blastn around 5-10 min 
./parser_best_hit_06_wo_description all_anemones_tr_vs_symABCtr_blastn_1e-20 > all_anemones_tr_vs_symABCtr_blastn_1e-20.bhtab

# Check the count is same as GC file and the line number should be same
wc -l SRA_tr_vs_symABCtr_blastn_1e-20.bhtab
wc -l all_anemones_tr.GC

# sort the coloum blast
$ awk -F"\t" '{ if (NR==1) {print "Q_Num\tQuery_Name\tQ_Len\tHit_Name\tE_val"} else if (NR>1) {printf"%s\t%s\t%s\t%s\t%s\n",NR-1,$2,$3,$4,$7}}' all_anemones_tr_vs_symABCtr_blastn_1e-20.bhtab > all_anemones_tr_vs_symABCtr_blastn_1e-20-ALL.bhtab.BLAST


# check the line number. Blast file has one more line than GC file
$ wc -l all_anemones_tr.GC
749635 all_anemones_tr.GC
$ wc -l all_anemones_tr_vs_symABCtr_blastn_1e-20-ALL.bhtab.BLAST
749636 all_anemones_tr_vs_symABCtr_blastn_1e-20-ALL.bhtab.BLAST


# vim delete first line on blast file
$ vi all_anemones_tr_vs_symABCtr_blastn_1e-20-ALL.bhtab.BLAST 
$ head all_anemones_tr_vs_symABCtr_blastn_1e-20-ALL.bhtab.BLAST 
1 Equa_10_S27_tr17623_c0_g1_i1 2196 no hits 0
2 Equa_10_S27_tr17623_c0_g2_i1 540 no hits 0
3 Equa_10_S27_tr17623_c0_g3_i1 595 no hits 0
4 Equa_10_S27_tr17623_c0_g4_i1 2501 no hits 0
5 Equa_10_S27_tr17623_c0_g4_i2 2416 no hits 0
6 Equa_10_S27_tr17663_c0_g1_i1 395 no hits 0
7 Equa_10_S27_tr17663_c0_g2_i1 348 no hits 0
8 Equa_10_S27_tr17628_c1_g1_i7 7552 no hits 0
9 Equa_10_S27_tr17628_c1_g1_i2 4462 no hits 0
10 Equa_10_S27_tr17628_c1_g1_i5 3572 no hits 0


# To conclude dataset I used Paste command
$ paste all_anemones_tr_vs_symABCtr_blastn_1e-20-ALL.bhtab.BLAST ../all_anemones_tr.GC > all_anemones_ABC-BLAST_GC.TAB


# To sort
./genome2tableNCBI.sh all_anemones_tr.fasta all_anemones_tr.fasta.TAB &
[1] 4941


# Based on GC content separation
＃Sea anemone
$ awk -F"\t" 'FNR==NR { a[FNR""] = $4; b[FNR""] =$7; next } {if (a[FNR""]=="no hits" && b[FNR""]<60) {print $0}}' all_anemones_ABC-BLAST_GC.TAB all_anemones_tr.fasta.TAB > all_anemones_tr.fasta.TAB.woHIT_60

# algae
awk -F"\t" 'FNR==NR { a[FNR""] = $4; b[FNR""] =$7; next } {if (a[FNR""]!="no hits" && b[FNR""]>=40) {print $0}}' all_anemones_ABC-BLAST_GC.TAB all_anemones_tr.fasta.TAB > all_algae_tr.fasta.TAB.wHIT_40

# Check identifier

# Anemone　around ６5-70％
$ cat all_anemones_tr.fasta.TAB.woHIT_60 | awk -F"_" '{print $1}' - | sort -s - | uniq -c - | head
 205062 >Equa
  87433 >G2
 428653 >Hcri
  47959 >Hhem
 105811 >MdorSRA
 158118 >Shad

# Alga　around ３０ -35％
$ cat all_algae_tr.fasta.TAB.wHIT_40| awk -F"_" '{print $1}' - | sort -s - | uniq -c - | head
  89439 >Equa
  41502 >G2
 215167 >Hcri
  40996 >Hhem
  77985 >MdorSRA
  77494 >Shad


#Unidentified 
$ cat all_anemones_tr.fasta.TAB | awk -F"_" '{print $1}' - | sort -s - | uniq -c - | head
317595 >Equa
 139256 >G2
 688981 >Hcri
  97220 >Hhem
 198134 >MdorSRA
 253812 >Shad

# create only Sea anemone dataset file
for SP in {Equa_9,Equa_10}; do grep ">${SP}" all_anemones_tr.fasta.TAB.woHIT_60 > ${SP}_anemone_nucle.TAB; done

# alga file
for SP in {Equa,Haur,Hcri,Hmag,Sgig,Shad,Smer}; do grep ">${SP}" all_algae_tr.fasta.TAB.wHIT_40 > ${SP}_alga_nucle.TAB; done

# The ID should be separate paragraph from the reads
for i in `cat ./list`;do fold -s ${i}_anemone_nucle.TAB > ${i}_anemone_nucle.fasta ; done
for i in `cat ./list`;do fold -s ${i}_alga_nucle.TAB > ${i}_alga_nucle.fasta ; done

**6 quarity check** 
#BUSCO 

#Longest transcript & Totall Transcripts (Maximum length of sequense= Longest Transcript  )

**7 Protein prediction**
#Step1
/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -S -m100  -t ${file}_anemone_nucle.fasta
#Step2
/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict --single_best_only -t ${file}_anemone_nucle.fasta

# Gene ID change
awk -F" " '{print $1}' ${file}_anemone_nucle.fasta.transdecoder.pep > ${file}_anemone.aa


**7 Orthofinder**
mv *{file}_anemone.aa Ortho_anemone

conda activate fasttree
/OrthoFinder-v2.4.0/orthofinder -f /Ortho_anemone -M msa

# Output location of the Orthofinder
/OrthoFinder/Results_date/MultipleSequenceAlignments/SpeciesTreeAlignment.fa

**8 Phylogenetic Tree reconstraction**
iqtree2 -s SpeciesTreeAlignment.fa -m MFP -bb 1000 -alrt 1000 -nt AUTO

**9 BUSCO**
busco --mode protein --in ${file}_anemone.aa --out BUSCO_metazoa_anemone_aa --lineage /metazoa_odb10
