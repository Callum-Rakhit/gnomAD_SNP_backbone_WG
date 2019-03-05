#!/usr/bin/env sh

# Filter on AF and remove all INFO except AF from gnomad whole genome dowload
bcftools filter -i "INFO/AF>0.45 && INFO/AF<0.7" wholegenome.vcf.bgz | bcftools annotate -x ^INFO/AF | bgzip > wholegenome_AF_filter.vcf.bgz

# Filter out fails
bcftools view -Oz -f .,PASS wholegenome_AF_filter.vcf.bgz > wholegenome_AF_PASS_filtered.vcf.bgz

# Filter out INDELs/non-SNPs
bcftools view --types snps wholegenome_AF_PASS_filtered.vcf.bgz > wholegenome_AF_PASS_SNP_filtered.vcf

###########################################
#  Now using the genome_bins.bed to make  #
#  a bunch of vcf files for each bin      #
###########################################

temp=$(wc -l <genome_bins.bed)

for i in $(seq 1 $temp); 
do head -$i genome_bins.bed | tail -1 > temp_tabix.bed; 
tabix -h -R temp_tabix.bed wholegenome_AF_PASS_SNP_filtered.vcf.gz > ./temp_folder/${i}_vcf.vcf;
done

# Remove VCF header
for i in *.vcf; do egrep -v "^#" $i > no_header_${i}; done

# Remove empty vcf files 
find . -size 0 -delete

# Pick first line of each file, write to .bed
for i in *.vcf; do line=$(head -n 1 $i); echo $line >> collated_vcf.vcf; done

# convert vcf into a bed
awk '{print $1, $2, $2+1, $3"_"$8}' collated_vcf.vcf >> collated_bed.bed

# Change whitespace to tab
sed -i.bak 's/ /\t/g' collated_bed.bed

# Sort the bed file
sortBed -chrThenSizeA -i collated_bed.bed > collated_sorted_bed.bed

################################################
#  Split each vcf into chucks, pick the        #
#  largest file (i.e. highest number of SNPs)  #
################################################

##### Create bigzip and index files for tabix #####

# Bigzip the vcfs
for i in *vcf; do bgzip -c ${i} > ${i}.gz; done

# Create tabix index
for i in *vcf.gz; do tabix -p vcf ${i}; done

##### Tabix makes ~5,000 regions, pick largest, delete rest #####

# Creat the regions table to pass to tabix
for i in *vcf; do Rscript ~/cgh2bed/SNP_maximiser.R ~/cgh2bed/temp_folder/no_empties/${i} ~/cgh2bed/temp_folder/no_empties/temp_for_the_5000/${i}_output.bed; done

# Tabix a vcf based on the regions created for it
biggest_vcf () {
  for vcf in *vcf;
  do temp=$(wc -l < ./temp_for_the_5000/${vcf}_output.bed);
  for i in $(seq 1 $temp); 
  do head -$i ./temp_for_the_5000/${vcf}_output.bed | tail -1 > temp_tabix.bed; 
  tabix -h -R temp_tabix.bed ${vcf}.gz > ./temp_for_the_5000/${i}_${vcf};
  done;
  cp `find -L temp_for_the_5000 -maxdepth 1 -name \*.vcf -printf '%s %p\n' | sort -nr | head -1 | awk '{print $2}'` ./biggest_vcfs
  rm temp_tabix.bed
  rm ./temp_for_the_5000/*.vcf
  done
  
}

# Paste all the vcfs into one file
for i in *.vcf; do cat ${i} >> collated_vcf.vcf; done

# convert vcf into a bed
awk '{print $1, $2, $2+1, $3"_"$8}' collated_vcf.vcf >> collated_bed.bed

# Change whitespace to tab
sed -i.bak 's/ /\t/g' collated_bed.bed

# Sort the bed file
sortBed -i collated_bed.bed > collated_sorted_bed.bed

# Merge the bedtools regions based on overlap, save the features in the fourth column
bedtools merge -i collated_sorted_bed.bed -d 100 -c 4 -o collapse -delim "|" > collated_sorted_merged_bed.bed
