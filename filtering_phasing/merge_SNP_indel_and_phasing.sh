# Commands used to filter out variants with certain FILTER values out of vcf, concat SNP and INDEL file of each sample to obtain one file for each sample, and merge them all togheter to obtain one multisample
# file to feed to SHAPEIT5 for phasing. Two directories are required containing SNP and INDEL files separately, and a txt file containing the suffinx of the files (sample names)

# select variants with PASS filter, do this for both SNP and INDEL files. (to execute from the directory containing all the vcf files. "indel_filtered" is the destination directory for filtered indel files and
# complete_samples is the file containing all sample names/file suffix)

# for indels
while IFS="" read -r p || [ -n "$p" ]; do ls | grep "$p" | xargs -I '{}' bcftools view -i 'FILTER="PASS" || FILTER="VQSRTrancheINDEL99.90to99.95" || FILTER="VQSRTrancheINDEL99.95to100.00"' '{}' > ../indel_filtered/"$p".indel.vcf ; done < ../complete_samples.txt
#for snp
while IFS="" read -r p || [ -n "$p" ]; do ls | grep "$p" | xargs -I '{}' bcftools view -i 'FILTER="PASS" || FILTER="VQSRTrancheSNP99.90to99.95" || FILTER="VQSRTrancheSNP99.95to100.00"' '{}' > ../snp_filtered/"$p".snp.vcf ; done < ../complete_samples.txt

# sort, compress, and remove uncompressed files(to be executed inside the directory containing the files)
for f in *; do bcftools sort $f -Oz -o "$f".gz; done && rm *vcf

# create indexes (to be executed inside the directory containing the files)
for f in *; do bcftools index -t $f; done


# combine SNP and INDEL files togheter for each patient (takes files from the respective directories, conbine them and saves them in /combned_vcf)
 while IFS="" read -r p || [ -n "$p" ]; do bcftools concat -a snp_filtered/"${p}".snp.vcf.gz indel_filtered/"${p}".indel.vcf.gz -Oz -o combined_vcf/"$p".vcf.gz ; done < complete_samples.txt

