# nanopore_r10
This readme file will run through the codes and snippetes needed to run nanopore whole genome analysis 

Needed_softwares
  artic
  clair3

Creation of primer trimmed and trimmed bam through artic pipeline

```bash script
output_dir=/hpcigib/varsha.r/Dengue_wholegenome/2023_runs/gridion/lib36/FASTQ/artic_res
input_bed=/hpcigib/varsha.r/tools/artic-ncov2019/primer_schemes/denv2/V1
#v_caller=/hpcigib/varsha.r/tools/Clair3

for i in $(ls *.fastq | rev | cut -c 7- | rev | uniq);

do

##alinging_indexing
minimap2 -a -x map-ont -t 40 $input_bed/denv2.reference.fasta ${i}.fastq | samtools view -bS -F 4 - | samtools sort -o $output_dir/${i}.sorted.bam
samtools index $output_dir/${i}.sorted.bam

##trimming_primers
align_trim --normalise 200 $input_bed/denv2.scheme.bed --start --remove-incorrect-pairs --report $output_dir/${i}.alignreport.txt < $output_dir/${i}.sorted.bam 2> $output_dir/${i}.alignreport.er | samtools sort -T $output_dir/${i} -o $output_dir/${i}.trimmed.rg.sorted.bam

align_trim --normalise 200 $input_bed/denv2.scheme.bed --remove-incorrect-pairs --report $output_dir/${i}.alignreport.txt < $output_dir/${i}.sorted.bam 2> $output_dir/${i}.alignreport.er | samtools sort -T $output_dir/${i} -o $output_dir/${i}.primertrimmed.rg.sorted.bam

##indexing_the_primers_trimmed_bams
samtools index $output_dir/${i}.trimmed.rg.sorted.bam
samtools index $output_dir/${i}.primertrimmed.rg.sorted.bam

##coverage_mask_file
artic_make_depth_mask --store-rg-depths $input_bed/denv2.reference.fasta $output_dir/${i}.primertrimmed.rg.sorted.bam $output_dir/${i}.coverage_mask.txt


done
```
Variant calling through Clair3

```bash script
output_dir=/hpcigib/varsha.r/Dengue_wholegenome/2023_runs/gridion/lib36/FASTQ/artic_res/trimmed_bam/VCF
input_bed=/hpcigib/varsha.r/tools/artic-ncov2019/primer_schemes/denv2/V1
v_caller=/hpcigib/varsha.r/tools/Clair3

for i in $(ls *.trimmed.rg.sorted.bam | rev | cut -c 23- | rev | uniq); 

do

$v_caller/run_clair3.sh --bam_fn ${i}.trimmed.rg.sorted.bam --ref_fn $input_bed/denv2.reference.fasta --threads=40 --platform "ont" --model_path $v_caller/models/r1041_e82_400bps_sup_v420 --output $output_dir/${i} --no_phasing_for_fa --include_all_ctgs --haploid_precise --var_pct_full=0.8

done
```

Renaming the VCF with the folder name
```bash script
##renaming the filename with folder name
find . -iname '*.vcf.gz' | while read fn; do name=$(basename "$fn") ; dir=$(dirname "$fn") ; mv "$fn" "$dir/$(basename "$dir")-$name" ;done
```

Post processing of VCF to remove low depth and low genotype quality

```bash script
for i in $(ls *-merge_output.vcf | rev | cut -c 18- | rev | uniq); do bcftools view --include 'MIN(FMT/DP)<=20 | MIN(FMT/GQ)<=3' ${i}-merge_output.vcf > ${i}_fail.vcf; done;
for i in $(ls *-merge_output.vcf | rev | cut -c 18- | rev | uniq); do bcftools view --include 'MIN(FMT/DP)>20 & MIN(FMT/GQ)>3' ${i}-merge_output.vcf > ${i}_pass.vcf; done;

```
Creation of preconsensus fata with the low covered region as N

```bash script
input_bed=/hpcigib/varsha.r/tools/artic-ncov2019/primer_schemes/denv2/V1
coverage_mask=/hpcigib/varsha.r/Dengue_wholegenome/2023_runs/gridion/lib34/FASTQ/artic_res/coverage_mask/
fail_vcf=/hpcigib/varsha.r/Dengue_wholegenome/2023_runs/gridion/lib34/FASTQ/artic_res/VCF/merged_VCF/fail_vcf

# Iterate over files in the coverage_mask directory
# Iterate over files in the coverage_mask directory
for coverage_mask_file in $coverage_mask/*; do
    if [ -f $coverage_mask_file ]; then
        # Extract base name until the first dot
        coverage_mask_base=$(basename $coverage_mask_file | cut -d '.' -f 1)

        # Iterate over files in the fail_vcf directory
        for fail_vcf_file in $fail_vcf/*; do
            if [ -f $fail_vcf_file ]; then
                # Extract base name from the last underscore to the end
                fail_vcf_base=$(basename $fail_vcf_file | cut -d '.' -f 1)

                # Check if base names match
                if [ $coverage_mask_base = $fail_vcf_base ]; then
                    # Run the command
                    artic_mask $input_bed/denv2.reference.fasta $coverage_mask_file $fail_vcf_file ${coverage_mask_base}.preconsensus.fasta
                fi
            fi
        done
    fi
done
```
Creation of consensus fasta
```bash script
coverage_mask=/hpcigib/varsha.r/Dengue_wholegenome/2023_runs/gridion/lib34/FASTQ/artic_res/coverage_mask/
pass_vcf=/hpcigib/varsha.r/Dengue_wholegenome/2023_runs/gridion/lib34/FASTQ/artic_res/VCF/merged_VCF/pass_vcf
preconsensus=/hpcigib/varsha.r/Dengue_wholegenome/2023_runs/gridion/lib34/FASTQ/artic_res/preconsensus

# Iterate over files in the coverage_mask directory
for coverage_mask_file in $coverage_mask/*; do
    if [ -f $coverage_mask_file ]; then
        # Extract base name until the first dot
        coverage_mask_base=$(basename $coverage_mask_file | cut -d '.' -f 1)
        
        # Iterate over files in the pass_vcf directory
        for pass_vcf_file in $pass_vcf/*pass.vcf.gz; do
            if [ -f $pass_vcf_file ]; then
                # Extract base name from the last underscore to the end
                pass_vcf_base=$(basename $pass_vcf_file | cut -d '.' -f 1)
                
                # Check if base names match
                if [ $coverage_mask_base = $pass_vcf_base ]; then
                    # Iterate over files in the preconsensus directory
                    for preconsensus_file in $preconsensus/*; do
                        if [ -f $preconsensus_file ]; then
                            # Extract base name from the last underscore to the end
                            preconsensus_base=$(basename $preconsensus_file | cut -d '.' -f 1)
                            
                            # Check if base names match
                            if [ $coverage_mask_base = $preconsensus_base ]; then
                                # Run the command
                                bcftools consensus -f $preconsensus_file $pass_vcf_file -m $coverage_mask_file -o ${coverage_mask_base}.consensus.fasta
                            fi
                        fi
                    done
                fi
            fi
        done
    fi
done
```
Sequencing metadata sheet preparation

```bash script
##DEPTH
for j in $(ls *.bam);do echo -n $j;samtools depth $j | awk '{sum+=$3} END { print " ",sum/NR}' ; done;

##READ
for i in $(ls *.fastq | rev | cut -c 7- | rev | uniq); do wc -l ${i}.fastq; done;

##COVERAGE THROUGH BAM
for i in *.bam; do $sam/samtools coverage $i; done;

##MAPPING RATE
for i in $(ls *.bam | rev | cut -c 5- | rev | uniq); do samtools flagstat ${i}.bam > ${i}_mappingrate.txt; done;
**TAKING SPECIFIC LINES FROM FILE AND WRITE IN ANOTHER
for file in *mappingrate.txt; do filename=$(basename "$file"); fifth_line=$(sed -n '5p' "$file");echo "$filename: $fifth_line" >> mappingrate.txt;done

#COVERAGE THROUGH FASTA
#for changing the name of file name to fasta header name
for file in *.fasta; do awk -v name="${file%%.*}" '/^>/{print ">" name; next} 1' "$file" > tmp && mv tmp "$file";done
cat *.fasta > run5.fasta
seqtk comp -u run5.fasta > run5_coverage.txt
```


