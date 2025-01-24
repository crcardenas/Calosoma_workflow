Data is Downloaded into RAW directory

QC_assembly.sh script
```
#!/bin/bash

###########################################################################
# QC and assemble the data
# 
# have you made fastpQC directory?
# have you made trimmed directory?
# have you made spades directory?
# 
# to actually run it do:
# nohup conda run -n anch-hyb-UCE bash QC_assembly.sh > QC_assembly.out &
###########################################################################

cd /home/cody/Calosoma_phylo/Calosoma_Lib5_Gae1/0_data

# fastp versoin fastp 0.19.5
for SampleID in $(cat ../assembly_taxa.list | cut -f 1); do
    fastp -i ./raw/${SampleID}_*_R1*.fastq.gz \
       -o ./trimmed/${SampleID}_r1.fastq.gz \
       -I ./raw/${SampleID}_*_R2*.fastq.gz \
       -O ./trimmed/${SampleID}_r2.fastq.gz \
       --thread 6 -h ./fastpQC/${SampleID}_QC.html;
done

# now do assembly
# spades version 3.15.5
for SampleID in $(cat ../assembly_taxa.list | cut -f 1); do
    spades.py --sc --careful -t 6 -m 20000 \
        -1 ./trimmed/${SampleID}_r1.fastq.gz \
        -2 ./trimmed/${SampleID}_r2.fastq.gz \
        -o ./spades/${SampleID};
done
```
