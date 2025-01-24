# work flow for 0_data
This data was cleaned and assembled in a previous step (by plates) and symlinked into the 0_data contigs dir like so:

(contig list lives in this directory)
```
for i in $(cat contigs.list); do ID=$(echo ${i} | cut -d "," -f 1) CONTIG_DIR=$(echo ${i} | cut -d "," -f 2); ln -s ${CONTIG_DIR} ${PWD}/${ID}.fasta; done
mv contigs.list ../ # dont want contigs in this directory afterwards
```

Using the phyluce processing genomes workflow.
Create directory for "2bit conversion, and create a list subset compatable with phyluce. Its best to not work with more than 30 taxa at a time, so programatically process your files.
```
mkdir 2bit
ls contigs/ | cut -d "." -f 1 > sampleIDs.list
awk 'NR%30==1 {file=substr(FILENAME, 1, length(FILENAME)-5) sprintf("%02d.list", (++count))} {printf "%s ", $0 >> file} END {print ""}' sampleIDs.list
```

Use this script to programatically create a configuration file. It depends on the output of the files you created (e.g., the second input variable for the shell script)
```
#!/bin/bash/

# run like make_2bit_config.sh ${PWD} 17
# where the *number of lists in your working directory*
WORKINGDIR=${1}
END=${2}

for i in $(seq -w 01 ${END}); do
touch taxa${i}.config
printf "[scaffolds]\n" > taxa${i}.config
for TAXA in $(cat sampleIDs${i}.list | tr " " "\n"); do
printf "${TAXA}:${WORKINGDIR}/2bit/${TAXA}/${TAXA}.2bit \n";
done >> taxa${i}.config;
done
```

Now use this bash script to process the sequence data using the flanking and non-flanking data with phyluce
```
#!/bin/bash
#source /local/anaconda3/bin/activate
#conda activate phyluce-1.7.1_py36
#run like
#nohup conda run -n phyluce-1.7.1_py36 bash 2bit_slice.sh ${PWD} 17 6 > 2bit_slice.out &!
# this script creates 2bit files for phyluce to process using the slice from genomes pipline

WORKINGDIR=${1}
NLISTS=${2}
CORES=${3}

cd ${WORKINGDIR}
for i in $(ls contigs); do
    # extract taxa info from contigs
    TAXA=$(echo ${i} | cut -d "." -f 1);
    # bash logic to skip 2bit step if already done
    if [ -f ${WORKINGDIR}/2bit/${TAXA}/${TAXA}.tab ]; then
        continue # skip to next itteration of file
    else
        mkdir -p  2bit/${TAXA}
        printf "convert ${TAXA} to 2bit\n"
        faToTwoBit contigs/${i}  2bit/${TAXA}/${TAXA}.2bit
        twoBitInfo  2bit/${TAXA}/${TAXA}.2bit  2bit/${TAXA}/${TAXA}.tab
        printf "finished converting ${TAXA}\n";
    fi
done

cd ${WORKINGDIR}/2bit
# run lastz
for i in $(seq -w 01 ${NLISTS}); do
    FILE="sampleIDs${i}.list";
    phyluce_probe_run_multiple_lastzs_sqlite \
        --db set${i}.sqlite \
        --output ../Gaedephaga_probes_${i}-lastz \
        --scaffoldlist $(cat ../${FILE}) \
        --genome-base-path ./ \
        --probefile ../Gaedephaga_probes.fasta \
        --coverage 50 \
        --identity 80 \
        --cores ${CORES}
done

cd ${WORKINGDIR}
# extract core probe regions
for i in $(seq -w 01 ${NLISTS}); do
    phyluce_probe_slice_sequence_from_genomes \
        --lastz Gaedephaga_probes_${i}-lastz \
        --conf taxa${i}.config \
        --flank 0 \
        --name-pattern "Gaedephaga_probes.fasta_v_{}.lastz.clean" \
        --output Gaedephaga_probes_${i}-noflank-fasta
done

# extract core + flanking of 100bp reads shouldn't be longer
# than aroud 320 bp (100 + 120 + 100; Flanking + Core + Flanking)
for i in $(seq -w 01 ${NLISTS}); do
    phyluce_probe_slice_sequence_from_genomes \
        --lastz Gaedephaga_probes_${i}-lastz \
        --conf taxa${i}.config \
        --flank 100 \
        --name-pattern "Gaedephaga_probes.fasta_v_{}.lastz.clean" \
        --output Gaedephaga_probes_${i}-flank100-fasta
done
```

Phyluce changes the names of the files. So we need to get all data into the appropriate direcotry while changing the names back to what we named them
```
for CURATION in noflank flank50; do
mkdir /home/cody/Calosoma_phylo/Calosoma_2025/2_find_uces/contigs/${CURATION}
    for i in $(ls Gaedephaga_probes_*-${CURATION}-fasta/*); do
        PATHWAY=$(echo ${i} | cut -d "/" -f 1);
        OLD=$(echo ${i} | cut -d "/" -f 2 | cut -d "." -f 1);
        NEW=$(echo ${i} | cut -d "/" -f 2 | cut -d "." -f 1 | tr '[:lower:]' '[:upper:]');
        ln -s /home/cody/Calosoma_phylo/Calosoma_2025/0_data/${PATHWAY}/${OLD}.fasta /home/cody/Calosoma_phylo/Calosoma_2025/2_find_uces/contigs/${CURATION}/${NEW}.fasta;
    done;
done
```

# workflow for 1_mitofinder

Mitochondrial bycatch and mitochondrial genome construction using 
MitoFinder by RemiAllio: https://github.com/RemiAllio/MitoFinder#run-mitofinder-with-singularity
Singularity: https://sylabs.io/guides/3.5/user-guide/environment_and_metadata.html?highlight=conda


need to create a directory of the contigs from 0_data/contigs

Get reference database for your taxa
```
conda activate sratoolkit
esearch -db nuccore -query "\"mitochondrion\"[All Fields] AND (Carabinae[Organism]) AND (refseq[filter] AND mitochondrion[filter] AND (\"10000\"[SLEN] : \"20000\"[SLEN]))" | efetch -format gbwithparts > reference.gb
```

cant get mitofinder to work like before... :'); so I'm cheesing it and using interactive shell with putty on the musuem computer and letting that go. 
```
singularity shell -H ${PWD} mitofinder_v1.4.2.sif
for CONTIGS in $(ls assemblies/*.fasta); do TAXA=$(echo ${CONTIGS} | cut -d "/" -f 2 | cut -d "." -f 1);  mitofinder -j ${TAXA} -a ${CONTIGS} -r Harpalinae.gb -m 64 -p 4 -o 5 --rename-contig no ; done
```


.... additional work to be done, ID the COI data, isolate the circular genomes
```
grep "Evidences of circularization were found!" *.log > ../circularization.txt
```
CBX0250_MitoFinder.log:Evidences of circularization were found!
CBX0252_MitoFinder.log:Evidences of circularization were found!
CBX0269_MitoFinder.log:Evidences of circularization were found!
CBX0302_MitoFinder.log:Evidences of circularization were found!
CBX0322_MitoFinder.log:Evidences of circularization were found!
CBX0406_MitoFinder.log:Evidences of circularization were found!
CBX0591_MitoFinder.log:Evidences of circularization were found!
CBX1086_MitoFinder.log:Evidences of circularization were found!

One of these is a Calisthenes 
# workflow for 2_find_uces

contigs were recovered from phyluce slice from genomes step and symlinked to this directory need to first make a config file and decide what taxa to keep. This requires the program seqkit which can be installed with conda.
https://bioinf.shenwei.me/seqkit/

```
# get remove list
cd contigs
for CURATION in flank50 noflank; do
 seqkit stats -T ${CURATION}/*.fasta > ${CURATION}_stats.tsv
done
awk 'FNR>=1 && $4 <= 30 {print $1}' *_stats.tsv > remove_samples.list

# remove those taxa with less than 20/30 sequences
for i in $(cat remove_samples.list); do rm ${i}; done

cd ../
for CURATION in flank50 noflank; do
    printf "[all]\n" > ${CURATION}_taxon-set.conf
    ls ${PWD}/contigs/${CURATION} | cut -d "." -f 1 >> ${CURATION}_taxon-set.conf;
done
```

Now that we have cleaned up the data a little, we need to begin processing the data. This is easy to do quickly.
you input your working dir, and provide the name of the input data set that is defined in the config file created previously
```
#!/bin/bash
# need to be in the 2_find_uces directory!
# nohup conda run -n phyluce-1.7.1_py36 bash get_incomplete_fasta.sh ${PWD} all > get_incomplete_fasta.out &
# here you input your working dir, and provide the name of the input data set that is defined in the config file created previously
WORKINGDIR=${1}
DATASET=${2}
for CURATION in flank50 noflank; do
cd ${WORKINGDIR}
    phyluce_assembly_match_contigs_to_probes \
        --contigs contigs/${CURATION} \
        --probes ../0_data/Gaedephaga_probes.fasta \
        --min-coverage 50 \
        --min-identity 80 \
        --output Gaedephaga_uce-loci_${CURATION}_search-results;

    mkdir -p ../3_taxon-sets/Gaedephaga_uce-loci_${CURATION}/${DATASET}/log;

    phyluce_assembly_get_match_counts \
        --locus-db Gaedephaga_uce-loci_${CURATION}_search-results/probe.matches.sqlite \
        --taxon-list-config ${CURATION}_taxon-set.conf \
        --taxon-group ${DATASET} \
        --incomplete-matrix \
        --output ../3_taxon-sets/Gaedephaga_uce-loci_${CURATION}/${DATASET}all-taxa-incomplete.conf;

    cd ../3_taxon-sets/Gaedephaga_uce-loci_${CURATION}/${DATASET}
    #phyluce_assembly_get_fastas_from_match_counts --contigs /home/cody/Calosoma_phylo/Calosoma_WORKING/Gaedephaga_uce-loci_flank50 --locus-db /home/cody/Calosoma_phylo/Calosoma_WORKING/Gaedephaga_uce-loci_flank50_search-results/probe.matches.sqlite --match-count-output all-taxa-incomplete.conf --output all-taxa-incomplete.fasta --incomplete-matrix all-taxa-incomplete.incomplete --log-path log;
# always be very explicit about your paths since we "jump" around from directories
    phyluce_assembly_get_fastas_from_match_counts \
        --contigs ${WORKINGDIR}/contigs/${CURATION}  \
        --locus-db ${WORKINGDIR}/Gaedephaga_uce-loci_${CURATION}_search-results/probe.matches.sqlite \
        --match-count-output all-taxa-incomplete.conf \
        --output all-taxa-incomplete.fasta \
        --incomplete-matrix all-taxa-incomplete.incomplete \
        --log-path log;

# ALWAYS breaks here... hopefully not this time...
    phyluce_assembly_explode_get_fastas_file \
        --input all-taxa-incomplete.fasta \
        --output taxa-exploded-fastas \
        --by-taxon

    phyluce_assembly_explode_get_fastas_file \
        --input all-taxa-incomplete.fasta \
        --output loci-exploded-fastas
done
```
# workflow for 3_taxon-sets

Align and trim data using an awk script to linearize the fasta files & clean up the headers the fast files are aligned using mafft auto, and trimmed with trimal auto.

```
#!/bin/bash

# once conda env is set up on pyrgus, run like:
# nohup conda run -n phyluce-1.7.1_py36 bash align_trim.sh ${PWD} all 4 94 > align_trim.out &!
#
# otherwise:
# nohup bash align_trim.sh ${PWD} 4 94 > align_trim.out &!
source /local/anaconda3/bin/activate
conda activate seqkit

WORKINGDIR=${1}
TAXONSET=${2} # taxonset directory (e.g., all from the config file)
CORES=${3} # number of threads/cores to use for analyses
LOCUSOCCUPANCY=${4} # number of taxa that must be present in a locus

# loop through directories
for CURATION in flank50 noflank; do
    # change to directories
    cd ${WORKINGDIR}/Gaedephaga_uce-loci_${CURATION}/${TAXONSET}
    
    # create list with at loci containing at least 40% taxon occupancy (165/413 = ~40%)
    seqkit stats -T loci-exploded-fastas/*.fasta | \
        awk -v LOCUSOCCUPANCY=${LOCUSOCCUPANCY} '$4 >= LOCUSOCCUPANCY {print}' > 40p.list

    # make an unaligned directory
    mkdir ${WORKINGDIR}/Gaedephaga_uce-loci_${CURATION}/${TAXONSET}/unaligned

    # linearize the fasta files and clean up headers with awk based on the 40p.list file,
    for FASTA in $(awk 'NR>1 {print $1}' 40p.list); do
        UCE_ID=$(echo ${FASTA} | cut -d "/" -f 2 | cut -d "." -f 1); # create a variable from the for loop
        seqkit seq -w 0 loci-exploded-fastas/${UCE_ID}.unaligned.fasta | \
        awk 'BEGIN {RS=">"; FS= "\n"} NEXT NR>1 {split($0,header,"|"); split(header[1],code,"_"); print ">" code[2] "\n" $2 }' \
            > ${WORKINGDIR}/Gaedephaga_uce-loci_${CURATION}/${TAXONSET}/unaligned/${UCE_ID}.fasta;
    done

    # deactivate teh seqkit conda environment
    conda deactivate
    conda activate phyluce-1.7.1_py36

    # using the mafft and trimal that comes with phyluce align and trim the fasta files
    mkdir ${WORKINGDIR}/Gaedephaga_uce-loci_${CURATION}/${TAXONSET}l/mafft
    for FASTA in $(ls unaligned); do
        UCE_ID=$(echo ${FASTA} | cut -d "." -f 1); # create a variable from the for loop
        mafft --auto --thread ${CORES} \
        ${WORKINGDIR}/Gaedephaga_uce-loci_${CURATION}/${TAXONSET}/unaligned/${UCE_ID}.fasta \
            > ${WORKINGDIR}/Gaedephaga_uce-loci_${CURATION}/${TAXONSET}/mafft/${UCE_ID}.aligned.fasta;
    done
    mkdir ${WORKINGDIR}/Gaedephaga_uce-loci_${CURATION}/${TAXONSET}/mafft_trimal
    for FASTA in $(ls mafft); do
        UCE_ID=$(echo ${FASTA} | cut -d "." -f 1); # create a variable from the for loop
        trimal -automated1 \
            -in ${WORKINGDIR}/Gaedephaga_uce-loci_${CURATION}/${TAXONSET}/mafft/${UCE_ID}.aligned.fasta \
            -out ${WORKINGDIR}/Gaedephaga_uce-loci_${CURATION}/${TAXONSET}/mafft_trimal/${UCE_ID}.trimmed.fasta;
    done
done
```
NEED TO RERUN TRIMAL, did not properlly call the aligned file
