Taxonomy and cuntional pathways for kraken2 sequencing data was assigned using HUMAnN3 pipelines

##Species-level functional profiling of metagenomes and metatranscriptomes
export PATH=/xtdisk/jiapl_group/zhuxl/conda/bin:$PATH
export PATH=/pnas/pmod/yuann/software/bowtie2-2.2.3:$PATH
export PATH=/usr/bin:$PATH

/xtdisk/jiapl_group/zhuxl/conda/bin/humann3
/pnas/pmod/yuann/software/bowtie2-2.2.3/bowtie2
/usr/bin/perl

export PATH=/software/biosoft/htop/bin:$PATH
which htop

pip install humann
/xtdisk/jiapl_group/zhuxl/conda/bin/humann_test # everything is ok.
/xtdisk/jiapl_group/zhuxl/conda/bin/metaphlan --version 

MetaCyc pathway functional profiling of stool metagenomes
Franzosa, E. A. et al. Species-level functional profiling of metagenomes and metatranscriptomes. Nat. Methods 15, 962–968 (2018).
Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3

https://huttenhower.sph.harvard.edu/humann2
https://github.com/biobakery/biobakery/wiki/humann2
https://github.com/biobakery/biobakery/wiki/humann3
/p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann/humann/data
/p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann
http://huttenhower.sph.harvard.edu/humann_data/uniprot/
https://pypi.python.org/pypi/humann2

############################

/xtdisk/jiapl_group/zhuxl/conda/bin/humann_databases --download chocophlan full /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann --update-config yes
#To download the full UniRef90 database (20.7GB, recommended):
/xtdisk/jiapl_group/zhuxl/conda/bin/humann_databases --download uniref uniref90_diamond /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann --update-config yes
#To download the EC-filtered UniRef90 database (0.9GB): 
/xtdisk/jiapl_group/zhuxl/conda/bin/humann_databases --download uniref uniref90_ec_filtered_diamond /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann --update-config yes
#To download the full UniRef50 database (6.9GB):
/xtdisk/jiapl_group/zhuxl/conda/bin/humann_databases --download uniref uniref50_diamond /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann --update-config yes
#To download the EC-filtered UniRef50 database (0.3GB):
/xtdisk/jiapl_group/zhuxl/conda/bin/humann_databases --download uniref uniref50_ec_filtered_diamond /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann --update-config yes


mkdir chocophlan_v296_201901
mkdir uniref90_v201901
mkdir mapping_v201901

tar -zxvf full_chocophlan.v296_201901.tar.gz -C ./chocophlan_v296_201901/
tar -zxvf uniref90_annotated_v201901.tar.gz -C uniref90_v201901
tar -zxvf full_mapping_v201901.tar.gz -C ./mapping_v201901/

/xtdisk/jiapl_group/zhuxl/conda/bin/humann_config --update database_folders nucleotide /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann/database_chocophlan
/xtdisk/jiapl_group/zhuxl/conda/bin/humann_config --update database_folders protein /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann/database_uniref
/xtdisk/jiapl_group/zhuxl/conda/bin/humann_config --update database_folders utility_mapping /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann/database_mapping

/xtdisk/jiapl_group/zhuxl/conda/bin/humann_config --print
/xtdisk/jiapl_group/zhuxl/conda/bin/humann_databases --available


/xtdisk/jiapl_group/zhuxl/conda/bin/humann3 --input /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann/examples/demo.fasta.gz --output /p300s/jiapl_group/zhuxl/Microbiome/humann --nucleotide-database /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann/database_chocophlan --protein-database /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann/database_uniref

ls ../Result_classfied > list_993
sed -i s/.classified.fasta//g list_993
/xtdisk/jiapl_group/zhuxl/conda/bin/humann3 --input /p300s/jiapl_group/zhuxl/Microbiome/Result_classfied/10125714.classified.fasta --nucleotide-database /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann/database_chocophlan --protein-database /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann/database_uniref --output /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/10125714 --threads 10
awk -F '\t' '{print "mkdir /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/"$1" && /xtdisk/jiapl_group/zhuxl/conda/bin/humann3 --input /p300s/jiapl_group/zhuxl/Microbiome/Result_classfied/"$1".classified.fasta --nucleotide-database /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann/database_chocophlan --protein-database /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann/database_uniref --output /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/"$1" --threads 10 1>/p300s/jiapl_group/zhuxl/Microbiome/humann/Result/"$1"/log 2>/p300s/jiapl_group/zhuxl/Microbiome/humann/Result/"$1"/err && touch /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/"$1".SUCCESS"}' list_993 >work.sh

############################


#####################
/xtdisk/jiapl_group/zhuxl/conda/bin/python3.9 /xtdisk/jiapl_group/zhuxl/conda/bin/humann3 --input /p300s/jiapl_group/zhuxl/Microbiome/Result_classfied/26463284.classified.fasta --nucleotide-database /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann/database_chocophlan --protein-database /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann/database_uniref --output /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/26463284 --threads 10
/xtdisk/jiapl_group/zhuxl/conda/bin/python /xtdisk/jiapl_group/zhuxl/conda/bin/metaphlan /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/26463284/26463284.classified_humann_temp/tmpn5jd4nfb/tmp935u0tvy -t rel_ab -o /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/26463284/26463284.classified_humann_temp/26463284.classified_metaphlan_bugs_list.tsv --input_type fasta --bowtie2out /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/26463284/26463284.classified_humann_temp/26463284.classified_metaphlan_bowtie2.t
/xtdisk/jiapl_group/zhuxl/conda/bin/python /xtdisk/jiapl_group/zhuxl/conda/bin/read_fastx.py -l 70 /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/26463284/26463284.classified_humann_temp/tmpn5jd4nfb/tmp935u0tvy
perl /pnas/pmod/yuann/software/bowtie2-2.2.3/bowtie2 --seed 1992 --quiet --no-unal --very-sensitive -S - -x /xtdisk/jiapl_group/zhuxl/conda/lib/python3.9/site-packages/metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103 -p 10 -U - -f
/pnas/pmod/yuann/software/bowtie2-2.2.3/bowtie2-align-l --wrapper basic-0 --seed 1992 --quiet --very-sensitive -x /xtdisk/jiapl_group/zhuxl/conda/lib/python3.9/site-packages/metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103 -p 10 -f --passthrough -U -

/xtdisk/jiapl_group/zhuxl/conda/bin/diamond blastx --query /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/18099492/18099492.classified_humann_temp/tmp4izje6r5/tmpa5gl_rte --evalue 1.0 --threads 10 --top 1 --outfmt 6 --db /p300s/jiapl_group/zhuxl/Microbiome/Pipeline/humann/database_uniref/uniref90_201901b_full --out /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/18099492/18099492.classified_humann_temp/tmp4izje6r5/diamond_m8_gs0ezmot --tmpdir /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/18099492/18099492.classified_humann_temp/tmp4izje6r5
/xtdisk/jiapl_group/zhuxl/conda/bin/python /xtdisk/jiapl_group/zhuxl/conda/bin/metaphlan /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/25383485/25383485.classified_humann_temp/tmpghzkebxr/tmp8slc0tyt -t rel_ab -o /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/25383485/25383485.classified_humann_temp/25383485.classified_metaphlan_bugs_list.tsv --input_type fasta --bowtie2out /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/25383485/25383485.classified_humann_temp/25383485.classified_metaphlan_bowtie2.txt


cd /p300s/jiapl_group/zhuxl/Microbiome/humann
/xtdisk/jiapl_group/zhuxl/conda/bin/humann_join_tables --input /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/ --search-subdirectories --file_name genefamilies --output Result_genefamilies.tsv
/xtdisk/jiapl_group/zhuxl/conda/bin/humann_join_tables --input /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/ --search-subdirectories --file_name pathabundance --output Result_pathabundance.tsv
/xtdisk/jiapl_group/zhuxl/conda/bin/humann_join_tables --input /p300s/jiapl_group/zhuxl/Microbiome/humann/Result/ --search-subdirectories --file_name pathcoverage --output Result_pathcoverage.tsv

/xtdisk/jiapl_group/zhuxl/conda/bin/humann_barplot --input $TABLE.tsv --focal-feature $FEATURE --outfile $FIGURE
/xtdisk/jiapl_group/zhuxl/conda/bin/humann_regroup_table --input $TABLE --groups $GROUPS --output $TABLE2
/xtdisk/jiapl_group/zhuxl/conda/bin/humann_rename_table --input $TABLE --names $NAMES --output $TABLE2 
/xtdisk/jiapl_group/zhuxl/conda/bin/humann_renorm_table --input $TABLE --units $CHOICE --output $TABLE2
/xtdisk/jiapl_group/zhuxl/conda/bin/humann_split_stratified_table --input $TABLE --output $OUTPUT_DIR
/xtdisk/jiapl_group/zhuxl/conda/bin/humann_split_table --input $TABLE --output $OUTPUT_DIR
/xtdisk/jiapl_group/zhuxl/conda/bin/humann_unpack_pathways --input-genes Result_genefamilies.tsv --input-pathways Result_pathabundance.tsv -o Result_gene_pathway


##################



head -1 Result_genefamilies.tsv | sed "s/\# //g" | sed s/.classified_Abundance-RPKs//g >Result_genefamilies.tsv.head
grep unclassified Result_genefamilies.tsv | sed s/\|unclassified//g > Result_genefamilies.tsv.res
cat Result_genefamilies.tsv.head Result_genefamilies.tsv.res >Result_genefamilies.tsv.result
sed -i s/UniRef90_//g Result_genefamilies.tsv.result

cut -f 1 Result_genefamilies.tsv.result |sed '1d' | sed "s/UniRef90_//g" >geneID
https://www.uniprot.org/tool-dashboard (ID mapping)
source database: UniProtKB AC/ID
Target database: UniProtKB

head -1 Result_pathabundance.tsv | sed "s/\# //g" | sed s/.classified_Abundance//g >Result_pathabundance.tsv.head
grep unclassified Result_pathabundance.tsv | grep PWY |sed s/\|unclassified//g > Result_pathabundance.tsv.res
cat Result_pathabundance.tsv.head Result_pathabundance.tsv.res >Result_pathabundance.tsv.result


head -1 Result_pathcoverage.tsv | sed "s/\# //g" | sed s/.classified_Coverage//g >Result_pathcoverage.tsv.head
grep unclassified Result_pathcoverage.tsv | grep PWY | sed s/\|unclassified//g > Result_pathcoverage.tsv.res
cat Result_pathcoverage.tsv.head Result_pathcoverage.tsv.res >Result_pathcoverage.tsv.result

