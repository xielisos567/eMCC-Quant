#step 1: Reads mapping with the eMCCDB database:
##Example: the 3-HP_4-HP pathway
for i in `less ../metagenome.txt`
do
        bbmap.sh in=../metagenome/${i}_1.fq.gz in2=../metagenome/${i}_2.fq.gz rpkm=./${i}_3-HP_4-HP.rpkm covstats=./${i}_3-HP_4-HP.cov idtag=t  ambiguous=best  minid=0.95 -Xmx50g usemodulo ref=/mnt/hpc/home/jiashulei/download/antarctic/database_pathways_nucl/3-HP_4-HP_all.fna nodisk
done

##Example: the CBB pathway
for i in `less ../metagenome.txt`
do
        bbmap.sh in=../metagenome/${i}_1.fq.gz in2=../metagenome/${i}_2.fq.gz rpkm=./${i}_CBB.rpkm covstats=./${i}_CBB.cov idtag=t  ambiguous=best  minid=0.95 -Xmx50g usemodulo ref=/mnt/hpc/home/jiashulei/download/antarctic/database_pathways_nucl/CBB_all.fna nodisk
done


#step 2: data processing
for i in `ls ./`
do 
	sed -i '1, 5d' ./$i && less ./$i | awk -v var="$i" '{print var"\t"$0}' >> ../combine_results_all.out
done

##Filter rows in column 8:
awk -F "\t" '{if($8!=0) print}'  combine_3-HP_4-HP_all.txt > 3-HP_4-HP.txt
awk -F "\t" '{if($8!=0) print}'  combine_CBB_all.txt > CBB.txt

#step 3: gene && pathway abundance (TPM)
python eMCCDB.py

