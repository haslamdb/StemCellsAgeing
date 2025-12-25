### Pathway abundance ####

#cd ~/Documents/Alignments/Humann3Alignments/PathwayAbundance/

  #mkdir GeigerFiles

#Files=$(cat ~/Desktop/GeigerSamples.txt)

#for f in ${Files}; do

#cp ${f}_pathabundance.tsv GeigerFiles/


done
# Now normalize pathway to cpm

cd ~/Documents/Alignments/Humann3Alignments/PathwayAbundance

samples=$(cat ~/Desktop/GeigerSamples.txt)

for s in ${samples}; do

humann2_renorm_table --input GeigerFiles/${s}_pathabundance.tsv --output GeigerFiles/${s}_pathabundance-cpm.tsv --units cpm --update-snames

done

# Normalize

#mkdir GeigerFiles/Normalized

mv GeigerFiles/*_pathabundance-cpm.tsv GeigerFiles/Normalized

humann2_join_tables -i GeigerFiles/Normalized/ -o GeigerFiles/GeigerFiles_pathabundance20220910-cpm.tsv

humann2_split_stratified_table -i GeigerFiles/GeigerFiles_pathabundance20220910-cpm.tsv -o GeigerFiles_pathabundance_Unstratified_20220910.tsv

humann2_rename_table --input Unstratified_20220910/GeigerFiles_pathabundance20220910-cpm_unstratified.tsv --output Geiger_genefamilies-uniref90-names.tsv --names uniref90

# Normalize normalized table - not necessary

#humann2_renorm_table -i GeigerFiles/GeigerFiles_pathabundance20200103-cpm.tsv  -o GeigerFiles/GeigerFiles_pathabundance20200103_relab-cpm.tsv --units relab

##### Genefamilies #####

cd ~/Documents/Alignments/Humann3Alignments/GeneFamilies/

  mkdir GeigerFiles

Files=$(cat ~/Desktop/GeigerSamples.txt)

for f in ${Files}; do

cp ${f}_genefamilies.tsv GeigerFiles/

  done

# Now normalize genefamilies to cpm

cd ~/Documents/Alignments/Humann3Alignments/GeneFamilies

samples=$(cat ~/Desktop/GeigerSamples.txt)

for s in ${samples}; do

humann2_renorm_table --input GeigerFiles/${s}_genefamilies.tsv --output GeigerFiles/${s}_genefamilies-cpm.tsv --units cpm --update-snames

done

# Normalize

mkdir GeigerFiles/Normalized

mv GeigerFiles/*_genefamilies-cpm.tsv GeigerFiles/Normalized

humann2_join_tables -i GeigerFiles/Normalized/ -o GeigerFiles/GeigerFiles_genefamilies20211018-cpm.tsv

humann2_split_stratified_table -i GeigerFiles_genefamilies20211018-cpm.tsv -o GeigerFiles_genefamilies_Unstratified20211018.tsv

