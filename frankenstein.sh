#this gives us the data only to read into a pandas dataframe
cat France.vcf | \awk -F$'\t' '{
        if (NR<34){print$0}    
        if (NR==34){print$q}
        if (NR>34) {
            for (i=1; i<10; i++){printf $i "\t"};
            for (i>9; i<NF; i++){
                if  ($i =="1/1") {printf $i"\t"};
                if  ($i=="1/0" ) {printf $i"\t"};
                if  ($i=="0/1" ) {printf $i"\t"};
                if  ($i=="0/0" ) {printf $i"\t"};
                };
            {printf "\n"}
            }
        }'  > frankensteins.vcf

#to go from the full frankenstein genomes vcf to just one per population
cat frankensteins.vcf| \
    awk -F$'\t' '{
        if (NR<35){print$0}

        if (NR>34) {
            for (i=1; i<10; i++){printf $i "\t"};
            for (i>9; i<11; i++){
                if  ($i =="1/1") {printf $i"\t"};
                if  ($i=="1/0" ) {printf $i"\t"};
                if  ($i=="0/1" ) {printf $i"\t"};
                if  ($i=="0/0" ) {printf $i"\t"};
                };
            {printf "\n"}
            }
        }' > frankengenome.vcf
 #get header from vcf
cat France.vcf | \
     awk -F$'\t' '{
         if (NR<34){print$0}  
         }' > header.txt

#in ipython
ipython3
import pandas as pd 
df=pd.read_table("~/frankengenome.vcf", header=33, low_memory=False)
df
dropcol=df.columns[10:]
df=df.drop(dropcol, axis=1)
df.to_csv("~/for_rupa/Czech_1.csv", index=False, sep="\t") 

#combine to vcf
cat header.txt Czech_1.csv > Czech_1.vcf

#replace any gaps with ./.
cat Estonia_1.vcf | awk -F$'\t' '{
if (NR<35){print$0}
if (NR>34){
for (i=1; i<10; i++){printf $i "\t"};
            for (i>9; i<11; i++){
                if  ($i =="1/1") {printf $i"\t"};
                if  ($i=="1/0" ) {printf $i"\t"};
                if  ($i=="0/1" ) {printf $i"\t"};
                if  ($i=="0/0" ) {printf $i"\t"};
                if ($i==""){printf "./.""\t"};
                {printf "\n"}
                }
            }
        }' > Frankengenomes/Estonia_frankengenome.vcf

# old version
# cat E1_diploid_all_loci.vcf | \
#     awk -F$'\t' '{
#         if (NR>3384) {
#             for (i=1; i<10; i++){printf $i "\t"};
#             for (i>9; i<NF; i++){
#                 if  ($i =="1/1") {printf $i"\t"};
#                 if  ($i=="1/0" ) {printf $i"\t"};
#                 if  ($i=="0/1" ) {printf $i"\t"};
#                 if  ($i=="0/0" ) {printf $i"\t"};
#                 };
#             {printf "\n"}
#             }
#         }' > frankenstein1.csv
# #check that everything looks normal
# cat frankenstein.csv | tail -n 1000 | awk -F$'\t' '{
#     for (i=1; i<4; i++){printf $i "\t"}; {printf "\n"}
#     }' | less -N

# #get header from vcf
# cat E1_diploid_all_loci.vcf | \
#     awk -F$'\t' '{
#         if (NR<3384){print$0}  
#         }' > header.txt

# #add the header back onto the frankenstein stitched vcf
# {
#     cat header.txt;
#     cat frankenstein.csv;
# } > frankenstein.vcf

# #check that the vcf reads in as a vcf 
# bcftools query -l frankenstein.vcf | less -N  

# #this includes the header from the original
# cat reich_aDNA.vcf | \
#     awk -F$'\t' '{
#         if (NR<29){
#             print $0
#         }
#         else if (NR>29) {
#             for (i=1; i<10; i++){printf $i "\t"};
#             for (i>9; i<NF; i++){
#                 if  ($i =="1/1") {printf $i"\t"};
#                 if  ($i=="1/0" ) {printf $i"\t"};
#                 if  ($i=="0/1" ) {printf $i"\t"};
#                 if  ($i=="0/0" ) {printf $i"\t"};
#             };
#             {printf "\n"}
#     }}' > stitched_frankenstein.csv
# #this includes the original sample names -- we don't want that
# cat reich_aDNA.vcf | \
#     awk -F$'\t' '{
#         if (NR<30){
#             print $0
#         }
#         else if (NR>29) {
#             for (i=1; i<NF; i++){
#                 {printf $1 "\t"}
#                 if  ($i =="1/1") {printf $i"\t"};
#                 if  ($i=="1/0" ) {printf $i"\t"};
#                 if  ($i=="0/1" ) {printf $i"\t"};
#                 if  ($i=="0/0" ) {printf $i"\t"};
#             };
#             {printf "\n"}
#     }}' | less -N 
