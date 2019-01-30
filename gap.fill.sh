ref1=$1
ref2=$2
iTools Fatools findN -InPut ${ref1}|awk '{print $1"\t"$2-1"\t"$3+1}' > ${ref1}.gap
bedtools merge -d 500 -i ${ref1}.gap > ${ref1}.gap
les ${ref1}.gap|awk '{print $1"\t"$2-500"\t"$2}' > temp.gap.1.bed
les ${ref1}.gap|awk '{print $1"\t"$3"\t"$3+500}' > temp.gap.2.bed
fastaDeal.pl -sub temp.gap.1.bed ${ref1} > temp.gap.1.bed.fa
fastaDeal.pl -sub temp.gap.2.bed ${ref1} > temp.gap.2.bed.fa
fastaDeal.pl -reform line500000 temp.gap.1.bed.fa > temp.gap.1.bed.line.fa
fastaDeal.pl -reform line500000 temp.gap.2.bed.fa > temp.gap.2.bed.line.fa
awk -F "" '{if(NR%2==1){print "@"$0}else{printf ($0"\n+\n");for(i=1;i<=NF;i++){printf "I"};printf "\n"}}' temp.gap.1.bed.line.fa   |sed 's/@>/@/'  > temp.gap.1.bed.line.fq
awk -F "" '{if(NR%2==1){print "@"$0}else{printf ($0"\n+\n");for(i=1;i<=NF;i++){printf "I"};printf "\n"}}'  temp.gap.2.bed.line.fa |sed 's/@>/@/'  > temp.gap.2.bed.line.fq
paste temp.gap.1.bed temp.gap.2.bed |awk '{t="@"$1"_"$2"_"$3"_"$5"_"$6; print "@"$1"_"$2"_"$3"\t"t"\n@"$1"_"$5"_"$6"\t"t}'  > temp.gap.bed.id
awk 'BEGIN{while (getline<"temp.gap.bed.id"){a[$1]=$2}}{if($1~/^@/){print a[$1]}else{print $0}}' temp.gap.1.bed.line.fq > temp.gap.1.bed.line.rename.fq
awk 'BEGIN{while (getline<"temp.gap.bed.id"){a[$1]=$2}}{if($1~/^@/){print a[$1]}else{print $0}}'  temp.gap.2.bed.line.fq > temp.gap.2.bed.line.rename.fq
bwa mem -a -t 10 ${ref2} temp.gap.1.bed.line.rename.fq temp.gap.2.bed.line.rename.fq > temp.gap.bed.line.rename.fq.sam
python gap.fill.step1.py temp.gap.bed.line.rename.fq.sam > temp.gap.bed.line.rename.fq.sam.gap
python gap.fill.step2.py temp.gap.bed.line.rename.fq.sam.gap  ${ref1}  > {ref1}.fasta 