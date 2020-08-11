
kmer=5
meta_p=../../bin/metacortex_k31

export R_ENV_PATH='../'

b=100
n=15
g=100
max_edges=8
delta=1000
kmers_cut=0
#no_walk='_no_walk'; M=''
no_walk=''; M='-M'

filename="reads/disconnected.fa"
#filename="reads/polymorphism_snp1.fa"
#filename="reads/polymorphism_snp2.fa"
#filename="reads/polymorphism_snp3.fa"
#filename="reads/polymorphism_indel.fa"
#filename="reads/polymorphism_indel2.fa"
#filename="reads/polymorphism_indel3.fa"
#filename="reads/polymorphism_none.fa"
#filename="reads/polymorphism_double.fa"
#filename="reads/polymorphism_two_snps.fa"
#filename="reads/polymorphism_allele.fa"
#filename="reads/polymorphism_allele2.fa"
#filename="reads/polymorphism_allele3.fa"
#filename="reads/polymorphism_allele_long.fa"
#filename="reads/two_paths.fa"
#filename="reads/polymorphism_reversecompl.fa"
#filename="reads/polymorphism_reversecompl2.fa"

name=`basename ${filename}`
name=`echo ${name%.*}`

cortex_file="${name}.ctx"
contig_file="${name}.fa"
log_file="${name}.txt"
file_list="allfiles.txt"
graph_file="${name}.gv"

if [ -f ${cortex_file} ] ; then
	rm ${cortex_file}
fi
if [ -f ${contig_file} ] ; then
        rm ${contig_file}
fi

mkdir graphs; chmod 755 graphs

echo ${filename} > ${file_list}

${meta_p} -k ${kmer} -n ${n} -b ${b} -i ${file_list} -t fasta -o ${cortex_file} -f ${contig_file} -g 5 -l ${log_file} -r ${max_edges} -R ${delta} -y ${kmers_cut}  -S ${M} -G ${graph_file}
dot -Tdot ${graph_file} -o ${name}.dot
dot -Tpdf ${name}.dot -o ${name}.pdf