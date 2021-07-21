mkdir tmp &&
cat data/reference/*.fasta > data/reference/tmp.fasta &&
mmseqs easy-linclust -e Inf -c 0.01 --spaced-kmer-mode 1 \
--spaced-kmer-pattern 110101 data/reference/tmp.fasta results tmp &&
rm data/reference/tmp.fasta