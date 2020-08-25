
perl split_fa.pl s 9 ../gene_catalog.fas split
blastn -word_size 16 -query split_1.fa -out split_1.bt -db NCBI-NT.fa -evalue 1e-10 -outfmt 6 -num_threads 20
perl get_length.pl gene_catalog.fas.pep gene_catalog.fas.pep.len
