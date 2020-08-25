
perl split_fa.pl s 9 ../gene_catalog.fas split
blastn -word_size 16 -query split_1.fa -out split_1.bt -db NCBI-NT.fa -evalue 1e-10 -outfmt 6 -num_threads 20
perl get_length.pl gene_catalog.fas.pep gene_catalog.fas.pep.len
filter_blast -i split_1.bt -o split_1.bt.f --qfile gene_catalog.fas.pep.len --qper 70 --tops 5
perl get_species_anno.v2019.pl split_1.bt.f split_1.bt.f.tax 
