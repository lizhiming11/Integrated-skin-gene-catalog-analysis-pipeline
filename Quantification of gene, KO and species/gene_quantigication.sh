2bwt-builder gene_catalog.fas && touch gene_catalog.fas.index
fasta_length -i gene_catalog.fas | awk '{print NR"\t"$0}'>  gene_catalog.fas.length

sample = $1
perl abundance.pl -a $sample_rmHost.1.fq.gz -b $sample_rmHost.2.fq.gz -d gene_catalog.fas.index -g  gene_catalog.fas.length -p ./abundance/$sample

profile_v1.0 -i skin.list -p gene_profile


