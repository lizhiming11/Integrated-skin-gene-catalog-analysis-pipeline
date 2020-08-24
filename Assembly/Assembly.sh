fq1 = $1
fq2 = $2
outfile = $3

spades.py --meta -1 $fq1 -2 $fq2  -o $outfile --threads 20 --memory 80 -k 21,33,55,77,99

