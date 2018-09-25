# input  snp file format: chrom position
# output bed file : standard bed

import sys, os

def proc_snp(in_fh, out_fh, bin_size):
    first_line = True
    for line in in_fh:
        if first_line: #skip the headers
            first_line = False
            continue
        fields = line.split('\t')
        chrom = fields[0]
        pos   = int(fields[1])
        start = int(pos - int((bin_size)/2))
        end   = start + bin_size
        assert(pos == int((start + end)/2)
        out_fh.write(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + line)

#in_fn  = sys.argv[1]
#out_fn = sys.argv[2]
#with open(in_fn, 'r') as in_fh, open(out_fn, 'w') as out_fh:

if len(sys.argv) == 1:
    bin_size = 1000
elif len(sys.argv) == 2:
    bin_size = int(sys.argv[1])
else:
    sys.stderr.write("Wrong arguments: " + argv[0] + " [bin_size]\n")
    quit()

proc_snp(sys.stdin, sys.stdout, bin_size)
