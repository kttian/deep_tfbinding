from __future__ import print_function, division

#field 0     1           2       3       4
#     93665  M2292_1.02  chr10   3485203 3485213 -   16.2182 2.05e-07    0.336   GGTGACTCATC JUND


#MAX_LOC = MMMKKKSSS
MAX_LOC  = 250000000L
def chrom_to_idx(chrom, loc):
    chrid = chrom[3:]
    if chrid == 'x' or chrid == 'X':
        chrnm = 23
    elif chrid == 'y' or chrid == 'Y':
        chrnm = 24
    else:
        try:
            chrnm = int(chrid)
            if chrnm > 22 or chrnm <= 0:
                return -1
        except:
            return -1
    return long(chrnm * MAX_LOC + loc)

import sys
in_fh = sys.stdin
in_fh.next()
for line in in_fh:
    f = line.split('\t')
    chrom = f[2]
    try:
        start = int(f[3])
        end   = int(f[4])
    except:
        continue
    pos   = int((start + end) / 2)
    idx   = chrom_to_idx(chrom, pos)
    if idx < 0:
        continue
    print(f[0], idx) # line_no, idx

