"""
Converts Cufflinks gtf to a sorted bed.
"""
#!/usr/bin/env python
# encoding: utf-8

import argparse
import os
import sys
import gzip
import bz2
import urllib
from itertools import izip
from re import findall
from collections import defaultdict


def nopen(f, mode="rb"):
    # author: github.com/brentp
    if not isinstance(f, basestring):
        return f
    if f.startswith("|"):
        p = Popen(f[1:], stdout=PIPE, stdin=PIPE, shell=True)
        if mode[0] == "r": return p.stdout
        return p
    return {"r": sys.stdin, "w": sys.stdout}[mode[0]] if f == "-" \
         else gzip.open(f, mode) if f.endswith((".gz", ".Z", ".z")) \
         else bz2.BZ2File(f, mode) if f.endswith((".bz", ".bz2", ".bzip2")) \
         else urllib.urlopen(f) if f.startswith(("http://", "https://",
             "ftp://")) \
        else open(f, mode)

 
def reader(fname, header=True, sep="\t"):
    # author: github.com/brentp
    line_gen = (l.rstrip("\r\n").split(sep) for l in nopen(fname))
    if header == True:
        header = line_gen.next()
        header[0] = header[0].lstrip("#")

    if header:
        for toks in line_gen:
            yield dict(izip(header, toks))
    else:
        for toks in line_gen:
            yield toks


def makedicts(gtf):
    """Returns transcript dictionary and exon defaultdict using transcript id
as the key
"""
    exons = defaultdict(list)
    transcripts = {}
    for line in reader(gtf, header=False):
        transid = findall(r'transcript_id \"([\w\.]+)\"',line[8])[0].strip()
        #fpkmval = findall(r'FPKM \"([\d\s\w\.]+)\"',line[8])[0].strip()
        # 0-chrome, 3-start, 4-end, 6-strand
        if(line[2] == 'transcript'):
            transcripts[transid] = [line[0], str(int(line[3])-1), line[4], line[6]]
        if(line[2] == 'exon'):
            exons[transid].append([str(int(line[3])-1), line[4]])
    return {'transcripts':transcripts, 'exons':exons}


def gtf2bed(gtf):
    data = makedicts(gtf)
    # full bed header
    header = ['chrom', 'start', 'end', 'name', 'score', 'strand', \
                'thickStart', 'thickEnd', 'itemRGB', 'blockCount', \
                'blockSizes', 'blockStarts']
    small_head = ['chrom','start','end','strand','thickStart','thickEnd']
    positions = [0, 1, 2, 3, 1, 2]
    for k, v in sorted(data['transcripts'].iteritems(), \
                        key=lambda (k,v): [v[0], int(v[1])]):
        print_line = {}
        sizes = []
        starts = []
        for t, p in izip(small_head, positions):
            print_line[t] = v[p]
        print_line['itemRGB'] = '0'
        print_line['name'] = k
        fpkm = 0

        print_line['score'] = '%d' % round(float(fpkm))
        for e in data['exons'][k]:
            sizes.append('%d' % (int(e[1]) - int(e[0])))
            starts.append('%d' % (int(e[0]) - int(print_line['start'])))
        if int(starts[0])!=0:
        	sizes=sizes[::-1]
        	starts=starts[::-1]
        print_line['blockCount'] = '%d' % len(sizes)
        print_line['blockSizes'] = ",".join(s for s in sizes)
        print_line['blockStarts'] = ",".join(s for s in starts)
        print "\t".join(print_line[h] for h in header)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-g", "--gtf", dest="gtf", \
                    help="Cufflinks .gtf to convert", type=str)
    args = p.parse_args()
    if not (args.gtf):
        sys.exit(p.print_help())
    gtf2bed(args.gtf)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
