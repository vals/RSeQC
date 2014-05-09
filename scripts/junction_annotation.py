#!/usr/bin/env python
'''-------------------------------------------------------------------------------------------------
Annotate splicing reads against gene model in two levels: splicing events level and splicing juncion
level. Note that a read (especially long read) can be spliced 2 or more times, each time is called a
splicing event; therefore, 100 spliced reads can produce more than 100 splicing events. On the other
hand, multiple splicing events spanning the same intron can be consolidated into one splicing junction.
-------------------------------------------------------------------------------------------------'''

#import built-in modules
import os,sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
	print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " RSeQC needs python2.7!\n"
	sys.exit()

import re
import string
from optparse import OptionParser
import warnings
import string
import collections
import math
import sets
from time import strftime
import subprocess

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *

#import my own modules
from qcmodule import SAM
#changes to the paths

#changing history to this module


__author__ = "Liguo Wang"
__copyright__ = "Copyright 2012. All rights reserved."
__credits__ = []
__license__ = "GPL"
__version__="2.3.9"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def printlog (mesg):
        '''print progress into stderr and log file'''
        mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
        LOG=open('class.log','a')
        print >>sys.stderr,mesg
        print >>LOG,mesg


def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Alignment file in BAM or SAM format.")
	parser.add_option("-r","--refgene",action="store",type="string",dest="ref_gene_model",help="Reference gene model in bed format. This file is better to be a pooled gene model as it will be used to annotate splicing junctions [required]")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files(s). [required]")
	parser.add_option("-m","--min-intron",action="store",type="int",dest="min_intron",default=50, help="Minimum intron length (bp). default=%default [optional]")
	parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality (phred scaled) for an alignment to be called \"uniquely mapped\". default=%default")

	(options,args)=parser.parse_args()
		
	if not (options.output_prefix and options.input_file and options.ref_gene_model):
		parser.print_help()
		sys.exit(0)
	if not os.path.exists(options.ref_gene_model):
		print >>sys.stderr, '\n\n' + options.ref_gene_model + " does NOT exists" + '\n'
		sys.exit(0)
	if os.path.exists(options.input_file):
		obj = SAM.ParseBAM(options.input_file)
		obj.annotate_junction(outfile=options.output_prefix,refgene=options.ref_gene_model,min_intron=options.min_intron, q_cut = options.map_qual)
		try:
			subprocess.call("Rscript " + options.output_prefix + '.junction_plot.r', shell=True)
		except:
			print >>sys.stderr, "Cannot generate pdf file from " + '.junction_plot.r'
			pass
	else:
		print >>sys.stderr, '\n\n' + options.input_file + " does NOT exists" + '\n'
		sys.exit(0)			


if __name__ == '__main__':
        main()
