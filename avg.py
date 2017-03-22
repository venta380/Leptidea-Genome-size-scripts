import sys
import string
import numpy
infile= str(sys.argv[1])
table = [int(line.strip()) for line in open(infile)]

print numpy.average(table)
