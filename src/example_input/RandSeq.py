import random
import sys

times = sys.argv[1]
s = sys.argv[2]

with open(s) as f:
        seq = f.read().splitlines()



#print "example input: python RandSeq.py 1000 sequence"

for i in range(int(times)):
	print ''.join(random.sample(seq[0],len(seq[0])))
