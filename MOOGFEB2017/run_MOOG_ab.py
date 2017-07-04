import os, sys

line = sys.argv
f = line[-1]

cmd = './MOOGSILENT << EOF\n%s\n\n\n\n\n\n\n\n\n\n\n\n\n\nEOF' % f 
#cmd = 'MOOGSILENT << EOF \n EOF' 
os.system(cmd)
