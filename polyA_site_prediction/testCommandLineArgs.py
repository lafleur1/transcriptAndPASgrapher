#test command line arguments 

import sys

n = len(sys.argv)
print ("TOTAL ARGS: ", n)


others = []
for arg in sys.argv[1:]:
	others.append(arg)

print (others)
