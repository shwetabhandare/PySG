import sys
import itertools

a = [[1,0,1,1], [0,0,0,0]]
b = [[0,1,0,1], [1,0,0,1]]

l = []
for i, j in itertools.izip(a,b):
	print (i,j)
	p = []
	for x,y in itertools.izip(i,j):
		print (x,y)
		if (x==y):
			p.append(1);
		if (x != y):
			p.append(0);
	print p
	l.append(p);
print l
