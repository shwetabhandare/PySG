import sys



def StemLoopSequence(seq, annotation):
	stems = []
	loops = []
	loopSeq = "";
	for i, j in zip(seq, annotation):	
		if j == '(':
			stems.append(i)
			if  loopSeq:
				print loopSeq;
				loops.append(loopSeq);
				loopSeq = ""
		elif j == ')':
			loops.append(loopSeq);
			loopSeq = "";
			print stems[-1], "-", loops[-1], "-", i;
			del stems[-1]
			del loops[-1]
		else:	
			loopSeq = loopSeq + i;

	print ''.join(loops)

if __name__ == "__main__":
	import sys
	seq = sys.argv[1]
	annotation = sys.argv[2]
	StemLoopSequence(seq, annotation);