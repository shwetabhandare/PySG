import sys



def StemLoopSequence(seq, annotation):
	stems = []
	loopSeq = "";
	for i, j in zip(seq, annotation):	
		if j == '(':
			stems.append(i)
		elif j == ')':
			print stems[-1], "-", loopSeq, "-", i;
			del stems[-1]
		else:
			loopSeq = loopSeq + i;




if __name__ == "__main__":
	import sys
	seq = sys.argv[1]
	annotation = sys.argv[2]
	StemLoopSequence(seq, annotation);