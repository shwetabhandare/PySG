from random import choice
def weightedchoice(items): # this doesn't require the numbers to add up to 100
	return choice("".join(x * y for x, y in items))

def GetRandomSequence(seqLen):
	seq=""
	for count in range(seqLen):
		seq+=weightedchoice([("C", 10), ("G", 20), ("A", 40), ("T", 30)]);
	return seq;


#seq = GetRandomSequence(50);
#print seq;