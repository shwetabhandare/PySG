from random import choice
def weightedchoice(items): # this doesn't require the numbers to add up to 100
	return choice("".join(x * y for x, y in items))



#seq = GetRandomSequence(50);
#print seq;
