from random import choice
import os, fnmatch
def weightedchoice(items): # this doesn't require the numbers to add up to 100
	return choice("".join(x * y for x, y in items))


def findFiles (path, filter):
   for root, dirs, files in os.walk(path):
      for file in fnmatch.filter(files, filter):
         yield os.path.join(root, file)
#seq = GetRandomSequence(50);
#print seq;
