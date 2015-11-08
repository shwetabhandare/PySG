#Python script for reading file, counting number of sequences, and concatenating +ve to -ve files
#Has no functions but I can do a proper program with main, etc.
import sys
import os

#Taking in filenames as arguments, +ve and then -ve




def CreateCombinedFile(posFile, negFile):
## This needs to be completed.

	filename_plus = posFile;
	filename_minus = negFile;
	
	#Open files as read-only
	file_plus = open(filename_plus, "r")
	file_minus = open(filename_minus, "r")
	
	#Creating "text" of the files to use
	plus_text = file_plus.read()
	minus_text = file_minus.read()
	
	plus_count = 0
	minus_count = 0
	
	#Split and search for the character '>', which denotes new sequence
	plus_split = plus_text.split()
	minus_split = minus_text.split()
	
	for word in plus_split:
		if '>' in word:
			plus_count = plus_count +1
	
	for word in minus_split:
		if '>' in word:
			minus_count = minus_count+1
	
	#Creating output file and writing other files to it
	just_plus_basename = os.path.basename(filename_plus)
	just_minus_basename = os.path.basename(filename_minus)
	just_plus = os.path.splitext(just_plus_basename)[0]		#splitting the filename from extension and only taking the filename
	just_minus = os.path.splitext(just_minus_basename)[0]		#splitting the filename from extension and only taking the filename
	#just_plus = os.path.splitext(filename_plus)[0]		#splitting the filename from extension and only taking the filename
	#just_minus = os.path.splitext(filename_minus)[0]
	pcount = str(plus_count)
	mcount = str(minus_count)
	filename_output = ""+just_plus+"_"+just_minus+"_"+pcount+"_"+mcount+"_CompleteSet.txt"	#creating the file name for the output file
	outputfile = open(filename_output, "w")	
	outputfile.write(plus_text)
	outputfile.write(minus_text)

	print "Combined Filename:", filename_output
	#Closing files
	file_plus.close()
	file_minus.close()
	outputfile.close()

	return (filename_output, plus_count, minus_count)

if __name__ == "__main__":
	import sys
	posFile = sys.argv[1]
	negFile = sys.argv[2]
	print posFile, negFile
	CreateCombinedFile(posFile, negFile);

