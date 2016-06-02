import sys
import SeqGenUtils


from tempfile import mkstemp
from shutil import move
from os import remove, close

def replace(file_path, pattern, subst):
	#Create temp file
	fh, abs_path = mkstemp()
	with open(abs_path,'w') as new_file:
		with open(file_path) as old_file:
			for line in old_file:
				if line[0] == ">":
					new_file.write(line)
					continue;
				else:
					new_file.write(line.replace(pattern, subst))
	close(fh)
	#Remove original file
	remove(file_path)
	#Move new file
	move(abs_path, file_path)

if __name__ == "__main__":
	import sys
	fastaFile = sys.argv[1]
	replace(fastaFile, "U", "T");


