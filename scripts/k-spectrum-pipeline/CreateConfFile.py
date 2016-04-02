import sys
import string
import yaml

def CreateConf(filename, output_prefix, confFile):
	filename_list = string.split(filename, "_")
	list_len = len(filename_list)

	print "FileName: ", filename
	print filename_list

	posLen = filename_list[list_len - 3]
	negLen = filename_list[list_len - 2]
	featureFileName =  output_prefix + "_Features.dat"
	resultsFileName = output_prefix + "_Results.yml"
	modelFileName = output_prefix + "_Model.dat"
	trainDataFileName = output_prefix + "_kspectrum.dat"
	paramFileName = output_prefix + "_params.dat"



	data = dict(
		data = dict (
			dataFile = filename,
			posLen = int(posLen),
			negLen = int(negLen),
			),
			kspectrum = dict(
				k1 = 8,
				k2 = 14,
				C = 10,
				folds = 5,
			),
			output = dict( 
				featureFile = featureFileName,
				resultsFile = resultsFileName,
				modelFile = modelFileName,
				trainDataFile = trainDataFileName,
				paramsFile = paramFileName,
				)
		)

	with open(confFile, 'w+') as outfile:
		outfile.write(yaml.dump(data, default_flow_style=True))
	outfile.close()

if __name__ == '__main__':
	import sys
	combinedFile = sys.argv[1]
	print "Combined File: ", combinedFile;
	prefix = sys.argv[2]
	confFile = sys.argv[3]
	CreateConf(combinedFile, prefix, confFile)
