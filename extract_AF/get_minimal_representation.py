#!/usr/bin/pythom
import pandas as pd
import argparse
########################################################
## this script is used to get the mini representation ##
## of VCF files ########################################
########################################################

def get_minimal_representation(pos, ref, alt): 
	# If it's a simple SNV, don't remap anything
	if len(ref) == 1 and len(alt) == 1: 
		return pos, ref, alt
	else:
		# strip off identical suffixes
		while(alt[-1] == ref[-1] and min(len(alt),len(ref)) > 1):
			alt = alt[:-1]
			ref = ref[:-1]
		# strip off identical prefixes and increment position
		while(alt[0] == ref[0] and min(len(alt),len(ref)) > 1):
			alt = alt[1:]
			ref = ref[1:]
			pos += 1
		#print 'returning: ', pos, ref, alt
		return pos, ref, alt 



# get the arguments
parser = argparse.ArgumentParser()
parser.add_argument("--INPUT",help="the location of input file",type=str)
parser.add_argument("--OUTPUT",help="the location of output file",type=str)
args = parser.parse_args()
input_file = args.INPUT
output_file = args.OUTPUT

# read in the data
vcf = pd.read_csv(input_file,sep='\t')

# get minimal represented pos ref and alt
for i in range(vcf.shape[0]):
    pos_updated, ref_updated, alt_updated = get_minimal_representation(vcf['POS'][i], vcf['REF'][i], vcf['ALT'][i])
    vcf.loc[i,'POS'] = pos_updated 
    vcf.loc[i,'REF'] = ref_updated
    vcf.loc[i,'ALT'] = alt_updated   
#    print(str(vcf['POS'][i])+'\t'+str(vcf['REF'][i])+'\t'+str(vcf['ALT'][i]))
# update the pos, ref and alt columns
#vcf['POS'] = pos_updated
#vcf['REF'] = ref_updated
#vcf['ALT'] = alt_updated

# output the file
vcf.to_csv(output_file,sep='\t',index=False)
