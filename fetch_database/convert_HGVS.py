import hgvs.parser
import pandas as pd
import argparse
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import re
import sys
import pickle

#hgvs_c = 'NM_001637.3:c.1582G>A'
# deal with parameters
parser = argparse.ArgumentParser()
parser.add_argument("--inpu",help="the location of the input file",type=str)
parser.add_argument("--out",help="the output destination for the search term",type=str)
args = parser.parse_args()
inpf = args.inpu
outf = args.out




# parse the cDNA variants to the genomic variants
def parse_hgvs_cdna (hgvs_c,am):
    # parse the genomic variant into a Python structure
    hp = hgvs.parser.Parser()
    var_c = hp.parse_hgvs_variant(hgvs_c)
    # convert the cDNA variants to the genomic variants
    try:
        g = am.c_to_g(var_c)
    except hgvs.exceptions.HGVSInvalidVariantError:
        print(str(var_c)+" Invalid Variant Error")
        return
    except hgvs.exceptions.HGVSDataNotAvailableError:
        print(str(var_c)+" Data Not Available")
        return
    except:
        print ('Unexpected error: '+ str(sys.exc_info()[0]))
        return
    start = g.posedit.pos.start.base
   # print(g.posedit.edit)
    ref = g.posedit.edit.ref
    if str(g.posedit.edit) == 'dup':
        alt = ref + ref
    elif str(g.posedit.edit) == 'del':
        alt = ''
    else:
        alt = g.posedit.edit.alt
    
    return start,ref,alt
    


hdp=hgvs.dataproviders.uta.connect()
# initialize the mapper for GRCh37 with splign-based alignments
am = hgvs.assemblymapper.AssemblyMapper(hdp,
    assembly_name='GRCh37', alt_aln_method='splign',
    replace_reference=True)
# store the two objects above for later usage
#filehandler = open('GRCh37_am', 'w') 
#pickle.dump([hdp,am], filehandler)

pos = []
ref = []
alt = []
# load the input file
df = pd.read_csv(inpf,sep='\t')
# remove artifacts
#df = df.drop(df[df['Classification']=='Artifact'].index)
#df = df.reindex(range(df.shape[0]))
df = df[df.Classification!='Artifact']

for i in df.index.tolist():
    # in case that the duplication length and refseq length doesn't match
    # just skip those variants
    #print(str(i)+'\n')
    if 'dup' in df['Nucleotide_Change'][i]:
        pat = ':c.(\d+_\d+)'+'dup'+'(\S+)'
        matched = re.search(pat,df['Nucleotide_Change'][i])
        # started position annotated
        if matched:
            pos_int = matched.group(1)
            length = int(re.split('_',pos_int)[1]) - int(re.split('_',pos_int)[0])
	    if length != len(matched.group(2)):
                pos.append(0)
                ref.append('')
                alt.append('')
                continue
    if 'del' in df['Nucleotide_Change'][i]:
        pat = ':c.(\d+_\d+)'+'del'+'(\S+)'
        matched = re.search(pat,df['Nucleotide_Change'][i])
        if matched:
            # started position annotated
            pos_int = matched.group(1)
            length = int(re.split('_',pos_int)[1]) - int(re.split('_',pos_int)[0])
            if length != len(matched.group(2)):
                pos.append(0)
                ref.append('')
                alt.append('')
                continue
 
    if parse_hgvs_cdna(df['Nucleotide_Change'][i],am) is not None:
        start,ref_tmp,alt_tmp = parse_hgvs_cdna(df['Nucleotide_Change'][i],am)
        pos.append(start)
        ref.append(ref_tmp)
        alt.append(alt_tmp)
    else:
        pos.append(0)
        ref.append('')
        alt.append('')

result = df
result['pos']=pos
result['ref']=ref
result['alt']=alt

result.to_csv(outf,sep='\t',index=False,header=True)
