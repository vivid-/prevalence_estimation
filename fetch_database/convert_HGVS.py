import hgvs.parser
import pandas as pd
import argparse
import hgvs.dataproviders.uta
import hgvs.assemblymapper

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
    g = am.c_to_g(var_c)
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

pos = []
ref = []
alt = []
# load the input file
df = pd.read_csv(inpf,sep='\t')
for i in range(df.shape[0]):
    start,ref_tmp,alt_tmp = parse_hgvs_cdna(df['Nucleotide_Change'][i],am)
    pos.append(start)
    ref.append(ref_tmp)
    alt.append(alt_tmp)

result = df
result['pos']=pos
result['ref']=ref
result['alt']=alt

result.to_csv(outf,sep='\t',index=False,header=True)
