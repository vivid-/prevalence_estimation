from bs4 import BeautifulSoup
import urllib2
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--gene",help="the gene symbol for search",type=str)
parser.add_argument("--out",help="the output destination for the search term",type=str)
args = parser.parse_args()
gene = args.gene
out = args.out

url_gene = 'http://www.egl-eurofins.com/emvclass/emvclass.php?approved_symbol=' + gene
# fetch the page
page = urllib2.urlopen(url_gene)
# parge the webpage html
soup = BeautifulSoup(page,'html.parser')

# firstly, find contents with td tags
tb = soup.find_all('td')

# get the result table 
df_result = pd.DataFrame(columns = ['Order','Gene','Exon','Nucleotide_Change','Protein_Change',
                                    'Alias_Listing','Classification','Last_Reviewed'])
# for each content tag
j = 1
i = 0
while i<len(tb):
    # find the starting line for the result table
    if tb[i].string == str(j):
        #i_tmp = i+1
        #print(i)
        #while tb[i_tmp].string != str(j+1):
        df_result.loc[j-1]=[j, tb[i+1].string, tb[i+2].string,
                           tb[i+3].string, tb[i+4].string, tb[i+5].string,
                           tb[i+6].string, tb[i+7].string]
        i = i+8
        j = j+1
    else:
        i = i+1


df_result.to_csv(out,sep='\t',index=False,header=True, encoding='utf-8')
