from os import listdir
from os.path import isfile, join
import re
import numpy as np
import pandas as pd
import json
import scipy.stats
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('sampleF', type=str, help='sample information file')
args = parser.parse_args()

# Reading the sample info file
subtypes = np.genfromtxt(args.sampleF, skip_header=1, usecols=[1], dtype=np.str)
subtypes = np.unique(subtypes)

# Reading the occurences of 32 pathways in PL database
dishes = ['Csf1', 'Egf', 'Hgf', 'Ifnab', 'Ifng', 'Igf1', 'IL1', 'IL2', 'IL4', 'IL6', 'IL12', 'Ins', 'Ngf', 'Pdgf', 'Tgfb1', 'Tnf', 'Vegf', 'Cd40', 'Serum', 'Lps', 'PolyIC', 'TLR9', 'Adriamycin', 'Aniso', 'Bleomycin', 'Etoposide', 'Hydroxyurea', 'NCS', 'PMA', 'Sorbitol', 'UV', 'IRad']

sym2ent = {}
with open('data/entities.json', 'r') as jsonF:
  ents = json.load(jsonF)
ents_prot = {key:val for key, val in ents.items() if 'etype' in val.keys() and val['etype']=='Protein'}
ent2sym = {key:val['metadata']['hugosym'] for key, val in ents_prot.items() if 'metadata' in val.keys() and 'hugosym' in val['metadata']}
sym2ent = {}
for k, v in ent2sym.items() :
  if isinstance(v, list) :
    for sym in v : sym2ent[sym]=k
  else : sym2ent[v]=k

dishPath = 'data/dishes'
dish2ent = {}
dishFs = [f for f in listdir(dishPath) if isfile(join(dishPath, f)) and f[-5:]=='.json' and (f.split('Dish')[0] in dishes)]
for dishF in dishFs :
  dish_name = dishF.split('Dish')[0]
  with open(join(dishPath, dishF), 'r') as dishF :
    dish_json = json.load(dishF)
    dish2ent[dish_name] = set(ents_prot.keys())  & set(map(lambda x : x.split('@')[0].split('-')[0], dish_json['init']))

# Enrichment test
enrich_res = []
for sub in subtypes :
  # Reading the list of DEGs
  if sub.upper()=="NORMAL" : continue
  lim_res = pd.read_csv('res/Limma_Normal-'+sub+'.txt', delimiter='\t')
  DEGs = lim_res.index[(lim_res['adj.P.Val']<0.05).values]
  DEGs_ent = set(map(lambda x : '' if x not in sym2ent else sym2ent[x], DEGs))
  DEGs_ent.discard('')
  
  # Enrichment test for 32 pathways
  for path, ents in dish2ent.items() : 
    n = [len(DEGs_ent & ents), len(ents), DEGs.shape[0], lim_res.shape[0]]
    o, p = scipy.stats.fisher_exact([[n[0], n[1]-n[0]], [n[2]-n[0], n[3]-n[1]-n[2]+n[0]]], alternative='greater')
    enrich_res.append([sub, path]+n+[p])
    #enrich_res.append(f'{sub}\t{path}\t{"\t".join(map(lambda x : str(x), n))}\t{p}')

#writing an output file
enrich_res = pd.DataFrame(enrich_res, columns=['Subtype','Pathway','DEG-Path', 'Gene-Path','DEG','Gene', 'P-value'])
enrich_res = enrich_res.sort_values(by=['Subtype', 'P-value'])
enrich_res.to_csv('res/Enrichment-test_results.txt', sep='\t', index=False)
