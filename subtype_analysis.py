from os import listdir
from os.path import isfile, join
import json
import pandas as pd
import numpy as np
import scipy.stats
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('subtype', type=str, help='subtype name')
parser.add_argument('pathway', type=str, help='pathway name')
args = parser.parse_args()

# reading the list of TFs
GRN = pd.read_csv('GRN_symbol.txt', sep='\t')
TFs = set(GRN['TF_sym'].unique())

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

# generating the list of goal entities (DEGs)
lim_res = pd.read_csv('res/Limma_Normal-'+args.subtype+'.txt', delimiter='\t')
DEGs = lim_res.index[(lim_res['adj.P.Val']<0.05).values]
DEGs_ent = set(map(lambda x : '' if x not in sym2ent else sym2ent[x], DEGs))
DEGs_ent.discard('')

ents = dish2ent[args.pathway]
with open('res/DEGs_in_'+args.subtype+'-'+args.pathway+'.txt', 'w') as outF : 
  outF.write('\n'.join(list(DEGs_ent & ents)))

# enrichment tests for TFs
lim_res2 = pd.read_csv('res/Limma_'+args.subtype+'-Others.txt', delimiter='\t')
DEGs2 = lim_res2.index[(lim_res2['adj.P.Val']<0.05).values]

TFs_ent = set(['' if t not in sym2ent else sym2ent[t] for t in TFs])
TFs_path = TFs_ent & ents

enrich_res = []
for TF_ent in TFs_path :
  TF_sym = ent2sym[TF_ent]
  TF_sym = TF_sym if isinstance(TF_sym, list) else [TF_sym]
  TGs = set(GRN.loc[GRN['TF_sym'].isin(TF_sym), 'Gene_sym'])
  n = [len(DEGs2 & TGs), len(TGs), DEGs2.shape[0], lim_res2.shape[0]]
  o, p = scipy.stats.fisher_exact([[n[0], n[1]-n[0]], [n[2]-n[0], n[3]-n[1]-n[2]+n[0]]], alternative='greater')
  enrich_res.append([TF_ent, ','.join(TF_sym)]+n+[p])
 
enrich_res = pd.DataFrame(enrich_res, columns=['TF_ent', 'TF_sym', 'DE-TG', 'TG','DEG','Gene', 'P-value'])
enrich_res = enrich_res.sort_values(by=['P-value'])
enrich_res.to_csv('res/Subtype-specific_TFs_in_'+args.subtype+'-'+args.pathway+'.txt', sep='\t', index=False)

