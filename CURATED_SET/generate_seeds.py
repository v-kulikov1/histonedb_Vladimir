from curated_set_services import CuratedSet
from curated_set_services import dict2tree, muscle_p2p_aln

# from Bio.Align.AlignInfo import SummaryInfo
# from Bio import AlignIO
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq

from ete3 import Tree
import json, os


cs = CuratedSet()
with open('classification.json') as json_file:
    data = json.load(json_file)
hist_tree=Tree()
dict2tree(hist_tree,data['tree'])
# print(hist_tree.get_ascii(show_internal=True))
draft_seeds_msa={}
folder = '../data/draft_seeds/'
if not os.path.exists(folder): os.makedirs(folder)
#generate draft seeds - needs debugging
for node in hist_tree.traverse("postorder"):
    print("Processing ",node.name)
    if(node.is_leaf()): # we get sequences for that variant and align them.
        draft_seeds_msa[node.name]=cs.muscle_aln(cs.data.query(f'type=="{node.name}" | variant.str.contains(r"(^|\\s){node.name}($|\\s)")', engine='python')['accession'],debug=False)
        print(node.name,"Alignment length:",len(draft_seeds_msa[node.name]))
        with open(f'{folder}{node.name}.fasta','w') as f:
            f.write(draft_seeds_msa[node.name].format('fasta'))
    elif not node.is_root(): # we will do profile to profile alignment
        #we should first check if there are seqs with this subrariant as the most specific one
        print(f"\t Node is internal, progressive alignment:")
        dq = cs.data.query(f'variant.str.contains(r"(^|\\s){node.name}($|\\s)")')['accession']
        msa=cs.muscle_aln(cs.data.query(f'variant.str.contains(r"(^|\\s){node.name}($|\\s)")')['accession'],debug=True)
        draft_seeds_msa[node.name+'_only']=msa
        print(f"\t\t For {node.name} aligned {len(msa)} sequences")
        ch=node.get_children() #get children
#         msa=draft_seeds_msa[ch[0].name]
        #progressively align
        for i in range(len(ch)):
            if(len(msa)==0):
                msa=draft_seeds_msa[ch[i].name]
                print(f"\t\t Adding child {node.name} aligned {len(draft_seeds_msa[ch[i].name])} sequences")
                continue
            elif(len(draft_seeds_msa[ch[i].name])!=0):
                msa=muscle_p2p_aln(msa,draft_seeds_msa[ch[i].name])
                print(f"\t\t Adding child {node.name} aligned {len(draft_seeds_msa[ch[i].name])} sequences")
            else:
                continue
        draft_seeds_msa[node.name]=msa
        print(node.name,"Alignment length:",len(draft_seeds_msa[node.name]))
        with open(f'{folder}{node.name}.fasta','w') as f:
            f.write(draft_seeds_msa[node.name].format('fasta'))
        with open(f'{folder}{node.name}_only.fasta','w') as f:
            f.write(draft_seeds_msa[node.name+'_only'].format('fasta'))
#         print(f"\t\t Final for {node.name} aligned {len(draft_seeds_msa[node.name])} sequences")
    else:
        continue