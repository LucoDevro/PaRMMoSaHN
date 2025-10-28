#! /bin/python

import sys
import pandas as pd

pan_reaction_table_file = sys.argv[1]
strain_blast_hits_file = sys.argv[2]
strain_reaction_table_file = sys.argv[3]

# pan_reaction_table_file = "no_singleton_protein_families-all-Reactions.tbl"
# strain_reaction_table_file = "no_singleton_protein_families-all-Reactions-filteredb100-fixed.tbl"
# strain_blast_hits_file = "matches_b100.tsv"

def format_number(num):
    if num.is_integer():
        return "{:.0f}".format(num)
    else:
        return str(num)

strain_blast_hits = pd.read_table(strain_blast_hits_file, header = None, 
                                  usecols=[0,1,11], names = ['Query', 'Target', 'Repr. bitscore'])
pan_reaction_table = pd.read_table(pan_reaction_table_file, comment = '#')

strain_reaction_table = pd.merge(strain_blast_hits, pan_reaction_table, 
                                 how = 'inner', left_on='Target', right_on='stitle')
strain_reaction_table.drop_duplicates(inplace = True, ignore_index = True,
                                      subset = ['Query', 'Target', 'rxn', 'pathway'])
strain_reaction_table.sort_values(by = "Repr. bitscore", ascending = False, inplace = True)
strain_reaction_table.drop_duplicates(subset = ['Query', 'pathway'], keep = "first", inplace = True)
strain_reaction_table.sort_index(inplace = True)
strain_reaction_table['stitle'] = strain_reaction_table['Query']
strain_reaction_table.drop(columns = ['Query', 'Target', 'Repr. bitscore'], inplace = True)
strain_reaction_table.to_csv(strain_reaction_table_file,
                              sep = '\t', index = False, na_rep = "NA", float_format = format_number)