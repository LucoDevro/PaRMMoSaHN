## TODO: replace Linux paths to more general ones using Pathlib

import cobra as cb
import pandas as pd
import json
import os
import re
import sys

root = os.getcwd()
report_dir = sys.argv[1]
model_dir = sys.argv[2]
output_dir = sys.argv[3]

#report_dir = '/mnt/DATA/PhD/WPs/WP1/clostridia_models/05-automatic_curations/memote'
#model_dir = '/mnt/DATA/PhD/WPs/WP1/clostridia_models/05-automatic_curations/00-input_models'
#output_dir = '/mnt/DATA/PhD/WPs/WP1/clostridia_models/05-automatic_curations/01-duplicates_imbalances'

def read_test_results(test):
    rxns = set()
    first_encountered = {}
    for model in filter(lambda x: '.json' in x, os.listdir(report_dir)):
        with open(os.path.join(report_dir, model), 'r') as handle:
            cont = json.load(handle)    
        encountered_rxns = cont['tests'][test]['data']
        try:
            if isinstance(encountered_rxns[0], list):
                encountered_rxns = [tuple(i) for i in encountered_rxns]
        except IndexError:
            continue
        new_rxns = set(encountered_rxns).difference(rxns)
        for r in new_rxns:
            first_encountered[r] = model.split('.')[0] # removes extension
        rxns = rxns.union(new_rxns)
    return first_encountered
    
first_encounter_rxns_charge_imbalance = read_test_results('test_reaction_charge_balance')
first_encounter_rxns_mass_imbalance = read_test_results('test_reaction_mass_balance')
first_encounter_dupl_rxns = read_test_results('test_find_duplicate_reactions')

def report_problem_balances(incidence_dict, save = False, save_name = None):
    df = []
    for r,m in incidence_dict.items():
        network = cb.io.read_sbml_model(os.path.join(model_dir, m + '.xml'))
        rxn_obj = network.reactions.get_by_id(r)
        new_record = {'ID': rxn_obj.id,
                      'Name': rxn_obj.name,
                      'Equation': rxn_obj.build_reaction_string(),
                      'Imbalance': str(rxn_obj.check_mass_balance()).strip('{}'),
                      'First_encounter': m}
        df.append(new_record)
    df = pd.DataFrame.from_dict(df)
    if save:
        df.to_csv(os.path.join(output_dir, save_name), sep = "\t", header = True, index = False)
    return df

report_problem_balances(first_encounter_rxns_charge_imbalance, True, "charge_imbalances.tsv")
report_problem_balances(first_encounter_rxns_mass_imbalance, True, "mass_imbalances.tsv")

def report_duplicate_reactions(incidence_dict, save = False, save_name = None):
    df = []
    for rp,m in incidence_dict.items():
        network = cb.io.read_sbml_model(os.path.join(model_dir, m + '.xml'))
        rxn_pair_obj = [network.reactions.get_by_id(r) for r in rp]
        id_regex = re.compile("[0-9]{5}")
        rxn_ids = [int(id_regex.findall(r)[0]) for r in rp]
        min_idx = rxn_ids.index(min(rxn_ids))
        default_rxn = rxn_pair_obj.pop(min_idx)
        non_default_rxn = rxn_pair_obj[0]
        new_record = {'First_encounter': m,
                      'Pair': str(rp).strip('()'),
                      'ID 1': default_rxn.id,
                      'Name 1': default_rxn.name,
                      'Equation 1': default_rxn.build_reaction_string(),
                      'Lower bound 1': default_rxn.lower_bound,
                      'Upper bound 1': default_rxn.upper_bound,
                      'GPR 1': default_rxn.gene_reaction_rule,
                      'ID 2': non_default_rxn.id,
                      'Name 2': non_default_rxn.name,
                      'Equation 2': non_default_rxn.build_reaction_string(),
                      'Lower bound 2': non_default_rxn.lower_bound,
                      'Upper bound 2': non_default_rxn.upper_bound,
                      'GPR 2': non_default_rxn.gene_reaction_rule}
        df.append(new_record)
    df = pd.DataFrame.from_dict(df)
    if save:
        df.to_csv(os.path.join(output_dir, save_name), sep = "\t", header = True, index = False)
    return df

report_duplicate_reactions(first_encounter_dupl_rxns, True, "duplicate_reactions.tsv")

