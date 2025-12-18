#! /bin/python

from os import listdir

from cobra.io import read_sbml_model, save_json_model, save_matlab_model

models = list(filter(lambda x: ".xml" in x, listdir()))
for m in models:
    basename = m[:-4]

    # Converting models
    print("Converting " + basename)
    model = read_sbml_model(basename + ".xml")
    save_matlab_model(model, basename + ".mat")
    save_json_model(model, basename + ".json")
