import analysis
import clean_files as cf
import find_triads as ft
import os
import configparser

import util

config = configparser.ConfigParser()
config.read(os.path.join(os.pardir, 'config.ini'))
config = config['residue']

transpath = config['transformed_location']
output = config['output_location']

output_analysis = config['analysis_output_location']

if __name__ == "__main__":
    # check if files have been transformed
    if not os.path.isdir(transpath) or not os.listdir(transpath):
        cf.clean_files()

    # find triades and make csv files
    ft.find_triads()

    # results analysis
    util.create_folder(output_analysis)
    analysis.store_triad_count(output, output_analysis, similar=True)
