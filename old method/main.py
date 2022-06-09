import clean_files as cf
import find_triades as ft
import os
import configparser

config = configparser.ConfigParser()
config.read(os.path.join(os.pardir, 'config.ini'))
config = config['default']

transpath = config['transformed_location']

if __name__ == "__main__":
    # check if files have been transformed
    if not os.path.isdir(transpath) or not os.listdir(transpath):
        cf.clean_files()

    # find triades and make csv files
    ft.find_triades()