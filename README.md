# Searching for triades in enzymes

### Project Description
This program aims to find triades, describe them and store in csv format.
Triades consist of:
* NUC - OG atom in SER or SG atom in CYS 
* BASE - CG atom in HIS, ASP or GLU
* ACID - OD1 or OD2 atom

Structures are analysed and information is collected using Biopython PDBParser.

Triades must fit given angle and side length criteria.
There are two kinds of triades:
* those that fit both criteria - accepted triades
* thos that only fit angle criteria - similar triangles


### Project structure and usage
##### config.ini
The `config.ini` file stores the locations of folders used in running this program.
You should match the paths in `config.ini` to the paths in your computer where you want to store those items.
* `location` - the location of PDB files 
* `transformed_location` - the location of transformed PDB files (after running `clean_files.py`)
* `output_location` - csv outputs
* `keywords` - keywords with triad candidates used to clear out files

##### main.py
This file automatically runs the entire program.
* If files haven't been transformed yet for search optimization by weeding out all atom rows that cannot be triad candidates, `clear_files()` does so.
* `find_triades` makes the search and stores output to 44 csv files (22 for accepted and 22 for similar).
