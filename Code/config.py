# tag name for the dataset of interest
tag='Mpop'

# file extention of the MuSiC files
ext='pop' 

# absolute path of the file containing the list of PDBs for which you want to generate the data
lst='/path/to/Main/List/Mpop_list.txt' 

# absolute path of CLUSTALW2 for the alignment step
clustalw2 = '/path/to/clustalw2'

# absolute path of the base directory
import os
base_dir = os.path.abspath(os.getcwd()).strip('Code') 

############################################################################

# Main directories
data_dir = base_dir + 'Data/'
list_dir = base_dir + 'List/'
code_dir = base_dir + 'Code/'

# Data directories
mus_dir = data_dir + ext + '/'
aas_dir = data_dir + 'aa_seq/'
nas_dir = data_dir + 'na_seq/'
ali_dir = data_dir + 'ali/'
mut_dir = data_dir + 'mut/'
mux_dir = data_dir + ext + 'x/'

# List paths
lst_pth = lst
map_pth = list_dir + tag + '_seqmap.txt'
