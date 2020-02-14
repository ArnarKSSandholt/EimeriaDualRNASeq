# Script that processes the read count data from htseq-count and splits the read counts by source as well
# as computing statistics about them, such as how lareg a fraction of the mapped reads map to the parasite
# Usage: python read_count_process.py /Path/to/metadata/file.tsv /Path/to/read/count/folder /Path/to/output/folder

import pandas as pd
import numpy as np
from os import listdir
import sys
import re

metadata_path = sys.argv[1]
read_count_path = sys.argv[2]
output_path = sys.argv[3]

read_count_filenames = listdir(metadata_path)
metadata_table = pd.read_table(metadata_path)

