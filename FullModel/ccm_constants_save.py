"""
This file saves all the biological constants of the Chlamy CCM 
to a file in /data_v5 called ccm_constants.csv .
"""
from ccm_constants import *
import csv

_ofile = csv.writer(open("data_v5/ccm_constants.csv", "w"))
_fields = []
_vals = []

for key in dir():

    if key[0] == "_" or key == "csv":
        continue

    _fields.append(key)
    _vals.append(globals()[key])

_ofile.writerow(_fields)
_ofile.writerow(_vals)
