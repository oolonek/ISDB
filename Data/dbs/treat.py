from __future__ import print_function

import csv
import sys

# Solving an issue on some machines (32 bits?) 
maxInt = sys.maxsize
csv.field_size_limit(maxInt)

if len(sys.argv) != 5:
    print(" Incorrect number of arguments")
    sys.exit()

# The tremolo file
result_file = sys.argv[1]
cytoscape_file = sys.argv[2]
db_file = sys.argv[3]
output_file = sys.argv[4]

# Lazy, we load everything in memory
print(" Treating the SMILES list: ", end='')
with open(db_file, 'rt') as f:
    reader = csv.reader(f)
    header = True
    smiles_store = {}
    cn_store = {}
    for row in reader:
        if header:
            unpd_id_pos = row.index('UNPD_ID')
            smiles_pos = row.index('SMILES')
            cn_pos = row.index('cn')
            header = False
        else:
            smiles_store[row[unpd_id_pos]] = row[smiles_pos]
            cn_store[row[unpd_id_pos]] = row[cn_pos]

print("treated {} compounds".format(len(smiles_store)))

print(" Treating the result file {}".format(result_file))

with open(result_file, 'rt') as f:
    reader = csv.reader(f, delimiter='\t')
    header = False
    result = {}
    result_id = {}
    result_cn= {}
    for row in reader:
        if not header:
            scan_id_pos = row.index('#Scan#')
            unpd_id_pos = row.index('CompoundName')
            header = row
        else:
            if row[scan_id_pos] in result:
                # Handle multiple molecules by node
                result[row[scan_id_pos]] += "," + smiles_store[row[unpd_id_pos]]
                result_cn[row[scan_id_pos]] += "," + cn_store[row[unpd_id_pos]]
                result_id[row[scan_id_pos]] += "," + row[unpd_id_pos]
            else:
                result[row[scan_id_pos]] = smiles_store[row[unpd_id_pos]]
                result_cn[row[scan_id_pos]] = cn_store[row[unpd_id_pos]]
                result_id[row[scan_id_pos]] = row[unpd_id_pos]

print(" Merging to the cytoscape result")

with open(cytoscape_file, 'rt') as f:
    reader = csv.reader(f, delimiter='\t')
    header = False
    output = []
    for row in reader:
        if not header:
            cluster_index = row.index('cluster index')
            header = row + ["SMILES"] + ["UNPD_IDs"] + ["ChemicalNames"]
        else:
            temp = row
            if row[cluster_index] in result:
                temp += [result[row[cluster_index]]]
            else:
                temp += []

            if row[cluster_index] in result_id:
                temp += [result_id[row[cluster_index]]]
            else:
                temp += []

            if row[cluster_index] in result_cn:
                temp += [result_cn[row[cluster_index]]]
            else:
                temp += []
            
            output += [temp]


print("Outputing")
with open(output_file, 'w') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t',
                        quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(header)
    for row in output:
        writer.writerow(row)
