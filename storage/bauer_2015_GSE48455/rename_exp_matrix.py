import sys

with open('probe_mapping.tsv', 'r') as fh:
    fh = fh.readlines()[1:]
    mappings = {}
    for line in fh:
        probe, gene = line.split('\t')
        mappings[probe.strip('\"')] = gene.strip().strip('\"')

with open('gene_expression_matrix', 'r') as fh:
    df = {}
    count = 1
    for line in fh:
        if count == 1:
            colnames = [i.strip('\"') for i in line.split('\t')]
        else:
            probe, exp = line.split('\t')[0].strip('\"'), line.split('\t')[1:]
            if probe in mappings.keys():
                if mappings[probe] not in "Nothing":
                    if mappings[probe] not in df.keys():
                        df[mappings[probe]] = exp

                    else:
                        for i in range(1, 10):
                            name = '{0}_{1}'.format(mappings[probe], i)
                            if name not in df.keys():
                                df[name] = exp
                                break
        count += 1
with open("gene_expression_matrix_mapped", "w") as fh:
    fh.write('\t{0}\n'.format('\t'.join(colnames).strip('"\n')))
    for name, exp in df.items():
        fh.write('{0}\t{1}'.format(name, '\t'.join(exp)))
