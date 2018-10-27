import sys

with open('/home/giulianonetto/windows/tcc/rmd-files/data/GSE48455/Bauer2015_features.tsv', 'r') as fh:
    fh = fh.readlines()[1:]
    features = {}
    for line in fh:
        line = line.split('\t')
        probe, gene = line[0], line[3]
        if len(gene) == 0:
            gene = "empty"
        features[probe] = gene

with open('/home/giulianonetto/windows/tcc/rmd-files/data/GSE48455/Bauer2015_vsn_cDNA_all.tsv', 'r') as fh:
    df = {}
    count = 1
    for line in fh:
        if count == 1:
            colnames = line
        else:
            line = line.split('\t')
            probe, exp = line[0], line[1:]
            probe = probe.split(".")[0]
            if probe in features.keys():
                if features[probe] not in df.keys():
                    df[features[probe]] = exp
                else:
                    for i in range(1, 100):
                        name = '{0}_{1}'.format(features[probe], i)
                        if name not in df.keys():
                            df[name] = exp
                            break
        count += 1
with open("/home/giulianonetto/windows/tcc/rmd-files/data/GSE48455/Bauer2015_vsn_cDNA_all_mapped.tsv", "w") as fh:
    fh.write('{0}'.format(colnames))
    for name, exp in df.items():
        fh.write('{0}\t{1}'.format(name, '\t'.join(exp)))
