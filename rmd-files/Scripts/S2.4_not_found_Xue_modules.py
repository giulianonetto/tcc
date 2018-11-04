import re

with open("/home/giulianonetto/windows/tcc/rmd-files/data/GSE48455/suppdata/not_found.gmt", "r") as fh:
    not_found = []
    for line in fh:
        line = line.split('\t')
        for gene in line[1:]:
            if gene != "NA":
                not_found.append(gene)

with open("/home/giulianonetto/windows/tcc/storage/gene_info/Aliases/aliases_humans_and_rats.tsv", "r") as fh:
    aliases = []
    for line in fh:
        names = re.split("\t| ", line.strip("\n"))
        aliases.append(names)

matches = {}
for gene in not_found:
    for entry in aliases:
        if gene in entry:
            if gene not in matches.keys():
                matches[gene] = entry
            else:
                if matches[gene] != entry:
                    for i in entry:
                        if i not in matches[gene]:
                            matches[gene].append(i)

with open("/home/giulianonetto/windows/tcc/rmd-files/data/GSE48455/suppdata/not_found_map.gmt", "w") as fh:
    for gene, alias_list in matches.items():
        alias_list = ' '.join(alias_list)
        fh.write(f'{gene} {alias_list}\n')
