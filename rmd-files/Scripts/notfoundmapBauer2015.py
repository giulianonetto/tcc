with open("/home/giulianonetto/windows/tcc/rmd-files/data/GSE48455/suppdata/not_found.gmt", "r") as fh:
    not_found = {}
    for line in fh:
        line = line.split('\t')
        not_found[line[0]] = line[1:]

with open("/home/giulianonetto/windows/tcc/storage/gene_info/Aliases/aliases_rats.tsv", "r") as fh:
    aliases = []
    for line in fh:
        line = line.strip('\n').split('\t')
        aliases.append(line)

matches = {}
for module, genes in not_found.items():
    for gene in genes:
        for entry in aliases:
            if gene in entry:
                if gene not in matches:
                    matches[gene] = entry
                else:
                    print(gene, "---", entry)

with open("/home/giulianonetto/windows/tcc/rmd-files/data/GSE48455/suppdata/not_found_map.gmt", "w") as fh:
    for gene, match in matches.items():
        match_as_string = ' '.join(match)
        fh.write(f"{gene}\t{match_as_string}\n")