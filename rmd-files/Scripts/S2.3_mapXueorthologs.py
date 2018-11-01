orthologs = {}
with open("/home/giulianonetto/windows/tcc/checklist/genesigs/human_mouse.map", "r") as fh:
	fh = fh.readlines()[1:]
	for line in fh:
		human = line.split("\t")[0]
		rat = line.split("\t")[1].strip('\n').upper()
		if rat != "NA":
			orthologs[human] = rat


geneSets = {}
with open("data/XueModules/genesets.gmt", "r") as fh:
	for line in fh:
		line = line.strip('\n').split('\t')
		module = '\t'.join([line[0], line[1]])
		genes = []
		for i in line[2:]:
			genes.append(i)
		geneSets[module] = genes


with open('data/XueModules/genesets.ortho.gmt', 'w') as fh:
	for module, genes in geneSets.items():
		mapped = []
		for gene in genes:
			if gene in orthologs:
				mapped.append(orthologs[gene])
			else:
				mapped.append(gene)
		mapped_string = '\t'.join(mapped)
		fh.write(f"{module}\t{mapped_string}\n")
