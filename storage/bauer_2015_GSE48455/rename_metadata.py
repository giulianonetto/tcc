with open("metadata.tsv", "r") as fh, open("names.tsv", "r") as fh2:
    content = fh.readlines()
    header = content[0]
    metadata = {}
    for line in content[1:]:
        line = line.strip().split("\t")
        gsm, data = line[0], line[1:]
        metadata[gsm] = data

    mapping = {}
    for entry in fh2:
        tag, name = entry.split("\t")[0].strip('"\n'), entry.split("\t")[1].strip('"\n')
        mapping[name] = tag
with open("metadata_renamed.tsv", "w") as fh:
    fh.write(header)
    for gsm in metadata.keys():
        fh.write('{0}\t{1}\n'.format(mapping[gsm], '\t'.join(metadata[gsm])))