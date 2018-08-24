import re


with open("meta-tmp.tsv", "r") as fh:
    fh = fh.readlines()[1:]
    meta = {}
    for line in fh:
        geo = re.findall("GSM[0-9]+", line)[0]
        treatment = re.findall("[a-z]+[-]*treated", line)[0]
        timen = re.findall("[0-9][DW]", line)[0][0]
        timel = re.findall("[0-9][DW]", line)[0][1]
        meta[geo] = [treatment, int(timen), timel]
counter = 1
for key, value in meta.items():
    n = 0
    if value[2] == "W":
        n = 7
    elif value[2] == "D":
        n = 1
    else:
        print(value[2], "WHAT THE FUCK")
        break
    meta[key] = [value[0], n*value[1]]
    if meta[key][0] == "bleomycin-treated":
        meta[key].extend([0, 1])
    elif meta[key][0] == "vehicle-treated":
        meta[key].extend([1, 0])
    elif meta[key][0] == "untreated":
        meta[key].extend([1, 1])
    else:
        print("WTF")
        break


def order_my_dict(mytuple):
    untreated = []
    for i in mytuple[1]:
        if i == 'untreated':
            untreated.append(0)
        else:
            untreated.append(1)
    return untreated, mytuple[1][1], mytuple[1][3], mytuple[1][2]


ordered_meta = dict()
counter = 0
track = []
for k, v in sorted(meta.items(), key=order_my_dict):
    ordered_meta[k] = v
    if v == track:
        ordered_meta[k].append(counter)
    else:
        counter += 1
        ordered_meta[k].append(counter)
    track = v[:-1]

with open("metadata.tsv", "w") as fh:
    fh.write("Geo\tTreatment\tTime\tControl\tBleomycin\tReplicate\n")
    for key, value in ordered_meta.items():
        fh.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(key, value[0],
                                                         value[1], value[2],
                                                         value[3], value[4]))
