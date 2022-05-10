import pandas as pd
import urllib.request
import numpy as np

data = pd.read_csv("~/splicing_project/Data/output_non_sQTL_full.csv")
start = None
data["DOMAIN"] = None
data["MOTIF"] = None
print(data.columns)
for ind, acc in enumerate(data["ALPHAFOLD.NAME"]):
    print(acc)
    try:
        start_need, stop_need = data["ALIGN.COORDS"][ind].split("-")
    except AttributeError:
        start_need = None
    url = "https://www.uniprot.org/uniprot/?query=" + str(acc) + "&columns=feature(REGION)&format=tab&compress=no"
    with urllib.request.urlopen(url) as r:
        fasta = r.read().decode('utf-8').strip()
    fasta = fasta.split("\n")
    for i, line in enumerate(fasta):
        spl = line.split("/")
        print(spl)
        domain = None
        for j, st in enumerate(spl):
            if "REGION" in st:
                l = st.split(" ").index("REGION")
                coord = st.split(" ")[l+1]
                #print(coord)
                try:
                    start, stop = coord.split("..")
                    stop = stop.replace(";", '')
                except ValueError:
                    stop = coord.split("..")[0]
                    stop = stop.replace(";", '')
                    start = stop
                #print(start, stop, start_need, stop_need)
                if start and start_need:
                    try:
                        inter = set(np.arange(int(start), int(stop)+1))
                        inter_need = set(np.arange(int(start_need), int(stop_need) + 1))
                        if inter.intersection(inter_need) :
                            print(j, spl[j + 1])
                            if not domain or len(spl[j + 1].split('=')[1].replace('"', '')) != "Disordered;":
                                domain = spl[j + 1].split('=')[1].replace('"', '')
                                domain = domain.replace(";", '')
                        start = None
                    except ValueError:
                        continue
                        #print(spl[j+1].split('=')[1])
        if domain:
            data["DOMAIN"][ind] = domain
for i, dom in enumerate(data["DOMAIN"]):
    data["DOMAIN"][i] = ";".join(dom)

data.to_csv("~/splicing_project/Data/output_non_sQTL_full.csv")
