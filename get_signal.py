import pandas as pd
import urllib.request
import numpy as np
import requests, sys
import json
from collections import defaultdict




data = pd.read_csv("~/splicing_project/Data/output_data_non_coloc.csv")
data["DOMAIN"] = None
data["SIGNAL"] = None
data["TOPOLOGY"] = None
data["TRANSMEMBRANE"] = None
data["MOTIF"] = None
for ind, acc in enumerate(data["ALPHAFOLD.NAME"]):
    start = None
    print(acc)
    try:
        start_need, stop_need = data["ALIGN.COORDS"][ind].split("-")
    except AttributeError:
        start_need = None
    if type(acc) == str:
        requestURL = 'https://www.ebi.ac.uk/proteins/api/features/' + str(acc) + '?categories=DOMAINS_AND_SITES%2CTOPOLOGY%2CPTM'
        r = requests.get(requestURL, headers={ "Accept" : "application/json"})

        if not r.ok:
          #r.raise_for_status()
          #sys.exit()
          continue

        responseBody = r.text

        json_acceptable_string = responseBody.replace("'", "\"")
        try:
            d = json.loads(json_acceptable_string)
        except json.decoder.JSONDecodeError:
            continue
        #print(d)
        #print(d["features"][1]['ptms'])
        #print(d.keys())
        #print(d["features"])
        d = d["features"]
        signal = []
        domain = []
        topology = []
        motif = []
        annot = defaultdict(list)
        transmembrane = None
        for i in d:
            print(i)
            try:
                inter_need = set(np.arange(int(start_need), int(stop_need) + 1))
            except TypeError:
                continue
            try:
                inter = set(np.arange(int(i["begin"]), int(i["end"]) + 1))
            except ValueError:
                print(i)
            if inter.intersection(inter_need):
                try:
                    annot[i["type"]].append(i["description"])
                except KeyError:
                    annot[i["type"]].append(i["type"])
            if inter.intersection(inter_need) and i["category"] == "TOPOLOGY" and i["type"] == "TRANSMEM":
                    transmembrane = True

    for key, values in annot.items():
        try:
            try :
                data[key][ind] = ";".join(values)
            except TypeError :
                data[key][ind] = None
        except KeyError:
            data.insert(len(data.columns), key, None)
            try :
                data[key][ind] = ";".join(values)
            except TypeError :
                data[key][ind] = None

    data["TRANSMEMBRANE"][ind] = transmembrane


print(data)
data.to_csv("~/splicing_project/Data/output_data_non_coloc.csv", index=False)
"""
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
print(data["DOMAIN"])
data.to_csv("~/splicing_project/Data/output_data_non_sQTL_try.csv")
"""