#import __main__
#__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

import sys, time, os
import pymol
import pandas as pd


data = pd.read_csv("~/splicing_project/Data/output_data_non_coloc_try.csv")

#pymol.finish_launching()

##
# Read User Input
#spath = os.path.abspath(sys.argv[1])
#sname = spath.split('/')[-1].split('.')[0]

# Load Structures
print(data.shape[0])
for i in range(data.shape[0]):
    try:
        print(i, data["ALPHAFOLD NAME"][i])
        pymol.cmd.load("~/splicing_project/UP000005640_9606_HUMAN/AF-"+data["ALPHAFOLD NAME"][i]+"-F1-model_v1.pdb", data["ALPHAFOLD NAME"][i])
        pymol.cmd.disable("all")
        pymol.cmd.enable(data["ALPHAFOLD NAME"][i])
        pymol.cmd.color("orange")
        pymol.cmd.color('green',"resi "+str(data["ALIGN COORDS"][i]))
        pymol.cmd.ray(1200, 1200)
        pymol.cmd.png("~/splicing_project/Data/visuals/proteins/non_coloc/"+ data["ALPHAFOLD NAME"][i]+".png")
    except pymol.CmdException:
        print("HERE")
        continue

    # Get out!
pymol.cmd.quit()