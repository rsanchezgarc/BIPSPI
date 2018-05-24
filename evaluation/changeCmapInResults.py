import os, sys
import pandas as pd

RESULTS_PATH="/home/rsanchez/Tesis/rriPredMethod/data/sppider/results/struct_2"
ALTERNATIVE_CMAPS_PATH="/home/rsanchez/Tesis/rriPredMethod/data/sppider/computedFeatures/common/contactMaps_homoFixed"
NEW_RESULTS_PATH="/home/rsanchez/Tesis/rriPredMethod/data/sppider/results/struct_2_homoFixed"

def fixOneFname(fname, inpath= RESULTS_PATH, cmapspath=ALTERNATIVE_CMAPS_PATH, outpath=NEW_RESULTS_PATH):
  prefix= fname.split(".")[0]
  full_fname_in= os.path.join(inpath, fname)
  full_fname_out= os.path.join(outpath, fname)
  data= pd.read_table(full_fname_in, sep='\s+', comment="#", dtype={"structResIdL":str, "structResIdR":str,
                                                            "chainIdL":str, "chainIdR":str})
  full_fname_cmap= os.path.join(cmapspath, prefix+".cMap.tab")
  cmap= pd.read_table(full_fname_cmap, sep='\s+', comment="#", dtype={"structResIdL":str, "structResIdR":str,
                                                            "chainIdL":str, "chainIdR":str})
  data_fixed= data.merge(cmap, on=["chainIdL", "structResIdL", "resNameL", "chainIdR", "structResIdR", "resNameR"])
#  print(data_fixed["categ_x"].equals(data_fixed["categ_y"]))
  data_fixed["categ_x"]= data_fixed["categ_y"]
  data_fixed= data_fixed.drop(["categ_y"], axis=1)
  data_fixed.rename(index=str, columns={"categ_x": "categ"} , inplace=True)
  data_fixed.to_csv(full_fname_out, index=False, sep=" ")

def main():
  for fname in os.listdir(RESULTS_PATH):
    if fname.endswith(".tab"):
      print(fname)
      fixOneFname(fname)
#      raw_input("press enter for next pdb result")
    
if __name__=="__main__":
  main()
  
