'''
Used to compare to other results pair predictions it they have different contact map criterium
'''

import os, sys
import pandas as pd


def fixOneFname(fname, inpath, cmapspath, outpath):

  prefix= fname.split(".")[0]
  full_fname_in= os.path.join(inpath, fname)
  full_fname_out= os.path.join(outpath, fname)

  data= pd.read_table(full_fname_in, sep='\s+', comment="#", dtype={"resIdL":str, "resIdR":str,
                                                            "chainIdL":str, "chainIdR":str})

  data.rename(index=str, columns={"structResIdL": "resIdL", "structResIdR": "resIdR"} , inplace=True)
  data= data.astype({colname: str for colname in data.columns[:6]})
  full_fname_cmap= os.path.join(cmapspath, prefix+"_.cMap.gz")

  cmap= pd.read_table(full_fname_cmap, sep='\s+', comment="#", dtype={"resIdL":str, "resIdR":str,
                                                            "chainIdL":str, "chainIdR":str})

  data_fixed= data.merge(cmap, on=["chainIdL", "resIdL", "resNameL", "chainIdR", "resIdR", "resNameR"])

  if "categ_x" in data_fixed:
    data_fixed["categ_x"]= data_fixed["categ_y"]
    data_fixed= data_fixed.drop(["categ_y"], axis=1)
    data_fixed.rename(index=str, columns={"categ_x": "categ"} , inplace=True)
  else:
    data_fixed= data_fixed[ ["chainIdL", "resIdL", "resNameL", "chainIdR", "resIdR", "resNameR", "categ", "prediction"] ]

  data_fixed.to_csv(full_fname_out, index=False, sep=" ", compression="gzip" if fname.endswith("gz") else None )

def main(inPath, cmapPath, resultsPath):
  for fname in os.listdir(inPath):
    if fname.endswith(".tab") or fname.endswith(".tab.gz"):
      print(fname)
      fixOneFname(fname, inPath, cmapPath, resultsPath)
    
if __name__=="__main__":
  inPath, cmapPath, resultsPath = sys.argv[1:]
  main( inPath, cmapPath, resultsPath )
  
