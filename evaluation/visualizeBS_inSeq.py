import sys, os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import transforms

N_AA_PER_LINE=50

BS_THR=0.5

resultsFnamePairs= sys.argv[1]
receptorInstead=len(sys.argv)>2


lastElem= resultsFnamePairs.split(".")[-1]
replacePatter= "%s.gz" if lastElem=="gz" else "tab.%s"

resultsFnameL= resultsFnamePairs.replace(lastElem, replacePatter%"lig")
resultsFnameR= resultsFnamePairs.replace(lastElem, replacePatter%"rec")

pairScores= pd.read_table(resultsFnamePairs, comment="#", sep="\s+",
              dtype={"resIdL":str, "resIdR":str, "chainIdL":str, "chainIdR":str})


if receptorInstead:
  pairScores=pairScores.rename({"chainIdR": "chainIdL", "resIdR":"resIdL", "resNameR":"resNameL",
                                "chainIdL": "chainIdR", "resIdL": "resIdR", "resNameL": "resNameR",
                                }, axis=1)
  print( pairScores.head() )

partnerSeq= pairScores[["chainIdL", "resIdL", "resNameL"]].drop_duplicates()

partnerSeq= partnerSeq.rename( {"chainIdL": "chainId", "resIdL":"resId", "resNameL": "resName"}, axis=1)


partnerScores= pd.read_table(resultsFnameR if receptorInstead else resultsFnameL, comment="#", sep="\s+",
                         dtype={"resId": str, "chainId": str})

partnerScores= partnerScores.merge(partnerSeq )


resIds=partnerScores["resId"].map(lambda x: int(x) if x[-1].isdigit() else int(x[:-1])).sort_values()


partnerScores= partnerScores.iloc[resIds.index, :]


for chainId, df in partnerScores.groupby("chainId"):
  print("chainId: ", chainId)
  # print(df)
  seq= []
  curFragment=""
  preds=[]
  lastCateg=df.iloc[0,:]["categ"]
  for i, (aa, categ, pred) in df[["resName", "categ", "prediction"]].iterrows():
    if categ!= lastCateg:
      seq.append( (lastCateg, preds, curFragment) )
      curFragment=""
      preds = []
      lastCateg=categ

    curFragment+=aa
    preds+=[pred]

  seq.append( (categ, preds, curFragment ) )


  newSeq=[]
  curLine=0
  nAAInLine=0
  prevCateg= seq[0][0]
  for categ, preds, fragment in seq:
    fragLen= len(fragment)
    curFragment = ""
    curPreds= []
    while fragLen>0:
      n2Add = min(fragLen, max(0, N_AA_PER_LINE - nAAInLine))

      curFragment += fragment[:n2Add]
      curPreds += preds[:n2Add]
      # print(curLine, fragment, curFragment, n2Add, fragLen, nAAInLine); raw_input("enter")

      fragment= fragment[n2Add: ]
      preds= preds[n2Add: ]

      fragLen -= n2Add
      nAAInLine += n2Add
      newSeq.append( (categ, curPreds, curLine, curFragment) )

      if nAAInLine>=N_AA_PER_LINE:
        curLine += 1
        curFragment = ""
        curPreds=[]
        nAAInLine = 0


  yposition= 0.9
  deltaY= 0.1
  xposition=0.02
  fontSize=10
  prevLine= newSeq[0][1]

  t = plt.gca().transData
  fig = plt.gcf()

  for categ, preds, line, fragment in newSeq:
    if (line!=prevLine):
      t = plt.gca().transData
      prevLine= line
    text= plt.text(xposition, yposition-line*deltaY, fragment, dict(size=fontSize, fontname='Courier New'), color="red" if categ==1 else "black", transform=t)

    text.draw(fig.canvas.get_renderer())
    ex = text.get_window_extent()
    t = transforms.offset_copy(text._transform, x=ex.width, units='dots')


  prevLine= newSeq[0][1]
  for categ, preds, line, fragment in newSeq:
    if (line!=prevLine):
      t = plt.gca().transData
      prevLine= line

    thr_preds= "".join([" " if pred<BS_THR else "*" for pred in preds])
    text= plt.text(xposition, yposition-line*deltaY - (deltaY*0.5), thr_preds,
                   dict(size=fontSize, fontname='Courier New'), color="green" , transform=t)

    text.draw(fig.canvas.get_renderer())
    ex = text.get_window_extent()
    t = transforms.offset_copy(text._transform, x=ex.width, units='dots')


  plt.show()