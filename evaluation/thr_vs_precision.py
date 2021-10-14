import sys, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.isotonic import IsotonicRegression
import json
import joblib


'''
input is a json {"thr": [thr0, thr1, ...], "prec":[precision0, precision1, ...] }

json is computed using evaluation/getPrecsAt_diffThr.py

python evaluation/thr_vs_precision.py ~/Tesis/rriPredMethod/data/bench5Data/newCodeData/results/thr_2_pretc_b5_seq.json

'''
DO_PLOT=True
if len(sys.argv)<2:
  raise ValueError("Error, 1 argument should be provided")
else:
  if len(sys.argv)==3:
    save_isotonic_modelFname=  os.path.expanduser( sys.argv[2])
  else:
    save_isotonic_modelFname= None
  jsonFname= os.path.expanduser( sys.argv[1])

with open(jsonFname) as f:
  data= json.load(f)
thr= data["thr"]
prec= data["prec"]


def logFit(thr, prec, maxThr= 4):
  thr= np.array(thr)
  prec= np.array(prec)
  prec= prec[ thr<=maxThr]
  thr= thr[ thr<=maxThr]
  a_b= np.polyfit( np.log(thr), prec, 1)
  formulaStr="prec= %f log(thr) + %f"%tuple(a_b)
  return lambda thr: a_b[0]*np.log(thr)+a_b[1], formulaStr

def logisticFit(thr, prec, minThr=0.01, maxThr= 4):
  thr= np.array(thr)
  prec= np.array(prec)
  prec= prec[ thr<=maxThr]
  thr= thr[ thr<=maxThr]

  prec= prec[ thr>minThr]
  thr= thr[ thr>minThr]
  
  objFun= lambda x, alpha, beta: 1./ (1+ np.exp(np.multiply(alpha,x)+ beta))  
  popt, pcov = curve_fit(objFun, thr, prec)
  print(popt)
  print( np.sqrt(np.diag(pcov)) )
#  formulaStr="prec= 1/( 1 + exp(%f thr + %f) )"%tuple(popt)
  formulaStr="sigmoid"
  return lambda x: objFun(x, *popt ), formulaStr  


def polyNomialFit(thr, prec, maxThr= 4):
  thr= np.array(thr)
  prec= np.array(prec)
  prec= prec[ thr<=maxThr]
  thr= thr[ thr<=maxThr]
  objFun= lambda thr, alpha, beta: alpha*thr**beta
  popt, pcov = curve_fit(objFun, thr, prec)
  print(popt)
  print( np.sqrt(np.diag(pcov)) )
  formulaStr=str(tuple(popt))
  return lambda thr: objFun(thr, *popt ), formulaStr
  
def isotonicFit(thr, prec, maxThr= 999):
  thr= np.array(thr)
  prec= np.array(prec)
  prec= prec[ thr<=maxThr]
  thr= thr[ thr<=maxThr]
  objFun= lambda thr, alpha, beta: alpha*thr**beta
  
  isoReg= IsotonicRegression(y_min=0, y_max=1)
  isoReg.fit(thr, prec)
  if save_isotonic_modelFname:
    joblib.dump(isoReg, save_isotonic_modelFname)
  return lambda x: isoReg.predict(x) , "isotonic"
  
#fitFun, formulaStr= logFit(thr, prec)
fitFun_1, formulaStr_1= logisticFit(thr, prec)
#fitFun, formulaStr= polyNomialFit(thr, prec)
fitFun_2, formulaStr_2= isotonicFit(thr, prec)

scores= np.linspace(0,3,200)
precision= fitFun_2(scores)
#for s, p in zip(scores, precision):
#  if abs(p-0.5)<0.01:
#    print(s,p)
#  if abs(p-0.6)<0.01:
#    print(s,p)
#  if abs(p-0.7)<0.01:
#    print(s,p)  
#raw_input()

'''
prec score
0.5  0.18090
0.6  0.58794
0.7  2.39698
'''

if DO_PLOT:
  fig, ax = plt.subplots()

  ax.plot(thr, prec, 'c.', label='precision data')
  ax.plot(thr,  fitFun_1(thr), 'b--', label=formulaStr_1+" fit")
  ax.plot(thr,  fitFun_2(thr), 'r--', label=formulaStr_2+" fit")
  ax.set_xlabel("score threshold")
  ax.set_ylabel("precision")
  ax.set_title("DBv5 scores seq model")
  legend = ax.legend(loc='lower right', shadow=True)

  #plt.savefig("/home/rsanchez/Tesis/papers/myPapers/2018/rriInteraction/major1Images/scores_vs_precision_bothFit_seq.png",
  #            format="png", dpi=350)
  plt.show()

