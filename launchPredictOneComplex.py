import os

ARE_FOR_TRAIN_AND_TEST=True
NORMALIZE_SCORES_IN_BFACTOR=True

def launcPredictOneComplex(workingDir=None,  lPdbId=None, lFname=None, lSequence=None,
                           rPdbId=None, rFname=None, rSequence=None, isHomoComplex=False, areForTrainAndTest=ARE_FOR_TRAIN_AND_TEST):

  from predictComplexes import predictOneComplex
  from computeFeatures.featuresComputer import getPrefix

  def myGetPrefix(oneStr):
    if oneStr is None:
      return oneStr
    else:
      return  getPrefix(oneStr)

  lPrefix, rPrefix=None, None
  if not lSequence:
    lPrefix= myGetPrefix(lFname) if lFname else myGetPrefix(lPdbId)

  if not rSequence:
    rPrefix= myGetPrefix(rFname) if rFname else myGetPrefix(rPdbId)

  if lPrefix and rPrefix:
    if lPrefix != rPrefix:
      areForTrainAndTest=False
    else:
      areForTrainAndTest=True
  else:
    areForTrainAndTest = False

  if areForTrainAndTest: print("WARNING: inputs are for test")
  print(isHomoComplex)

  resultsDir, (p, l, r) = predictOneComplex(allRootDir=workingDir,
                                         lPdbId=lPdbId, lFname=lFname, lSequence=lSequence,
                                         rPdbId=rPdbId, rFname=rFname, rSequence=rSequence, isHomoComplex=isHomoComplex,
                                         areForTrainAndTest=areForTrainAndTest, removeInputs=False)


  from pythonTools.insertPredsInBfactor import insertPredsInBfactor

  insertPredsInBfactor(lFname, os.path.join(resultsDir, "PRED.res.tab.lig.gz"), os.path.join(resultsDir, os.path.basename(lFname) ) , min_max_norm=NORMALIZE_SCORES_IN_BFACTOR, verbose=False)

  if rFname:
    insertPredsInBfactor(rFname, os.path.join(resultsDir, "PRED.res.tab.rec.gz"), os.path.join(resultsDir, os.path.basename(rFname) ), min_max_norm=NORMALIZE_SCORES_IN_BFACTOR, verbose=False )


if __name__=="__main__":
  import argparse
  parser = argparse.ArgumentParser(description='Predicts binding site from input files, ids or sequences using the model specified in config file')
  parser.add_argument('workingDir', type=str,  help='path where computations will be saved')
  parser.add_argument('--lFname', type=str, nargs='?',  help='path to ligand protein filename', required=False)
  parser.add_argument('--lSequence', type=str, nargs='?',    help='sequence of ligand protein', required=False)
  parser.add_argument('--lPdbId', type=str, nargs='?',  help='pdb id of ligand protein', required=False)


  parser.add_argument('--rFname', type=str, nargs='?',  help='path to receptor protein filename. Skip receptor if predicting homocomplex', required=False)
  parser.add_argument('--rSequence', type=str, nargs='?',    help='sequence of receptor protein. Skip receptor if predicting homocomplex', required=False)
  parser.add_argument('--rPdbId', type=str, nargs='?',  help='pdb id of receptor protein. Skip receptor if predicting homocomplex', required=False)

  args = vars(parser.parse_args())

  nLArgs= len([ varName for varName in args if args[varName] is not None and varName.startswith("l")] )
  nRArgs= len([ varName for varName in args if args[varName] is not None and varName.startswith("r")] )

  assert nLArgs<=1, "Error, only one filename, sequence or PDB ID must be provided for ligand partner."
  assert nLArgs==1, "Error you always needs to provide --lFname or --rSequence or  --lPdbId"

  assert nRArgs<=1, "Error, only one filename, sequence or PDB ID must be provided for receptor partner."

  if nLArgs== 0 or nRArgs==0:
    args["isHomoComplex"]=True
  elif nLArgs== 0 and nRArgs==0:
    raise AssertionError("Error, at least on partner must be provided")
  else:
    args["isHomoComplex"]=False

  print(args)
  launcPredictOneComplex( **args)