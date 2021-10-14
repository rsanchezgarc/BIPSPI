import sys, os
import joblib
from codifyComplexes.ComplexCodified import ComplexCodified, ComplexSeqStructCodified #required for joblib.load
from joblib import Parallel, delayed
from trainAndTest.crossValidationSplitter import N_SUB_FOLDS

N_JOBS=8

def main(rootPath):

  Parallel(n_jobs=N_JOBS)( delayed(checkOnePrefix)(fname, rootPath) for fname in
                                                        os.listdir( os.path.join(rootPath, "sampledInputs")) )


def removePrefix(prefix, rootPath):
  fnameSampled = os.path.join(rootPath, "sampledInputs", prefix + ".train.pkl.gz")
  fnameFull = os.path.join(rootPath, "allInputs", prefix + ".predict.pkl.gz")
  print("Removing: %s\n%s" % (fnameSampled, fnameFull))
  try:
    os.remove(fnameSampled)
  except OSError:
    pass
  try:
    os.remove(fnameFull)
  except OSError:
    pass

def getRelatedPrefixes(prefix):
  if "@" in prefix:
    raw_prefix= prefix.split("@")[0]
    return [ raw_prefix+"@%0.2d"%i for i in range(N_SUB_FOLDS)]
  else:
    return []

def errorLoading(prefix, rootPath):
  fnameSampled = os.path.join(rootPath, "sampledInputs", prefix + ".train.pkl.gz")
  fnameFull = os.path.join(rootPath, "allInputs", prefix + ".predict.pkl.gz")
  try:
    joblib.load(fnameSampled)
    joblib.load(fnameFull)
    return False
  except (ValueError, EOFError, KeyError, TypeError, IOError) as e:
    print("Error with: %s\n%s"%(fnameSampled, fnameFull))
    print(e)
    return True

def checkOnePrefix(baseName, rootPath):
  prefix = baseName.split(".")[0]
  if errorLoading(prefix, rootPath):
    for prefixToRemove in getRelatedPrefixes(prefix):
      removePrefix(prefixToRemove, rootPath)

if __name__=="__main__":
  main( os.path.abspath(os.path.expanduser(sys.argv[1])))
