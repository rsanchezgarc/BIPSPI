import sys, os
reverseSelection=True
def separate(inPath, outPath, listFname):
  readPrefixes=[]
  with open(listFname) as f:
    for line in f:
      readPrefixes.append( line.split()[0].upper() )
  for fname in os.listdir(inPath):
    isInList= fname.split(".")[0] in readPrefixes
    if (not reverseSelection and isInList) or (reverseSelection and not isInList):
      print(fname)
      fnameOld= os.path.join(inPath, fname)
      newFname= os.path.join(outPath, fname)
      os.symlink(fnameOld, newFname)

if __name__=="__main__":
  if len(sys.argv)!=4:
    raise ValueError("Incorrect number of arguments")
  inPath= os.path.abspath(os.path.expanduser(sys.argv[1]))
  outPath= os.path.abspath(os.path.expanduser(sys.argv[2]))
  listFname= os.path.abspath(os.path.expanduser(sys.argv[3]))
  separate(inPath, outPath, listFname)
