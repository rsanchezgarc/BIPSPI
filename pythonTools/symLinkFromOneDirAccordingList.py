import sys, os

inDir, restrictionFile, outDir= sys.argv[1:]

prefixes=[]
with open(restrictionFile) as f:
  for line in f:
    lineArray= line.split()
    if len(lineArray)>0:
      prefixes.append( lineArray[0].strip())

for fname in os.listdir(inDir):
  for prefix in prefixes:
    if fname.startswith(prefix):
      baseFname= os.path.basename(fname)
      os.symlink(os.path.join(inDir, baseFname), os.path.join(outDir, baseFname))
      break
