import sys, os

def makeFakeBoundFromUnboud(inPath):

  for fname in os.listdir(inPath):
    if fname[:4].isupper():
      os.remove(os.path.join(inPath, fname)); continue

    old_name= os.path.join(inPath, fname)
    new_name= os.path.join(inPath, fname.replace("_u.", "_b.") )
    os.symlink(old_name, new_name)


def copyCapitalExamples(inPath, outPath):
  for fname in os.listdir(inPath):
    if fname[:4].isupper():
      old_name= os.path.join(inPath, fname )
      new_name= os.path.join(outPath, fname )
      os.symlink(old_name, new_name)


if __name__=="__main__":
  inPath= os.path.abspath(os.path.expanduser(sys.argv[1]))
  makeFakeBoundFromUnboud(inPath)

  # outPath= os.path.abspath(os.path.expanduser(sys.argv[2]))
  # copyCapitalExamples(inPath, outPath)