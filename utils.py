import os, shutil
from operator import itemgetter
import gzip
from Config import Configuration
conf=Configuration()

def myMakeDir(dirPath, pathTail=None):
  '''
    Equivalent to os.path.join but if the resultant path does not exists, then it is created as a directory.
    :param dirPath: str. prefix  of the path.
    :param pathTail: str. tail of the path. Optional
    :return the result of path concatenate dirPath and pathTail (which will be a directory)
  '''  
  if pathTail != None:
    dirPath= os.path.join(dirPath, pathTail)
  if not os.path.isdir(dirPath):
    try:
      os.makedirs(dirPath)
    except OSError as e:
      if e.errno== 17: # Directory already exists
        pass
  return dirPath

def myMakeDirUnique(dirPath, pathTail=None):
  '''
    Equivalent to os.path.join but if the resultant path exists, it adds a suffix
    of form _#Nun
    P.e. 
    myMakeDirUnique("/hola/path", "foo")
    will create dir "/hola/path/foo"
    If executed again,     myMakeDirUnique("/hola/path", "foo")
    will create dir "/hola/path/foo_1"    
    If executed again,     myMakeDirUnique("/hola/path", "foo")
    will create dir "/hola/path/foo_2"
    and so on
    and will return the new name
    
    :param dirPath: str. prefix  of the path.
    :param pathTail: str. tail of the path. Optional
    :return the result of path concatenation of dirPath and pathTail
  ''' 

  if pathTail != None:
    if pathTail.endswith("/"):
      pathTail= pathTail[:-1]
    dirPath= os.path.join(dirPath, pathTail)
    
  origDirPath=  dirPath    
  i=1
    
  if os.path.isdir(origDirPath):
    dirPath= origDirPath+"_"+str(i)
    i+=1
  else:
    os.makedirs(dirPath)
    return dirPath
  while os.path.isdir(dirPath):
    dirPath= origDirPath+"_"+str(i)
    i+=1
  os.makedirs(dirPath)
  return dirPath
  
def tryToRemove( fname):
  '''
    Try to remove one file. If it is not possible, it does not do anything.
    :param fname: str. file to be removed
  '''    
  try:
    os.remove( fname)
  except OSError:
    pass

def tryToMove(source, dest):
  '''
    Try to move source file to dest file. If it is not possible, it does not do anything.
    :param source: str. path to source file
    :param dest: str. path where source will be copied
  '''    
  try:
    os.rename( source, dest)
  except IOError:
    pass

def tryToCopy( source, dest):
  '''
    Try to copy source file to dest file. If it is not possible, it does not do anything.
    :param source: str. path to source file
    :param dest: str. path where source will be copied
  '''    
  try:
    shutil.copyfile( source, dest)
  except IOError:
    pass
    
def tryToSymlink( source, dest):
  '''
    Try to copy source file to dest file. If it is not possible, it does not do anything.
    :param source: str. path to source file
    :param dest: str. path where source will be copied
  '''    
  try:
    os.symlink( source, dest)
  except (IOError, OSError):
    pass 
    

def tryToCleanDir(dirName, substr="_", rootDataDir=conf.computedFeatsRootDir):
  for name in os.listdir(dirName):
    if not substr or substr in name:
      nameToRemove= os.path.join(dirName, name)
      assert nameToRemove.startswith(rootDataDir ), "Error, trying to remove not allowed file %s"%(nameToRemove)
      os.remove(os.path.join(dirName, name))
      
def openForReadingFnameOrGz( fname):
  if fname.endswith(".gz"):
    return gzip.open(fname)  
  elif os.path.isfile(fname+".gz"):
    return gzip.open(fname+".gz")
  else:
    return open(fname)



def getItemsFromList(idxs, l):
  l=  itemgetter( *idxs )( l)
  if len(idxs)==1:
    l= [l]
  return l
  
def checkFreeMemory():
  '''
  return free memory in GBytes
  '''
  import psutil
  x=psutil.virtual_memory()
  return x.free/(1024.0 ** 3)

def getTotalMemory():
  '''
  return free memory in GBytes
  '''
  import psutil
  x=psutil.virtual_memory()
  return x.total/(1024.0 ** 3)

def getFileSize(fname):
  '''
  return file disk usage in GBytes
  '''
  statinfo = os.stat(fname)
  return statinfo.st_size/(1024.0 ** 3)