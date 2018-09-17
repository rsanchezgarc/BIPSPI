import sys, os, shutil
from subprocess import call
import requests
import re
from bigExceptions import NoAvailableForDownloadPDB

def myMakeDir(dirPath, pathTail=None):
  '''
    Equivalent to os.path.join but if the resultant path does not exists, then it is created as a directory.
    @param dirPath: str. prefix  of the path.
    @param pathTail: str. tail of the path. Optional
    @return the result of path concatenate dirPath and pathTail (which will be a directory)
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
    
    @param dirPath: str. prefix  of the path.
    @param pathTail: str. tail of the path. Optional
    @return the result of path concatenation of dirPath and pathTail
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
    @param fname: str. file to be removed
  '''    
  try:
    os.remove( fname)
  except OSError:
    pass

def tryToCopy( source, dest):
  '''
    Try to copy source file to dest file. If it is not possible, it does not do anything.
    @param source: str. path to source file
    @param dest: str. path where source will be copied
  '''    
  try:
    shutil.copyfile( source, dest)
  except IOError:
    pass
    
def tryToSymlink( source, dest):
  '''
    Try to copy source file to dest file. If it is not possible, it does not do anything.
    @param source: str. path to source file
    @param dest: str. path where source will be copied
  '''    
  try:
    os.symlink( source, dest)
  except (IOError, OSError):
    pass 

def openForReadingFnameOrGz( fname):
  if fname.endswith(".gz"):
    return gzip.open(fname)  
  elif os.path.isfile(fname+".gz"):
    return gzip.open(fname+".gz")
  else:
    return open(fname)
