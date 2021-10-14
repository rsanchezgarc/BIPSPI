from __future__ import absolute_import, print_function
import os
from Config import Configuration
from subprocess import check_call


class FeaturesComputer(Configuration):
  '''
  Abstract class. It will be extended with different features computer classes, p.e, PSAIA computer, DSSP computer...
  Intended to be used for computing one type of features each
  '''
  def __init__(self, prefix, computedFeatsRootDir= None, statusManager=None):
    '''
      @prefix. An id for a complex. Example: 1A2K
      :param computedFeatsRootDir: str. path where features will be stored. If None, read from Confinguration
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''
    Configuration.__init__(self)  # Load configuration parameters such as path to programs
    self.prefix= prefix
    self.statusManager= statusManager
    if computedFeatsRootDir!= None:
      self.computedFeatsRootDir= computedFeatsRootDir
    self.computedFeatsRootDir= self.computedFeatsRootDir #Creates root path where features will be saved
    
  def reportStatus(self, msg):
    if not self.statusManager is None:
      self.statusManager.appendToStatus(msg)

  def getExtendedPrefix(self, fname, splitTag="."):
    '''
      Given a filename, obtains its unambiguous id
      :param fname: str. A filename. pe. "/path/to/file/1A2K_l_.updb"
      :return unambiguous id: str. pe "1A2K_l_"  This id will be the prefix of output names
    '''
    return getExtendedPrefix(fname, splitTag= splitTag)
    
  def getPrefix(self, fname):
    '''
      Given a filename, obtains its unambiguous id
      :param fname: str. A filename. e.g. "/path/to/file/1A2K_l_.pdb"
      :return id at complex level: str. e.g. "1A2K"
    '''
    return getPrefix(fname)
  
  def splitExtendedPrefix(self, extendedPrefix):
    '''
      Given an extendedPrefix, splits it into: prefix, chainType ,[chainId]
      :param fname: str. An extendedPrefix. pe. "1A2K_l_"
      :return    prefix,  chainType, ""
             or  prefix,  chainType, chainId
    '''
    return splitExtendedPrefix(extendedPrefix)
    
  def compressFname(self, fname):
    if fname.split(".")[-1]!="gz":
      cmd= ["gzip", "-f", fname]
      check_call(cmd)
  
def getExtendedPrefix( fname, splitTag="."):
  '''
    Given a filename, obtains its unambiguous id
    :param fname: str. A filename. e.g. "/path/to/file/1A2K_l.pdb"
    :return unambiguous id: str. e.g "1A2K_l"  This id will be the prefix of output names
  ''' 
  extendedPrefix= os.path.basename(fname).split(splitTag)[0]
  if extendedPrefix[-1]!="_":
    extendedPrefix+="_"
  return extendedPrefix
  
def getPrefix( fname):
  '''
    Given a filename, obtains its unambiguous id
    :param fname: str. A filename. e.g. "/path/to/file/1A2K_l_.pdb"
    :return id at complex level: str. e.g. "1A2K"
  '''   
  return os.path.basename(fname).split(".")[0].split("_")[0]
  
def splitExtendedPrefix( extendedPrefix):
  '''
    Given an extendedPrefix, splits it into: prefix, chainType ,[chainId]
    :param fname: str. An extendedPrefix. e.g. "1A2K_l_"
    :return    prefix,  chainType, ""
           or  prefix,  chainType, chainId
  '''
  contain= extendedPrefix.split("_")
  assert contain[-1]==""
  contain= contain[:-1]
  return contain

         
class FeatureComputerException(Exception):
  pass
  

