from __future__ import absolute_import
import os
import socket

PROYECT_ROOT= os.path.dirname(os.path.realpath(__file__))

hostname= socket.gethostname()


def loadConfig(configFile):
  paramsDict= {}
  with open(configFile) as f:
    for lNum, line in enumerate(f):
      if line.startswith("#"):
        continue
      else:
        lineArray= line.split()
        nElems= len(lineArray)
        if nElems==0:
          continue
        elif nElems!=2:
          raise ValueError("Bad format in configFile %s line %d"%( configFile, lNum))
        else:
          if lineArray[1].strip() in ["True", "False"]:
            paramsDict[lineArray[0]]= "True"==lineArray[1]
          else:
            try:
              paramsDict[lineArray[0]]= int(lineArray[1])
            except ValueError:
              try:
                paramsDict[lineArray[0]]= float(lineArray[1])
              except ValueError:
                paramsDict[lineArray[0]]= os.path.expanduser(lineArray[1])
            
  return paramsDict

if hostname.startswith( "xgbpred"):
  DEFAULT_CONFIG_FILE= os.path.join(PROYECT_ROOT, "configFiles", "configFileWebServer.cfg")
else:
  DEFAULT_CONFIG_FILE= os.path.join(PROYECT_ROOT, "configFiles", "configFile.cfg")

CONFIG_DICT= loadConfig(DEFAULT_CONFIG_FILE)

class Configuration(object):
  '''
   This class contains basic information about paths, parameters and tools that will be used by differente module elements.
   It also contains several utility funcitons such as tryToCopy...
  '''
  def __init__(self, configDict= CONFIG_DICT):
    self.__dict__.update( configDict.copy()) # Update object parameters to include DEFAULT_PARAMETERS elements
    
    assert self.modelType in ["mixed", "struct", "seq", "seq_train"], \
                                  "Error modelType in configFile must be mixed or struct or seq"
    if self.modelType.startswith("seq"): 
      self.modelType= self.modelType.replace("seq", "seq_train")

    assert self.ncpu >=1 ,  "Error ncpu in configFile must be >=1"
    
if __name__== "__main__":
  conf=Configuration()
  print (conf.__dict__)


