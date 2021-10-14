from __future__ import absolute_import
import os,sys
import socket
from argparse import ArgumentParser, SUPPRESS

PROYECT_ROOT= os.path.dirname(os.path.realpath(__file__))

def loadConfig(configFilesList):
  paramsDict= {}
  for configFile in configFilesList:
    if not os.path.isfile(configFile):
      raise ValueError("Error, config file (%s) does not exits" % (configFile))
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
                  key = lineArray[0]
                  str_val = os.path.expanduser(lineArray[1])
                  if key.endswith("_path"):
                    key = key.replace("_path", "")
                    paramsDict[key] =  os.path.expanduser(str_val)
                    if  not str_val.startswith("/") and  not str_val.startswith("%") and not "@" in str_val:
                      paramsDict[key] = os.path.abspath(str_val)
                  else:
                    paramsDict[key] = str_val

  return paramsDict

hostname= socket.gethostname()
if hostname.startswith( "xgbpred") or hostname.startswith("bipspi"):  #TODO: Change it to make the same config file
  DEFAULT_CONFIG_DIR = os.path.join(PROYECT_ROOT, "configFiles", "webServer")
  DEFAULT_CONFIG_FILE= os.path.join(DEFAULT_CONFIG_DIR, "configFileWebServer.cfg")
else:
  DEFAULT_CONFIG_DIR = os.path.join(PROYECT_ROOT, "configFiles", "cmdTool")
  DEFAULT_CONFIG_FILE= os.path.join(DEFAULT_CONFIG_DIR, "configFile.cfg")

DEPENDENCIES_CONFIG_FILE= os.path.join(DEFAULT_CONFIG_DIR, "dependencies.cfg")


CONFIG_FILES= [DEFAULT_CONFIG_FILE, DEPENDENCIES_CONFIG_FILE]


CONFIG_DICT = loadConfig(CONFIG_FILES)

class Configuration(object):
  '''
   This class contains basic information about paths, parameters and tools that will be used by differente module elements.
   It also contains several utility funcitons such as tryToCopy...
  '''

  CONFIG_FLAG = "configFile"
  DEPENDENCIES_FLAG = "configDependencies"
  def __init__(self):

    self.already_set_attributes = set([])

    self.check_args()

  def __getattribute__(self,  name):
    if name == "already_set_attributes":
      return super(Configuration, self).__getattribute__("already_set_attributes")

    if name in CONFIG_DICT and name not in super(Configuration, self).__getattribute__("already_set_attributes") :
        val =   CONFIG_DICT.get(name, None)
        if isinstance(val, str) and "%" in val:
          val = val % CONFIG_DICT
        return val
    else:
        return super(Configuration, self).__getattribute__(name)


  def __setattr__(self, key, value):
    if key != "already_set_attributes":
      self.already_set_attributes.add(key)
    super(Configuration, self).__setattr__( key, value)

  def getAllAttrDict(self):
    dict1 = self.__dict__.copy()
    dict1.update( CONFIG_DICT)
    return dict1

  @classmethod
  def getArgParser(cls, *args, **kwargs):
    parser = MyArgParser(*args, **kwargs)

    parser.add_argument("-c", "--" + cls.CONFIG_FLAG, type=str, help="Filename for the per-project config file", default=DEFAULT_CONFIG_FILE)
    parser.add_argument("-d", "--" + cls.DEPENDENCIES_FLAG, type=str, help="Filename for the dependencies config file", default=DEPENDENCIES_CONFIG_FILE)

    for key, val in CONFIG_DICT.iteritems():
      parser.add_argument("--" + key, type=type(val), default=val,  help=SUPPRESS)

    parser.modify_field("psiBlastNThrs",  _type=int, help="Number of threads for psi-BLAST", default=CONFIG_DICT.get("psiBlastNThrs"))
    return parser

  @classmethod
  def getRawConfig(self):
    return CONFIG_DICT

  def check_args(self):
    assert self.ncpu >= 1, "Error ncpu in configFile must be >=1"
    # assert os.path.isdir(self.savedModelsPath), "Error, savedModelsPath (%s) does not exists or is not a directory" % self.savedModelsPath
    assert self.minNumResiduesPartner>0, "Error, minNumResiduesPartner must be >0"
    if not hasattr(self, "PREDICT_ONLY"):
      assert self.modelType in ["mixed", "struct", "seq"], \
                                  "Error modelType in configFile must be mixed or struct or seq"

  @classmethod
  def int_or_filePath(cls, value):
    try:
      return int(value)
    except ValueError:
      return cls.file_path(value)

  @classmethod
  def file_path(cls, str_val):
    val = str_val
    if not str_val.startswith("/") and not str_val.startswith("%") and not "@" in str_val:
     val = os.path.abspath(str_val)
    return val

  @classmethod
  def update_config_dict(cls, project_config=None, depend_confi=None, other_dict={}):

    if project_config is not None:
      os.path.abspath(os.path.expanduser(project_config))
      assert os.path.isfile(project_config), "Error, %s does not exists" % project_config
      CONFIG_FILES[0] = project_config
    if depend_confi is not None:
      os.path.abspath(os.path.expanduser(depend_confi))
      assert os.path.isfile(depend_confi), "Error, %s does not exists" % depend_confi
      CONFIG_FILES[1] = depend_confi

    if project_config is not None or depend_confi is not None:
      global CONFIG_DICT
      CONFIG_DICT = loadConfig(CONFIG_FILES)

    if other_dict is not None:
      CONFIG_DICT.update( other_dict)
      Configuration().check_args()

class MyArgParser(ArgumentParser):

    def parse_args(self, args=None, namespace=None):
      out = vars(super(MyArgParser, self).parse_args(args, namespace))
      project_config =  out.get(Configuration.CONFIG_FLAG, None)
      depend_confi = out.get(Configuration.DEPENDENCIES_FLAG, None)
      # if  out.get(Configuration.CONFIG_FLAG, None):
      #   fname = os.path.abspath(os.path.expanduser(out[Configuration.CONFIG_FLAG]))
      #   assert  os.path.isfile(fname), "Error, %s does not exists"%fname
      #   CONFIG_FILES[0] = fname
      # if out.get(Configuration.DEPENDENCIES_FLAG, None):
      #   fname = os.path.abspath(os.path.expanduser(out[Configuration.DEPENDENCIES_FLAG]))
      #   CONFIG_FILES[1] =fname
      #   assert  os.path.isfile(fname), "Error, %s does not exists"%fname

      Configuration.update_config_dict(project_config, depend_confi, out )
      # global CONFIG_DICT
      # CONFIG_DICT = loadConfig(CONFIG_FILES)
      # CONFIG_DICT.update( out)
      # Configuration().check_args()

      return out

    def modify_field(self, name, help=None, _type=None, default=None, choices=None):
      for action in self._actions:
        if action.dest == name:
          if help:
            action.help= "Overwrites config file content. " +help
          if _type:
            action.type= _type
          if default is not None:
            action.default= default
          if choices is not None:
            action.choices = choices


if __name__== "__main__":
  conf=Configuration()
  print (conf.getAllAttrDict())

  parser = Configuration.getArgParser(prog="BIPSPI v2")
  args = parser.parse_args()
  print("--------------------------------------")
  print (conf.getAllAttrDict())
  print(conf.wdir)
  conf.wdir = "kk"
  print(conf.wdir)

  # print(CONFIG_DICT["wdir"])
  # print(CONFIG_DICT["computedFeatsRootDir"])
  # print(conf.computedFeatsRootDir)
