from __future__ import print_function
import getopt
import sys, os
from multiprocessing import cpu_count
from .trainAndTestWorker import TrainAndTestWorker

class TrainAndTest(object):
  '''
    This class performs training and evaluation of a data set by cross validation
  '''
  def __init__(self):
    '''
      Builder. No parameters. See self.setParameters() to know about needed parameters.
      self.setParameters() should be called before main if this program is not used as
      the main program.
    '''
    self.sampledComplexesPath= None
    self.wholeComplexesPath= None
    self.predictOutputPath= None
    self.numProc= 2
    self.nFolds= -1 #Leave-one-out if nFolds==-1. Testing independent testSet if nFolds==0
    self.verbose= False
    self.saveModelFname= None
    self.doGridSearch=False
    self.areParamSet= False
    
  def setParameters(self, sampledComplexesPath, wholeComplexesPath, predictOutputPath, numProc,
                    nFolds=-1, verbose=True, saveModelFname=None, doGridSearch=False):
    '''
       Sets object attibutes needed to perform training and testing
       
       @param sampledComplexesPath: str. Path to a dir where training data files are stored
       @param wholeComplexesPath: str. Path to a dir where testing data files are stored. None if doGridSearch==True
       @param predictOutputPath: str. Path to a dir where predictions will be stored
       @param predictOutputPath: str. Path to a dir where predictions will be stored. None if they won't be saved
       @param numProc: int. Number of processes to use in parallel. If -1 all cpus will be used
       @param nFolds: int. Number of folds for k-fold cross-validation. If -1, leave-one-out will be performed.
                           If 0, testing set will be used as if it were independent. Cannot be 1
       @param verbose: boolean. Whether or not print to stdout info
       @param saveModelFname: str. A path where the final model, trained with all data will be saved. If None,
                                  model won't be saved
       @param doGridSearch: boolean. Whether or not to perform a grid search for hyperparameter tunning
    '''                    
    self.sampledComplexesPath= sampledComplexesPath
    self.wholeComplexesPath= wholeComplexesPath
    self.predictOutputPath= predictOutputPath
    self.numProc= numProc
    self.nFolds= nFolds
    self.verbose= verbose
    self.saveModelFname= saveModelFname
    self.doGridSearch= doGridSearch
    self.areParamSet= True
    
  def main(self, cmds= sys.argv ):
    '''
       Performs train and test operations calling self.doTrainAndTest()
       
       @param cmds: list whith arguments or None. By default it uses console args if
                    attributes were not set by self.setParameters(). List of arguments 
                    must follow de instructions showed in self.usage() 
    '''
    if not self.areParamSet and not cmds is None:
      self.parseCmdArgs(cmds)
    self.checkIfParametersAreValid()
    return self.doTrainAndTest()
      
  def doTrainAndTest(self):
    '''
       Performs train and test by doing either cross validation or train-test split. Returns a pandas.DataFrame which is
       a summary of the performance evaluation
       
       @return crossValidationSummaryResults: pandas.DataFrame. A summary of the evaluation results or 
                                              None if self.doGridSearch==True
    '''  
    print ("Training:\nInput %s"%(self.sampledComplexesPath))
    print ("Test %s"%(self.wholeComplexesPath))
    if self.predictOutputPath!=None:
      print ("Output %s"%(self.predictOutputPath))
    if not self.saveModelFname is None:
      print ("Model save fname: %s"%self.saveModelFname)
    print ( "running in %d cpus"%(self.numProc))

    trainerAndTester= TrainAndTestWorker(self.sampledComplexesPath, self.wholeComplexesPath, self.predictOutputPath, 
                                           self.nFolds, self.saveModelFname, self.verbose, self.numProc)
    if self.doGridSearch:
      trainerAndTester.doGridSearch()
    else:
      return trainerAndTester.computeTrainAndTest()

  def parseCmdArgs(self,cmds):
    '''
       Parse cmds to set object attributes. cmds should be like sys.argv
       
       @param cmds: list whith arguments. For full description see self.usage()
    '''  
    options, remainder= getopt.getopt(cmds[1:], 'i:t:o:n:j:k:vs:g',
                                     ['input=',
                                      'test=',
                                      'output=',
                                      'numProc=',
                                      'nFolds=',
                                      'verbose',
                                      'saveModelFname=',
                                      'gridSearch'
                                     ])

#    print("options:",cmds[1:], options, remainder)
#    raw_input()
    for opt, arg in options:
      if opt in   ('-i', '--input'):
        self.sampledComplexesPath = os.path.abspath(os.path.expanduser(arg))
      elif opt in ('-t', '--test'):
        self.wholeComplexesPath = os.path.abspath(os.path.expanduser(arg))
      elif opt in ('-o', '--output'):
        self.predictOutputPath = os.path.abspath(os.path.expanduser(arg))
      elif opt in ('-j', '--numProc'):
        self.numProc = int(arg)
      elif opt in ('-k', '--numProc'):
        self.nFolds = int(arg)        
      elif opt in ('-v', '--verbose'):
        self.verbose = True
      elif opt in ('-s', '--saveModelFname'):
        self.saveModelFname = True
      elif opt in ('-g', '--gridSearch'):
        self.doGridSearch = True
        
  def usage(self):
    '''
      returns an str with the usage of this program
       
       @return usage: str. Usage of this program.
    '''  
      
    descriptionStr='''
python -m trainAndTest.trainAndTest  [options]
-i, --input.         input dirPath to train data. Required. It must contain pickled data (joblib.dump).

-t --test.           test dirPath to test data. Required if -g is not selected. It must contain pickled data (joblib.dump).

-o --output.         output dirPath where results want to be saved. If not specified, results will not be saved and just
                      performance evaluation will be printed out to stdout. If output dir already exists the program
                      will not override previous results and continue evaluation

-j --numProc.        number of processors to use (Parallelization). If -1, all cpu's will be used. Default 2

-k --nFolds          number of folds to performe k-fold cross-validation. If no specified, leave-one-out  will be performed.
                     If 0, testing set will be used as if it where independent. Cannot be 1

-v --verbose.        print out information about training and evaluation process.

-s --saveModelFname  path where model will be saved. If no specified, no model will be saved

-g --gridSearch      To perform gridSearch for hyperparameter tunning

'''
    return descriptionStr
    
  def checkIfParametersAreValid(self):
    '''
       Checks if parameters are correct.
    '''
    if self.sampledComplexesPath== None or not os.path.isdir(self.sampledComplexesPath):
      print (self.usage())
      print ("Error, there is no train input path or it is not valid: %s"%(self.sampledComplexesPath))
      sys.exit(1)

    if (self.doGridSearch==False and (self.wholeComplexesPath== None or not os.path.isdir(self.wholeComplexesPath))):
      print (self.usage())
      print ("Error, there is no test path or it is not valid: %s"%(self.wholeComplexesPath))
      sys.exit(1)    

    if self.predictOutputPath!= None and not os.path.isdir(self.predictOutputPath):
      print (self.usage())
      print ("Error output path does not exists it is not valid:%s"%(self.predictOutputPath))
      sys.exit(1)
    elif self.predictOutputPath!= None and len(os.listdir(self.predictOutputPath))>0:
      print ("There are previous results in outPath, evaluation will continue")
      
    if not self.saveModelFname is None:
      if os.path.isdir(self.saveModelFname):
        self.saveModelFname= os.path.join(self.saveModelFname,"rriPred.model.pkl")
      if os.path.isfile(self.saveModelFname):
        print (self.usage())
        print ("Error saveModelFname already exists: %s"%(self.saveModelFname))
        sys.exit(1)
        
    if self.numProc==-1:
      self.numProc= cpu_count()
    elif self.numProc<0:
      self.numProc= 1

    if self.nFolds<-1:
      print (self.usage())
      print ("Error nFolds must be >=-1")
      sys.exit(1)
    if self.nFolds==1:
      print (self.usage())
      print ("Error nFolds cannot be ==1")
      sys.exit(1)
    return True
    
def testModule():
  trainAndTest= TrainAndTest()
  sampledComplexesPath=os.path.expanduser("~/Tesis/rriPredMethod/data/ppdockingBenchData/newCodeData/codifiedInput/struct_2/sampledInputs")
  wholeComplexesPath=os.path.expanduser("~/Tesis/rriPredMethod/data/ppdockingBenchData/newCodeData/codifiedInput/struct_2/allInputs/")
#  predictOutputPath=os.path.expanduser("~/Tesis/rriPredMethod/data/ppdockingBenchData/tmpResults")
  predictOutputPath=None
  numProc=4
  nFolds=-1
  verbose=True
  saveModelFname=None
  doGridSearch= True
  trainAndTest.setParameters(sampledComplexesPath, wholeComplexesPath, predictOutputPath, numProc,
                             nFolds, verbose, saveModelFname, doGridSearch=True)
  trainAndTest.main()
if __name__== "__main__":
#  testModule()
  TrainAndTest().main()
  
