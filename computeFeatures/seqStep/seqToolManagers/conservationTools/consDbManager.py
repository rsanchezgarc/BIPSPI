from __future__ import absolute_import
import sys, os
from subprocess import Popen, PIPE, check_output
import joblib
from Config import Configuration
import sqlite3

consDbSqlite= "/home/rsanchez/Tesis/pssmDB/data/unirefData100/smallPdbOut/pssmDB_v11_u90and100.db"
consDbFilesPath="rsanchez@campins:~/Tesis/pssmDB/data/databaseUpdate2/uniref90"

class ConsDbManager(Configuration):
  def __init__(self, consDbSqlite=None, consDbFilesPath=None):
    Configuration.__init__(self)
    self.isReady= True
    if consDbFilesPath:
      self.consDbFilesPath= consDbFilesPath
    if consDbSqlite:
      self.consDbSqlite= consDbSqlite
      
    self.unirefType=None
    if not os.path.isfile(self.consDbSqlite) or not self.checkIfDbFilesAvailable():
      self.isReady=False
    else:      
      self.unirefType= os.path.basename(self.consDbFilesPath )
      assert self.unirefType in ["uniref90", "uniref100"], ("Error, consDbFilesPath %s "+\
                                              ":must be path/to/data/[uniref90|uniref100]")%(self.consDbFilesPath)
                       
    self.sqliteConn= sqlite3.connect(self.consDbSqlite)
    self.sqliteCursor= self.sqliteConn.cursor()
    try: #check if sqlite was correctly opened
      self.sqliteCursor.execute("SELECT seqId FROM sequencesTable where sequence== 0").fetchone()
    except sqlite3.OperationalError:
      self.isReady=False

  def close(self):
    try:
      self.sqliteConn.close()
    except Exception:
      pass
    
  def __del__(self):
    self.close()
  def __exit__(self, exc_type, exc_value, traceback):
    self.close()

  def consDbIsAvailable(self): 
    return self.isReady

  def checkIfDbFilesAvailable(self):
    
    try:
      host, filesPath= self.consDbFilesPath.split(":")
      assert not ";" in filesPath, "Error, consDbFilesPath contain ';' "
      cmd= "ssh -o BatchMode=yes -o ConnectTimeout=10 %s 'ls %s '"%( host, filesPath)
      print(cmd)      
      process= Popen( cmd, shell=True, stdout=PIPE, stderr=PIPE)
      processOut= process.communicate()
      if len(processOut[1])>0:
        return False
      else:
        return True
        
    except ValueError:
      return os.path.isdir( self.consDbFilesPath )

  def transferZip(self, relatPathInDdFiles, outName):
    try:
      host, filesPath= self.consDbFilesPath.split(":")
      cmd= "ssh -o BatchMode=yes -o ConnectTimeout=10 %s 'zcat %s ' > %s "%( host, os.path.join(filesPath,relatPathInDdFiles), outName)        
    except ValueError:
      cmd= "zcat %s > %s "%(  os.path.join(self.consDbFilesPath, relatPathInDdFiles), outName)    
    print(cmd)
    process= Popen( cmd, shell=True, stdout=PIPE, stderr=PIPE)
    processOut= process.communicate()
    if len(processOut[1])>0:
      return False
    else:
      if outName.endswith(".gz"):
        import gzip
        openFun= gzip.open
      else:
        openFun= open
      if not os.path.isfile(outName):
        return False
      okWords= ["PSI Gapped", "Window for multiple hits:"]
      with openFun(outName) as f:
        for line in f:
          for okword in okWords:
            if line.startswith(okword):
              return True
      os.remove(outName)
      return False
    
  def retrieve3DConsFromSeq(self, seqStr, pssmName, psioutName=None):
    '''
      It returns values according statusTable codes in database. 0 if ok.
    '''
    print("Trying to recover 3dcons pssm from sequence search")
    seqId= self.sqliteCursor.execute("SELECT seqId FROM sequencesTable where sequence== ?" ,(seqStr,)).fetchone()
    if seqId is None:
      return -1
    else:
      seqId= seqId[0]
      if self.unirefType=="uniref90":
        status= self.sqliteCursor.execute("SELECT status_uniref90 FROM pssmsTableUniref where seqId== ?" ,
                                                                                          (seqId,)).fetchone()
      elif self.unirefType=="uniref100":
        status= self.sqliteCursor.execute("SELECT status_uniref100 FROM pssmsTableUniref where seqId== ?" ,
                                                                                          (seqId,)).fetchone()
      else:
        raise ValueError("unirefType must be uniref90 or uniref100: now it is %s"%(self.unirefType))
      if status is None:
        return -1
      elif status[0]!= 0 : 
        return status[0]
        
      if psioutName:
        wasOk= self.transferZip(os.path.join("iterNum3","%d.step3.psi_out.zip"%(seqId)), psioutName)
        if not wasOk:
          return -1
      wasOk= self.transferZip(os.path.join("iterNum3","%d.step3.pssm.zip"%(seqId)), pssmName)
      if not wasOk:
        return -1
      else:
        return 0

  def retrieve3DConsFromPDBChain(self, pdbCode, chainId, pssmName, psioutName=None):
    '''
      It returns values according statusTable codes in database. 0 if ok.
    '''
    print("Trying to recover pssm for pdb %s:%s"%(pdbCode, chainId))
    pdbCode= pdbCode.lower()
    seqId= self.sqliteCursor.execute("SELECT seqId FROM pdbsTable where pdb== ? and chain== ?" ,(
                                      pdbCode,chainId)).fetchone()
    if seqId is None:
      return -1
    else:
      seqId= seqId[0]                                         
      if self.unirefType=="uniref90":
        status= self.sqliteCursor.execute("SELECT status_uniref90 FROM pdbChainIterStatus where pdb== ? and chain== ?" ,(
                                          pdbCode,chainId)).fetchone()
      elif self.unirefType=="uniref100":
        status= self.sqliteCursor.execute("SELECT status_uniref100 FROM pdbChainIterStatus where pdb== ? and chain== ?" ,(
                                          pdbCode,chainId)).fetchone()
      else:
        raise ValueError("unirefType must be uniref90 or uniref100: now it is %s"%(self.unirefType))
                     
      if status is None:
        return -1
      elif status[0]!= 0 : 
        return status[0]
        
      if psioutName:
        wasOk= self.transferZip(os.path.join("iterNum3","%d.step3.psi_out.zip"%(seqId)), psioutName)
        if not wasOk:
          return -1
      wasOk= self.transferZip(os.path.join("pssmInZipWithStructId","%s_%s_3.pssm.zip"%(pdbCode, chainId)), pssmName)
      if not wasOk:
        return 2
      else:
        return 0

def testSeq():

  consManager= ConsDbManager()
  seq="TEFGSELKSFPEVVGKTVDQAREYFTLHYPQYDVYFLPEGSPVTLDLRYNRVRVFYNPGTNVVNHVPHVG"
  pssmName= os.path.expanduser("~/tmp/tmp3dcons/pssmFromSeq.txt")
  psiName= os.path.expanduser("~/tmp/tmp3dcons/aligsFromSeq.txt")
  print(seq)
  returnCode=consManager.retrieve3DConsFromSeq(seq, pssmName, psiName)
  print("return code %d"%(returnCode))
  seq="CACA"
  print(seq)
  returnCode=consManager.retrieve3DConsFromSeq(seq, pssmName, psiName)
  print("return code %d"%(returnCode))

def testPdbId():
  consManager= ConsDbManager()
  print("1acb","E")
  returnCode=consManager.retrieve3DConsFromPDBChain("1acb","E", "/home/rsanchez/tmp/tmp3dcons/pssm1.txt", "/home/rsanchez/tmp/tmp3dcons/aligs1.txt")
  print("return code %d"%(returnCode))
  print("1acb","W")
  returnCode=consManager.retrieve3DConsFromPDBChain("1acb","W", "/home/rsanchez/tmp/tmp3dcons/pssm2.txt", "/home/rsanchez/tmp/tmp3dcons/aligs2.txt")
  print("return code %d"%(returnCode))
  print("4zez","E")
  returnCode=consManager.retrieve3DConsFromPDBChain("4zez","E", "/home/rsanchez/tmp/tmp3dcons/pssm3.txt", "/home/rsanchez/tmp/tmp3dcons/aligs3.txt")
  print("return code %d"%(returnCode))
  
  print("1A2K","A")
  returnCode=consManager.retrieve3DConsFromPDBChain("1A2K","A", "/home/rsanchez/tmp/tmp3dcons/pssm4.txt", "/home/rsanchez/tmp/tmp3dcons/aligs3.txt")
  print("return code %d"%(returnCode))

if __name__=="__main__":
  testSeq()
#  testPdbId()
  

