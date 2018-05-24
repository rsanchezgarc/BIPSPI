from Config import Configuration

class MyException(Exception):
  def __init__(self, msg):
    Exception.__init__(self, msg)
    self.msg= msg
  
class NoAvailableForDownloadPDB(MyException):
  def __init__(self, msg):
    MyException.__init__(self, msg)
  
class NoValidPDBFile(MyException):
  def __init__(self, msg):
    MyException.__init__(self, msg)
    
class BadNumberOfResidues(Configuration, MyException):
  def __init__(self, nResidues, partnerId):
    Configuration.__init__(self)
    MyException.__init__(self, "Bad number of residues for partner %s: %d. Number of residues must be %d < nResidues < %d"%(
                                 partnerId, nResidues, self.minNumResiduesPartner , self.maxNumResiduesPartner))
    
class BadSequence(Configuration, MyException):
  def __init__(self, msg):
    Configuration.__init__(self)
    MyException.__init__(self, msg)
    
