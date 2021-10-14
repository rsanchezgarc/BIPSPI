import sys, os, shutil
from subprocess import Popen, PIPE, check_call
import gzip

from Config import Configuration
from .patchDockParamsTemplate import patchDockTemplate
from pythonTools.rotTransPdb import rotateTranslatePdb
from utils import myMakeDir

BINDING_SITE_FNAME_TEMPLATE= "%(patchDockWdir)s/%(chainType)sBindingSite.tab"
N_MODELS_TO_EXTRACT=10

def writeBindingSite(bindingSiteResidues, patchDockWdir, chainType):
  fname= BINDING_SITE_FNAME_TEMPLATE%locals()

  if isinstance(bindingSiteResidues, str) and os.path.isfile(bindingSiteResidues):
    shutil.copyfile(bindingSiteResidues, fname)
  else:
    with open(fname, "w") as f:
      for res in bindingSiteResidues:
        f.write( " ".join(res)+"\n")
  return fname

def uncompressPdbIfGz(pdbFname):
  if pdbFname.endswith(".gz"):
    with gzip.open(pdbFname, 'rb') as f_in:
      with open(pdbFname[:-3], 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out) 
    return pdbFname[:-3]
  else:
    return pdbFname

def cleanDirectory(lPdbFname, rPdbFname, patchDockWdir, hideFnames=True):

  if lPdbFname.endswith(".pdb"):
    try:
      os.remove(lPdbFname)
    except OSError:
      pass
  if rPdbFname.endswith(".pdb"):
    try:
      os.remove(rPdbFname)
    except OSError:
      pass
  for fname in os.listdir( patchDockWdir):
    fname_full= os.path.join(patchDockWdir, fname)
    try:
      if fname.startswith("results.patchdock"):
        if hideFnames and fname=="results.patchdock":
          cmd='cat %s | grep -Pv "ligandActiveSite|ligandPdb|log-file|protLib|receptorActiveSite|receptorPdb" | gzip > %s'%(fname_full, fname_full+".gz")
          check_call( cmd , shell=True)
        else:
          check_call(['gzip', fname_full])
        if os.path.isfile(fname_full):
          os.remove( fname_full )
      else:
        os.remove( fname_full )
    except IOError:
      continue

def excutePatchDock(lPdbFname, rPdbFname, lBindingSite, rBindingSite, patchDockWdir,
                    writeOnlyLigand=False, cleanWorkingDir=False):
  conf= Configuration()
  patchDockRootDir= conf.patchDockRootDir
  lPdbFname= uncompressPdbIfGz(lPdbFname)
  rPdbFname= uncompressPdbIfGz(rPdbFname)
  
  configStr= patchDockTemplate%{ "lPdbFname":lPdbFname, "rPdbFname":rPdbFname, "patchDockRootDir":patchDockRootDir,
                               "patchDockWdir":patchDockWdir }                   
  configFname= "%(patchDockWdir)s/config.txt"%{"patchDockWdir":patchDockWdir }
  resultsFname= "%(patchDockWdir)s/results.patchdock"%{"patchDockWdir":patchDockWdir }
  myMakeDir(patchDockWdir)
  with open( configFname, "w") as f:
    f.write(configStr)

  writeBindingSite(lBindingSite, patchDockWdir, chainType="l")
  writeBindingSite(rBindingSite, patchDockWdir, chainType="r")
  cmd= [os.path.join(patchDockRootDir,"patch_dock.Linux"),configFname, resultsFname]
  proc = Popen(cmd, stdin= PIPE, stdout=PIPE, stderr=PIPE, cwd=patchDockWdir)
  print( " ".join(cmd) )
  output= proc.communicate()
  print(output[0])
  print("\n?????????????????????????????????????????????????????????????????????\n")
  print(output[1])
  if "error" in output[1]:
    raise Exception("Error executing patchDock")

  listOfSelectedModels=[]
  with open(resultsFname) as f:
    for line in f:
      if "# | score | pen." in line:
        break
    for i, line in enumerate(f):
      if i>=N_MODELS_TO_EXTRACT:
        break
      lineArray= line.split("|")
      if len(lineArray)>0:
        score= float( lineArray[1].strip())
        transformations=  lineArray[-1].split()
        rots= [float(elem) for elem in transformations[:3]]
        trans= [float(elem) for elem in transformations[3:]]
      listOfSelectedModels.append(  (i, score, rots, trans)  )


  for i, score, rots, trans  in listOfSelectedModels:
    rotX, rotY, rotZ= rots
    transX, transY, transZ= trans
    fnameOut= os.path.join(patchDockWdir, "results.patchdock.%d.pdb"%(i+1))
    rotateTranslatePdb(lPdbFname, rotX, rotY, rotZ, transX, transY, transZ,
                       fnameOut=fnameOut)
    if not writeOnlyLigand:
      cmd = "cat %(rPdbFname)s %(fnameOut)s > %(fnameOut)s.tmp && mv %(fnameOut)s.tmp %(fnameOut)s "%locals()
      proc = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True, cwd=patchDockWdir)
      output = proc.communicate()
      print(output[0])
      print(output[1])
      if "error" in output[1]:
        raise Exception("Error concatenating ligand and receptor pdb files")

  if cleanWorkingDir:
    cleanDirectory(lPdbFname, rPdbFname, patchDockWdir)
  print("patchDock DONE")
  return listOfSelectedModels
  
if __name__=="__main__":
  lPdbFname, rPdbFname, lBindingSite, rBindingSite, patchDockWdir= sys.argv[1:]
  excutePatchDock(lPdbFname, rPdbFname, lBindingSite, rBindingSite, patchDockWdir,
                  writeOnlyLigand=True, cleanWorkingDir=False)

  '''
python -m patchDock.launchPatchDock /home/rsanchez/side_projects/use_bipspi/yaoMcGill/docking/tnsCtnsE/inputs/tnscTnse_l_u.pdb  /home/rsanchez/side_projects/use_bipspi/yaoMcGill/docking/tnsCtnsE/inputs/tnscTnse_r_u.pdb /home/rsanchez/side_projects/use_bipspi/yaoMcGill/docking/tnsCtnsE/inputs/ligandBindingSite_1.txt  /home/rsanchez/side_projects/use_bipspi/yaoMcGill/docking/tnsCtnsE/inputs/recptorBindingSite.txt /home/rsanchez/side_projects/use_bipspi/yaoMcGill/docking/tnsCtnsE/patchdock_patch1


pymol for idx in range(1,1001):cmd.load("results.patchdock.%d.pdb"%idx,"mov")

  '''
