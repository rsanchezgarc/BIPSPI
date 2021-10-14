import os
from unittest import TestCase

import shutil

RESULTS_PATH="/home/rsanchez/tmp/BIPSPI/output/trial_ori_1"
class TestPredictOneComplex(TestCase):

  def _testPredictOneComplex(self, lPdbId=None, lFname=None, lSequence=None,
                                   rPdbId=None, rFname=None, rSequence=None, isHomo=False):
    shutil.rmtree(RESULTS_PATH)
    os.mkdir(RESULTS_PATH)
    print(lPdbId, lFname, lSequence, rPdbId, rFname, rSequence, isHomo)
    from predictComplexes import predictOneComplex
    outName, (p, l, r) = predictOneComplex(allRootDir=RESULTS_PATH,
                                           lPdbId=lPdbId, lFname=lFname, lSequence=lSequence,
                                           rPdbId=rPdbId, rFname=rFname, rSequence=rSequence, isHomoComplex=isHomo,
                                           areForTrainAndTest=False, removeInputs=False)



  def testPredictOneComplexPdb_None(self):
    lFname = "/home/rsanchez/Tesis/rriPredMethod/data/develData/pdbFiles/1ACB_l_u.pdb"
    self._testPredictOneComplex(lFname=lFname, isHomo=True)

  def testPredictOneComplexPdb_Pdb(self):
    lFname = "/home/ruben/Tesis/rriPredMethod/docs_v2/useCases/local_4ov6/input/ED/4ov6:E_ED.pdb"
    rFname = "/home/ruben/Tesis/rriPredMethod/docs_v2/useCases/local_4ov6/input/ED/4ov6:D_ED.pdb"

    self._testPredictOneComplex(lFname=lFname, rFname=rFname, isHomo=False)

  def testPredictUseCase(self):
    lFname = "/home/rsanchez/Tesis/rriPredMethod/data/develData/pdbFiles/1ACB_l_u.pdb"
    rFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/pdbFiles/1ACB_r_u.pdb"

    self._testPredictOneComplex(lFname=lFname, rFname=rFname, isHomo=False)


  def _other(self):

    lFname = "/home/rsanchez/Tesis/rriPredMethod/data/develData/pdbFiles/1ACB_l_u.pdb"
    # rFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/pdbFiles/1ACB_r_u.pdb"
    # rFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/inputFasta/seq1_r_u.fasta"
    # lFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/inputFasta/seq1_l_u.fasta"
    # rFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/inputFasta/1ACB_r_u.fasta"
    # lFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/inputFasta/1ACB_l_u.fasta"
    # lFname="/home/rsanchez/tmp/tmpRRI/inPdb/6f1y-f:j_l_u.pdb"
    # rFname="/home/rsanchez/tmp/tmpRRI/inPdb/6f1y-f:j_r_u.pdb"
    # lFname="/home/rsanchez/tmp/BIPSPI/input/inputSeq_l_u.fasta"
    # rFname="/home/rsanchez/tmp/BIPSPI/input/inputSeq_l_u.fasta"

    #  lPdbId="1acb"
    #  rPdbId="1acb"
    #  lPdbId="4OV6:E"
    #  rPdbId="4OV6:D"
    #  lSequence=">1ACB:E|PDBID|CHAIN|SEQUENCE\nCGVPAIQPVLSGLSRIVNGEEAVPGSWPWQVSLQDKTGFHFCGGSLINENWVVTAAHCGVTTSDVVVAGEFDQGSSSEKIQKLKIAKVFKNSKYNSLTINNDITLLKLSTAASFSQTVSAVCLPSASDDFAAGTTCVTTGWGLTRYTNANTPDRLQQASLPLLSNTNCKKYWGTKIKDAMICAGASGVSSCMGDSGGPLVCKKNGAWTLVGIVSWGSSTCSTSTPGVYARVTALVNWVQQTLAAN"
    #  rSequence=">1ACB:I|PDBID|CHAIN|SEQUENCE\nTEFGSELKSFPEVVGKTVDQAREYFTLHYPQYDVYFLPEGSPVTLDLRYNRVRVFYNPGTNVVNHVPHVG"
    #  lPdbId="5by8:B"
    #  rPdbId="5by8:A"

    # isHomo=False
    isHomo = True

    outName, (p, l, r) = predictOneComplex(allRootDir="/home/rsanchez/tmp/BIPSPI/output/trial_ori_1",
                                           lPdbId=lPdbId, lFname=lFname, lSequence=lSequence,
                                           rPdbId=rPdbId, rFname=rFname, rSequence=rSequence, isHomoComplex=isHomo,
                                           areForTrainAndTest=False, removeInputs=False)

