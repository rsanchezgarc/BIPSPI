import numpy as np
import itertools

def getBasePrefix(prefix):
  return prefix.split(":")[0].split("-")[0].split("@")[0].split("#")[0].upper()
  
def getScopeGroups(trainPrefixes, scopeFname, ligandReceptorOrder=True, considerPairsInsteadMonomer=True):
  '''

  scope independency will be carried out as pairs of families, e.g (A-B) (A-B') are independent (A-B) (A'-B') not
  
  given a path to a file that contains scope info as bellow, splits train and test prefixes
  as leave-one-scope-family-out

  <<<< scopeFname example >>>
  1qfw	IM	AB	b.1.1.1:l.1.1.1;b.1.1.1	g.17.1.4;g.17.1.4	2
  2jel	HL	P	b.1.1.1:b.1.1.2;b.1.1.1:b.1.1.2	d.94.1.1	2
  1avx	A	B	b.47.1.2	b.42.4.1	2
  1ay7	A	B	d.1.1.2	c.9.1.1	2
  1buh	A	B	d.144.1.7	d.97.1.1	2
  1bvn	P	T	b.71.1.1:c.1.8.1	b.5.1.1	2
  1clv	A	I	b.71.1.1:c.1.8.1	g.3.2.1	4
  1d6r	A	I	b.47.1.2	g.3.13.1	2
  1dfj	E	I	d.5.1.1	c.10.1.1	2
  1e6e	A	B	c.3.1.1:c.4.1.1	d.15.4.1	2

  or json file like
  {"1AK4": {"r": [["A", "c-134"]], "l": [["D", "c-1694"]]}, "3K58": {"r": [["A", "c-545"]], "l": [["B", "c-719"]]}, "4H03": {"r": [["A", "c-33"]], "l": [["A", "c-136"]]}
  
  :param ligandReceptorOrder: True if first scope column is for ligand and second is for receptor.
                              if False first scope column is for receptor and second is for ligand.
                              
  :param considerPairsInsteadMonomer: Ensure independency at the scopes pairs levels instead at the monomers level
  '''

  print("ensuring scope independency")
  if scopeFname.endswith(".json"):
    prefixToScope= loadFromJson(scopeFname)
  else:
    prefixToScope= loadFromTable(scopeFname, ligandReceptorOrder)

  scopeGroups={ prefix: (set([prefix]),i) for i, prefix in enumerate(trainPrefixes)}
  for i, prefix_i in enumerate(trainPrefixes):
    test_ix=[]
    train_ix=[]
    prefix_i_base=  getBasePrefix(prefix_i)
    if prefix_i_base in prefixToScope:
      scope_L, scope_R= prefixToScope[prefix_i_base]
    for j, prefix_other in enumerate(trainPrefixes):
      prefix_other_base= getBasePrefix(prefix_other)
      if prefix_other_base in prefixToScope:
        scope_L_other, scope_R_other= prefixToScope[prefix_other_base]
        if considerPairsInsteadMonomer:
          condition= (scope_L.intersection(scope_L_other) and scope_R.intersection(scope_R_other) or
            scope_L.intersection(scope_R_other) and scope_R.intersection(scope_L_other) )
        else:
          scopes1= scope_L.union(scope_R)
          scopes2= scope_L_other.union(scope_R_other)
          condition= scopes1.intersection(scopes2)
        if condition:
          #Then merge scope groups
          group_i, idx_i= scopeGroups[prefix_i]
          group_other, idx_others= scopeGroups[prefix_other]
          merged_group= group_i.union(group_other)
          merged_idx= min(idx_i, idx_others)
          scopeGroups[prefix_i]= (merged_group, merged_idx)
          scopeGroups[prefix_other]= (merged_group, merged_idx)
    

  prefix_to_idx= { prefix:idx for idx, prefix in enumerate(trainPrefixes)}
  # print(prefix_to_idx)
  groups=[ 0 for elem in trainPrefixes]
  last_g_id=0
  for i, (prefixes, g_id) in enumerate(sorted([ (prefixes, p_id) for prefixes, p_id in scopeGroups.values() ], 
                                            key=lambda x: x[1])):
    if g_id!= last_g_id:
      last_g_id= g_id
      for prefix in prefixes:
        groups[ prefix_to_idx[prefix] ] = g_id

  return groups
  
def loadFromTable(fname, ligandReceptorOrder):
  prefixToScope={}
  with open(fname) as f:
    for line in f:
      if line.startswith("#"): continue
      lineArray= line.split()
      if ligandReceptorOrder:
        prefix, chainsL, chainsR, scopesR, scopesL= lineArray[:5]
      else:
        prefix, chainsR, chainsL, scopesR, scopesL= lineArray[:5]
      prefix= prefix.upper()
      scopesL= set( itertools.chain.from_iterable([elem.split(":") for elem in scopesL.split(";")] ))
      scopesR= set( itertools.chain.from_iterable([elem.split(":") for elem in scopesR.split(";")] ))
      prefixToScope[prefix]= ( scopesL, scopesR )
  return prefixToScope
  
def loadFromJson(fname):
  import json
  prefixToScope={}
  with open(fname) as f:
    dataDict= json.load(f)
    for prefix in dataDict:
      prefix= prefix.upper()
      scopesL= set( itertools.chain.from_iterable([ dataDict[prefix]["l"][chain] for chain in dataDict[prefix]["l"] ] ))
      scopesR= set( itertools.chain.from_iterable([ dataDict[prefix]["r"][chain] for chain in dataDict[prefix]["r"] ] ))
      prefixToScope[prefix]= ( scopesL, scopesR )
  return prefixToScope
  

if __name__=="__main__":
  import sys, json
  scopeFname = sys.argv[1]
  trainTestOutName = sys.argv[2]
  nGroupsToSample=-1
  if len(sys.argv)==4:
    nGroupsToSample= int(sys.argv[3])

  testPrefixes=[]
  trainPrefixes=[]
  with open(scopeFname) as f:
    for line in f:
      lineArray= line.split()
      if len(lineArray)>0:
        if lineArray[0].islower():
          prefix = lineArray[0].strip() + "-" + lineArray[1] + lineArray[2]
          trainPrefixes.append( prefix )
        else:
          prefix = lineArray[0].strip() # + "-" + lineArray[1] + lineArray[2]
          testPrefixes.append( prefix )

  scopeGroups=getScopeGroups(trainPrefixes+testPrefixes, scopeFname, ligandReceptorOrder=True, considerPairsInsteadMonomer=True)
  nTrain= len(trainPrefixes)
  trainGroups= scopeGroups[:nTrain]
  testGropus= scopeGroups[nTrain:]
  independentGroups= set(trainGroups).difference(testGropus)
  nItems= [len(trainGroups), len(testGropus), len(independentGroups)]
  print( nItems )

  independentTrainSet={}
  for prefix, group in zip(trainPrefixes, trainGroups):
    if group in independentGroups:
      if group not in independentTrainSet:
        independentTrainSet[group]=[]
      independentTrainSet[group].append(prefix)

  if nGroupsToSample>1:
    samplingIndex= np.random.choice(independentTrainSet.keys(), size= nGroupsToSample,replace=False)
    independentTestSet_sampled= {key:independentTrainSet[key] for key in samplingIndex}
    independentTrainSet=independentTestSet_sampled

  data={"train":[], "test":[prefix for prefix in testPrefixes]}
  for pdbIds in independentTrainSet.values():
    data["train"].append( np.random.choice(pdbIds))

  with open(trainTestOutName, "w") as f:
    json.dump([data], f)