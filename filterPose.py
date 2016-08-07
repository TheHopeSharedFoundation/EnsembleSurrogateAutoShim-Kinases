_version = "0.1"
_usage="""
filterPose.py [options] <sdLigandFilename>

"""

from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
import sys,os,math
from rdkit.RDLogger import logger

# for now we are using the predefined definitions - however, it might be better to 
# reduce to the necessary features (Don/Acc) to improve speed
#fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
#fdefName = ('/home/sprincl1/autoshim/ligands/autoshimFilter/BaseFeatures.fdef')
fdefName = os.path.join(os.path.dirname(__file__),"BaseFeatures.fdef")

#############
# functions #
#############
def checkConstraints(mol,recConf,feat,filt):
  # identify what constraint we have to test:
  mC = mol.GetConformer()
  matchId = feat.GetAtomIds()[0]
  # all values in the filters are "anded" so we need to loop over them
  for i in range(2,len(filt),3):
    if filt[i]=="Distance":
      # get the partners
      dist = mC.GetAtomPosition(matchId).Distance(recConf.GetAtomPosition(filt[1][0]))
#!!!POTENTIAL ERROR!!! Identified by SamudBe1 on 04AUG2014 18:44 PDT
#!!!POTENTIAL ERROR!!! In the following line, the code should be "dist>=filt[i+1]" AND "dist<filt[i+2]"
      if (dist<=filt[i+1]) or (dist>filt[i+2]):
        return False
    elif filt[i]=="Angle":
      # changed to be able to handle multiple neighbour atoms for rec
      # the first entry in filt is the matching atom - then a list of all others are neighbours
      # and the same is now also added for the ligand (ie multiple neighbour atoms)
      # get the matching atom neighbours - the ligans don't get the H's read in - so we can use
      # that list directly
      fitsRecAngle = False
      fitsLigAngle = False
      # first the receptor angles - there are usually less nieghbours for rec atoms
      for neighIdx in filt[1][1]:
        # here we are looping over all possible combinations and check at the very end if
        # both angles are ok
        # if an angle is already ok we can skip the second calculation
        # get the vectors
        l1 = mC.GetAtomPosition(matchId)-recConf.GetAtomPosition(filt[1][0])
        l2 = recConf.GetAtomPosition(neighIdx)-recConf.GetAtomPosition(filt[1][0]) 
        angle = math.degrees(l1.AngleTo(l2))
        # old version: angle = math.degrees(math.acos(l1.DotProduct(l2)/(l1.Length()*l2.Length()))) 
        if (angle>filt[i+1]) and (angle<filt[i+2]):
          fitsRecAngle = True
          break
      if not fitsRecAngle:
        return False  

      # now we check on the ligands  
      neighbAtmIdx = [a.GetIdx() for a in mol.GetAtomWithIdx(matchId).GetNeighbors()]
      for idx in neighbAtmIdx:
        l1 = mC.GetAtomPosition(idx)-mC.GetAtomPosition(matchId)
        l2 = recConf.GetAtomPosition(filt[1][0])-mC.GetAtomPosition(matchId)
        angle = math.degrees(l1.AngleTo(l2))
        # old version: angle = math.degrees(math.acos(l1.DotProduct(l2)/(l1.Length()*l2.Length())))
        if (angle>filt[i+1]) and (angle<filt[i+2]):
          fitsLigAngle = True
          break
      if not fitsLigAngle:
        return False

    else:
      logger.error("Requesting a constraint that is not defined %s" % filt[2])
      return False
  # we only reach this position if the ligand matches all queries
  return True

#-----------------------------------------------------------------------------

_usage="""
  filterPose.py [optional arguments] sdFileToFilter
"""
  

from optparse import OptionParser
parser=OptionParser(_usage,version='%prog '+_version)
parser.add_option('--outF','--outFile',default='-',
    help='The name of the output file. The default is the console (stdout).')
parser.add_option('--recF','--receptorFile',
    help='The name of the receptor mol2 file.')
parser.add_option('--scoreFilter',default=None,type='float',
    help='Upper boundary for score filter - if not set, no score filtering is used')
parser.add_option('--scoreF','--scoreFile',default='',
    help='The name of the score file used for filtering. Contains Name - Score combinations')
parser.add_option('--ph4F','--ph4FiltFile',
    help='The name of the file with the ph4 definitions')

if __name__=='__main__':
  # input option handling
  options,args = parser.parse_args()
  if len(args)!=1:
    parser.error('Please provide an input sd file')
  if options.recF is None:
    parser.error('You have to specify a receptor file')
  if options.ph4F is None:
    parser.error('You have to specify a ph4 constraint file')

  # output file or stdout
  if options.outF!='-':
    options.outF = file(options.outF,'w+')
  else:
    options.outF = sys.stdout

  upperScore = None
  if options.scoreFilter and options.scoreF == '':
    parser.error('If you want to use score filtering you have to provide a score file') 
  elif options.scoreFilter:
    upperScore = options.scoreFilter

  # get the input data
  # the receptor
  try:
    rec = Chem.MolFromMol2File(options.recF,removeHs=False)
    # setup the conformer for the receptor
    recConf = rec.GetConformer()
  except:
    print ("Problem reading receptor from %s" % options.recF)
    print "Error:", sys.exc_info()[0]
    quit()

  try:
    f = open(options.ph4F,'r')
  except IOError, err:
    logger.error(err)
  ph4Filters=[]
  for line in f:
    line=line.rstrip("\n")
    splitL=line.split(" ")
    # modify the types from str to the relevant types
    splitL[1]=int(splitL[1])
    # I am not going to worry about "Represntation Errors"!
    for i in range(2,len(splitL),3):
      splitL[i+1]=float(splitL[i+1])
      splitL[i+2]=float(splitL[i+2])
    ph4Filters.append(splitL) 
      
  if options.scoreFilter:
    try:
      f = open(options.scoreF,'r')
    except IOError,err:
      #logger.error(err)
      print err
    scores={}
    for line in f:
      molId,score = line.split(" ")
      scores[molId]=float(score)

  # we create the set of filters based on the input - it is actually quite simple - 
  # for each entry the receptor atoms is expanded automatically to its first heavy atom neighbour
  # this tuple is stored in the respective filters
  for filt in ph4Filters:
    # keep the non-H rec atoms
    neigh=[]
    for a in rec.GetAtomWithIdx(filt[1]-1).GetNeighbors():
      if a.GetAtomicNum() != 1:
        neigh.append(a.GetIdx())
 
    # there should be at least one heavy atom left - if not - we throw an error and break
    if len(neigh)<1:
      logger.error("Problem parsing filter %s - abort" % filt)
      quit()
    # assign into original filter set - this should speed up things
    # I am subtracting the original receptor number by 1 to account for the one-off in atomIdx
    filt[1] = [ filt[1]-1 , neigh ]

  # setup the feature Factory
  featFact = ChemicalFeatures.BuildFeatureFactory(fdefName) 

  # here we maybe should switch to an early switch to avoid multiple if statements for 
  # score filtering - however, I don't think this will really make a big difference
  # so I leave it for the sake of readability
  # loop over all ligands and assign the types
  inFile = args[0]
  supplier=Chem.SDMolSupplier(inFile)
  for i,mol in enumerate(supplier):
    if not mol:
      print  >>sys.stderr,"Error> Bad mol - entry %d" %(i+1)
# parenthesis added around i+1 on 31JAN2013PST1706 by SamudBe1
      continue
    if not mol.HasProp('_Name'):
      molName = 'Mol_%d'%(i+1)
    else:
      molName = mol.GetProp('_Name')

    # check the score - we skip the mol if the score is not ok
    if upperScore and scores.get(molName,upperScore+1)>upperScore:
        continue
 
    # now we have the mol - assign the types
    assignedFeatures = featFact.GetFeaturesForMol(mol)
    # loop over all features and check if the fulfill the constraints 
    meetsConstraints = False
    # I got this of google - I like the idea to break out of deeply nested loops using an exception
    try:
      for feat in assignedFeatures:
        for filt in ph4Filters:
          if feat.GetType()==filt[0]:
            #print "Handling %s for %s on %s atomIdx %d atomType %s recType %s" % (filt[0],molName,feat.GetType(),\
            #      feat.GetAtomIds()[0],mol.GetAtomWithIdx(feat.GetAtomIds()[0]).GetSymbol(),rec.GetAtomWithIdx(filt[1][0]).GetSymbol())
            # check the constraints - the filters are ored but the constraints per filter are anded 
            meetsConstraints = checkConstraints(mol,recConf,feat,filt)
          if meetsConstraints==True:
            raise StopIteration()
    except StopIteration:
      pass

    # at this level meetsConstraints is either true or fals - if it is true - write the molecule
    if meetsConstraints==True:
      #print >>options.outF,mol
      #print >>options.outF write(mol)
      options.outF.write(supplier.GetItemText(i))

