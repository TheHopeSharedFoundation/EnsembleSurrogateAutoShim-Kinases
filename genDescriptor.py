_version = "0.1"
_usage="""
genDescriptor.py [options] <sdLigandFilename>

"""

from rdkit import Chem
from rdkit import Geometry
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
import sys,os,math
from rdkit.RDLogger import logger

# for now we are using the predefined definitions - however, it might be better to 
# reduce to the necessary features (Don/Acc) to improve speed
#fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
#fdefName = ('/home/samudbe1/201211NOV02PDT1222_MakeRDkitFeatGenSimilarToFloAS/AutoShimBaseFeatures.fdef')
fdefName = os.path.join(os.path.dirname(__file__),"AutoShimBaseFeatures.fdef")
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

def countConstraints(mol,recConf,feat,filt):
  # identify what constraint we have to test:
  mC = mol.GetConformer()
  matchId = feat.GetAtomIds()[0]
  # all values in the filters are "anded" so we need to loop over them
  for i in range(2,len(filt),3):
    if filt[i]=="Distance":
      # get the partners
      dist = mC.GetAtomPosition(matchId).Distance(recConf.GetAtomPosition(filt[1][0]))
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
parser.add_option('--ph4desc',
    help='List of ph4 info for scoring Name Initial_Value Type radius x y z ph4_type')
parser.add_option('--ph4F','--ph4FiltFile',
    help='The name of the file with the ph4 definitions')
parser.add_option('--ph4DescOutFile',
    help='The name of the file with the ph4 definitions')

if __name__=='__main__':
  # input option handling
  options,args = parser.parse_args()
  if len(args)!=1:
    parser.error('Please provide an input sd file')
  if options.ph4desc is None:
    parser.error('You need to have a ph4 desc file')

  # output file or stdout
  if options.outF!='-':
    options.outF = file(options.outF,'w+')
  else:
    options.outF = sys.stdout

  # get the ph4 Scoring parameters

  try:
    options.ph4DescOutFile = open(options.ph4DescOutFile,'w')
  except:
    parser.error('Error opening You need to have a ph4DescOutFile file')

  try:
    f = open(options.ph4desc,'r')
  except IOError, err:
    logger.error(err)


  featureOrder=[]  
  splitL = []
  Ph4_Init_value={}
  Ph4_Descriptors={}
  Ph4_location={}
  Ph4_Type={}
  Ph4_radius={}
  Ph4_FeatureType={}

  for line in f:
    line=line.rstrip("\n")
    #(Name,Init_value,Type,radius,x,y,z,Feature_Type)=line.split("\t")
    splitL = line.split(" ")
    print len (splitL)
    print splitL
    if (len (splitL) < 8):
      print "skipping line in options.ph4desc"
      continue
    Name = splitL[0]
    featureOrder.extend([Name])
    Init_value = splitL[1]
    Type = splitL[2]
    radius = splitL[3]
    x = splitL[4]
    y = splitL[5]
    z = splitL[6]
    Feature_Type = splitL[7]


#  data structure for ph4_Descriptors needs to be decided
    if (Ph4_Descriptors.has_key(Name)):
      print "Caution repeating " + Name + "    "

    
    Ph4_Descriptors[Name] = float (Init_value)
    Ph4_Init_value[Name] = float (Init_value)
    #print Ph4_Descriptors[Name]
    Ph4_location[Name]    = Geometry.rdGeometry.Point3D(float(x),float(y),float(z))
    Ph4_Type[Name]        = Type
    Ph4_radius[Name]      = float (radius)
    Ph4_FeatureType[Name] = Feature_Type
    

  print "featureOrder\n" +  str (featureOrder)

#The following lines are what print the feature headers. -- BMS 26OCT2012PDT1922

  #print >>options.ph4DescOutFile, "molName" + "\t" + "confNum" + "\t" + "protein" + "\t" ,
  print >>options.ph4DescOutFile, "CHIRONID" + "\t" + "ID" + "\t" + "Vscore" + "\t" + "Source" + "\t" + "VAtmCnt" + "\t",

  for key in featureOrder:
    print >>options.ph4DescOutFile, key  + "\t",

  print >>options.ph4DescOutFile, ""
  


  #for key in featureOrder:
      #print key
      #print type (key)
      #print type (Ph4_Descriptors)
    #print key + " " + str(Ph4_Descriptors[key]) + " " + str(Ph4_location[key].x) + " " + str(Ph4_location[key].y) + " " + str(Ph4_location[key].z) + " " + str(Ph4_Type[key]) + " " + str(Ph4_radius[key]) + " " + str(Ph4_FeatureType[key])
      #print key Ph4_Descriptors[key] Ph4_location[key] Ph4_Type[key]  Ph4_radius[key]  Ph4_FeatureType[key]
    #print key Ph_location[key]


    
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
      print  >>sys.stderr,"Error> Bad mol - entry %d" % i+1
      continue
   
    #molName = mol.GetProp('NAME')
    #confNum = mol.GetProp('IX')
    #ensMmb = mol.GetProp('EnsembleMember')
    #icmScore = mol.GetProp('ScoreAfterMin')
    title = mol.GetProp('title')
    numHeavy = str (float (float(mol.GetNumAtoms()) / (10)))

    
#The following line is what needs to be modified to match David's feature file format.  The name of the protein ensemble
#member is derived from the path specified by the SDF property <FILE>.  The path is broken into elements of the array
#named confFileDir -- BMS 26OCT2012PDT1849

    print >>options.ph4DescOutFile, title + "_1" + "\t" + numHeavy + "\t",  
  # now we have the mol - assign the types
    assignedFeatures = featFact.GetFeaturesForMol(mol)
    # loop over all features and check if the fulfill the constraints 
    meetsConstraints = False

    for feat in assignedFeatures:
        pass
        #print "\n\n\n"
#Changed feat.GetType() to feat.GetFamily() on the following line -- BMS 02NOV2012PDT1239 
        #print "Handling %s CONF %s on %s %s atomIdx %d atomType %s recType %.6f %.6f %.6f" % (molName,confNum,feat.GetFamily(),feat.GetFamily(),feat.GetAtomIds()[0],mol.GetAtomWithIdx(feat.GetAtomIds()[0]).GetSymbol(),feat.GetPos().x ,feat.GetPos().y  , feat.GetPos().z)
        #for descriptor in ph4_Descriptors:
#The following calculates the features -- BMS 30OCT2012PDT1337
        for key in featureOrder:
          print "\n"
          print key + " " + str(Ph4_Descriptors[key]) + " " + str(Ph4_location[key].x) + " " + str(Ph4_location[key].y) + " " + str(Ph4_location[key].z)  + " " + str(Ph4_Type[key]) + " " + str(Ph4_radius[key]) + " " + str(Ph4_FeatureType[key])
#Changed following from feat.GetType() to feat.GetFamily() -- BMS 02NOV2012PDT1240           
          print str(feat.GetFamily()) + "?=" + str (Ph4_FeatureType[key])
#Changed following from feat.GetType() to feat.GetFamily() -- BMS 02NOV2012PDT1237 
          if ( feat.GetFamily() == Ph4_FeatureType[key]):
#Changed following from feat.GetType() to feat.GetFamily() -- BMS 02NOV2012PDT1238
            print  str (feat.GetFamily()) + " == " + str (Ph4_FeatureType[key]) + " Match"
            if (Ph4_Type[key] == "Counted"):
              # keep track of the number of times the feat occurs with a certain radius.
              print "feat.ifGetPos().Distance(Ph4_location[key] " + str ( feat.GetPos().Distance(Ph4_location[key]))
              print str (feat.GetPos().Distance(Ph4_location[key])) + " <? " +  str(Ph4_radius[key])
              print str (type (feat.GetPos().Distance(Ph4_location[key]))) + " <? " +  str (type(Ph4_radius[key]))
              if (feat.GetPos().Distance(Ph4_location[key]) < Ph4_radius[key]):
                Ph4_Descriptors[key] = Ph4_Descriptors[key] + 1
                print "Ph4_Descriptors[" + key + "] is bigger " + str (Ph4_Descriptors[key]) + "  " + str ( Ph4_radius[key])
              else:
                print "Ph4_Descriptors[" + key + "] is unchanged " + str (Ph4_Descriptors[key])
                continue
            elif (Ph4_Type[key] == "Closest"):
              pass
              print "feat.GetPos().Distance(Ph4_location[key]) " + str ( feat.GetPos().Distance(Ph4_location[key]) )
              print str (feat.GetPos().Distance(Ph4_location[key])) + " <? " +  str(Ph4_radius[key])
              print str (type (feat.GetPos().Distance(Ph4_location[key]))) + " <? " +  str (type(Ph4_radius[key]))
              if (feat.GetPos().Distance(Ph4_location[key]) < Ph4_Descriptors[key] ):
                Ph4_Descriptors[key] = feat.GetPos().Distance(Ph4_location[key])
                print "Ph4_Descriptors[" + key + "] is smaller " + str (Ph4_Descriptors[key])
            else:
              pass
              print "Ph4_Type[key] " + str ( Ph4_Type[key]) + " did not match"
          else:
            pass
#Changed following from feat.GetType() to feat.GetFamily() -- BMS 02NOV2012PDT1240
            print  str (feat.GetFamily()) + " == " + str (Ph4_FeatureType[key]) + " are not the same"
        
    for key in featureOrder:

#The following line prints out the feature values for each compound.  -- BMS 26OCT2012PDT1852

      print >>options.ph4DescOutFile, str (Ph4_Descriptors[key]) + "\t",
      Ph4_Descriptors[key] = Ph4_Init_value[key]
      #print "resetting Ph4_Descriptors[" + key + "] to " + str (Ph4_Init_value[key])
      

    print >>options.ph4DescOutFile, ""
