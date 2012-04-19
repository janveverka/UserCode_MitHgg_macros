import ROOT
import array
ROOT.gROOT.SetStyle("Plain")
# Get The ROOT files which contain the limit plotting.
# {"Mass":File}
intlumi = str(4763.0)
doratio=False
dofp=False
usegrid=True

if not usegrid:
  ROOT.gROOT.ProcessLine(".L $CMSSW_BASE/src/MitHgg/macros/medianCalc.C++")
  from ROOT import medianCalc


#EXPfile = ROOT.TFile("unbinnedsmnosystexpected.root")
#EXPfile = ROOT.TFile("unbinnedsmabsJun7expected.root")
#OBSfile = ROOT.TFile("unbinnedsmabsJun7observed.root")
#EXPfile = ROOT.TFile("unbinnedsmrelJun7aexpected.root")
#OBSfile = ROOT.TFile("unbinnedsmrelJun7aobserved.root")
#EXPfile = ROOT.TFile("smabsJul7clsobservedexpected.root")
#OBSfile = ROOT.TFile("smabsJul7clsobservedexpected.root")
#BOBSfile = ROOT.TFile("smabsJul7clsobservedexpected.root")
EXPfile = ROOT.TFile("smrelpobs.root")
OBSfile = ROOT.TFile("smrelpobs.root")

BOBSfile = ROOT.TFile("smrelpobscat4.root")
#COBSfile = ROOT.TFile("smabsJul10PromptSmearpl.root")

Obscfiles = []
Obscnames = []
ObscCols = []

allcats = False

if allcats:

  Obscfiles.append(ROOT.TFile("smrelpobscat0.root"))
  Obscfiles.append(ROOT.TFile("smrelpobscat1.root"))
  Obscfiles.append(ROOT.TFile("smrelpobscat2.root"))
  Obscfiles.append(ROOT.TFile("smrelpobscat3.root"))
  Obscfiles.append(ROOT.TFile("smrelpobscat4.root"))

  Obscnames.append("cat0")
  Obscnames.append("cat1")
  Obscnames.append("cat2")
  Obscnames.append("cat3")
  Obscnames.append("cat4 (VBFTag)")


  ObscCols.append(ROOT.kBlue)
  ObscCols.append(ROOT.kGreen)
  ObscCols.append(ROOT.kOrange)
  ObscCols.append(ROOT.kMagenta)
  ObscCols.append(ROOT.kRed)
  
else:
  Obscfiles.append(ROOT.TFile("smrelpobsinclonly.root"))
  Obscfiles.append(ROOT.TFile("smrelpobscat4.root"))
  
  Obscnames.append("cat0-3 (Non-VBFTag)")
  Obscnames.append("cat4 (VBFTag)")  

  ObscCols.append(ROOT.kBlue)
  ObscCols.append(ROOT.kRed)



Method = "ProfileLikelihood"

EXPName = Method+"/expected"+Method
OBSName = Method+"/higgsCombineTest."+Method
EXPmasses = [110
       ,115
       ,120
       ,125
       ,130
       ,135
       ,140]
#OBSmasses = [110,110.25,110.5,110.75,111,111.5,112,112.5,113,113.5,114,114.5,115,115.5,116,116.5,117,117.5,118,118.5,119,119.5,120,120.5,121,121.5,122,122.5,123,123.5,124,124.5,125,125.5,126,126.5,127,127.5,128,128.5,129,129.5,130,130.5,131,131.5,132,132.5,133,133.5,134,134.5,135,135.5,136,136.5,137,137.5,138,138.5,139,139.5,140]
OBSmasses = []
#newmass=110.
#mcount=0
for i in range(0,81):
  OBSmasses.append(110+i*0.5)
  
EXPmasses = OBSmasses
#OBSmasses = EXPmasses
#-------------------------------------------------------------------------

allMasses = array.array('d',[110.0,110.5,111.0,111.5,112.0,112.5,113.0,113.5,114.0,114.5,115.0,115.5,116.0,116.5,117.0,117.5,118.0,118.5,119.0,119.5,120.0,120.5,121.0,121.5,122.0,122.5,123.0,123.5,124.0,124.5,125.0,125.5,126.0,126.5,127.0,127.5,128.0,128.5,129.0,129.5,130.0,130.5,131.0,131.5,132.0,132.5,133.0,133.5,134.0,134.5,135.0,135.5,136.0,136.5,137.0,137.5,138.0,138.5,139.0,139.5,140.0,140.5,141.0,141.5,142.0,142.5,143.0,143.5,144.0,144.5,145.0,145.5,146.0,146.5,147.0,147.5,148.0,148.5,149.0,149.5,150.0])

xSec = array.array('d',[0.04474106,0.04476087,0.04458980,0.04461816,0.04465764,0.04445241,0.04450231,0.04452202,0.04432575,0.04435697,0.04417173,0.04419175,0.04402058,0.04382749,0.04365297,0.04367220,0.04349037,0.04330598,0.04314091,0.04297301,0.04277804,0.04260475,0.04226444,0.04208721,0.04190975,0.04154534,0.04122657,0.04106506,0.04072321,0.04038161,0.04006593,0.03972509,0.03940765,0.03909069,0.03877650,0.03846030,0.03797801,0.03766622,0.03721147,0.03689953,0.03645195,0.03600491,0.03556311,0.03514415,0.03470779,0.03427222,0.03370990,0.03328123,0.03290089,0.03235105,0.03180601,0.03141397,0.03087743,0.03034568,0.02981916,0.02933985,0.02882227,0.02830964,0.02782357,0.02733986,0.02672205,0.02624931,0.02578068,0.02525736,0.02473893,0.02414999,0.02356720,0.02307409,0.02258560,0.02202963,0.02147946,0.02101546,0.02055579,0.02003015,0.01950998,0.01893346,0.01836331,0.01786838,0.01737859,0.01683392,0.01629523])

#xSecErrPlus = array.array('d',[23.8874,23.6736,23.46,23.2465,23.0332,22.82,22.6069,22.394,22.1812,21.9685,21.756,21.5706,21.3853,21.2001,21.015,20.8299,20.645,20.4601,20.2754,20.0907,19.9061,19.7448,19.5836,19.4224,19.2612,19.1001,18.9391,18.7781,18.6172,18.4563,18.2955,18.1487,18.002,17.8554,17.7088,17.5624,17.4159,17.2696,17.1234,16.9772,16.831,16.7029,16.5748,16.4467,16.3187,16.1908,16.0629,15.9351,15.8074,15.6797,15.5521,15.4379,15.3236,15.2094,15.0952,14.981,14.8669,14.7527,14.6386,14.5245,14.5245])

#xSecErrMinus = array.array('d',[16.8045,16.6596,16.5148,16.37,16.2251,16.0803,15.9355,15.7906,15.6458,15.5009,15.3561,15.2327,15.1091,14.9856,14.8619,14.7382,14.6145,14.4907,14.3668,14.2429,14.1189,14.0068,13.8947,13.7827,13.6706,13.5585,13.4465,13.3344,13.2223,13.1103,12.9982,12.8972,12.7961,12.6951,12.5941,12.493,12.392,12.291,12.1899,12.0889,11.9879,11.901,11.8141,11.7271,11.6402,11.5532,11.4662,11.3792,11.2921,11.2051,11.118,11.0385,10.9591,10.8796,10.8001,10.7206,10.641,10.5614,10.4819,10.4023,10.4023])
xSecErrPlus = array.array('d',[0.0]*len(allMasses))
xSecErrMinus = array.array('d',[0.0]*len(allMasses))
xSecErrPlusPercent = array.array('d',[0.0]*len(allMasses))
xSecErrMinusPercent = array.array('d',[0.0]*len(allMasses))

#xSecErrPlusPercent = array.array('d',[20.4,20.36,20.32,20.28,20.24,20.2,20.16,20.12,20.08,20.04,20,19.97,19.94,19.91,19.88,19.85,19.82,19.79,19.76,19.73,19.7,19.68,19.66,19.64,19.62,19.6,19.58,19.56,19.54,19.52,19.5,19.47,19.44,19.41,19.38,19.35,19.32,19.29,19.26,19.23,19.2,19.17,19.14,19.11,19.08,19.05,19.02,18.99,18.96,18.93,18.9,18.89,18.88,18.87,18.86,18.85,18.84,18.83,18.82,18.81,18.81])

#xSecErrMinusPercent = array.array('d',[15.3,15.3,15.3,15.3,15.3,15.3,15.3,15.3,15.3,15.3,15.3,15.28,15.26,15.24,15.22,15.2,15.18,15.16,15.14,15.12,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.09,15.08,15.07,15.06,15.05,15.04,15.03,15.02,15.01,15,14.99,14.98,14.97,14.96,14.95,14.94,14.93,14.92,14.91,14.91])

br = array.array('d',[0.00197,0.00198737,0.00200441,0.00202113,0.00203754,0.00205365,0.00206948,0.00208501,0.00210028,0.00211527,0.00213,0.00214331,0.0021563,0.00216899,0.00218138,0.00219348,0.0022053,0.00221685,0.00222815,0.0022392,0.00225,0.00225571,0.00226125,0.00226662,0.00227182,0.00227687,0.00228177,0.00228652,0.00229114,0.00229563,0.0023,0.00229526,0.00229072,0.00228635,0.00228215,0.00227811,0.00227422,0.00227047,0.00226686,0.00226337,0.00226,0.00224526,0.00223124,0.00221791,0.0022052,0.00219308,0.00218151,0.00217044,0.00215986,0.00214972,0.00214,0.00211438,0.00209031,0.00206765,0.00204629,0.00202612,0.00200703,0.00198895,0.0019718,0.0019555,0.00193])

if dofp:
  tmpallMasses = array.array('d',[110,115,120,125,130,135,140])
  tmpxSec = array.array('d',[44.49/1000.*3.567,43.94/1000.*1.987,42.56/1000.*1.198,40.04/1000.*0.783,36.32/1000.*0.545,31.80/1000.*0.406,26.74/1000.*0.320])

  for j in range(0,6):
    for i in range(0,10):
      xSec[10*j+i] = tmpxSec[j] + (tmpxSec[j+1]-tmpxSec[j])*float(i)/10.0
  xSec[60]=tmpxSec[6]
  
  for j in range(0,81):
    xSecErrPlusPercent[j]=0
    xSecErrMinusPercent[j]=0
    #br[j]=1.0
  

  xSecErrPlus = xSec
  xSecErrMinus = xSec

  #xSecErrPlusPercent = array.array('d',[0.0,0.0,0.0,0.0,0.0,0.0,0.0])

  #xSecErrMinusPercent = array.array('d',[0.0,0.0,0.0,0.0,0.0,0.0,0.0])

  #br = array.array('d',[1.0,1.0,1.0,1.0,1.0,1.0,1.0])

for j in range(0,81):
    xSecErrPlusPercent[j]=0
    xSecErrMinusPercent[j]=0

## Nick's python for limits plus additional stuff to plot SM lines
ROOT.gROOT.ProcessLine( \
   "struct Entry{	\
    double r;		\
   };"
)
from ROOT import Entry
def getOBSERVED(file,m,quantile):
  tree = file.Get("limit")
  #br = tree.GetBranch("limit")
  #c = Entry()
  #br.SetAddress(ROOT.AddressOf(c,'r'))
  ientries = tree.GetEntries()
  for ientry in xrange(ientries):
    tree.GetEntry(ientry)
    if (abs(tree.mh-m)<0.001 and abs(tree.quantileExpected-quantile)<0.01): return tree.limit
    #if (abs(tree.mh-m)<0.001): return tree.limit

  #print "recovering"

  return 0.0
  mindiff=999.
  limit=0.
  for ientry in xrange(ientries):
    tree.GetEntry(ientry)
   # print abs(tree.mh-m)
   # print mindiff
    if (abs(tree.mh-m)<mindiff and abs(tree.quantileExpected-quantile)<0.01 and tree.limit>0):
    #if (abs(tree.mh-m)<mindiff and tree.limit>0):
      mindiff=abs(tree.mh-m)
      limit=tree.limit
    #  print limit
  return limit
  

def getObservedGraph(file,quantile):
  
  sobs = [getOBSERVED(file,m,-1.0) for m in OBSmasses]

  isobs=0
  
  graphSObs = ROOT.TGraphAsymmErrors()
  
  for i,mass in zip(range(len(OBSmasses)),OBSmasses):
    if sobs[i]>0.0:
      graphSObs.SetPoint(isobs,float(mass),sobs[i])
      graphSObs.SetPointError(isobs,0,0,0,0)
      isobs+=1 

  return graphSObs

#EXPfiles = [ROOT.TFile(EXPName+".mH%d.root"%m) for m in EXPmasses]
#OBSfiles = [ROOT.TFile(OBSName+".mH%d.root"%m) for m in OBSmasses]


obs = [getOBSERVED(OBSfile,m,-1.0) for m in OBSmasses]
eobs = [getOBSERVED(EXPfile,m,-1.0) for m in OBSmasses]
bobs = [getOBSERVED(BOBSfile,m,-1.0) for m in OBSmasses]


print "observed limits:"
print obs

leg=ROOT.TLegend(0.42,0.12,0.87,0.35)
#leg=ROOT.TLegend(0.52,0.10,0.87,0.25)
leg.SetFillColor(0)
leg.SetBorderSize(0)

graph68  = ROOT.TGraphAsymmErrors()
graph95  = ROOT.TGraphAsymmErrors()
graphMed = ROOT.TGraphAsymmErrors()
graphObs = ROOT.TGraphAsymmErrors()
graphBObs = ROOT.TGraphAsymmErrors()
graphCObs = ROOT.TGraphAsymmErrors()
graphEObs = ROOT.TGraphAsymmErrors()

graphOne = ROOT.TGraphAsymmErrors()

ROOT.gStyle.SetEndErrorSize(7)


#0.989458 +/- 0.00198175
graphT = ROOT.TGraphErrors()
graphT.SetPoint(0,125,1.0-0.997975)
graphT.SetPointError(0,0.,0.000764777)

#graphT.SetPoint(0,123.5,1.0-0.989458)
#graphT.SetPointError(0,0.,0.00198)

graphT.SetMarkerSize(2)
graphT.SetMarkerStyle(8)
graphT.SetMarkerColor(9)
graphT.SetLineWidth(2)


leg.AddEntry(graphObs,"Observed Asymptotic","L")
leg.AddEntry(graphT,"Observed Ensemble","LPE")
#leg.AddEntry(graphBObs,"Observed (Fit 100-150 GeV)","L")


#leg.AddEntry(graphEObs,"1xSM Higgs Median Expected","L")
#leg.AddEntry(graphBObs,"1xSM Higgs Single Mass 123.5 GeV","L")

#leg.AddEntry(graphEObs,"1xFP Higgs Median Expected","L")
#leg.AddEntry(graphBObs,"1xFP Higgs Single Mass 113.5 GeV","L")

#leg.AddEntry(graphBObs,"May10+Prompt CL","L")
#leg.AddEntry(graphCObs,"05July ReReco (Prompt Res) PL","L")

#leg.AddEntry(graphMed,"Expected CLs Limit","L")
#leg.AddEntry(graph68,"#pm 1#sigma Expected CLs","F")
#leg.AddEntry(graph95,"#pm 2#sigma Expected CLs","F")

MG = ROOT.TMultiGraph()


sm = 1
#EXPECTED
iexp=0
for i,mass in zip(range(len(EXPmasses)),EXPmasses):
  


  for j,mm in enumerate(allMasses): 
	if mm==mass: sm = 1.0/(xSec[j])
  
  if not doratio: sm=1.0
  
  tree = EXPfile.Get("limit")
  
  median = array.array('d',[0])
  up68   = array.array('d',[0])
  dn68   = array.array('d',[0])
  up95   = array.array('d',[0])
  dn95   = array.array('d',[0])
  
  if usegrid:
    median[0] = getOBSERVED(EXPfile,mass,0.5)
    dn68[0] = getOBSERVED(EXPfile,mass,0.16)
    up68[0] = getOBSERVED(EXPfile,mass,0.84)
    dn95[0] = getOBSERVED(EXPfile,mass,0.025)
    up95[0] = getOBSERVED(EXPfile,mass,0.975)
    #print median[0]
    #print up68[0]
    #print dn68[0]
    #print up95[0]
    #print dn95[0]
  else:
    medianCalc("r_mH"+str(mass),tree,median,up68,dn68,up95,dn95,mass)


  print "mass = %d, median = %f, dn68 = %f, up68 = %f, dn95 = %f, up95 = %f" % (mass,median[0],dn68[0],up68[0],dn95[0],up95[0])

  #print "median"
  #print median
  #print "dn68"
  #print dn68
  #print "up68"
  #print up68
  #print "dn95"
  #print dn95
  #print "up95"
  #print up95
    
    
  if (median[0]>0.0):
    graph68.SetPoint(iexp,float(mass),median[0]*sm)
    graph95.SetPoint(iexp,float(mass),median[0]*sm)
    graphMed.SetPoint(iexp,float(mass),median[0]*sm)
    graphOne.SetPoint(iexp,float(mass),1.*sm)

    diff95_up = abs(median[0] - up95[0])*sm
    diff95_dn = abs(median[0] - dn95[0])*sm
    diff68_up = abs(median[0] - up68[0])*sm 
    diff68_dn = abs(median[0] - dn68[0])*sm 

    
    graph68.SetPointError(iexp,0,0,diff68_dn,diff68_up)
    graph95.SetPointError(iexp,0,0,diff95_dn,diff95_up)
    graphMed.SetPointError(iexp,0,0,0,0)
    graphOne.SetPointError(iexp,0,0,0,0)
    iexp+=1

#OBSERVED
iobs=0
ibobs=0
icobs=0
ieobs=0
for i,mass in zip(range(len(OBSmasses)),OBSmasses):

  for j,mm in enumerate(allMasses): 
	if mm==mass: sm = 1.0/(xSec[j])

  if not doratio: sm=1.0

  if obs[i]>0.0:
    graphObs.SetPoint(iobs,float(mass),obs[i]*sm)
    graphObs.SetPointError(iobs,0,0,0,0)
    iobs+=1
    
  if bobs[i]>0.0:
    graphBObs.SetPoint(ibobs,float(mass),bobs[i]*sm)
    graphBObs.SetPointError(ibobs,0,0,0,0)
    ibobs+=1 

  if eobs[i]>0.0:
    graphEObs.SetPoint(ieobs,float(mass),eobs[i]*sm)
    graphEObs.SetPointError(ieobs,0,0,0,0)
    ieobs+=1   

#-------------------------------------------------------------------------
xSec2 = array.array('d',[0.0]*len(allMasses))
xSec5 = array.array('d',[0.0]*len(allMasses))
xSec10 = array.array('d',[0.0]*len(allMasses))
xSec20 = array.array('d',[0.0]*len(allMasses))
xSec30 = array.array('d',[0.0]*len(allMasses))
xSec40 = array.array('d',[0.0]*len(allMasses))

xSecErrPlus2 = array.array('d',[0.0]*len(allMasses))
xSecErrPlus5 = array.array('d',[0.0]*len(allMasses))
xSecErrPlus10 = array.array('d',[0.0]*len(allMasses))
xSecErrPlus20 = array.array('d',[0.0]*len(allMasses))
xSecErrPlus30 = array.array('d',[0.0]*len(allMasses))
xSecErrPlus40 = array.array('d',[0.0]*len(allMasses))
xSecErrMinus2 = array.array('d',[0.0]*len(allMasses))
xSecErrMinus5 = array.array('d',[0.0]*len(allMasses))
xSecErrMinus10 = array.array('d',[0.0]*len(allMasses))
xSecErrMinus20 = array.array('d',[0.0]*len(allMasses))
xSecErrMinus30 = array.array('d',[0.0]*len(allMasses))
xSecErrMinus40 = array.array('d',[0.0]*len(allMasses))
dPlus = array.array('d',[0.0]*len(allMasses))
dPlus2 = array.array('d',[0.0]*len(allMasses))
dPlus5 = array.array('d',[0.0]*len(allMasses))
dPlus10 = array.array('d',[0.0]*len(allMasses))
dPlus20 = array.array('d',[0.0]*len(allMasses))
dPlus30 = array.array('d',[0.0]*len(allMasses))
dPlus40 = array.array('d',[0.0]*len(allMasses))
dMinus = array.array('d',[0.0]*len(allMasses))
dMinus2 = array.array('d',[0.0]*len(allMasses))
dMinus5 = array.array('d',[0.0]*len(allMasses))
dMinus10 = array.array('d',[0.0]*len(allMasses))
dMinus20 = array.array('d',[0.0]*len(allMasses))
dMinus30 = array.array('d',[0.0]*len(allMasses))
dMinus40 = array.array('d',[0.0]*len(allMasses))



for i in range(len(allMasses)):
  
  xSecErrPlus[i] =  xSec[i] +  xSec[i]*xSecErrPlusPercent[i]/100.;
  xSecErrMinus[i] = xSec[i] -  xSec[i]*xSecErrMinusPercent[i]/100.;  
  xSec[i]=xSec[i];
  xSecErrPlus[i] = xSecErrPlus[i];
  xSecErrMinus[i]= xSecErrMinus[i];
  xSec2[i]=2.*xSec[i];
  xSec5[i]=5.*xSec[i];
  xSec10[i]=10.*xSec[i];
  xSec20[i]=20.*xSec[i];
  xSec30[i]=30.*xSec[i];
  xSec40[i]=40.*xSec[i];
  
  xSecErrPlus2[i]=2.*xSecErrPlus[i]
  xSecErrPlus5[i]=5.*xSecErrPlus[i]
  xSecErrPlus10[i]=10.*xSecErrPlus[i]
  xSecErrPlus20[i]=20.*xSecErrPlus[i]
  xSecErrPlus30[i]=30.*xSecErrPlus[i]
  xSecErrPlus40[i]=40.*xSecErrPlus[i]
  xSecErrMinus2[i]=2.*xSecErrMinus[i]
  xSecErrMinus5[i]=5.*xSecErrMinus[i]
  xSecErrMinus10[i]=10.*xSecErrMinus[i]
  xSecErrMinus20[i]=20.*xSecErrMinus[i]
  xSecErrMinus30[i]=30.*xSecErrMinus[i]
  xSecErrMinus40[i]=40.*xSecErrMinus[i]
  dPlus[i]  =  abs(xSecErrPlus[i]-xSec[i])
  dMinus[i] =  abs(xSec[i] - xSecErrMinus[i])
  dPlus2[i]  =  abs(xSecErrPlus2[i]-xSec2[i])
  dMinus2[i] =  abs(xSec2[i] - xSecErrMinus2[i])
  dPlus5[i]  =  abs(xSecErrPlus5[i]-xSec5[i])
  dMinus5[i] =  abs(xSec5[i] - xSecErrMinus5[i])
  dPlus10[i]  =  abs(xSecErrPlus10[i]-xSec10[i])
  dMinus10[i] =  abs(xSec10[i] - xSecErrMinus10[i])
#  print "debug %f %f %f" %(,xSec10[i],dPlus10,dMinus10)  
  dPlus20[i]  =  xSecErrPlus20[i]-xSec20[i]
  dMinus20[i] =  xSec20[i] - xSecErrMinus20[i]
  dPlus30[i]  =  xSecErrPlus30[i]-xSec30[i]
  dMinus30[i] =  xSec30[i] - xSecErrMinus30[i]
  dPlus40[i]  =  xSecErrPlus40[i]-xSec40[i]
  dMinus40[i] =  xSec40[i] - xSecErrMinus40[i]


  
  myGraphXSecSM   = ROOT.TGraphAsymmErrors()
  myGraphXSec2SM = ROOT.TGraphAsymmErrors()
  myGraphXSec5SM = ROOT.TGraphAsymmErrors()
  myGraphXSec10SM = ROOT.TGraphAsymmErrors()
  myGraphXSec20SM = ROOT.TGraphAsymmErrors()
  myGraphXSec30SM = ROOT.TGraphAsymmErrors()
  myGraphXSec40SM = ROOT.TGraphAsymmErrors()
  
  
  
  for i,mass in zip(range(len(allMasses)),allMasses):
    for j,mm in enumerate(allMasses): 
      if mm==mass: sm = xSec[j]
  
    if not doratio: sm=1.0

    myGraphXSecSM.SetPoint(i,allMasses[i],xSec[i]/sm)
    myGraphXSecSM.SetPointError(i,0,0,dMinus[i]/sm,dPlus[i]/sm)
    myGraphXSec2SM.SetPoint(i,allMasses[i],xSec2[i]/sm)
    myGraphXSec2SM.SetPointError(i,0,0,  dMinus2[i]/sm,dPlus2[i]/sm)
    myGraphXSec5SM.SetPoint(i,allMasses[i],xSec5[i]/sm)
    myGraphXSec5SM.SetPointError(i,0,0,  dMinus5[i]/sm,dPlus5[i]/sm)
    myGraphXSec10SM.SetPoint(i,allMasses[i],xSec10[i]/sm)
    myGraphXSec10SM.SetPointError(i,0,0,  dMinus10[i]/sm,dPlus10[i]/sm)
    myGraphXSec20SM.SetPoint(i,allMasses[i],xSec20[i]/sm)
    myGraphXSec20SM.SetPointError(i,0,0,  dMinus20[i]/sm,dPlus20[i]/sm)
    myGraphXSec30SM.SetPoint(i,allMasses[i],xSec30[i]/sm)
    myGraphXSec30SM.SetPointError(i,0,0,  dMinus30[i]/sm,dPlus30[i]/sm)
    myGraphXSec40SM.SetPoint(i,allMasses[i],xSec40[i]/sm)
    myGraphXSec40SM.SetPointError(i,0,0,  dMinus40[i]/sm,dPlus40[i]/sm)






#---------------------------------------------------------------------------
#graph95.SetFillColor(ROOT.kYellow-4)
graph95.SetFillColor(ROOT.kYellow-4)
#graph95.SetFillStyle(3001)
#graph68.SetFillColor(ROOT.kGreen+2)
graph68.SetFillColor(ROOT.kGreen)
#graph68.SetFillStyle(3001)
graphMed.SetLineStyle(2)
graphMed.SetLineColor(2)
graphMed.SetLineWidth(3)
#graphObs.SetMarkerStyle(2)

graphOne.SetLineWidth(3)

xsecbandfill=3244
#xsecbandfill=0

myGraphXSecSM.SetLineStyle(2)
myGraphXSecSM.SetLineColor(ROOT.kAzure+7)
myGraphXSecSM.SetLineWidth(4)
myGraphXSecSM.SetFillColor(ROOT.kAzure+7)
myGraphXSecSM.SetFillStyle(xsecbandfill)
myGraphXSec2SM.SetLineStyle(2)
myGraphXSec2SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec2SM.SetLineWidth(4)
myGraphXSec2SM.SetFillColor(ROOT.kAzure+7)
myGraphXSec2SM.SetFillStyle(xsecbandfill)
myGraphXSec5SM.SetLineStyle(2)
myGraphXSec5SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec5SM.SetLineWidth(4)
myGraphXSec5SM.SetFillColor(ROOT.kAzure+7)
myGraphXSec5SM.SetFillStyle(xsecbandfill)
myGraphXSec10SM.SetLineStyle(2)
myGraphXSec10SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec10SM.SetLineWidth(4)
myGraphXSec10SM.SetFillColor(ROOT.kAzure+7)
myGraphXSec10SM.SetFillStyle(xsecbandfill)
myGraphXSec20SM.SetLineStyle(2)
myGraphXSec20SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec20SM.SetLineWidth(4)
myGraphXSec20SM.SetFillColor(ROOT.kAzure+7)
myGraphXSec20SM.SetFillStyle(xsecbandfill)
myGraphXSec30SM.SetLineStyle(2)
myGraphXSec30SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec30SM.SetLineWidth(4)
myGraphXSec30SM.SetFillColor(ROOT.kAzure+7)
myGraphXSec30SM.SetFillStyle(xsecbandfill)
myGraphXSec40SM.SetLineStyle(2)
myGraphXSec40SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec40SM.SetLineWidth(4)
myGraphXSec40SM.SetFillColor(ROOT.kAzure+7)
myGraphXSec40SM.SetFillStyle(xsecbandfill)



# use 95 as the default guy

#MG.Add(graph95)
#MG.Add(graph68)
##MG.Add(graphOne)
#MG.Add(graphMed)
#MG.Add(myGraphXSecSM)
#if not doratio:
  #if dofp: MG.Add(myGraphXSec2SM)
  #MG.Add(myGraphXSec20SM)
  #MG.Add(myGraphXSec30SM)
#MG.Add(myGraphXSec5SM)
#MG.Add(myGraphXSec10SM)
  
#MG.Add(myGraphXSec40SM)
#graphObs.SetMarkerStyle(22)
#graphObs.SetMarkerSize(1.5)

graphBObs.SetLineStyle(7)
graphBObs.SetLineColor(ROOT.kRed)

#graphBObs.SetMarkerStyle(20)
#graphBObs.SetMarkerSize(1.5)

graphEObs.SetLineStyle(7)
graphEObs.SetLineWidth(4)
graphEObs.SetLineColor(ROOT.kBlue)


graphObs.SetLineWidth(3)
graphBObs.SetLineWidth(3)

MG.Add(graphObs)
#MG.Add(graphBObs)
#MG.Add(graphCObs)
#MG.Add(graphEObs)

for i in range(len(Obscfiles)):
  graph = getObservedGraph(Obscfiles[i],-1.0)
  graph.SetLineColor(ObscCols[i])
  graph.SetLineWidth(2)
  MG.Add(graph)
  leg.AddEntry(graph,Obscnames[i],"L")


MG.Add(graphT,"PE1")

# -------------------------------------
C = ROOT.TCanvas("#int L = %s"%intlumi,"#int L = %s"%intlumi,1600,1100)
#C.SetLogy()
C.SetGrid(True)
#MG.Draw("AL3")
#if not dofp:
dummy = ROOT.TH1D("dummy","",1,110.,150.);
dummy.SetStats(False)

dummy.GetYaxis().SetTitle("p-value")

dummy.Draw()
MG.Draw("L3")
dummy.GetXaxis().SetTitle("m_{H} (GeV/c^{2})")


dummy.SetMinimum(0.0001)
dummy.SetMaximum(3.0)
#maxy=170
C.SetLogy()



onesigma = 1.58655253931457074e-01
twosigma = 2.27501319481792155e-02
threesigma = 1.34989803163009588e-03
foursigma = 3.16712418331199785e-05


onesigmaline = ROOT.TLine(110.0,onesigma,150.0,onesigma)
onesigmaline.SetLineColor(ROOT.kRed)
onesigmaline.SetLineWidth(2)
onesigmaline.Draw()
onesigmalabel = ROOT.TLatex(150.5,onesigma,"1\sigma")
onesigmalabel.Draw()

twosigmaline = ROOT.TLine(110.0,twosigma,150.0,twosigma)
twosigmaline.SetLineColor(ROOT.kRed)
twosigmaline.SetLineWidth(2)
twosigmaline.Draw()
twosigmalabel = ROOT.TLatex(150.5,twosigma,"2\sigma")
twosigmalabel.Draw()

threesigmaline = ROOT.TLine(110.0,threesigma,150.0,threesigma)
threesigmaline.SetLineColor(ROOT.kRed)
threesigmaline.SetLineWidth(2)
threesigmaline.Draw()
threesigmalabel = ROOT.TLatex(150.5,threesigma,"3\sigma")
threesigmalabel.Draw()

#foursigmaline = ROOT.TLine(110.0,foursigma,150.0,foursigma)
#foursigmaline.SetLineColor(ROOT.kRed)
#foursigmaline.SetLineWidth(2)
#foursigmaline.Draw()
#foursigmalabel = ROOT.TLatex(150.5,foursigma,"4\sigma")
#foursigmalabel.Draw()

leg.Draw()

#box=ROOT.TLegend(0.11,0.75,0.89,0.89)
box=ROOT.TLegend(0.11,0.82,0.89,0.89)
box.SetFillColor(ROOT.kWhite)
box.SetFillStyle(1001)
box.SetBorderSize(0)
box.Draw()

mytext = ROOT.TLatex()
mytext.SetTextSize(0.04)
mytext.DrawLatex(133,1.55,"#splitline{CMS preliminary}{#sqrt{s} = 7 TeV L = 4.76 fb^{-1}}")
#mytext.DrawLatex(111,1.75,"#scale[1.0]{Trials Factor of ~20 Not Included}")
#mytext.DrawLatex(111,1.3,"#scale[1.0]{Local  Significance: 2.3#sigma}")
#mytext.DrawLatex(111,0.9,"#scale[1.0]{Global Significance: 0.8#sigma}")
mytext.DrawLatex(111,1.75,"#scale[1.0]{Interpretation Requires LEE}")


graphObs.Print()
graphBObs.Print()

#C.SaveAs("LimitBayesianFFRel.eps")
C.SaveAs("pvaluesmcat.pdf")
C.SaveAs("pvaluesmcat.eps")
#C.SaveAs("test.eps")

