import ROOT
import array
ROOT.gROOT.SetStyle("Plain")
# Get The ROOT files which contain the limit plotting.
# {"Mass":File}
intlumi = str(4763.0)
doratio=True
dofp=False
usegrid=True
divxsec = False

if not usegrid:
  ROOT.gROOT.ProcessLine(".L $CMSSW_BASE/src/MitHgg/macros/medianCalc.C++")
  from ROOT import medianCalc

BOBSfile = ROOT.TFile("smrel.root")
EXPfile = ROOT.TFile("smrel.root")
OBSfile = ROOT.TFile("smrel.root")

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

       OBSmasses = []
for i in range(0,81):
  OBSmasses.append(110+i*0.5)
  
EXPmasses = OBSmasses

allMasses = array.array('d',[110.0,110.5,111.0,111.5,112.0,112.5,113.0,113.5,114.0,114.5,115.0,115.5,116.0,116.5,117.0,117.5,118.0,118.5,119.0,119.5,120.0,120.5,121.0,121.5,122.0,122.5,123.0,123.5,124.0,124.5,125.0,125.5,126.0,126.5,127.0,127.5,128.0,128.5,129.0,129.5,130.0,130.5,131.0,131.5,132.0,132.5,133.0,133.5,134.0,134.5,135.0,135.5,136.0,136.5,137.0,137.5,138.0,138.5,139.0,139.5,140.0,140.5,141.0,141.5,142.0,142.5,143.0,143.5,144.0,144.5,145.0,145.5,146.0,146.5,147.0,147.5,148.0,148.5,149.0,149.5,150.0])
xSec = array.array('d',[0.04474106,0.04476087,0.04458980,0.04461816,0.04465764,0.04445241,0.04450231,0.04452202,0.04432575,0.04435697,0.04417173,0.04419175,0.04402058,0.04382749,0.04365297,0.04367220,0.04349037,0.04330598,0.04314091,0.04297301,0.04277804,0.04260475,0.04226444,0.04208721,0.04190975,0.04154534,0.04122657,0.04106506,0.04072321,0.04038161,0.04006593,0.03972509,0.03940765,0.03909069,0.03877650,0.03846030,0.03797801,0.03766622,0.03721147,0.03689953,0.03645195,0.03600491,0.03556311,0.03514415,0.03470779,0.03427222,0.03370990,0.03328123,0.03290089,0.03235105,0.03180601,0.03141397,0.03087743,0.03034568,0.02981916,0.02933985,0.02882227,0.02830964,0.02782357,0.02733986,0.02672205,0.02624931,0.02578068,0.02525736,0.02473893,0.02414999,0.02356720,0.02307409,0.02258560,0.02202963,0.02147946,0.02101546,0.02055579,0.02003015,0.01950998,0.01893346,0.01836331,0.01786838,0.01737859,0.01683392,0.01629523])
ffxsec = array.array('d',[0.16346707,0.153639388,0.144474184,0.135978954,0.128011936,0.120580658,0.113698278,0.1072188,0.101204298,0.095574901,0.090277395,0.085366199,0.080723132,0.076390615,0.072348276,0.068536512,0.064981344,0.061626793,0.058491355,0.055518848,0.052712643,0.050108436,0.047635539,0.045289044,0.043111872,0.041038726,0.0390816,0.03724344,0.035498398,0.033842835,0.032293386,0.03081212,0.029395632,0.028062999,0.026825547,0.025612992,0.024479542,0.023414238,0.02239302,0.021406322,0.020504886,0.01962496,0.018787808,0.0179886462,0.0172365732,0.0165092732,0.0158259735,0.0151657605,0.014543112,0.0139496598,0.01338225,0.0128386206,0.012320084,0.0118226735,0.0113473074,0.010898992,0.0104629713,0.010045404,0.0096510688,0.0092662975,0.0088983973,0.0085471461,0.0082114608,0.007885017,0.007570577,0.0072704658,0.0069803288,0.0067006628,0.0064322028,0.00617446755,0.005927448,0.005686585,0.00545616,0.005232885,0.0050181952,0.00480948195,0.004609108,0.0044133128,0.0042256148,0.0040421231,0.0038650707])

xSecErrPlus = array.array('d',[0.0]*len(allMasses))
xSecErrMinus = array.array('d',[0.0]*len(allMasses))
xSecErrPlusPercent = array.array('d',[0.0]*len(allMasses))
xSecErrMinusPercent = array.array('d',[0.0]*len(allMasses))


if dofp:  
  xSec = ffxsec
  
  xSecErrPlus = xSec
  xSecErrMinus = xSec



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
    #if (abs(tree.mh-m)<0.001 and abs(tree.quantileExpected-quantile)<0.02) and tree.limit>0.005 and tree.limit<0.2: return [tree.limit,tree.limitErr]
    if (abs(tree.mh-m)<0.001 and abs(tree.quantileExpected-quantile)<0.04): return [tree.limit,tree.limitErr]
    #if (abs(tree.mh-m)<0.001): return [tree.limit,tree.limitErr]

  #print "recovering"

  return [0.0,0.0]
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
  


#EXPfiles = [ROOT.TFile(EXPName+".mH%d.root"%m) for m in EXPmasses]
#OBSfiles = [ROOT.TFile(OBSName+".mH%d.root"%m) for m in OBSmasses]


obs = [getOBSERVED(OBSfile,m,-1.0) for m in OBSmasses]
bobs = [getOBSERVED(BOBSfile,m,0.5) for m in OBSmasses]
bobsmedian = [getOBSERVED(BOBSfile,m,0.5) for m in OBSmasses]
bobs68dn = [getOBSERVED(BOBSfile,m,0.16) for m in OBSmasses]
bobs68up = [getOBSERVED(BOBSfile,m,0.84) for m in OBSmasses]
bobs95dn = [getOBSERVED(BOBSfile,m,0.025) for m in OBSmasses]
bobs95up = [getOBSERVED(BOBSfile,m,0.975) for m in OBSmasses]

print "observed limits:"
print obs

#leg=ROOT.TLegend(0.46,0.56,0.79,0.89)
leg=ROOT.TLegend(0.12,0.70,0.58,0.89)
#leg=ROOT.TLegend(0.14,0.50,0.49,0.69)
leg.SetFillColor(0)
leg.SetBorderSize(0)

graph68  = ROOT.TGraphAsymmErrors()
graph95  = ROOT.TGraphAsymmErrors()
#graphMed = ROOT.TGraphAsymmErrors()
graphMed = ROOT.TGraph()
graphObs = ROOT.TGraphAsymmErrors()
graphBObs = ROOT.TGraphAsymmErrors()
graphOne = ROOT.TGraphAsymmErrors()

graph68up = ROOT.TGraphErrors()
graph68dn = ROOT.TGraphErrors()
graph95up = ROOT.TGraphErrors()
graph95dn = ROOT.TGraphErrors()
graphmede = ROOT.TGraphErrors()

graphb68up = ROOT.TGraph()
graphb68dn = ROOT.TGraph()
graphb95up = ROOT.TGraph()
graphb95dn = ROOT.TGraph()
graphbmede = ROOT.TGraph()



leg.AddEntry(graphObs,"Observed","L")

leg.AddEntry(graphMed,"Median Expected (MET)","L")
leg.AddEntry(graph68,"#pm 1#sigma Expected","F")
leg.AddEntry(graph95,"#pm 2#sigma Expected","F")
MG = ROOT.TMultiGraph()

sm = 1
#EXPECTED
iexp=0
for i,mass in zip(range(len(EXPmasses)),EXPmasses):
  


  for j,mm in enumerate(allMasses): 
	if mm==mass: sm = 1.0/(xSec[j])
  
  if not divxsec: sm=1.0
  
  tree = EXPfile.Get("limit")
  

  
  if usegrid:    
    median = getOBSERVED(EXPfile,mass,0.5)
    dn68 = getOBSERVED(EXPfile,mass,0.16)
    up68 = getOBSERVED(EXPfile,mass,0.84)
    dn95 = getOBSERVED(EXPfile,mass,0.025)
    up95 = getOBSERVED(EXPfile,mass,0.975)
    #print median[0]
    #print up68[0]
    #print dn68[0]
    #print up95[0]
    #print dn95[0]
  else:
    median = array.array('d',[0])
    up68   = array.array('d',[0])
    dn68   = array.array('d',[0])
    up95   = array.array('d',[0])
    dn95   = array.array('d',[0])
    medianCalc("r_mH"+str(mass),tree,median,up68,dn68,up95,dn95,mass)


  print "mass = %d, median = %f, dn68 = %f, up68 = %f, dn95 = %f, up95 = %f" % (mass,median[0],dn68[0],up68[0],dn95[0],up95[0])
  print "medianerror = %f" % (median[1])

    
    
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
    #graphOne.SetPointError(iexp,0,0,0,0)

    graphmede.SetPoint(iexp,float(mass),median[0])
    graph68up.SetPoint(iexp,float(mass),up68[0])
    graph68dn.SetPoint(iexp,float(mass),dn68[0])
    graph95up.SetPoint(iexp,float(mass),up95[0])
    graph95dn.SetPoint(iexp,float(mass),dn95[0])

    if (len(median)>1):
      graphmede.SetPointError(iexp,0.0,median[1])
      graph68up.SetPointError(iexp,0.0,up68[1])
      graph68dn.SetPointError(iexp,0.0,dn68[1])
      graph95up.SetPointError(iexp,0.0,up95[1])
      graph95dn.SetPointError(iexp,0.0,dn95[1])

    iexp+=1


    print "%d: Median = %f" % (mass,median[0]*sm)


fitstring = "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x"

medfunc = ROOT.TF1("medfunc",fitstring,109.75,150.25);
up68func = ROOT.TF1("up68func",fitstring,109.75,150.25);
dn68func = ROOT.TF1("dn68func",fitstring,109.75,150.25);
up95func = ROOT.TF1("up95func",fitstring,109.75,150.25);
dn95func = ROOT.TF1("dn95func",fitstring,109.75,150.25);

graphmede.Fit(medfunc,"R,0,EX0","")
graph68up.Fit(up68func,"R,0,EX0","")
graph68dn.Fit(dn68func,"R,0,EX0","")
graph95up.Fit(up95func,"R,0,EX0","")
graph95dn.Fit(dn95func,"R,0,EX0","")


#for i,mass in zip(range(len(EXPmasses)),EXPmasses):
  #for j,mm in enumerate(allMasses): 
    #if mm==mass: sm = 1.0/(xSec[j])
  
  #if not doratio: sm=1.0
  
  #mediansmooth = medfunc.Eval(mass)
  
  #graphMed.SetPoint(i,mass,mediansmooth*sm)
  #graph68.SetPoint(i,mass,mediansmooth*sm)
  #graph95.SetPoint(i,mass,mediansmooth*sm)

  #diff95_up = abs(mediansmooth - up95func.Eval(mass))*sm
  #diff95_dn = abs(mediansmooth - dn95func.Eval(mass))*sm
  #diff68_up = abs(mediansmooth - up68func.Eval(mass))*sm 
  #diff68_dn = abs(mediansmooth - dn68func.Eval(mass))*sm 

  #graph68.SetPointError(i,0,0,diff68_dn,diff68_up)
  #graph95.SetPointError(i,0,0,diff95_dn,diff95_up)

 ## print "i = %d, sm = %f" % (i,sm)
  #print "med= %f, 95up = %f, 95dn = %f, 68up = %f, 68dn = %f" % (mediansmooth,up95func.Eval(mass),dn95func.Eval(mass),up68func.Eval(mass),dn68func.Eval(mass))



#OBSERVED
iobs=0
ibobs=0
for i,mass in zip(range(len(OBSmasses)),OBSmasses):

  for j,mm in enumerate(allMasses): 
	if mm==mass: sm = 1.0/(xSec[j])

  bsm = sm

  if not divxsec: sm=1.0

  #if obs[i][0]>0.0:
    #graphObs.SetPoint(iobs,float(mass),obs[i][0]*sm)
    #graphObs.SetPointError(iobs,0,0,0,0)
    ##if (mass==115 or mass==147): graphObs.SetPoint(iobs,float(mass),bobs[i][0]*sm)
    #print "%f: Observed CLs = %f" % (mass,obs[i][0]*sm)
    #iobs+=1
    
  ##if (bobs[i][0]>0.0 and bobs[i][0]>obs[i][0]):
  if bobs[i][0]>0.0:
    graphBObs.SetPoint(ibobs,float(mass),bobs[i][0]*bsm)
    graphBObs.SetPointError(ibobs,0,0,0,0)
    
    graphbmede.SetPoint(ibobs,float(mass),bobsmedian[i][0]*sm)
    graphb68dn.SetPoint(ibobs,float(mass),bobs68dn[i][0]*sm)
    graphb68up.SetPoint(ibobs,float(mass),bobs68up[i][0]*sm)
    graphb95dn.SetPoint(ibobs,float(mass),bobs95dn[i][0]*sm)
    graphb95up.SetPoint(ibobs,float(mass),bobs95up[i][0]*sm)
    
    print "%f: Observed Bay = %f" % (mass,bobs[i][0]*sm)
    ibobs+=1    
    
  graphObs.SetPoint(iobs,float(mass),obs[i][0]*sm)
  graphObs.SetPointError(iobs,0,0,0,0)
  #if (mass==115 or mass==147): graphObs.SetPoint(iobs,float(mass),bobs[i][0]*sm)
  print "%f: Observed CLs = %f" % (mass,obs[i][0]*sm)
  iobs+=1
  

  #graphBObs.SetPoint(ibobs,float(mass),bobs[i][0]*sm)
  #graphBObs.SetPointError(ibobs,0,0,0,0)
  
  #graphbmede.SetPoint(ibobs,float(mass),bobsmedian[i][0]*sm)
  #graphb68dn.SetPoint(ibobs,float(mass),bobs68dn[i][0]*sm)
  #graphb68up.SetPoint(ibobs,float(mass),bobs68up[i][0]*sm)
  #graphb95dn.SetPoint(ibobs,float(mass),bobs95dn[i][0]*sm)
  #graphb95up.SetPoint(ibobs,float(mass),bobs95up[i][0]*sm)
  
  #print "%f: Observed Bay = %f" % (mass,bobs[i][0]*sm)
  #ibobs+=1        



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
#graph95.SetFillColor(ROOT.kYellow-4)
#graph95.SetFillColor(ROOT.kYellow-9)
#graph95.SetFillStyle(3001)
#graph95.SetFillColor(ROOT.kGreen-10)
#graph95.SetFillStyle(3001)
#graph68.SetFillColor(ROOT.kGreen+2)
#graph68.SetFillColor(ROOT.kGreen)
#graph68.SetFillColor(ROOT.kGreen-10)
#graph68.SetFillStyle(3001)
#graph68.SetFillColor(ROOT.kYellow-9)
#graph68.SetFillStyle(3001)

#graph95.SetFillColor(ROOT.kGreen)
#graph68.SetFillColor(ROOT.kYellow-4)

graph68.SetFillColor(ROOT.kYellow-4)
graph95.SetFillColor(ROOT.kGreen)


graphMed.SetLineStyle(2)
graphMed.SetLineColor(2)
graphMed.SetLineWidth(3)
graphObs.SetLineWidth(3)
graphBObs.SetLineWidth(4)
graphCObs.SetLineWidth(3)
graphDObs.SetLineWidth(3)


#graphObs.SetMarkerStyle(2)

graphOne.SetLineWidth(4)

xsecbandfill=3244
#xsecbandfill=0

#myGraphXSecSM.SetLineColor(ROOT.kAzure+7)

myGraphXSecSM.SetLineStyle(1)
myGraphXSecSM.SetLineColor(ROOT.kAzure+7)
myGraphXSecSM.SetLineWidth(2)
myGraphXSecSM.SetFillColor(ROOT.kAzure+7)
myGraphXSecSM.SetFillStyle(xsecbandfill)
myGraphXSec2SM.SetLineStyle(1)
myGraphXSec2SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec2SM.SetLineWidth(2)
myGraphXSec2SM.SetFillColor(ROOT.kAzure+7)
myGraphXSec2SM.SetFillStyle(xsecbandfill)
myGraphXSec5SM.SetLineStyle(1)
myGraphXSec5SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec5SM.SetLineWidth(2)
myGraphXSec5SM.SetFillColor(ROOT.kAzure+7)
myGraphXSec5SM.SetFillStyle(xsecbandfill)
myGraphXSec10SM.SetLineStyle(1)
myGraphXSec10SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec10SM.SetLineWidth(2)
myGraphXSec10SM.SetFillColor(ROOT.kAzure+7)
myGraphXSec10SM.SetFillStyle(xsecbandfill)
myGraphXSec20SM.SetLineStyle(1)
myGraphXSec20SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec20SM.SetLineWidth(2)
myGraphXSec20SM.SetFillColor(ROOT.kAzure+7)
myGraphXSec20SM.SetFillStyle(xsecbandfill)
myGraphXSec30SM.SetLineStyle(1)
myGraphXSec30SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec30SM.SetLineWidth(2)
myGraphXSec30SM.SetFillColor(ROOT.kAzure+7)
myGraphXSec30SM.SetFillStyle(xsecbandfill)
myGraphXSec40SM.SetLineStyle(1)
myGraphXSec40SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec40SM.SetLineWidth(2)
myGraphXSec40SM.SetFillColor(ROOT.kAzure+7)
myGraphXSec40SM.SetFillStyle(xsecbandfill)



# use 95 as the default guy

MG.Add(graph95)
MG.Add(graph68)

#MG.Add(graphCObs)
#MG.Add(graphDObs)


#MG.Add(graphOne)
#MG.Add(graphmede)
MG.Add(graphMed)
MG.Add(myGraphXSecSM)
if not doratio:
  #if dofp: MG.Add(myGraphXSec2SM)
  #MG.Add(myGraphXSec20SM)
  #MG.Add(myGraphXSec30SM)
  if not dofp: MG.Add(myGraphXSec5SM)
  #MG.Add(myGraphXSec10SM)
  
#MG.Add(myGraphXSec40SM)
#graphObs.SetMarkerStyle(22)
#graphObs.SetMarkerSize(1.5)

graphBObs.SetLineStyle(7)
#graphBObs.SetLineColor(4)
graphBObs.SetLineColor(ROOT.kBlue-4)

graphCObs.SetLineStyle(7)
#graphBObs.SetLineColor(4)
graphCObs.SetLineColor(ROOT.kBlue-4)

graphDObs.SetLineStyle(7)
#graphBObs.SetLineColor(4)
#graphDObs.SetLineColor(ROOT.kSpring-7)
graphDObs.SetLineColor(ROOT.kBlue-4)

#graphBObs.SetMarkerStyle(20)
#graphBObs.SetMarkerSize(1.5)


MG.Add(graphObs)



#MG.Add(graphBObs)

#MG.Add(graphbmede)
#MG.Add(graphb68up)
#MG.Add(graphb68dn)
#MG.Add(graphb95up)
#MG.Add(graphb95dn)

# -------------------------------------
C = ROOT.TCanvas("#int L = %s"%intlumi,"#int L = %s"%intlumi,1600,1100)

C.SetGrid(True)
#MG.Draw("AL3")

#if not dofp:
if not (dofp and doratio):
  #dummy = ROOT.TH1D("dummy","",1,109.75,154.5);
  dummy = ROOT.TH1D("dummy","",1,109.75,150.);
else:
  dummy = ROOT.TH1D("dummy","",1,109.75,150.);  
#else:
#  dummy = ROOT.TH1D("dummy","",1,106.5,140.);
dummy.SetStats(False)

dummy.Draw("AXIS")

if (dofp and doratio):
  MG.Draw("L3")  
  dummy.Draw("AXIGSAME")
  #MG.Draw("L3")
else:
  MG.Draw("L3")
  dummy.Draw("AXIGSAME")

#MG.Draw("C4")

#MG.GetXaxis().SetRangeUser(110,145)
dummy.GetXaxis().SetTitle("m_{H} (GeV/c^{2})")
#MG.GetXaxis().SetRangeUser(110,140)
maxy=0.23
if dofp:
  maxy=0.20


#dummy.SetMinimum(0.0);
dummy.SetMaximum(maxy);
if doratio:
#if True:
  #MG.GetYaxis().SetRangeUser(0.0,20.0)
  #MG.GetYaxis().SetMaximum(20.0)
  #maxy=4.0
  maxy=40.0
  dummy.SetMaximum(maxy)
  #maxy=35.0
  

if dofp and doratio:
  maxy=20
  dummy.SetMinimum(0.1)
  dummy.SetMaximum(30)
  C.SetLogy()
else:
  dummy.SetMinimum(0)

#MG.GetYaxis().SetRangeUser(0.0,maxy)  
#MG.GetYaxis().SetMaximum(maxy)  

#MG.GetYaxis().SetTitle("\sigma #times BR(H#rightarrow #gamma #gamma)_{95%CL} (pb)")
dummy.GetYaxis().SetTitle("#sigma#timesBR(H#rightarrow#gamma#gamma)_{95%CL}(pb)")
if doratio: 
  dummy.GetYaxis().SetTitle("\sigma(H#rightarrow #gamma #gamma)_{95%CL}/\sigma(H#rightarrow #gamma #gamma)_{SM}")
  if dofp:
    dummy.GetYaxis().SetTitle("\sigma(H#rightarrow #gamma #gamma)_{95%CL}/\sigma(H#rightarrow #gamma #gamma)_{FP}")
MG.SetTitle("#int L = %s"%intlumi)
mytext = ROOT.TLatex()
mytext.SetTextSize(0.04)

labely=1.0
if not doratio:
  if dofp: labely=labely*0.00387
  else: labely=labely*0.0163

if dofp:
  pass;
  if doratio:
    mytext.DrawLatex(150.5,0.8*labely,"1#times#sigma_{FP}")
  if not doratio:
    mytext.DrawLatex(150.5,1.0*labely,"1#times#sigma_{FP}")

else:
  mytext.DrawLatex(150.5,0.8*labely,"1#times#sigma_{SM}")
  if not doratio:
    mytext.DrawLatex(150.5,4.8*labely,"5#times#sigma_{SM}")

leg.Draw()
mytext.DrawLatex(135,0.90*maxy,"#splitline{CMS preliminary}{#sqrt{s} = 7 TeV L = 4.76 fb^{-1}}")

C.SaveAs("LimitSMRel.pdf")


graphObs.Print()
graphMed.Print()

ox = ROOT.Double(0.0)
oy = ROOT.Double(0.0)
oex = ROOT.Double(0.0)
oey = ROOT.Double(0.0)
print "mh, observed, expected, plusone, minusone, plustwo, minustwo"
for i in range(graphObs.GetN()):
 graphObs.GetPoint(i,ox,oy)
 graphMed.GetPoint(i,oex,oey)
 print "%f, %f, %f, %f, %f, %f, %f" % (ox, oy, oey, oey+graph68.GetErrorYhigh(i),oey-graph68.GetErrorYlow(i),oey+graph95.GetErrorYhigh(i),oey-graph95.GetErrorYlow(i))


