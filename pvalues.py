import os

from multiprocessing import Pool


mmin=110.0
mmax=150.0
nm=81

card="/home/bendavid/cms/cmssw/021/CMSSW_4_2_3_patch2/src/MitHgg/macros/cards/feb16/hggcardfeb16smmvavbf.txt"

cmdlist=[]

for i in range(nm):
  mass = mmin + i*(mmax-mmin)/float(nm-1)
  obsexec = "combine -d %s -m %g -M ProfileLikelihood --significance --pvalue -n SigFromGridObs" % (card,mass)
  print obsexec
  cmdlist.append(obsexec) 
 
pool = Pool(processes=20)
pool.map_async(os.system,cmdlist)
pool.close()
pool.join()    
 
os.system("hadd smpvalues.root higgsCombineSigFromGridObs*.root")
os.system("rm higgsCombineSigFromGrid*.root")
