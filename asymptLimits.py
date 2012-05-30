import os

from multiprocessing import Pool


mmin=110.0
mmax=150.0
nm=81

card="testcard.txt"



cmdlist=[]

for i in range(nm):
  mass = mmin + i*(mmax-mmin)/float(nm-1)
  obsexec = "combine -d %s -m %g -U -M Asymptotic --rRelAcc=0.001 --rAbsAcc=0.001 --minimizerStrategy=0 --rMax=30 --run=expected -n LimitsFromGridObs" % (card,mass)
  print obsexec
  cmdlist.append(obsexec)

 
 
pool = Pool(processes=20)
pool.map_async(os.system,cmdlist)
pool.close()
pool.join()    
 
os.system("hadd smrel.root higgsCombineLimitsFromGridObs*.root")
os.system("rm higgsCombineLimitsFromGrid*.root")
