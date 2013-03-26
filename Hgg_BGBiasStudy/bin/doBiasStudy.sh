#! /bin/bash

root -l -q -b "macros/doBiasStudy.C+(\"hgg_8TeV_2013moriond_bdt0\",\"BDT Class 0\")"
root -l -q -b "macros/doBiasStudy.C+(\"hgg_8TeV_2013moriond_bdt1\",\"BDT Class 1\")"
root -l -q -b "macros/doBiasStudy.C+(\"hgg_8TeV_2013moriond_bdt2\",\"BDT Class 2\")"
root -l -q -b "macros/doBiasStudy.C+(\"hgg_8TeV_2013moriond_bdt3\",\"BDT Class 3\")"
root -l -q -b "macros/doBiasStudy.C+( \"hgg_8TeV_2013moriond_dijettight\"  ,\"Dijet Tight\")"
root -l -q -b "macros/doBiasStudy.C+( \"hgg_8TeV_2013moriond_dijetloose\"  ,\"Dijet Loose\")"
root -l -q -b "macros/doBiasStudy.C+( \"hgg_8TeV_2013moriond_mtag\"	     ,\"Muon Tag\")"
root -l -q -b "macros/doBiasStudy.C+( \"hgg_8TeV_2013moriond_etag\"	     ,\"Electron Tag\")"
root -l -q -b "macros/doBiasStudy.C+( \"hgg_8TeV_2013moriond_mettag\"      ,\"MET Tag\")"
