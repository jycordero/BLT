#! /bin/sh


# Don't try to do anything if the wrong version of ROOT is loaded
myROOTSYS="/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.10.02-19565/x86_64-slc6-gcc62-opt"
if [ "$ROOTSYS" != "$myROOTSYS" ]
then
    echo "ROOT 6.10 not loaded"
else

macroPath="/uscms/home/jrainbol/nobackup/CMSSW_8_0_26_patch1/src/BLT/BLTAnalysis/macros"
eosPath="/store/user/jrainbol/Trees/2016_12d"


##  MUON  ##

#   #  Run B

#   if eval "hadd muon_2016_1.root output_muon_2016B_v2_{1..25}.root"
#   then
#       if eval "rm output_muon_2016B_v2_{1..25}.root"
#       then
#           echo "Deleted output_muon_2016B_v2_{1..25}.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"muon_2016_1.root\", \"tree_muon_2016B\", \"tree_muon_2016\")"
#       rootrm  "muon_2016_1.root:tree_muon_2016B"
#       rootmv  "muon_2016_1.root:TotalEvents_muon_2016B" "muon_2016_1.root:TotalEvents_muon_2016"
#       rootmv  "muon_2016_1.root:AcceptedEvents_muon_2016B" "muon_2016_1.root:AcceptedEvents_muon_2016"
#   fi
#   if eval "hadd muon_2016_2.root output_muon_2016B_v2_{26..51}.root"
#   then
#       if eval "rm output_muon_2016B_v2_{26..51}.root"
#       then
#           echo "Deleted output_muon_2016B_v2_{26..51}.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"muon_2016_2.root\", \"tree_muon_2016B\", \"tree_muon_2016\")"
#       rootrm  "muon_2016_2.root:tree_muon_2016B"
#       rootmv  "muon_2016_2.root:TotalEvents_muon_2016B" "muon_2016_2.root:TotalEvents_muon_2016"
#       rootmv  "muon_2016_2.root:AcceptedEvents_muon_2016B" "muon_2016_2.root:AcceptedEvents_muon_2016"
#   fi


#   # Run C

#   if eval "hadd muon_2016_3.root output_muon_2016C_v1_*.root"
#   then
#       if eval "rm output_muon_2016C_v1_*.root"
#       then
#           echo "Deleted output_muon_2016C_v1_*.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"muon_2016_3.root\", \"tree_muon_2016C\", \"tree_muon_2016\")"
#       rootrm  "muon_2016_3.root:tree_muon_2016C"
#       rootmv  "muon_2016_3.root:TotalEvents_muon_2016C" "muon_2016_3.root:TotalEvents_muon_2016"
#       rootmv  "muon_2016_3.root:AcceptedEvents_muon_2016C" "muon_2016_3.root:AcceptedEvents_muon_2016"
#   fi


#   # Run D

#   if eval "hadd muon_2016_4.root output_muon_2016D_v1_*.root"
#   then
#       if eval "rm output_muon_2016D_v1_*.root"
#       then
#           echo "Deleted output_muon_2016D_v1_*.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"muon_2016_4.root\", \"tree_muon_2016D\", \"tree_muon_2016\")"
#       rootrm  "muon_2016_4.root:tree_muon_2016D"
#       rootmv  "muon_2016_4.root:TotalEvents_muon_2016D" "muon_2016_4.root:TotalEvents_muon_2016"
#       rootmv  "muon_2016_4.root:AcceptedEvents_muon_2016D" "muon_2016_4.root:AcceptedEvents_muon_2016"
#   fi


#   # Run E

#   if eval "hadd muon_2016_5.root output_muon_2016E_v1_*.root"
#   then
#       if eval "rm output_muon_2016E_v1_*.root"
#       then
#           echo "Deleted output_muon_2016E_v1_*.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"muon_2016_5.root\", \"tree_muon_2016E\", \"tree_muon_2016\")"
#       rootrm  "muon_2016_5.root:tree_muon_2016E"
#       rootmv  "muon_2016_5.root:TotalEvents_muon_2016E" "muon_2016_5.root:TotalEvents_muon_2016"
#       rootmv  "muon_2016_5.root:AcceptedEvents_muon_2016E" "muon_2016_5.root:AcceptedEvents_muon_2016"
#   fi


#   # Run F

#   if eval "hadd muon_2016_6.root output_muon_2016F_v1_*.root"
#   then
#       if eval "rm output_muon_2016F_v1_*.root"
#       then
#           echo "Deleted output_muon_2016F_v1_*.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"muon_2016_6.root\", \"tree_muon_2016F\", \"tree_muon_2016\")"
#       rootrm  "muon_2016_6.root:tree_muon_2016F"
#       rootmv  "muon_2016_6.root:TotalEvents_muon_2016F" "muon_2016_6.root:TotalEvents_muon_2016"
#       rootmv  "muon_2016_6.root:AcceptedEvents_muon_2016F" "muon_2016_6.root:AcceptedEvents_muon_2016"
#   fi


#   # Run G

#   if eval "hadd muon_2016_7.root output_muon_2016G_v1_{1..25}.root"
#   then
#       if eval "rm output_muon_2016G_v1_{1..25}.root"
#       then
#           echo "Deleted output_muon_2016G_v1_{1..25}.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"muon_2016_7.root\", \"tree_muon_2016G\", \"tree_muon_2016\")"
#       rootrm  "muon_2016_7.root:tree_muon_2016G"
#       rootmv  "muon_2016_7.root:TotalEvents_muon_2016G" "muon_2016_7.root:TotalEvents_muon_2016"
#       rootmv  "muon_2016_7.root:AcceptedEvents_muon_2016G" "muon_2016_7.root:AcceptedEvents_muon_2016"
#   fi
#   if eval "hadd muon_2016_8.root output_muon_2016G_v1_{26..51}.root"
#   then
#       if eval "rm output_muon_2016G_v1_{26..51}.root"
#       then
#           echo "Deleted output_muon_2016G_v1_{26..51}.root"
#       fi
#       root.exe -q -b "../../renameTree.C(\"muon_2016_8.root\", \"tree_muon_2016G\", \"tree_muon_2016\")"
#       rootrm  "muon_2016_8.root:tree_muon_2016G"
#       rootmv  "muon_2016_8.root:TotalEvents_muon_2016G" "muon_2016_8.root:TotalEvents_muon_2016"
#       rootmv  "muon_2016_8.root:AcceptedEvents_muon_2016G" "muon_2016_8.root:AcceptedEvents_muon_2016"
#   fi


#   # Run H

#   if eval "hadd muon_2016_9.root output_muon_2016H_v2_{1..26}.root"
#   then
#       if eval "rm output_muon_2016H_v1_{1..26}.root"
#       then
#           echo "Deleted output_muon_2016H_v1_{1..26}.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"muon_2016_9.root\", \"tree_muon_2016H\", \"tree_muon_2016\")"
#       rootrm  "muon_2016_9.root:tree_muon_2016H"
#       rootmv  "muon_2016_9.root:TotalEvents_muon_2016H" "muon_2016_9.root:TotalEvents_muon_2016"
#       rootmv  "muon_2016_9.root:AcceptedEvents_muon_2016H" "muon_2016_9.root:AcceptedEvents_muon_2016"
#   fi
#   if eval "hadd muon_2016_10.root output_muon_2016H_v2_{27..51}.root output_muon_2016H_v3_*.root"
#   then
#       if eval "rm output_muon_2016H_v2_{27..51}.root"
#       then
#           echo "Deleted output_muon_2016H_v2_{27..51}.root"
#       fi
#       if eval "rm output_muon_2016H_v3_*.root"
#       then
#           echo "Deleted output_muon_2016H_v3_*.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"muon_2016_10.root\", \"tree_muon_2016H\", \"tree_muon_2016\")"
#       rootrm  "muon_2016_10.root:tree_muon_2016H"
#       rootmv  "muon_2016_10.root:TotalEvents_muon_2016H" "muon_2016_10.root:TotalEvents_muon_2016"
#       rootmv  "muon_2016_10.root:AcceptedEvents_muon_2016H" "muon_2016_10.root:AcceptedEvents_muon_2016"
#   fi

#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/muon_2016
#   xrdcp -f muon_2016_*.root root://cmseos.fnal.gov/"$eosPath"/muon_2016/.



##  ELECTRON  ##

#   # Run B

#   if eval "hadd electron_2016_1.root output_electron_2016B_v2_{1..25}.root"
#   then
#       if eval "rm output_electron_2016B_v2_{1..25}.root"
#       then
#           echo "Deleted output_electron_2016B_v2_{1..25}.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"electron_2016_1.root\", \"tree_electron_2016B\", \"tree_electron_2016\")"
#       rootrm  "electron_2016_1.root:tree_electron_2016B"
#       rootmv  "electron_2016_1.root:TotalEvents_electron_2016B" "electron_2016_1.root:TotalEvents_electron_2016"
#       rootmv  "electron_2016_1.root:AcceptedEvents_electron_2016B" "electron_2016_1.root:AcceptedEvents_electron_2016"
#   fi
#   if eval "hadd electron_2016_2.root output_electron_2016B_v2_{26..51}.root"
#   then
#       if eval "rm output_electron_2016B_v2_{26..51}.root"
#       then
#           echo "Deleted output_electron_2016B_v2_{26..51}.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"electron_2016_2.root\", \"tree_electron_2016B\", \"tree_electron_2016\")"
#       rootrm  "electron_2016_2.root:tree_electron_2016B"
#       rootmv  "electron_2016_2.root:TotalEvents_electron_2016B" "electron_2016_2.root:TotalEvents_electron_2016"
#       rootmv  "electron_2016_2.root:AcceptedEvents_electron_2016B" "electron_2016_2.root:AcceptedEvents_electron_2016"
#   fi


#   # Run C

#   if eval "hadd electron_2016_3.root output_electron_2016C_v1_*.root"
#   then
#       if eval "rm output_electron_2016C_v1_*.root"
#       then
#           echo "Deleted output_electron_2016C_v1_*.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"electron_2016_3.root\", \"tree_electron_2016C\", \"tree_electron_2016\")"
#       rootrm  "electron_2016_3.root:tree_electron_2016C"
#       rootmv  "electron_2016_3.root:TotalEvents_electron_2016C" "electron_2016_3.root:TotalEvents_electron_2016"
#       rootmv  "electron_2016_3.root:AcceptedEvents_electron_2016C" "electron_2016_3.root:AcceptedEvents_electron_2016"
#   fi


#   # Run D

#   if eval "hadd electron_2016_4.root output_electron_2016D_v1_*.root"
#   then
#       if eval "rm output_electron_2016D_v1_*.root"
#       then
#           echo "Deleted output_electron_2016D_v1_*.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"electron_2016_4.root\", \"tree_electron_2016D\", \"tree_electron_2016\")"
#       rootrm  "electron_2016_4.root:tree_electron_2016D"
#       rootmv  "electron_2016_4.root:TotalEvents_electron_2016D" "electron_2016_4.root:TotalEvents_electron_2016"
#       rootmv  "electron_2016_4.root:AcceptedEvents_electron_2016D" "electron_2016_4.root:AcceptedEvents_electron_2016"
#   fi


#   # Run E

#   if eval "hadd electron_2016_5.root output_electron_2016E_v1_*.root"
#   then
#       if eval "rm output_electron_2016E_v1_*.root"
#       then
#           echo "Deleted output_electron_2016E_v1_*.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"electron_2016_5.root\", \"tree_electron_2016E\", \"tree_electron_2016\")"
#       rootrm  "electron_2016_5.root:tree_electron_2016E"
#       rootmv  "electron_2016_5.root:TotalEvents_electron_2016E" "electron_2016_5.root:TotalEvents_electron_2016"
#       rootmv  "electron_2016_5.root:AcceptedEvents_electron_2016E" "electron_2016_5.root:AcceptedEvents_electron_2016"
#   fi


#   # Run F

#   if eval "hadd electron_2016_6.root output_electron_2016F_v1_*.root"
#   then
#       if eval "rm output_electron_2016F_v1_*.root"
#       then
#           echo "Deleted output_electron_2016F_v1_*.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"electron_2016_6.root\", \"tree_electron_2016F\", \"tree_electron_2016\")"
#       rootrm  "electron_2016_6.root:tree_electron_2016F"
#       rootmv  "electron_2016_6.root:TotalEvents_electron_2016F" "electron_2016_6.root:TotalEvents_electron_2016"
#       rootmv  "electron_2016_6.root:AcceptedEvents_electron_2016F" "electron_2016_6.root:AcceptedEvents_electron_2016"
#   fi


#   # Run G

#   if eval "hadd electron_2016_7.root output_electron_2016G_v1_{1..25}.root"
#   then
#       if eval "rm output_electron_2016G_v1_{1..25}.root"
#       then
#           echo "Deleted output_electron_2016G_v1_{1..25}.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"electron_2016_7.root\", \"tree_electron_2016G\", \"tree_electron_2016\")"
#       rootrm  "electron_2016_7.root:tree_electron_2016G"
#       rootmv  "electron_2016_7.root:TotalEvents_electron_2016G" "electron_2016_7.root:TotalEvents_electron_2016"
#       rootmv  "electron_2016_7.root:AcceptedEvents_electron_2016G" "electron_2016_7.root:AcceptedEvents_electron_2016"
#   fi
#   if eval "hadd electron_2016_8.root output_electron_2016G_v1_{26..51}.root"
#   then
#       if eval "rm output_electron_2016G_v1_{26..51}.root"
#       then
#           echo "Deleted output_electron_2016G_v1_{26..51}.root"
#       fi
#       root.exe -q -b "../../renameTree.C(\"electron_2016_8.root\", \"tree_electron_2016G\", \"tree_electron_2016\")"
#       rootrm  "electron_2016_8.root:tree_electron_2016G"
#       rootmv  "electron_2016_8.root:TotalEvents_electron_2016G" "electron_2016_8.root:TotalEvents_electron_2016"
#       rootmv  "electron_2016_8.root:AcceptedEvents_electron_2016G" "electron_2016_8.root:AcceptedEvents_electron_2016"
#   fi


#   # Run H

#   if eval "hadd electron_2016_9.root output_electron_2016H_v2_{1..26}.root"
#   then
#       if eval "rm output_electron_2016H_v1_{1..26}.root"
#       then
#           echo "Deleted output_electron_2016H_v1_{1..26}.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"electron_2016_9.root\", \"tree_electron_2016H\", \"tree_electron_2016\")"
#       rootrm  "electron_2016_9.root:tree_electron_2016H"
#       rootmv  "electron_2016_9.root:TotalEvents_electron_2016H" "electron_2016_9.root:TotalEvents_electron_2016"
#       rootmv  "electron_2016_9.root:AcceptedEvents_electron_2016H" "electron_2016_9.root:AcceptedEvents_electron_2016"
#   fi
#   if eval "hadd electron_2016_10.root output_electron_2016H_v2_{27..51}.root output_electron_2016H_v3_*.root"
#   then
#       if eval "rm output_electron_2016H_v2_{27..51}.root"
#       then
#           echo "Deleted output_electron_2016H_v2_{27..51}.root"
#       fi
#       if eval "rm output_electron_2016H_v3_*.root"
#       then
#           echo "Deleted output_electron_2016H_v3_*.root"
#       fi
#       root.exe -q -b "$macroPath""/renameTree.C(\"electron_2016_10.root\", \"tree_electron_2016H\", \"tree_electron_2016\")"
#       rootrm  "electron_2016_10.root:tree_electron_2016H"
#       rootmv  "electron_2016_10.root:TotalEvents_electron_2016H" "electron_2016_10.root:TotalEvents_electron_2016"
#       rootmv  "electron_2016_10.root:AcceptedEvents_electron_2016H" "electron_2016_10.root:AcceptedEvents_electron_2016"
#   fi

#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/electron_2016
#   xrdcp -f electron_2016_*.root root://cmseos.fnal.gov/"$eosPath"/electron_2016/.




##  DRELL-YAN  ##

#   # Inclusive

#   if eval "hadd zjets_10to50_1.root output_DYJetsToLL_M-10to50_*.root"
#   then
#       if eval "rm output_DYJetsToLL_M-10to50_*.root"
#       then
#           echo "Deleted output_DYJetsToLL_M-10to50_*.root"
#       fi
#   fi
#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/zjets_10to50
#   xrdcp -f zjets_10to50_*.root root://cmseos.fnal.gov/"$eosPath"/zjets_10to50/.

#   if eval "hadd zjets_50_1.root output_DYJetsToLL_M-50_{1..2}.root"
#   then
#       if eval "rm output_DYJetsToLL_M-50_{1..2}.root"
#       then
#           echo "Deleted output_DYJetsToLL_M-50_{1..2}.root"
#       fi
#   fi
#   if eval "hadd zjets_50_2.root output_DYJetsToLL_M-50_{3..4}.root"
#   then
#       if eval "rm output_DYJetsToLL_M-50_{3..4}.root"
#       then
#           echo "Deleted output_DYJetsToLL_M-50_{3..4}.root"
#       fi
#   fi
#   if eval "hadd zjets_50_3.root output_DYJetsToLL_M-50_{5..6}.root"
#   then
#       if eval "rm output_DYJetsToLL_M-50_{5..6}.root"
#       then
#           echo "Deleted output_DYJetsToLL_M-50_{5..6}.root"
#       fi
#   fi
#   if eval "hadd zjets_50_4.root output_DYJetsToLL_M-50_{7..8}.root"
#   then
#       if eval "rm output_DYJetsToLL_M-50_{7..8}.root"
#       then
#           echo "Deleted output_DYJetsToLL_M-50_{7..8}.root"
#       fi
#   fi
#   if eval "hadd zjets_50_5.root output_DYJetsToLL_M-50_{9..10}.root"
#   then
#       if eval "rm output_DYJetsToLL_M-50_{9..10}.root"
#       then
#           echo "Deleted output_DYJetsToLL_M-50_{9..10}.root"
#       fi
#   fi
#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/zjets_50
#   xrdcp -f zjets_50_*.root root://cmseos.fnal.gov/"$eosPath"/zjets_50/.


#   # 1 parton

#   if eval "hadd z1jets_10to50_1.root output_DY1JetsToLL_M-10to50_*.root"
#   then
#       if eval "rm output_DY1JetsToLL_M-10to50_*.root"
#       then
#           echo "Deleted output_DY1JetsToLL_M-10to50_*.root"
#       fi
#   fi
#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/z1jets_10to50
#   xrdcp -f z1jets_10to50_*.root root://cmseos.fnal.gov/"$eosPath"/z1jets_10to50/.

#   for i in {1..6}
#   do
#       if mv "output_DY1JetsToLL_M-50_""$i"".root" "z1jets_50_""$i"".root"
#       then
#           echo "Moved output_DY1JetsToLL_M-50_""$i"".root to z1jets_50_""$i"".root"
#       fi
#   done
#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/z1jets_50
#   xrdcp -f z1jets_50_*.root root://cmseos.fnal.gov/"$eosPath"/z1jets_50/.


#   # 2 parton

#   if eval "hadd z2jets_10to50_1.root output_DY2JetsToLL_M-10to50_*.root"
#   then
#       if eval "rm output_DY2JetsToLL_M-10to50_*.root"
#       then
#           echo "Deleted output_DY2JetsToLL_M-10to50_*.root"
#       fi
#   fi
#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/z2jets_10to50
#   xrdcp -f z2jets_10to50_*.root root://cmseos.fnal.gov/"$eosPath"/z2jets_10to50/.

#   if eval "hadd z2jets_50_1.root output_DY2JetsToLL_M-50_{1,4,5}.root"
#   then
#       if eval "rm output_DY2JetsToLL_M-50_{1,4,5}.root"
#       then
#           echo "Deleted output_DY2JetsToLL_M-50_{1,4,5}.root"
#       fi
#   fi
#   if eval "hadd z2jets_50_2.root output_DY2JetsToLL_M-50_{2,6,7}.root"
#   then
#       if eval "rm output_DY2JetsToLL_M-50_{2,6,7}.root"
#       then
#           echo "Deleted output_DY2JetsToLL_M-50_{2,6,7}.root"
#       fi
#   fi
#   if eval "hadd z2jets_50_3.root output_DY2JetsToLL_M-50_{3,8,9}.root"
#   then
#       if eval "rm output_DY2JetsToLL_M-50_{3,8,9}.root"
#       then
#           echo "Deleted output_DY2JetsToLL_M-50_{3,8,9}.root"
#       fi
#   fi
#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/z2jets_50
#   xrdcp -f z2jets_50_*.root root://cmseos.fnal.gov/"$eosPath"/z2jets_50/.


#   # 3 parton

#   if eval "hadd z3jets_10to50_1.root output_DY3JetsToLL_M-10to50_*.root"
#   then
#       if eval "rm output_DY3JetsToLL_M-10to50_*.root"
#       then
#           echo "Deleted output_DY3JetsToLL_M-10to50_*.root"
#       fi
#   fi
#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/z3jets_10to50
#   xrdcp -f z3jets_10to50_*.root root://cmseos.fnal.gov/"$eosPath"/z3jets_10to50/.

#   if eval "hadd z3jets_50_1.root output_DY3JetsToLL_M-50_*.root"
#   then
#       if eval "rm output_DY3JetsToLL_M-50_*.root"
#       then
#           echo "Deleted output_DY3JetsToLL_M-50_*.root"
#       fi
#   fi
#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/z3jets_50
#   xrdcp -f z3jets_50_*.root root://cmseos.fnal.gov/"$eosPath"/z3jets_50/.


#   # 4 parton

#   if eval "hadd z4jets_10to50_1.root output_DY4JetsToLL_M-10to50_*.root"
#   then
#       if eval "rm output_DY4JetsToLL_M-10to50_*.root"
#       then
#           echo "Deleted output_DY4JetsToLL_M-10to50_*.root"
#       fi
#   fi
#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/z4jets_10to50
#   xrdcp -f z4jets_10to50_*.root root://cmseos.fnal.gov/"$eosPath"/z4jets_10to50/.

#   if eval "hadd z4jets_50_1.root output_DY4JetsToLL_M-50_*.root"
#   then
#       if eval "rm output_DY4JetsToLL_M-50_*.root"
#       then
#           echo "Deleted output_DY4JetsToLL_M-50_*.root"
#       fi
#   fi
#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/z4jets_50
#   xrdcp -f z4jets_50_*.root root://cmseos.fnal.gov/"$eosPath"/z4jets_50/.



##  TTBAR  ##

#   if eval "hadd ttbar_1.root output_ttbar_inclusive_*.root"
#   then
#       if eval "rm output_ttbar_inclusive_*.root"
#       then
#           echo "Deleted output_ttbar_inclusive_*.root"
#       fi
#   fi
#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/ttbar
#   xrdcp -f ttbar_*.root root://cmseos.fnal.gov/"$eosPath"/ttbar/.
#   

 
##  W+ DIBOSON  ##

#   # WW to 2l 2nu
#   if eval "hadd ww_2l2nu_1.root output_WW_*.root"
#   then
#       if eval "rm output_WW_*.root"
#       then
#           echo "Deleted output_WW_*.root"
#       fi
        root.exe -q -b "$macroPath""/renameTree.C(\"ww_2l2nu_1.root\", \"tree_ww\", \"tree_ww_2l2nu\")"
        rootrm  "ww_2l2nu_1.root:tree_ww"
        rootmv  "ww_2l2nu_1.root:TotalEvents_ww" "ww_2l2nu_1.root:TotalEvents_ww_2l2nu"
        rootmv  "ww_2l2nu_1.root:AcceptedEvents_ww" "ww_2l2nu_1.root:AcceptedEvents_ww_2l2nu"
#   fi
    eos root://cmseos.fnal.gov mkdir -p "$eosPath"/ww_2l2nu
    xrdcp -f ww_2l2nu_*.root root://cmseos.fnal.gov/"$eosPath"/ww_2l2nu/.

#   # WZ to 2l 2q
#   if eval "hadd wz_2l2q_1.root output_WZJetsTo2L2Q_{2..8}.root"
#   then
#       if eval "rm output_WZJetsTo2L2Q_{2..8}.root"
#       then
#           echo "Deleted output_WZJetsTo2L2Q_{2..8}.root"
#       fi
#   fi
#   if eval "hadd wz_2l2q_2.root output_WZJetsTo2L2Q_{9..15}.root"
#   then
#       if eval "rm output_WZJetsTo2L2Q_{9..15}.root"
#       then
#           echo "Deleted output_WZJetsTo2L2Q_{9..15}.root"
#       fi
#   fi
#   if eval "hadd wz_2l2q_3.root output_WZJetsTo2L2Q_{16..20}.root"
#   then
#       if eval "rm output_WZJetsTo2L2Q_{16..20}.root"
#       then
#           echo "Deleted output_WZJetsTo2L2Q_{16..20}.root"
#       fi
#   fi
#   if eval "hadd wz_2l2q_4.root output_WZJetsTo2L2Q_1.root output_WZJetsTo2L2Q_{21..25}.root"
#   then
#       if eval "rm output_WZJetsTo2L2Q_1.root output_WZJetsTo2L2Q_{21..25}.root"
#       then
#           echo "Deleted output_WZJetsTo2L2Q_{1,21..25}.root"
#       fi
#   fi
#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/wz_2l2q
#   xrdcp -f wz_2l2q_*.root root://cmseos.fnal.gov/"$eosPath"/wz_2l2q/.

#   # WZ to 3l nu
#   if eval "hadd wz_3lnu_1.root output_WZJetsTo3LNu_*.root"
#   then
#       if eval "rm output_WZJetsTo3LNu_*.root"
#       then
#           echo "Deleted output_WZJetsTo3LNu_*.root"
#       fi
#   fi
#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/wz_3lnu
#   xrdcp -f wz_3lnu_*.root root://cmseos.fnal.gov/"$eosPath"/wz_3lnu/.



##  ZZ DIBOSON  ##

#   # ZZ to 2l 2q
#   if mv "output_ZZJetsTo2L2Q_1.root" "zz_2l2q_1.root"
#   then
#       echo "Moved output_ZZJetsTo2L2Q_1.root to zz_2l2q_1.root" 
#   fi
#   if eval "hadd zz_2l2q_2.root output_ZZJetsTo2L2Q_{2..9}.root"
#   then
#       if eval "rm output_ZZJetsTo2L2Q_{2..9}.root"
#       then
#           echo "Deleted output_ZZJetsTo2L2Q_{2..9}.root"
#       fi
#   fi
#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/zz_2l2q
#   xrdcp -f zz_2l2q_*.root root://cmseos.fnal.gov/"$eosPath"/zz_2l2q/.

#   # ZZ to 2l 2nu
#   if mv "output_ZZJetsTo2L2Nu_1.root" "zz_2l2nu_1.root"
#   then
#       echo "Moved output_ZZJetsTo2L2Nu_1.root to zz_2l2nu_1.root" 
#   fi
#   if eval "hadd zz_2l2nu_2.root output_ZZJetsTo2L2Nu_{2..5}.root output_ZZJetsTo2L2Nu_10.root"
#   then
#       if eval "rm output_ZZJetsTo2L2Nu_{2..5}.root output_ZZJetsTo2L2Nu_10.root"
#       then
#           echo "Deleted output_ZZJetsTo2L2Nu_{2..5,10}.root"
#       fi
#   fi
#   for i in {6..9}
#   do
#       j=$((i - 3))
#       if mv "output_ZZJetsTo2L2Nu_""$i"".root" "zz_2l2nu_""$j"".root"
#       then
#           echo "Moved output_ZZJetsTo2L2Nu_""$i"".root to zz_2l2nu_""$j"".root"
#       fi
#   done
#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/zz_2l2nu
#   xrdcp -f zz_2l2nu_*.root root://cmseos.fnal.gov/"$eosPath"/zz_2l2nu/.

#   # ZZ to 4l
#   if mv "output_ZZJetsTo4L_1.root" "zz_4l_1.root"
#   then
#       echo "Moved output_ZZJetsTo4L_1.root to zz_4l_1.root" 
#   fi
#   if eval "hadd zz_4l_2.root output_ZZJetsTo4L_{2..9}.root"
#   then
#       if eval "rm output_ZZJetsTo4L_{2..9}.root"
#       then
#           echo "Deleted output_ZZJetsTo4L_{2..9}.root"
#       fi
#   fi
#   eos root://cmseos.fnal.gov mkdir -p "$eosPath"/zz_4l
#   xrdcp -f zz_4l_*.root root://cmseos.fnal.gov/"$eosPath"/zz_4l/.

#   xrdcp -f source.tar.gz "$eosPath"/.
#   echo "Remember to delete!"
fi
