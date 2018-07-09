#! /bin/sh


# Don't try to do anything if the wrong version of ROOT is loaded
myROOTSYS="/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.10.02-19565/x86_64-slc6-gcc62-opt"
if [ "$ROOTSYS" != "$myROOTSYS" ]
then
    echo "ROOT 6.10 not loaded"
else


## Rename data trees and histograms
#   for suffix in "electron_2016" "muon_2016"
#   do
#       for tag in B C D E F G H
#       do
#           hadd    "$suffix$tag".root    output_"$suffix$tag"_*.root
#           root.exe -q -b "../../renameTree.C(\""$suffix$tag".root\", \"tree_"$suffix$tag"\", \"tree_"$suffix"\")"
#           if rootrm  "$suffix$tag"".root:tree_""$suffix$tag"
#           then
#               echo Deleted tree_"$suffix$tag".
#           fi
#           if rootmv  "$suffix$tag"".root:TotalEvents_""$suffix$tag"    "$suffix$tag"".root:TotalEvents_""$suffix"
#           then
#               echo Moved TotalEvents_"$suffix$tag" to TotalEvents_"$suffix".
#               echo ""
#               echo ""
#               echo ""
#           fi
#       done
#   done

## Rename files
#   mv      electron_2016B.root     electron_2016_1.root
#   mv      electron_2016C.root     electron_2016_2.root
#   mv      electron_2016D.root     electron_2016_3.root
#   mv      electron_2016E.root     electron_2016_4.root
#   mv      electron_2016F.root     electron_2016_5.root
#   mv      electron_2016G.root     electron_2016_6.root
#   mv      electron_2016H.root     electron_2016_7.root

#   mv      muon_2016B.root         muon_2016_1.root
#   mv      muon_2016C.root         muon_2016_2.root
#   mv      muon_2016D.root         muon_2016_3.root
#   mv      muon_2016E.root         muon_2016_4.root
#   mv      muon_2016F.root         muon_2016_5.root
#   mv      muon_2016G.root         muon_2016_6.root
#   mv      muon_2016H.root         muon_2016_7.root


 # Merge remaining files
    hadd    dy_m-10to50_1.root      output_DYJetsToLL_M-10to50_*.root

    hadd    dy_m-50_1.root          output_DYJetsToLL_M-50_1.root   output_DYJetsToLL_M-50_?1.root
    hadd    dy_m-50_2.root          output_DYJetsToLL_M-50_2.root   output_DYJetsToLL_M-50_?2.root
    hadd    dy_m-50_3.root          output_DYJetsToLL_M-50_3.root   output_DYJetsToLL_M-50_?3.root
    hadd    dy_m-50_4.root          output_DYJetsToLL_M-50_4.root   output_DYJetsToLL_M-50_?4.root
    hadd    dy_m-50_5.root          output_DYJetsToLL_M-50_5.root   output_DYJetsToLL_M-50_?5.root
    hadd    dy_m-50_6.root          output_DYJetsToLL_M-50_6.root   output_DYJetsToLL_M-50_?6.root
    hadd    dy_m-50_7.root          output_DYJetsToLL_M-50_7.root   output_DYJetsToLL_M-50_?7.root
    hadd    dy_m-50_8.root          output_DYJetsToLL_M-50_8.root   output_DYJetsToLL_M-50_?8.root
    hadd    dy_m-50_9.root          output_DYJetsToLL_M-50_9.root   output_DYJetsToLL_M-50_?9.root

    hadd    ggH_zz_4l_1.root        output_GluGluHToZZTo4L_*.root
    hadd    H_zz_4l_1.root          output_HToZZTo4L_*.root

    hadd    ttbar_1.root            output_ttbar_inclusive_?.root
    hadd    ttbar_2.root            output_ttbar_inclusive_1?.root
    hadd    ttbar_3.root            output_ttbar_inclusive_2?.root
    hadd    ttbar_4.root            output_ttbar_inclusive_3?.root

    hadd    ttz_2l2nu_1.root        output_TTZToLLNuNu_M-10_*.root
    hadd    ww_2l2nu_1.root         output_WW_*.root
    hadd    wz_3lnu_1.root          output_WZJetsTo3LNu_*.root
    hadd    zz_4l_1.root            output_ZZJetsTo4L_*.root


 # Move to EOS
#   for suffix in "electron_2016" "muon_2016" "dy_m-10to50" "dy_m-50" "ggH_zz_4l" "H_zz_4l" "ttbar" "ttz_2l2nu" "ww_2l2nu" "wz_3lnu" "zz_4l"
    for suffix in "dy_m-10to50" "dy_m-50" "ggH_zz_4l" "H_zz_4l" "ttbar" "ttz_2l2nu" "ww_2l2nu" "wz_3lnu" "zz_4l"
    do
        eos root://cmseos.fnal.gov mkdir -p /store/user/jrainbol/Trees/2016/"$suffix"
        xrdcp       -f  "$suffix"_*.root    root://cmseos.fnal.gov//store/user/jrainbol/Trees/2016/"$suffix"/.
    done

    xrdcp   -f  source.tar.gz   root://cmseos.fnal.gov//store/user/jrainbol/Trees/2016/.

    echo "Remember to delete!"
fi
