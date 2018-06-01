#! /bin/sh


# Rename data trees and histograms
    for suffix in "electron_2016" "muon_2016"
    do
        for tag in B C D E F G H
        do
            hadd    "$suffix$tag".root    output_"$suffix$tag"_*.root
            root.exe -q -b "../../renameTree.C(\""$suffix$tag".root\", \"tree_"$suffix$tag"\", \"tree_"$suffix"\")"
            if rootrm  "$suffix$tag".root:tree_"$suffix$tag"
            then
                echo Deleted "$suffix$tag".
            fi
            if rootmv  "$suffix$tag".root:TotalEvents_"$suffix$tag"    "$suffix$tag".root:TotalEvents_"$suffix"
            then
                echo Moved TotalEvents_"$suffix$tag" to TotalEvents_"$suffix".
            fi
        done
    done


# Merge files
    hadd    muon_2016_1.root        muon_2016B.root
    hadd    muon_2016_2.root        muon_2016C.root muon_2016D.root
    hadd    muon_2016_3.root        muon_2016E.root muon_2016F.root
    hadd    muon_2016_4.root        muon_2016G.root
    hadd    muon_2016_5.root        muon_2016H.root

    hadd    electron_2016_1.root    electron_2016B.root electron_2016C.root electron_2016D.root electron_2016E.root electron_2016F.root
    hadd    electron_2016_2.root    electron_2016G.root electron_2016H.root

#   hadd    dy_m-10to50_1.root      output_DYJetsToLL_M-10to50_*.root
#   hadd    dy_m-50_1.root          output_DYJetsToLL_M-50_?.root
#   hadd    dy_m-50_2.root          output_DYJetsToLL_M-50_1?.root
#   hadd    dy_m-50_3.root          output_DYJetsToLL_M-50_2?.root
#   hadd    dy_m-50_4.root          output_DYJetsToLL_M-50_3?.root
#   hadd    dy_m-50_5.root          output_DYJetsToLL_M-50_4?.root

    hadd    zjets_m-10to50_1.root   output_DYJetsToLL_M-10to50_*.root

    hadd    zjets_m-50_1.root       output_DYJetsToLL_M-50_1.root output_DYJetsToLL_M-50_2.root output_DYJetsToLL_M-50_3.root output_DYJetsToLL_M-50_4.root output_DYJetsToLL_M-50_5.root
    hadd    zjets_m-50_2.root       output_DYJetsToLL_M-50_6.root output_DYJetsToLL_M-50_7.root output_DYJetsToLL_M-50_8.root output_DYJetsToLL_M-50_9.root output_DYJetsToLL_M-50_10.root

    hadd    ttbar_1.root            output_ttbar_inclusive_*.root
    hadd    ttz_2l2nu_1.root        output_TTZToLLNuNu_M-10_*.root
    hadd    ww_2l2nu_1.root         output_WW_*.root
    hadd    wz_3lnu_1.root          output_WZJetsTo3LNu_*.root
    hadd    zz_4l_1.root            output_ZZJetsTo4L_*.root


# Move to EOS
#   for suffix in "dy_m-10to50" "dy_m-50" "ttbar" "ttz_2l2nu" "ww_2l2nu" "wz_3lnu" "zz_4l" "electron_2016" "muon_2016"
    for suffix in "zjets_m-10to50" "zjets_m-50" "ttbar" "ttz_2l2nu" "ww_2l2nu" "wz_3lnu" "zz_4l" "electron_2016" #"muon_2016"
    do
        xrdcp -f "$suffix"_*.root root://cmseos.fnal.gov//store/user/jrainbol/Trees/2016/"$suffix"/.
    done

    echo "Remember to delete!"
