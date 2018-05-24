#! /bin/sh

# Merge files
hadd    dy_m-10to50_1.root      output_DYJetsToLL_M-10to50_*.root

#hadd    dy_m-50_1.root          output_DYJetsToLL_M-50_?.root
#hadd    dy_m-50_2.root          output_DYJetsToLL_M-50_1?.root
#hadd    dy_m-50_3.root          output_DYJetsToLL_M-50_2?.root
#hadd    dy_m-50_4.root          output_DYJetsToLL_M-50_3?.root
#hadd    dy_m-50_5.root          output_DYJetsToLL_M-50_4?.root

#hadd    ttbar_1.root            output_ttbar_inclusive_*.root
#hadd    ttz_2l2nu_1.root        output_TTZToLLNuNu_M-10_*.root
#hadd    ww_2l2nu_1.root         output_WW_*.root
#hadd    wz_3lnu_1.root          output_WZJetsTo3LNu_*.root
#hadd    zz_4l_1.root            output_ZZJetsTo4L_*.root


#hadd    electron_2016_1.root    output_electron_2016B_*.root output_electron_2016C_*.root output_electron_2016D_*.root output_electron_2016E_*.root output_electron_2016F_*.root
#hadd    electron_2016_2.root    output_electron_2016G_*.root output_electron_2016H_*.root
#hadd    muon_2016_1.root        output_muon_2016B_*.root
#hadd    muon_2016_2.root        output_muon_2016C_*.root output_muon_2016D_*.root
#hadd    muon_2016_3.root        output_muon_2016E_*.root output_muon_2016F_*.root
#hadd    muon_2016_4.root        output_muon_2016G_*.root
#hadd    muon_2016_5.root        output_muon_2016H_*.root


# Move to EOS
#for suffix in "dy_m-10to50" "dy_m-50" "ttbar" "ttz_2l2nu" "ww_2l2nu" "wz_3lnu" "zz_4l" "electron_2016" "muon_2016"
for suffix in "dy_m-10to50"
do
    xrdcp "$suffix"_*.root root://cmseos.fnal.gov//store/user/jrainbol/Trees/2016/"$suffix"/.
done
