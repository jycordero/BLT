#! /bin/sh

set -e

batchDir=$1
year="2017"
srcPref="output"
batchPath="../test/batch/"
eosPath="/store/user/jrainbol/Trees/${year}/"

 # DATA
#   dataTag=("B" "C" "D" "E" "F")
#   nTarget=(2 4 2 4 5)

#   for dataName in "muon" "electron"
#   do
#       suffix=${dataName}_${year}

#       for ((i=0; i<${#dataTag[@]}; i+=1))
#       do
#           ./hadd_data_run.sh ${nTarget[i]} ${srcPref} ${suffix} ${dataTag[i]} ${batchPath}${batchDir}
#           echo ""
#           echo ""
#       done

#       ./move_data.sh ${suffix} ${suffix} 0 ${eosPath}${suffix}
#   done



 # MONTE CARLO
#   mcName=("DYJetsToLL_M-50" "TTJets" "WWTo2L2Nu" "WZTo2L2Q" "WZTo3LNu" "ZZTo2L2Q" "ZZTo4L" "GluGluHToZZTo4L" "VBF_HToZZTo4L")
#   nTarget=(5 3 1 4 2 4 2 1 1)
#   suffix=("zjets_m-50" "ttbar" "ww_2l2nu" "wz_2l2q" "wz_3lnu" "zz_2l2q" "zz_4l" "ggH_zz_4l" "vbfH_zz_4l")
    mcName=("TTJets")
    nTarget=(3)
    suffix=("ttbar")

    for ((i=0; i<${#mcName[@]}; i+=1))
    do
        ./hadd_move_mc.sh ${nTarget[i]} ${srcPref}_${mcName[i]} ${suffix[i]} ${batchPath}${batchDir} ${eosPath}${suffix[i]}
    done


echo ""
echo ""
echo "Remember to delete files in ${batchDir}!"
