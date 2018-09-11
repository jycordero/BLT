#! /bin/sh

set -e

batchDir=$1
year="2017"
srcPref="output"
batchPath="../test/batch/"
eosPath="/store/user/jrainbol/Trees/${year}/"

 # DATA
   dataName=("muon" "electron")
    dataTag=("B" "C" "D" "E" "F")
    nTarget=(2 4 2 4 5)

    for name in ${dataName}
    do
        suffix=${name}_${year}

        for ((i=0; i<${#dataTag[@]}; i+=1))
        do
            ./hadd_data_run.sh ${nTarget[i]} ${srcPref} ${suffix} ${dataTag[i]} ${batchPath}${batchDir}
            echo ""
            echo ""
        done

        ./move_data.sh ${suffix} ${suffix} 0 ${eosPath}${suffix}
    done



 # MONTE CARLO
    mcName=("TTJets" "WWTo2L2Nu" "WZTo2L2Q" "WZTo3LNu" "ZZTo2L2Q" "ZZTo4L" "GluGluHToZZTo4L" "VBF_HToZZTo4L")
    nTarget=(3 1 4 2 4 2 1 1)
    suffix=("ttbar" "ww_2l2nu" "wz_2l2q" "wz_3lnu" "zz_2l2q" "zz_4l" "ggH_zz_4l" "H_zz_4l")

    for ((i=0; i<${#mcName[@]}; i+=1))
    do
        ./hadd_move_mc.sh ${nTarget[i]} ${srcPref}_${mcName[i]} ${suffix[i]} ${batchPath}${batchDir} ${eosPath}${suffix[i]}
    done


echo ""
echo ""
echo "Remember to delete files in ${batchDir}!"
