#! /bin/sh

set -e

batchDir=$1
year="2017"
srcPref="output"
batchPath="../test/batch/"
eosPath="/store/user/jrainbol/Trees/${year}/"

 # DATA
#   dataName=("muon", "eg")
    dataName=("electron")
    dataTag=("B" "C" "D" "E" "F")
    nTarget=(1 4 2 4 5)

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



## MONTE CARLO
#   mcName=("DYJetsToLL_M-50")
#   nTarget=(5)
#   suffix=("zjets_m-50")

#   for ((i=0; i<${#mcName[@]}; i+=1))
#   do
#       ./hadd_move_mc.sh ${nTarget[i]} ${srcPref}_${mcName} ${suffix} ${batchPath}${batchDir} ${eosPath}${suffix}
#   done


echo ""
echo ""
echo "Remember to delete files in ${batchDir}!"
