#1:channel{e,mu},2:dataset{DATA,MC},3:max_thread_number,4:output_file path and name

channel="0"
dataset="0"

if test $1 = "e"
then
    channel="Zee"
elif test $1 = "mu"
then
    channel="MuMu"
else
    echo "wrong channel"
    exit
fi

if test $2 = "MC"
then
    dataset=$2
elif test $2 = "DATA"
then
    dataset=$2
else
    echo "wrong dataset"
    exit
fi

rm -rf ../$channel/$dataset/out
mkdir ../$channel/$dataset/out/


cd ../$channel
./SR_prepare.sh
python3 Compiler.py $dataset

cd ./$dataset
parallel -j $3 --progress -a ../../parallelRun/input/$1_$2.txt python3 RunAnalysis.py

cd ../../parallelRun

hadd ./root_files/$4 ../$channel/$dataset/out/*.root
