n_replicas=20
name="9_10_new_20x4E5" #"9_10_new_20x1E5" #output folder name
cwd=$(pwd)

for i in `seq 0 $n_replicas`; do
    folder=../re-submission/results_bootstrap/$name/rep_$i
    mkdir $folder
    cp *py *x $folder
    cd $folder
    nohup python read_fastq.py > logfile.$i &
    cd $cwd
done
