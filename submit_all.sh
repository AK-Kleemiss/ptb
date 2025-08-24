wd=`pwd`
rundir="runs"
for filter in 0.001 0.0001 0.00001 0.000001 0.0;
#for filter in 0.0001;
do
  #for col in full_matrix kmeans_r_heuristic atoms_columns;
  #for col in full_matrix kmeans_r_heuristic atoms_columns conventional;
  for col in full_matrix kmeans_r_heuristic;
  do
    cd $wd
    for sys in `ls systems/*.xyz | sed "s/\.xyz//g" | sed "s/systems\///g"`;
    #for sys in h2o_512;
    #for sys in h2o_512 nylon 1nx2 2b59;
    do
        cd $wd
        charge=`cat systems/$sys.charge`
        dir="$rundir/${sys}_$charge/filter_$filter/$col"
        echo $dir
        rm -rf $dir
        mkdir -p $dir
        cp job.sh $dir
        cd $dir
        sed -i "s/_SYS_/$sys/g" job.sh
        sed -i "s/_CHRG_/$chrg/g" job.sh
        sed -i "s/_FILTER_/$filter/g" job.sh
        if [ "$col" != "conventional" ];
        then
          echo "mode = submatrix" > .PUR
          #echo "check" >> .PUR
          echo "submatrix_mode = submatrix_sygv" >> .PUR
          echo "submatrix_columns = $col " >> .PUR
          echo "check " >> .PUR
          sbatch job.sh -check -purify
        else
          sbatch job.sh -check
        fi

    done  
  done
done
