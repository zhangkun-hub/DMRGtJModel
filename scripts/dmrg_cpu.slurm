#!/bin/bash
#SBATCH --job-name tJModelJ0.1Mu2D15000
#SBATCH --output=slurm-%j-%x.log
#SBATCH --mail-user=wanghaoxin1996@gmail.com
#SBATCH --mail-type=END
#SBATCH --partition 256G56c
#SBATCH --exclusive   
#SBATCH --ntasks 7
#SBATCH --cpus-per-task=8
#SBATCH --dependency=369159

ulimit -s unlimited 
ulimit -l unlimited

source ${HOME}/myenv.sh


rundir=$(pwd)
rootdir="${rundir}/../../"
exedir="${rootdir}/build_real_nou1/"
datadir="${rootdir}/data/"

for i in {1..4}; do
   paramsfile="params${i}.json"

   Ly=`cat ${paramsfile} | grep \"Ly\" | grep -o ':.*,' |grep -o [0-9].* | sed 's/,//g'`
   Lx=`cat ${paramsfile} | grep \"Lx\" | grep -o ':.*,' |grep -o [0-9].* | sed 's/,//g'`
   t2=`cat ${paramsfile} | grep \"t2\" | grep -o ':.*,' |grep -o [-,0-9].* | sed 's/,//g'`
   J=`cat ${paramsfile} | grep \"J\" | grep -o ':.*,' |grep -o [0-9].* | sed 's/,//g'`
   mu=`cat ${paramsfile} | grep mu | grep -o ':.*,' |grep -o [0-9].* | sed 's/,//g'`
   D=`cat ${paramsfile} | grep Dmax | grep -o ':.*,' |grep -o [0-9].* | sed 's/,//g'`

   #if [[ "${i}" == "1" ]]; then
   #   mpirun -n ${SLURM_NTASKS} ${exedir}/dmrg "${paramsfile}" --D=10,100,300,1000,2000,3000
   #else
      mpirun -n ${SLURM_NTASKS} ${exedir}/dmrg "${paramsfile}" 
   #fi
   cp -r mps "mpsD${D}"
   ${exedir}/measure1 ${paramsfile} 24 
   mpirun -n ${SLURM_NTASKS} ${exedir}/measureSC ${paramsfile} --start=24 --end=200 
   for name in $(ls *.json | grep -v "params"); do
     mv $name "${datadir}/${name%.*}${Ly}x${Lx}t2=${t2}J=${J}Mu${mu}D${D}Cylinder.json"
   done
done
