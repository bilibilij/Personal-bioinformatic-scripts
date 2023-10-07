#module load openmpi/1.6.5
MPI=$4
if [ $MPI == "True" ]; then
  CMD=$(echo mpirun -np $1 $2 $3)
else
 CMD=$(echo $2 $3)
fi
echo $CMD
$CMD