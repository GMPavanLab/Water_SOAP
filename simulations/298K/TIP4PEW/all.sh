source /Users/cape/programs/gromacs-2021.2/bin/GMXRC
#source /home/pavan/Riccardo/GMX212/bin/GMXRC
#source /home/cape/programs/GMX21_2/bin/GMXRC

rm */*.tpr

cd minimization

gmx_mpi grompp -f min.mdp -c ../1024_wat.gro -p ../topol.top -o min.tpr -maxwarn 3
mpirun -n 2 gmx_mpi mdrun  -v -deffnm min -nsteps -1
rm \#*
rm min.trr

cd ..

cd eq_nvt 

gmx_mpi grompp -f eq_nvt.mdp -c ../minimization/min.gro -p ../topol.top -o eq_nvt.tpr
mpirun -n 2 gmx_mpi mdrun  -v -deffnm eq_nvt
rm \#*

cd ..

cd eq_npt

gmx_mpi grompp -f eq_npt.mdp -c ../eq_nvt/eq_nvt.gro -p ../topol.top -o eq_npt.tpr
mpirun -n 2 gmx_mpi mdrun  -v -deffnm eq_npt
rm \#*

cd ..

cd run

gmx_mpi grompp -f run.mdp -c ../eq_npt/eq_npt.gro -p ../topol.top -o run.tpr
mpirun -n 2 gmx_mpi mdrun  -v -deffnm run
rm \#*

cd ..
