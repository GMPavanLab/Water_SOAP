#source /Users/cape/programs/gromacs-2021.2/bin/GMXRC
source /home/pavan/Riccardo/GMX212/bin/GMXRC

rm */*.tpr

cd minimization

gmx grompp -f min.mdp -c ../1024_wat.gro -p ../topol.top -o min.tpr -maxwarn 3
gmx mdrun -nt 6 -v -deffnm min -nsteps -1
rm \#*

cd ..

cd eq_nvt 

gmx grompp -f eq_nvt.mdp -c ../minimization/min.gro -p ../topol.top -o eq_nvt.tpr
gmx mdrun -nt 6 -v -deffnm eq_nvt
rm \#*

cd ..

cd eq_npt

gmx grompp -f eq_npt.mdp -c ../eq_nvt/eq_nvt.gro -p ../topol.top -o eq_npt.tpr
gmx mdrun -nt 6 -v -deffnm eq_npt
rm \#*

cd ..

cd run

gmx grompp -f run.mdp -c ../eq_npt/eq_npt.gro -p ../topol.top -o run.tpr
gmx mdrun -nt 6 -v -deffnm run
rm \#*

cd ..
