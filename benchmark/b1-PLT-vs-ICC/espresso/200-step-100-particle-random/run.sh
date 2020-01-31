
#for var in 10
for var in 2 4 8 16 32
do
echo $var
/home/yougi/codes/espressomd/espresso4-1-2/build/pypresso icc-new.py $var
done
