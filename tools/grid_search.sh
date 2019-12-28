#!/bin/bash
#echo '' > results.dat
echo 'COMPS'  'RMSE_E' 'RMSE_F' 'RMSE_S'  | awk '{printf "%10s %15s %15s %15s\n", $1, $2, $3,$4}' > results.dat
for a in {-3..3}
do
for b in {-3..3}
do
for c in {-3..3}
do
if [ $a -eq 0  ] 
then
    continue
fi

if [ $b -eq 0  ]
then
    continue
fi

if [ $c -eq 0  ] 
then
    continue
fi
if [ $a -eq $b ]
then
    continue
fi
if [ $a -eq $c ]
then
    continue
fi
if [ $c -eq $b ]
then
    continue
fi

echo ${a}_${b}_$c
mkdir A_${a}_${b}_$c
cp config test control neural.in  A_${a}_${b}_$c
cd A_${a}_${b}_$c
sed  -i  "s/WeightOfAtoms=/WeightOfAtoms= $a $b $c/" control
~/NGAP_RELEASE/ngap/ngap.x > log
e=`grep 'ENERGY' log | awk '{print $3}'`
f=`grep 'FORCE' log | awk '{print $3}'`
s=`grep 'STRESS' log | awk '{print $5}'`
echo ${a}_${b}_$c  $e $f $s  | awk '{printf "%10s %15.3f %15.2f %15.1f\n", $1, $2, $3,$4}' >> ../results.dat
cd ..
done
done
done
