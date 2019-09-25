#!/bin/bash
echo '' > results.dat
for a in {-5..5}
do
for b in {-5..5}
do
for c in {-5..5}
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

echo ${a}_${b}_$c
mkdir A_${a}_${b}_$c
cp config test contral neural.in  A_${a}_${b}_$c
cd A_${a}_${b}_$c
sed  -i  "s/ELEMENTS_WEIGHT=/ELEMENTS_WEIGHT= $a $b ${c}/" contral
~/nGAP/src/ngap.x > log
e=`grep 'ENERGY' log | awk '{print $3}'`
f=`grep 'FORCE' log | awk '{print $3}'`
s=`grep 'STRESS' log | awk '{print $5}'`
echo ${a}_${b}_$c  $e $f $s >> ../results.dat
cd ..
done
done
done
