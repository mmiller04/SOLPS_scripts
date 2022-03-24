#!/bin/bash

if [ $# -ne 2 ]
then
        read -p 'Shot Number: ' SN

        read -p 'Attempt Number: ' -a ARR
else

        SN=$1
        ARR=$2
fi

echo ${#ARR[@]} Attempts Entered...

for AN in ${ARR[@]}
do

echo Beginning data collection for Attempt$AN ...

cd ../$SN/attempt$AN

mkdir -p Output/

rm -rfv Output/*

2d_profiles

cp {an3da*,dn3da*,ke3da*,ki3da*,ne3da*,ne3dr*,te3da*,te3dr*,ti3da*,ti3dr*} Output/

echo "cr 'RadLoc${AN}' mprt" | b2plot

mv b2pl.exe.dir/RadLoc${AN} Output/

echo "cz 'VertLoc${AN}' mprt" | b2plot

mv b2pl.exe.dir/VertLoc${AN} Output/

echo "vol 'VOL${AN}' mprt" | b2plot

mv b2pl.exe.dir/VOL${AN} Output/

echo "d 'DN${AN}' mprt" | b2plot

mv b2pl.exe.dir/DN${AN} Output/

echo "kye 'KYE${AN}' mprt" | b2plot

mv b2pl.exe.dir/KYE${AN} Output/

echo "kyi 'KYI${AN}' mprt" | b2plot

mv b2pl.exe.dir/KYI${AN} Output/

echo "fnay 'IonFlx${AN}' mprt" | b2plot

mv b2pl.exe.dir/IonFlx${AN} Output/

echo "fnax 'IonPol${AN}' mprt" | b2plot

mv b2pl.exe.dir/IonPol${AN} Output/

echo "eirc dab2 'NeuDen${AN}' mprt" | b2plot

mv b2pl.exe.dir/NeuDen${AN} Output/

echo "eirc dmb2 'MolDen${AN}' mprt" | b2plot

mv b2pl.exe.dir/MolDen${AN} Output/

echo "eirc tab2 'NeuTemp${AN}' mprt" | b2plot

mv b2pl.exe.dir/NeuTemp${AN} Output/

echo "eirc tmb2 'MolTemp${AN}' mprt" | b2plot

mv b2pl.exe.dir/MolTemp${AN} Output/

echo "te 'Te${AN}' mprt" | b2plot

mv b2pl.exe.dir/Te${AN} Output/

echo "ti 'Ti${AN}' mprt" | b2plot

mv b2pl.exe.dir/Ti${AN} Output/

echo "ne 'Ne${AN}' mprt" | b2plot

mv b2pl.exe.dir/Ne${AN} Output/

echo "vlay 'RadPinch${AN}' mprt" | b2plot

mv b2pl.exe.dir/RadPinch${AN} Output/

2da timesa > Output/TimeStamps

2da nesepm tesepm tisepm nesepa tesepa tisepa nesepi tesepi tisepi> Output/TimeTraces

if [ -f "Note" ]; then
	cp Note Output/Note
fi

echo Output files created for Attempt${AN}



echo Now starting GeoMeshGen for Output2

mkdir -p Output2/

rm -rfv Output2/*

echo "0 crx 'Rad0Cor${AN}' mprt" | b2plot

mv b2pl.exe.dir/Rad0Cor${AN} Output2/

echo "0 cry 'Vert0Cor${AN}' mprt" | b2plot

mv b2pl.exe.dir/Vert0Cor${AN} Output2/

echo "1 crx 'Rad1Cor${AN}' mprt" | b2plot

mv b2pl.exe.dir/Rad1Cor${AN} Output2/

echo "1 cry 'Vert1Cor${AN}' mprt" | b2plot

mv b2pl.exe.dir/Vert1Cor${AN} Output2/

echo "2 crx 'Rad2Cor${AN}' mprt" | b2plot

mv b2pl.exe.dir/Rad2Cor${AN} Output2/

echo "2 cry 'Vert2Cor${AN}' mprt" | b2plot

mv b2pl.exe.dir/Vert2Cor${AN} Output2/

echo "3 crx 'Rad3Cor${AN}' mprt" | b2plot

mv b2pl.exe.dir/Rad3Cor${AN} Output2/

echo "3 cry 'Vert3Cor${AN}' mprt" | b2plot

mv b2pl.exe.dir/Vert3Cor${AN} Output2/

echo "hx 'HX${AN}' mprt" | b2plot

mv b2pl.exe.dir/HX${AN} Output2/

echo "hy 'HY${AN}' mprt" | b2plot

mv b2pl.exe.dir/HY${AN} Output2/

echo "sx 'SX${AN}' mprt" | b2plot

mv b2pl.exe.dir/SX${AN} Output2/

echo "sy 'SY${AN}' mprt" | b2plot

mv b2pl.exe.dir/SY${AN} Output2/

echo "sz 'SZ${AN}' mprt" | b2plot

mv b2pl.exe.dir/SZ${AN} Output2/

echo Output2 files created for Attempt${AN}

cd ..

done

echo ALL REQUESTED OUTPUT DATA COLLECTED AND UPLOADED
