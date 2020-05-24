#!/bin/bash
path2rosetta='/projects/pw8/mc70/rosetta_bin_linux_2019.35.60890_bundle/'

## get native sequence from input PDB;
cat $1 | awk '/ATOM/ && $3 == "CA"  {print $4}' | tr '\n' ' ' | sed 's/ALA/A/g;s/CYS/C/g;s/ASP/D/g;s/GLU/E/g;s/PHE/F/g;s/GLY/G/g;s/HIS/H/g;s/ILE/I/g;s/LYS/K/g;s/LEU/L/g;s/MET/M/g;s/ASN/N/g;s/PRO/P/g;s/GLN/Q/g;s/ARG/R/g;s/SER/S/g;s/THR/T/g;s/VAL/V/g;s/TRP/W/g;s/TYR/Y/g' | sed 's/ //g' > native.seq

## randomize the native sequence and generate decoys;
python RandSeq.py $2 native.seq > peptide

## copy the input PDB into a dummy PDB id (for now, we us 3gso);
cp $1 3gso.pdb
cp $1 native.pdb

## evaluate energetics in the input pose with native sequence;
$path2rosetta/main/source/bin/rosetta_scripts.static.linuxgccrelease -parser:protocol  native.xml -s 3gso.pdb  -overwrite -ignore_zero_occupancy false > native.log
cat native.log | grep "ResResE" > temp
mv temp native.log

## evaluate energetics in the input pose threaded by decoy sequences;
let "j=0"
while read line; do
	let "j=j+1"
	sed -i "s/GKRSNTTGK/$line/" test.xml
	$path2rosetta/main/source/bin/rosetta_scripts.static.linuxgccrelease -parser:protocol  test.xml -s 3gso.pdb  -overwrite -ignore_zero_occupancy false > $j.log
	cat $j.log | grep "ResResE" > temp
	mv temp  $j.log
	mv 3gso_0001.pdb $j.pdb
	sed -i "s/$line/GKRSNTTGK/" test.xml
done < peptide

## Analyze the energetics and output frustration patterns of the input pose;
python Frust_Post_public.py 100 -2.5 0.5 $2 9 0 Function1 -1.5 0.5
