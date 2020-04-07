cp 3gso.pdb native.pdb 
$path2rosetta/main/source/bin/rosetta_scripts.static.linuxgccrelease -parser:protocol  native.xml -s 3gso.pdb  -extra_res_fa UNK.params   -overwrite -ignore_zero_occupancy false > native.log
cat native.log | grep "ResResE" > temp
mv temp native.log

let "j=0"
while read line; do
	let "j=j+1"
	sed -i "s/GKRSNTTGK/$line/" test.xml
	$path2rosetta/main/source/bin/rosetta_scripts.static.linuxgccrelease -parser:protocol  test.xml -s 3gso.pdb  -extra_res_fa UNK.params   -overwrite -ignore_zero_occupancy false > $j.log
	cat $j.log | grep "ResResE" > temp
	mv temp  $j.log
	mv 3gso_0001.pdb $j.pdb
	sed -i "s/$line/GKRSNTTGK/" test.xml
done < peptide
python Frust_Post.py 100 -2.5 0.5 300 9 0 Function1 -1.5 0.5
