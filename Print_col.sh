#!/bin/bash
result=0
for ii in {2..55}
do
	infile="$(cat $1 | head -n 1 | tr -s ' '| cut -d ' ' -f $ii)"
	if [[ "$infile" = "$2_real" ]]
	then
		result=$ii
		break
	fi
done
if [[ "$result" -eq 0 ]]
then
	echo "Error: Index not found!"
	exit 1
fi
next="$(( $result + 1 ))"
cat $1 | tr -s ' '| cut -d ' ' -f 1,$result,$next
exit 0

