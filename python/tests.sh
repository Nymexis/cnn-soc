#!/bin/bash 
FILES="./tests/*"
RESULT="./results.txt"
EXCLUDED="fig1.png dl.sh fig2.png"

if [ -f $RESULT ]; then
	rm $RESULT
	echo "Deleted previous results file."
else
	echo "No previous results file found."
fi

touch results.txt

for ((i = 1; i <= $1; i++ )); do
  python3 main.py $i
done

RES=`cat results.txt | grep 1 | wc -l`
TOTAL=`cat results.txt | wc -l`
echo "$RES"
