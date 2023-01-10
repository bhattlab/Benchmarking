#!/bin/bash
FILE="$1"

info=$(cat $FILE | cut -f 3,4,5,6 | awk '{for(i=1;i<=NF;i++)$i=(a[i]+=$i)}END{print}')
echo -e $FILE"\t"$info

