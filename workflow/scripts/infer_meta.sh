#!/bin/bash

if [ "${1}" == "endness" ]
then
    awk '/Data/ {print}' "${2}" | sed -e 's/^This is //' -e 's/ Data$//'
elif [ "${1}" == "fail" ]
then
    awk '/Fraction of reads failed/ {print}' "${2}" | sed -e 's/^Fraction of reads failed to determine: //'
elif [ "${1}" == "sef" ]
then
    awk '/\+\+,--/ {print}' "${2}" | sed -e 's/^Fraction of reads explained by "++,--": //'
elif [ "${1}" == "ser" ]
then
    awk '/\+-,-\+/ {print}' "${2}" | sed -e 's/^Fraction of reads explained by "+-,-+": //'
elif [ "${1}" == "pef" ]
then
    awk '/1\+\+,1--,2\+-,2-\+/ {print}' "${2}" | sed -e 's/^Fraction of reads explained by "1++,1--,2+-,2-+": //'
elif [ "${1}" == "per" ]
then
    awk '/1\+-,1-\+,2\+\+,2--/ {print}' "${2}" | sed -e 's/^Fraction of reads explained by "1+-,1-+,2++,2--": //'
fi