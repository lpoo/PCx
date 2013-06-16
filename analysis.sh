#!/bin/bash
stdlib="netlib netlibqap kenlib meslib pdslib raillib fctplib"
while getopts h name
do
    case $name in
        h)
            echo "Give some statistics about some library based on."
            echo "the output of filter.sh.";&
        ?)
            echo 'Usage: analysis.sh lib ...'
            echo ''
            echo 'Avaliable test library:'
            echo ''
            echo '    * netlib'
            echo '    * netlibqap'
            echo '    * kenlib'
            echo '    * meslib'
            echo '    * pdslib'
            echo '    * raillib'
            echo '    * fctplib'
            echo ''
            echo 'If any argument is given it will run for all the library'
            echo 'and specifications avaliable.'
            exit 2;;
    esac
done

if $(test ! -v lib)
then
    lib=$stdlib
elif $(test -z $lib)
then
    lib=$stdlib
fi

for l in $lib
do
    echo "Check if file \`bench-${l}.join\` exist."
    if $(test ! -e bench-${l}.join)
    then
        echo "Can\'t find file \`bench-${l}.join\`. Aborting."
        exit 2;
    fi
    echo "Found file \`bench-${l}.join\`. Starting analysis."
    tail -n +3  bench-${l}.join | awk -F ',' -f analysis.awk
done
