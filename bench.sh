#!/bin/bash
stdlib="netlib netlibqap kenlib meslib pdslib raillib fctplib"
stdspecs="mmd rcm"
while getopts l:s:h name
do
    case $name in
        (l)
            lib="$OPTARG";;
        (s)
            specs="$OPTARG";;
        (h|?)
            echo 'Call `bench.py` for test libray and specifications given'
            echo 'as arguments.'
            echo ''
            echo 'Usage: bench.sh [-l lib] [-s specs]'
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
            echo 'Avaliable specifications:'
            echo ''
            echo '    * mmd'
            echo '    * rcm'
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

if $(test ! -v specs)
then
    specs=$stdspecs
elif $(test -z $specs)
then
    specs=$stdspecs
fi

for l in $lib
do
    for s in $specs
    do
        python3 bench.py -K -o bench-$l-$s.out -S $s.specs -c --$l
    done
done
