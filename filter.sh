#!/bin/bash
# Need for the sort.
LC_ALL=C
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
            echo 'Filter benchmakr log get from `bench.py`'
            echo ''
            echo 'Usage: filter.sh [-l lib] [-s specs]'
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

function specs_join(){
    echo ","$(for s in ${specs};do printf "%s,,,," $s; done;) > bench-${1}.join
    join -j 1 -t , $(for s in ${specs}; do echo bench-${1}-${s}.out-sort; done) >> bench-${1}.join
}

function bench_sort(){
    head -n 1 ${1}-select > ${1}-sort;
    tail -n +2 ${1}-select | sort >> ${1}-sort;
}

function specs_sort(){
    for s in ${specs}
    do
        bench_sort bench-${1}-${s}.out
    done
}

function bench_select(){
    awk -F , '{printf "%s,%s,%.4E,%s,%.4E\n", $1, $2, $7, $14, $15}' ${1} > ${1}-select
}

function specs_select(){
    for s in ${specs}
    do
        bench_select bench-${1}-${s}.out
    done
}

for l in ${lib}
do
    # Select only part of the fields from the log
    specs_select ${l}
    # Sort the problems by name
    specs_sort ${l}
    # Join the specs
    specs_join ${l}
done
