#!/bin/bash
stdlib="netlib netlibqap kenlib meslib pdslib raillib fctplib"
stdlang='en'
while getopts l:h name
do
    case $name in
        l)
            lang="$OPTARG";;
        h)
            echo "Convert \`bench-*.out\` to LaTeX table format.";&
        ?)
            echo 'Usage: analysis.sh [-l lang] lib ...'
            echo ''
            echo 'Avaliable languages:'
            echo ''
            echo '    * en (default)'
            echo '    * pt'
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

if $(test ! -v lang)
then
    lang=$stdlang
elif $(test -z $lang)
then
    lang=$stdlang
fi

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

    echo "Generating LaTeX table(s) from \`bench-${l}.join\`."
    ifilename=bench-${l}.join
    ofilename=bench-${l}.join.tex

    if $(test $lang = 'en')
    then
        echo '\begin{tabular}{|l|r|r|r|r|r|r|r|r|}' > $ofilename
        echo '\hline' >> $ofilename
        head -n 1 $ifilename | sed 's/,/ /g' | awk '{print \
        "\\multicolumn{1}{|c|}{Problem} & \\multicolumn{4}{|c|}{" toupper($1) "} & \
        \\multicolumn{4}{|c|}{" toupper($2) "} \\\\ \\hline"}' >> $ofilename
        echo '\multicolumn{1}{|c|}{Name} & \multicolumn{1}{|c|}{R} &
        \multicolumn{1}{|c|}{NNZ} & \multicolumn{1}{|c|}{IT} &
        \multicolumn{1}{|c|}{T} & \multicolumn{1}{|c|}{R} &
        \multicolumn{1}{|c|}{NNZ} & \multicolumn{1}{|c|}{IT} &
        \multicolumn{1}{|c|}{T} \\ \hline' >> $ofilename
        tail -n +3 $ifilename | sed 's/^[a-z]*-//; s/\.mps//; s/_/-/g; s/,/ \& /g;
        s/$/ \\\\ \\hline/' >> $ofilename
        echo '\end{tabular}' >> $ofilename
    fi

    if $(test $lang = 'pt')
    then
        echo '\begin{tabular}{|l|r|r|r|r|r|r|r|r|}' > $ofilename
        echo '\hline' >> $ofilename
        head -n 1 $ifilename | sed 's/,/ /g' | awk '{print \
        "\\multicolumn{1}{|c|}{Problema} & \\multicolumn{4}{|c|}{" toupper($1) "} & \
        \\multicolumn{4}{|c|}{" toupper($2) "} \\\\ \\hline"}' >> $ofilename
        echo '\multicolumn{1}{|c|}{Nome} & \multicolumn{1}{|c|}{R} &
        \multicolumn{1}{|c|}{NNZ} & \multicolumn{1}{|c|}{IT} &
        \multicolumn{1}{|c|}{T} & \multicolumn{1}{|c|}{R} &
        \multicolumn{1}{|c|}{NNZ} & \multicolumn{1}{|c|}{IT} &
        \multicolumn{1}{|c|}{T} \\ \hline' >> $ofilename
        tail -n +3 $ifilename | sed 's/^[a-z]*-//; s/\.mps//; s/_/-/g; s/,/ \& /g; s/\./,/g;
        s/$/ \\\\ \\hline/' >> $ofilename
        echo '\end{tabular}' >> $ofilename
    fi

    # Split file if it has more than 40 lines
    if $(test $(( ($(wc -l < $ofilename) - 9) / 40 + 1)) -gt 1)
    then
        header=$(head -n 8 $ofilename)
        footer=$(tail -n 1 $ofilename)
        head -n -1 $ofilename | tail -n +5 | split -a 1 -d -l 40 - x
        for i in $(seq 0 $(( ($(wc -l < $ofilename) - 9) / 40)) )
        do
            echo $header > ${ofilename}$i
            cat x$i >> ${ofilename}$i
            echo $footer >> ${ofilename}$i
        done
    fi
    # Remove temporary files
    rm -f x[0-9]*
done
