A script to run a benchmark with PCx for the given problems can be found in
`bench.py`. For informations of how to use it:

    $ bench.py --help

Some files, see `bench-*.out` have the output of `bench.py` when run it in a
machine with the hardware describe in `bench.hw`. For a pretty print of the
information in `bench-*.out` you can use

    $ column -s , -t bench-file.out

or a spreedsheet program.

Some script run a benchmark, filter and analysis the output for the problems in
the folder `mps` are given:

`bench.sh`
    Run the `bench.py` for some problems of the folder `mps`.
`filter.sh`
    Filter the output of `bench.sh` to get only some information.
`filter2tex.sh`
    Convert the `bench-*.join` to LaTeX table format.
`analysis.sh`
    Make a shallow analysis of the output of `filter.sh`.

    It depends of `analysis.awk`.
