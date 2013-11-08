#!/usr/bin/python
import os
import re

class PTInfo():
    def __init__(self, name):
        self.name = str(name)
        self.brows = 0
        self.bcols = 0
        self.arows = 0
        self.acols = 0
        self.fact_nonzeros = 0
        self.fact_density = 0
        self.rel_infeas = 0
        self.rel_compl = 0
        self.prim_obj = 0
        self.dual_obj = 0
        self.max_corr = 0
        self.it = 0
        self.cpu = 0

    def add_status(self, status):
        self.status = int(status)

    def add_brows(self, rows):
        self.brows = int(rows)

    def add_bcols(self, cols):
        self.bcols = int(cols)

    def add_arows(self, rows):
        self.arows = int(rows)

    def add_acols(self, cols):
        self.acols = int(cols)

    def add_fact_nonzeros(self, nonzeros):
        self.fact_nonzeros = int(nonzeros)

    def add_fact_density(self, density):
        self.fact_density = float(density)

    def add_rel_infeas(self, rel_infeas):
        self.rel_infeas = float(rel_infeas)

    def add_rel_compl(self, rel_compl):
        self.rel_compl = float(rel_compl)

    def add_prim_obj(self, prim_obj):
        self.prim_obj = float(prim_obj)

    def add_dual_obj(self, dual_obj):
        self.dual_obj = float(dual_obj)

    def add_max_corr(self, max_corr):
        self.max_corr = int(max_corr)

    def add_it(self, it):
        self.it = int(it)

    def add_cpu(self, cpu):
        self.cpu = float(cpu)

    def csv_header(self):
        return ("Name,Status,Rows b.p.,Cols b.p.,Rows a.p.,Cols a.p.,Nonzeros,Density,"
                "Relative Infeas,Relative Compl,Primal Objective,Dual Objetive,"
                "Max Add. Corr.,Iters,CPU Time [secs]\n")

    def csv(self):
        return "{},{},{},{},{},{},{},{:e},{:e},{:e},{:e},{:e},{},{},{:e}\n".format(self.name,
                self.status, self.brows, self.bcols, self.arows, self.acols,
                self.fact_nonzeros, self.fact_density, self.rel_infeas,
                self.rel_compl, self.prim_obj, self.dual_obj, self.max_corr, self.it, self.cpu)

def bench(S, p, o, mf, s, g, t, k):
    """
    Run the benchmark.

    :param S: Specification file to be used.
    :type S: str
    :param p: Path of the files to be used.
    :type p: str
    :param o: Output file name.
    :type o: str
    :param mf: Maximum number of files to be used.
    :type mf: int
    :param s: Sort method of files to be used.
    :type s: int
    :param g: Group of files to be tested.
    :type g: list
    :param t: time (in seconds) to keeping search for the solution of a file.
    :type t: int
    :param k: keep the log and stdout
    :type k: bool
    """
    import subprocess
    import fnmatch

    # Check if path exist.
    if not os.path.exists(p):
        raise FileNotFoundError("can't find directory '{0}'".format(p))

    of = open(o, 'w')
    of.write(PTInfo('').csv_header())

    # Process list of files to be used in the benchmark.
    if g[0] == 1:
        fl = g[1]
    else:
        fl = os.listdir(p)
    if g[0] == 2:
        for r in g[1]:
            try:
                fl.remove(r)
            except ValueError:
                pass
    elif g[0] == 3:
        fl = fnmatch.filter(fl, 'netlib-*')
    elif g[0] == 4:
        fl = fnmatch.filter(fl, 'netlibqap-*')
    elif g[0] == 5:
        fl = fnmatch.filter(fl, 'kennington-*')
    elif g[0] == 6:
        fl = fnmatch.filter(fl, 'meszaros-*')
    elif g[0] == 7:
        fl = fnmatch.filter(fl, 'pds-*')
    elif g[0] == 8:
        fl = fnmatch.filter(fl, 'rail-*')
    elif g[0] == 9:
        fl = fnmatch.filter(fl, 'fctp-*')

    # Sort the list of files to be used in the benchmark.
    if s == 1:
        fl.sort(key=lambda fname: os.path.getsize('{0}/{1}'.format(p,fname)),
                reverse=True)
    elif s == 2:
        fl.sort(key=lambda fname: os.path.getsize('{0}/{1}'.format(p,fname)))
    elif s == 3:
        fl.sort()

    # Limit the size of the list of files to be used in the benchmark.
    if mf:
        fl = fl[0:mf]

    # Run PCx for every file in ``fl``.
    for f in fl:
        print("Running PCx for {0}. It can take some minutes.".format(f))
        with open(f.replace(".mps", ".stdout"), 'w') as log:
            if S:
                PCx_instance = subprocess.Popen(["./PCx", "-s", S, f], stdout=log)
            else:
                PCx_instance = subprocess.Popen(["./PCx", f], stdout=log)
            try:
                retcod = PCx_instance.wait(t)
            except subprocess.TimeoutExpired:
                print("PCx is taking too long for find a solution to {0}."
                      " Stoping the search and killing PCx.".format(f))
                PCx_instance.kill()
                retcod = -1
        print("End running PCx.")
        info = PTInfo(f)
        if retcod in [0, 11, 12, 13, 14, 15]:
            with open(f.replace(".mps", ".log"), 'r') as log:
                step = 0  # step is responsible to select the regex to be used.
                for l in log:
                    if step == 0:
                        m = re.match('Iterations=(?P<it>[0-9]*), Termination Code=(?P<status>[0-9]*)', l)
                        if m:
                            info.add_status(m.group('status'))
                            info.add_it(m.group('it'))
                            step += 1
                    if step == 1:
                        m = re.match('.*Maximum Gondzio corrections = *(?P<cor>[0-9]*)', l)
                        if m:
                            info.add_max_corr(m.group('cor'))
                            step += 1
                    elif step == 2:
                        m = re.match('.*Before Presolving: *(?P<rows>[0-9]*) rows, *(?P<cols>[0-9]*) columns', l)
                        if m:
                            info.add_brows(m.group('rows'))
                            info.add_bcols(m.group('cols'))
                            step += 1
                    elif step == 3:
                        m = re.match('.*After  Presolving: *(?P<rows>[0-9]*) rows, *(?P<cols>[0-9]*) columns', l)
                        if m:
                            info.add_arows(m.group('rows'))
                            info.add_acols(m.group('cols'))
                            step += 1
                    elif step == 4:
                        m = re.match('.*Nonzeros in L=(?P<nonzeros>[0-9]*); *Density of L=(?P<density>-?[0-9].[0-9]*).*', l)
                        if m:
                            info.add_fact_nonzeros(m.group('nonzeros'))
                            info.add_fact_density(m.group('density'))
                            step += 1
                    elif step == 5:
                        m = re.match('Primal Objective = *(?P<po>-?[0-9]*.[0-9]*e[+-][0-9]*).*', l)
                        if m:
                            info.add_prim_obj(m.group('po'))
                            step += 1
                    elif step == 6:
                        m = re.match('Dual   Objective = *(?P<po>-?[0-9]*.[0-9]*e[+-][0-9]*).*', l)
                        if m:
                            info.add_dual_obj(m.group('po'))
                            step += 1
                    elif step == 7:
                        m = re.match('Relative Complementarity = *(?P<rc>-?[0-9]*.[0-9]*e[+-][0-9]*).*', l)
                        if m:
                            info.add_rel_compl(m.group('rc'))
                            step += 1
                    elif step == 8:
                        m = re.match('Primal = *(?P<rip>-?[0-9]*.[0-9]*e[+-][0-9]*).*', l)
                        if m:
                            info.add_rel_infeas(m.group('rip'))
                            step += 1
                    elif step == 9:
                        m = re.match('Time to solve.*: *(?P<cpu>[0-9]*.\.[0-9]*).*', l)
                        if m:
                            info.add_cpu(m.group('cpu'))
                            step += 1
        else:
            info.add_status(retcod)
        of.write(info.csv())
        # Remove log and stdout
        if not k:
            subprocess.call(["rm", "-f", f.replace(".mps", ".log")])
            subprocess.call(["rm", "-f", f.replace(".mps", ".stdout")])
    of.close()

def check_spc_file(f2check, p):
    """
    This function check if a given specification file can be used.

    :param f2check: Name of the file to be checked.
    :type f2check: str
    :param p: Path to be used when check the file.
    :type p: str
    """
    with open(f2check, 'r') as f:
        for l in f:
            m = re.match('history.*(?P<toggle>(yes|no|))', l)
            if m:
                if m.group('toggle') in ['yes', '']:
                    history = True
                elif m.group('toggle') == 'no':
                    history == False
            m = re.match('inputdirectory (?P<dir>.*)', l)
            if m:
                if m.group('dir') == p:
                    path = True
                else:
                    path = False
    return history and path

def check_spc(path):
    """
    This function check if the specifications file can be used.

    :param p: Path to be used when check the file.
    :type p: str
    """
    spc_correct = False
    if os.path.isfile('spc'):
        spc_correct = check_spc_file('spc', path)
    elif os.path.isfile('specs'):
        spc_correct = check_spc_file('sepcs', path)
    elif os.path.isfile('PCx.specs'):
        spc_correct = check_spc_file('PCx.specs', path)
    return spc_correct

if __name__ == "__main__":
    """Benchmark for PCx. ::

        $ bench.py --help
    """
    import sys
    import argparse
    from argparse import RawTextHelpFormatter

    if sys.version_info[0] >= 3:
        if sys.version_info[1] < 3:
            raise Exception('It require Python 3.3 or later.')
    else:
        raise Exception('It require Python 3.3 or later.')

    # Parse of flags.
    parser = argparse.ArgumentParser(description='Benchmark for PCx.',
            formatter_class=RawTextHelpFormatter)
    parser.add_argument('-p', '--path', type=str, default='./mps',
            help='Path of the files to be used. [default: `./mps`]')
    parser.add_argument('-o', '--output', type=str, default='bench.out',
            help='Output file name. [default: `bench.out`]')
    parser.add_argument('--max', type=int, default=None,
            help="Maximum number of files to be used. [default: None."
            " All files will be used]")
    parser.add_argument('-t', '--time', type=int, default=7200,
            help="The maximum time (in seconds) keeping search for the"
            " solution of a file. [default: 7200s]")
    parser.add_argument('-S', '--specs', type=str, default=None,
            help='Specification file to be used. [default: None]')
    parser.add_argument('-K', '--keep-all-log', action='store_true',
            help='Keep the log and stdout.')
    bench_sort = parser.add_mutually_exclusive_group()
    bench_sort.add_argument('-n', '--name', action='store_true',
            help='Sort list of files to be used by name.')
    bench_sort.add_argument('-d', '--decrescent-size', action='store_true',
            help='Sort list of files to be used by size (in decrescent order).')
    bench_sort.add_argument('-c', '--crescent-size', action='store_true',
            help='Sort list of files to be used by size (in crescent order).')
    bench_file = parser.add_mutually_exclusive_group()
    bench_file.add_argument('-i', '--input', nargs='+', type=str, default=[],
            help="List of files to be used for the benchmark.")
    bench_file.add_argument('-x', '--exclude', nargs='+', type=str, default=[],
            help="List of files to be exclude from the benchmark.")
    bench_file.add_argument('--netlib', action='store_true',
            help="Use only problem tests from NETLIB.")
    bench_file.add_argument('--netlibqap', action='store_true',
            help="Use only problem tests from NETLIB QAP.")
    bench_file.add_argument('--kenlib', action='store_true',
            help="Use only problem tests from KENNINGTON.")
    bench_file.add_argument('--meslib', action='store_true',
            help="Use only problem tests from MESZAROS.")
    bench_file.add_argument('--pdslib', action='store_true',
            help="Use only problem tests from PSD.")
    bench_file.add_argument('--raillib', action='store_true',
            help="Use only problem tests from RAIL.")
    bench_file.add_argument('--fctplib', action='store_true',
            help="Use only problem tests from FCTP.")

    args = parser.parse_args()

    # Reduce group bench_sort
    # 0 for no sort
    # 1 for descrescent size order
    # 2 for crescent size order
    # 3 for alphabetic order
    if args.decrescent_size:
        s = 1
    elif args.crescent_size:
        s = 2
    elif args.name:
        s = 3
    else:
        s = 0

    # Reduce group bench_file
    # 0 for none
    # 1 for include
    # 2 for exclude
    # 3 for netlib
    # 4 for netlibqap
    # 5 for kenlib
    # 6 for meslib
    # 7 for pdslib
    # 8 for raillib
    # 9 for fctplib
    if args.input:
        f = [1, args.input]
    elif args.exclude:
        f = [2, args.exclude]
    elif args.netlib:
        f = [3]
    elif args.netlibqap:
        f = [4]
    elif args.kenlib:
        f = [5]
    elif args.meslib:
        f = [6]
    elif args.pdslib:
        f = [7]
    elif args.raillib:
        f = [8]
    elif args.fctplib:
        f = [9]
    else:
        f = [0]

    if check_spc(args.path):
        bench(args.specs, args.path, args.output, args.max, s,
                f, args.time, args.keep_all_log)
    else:
        print("The specification file must contain 'history yes'.\n")
