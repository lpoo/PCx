#!/usr/bin/python
import os
import re

class PTInfo():
    def __init__(self, name):
        self.name = name
        self.brows = 0
        self.bcols = 0
        self.arows = 0
        self.acols = 0
        self.rel_infeas = 0
        self.rel_compl = 0
        self.prim_obj = 0
        self.max_corr = 0
        self.it = 0
        self.cpu = 0

    def add_status(self, status):
        self.status = status

    def add_brows(self, rows):
        self.brows = rows

    def add_bcols(self, cols):
        self.bcols = cols

    def add_arows(self, rows):
        self.arows = rows

    def add_acols(self, cols):
        self.acols = cols

    def add_rel_infeas(self, rel_infeas):
        self.rel_infeas = rel_infeas

    def add_rel_compl(self, rel_compl):
        self.rel_compl = rel_compl

    def add_prim_obj(self, prim_obj):
        self.prim_obj = prim_obj

    def add_max_corr(self, max_corr):
        self.max_corr = max_corr

    def add_it(self, it):
        self.it = it

    def add_cpu(self, cpu):
        self.cpu = cpu

    def csv_header(self):
        return ("Name,Rows b.p.,Cols b.p.,Rows a.p.,Cols a.p.,Relative Infeas,"
                "Relative Compl,Primal Objective,Max Add. Corr.,Iters,CPU Time"
                " [secs]")

    def csv(self):
        return "{},{},{},{},{},{},{},{},{},{},{}".format(self.name, self.brows,
                self.bcols, self.arows, self.acols, self.rel_infeas,
                self.rel_compl, self.prim_obj, self.max_corr, self.it, self.cpu)

def bench(out_file, path, exclude):
    import subprocess

    of = open(out_file, 'w')
    of.write(PTInfo('').csv_header())

    for f in os.listdir(path):
        subprocess.call(["./PCx", f])
        with open(f.replace("mps", "log"), 'r') as log:
            info = PTInfo(f)
            step = 0  # step is responsible to select the regex to be used.
            for l in log:
                if step == 0:
                    m = re.match('Iterations=(?P<it>[0-9]*), Termination Code=(?P<status>[0-9])', l)
                    if m:
                        print(l)
                        info.add_status(m.group('status'))
                        info.add_it(m.group('it'))
                        step += 1
                elif step == 1:
                    m = re.match('.*Before Presolving:.* (?P<rows>[0-9]*) rows,.* (?P<cols>[0-9]*) columns', l)
                    if m:
                        info.add_brows(m.group('rows'))
                        info.add_bcols(m.group('cols'))
                        step += 1
                elif step == 2:
                    print(l)
                    m = re.match('.*After  Presolving:.* (?P<rows>[0-9]*) rows,.* (?P<cols>[0-9]*) columns', l)
                    if m:
                        info.add_arows(m.group('rows'))
                        info.add_acols(m.group('cols'))
                        step += 1
                elif step == 3:
                    m = re.match('Primal Objective = (?P<po>[ -][0-9]*.[0-9]*e[+-][0-9]*).*', l)
                    if m:
                        info.add_prim_obj(m.group('po'))
                        step += 1
                elif step == 4:
                    m = re.match('Relative Complementarity = (?P<rc>[ -][0-9]*.[0-9]*e[+-][0-9]*).*', l)
                    if m:
                        info.add_rel_compl(m.group('rc'))
                        step += 1
                elif step == 5:
                    m = re.match('Primal =(?P<rip>[ -][0-9]*.[0-9]*e[+-][0-9]*).*', l)
                    if m:
                        info.add_rel_infeas(m.group('rip'))
                        step += 1
                elif step == 6:
                    m = re.match('Time to solve.*:= (?P<cpu>[0-9]*.\.[0-9]*).*', l)
                    if m:
                        info.add_cpu(m.group('cpu'))
                        step += 1
            of.write(info.csv())
        subprocess.call(["rm", f.replace(".mps", "*")])
    of.close()

def check_spc_file(f2check):
    """This function check if a given specification file can be used."""
    f_correct = False
    with open(f2check, 'r') as f:
        for l in f:
            m = re.match('history.*(?P<toggle>(yes|no|))', l)
            if m:
                if m.group('toggle') in ['yes', '']:
                    f_correct = True
                elif m.group('toggle') == 'no':
                    f_correct == False
    return f_correct

def check_spc():
    """This function check if the specifications file can be used."""
    spc_correct = False
    if os.path.isfile('spc'):
        spc_correct = check_spc_file('spc')
    elif os.path.isfile('specs'):
        spc_correct = check_spc_file('sepcs')
    elif os.path.isfile('PCx.specs'):
        spc_correct = check_spc_file('PCx.specs')
    return spc_correct

if __name__ == "__main__":
    """Benchmark for PCx. ::

        $ bench.py --help
    """
    import sys
    import argparse
    from argparse import RawTextHelpFormatter

    if sys.version_info[0] < 3:
        raise Exception('It require Python 3.0 or later.')

    # Parse of flags.
    parser = argparse.ArgumentParser(description='Benchmark for PCx.',
            formatter_class=RawTextHelpFormatter)
    parser.add_argument('-o', type=str, default='bench.out',
            help='Output file.')
    parser.add_argument('-p', type=str, default='mps',
            help='Path of the files to test.')
    parser.add_argument('-x', nargs='+', type=str, default=[],
            help="List of files to be exclude from the benchmark.\n")

    args = parser.parse_args()

    if check_spc():
        bench(args.o, args.p, args.x)
    else:
        print("The specification file must contain 'history yes'.\n")
