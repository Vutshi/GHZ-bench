#!/usr/bin/env python3

import os
import sys
import csv
import argparse
from glob import glob

class Test:
    def __init__(self, cpu, frequency, compiler):
        self.cpu = cpu
        self.frequency = frequency
        self.compiler = compiler
        self.naive = 0.0
        self.unroll = 0.0
        self.unroll_soa = 0.0
        self.unroll_soa_split = 0.0

def load_tests(sort=False):
    logs = glob('logs/*.log')
    tests = []
    for path in logs:
        with open(path, 'r', encoding='utf-8') as f:
            components = os.path.basename(path)[:-4].split('_')
            cpu = components[0]
            freq = components[1]
            freq = '{0} {1}'.format(freq[:-3], freq[-3:])
            compiler = components[2]
            test = Test(cpu, freq, compiler)
            while True:
                line = ''
                while True:
                    line = f.readline()
                    if line == '' or line.startswith('#'):
                        break
                if line == '':
                    break
                testid = line.split(' ')[0]
                f.readline()
                gflops = float(f.readline().split('=')[1].strip())
                f.readline()
                if testid == '#1':
                    test.naive = gflops
                elif testid == '#2':
                    test.unroll = gflops
                elif testid == '#3':
                    test.unroll_soa = gflops
                elif testid == '#4':
                    test.unroll_soa_split = gflops
            tests.append(test)
    if sort:
        tests.sort(
            key=lambda test: max(test.naive, test.unroll, test.unroll_soa, test.unroll_soa_split),
            reverse=True
        )
    return tests

def print_csv(tests, outfile):
    fieldnames = ['CPU', 'Frequency', 'Compiler', 'Naive', 'Unroll', 'Unroll+SoA', 'Unroll+SoA+split']
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()
    for test in tests:
        tmp = dict()
        tmp['CPU'] = test.cpu
        tmp['Frequency'] = test.frequency
        tmp['Compiler'] = test.compiler
        tmp['Naive'] = '%.2f' % test.naive
        tmp['Unroll'] = '%.2f' % test.unroll
        tmp['Unroll+SoA'] = '%.2f' % test.unroll_soa
        tmp['Unroll+SoA+split'] = '%.2f' % test.unroll_soa_split
        writer.writerow(tmp)

def print_md(test, outfile):
    def p(*args, **kwargs):
        print(*args, **kwargs, file=outfile)

    p('# Single thread, GFLOPS')
    p()

    p('| {0:<18} '.format('CPU'), end='')
    p('| {0:>9} '.format('Frequency'), end='')
    p('| {0:<18} '.format('Compiler'), end='')
    p('| {0:>8} '.format('Naive'), end='')
    p('| {0:>8} '.format('Unroll'), end='')
    p('| {0:>10} '.format('Unroll+SoA'), end='')
    p('| {0:>16} |'.format('Unroll+SoA+split'))

    p('|-{0:-<18}-'.format('-'), end='')
    p('|-{0:->9}:'.format('-'), end='')
    p('|-{0:-<18}-'.format('-'), end='')
    p('|-{0:-<8}:'.format('-'), end='')
    p('|-{0:-<8}:'.format('-'), end='')
    p('|-{0:-<10}:'.format('-'), end='')
    p('|-{0:-<16}:|'.format('-'))

    for test in tests:
        p('| {0:<18} '.format(test.cpu), end='')
        p('| {0:>9} '.format(test.frequency), end='')
        p('| {0:<18} '.format(test.compiler), end='')
        p('| {0:>8} '.format('%.2f' % test.naive), end='')
        p('| {0:>8} '.format('%.2f' % test.unroll), end='')
        p('| {0:>10} '.format('%.2f' % test.unroll_soa), end='')
        p('| {0:>16} |'.format('%.2f' % test.unroll_soa_split))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process logs and generate table')
    parser.add_argument('file', metavar='FILE', type=str, nargs='?',
                        help='output to file (default: stdout)')
    parser.add_argument('--csv', dest='print', action='store_const',
                        const=print_csv, default=print_md,
                        help='Output CSV (default: Markdown)')
    parser.add_argument('--no-sort', dest='sort', action='store_false',
                        help='Sort by GFLOPS')
    args = parser.parse_args()
    tests = load_tests(args.sort)
    if args.file:
        with open(args.file, 'w', newline='') as outfile:
            args.print(tests, outfile)
    else:
        args.print(tests, sys.stdout)
