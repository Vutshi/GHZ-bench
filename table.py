#!/usr/bin/env python3

import os
import sys
import csv
import argparse
from glob import glob

def load_tests():
    logs = glob('logs/*.log')
    logs.sort()
    tests = []
    for path in logs:
        with open(path, 'r', encoding='utf-8') as f:
            test = dict()
            components = os.path.basename(path)[:-4].split('_')
            test['CPU'] = components[0]
            freq = components[1]
            test['Frequency'] = '{0} {1}'.format(freq[:-3], freq[-3:])
            test['Compiler'] = components[2]
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
                gflops = '%.2f' % float(f.readline().split('=')[1].strip())
                f.readline()
                if testid == '#1':
                    test['Naive'] = gflops
                elif testid == '#2':
                    test['Unroll'] = gflops
                elif testid == '#3':
                    test['Unroll+SoA'] = gflops
                elif testid == '#4':
                    test['Unroll+SoA+split'] = gflops
            tests.append(test)
    return tests

def print_csv(tests, outfile):
    fieldnames = ['CPU', 'Frequency', 'Compiler', 'Naive', 'Unroll', 'Unroll+SoA', 'Unroll+SoA+split']
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()
    for test in tests:
        writer.writerow(test)

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
        p('| {0:<18} '.format(test['CPU']), end='')
        p('| {0:>9} '.format(test['Frequency']), end='')
        p('| {0:<18} '.format(test['Compiler']), end='')
        p('| {0:>8} '.format(test['Naive']), end='')
        p('| {0:>8} '.format(test['Unroll']), end='')
        p('| {0:>10} '.format(test['Unroll+SoA']), end='')
        p('| {0:>16} |'.format(test['Unroll+SoA+split']))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process logs and generate table')
    parser.add_argument('file', metavar='FILE', type=str, nargs='?',
                        help='output to file (default: stdout)')
    parser.add_argument('--csv', dest='print', action='store_const',
                        const=print_csv, default=print_md,
                        help='Output CSV (default: Markdown)')
    args = parser.parse_args()
    tests = load_tests()
    if args.file:
        with open(args.file, 'w', newline='') as outfile:
            args.print(tests, outfile)
    else:
        args.print(tests, sys.stdout)
