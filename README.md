# GHZ-bench

This is a simple benchmark created while studying different ways of optimization for x86 and e2k architechtures. It features the main loop of a program evaluating noise properties of a quantum interferometer employing a cascade of GHZ states.

Results for various CPUs are shown [here](results.md).

## Build

Build
`make`

Run
`./ghz-bench`

Alternatively, one can compile as follows

```clang -Ofast -march=native ghz-bench.c -o ghz-bench```. 

In some cases, it might be beneficial to point out the target architechture explicetly, e.g., `-march=icelake-server`. It is known that `gcc` does not generate efficent binary with this code, so `clang` is recommended on x86.

## About

This program was created by [Vutshi](https://github.com/Vutshi) and a friend.
