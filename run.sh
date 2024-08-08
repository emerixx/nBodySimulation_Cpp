#!/bin/bash
g++ -c -o vctr.o vctr.cpp
g++ -c -o main.o main.cpp
g++ -o program.out main.o vctr.o -lraylib
./program.out
