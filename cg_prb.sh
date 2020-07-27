#!/bin/bash
#
gfortran -c cg_prb.f
if [ $? -ne 0 ]; then
  echo "Errors compiling cg_prb.f"
  exit
fi
#
gfortran cg_prb.o -L$HOME/libf77 -lcg
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cg_prb.o"
  exit
fi
rm cg_prb.o
#
mv a.out cg_prb
./cg_prb > cg_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cg_prb"
  exit
fi
rm cg_prb
#
echo "Test program output written to cg_prb_output.txt."
