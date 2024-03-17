#! /bin/bash
# test/testHepMCIteration.sh.  Generated from testHepMCIteration.sh.in by configure.

rm -f testHepMCIteration*out

./testHepMCIteration >& testHepMCIteration.cout

OS=`uname`
case "$OS" in
CYGWIN*)
  #don't compare these on Windows (hopelessly different default output)
  cmd1=
  cmd3=`diff -q -b testHepMCIteration2.out testHepMCIteration.out`
  cmd4=`diff -q -b testHepMCIteration3.out testHepMCIteration.out`
  ;;
Darwin*)
  # MacOSX fix
  cmd1=`sed 's/(-0,0)/(0,0)/g'  testHepMCIteration.out  | \
  	diff -q -b - ./testHepMCIteration.dat`
  cmd3=`diff -q -b testHepMCIteration2.out testHepMCIteration.out`
  cmd4=`diff -q -b testHepMCIteration3.out testHepMCIteration.out`
  ;;
*)
  cmd1=`sed 's/(-0,0)/(0,0)/g'  testHepMCIteration.out  | \
  	diff -q -b - ./testHepMCIteration.dat`
  cmd3=`diff -q -b testHepMCIteration2.out testHepMCIteration.out`
  cmd4=`diff -q -b testHepMCIteration3.out testHepMCIteration.out`
esac

if [ -n "$cmd1" ]
then
  echo $cmd1
  exit 1;
fi

if [ -n "$cmd3" ]
then
  echo $cmd3
  exit 1;
fi

if [ -n "$cmd4" ]
then
  echo $cmd4
  exit 1;
fi

exit 0;
