#! /bin/bash
# test/testMass.sh.  Generated from testMass.sh.in by configure.

# want to see this output
./testMass

OS=`uname`
case "$OS" in
CYGWIN*)
  cmd1=`sed 's/e-0/e-/g' testMass1.out | \
        sed 's/e+0/e+/g' | sed 's/MEV/GEV/g' | \
	sed 's/GEV/GEV/g' | \
        sed 's/MM/MM/g' | \
         diff -q -b - ./testMass1.dat`
  cmd2=`diff -q -b testMass2.out testMass1.out`
  ;;
*)
  cmd1=`sed 's/GEV/GEV/g' testMass1.out | \
        sed 's/MM/MM/g' | \
        diff -q -b - ./testMass1.dat`
  cmd2=`diff -q -b testMass2.out testMass1.out`
esac

if [ -n "$cmd1" ]
then
  echo $cmd1
  exit 1;
fi

if [ -n "$cmd2" ]
then
  echo $cmd2
  exit 1;
fi

exit 0;

