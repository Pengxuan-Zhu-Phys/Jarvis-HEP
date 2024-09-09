

LIB = ../nngg.a



VPATH := $(VPATH) ./kin ./diag ./univariate ./system

OBJS = \
  logfile.o \
  VecSet.o


$(LIB): $(LIB)($(OBJS))

$(LIB)(Gauss.o Patterson.o): univariate/config.h

$(LIB)(Eigensystem.o SingularValues.o PseudoEigensystem.o \
  ApplyEigen.o ApplySVD.o): diag/config.h



clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) $(LIB)

