LIB = ../nngg.a

OBJS := \
  abbr1s1.o \
  abbr1s2.o \
  abbr1s3.o \
  abbr1a.o \
  abbr1h1.o \
  abbr1h2.o \
  abbr1h3.o \
  abbr1h4.o \
  abbr1h5.o \
  abbr1h6.o \
  vert_cha5_cha6_hia5.o \
  vert_hia5_lpd5.o \
  vert_hia5_qud5.o \
  vert_hia5_quu5.o \
  vert_cha5_cha6_his5.o \
  vert_his5_lpd5.o \
  vert_his5_qud5.o \
  vert_his5_quu5.o \
  vert_his5_sld5.o \
  vert_his5_sle5.o \
  vert_his5_sqd5.o \
  vert_his5_sqe5.o \
  vert_his5_squ5.o \
  vert_his5_sqv5.o \
  vert_his5.o \
  vert_cha5_cha6.o \
  vert_lpd5.o \
  vert_qud5.o \
  vert_quu5.o \
  vert_sld5.o \
  vert_sle5.o \
  vert_sqd5.o \
  vert_sqe5.o \
  vert_squ5.o \
  vert_sqv5.o \
  vert.o \
  vert_cha5_neu5.o \
  vert_lpd5_neu5.o \
  vert_neu5_qud5.o \
  vert_neu5_quu5.o \
  vert_neu5.o \
  box_cha5_cha6.o \
  box_cha5.o \
  box_lpd5.o \
  box_qud5.o \
  box_quu5.o \
  gzSquaredME.o

$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): vars.h renconst.h model.h util.h looptools.h

clean distclean:
	$(RM) $(OBJS) $(LIB)

