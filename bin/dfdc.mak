#===========================================
# Makefile for DFDC ducted fan code
# Version 070-ES3 with ESLOFT and MODI
#===========================================

PLTLIB = ../plotlib/libPlt.a
#PLTLIB = ../plotlib/libPltDP.a

SRC1 = ../src-lib
SRC2 = ../src-driver

##--------------------------
## mingw  gfortran
BDIR = g:\gcc452\bin
FC = $(BDIR)\gfortran.exe
CC  = $(BDIR)\gfortran.exe 
DP =
FFLAGS  = -O -pipe $(DP)
CFLAGS  = -O -pipe
AR = $(BDIR)\ar.exe r
RANLIB = $(BDIR)\ranlib.exe 
LINKLIB = -lgdi32
LFLAGS = -Wl,--strip-debug  -static-libgfortran -static-libgcc 
LIBS =
##--------------------------
OBJSLIB = \
adjpanl.o aero.o \
aic.o airio.o \
atmo.o blcalc.o \
dfdcsubs.o forces.o \
gauss.o geom.o geutil.o \
grdutils.o  inigrd.o \
lamp.o plutil.o pnsubs.o \
qaic.o rotoper.o \
sgutil.o sgutil2.o \
solve.o spline.o system.o \
userio.o \
vels.o viscvel.o \
wakesubs.o xutils.o

OBJ2 = debug.o \
dfdc.o dfdc2.o \
gdes.o gplots.o \
gpoint.o gui.o \
list.o modify.o \
oper.o oplots.o oplset.o \
qdes.o qplots.o \
rplots.o sl.o \
esloft.o lplots.o xfsubs.o \
xmodify.o modi.o

PROGS = dfdc

#===========================================================

all:	 $(PROGS)

install: 
	$(INSTALLCMD) $(PROGS) $(BINDIR)

clean:
	-/bin/rm *.o
	-/bin/rm $(PROGS)
	-/bin/rm plot.ps

cleanall:
	-/bin/rm *.o
	-/bin/rm *.a
	-/bin/rm $(PROGS)
	-/bin/rm plot.ps

#===========================================================

dfdc: $(OBJSLIB)  $(OBJ2) 
	$(FC) -v  -o dfdc.exe  $(LFLAGS) $(OBJ2) $(OBJSLIB) $(PLTLIB) $(LINKLIB) $(FTNLIB)

#dfdclib.a: $(OBJSLIB) 
#	$(AR) -r dfdclib.a $(OBJSLIB) 

#===========================================================

$(SRC1)/DFDC.INC: 


#=========================================================================

debug.o: $(SRC2)/debug.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/debug.f

dfdc.o: $(SRC2)/dfdc.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/dfdc.f

dfdc2.o: $(SRC2)/dfdc2.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/dfdc2.f

gdes.o: $(SRC2)/gdes.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC $(SRC1)/MASKS.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/gdes.f

gplots.o: $(SRC2)/gplots.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC $(SRC1)/MASKS.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/gplots.f

gpoint.o: $(SRC2)/gpoint.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/gpoint.f

gui.o: $(SRC2)/gui.f
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/gui.f

list.o: $(SRC2)/list.f $(SRC1)/DFDC.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/list.f

modify.o: $(SRC2)/modify.f
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/modify.f

oper.o: $(SRC2)/oper.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/oper.f

oplots.o: $(SRC2)/oplots.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC $(SRC1)/MASKS.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/oplots.f

oplset.o: $(SRC2)/oplset.f $(SRC1)/PLOT.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/oplset.f

qdes.o: $(SRC2)/qdes.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC $(SRC1)/MASKS.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/qdes.f

qplots.o: $(SRC2)/qplots.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC $(SRC1)/MASKS.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/qplots.f

rplots.o: $(SRC2)/rplots.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/rplots.f

sl.o: $(SRC2)/sl.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/sl.f

esloft.o: $(SRC2)/esloft.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/esloft.f

lplots.o: $(SRC2)/lplots.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/lplots.f

xfsubs.o: $(SRC2)/xfsubs.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/xfsubs.f

xmodify.o: $(SRC2)/xmodify.f
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/xmodify.f

modi.o: $(SRC2)/modi.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/modi.f
#=========================================================================

adjpanl.o: $(SRC1)/adjpanl.f $(SRC1)/DFDC.INC 
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/adjpanl.f

aero.o: $(SRC1)/aero.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/aero.f

aic.o: $(SRC1)/aic.f
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/aic.f

airio.o: $(SRC1)/airio.f
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/airio.f

atmo.o: $(SRC1)/atmo.f 
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/atmo.f

blcalc.o: $(SRC1)/blcalc.f  $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/blcalc.f

dfdcsubs.o: $(SRC1)/dfdcsubs.f $(SRC1)/DFDC.INC $(SRC1)/PLOT.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/dfdcsubs.f

forces.o: $(SRC1)/forces.f $(SRC1)/DFDC.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/forces.f

gauss.o: $(SRC1)/gauss.f
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/gauss.f

geom.o: $(SRC1)/geom.f $(SRC1)/DFDC.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/geom.f

geutil.o: $(SRC1)/geutil.f
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/geutil.f

grdutils.o: $(SRC1)/grdutils.f
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/grdutils.f

inigrd.o: $(SRC1)/inigrd.f $(SRC1)/DFDC.INC $(SRC1)/MASKS.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/inigrd.f

lamp.o: $(SRC1)/lamp.f
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/lamp.f

pane.o: $(SRC1)/pane.f
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/pane.f

plutil.o: $(SRC1)/plutil.f $(SRC1)/PLOT.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/plutil.f

pnsubs.o: $(SRC1)/pnsubs.f $(SRC1)/DFDC.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/pnsubs.f

qaic.o: $(SRC1)/qaic.f $(SRC1)/DFDC.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/qaic.f

rotoper.o: $(SRC1)/rotoper.f $(SRC1)/DFDC.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/rotoper.f

sgutil.o: $(SRC1)/sgutil.f $(SRC1)/PLOT.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/sgutil.f

sgutil2.o: $(SRC1)/sgutil2.f
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/sgutil2.f

solve.o: $(SRC1)/solve.f $(SRC1)/DFDC.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/solve.f

spline.o: $(SRC1)/spline.f
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/spline.f

system.o: $(SRC1)/system.f $(SRC1)/DFDC.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/system.f

userio.o: $(SRC1)/userio.f
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/userio.f

vels.o: $(SRC1)/vels.f $(SRC1)/DFDC.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/vels.f

viscvel.o: $(SRC1)/viscvel.f $(SRC1)/DFDC.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/viscvel.f

wakesubs.o: $(SRC1)/wakesubs.f $(SRC1)/DFDC.INC
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/wakesubs.f

xutils.o: $(SRC1)/xutils.f
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/xutils.f

#===========================================================
