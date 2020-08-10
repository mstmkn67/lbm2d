CC = g++
CFLAG = -c -O3 -I $(PF_FILES)/include
CCFLAG = -o lbm2d  -L $(PF_FILES)/lib/linux_64 -lplatform
OBJECT= Timer.o gourmain.o cmdlineopts.o \
	LBM2dSimulator.o LBM2dPlates.o LBM2dLeesEdwards.o LBM2dLattice.o  \
	LBM2d.o lbm2d_driver.o

all :$(OBJECT)
	$(CC) $(OBJECT) $(CCFLAG)

Timer.o : Timer.cpp Timer.h
	$(CC) Timer.cpp $(CFLAG)

gourmain.o : udf/gourmain.cpp udf/gourmain.h
	$(CC) udf/gourmain.cpp $(CFLAG)

cmdlineopts.o : udf/cmdlineopts.cpp udf/cmdlineopts.h
	$(CC) udf/cmdlineopts.cpp $(CFLAG)

LBM2dSimulator.o : lbm/LBM2dSimulator.cpp lbm/LBM2dSimulator.h
	$(CC) lbm/LBM2dSimulator.cpp $(CFLAG)

LBM2dPlates.o : lbm/LBM2dPlates.cpp lbm/LBM2dPlates.h
	$(CC) lbm/LBM2dPlates.cpp $(CFLAG)

LBM2dLeesEdwards.o : lbm/LBM2dLeesEdwards.cpp lbm/LBM2dLeesEdwards.h
	$(CC) lbm/LBM2dLeesEdwards.cpp $(CFLAG)

LBM2dLattice.o : lbm/LBM2dLattice.cpp lbm/LBM2dLattice.h
	$(CC) lbm/LBM2dLattice.cpp $(CFLAG)

LBM2d.o : LBM2d.cpp LBM2d.h
	$(CC) LBM2d.cpp $(CFLAG)

lbm2d_driver.o : lbm2d_driver.cpp
	$(CC) lbm2d_driver.cpp $(CFLAG)

clean : 
	rm -f *.o lbm2d

install :
	strip lbm2d
	cp lbm2d ../bin
