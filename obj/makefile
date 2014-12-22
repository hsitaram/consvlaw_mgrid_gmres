FILES=../gmres.o
OBJFILES = gmres.o
CXXFLAGS = -g

FORTFILES = ../gmres.f90
FORTPRGFILE = ../testgmres.f90
FILECONVDIFF = ../convdiffmain.f90
CONVDIFFOBJ = convdiff.o

FORT = gfortran

fortranconvdiff:
	${FORT} -c -fbounds-check ${FORTFILES} ../convdiff.f90
	${FORT} -fbounds-check ${OBJFILES} ${CONVDIFFOBJ} ${FILECONVDIFF} -o f90convdiff.exe

clean:
	rm -rf *.o *.mod *.exe *.dat
