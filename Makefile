#  Declaration
CXXFLAGS = -Wall -v -g -O3
PARFLAGS = -fopenmp -g
CC	=	gcc
C++	=	g++
F77    	=	f77

.SUFFIXES:

.SUFFIXES: .f .F .cpp .c .o .h

%.o: %.f
	${F77} ${FFLAGS} $< -c
	
%.o: %.F
	${F77} ${FFLAGS} $< -c

%.o: %.cpp %.h
	${C++} ${CXXFLAGS} $< -c

%.o: %.c %.h
	${CC} ${CFLAGS} $< -c


objects	= App_Damage.o App_Fracture.o Fem_3D.o Gauss.o Geometry_3D.o \
	  Geometry_2D.o Global_Stiff_Matrix.o Global_Load_Vector.o Hns.o \
	  Input_Reader.o MatBase.o MathMatrix.o MatPro.o Mesher.o \
	  Postprocessor.o SolveEqu.o MainProg.o WeightFunc.o\

neca : $(objects)        
	${C++} ${PARFLAGS} -o neca $(objects)
	              
MainProg.o : MainProg.cpp
	${C++} ${CXXFLAGS} MainProg.cpp -c

Global_Stiff_Matrix.o : Global_Stiff_Matrix.cpp
	${C++} ${PARFLAGS} ${CXXFLAGS} $< -c 

SolveEqu.o : SolveEqu.cpp
	${C++} ${PARFLAGS} ${CXXFLAGS} $< -c
		
# 伪目标文件
.PHONY : clean
clean :
	-rm neca $(objects)
