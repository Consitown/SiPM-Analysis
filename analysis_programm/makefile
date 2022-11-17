CXX             =g++
FC              =g++
LD              =g++
MKDEP           =makedepend

RKDIR           =$(shell pwd)
OBJDIR          =$(WRKDIR)/obj
INCDIR          =$(WRKDIR)/inc
SRCDIR          =$(WRKDIR)/src
INCROOT         =$(shell root-config --incdir)
CFLAGS  = -Dextname -Df2cFortran -I/usr/local/include 
#-I/home/chscharf/root/include/

CPPFLAGS        =-I $(INCROOT)/ -I $(INCDIR)/
CXXFLAGS        =-fPIC -g -O --std=c++17
MKDEPFLAGS      =-Y ${INCL} -m -w 110

DICTB           =ReadRunDictUX
DICTH           =${DICTB}.h
DICT            =${DICTB}.cc
DICTO           =${DICTB}.o

HDRS            =ReadRun.h
DICTHDRS        = $(HDRS) LinkDef.h
SRCS            =ReadRun.cc
OBJS            =ReadRun.o $(DICTO)

ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs) 
LIBSLIN       = $(ROOTGLIBS) -lpthread -lm -ldl

SLL             =ReadRunLib.sl

#__________________________________________________________

all:		${OBJS} ${DICTO}
		${LD} -shared ${CXXFLAGS} -o ${SLL} ${OBJS} ${LIBSLIN}	

depend:         ${HDRS} ${SRCS}
		${MKDEP} ${MKDEPFLAGS} ${SRCS}

${DICT}:        ${DICTHDRS}
		rootcint -f ${DICT} -c ${DICTHDRS}

clean:
		@rm ${SL} ${OBJS} ${DICTO} ${DICTB}.* core

%.o             : %.cc
		${LD} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@          





 








