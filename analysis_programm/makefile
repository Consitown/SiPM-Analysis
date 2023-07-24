CXX             =g++
LD              =g++

INCROOT         =$(shell root-config --incdir)

CPPFLAGS        =-I $(INCROOT)/ -I $(INCDIR)/
CXXFLAGS        =-fPIC -g -O2 -Wall -Wextra

# Check the version of the C++ compiler for older systems
CXXVER := $(shell $(CXX) -dumpversion)
# $(info g++ compiler version $(CXXVER))
ifeq ($(shell expr $(word 1,$(subst ., ,$(CXXVER))) \< 7), 1)
	ifeq ($(shell expr $(word 1,$(subst ., ,$(CXXVER))) \< 5), 1)
		CXXFLAGS += -std=c++11
	else
		CXXFLAGS += -std=c++14
	endif
else 
	CXXFLAGS += -std=c++17
endif


DICTB           =ReadRunDictUX
DICTH           =${DICTB}.h
DICT            =${DICTB}.cc
DICTO           =${DICTB}.o

HDRS            =src/PMT.h src/CosmicsBox.h src/FFT_WF.h src/Experimental.h src/ReadRun.h 
DICTHDRS        =$(HDRS) 
OBJS            =src/PMT.o src/CosmicsBox.o src/FFT_WF.o src/Experimental.o src/ReadRun.o $(DICTO)

LIBSLIN			=$(shell root-config --glibs)

SLL             =ReadRunLib.sl

#__________________________________________________________

all:		${OBJS} ${DICTO}
		${LD} -shared ${CXXFLAGS} -o ${SLL} ${OBJS} ${LIBSLIN}	
		$(MAKE) clean-intermediate

${DICT}:        ${DICTHDRS}
		rootcint -rml=ReadRun -f ${DICT} -c ${DICTHDRS} -I. misc/LinkDef.h

clean-intermediate:
		@rm ${OBJS} ${DICTB}.cc

clean:
		@rm ${SLL}

%.o             : %.cc
		${LD} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@          





 








