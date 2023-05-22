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
endif


DICTB           =ReadRunDictUX
DICTH           =${DICTB}.h
DICT            =${DICTB}.cc
DICTO           =${DICTB}.o

HDRS            =src/ReadRun.h src/LinkDef.h
DICTHDRS        =$(HDRS) 
OBJS            =src/ReadRun.o $(DICTO)

LIBSLIN			=$(shell root-config --glibs)

SLL             =ReadRunLib.sl

#__________________________________________________________

all:		${OBJS} ${DICTO}
		${LD} -shared ${CXXFLAGS} -o ${SLL} ${OBJS} ${LIBSLIN}	
		$(MAKE) clean-intermediate

${DICT}:        ${DICTHDRS}
		rootcint -f ${DICT} -c ${DICTHDRS}

clean-intermediate:
		@rm ${OBJS} ${DICTB}.cc

clean:
		@rm ${SLL}

%.o             : %.cc
		${LD} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@          





 








