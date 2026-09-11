CXX=g++
LINKER=g++
EXEC=analysis
EXEC_FRS_PID=plot_frs_pid_cocktailbeam
EXEC_GLAD_PID=plot_glad_pid

ROOTCONFIG := root-config

CFLAGS := $(shell $(ROOTCONFIG) --cflags)
CFLAGS += -I${FAIRROOTPATH}/include
CFLAGS += -I$(SIMPATH)/include
CFLAGS += -I$(ROOT_INCLUDE_PATH)
CFLAGS += -I$(ROOT_INCLUDE_DIR)
CFLAGS += -I$(VMCWORKDIR)
CFLAGS += -I$(VMCWORKDIR)/r3bdata/footData
CFLAGS += -I$(VMCWORKDIR)/r3bdata/wrData
CFLAGS += -I$(VMCWORKDIR)/r3bdata/frsData
CFLAGS += -I$(VMCWORKDIR)/r3bdata/tofData
CFLAGS += -I$(VMCWORKDIR)/r3bdata/califaData
CFLAGS += -I$(VMCWORKDIR)/r3bdata/neulandData
CFLAGS += -I$(VMCWORKDIR)/r3bdata
CFLAGS += -I$(VMCWORKDIR)/r3bsource/base
CFLAGS += -I$(VMCWORKDIR)/r3bbase
CFLAGS += -I$(VMCWORKDIR)/tracking

#CFLAGS += -I$(UCESB_DIR)/hbook

CFLAGS += --std=c++17 -g -O0 -fexceptions

#CFLAGS += -g -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings 

LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
LDFLAGS += -lEG $(shell $(ROOTCONFIG) --glibs)
LDFLAGS += -L$(ROOT_LIBRARY_DIR) -L$(FAIRROOTPATH)/lib
LDFLAGS += -g

#LDFLAGS += -L$(VMCWORKDIR)/../build/lib -lR3BSsd -lR3Bsource -lR3BBase

LDFLAGS += -L$(VMCWORKDIR)/../build/lib -lR3BSsd -lR3BBase -lR3BData -lR3BTracking
LDFLAGS += -L$(FAIRROOTPATH)/lib -lBase -lParBase -lFairTools
LDFLAGS += -lR3BData

INCLUDEDIR=include

DIR_INC=-I$(INCLUDEDIR)

SRC=elastic_proton_scattering.C
SRC_FRS_PID=plot_frs_pid_cocktailbeam.C
SRC_GLAD_PID=plot_glad_pid.C
#SRC=delta_production_mul2.C

OBJ=$(SRC:.C=.o)
OBJ:=$(OBJ:.cxx=.o)
OBJ_FRS_PID=$(SRC_FRS_PID:.C=.o)
OBJ_GLAD_PID=$(SRC_GLAD_PID:.C=.o)

all: $(EXEC) $(EXEC_FRS_PID) $(EXEC_GLAD_PID)

$(EXEC): $(OBJ)
	${LINKER} -o $@ $^ ${LDFLAGS}
	echo " COMP $@"
%.o : %.C
	$(MAKEDEPEND)
	${CXX} ${CFLAGS} $(DIR_INC) -c $< -o $@
	echo "	CXX $@"

$(EXEC_FRS_PID): $(OBJ_FRS_PID)
	${LINKER} -o $@ $^ ${LDFLAGS}
	@echo " COMP $@"

$(EXEC_GLAD_PID): $(OBJ_GLAD_PID)
	${LINKER} -o $@ $^ ${LDFLAGS}
	@echo " COMP $@"

%.o : %.cpp
	${CXX} ${CFLAGS} $(DIR_INC) -c $< -o $@
	@echo "   CXX $@"

clean:
	rm -f *.o $(EXEC) $(EXEC_FRS_PID) $(EXEC_GLAD_PID)
