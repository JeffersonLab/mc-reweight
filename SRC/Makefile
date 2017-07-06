my_objs = elastic.o gauss1.o mt19937.o f2allm.o r1990.o rcint.o rcext.o \
          rc_mod.o csb_jan05.o pos.o yp_optcor.o samp_eloss.o ran3.o \
          model_new.o f1f209.o 

INC= hms.cmn

RM        = rm -f 
SHELL     = /bin/sh

MYOS := $(subst -,,$(shell uname))

F77=g77
FC=gfortran

ifeq ($(MYOS),SunOS)
  FFLAGS= -e -O
  CERN_ROOT= /apps/cernlib/sun4_solaris2/95a
  CERNLIBS= -L$(CERN_ROOT)/lib -lpacklib
endif

ifeq ($(MYOS),HPUX)
  FFLAGS= +U77 +ppu -C +es -O +Onolimit +FPVZOU
  CERNLIBS= -Wl,-L$(CERN_ROOT)/lib -lpacklib
endif

ifeq ($(MYOS),ULTRIX)
  FFLAGS=-check_bounds
  CERNLIBS= -Wl,-L$(CERN_ROOT)/lib -lpacklib
endif

ifeq ($(MYOS),OSF1)
  FFLAGS=-O -C -extend_source -fpe
  CERNLIBS= -L$(CERN_ROOT)/lib -lpacklib
endif


ifeq ($(MYOS),Linux)
  FFLAGS=-O -C  -fno-range-check 
#  CERN_ROOT = /u/site/cernlib/i386_rhel4/2005
#  CERN_ROOT = /u/site/cernlib/x86_64_rhel7/2005
#  CERN_ROOT = /opt/cernlib/2005
  CERNLIBS= -L $(CERN_ROOT)/lib -lmathlib -lpacklib
endif
none: recon_mc

all: recon_mc

recon_mc: $(my_objs) recon_mc.o Makefile
	$(FC) -g -o  $@ $(FFLAGS) recon_mc.o $(my_objs) $(CERNLIBS)

clean: 
	$(RM) *.o recon_mc















