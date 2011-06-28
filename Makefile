# Build type
#   private:
#     Cephes library from canonical source.  Original cpoly and rpoly Fortran
#     functions from TOMS under ACM Software License Agreement, which is for
#     non-commercial use only.
#
#   gpl:
#     Based on GPL software: Cephes from LabPlot distribution under GPL, cpoly
#     from R distribution under GPL.  However, the result may not be
#     distributed under the GPL because the q library is neither compatible,
#     nor a system library (it's not normally distributed with KDB+).
#
#   bsd:
#     Cephes from Netlib, which has a note allowing free use.  poly implemented
#     using LAPACK.  Result is distributable under a BSD-style license.

#BUILD=private
#BUILD=gpl
BUILD=bsd


# Platform
PLATFORM=w32
#PLATFORM=w64


# Build tools
ifeq "$(PLATFORM)" "w64"
    MSSDK_BIN=$$MSSDK/Bin/win64/x86/AMD64
else
    MSSDK_BIN=$$MSSDK/vc/bin
endif
CL="$(MSSDK_BIN)/cl"
LINK="$(MSSDK_BIN)/link"


# Compilation flags
CFLAGS=-Ox
DEFINES=-D__STDC__ -D_IEEE_LIBM -D__LITTLE_ENDIAN -DWIN32


# Dependecy locations
SOURCEFORGE_DIR=http://prdownloads.sourceforge.net
#   FDLIBM
FDLIBM_SRC=$(SOURCEFORGE_DIR)/gnuwin32/fdlibm-5.2-2-src.zip
#   Cephes
CEPHES_SRC=http://www.moshier.net/double.zip
CPROB_SRC=http://www.netlib.org/cephes/cprob.tgz
LABPLOT_SRC=$(SOURCEFORGE_DIR)/labplot/LabPlot-1.6.0.2.tar.gz
#   CLAPACK
CLAPACK_DIR=http://www.netlib.org/clapack
CLAPACK_H=$(CLAPACK_DIR)/clapack.h
ifeq "$(PLATFORM)" "w64"
    CLAPACK_LIB_SUBDIR=LIB_WINDOWS/x64
else
    CLAPACK_LIB_SUBDIR=LIB_WINDOWS/Win32
endif
F2C_LIB=$(CLAPACK_DIR)/$(CLAPACK_LIB_SUBDIR)/libf2c.lib
BLAS_LIB=$(CLAPACK_DIR)/$(CLAPACK_LIB_SUBDIR)/BLAS.lib
CLAPACK_LIB=$(CLAPACK_DIR)/$(CLAPACK_LIB_SUBDIR)/clapack.lib
#   f2c program
F2C_BIN_GZ=http://www.netlib.org/f2c/mswin/f2c.exe.gz
F2C_H=http://www.netlib.org/f2c/f2c.h
#   cpoly and rpoly
CPOLY_F=http://www.netlib.org/toms/419
RPOLY_F=http://www.netlib.org/toms/493
CPOLY_C=http://svn.r-project.org/R/tags/R-2-9-2/src/appl/cpoly.c
#   q stuff
Q_DIR=http://kx.com/q
K_H=$(Q_DIR)/c/c/k.h
Q_LIB=$(Q_DIR)/$(PLATFORM)/q.lib


# Main build target
.PHONY: all clean distclean dist install
all: $(PLATFORM)/qml.dll


# Download targets
download/fdlibm-src.zip:
	mkdir -p download && wget -O '$@' '$(FDLIBM_SRC)'

download/cephes.zip:
	mkdir -p download && wget -O '$@' '$(CEPHES_SRC)'

download/cprob.tgz:
	mkdir -p download && wget -O '$@' '$(CPROB_SRC)'

download/labplot.tar.gz:
	mkdir -p download && wget -O '$@' '$(LABPLOT_SRC)'

download/k.h:
	mkdir -p download && wget -O '$@' '$(K_H)'

download/$(PLATFORM)/q.lib:
	mkdir -p 'download/$(PLATFORM)' && wget -O '$@' '$(Q_LIB)'

download/clapack.h:
	mkdir -p download && wget -O '$@' '$(CLAPACK_H)'

download/$(PLATFORM)/f2c.lib:
	mkdir -p 'download/$(PLATFORM)' && wget -O '$@' '$(F2C_LIB)'

download/$(PLATFORM)/blas.lib:
	mkdir -p 'download/$(PLATFORM)' && wget -O '$@' '$(BLAS_LIB)'

download/$(PLATFORM)/clapack.lib:
	mkdir -p 'download/$(PLATFORM)' && wget -O '$@' '$(CLAPACK_LIB)'

download/f2c.exe.gz:
	mkdir -p download && wget -O '$@' '$(F2C_BIN_GZ)'

download/f2c.h:
	mkdir -p download && wget -O '$@' '$(F2C_H)'

download/cpoly.f:
	mkdir -p download && wget -O '$@' '$(CPOLY_F)'

download/rpoly.f:
	mkdir -p download && wget -O '$@' '$(RPOLY_F)'

download/cpoly.c:
	mkdir -p download && wget -O '$@' '$(CPOLY_C)'


# Build FDLIBM
fdlibm/.extracted: download/fdlibm-src.zip
	mkdir -p fdlibm && unzip -o -j '$<' -d fdlibm
	touch '$@'

fdlibm/.patched: fdlibm/.extracted
	rm -f fdlibm/libm-dllversion.c fdlibm/w_gamma.c
	touch '$@'

fdlibm/fdlibm.lib: fdlibm/.patched
	cd fdlibm && $(CL) $(DEFINES) $(CFLAGS) -c -Oi- '*.c'
	cd fdlibm && $(LINK) -lib -out:fdlibm.lib '*.obj'

include/fdlibm.h: fdlibm/.patched
	mkdir -p include && cp -f 'fdlibm/fdlibm.h' '$@'

lib/fdlibm.lib: fdlibm/fdlibm.lib
	mkdir -p lib && cp -f '$<' '$@'


# Patch and build Cephes library
ifeq "$(BUILD)" "gpl"
    cephes/.extracted: download/labplot.tar.gz
	tar xf '$<' --strip-components 1 --wildcards LabPlot-*/cephes
	touch '$@'
else ifeq "$(BUILD)" "bsd"
    cephes/.extracted: download/cprob.tgz
	mkdir -p cephes
	tar xf download/cprob.tgz -C cephes
	touch '$@'
else
    cephes/.extracted: download/cephes.zip
	mkdir -p cephes && unzip -o '$<' -d cephes
	touch '$@'
endif

cephes/.patched: cephes/.extracted
	sed -i -e '$$a\' -e '#define IBMPC 1' -e '/# *define  *UNK/d' \
	    cephes/mconf.h
	sed -i -e '/( *!isfinite(x) *)/s//(x==INFINITY||x==-INFINITY)/' \
	    cephes/gamma.c
	sed -i -e 's/erf/c_erf/g;1i\' \
	    -e 'extern double c_erf(double), c_erfc(double);' cephes/ndtr.c
	sed -i -e 's/y01/y0/g;/double igami(/{N;N;N;N;N;' \
	    -e 's/$$/if(y0>.5||y0<=0)return(NAN);/}' cephes/igami.c
	mv -f cephes/const.c cephes/const.c_
	sed -i -e 's/true_gamma/gamma/g' cephes/*.h
	sed -i -e 's/true_gamma/gamma/g' \
	    -e 's/char st\{0,1\}\[\]="[[:alnum:]][[:alnum:]]*";//' \
	    -e 's/mtherr( *\("\|st\{0,1\} *,\)[^)]*)/return(NAN)/;1i\' \
	    -e 'extern double NAN;' cephes/*.c
	mv -f cephes/const.c_ cephes/const.c
	if test -e cephes/simq.c; then \
	    sed -i -e 's/printf([^)]*) *;//' cephes/simq.c; fi
	rm -f cephes/{floor,log,atan,pow,exp,lmdif,mtherr,mod2pi,dtestvec}.c
    ifeq "$(PLATFORM)" "w32"
	rm -f cephes/setprec.c
    else
	rm -f cephes/setpmsvc.c
    endif
	touch '$@'

cephes/cephes.lib: cephes/.patched
	cd cephes && $(CL) $(DEFINES) $(CFLAGS) -c -Oi- '*.c'
	cd cephes && $(LINK) -lib -out:cephes.lib '*.obj'

lib/cephes.lib: cephes/cephes.lib
	mkdir -p lib && cp -f '$<' '$@'


# Install q glue
include/k.h: download/k.h
	mkdir -p include && cp -f '$<' '$@'

lib/q.lib: download/$(PLATFORM)/q.lib
	mkdir -p lib && cp -f '$<' '$@'


# Install CLAPACK libraries
include/clapack.h: download/clapack.h
	mkdir -p include && cp -f '$<' '$@'

lib/f2c.lib: download/$(PLATFORM)/f2c.lib
	mkdir -p lib && cp -f '$<' '$@'

lib/blas.lib: download/$(PLATFORM)/blas.lib
	mkdir -p lib && cp -f '$<' '$@'

lib/clapack.lib: download/$(PLATFORM)/clapack.lib
	mkdir -p lib && cp -f '$<' '$@'


# Build cpoly
ifeq "$(BUILD)" "gpl"
    cpoly/cpoly.c: download/cpoly.c
	mkdir -p cpoly && cp -f '$<' '$@'
	sed -i -e 	's/Rboolean *\*fail/int *fail, double* work/' \
	    -e 's/Rboolean/int/;s/FALSE/0/;s/TRUE/1/;s/R_PosInf/INFINITY/' \
	    -e 's/<R_ext\/Arith\.h>/<fdlibm.h>/;/<R_ext\/Memory\.h>/d' \
	    -e 's/R_alloc[^;]*/work/;s/<Rmath\.h>/<math.h>/;s/R_pow_di/pow/' \
	    -e '/<R_ext\/Applic\.h>/d;1i\' -e 'extern double INFINITY;' \
	    cpoly/cpoly.c

    cpoly/cpoly.lib: cpoly/cpoly.c include/fdlibm.h lib/cephes.lib
	cd cpoly && $(CL) $(DEFINES) -I../include \
	    -DHAVE_HYPOT -D_USE_MATH_DEFINES $(CFLAGS) -c -Oi- cpoly.c
	cd cpoly && $(LINK) -lib -out:cpoly.lib cpoly.obj ../lib/cephes.lib

    lib/cpoly.lib: cpoly/cpoly.lib
	mkdir -p lib && cp -f '$<' '$@'
else ifeq "$(BUILD)" "bsd"
    # only need LAPACK
else
    cpoly/f2c.exe: download/f2c.exe.gz
	mkdir -p cpoly && gunzip -c '$<' >'$@' && chmod +x '$@'

    cpoly/f2c.h: download/f2c.h
	mkdir -p cpoly && cp -f '$<' '$@'

    cpoly/cpoly.orig.f cpoly/rpoly.orig.f: cpoly/%.orig.f: download/%.f
	mkdir -p cpoly && cp -f '$<' '$@'

    cpoly/cpoly.orig.c cpoly/rpoly.orig.c: %.c: %.f cpoly/f2c.exe
	cd cpoly && ./f2c.exe -a '$(notdir $<)'

    cpoly/cpoly.c cpoly/rpoly.c: %.c: %.orig.c
	sed -e 's/\([0-9][0-9]*e[0-9][0-9]*\)f/\1/;s/\(\.[0-9][0-9]*\)f/\1/' \
	    -e '/extern \/\* Subr/{/^[^;]*$$/{:n' -e 'N;s/^[^;]*$$/\0/;t n' \
	    -e '};};s/^\(struct \){/\1Global {/;s/^} global_;/};/' \
	    -e 's/_(/\0global_,/g' \
	    -e 's/\(int\|doublereal\) [^ ]*_(/\0struct Global* /g' \
	    -e '/extern \/\* Subr/s/,[^,_]*_(/\0struct Global* /g' \
	    -e 's/\([cr]poly_(\)[^,]*,/\1/g;/int [cr]poly_/{N;N' \
	    -e 's/$$/struct Global global,*global_=\&global;/;}' \
	    -e 's/^\(#define global_1 \)\(global_\)/\1(*\2)/;s/,void//' \
	    -e '/\/\* Main program alias/d' '$<' > '$@'

    cpoly/cpoly.lib: cpoly/cpoly.c cpoly/rpoly.c cpoly/f2c.h
	cd cpoly && $(CL) $(DEFINES) -DMSDOS $(CFLAGS) -c -Oi- cpoly.c rpoly.c
	cd cpoly && $(LINK) -lib -out:cpoly.lib cpoly.obj rpoly.obj

    lib/cpoly.lib: cpoly/cpoly.lib
	mkdir -p lib && cp -f '$<' '$@'
endif


# Build QML
VERSION=0.1.1
CONFIG=-DVERSION=$(VERSION)
ifeq "$(BUILD)" "gpl"
    CONFIG+= -DUSE_R_POLY
else ifeq "$(BUILD)" "bsd"
    CONFIG+= -DUSE_LAPACK_POLY
endif

DEFINES+= -D_USE_MATH_DEFINES

SOURCES=qml.c include/fdlibm.h include/k.h include/clapack.h
LIBS=lib/fdlibm.lib lib/cephes.lib lib/q.lib lib/f2c.lib \
    lib/blas.lib lib/clapack.lib

ifneq "$(BUILD)" "bsd"
    LIBS+= lib/cpoly.lib
endif

ifeq "$(PLATFORM)" "w64"
    lib/compat.c:
	mkdir -p lib && echo "void __fastcall __GSHandlerCheck() {}\
	    void __fastcall __security_check_cookie(unsigned* p) {}\
	    unsigned* __security_cookie;" >"$@"
    SOURCES+= lib/compat.c
endif

build/qml.dll: $(SOURCES) $(LIBS)
	mkdir -p build
	cd build && $(CL) $(DEFINES) $(CONFIG) $(CFLAGS) -LD -I../include \
	    $(patsubst %,../%,$(filter-out %.h,$(SOURCES))) \
	    -link $(patsubst %,../%,$(LIBS))

$(PLATFORM)/qml.dll: build/qml.dll
	mkdir -p '$(PLATFORM)' && cp -pf '$<' '$@'


# Create distributable archive
DIST=LICENSE.txt LICENSE_CLAPACK.txt LICENSE_LIBF2C.txt LICENSE_Q.txt \
    README.txt Makefile qml.c qml.q test.q w32/qml.dll w64/qml.dll

dist: $(DIST)
	rm -f 'qml-$(VERSION).zip'
	7z a -tzip -mx=9 'qml-$(VERSION).zip' $(DIST)


# Clean up
clean:
	rm -rf cephes fdlibm cpoly include lib build

distclean: clean
	rm -rf download '$(PLATFORM)'
	rm -f 'qml-$(VERSION).zip'
