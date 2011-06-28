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

ifeq "$(BUILD)" "private"
else ifeq "$(BUILD)" "gpl"
else
    BUILD=bsd
endif
$(info Build type: $(BUILD))


# Platform
#   w32
#   w64

ifeq "$(PLATFORM)" "w32"
else ifeq "$(PLATFORM)" "w64"
else
    ifeq "$(patsubst CYGWIN%,,$(shell uname -s).)" ""
        ifeq "$(patsubst %WOW64,,.$(shell uname -s))" ""
            PLATFORM=w64
        else
            PLATFORM=w32
        endif
    else ifeq "$(shell uname -s)" "Linux"
        ifeq "$(shell uname -m)" "x86_64"
            $(error l64 platform is currently not supported)
	    #PLATFORM=l64
        else
            PLATFORM=l32
        endif
    else
        $(error couldn't determine platform, please set PLATFORM variable)
    endif
endif
$(info Platform:   $(PLATFORM))


# Toolchain
#   mssdk:
#     Microsoft SDK
#
#   gnu:
#     GNU

ifeq "$(TOOLCHAIN)" "mssdk"
else
    TOOLCHAIN=gnu
endif
$(info Toolchain:  $(TOOLCHAIN))

ifeq "$(TOOLCHAIN)" "gnu"
    ifneq "$(TOOLPREFIX)" ""
        $(info Tool prefix: $(TOOLPREFIX))
    endif
endif


# Compiler and flags
DEFINES=-D_IEEE_LIBM -D__LITTLE_ENDIAN
CFLAGS=

ifeq "$(TOOLCHAIN)" "mssdk"
    ifeq "$(PLATFORM)" "w64"
        MSSDK_BIN=$$MSSDK/Bin/win64/x86/AMD64
    else
        MSSDK_BIN=$$MSSDK/vc/bin
    endif
    CL="$(MSSDK_BIN)/cl"
    LINK="$(MSSDK_BIN)/link"
    OBJEXT=obj
    LIBEXT=lib
    DLLEXT=dll
    CFLAGS+= -Ox -Oi-
    DEFINES+= -D__STDC__ -DWIN32 -D_USE_MATH_DEFINES
    cc=$(CL) $(DEFINES) $(2) $(CFLAGS) -c $(patsubst %,'%',$(1))
    ar=$(LINK) -lib -out:$(1).$(LIBEXT) $(patsubst %,'%',$(2))
    ccdll=$(CL) $(DEFINES) $(4) $(CFLAGS) -LD $(patsubst %,'%',$(2)) \
          -link $(patsubst %,'%',$(3)) -out:$(1).$(DLLEXT)
else
    ifneq "$(TOOLPREFIX)" ""
        TOOLPREPEND=$(TOOLPREFIX)-
    else
        TOOLPREPEND=
    endif
    CC=$(TOOLPREPEND)gcc
    LD=$(TOOLPREPEND)ld
    AR=$(TOOLPREPEND)ar
    RANLIB=$(TOOLPREPEND)ranlib
    DLLTOOL=$(TOOLPREPEND)dlltool
    NM=$(TOOLPREPEND)nm
    AS=$(TOOLPREPEND)as
        OBJEXT=o
        LIBEXT=a
    ifeq "$(patsubst w%,,$(PLATFORM))" ""
        DLLEXT=dll
        CFLAGS+= -mno-cygwin
    else
        CFLAGS+= -fPIC -fvisibility=hidden
        DLLEXT=so
    endif
    CFLAGS+= -std=gnu99 -O2 -fno-strict-aliasing -pipe
    CFLAGS+= -Wall -Wno-parentheses -Wno-uninitialized
    # -fstrict-aliasing breaks fdlibm
    cc=$(CC) $(DEFINES) $(2) $(CFLAGS) -c $(1)
    ar=$(AR) r $(1).$(LIBEXT) $(2)
    ccdll=$(CC) $(DEFINES) $(4) $(CFLAGS) -shared -Wl,-s \
          -o $(1).$(DLLEXT) $(2) $(3)
endif


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
CLAPACK_SRC=$(CLAPACK_DIR)/clapack.tgz
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
.PHONY: all cleanpart clean distclean dist install
all: $(PLATFORM)/qml.$(DLLEXT)


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

download/clapack.tgz:
	mkdir -p download && wget -O '$@' '$(CLAPACK_SRC)'

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
	sed -i $(foreach const,DOMAIN SING OVERFLOW UNDERFLOW TLOSS PLOSS __P,\
	    -e '/#define[[:space:]]\{1,\}$(const)\b/i\' \
	    -e '#undef $(const)') fdlibm/fdlibm.h
	rm -f fdlibm/libm-dllversion.c fdlibm/w_gamma.c fdlibm/*.o
	touch '$@'

fdlibm/fdlibm.$(LIBEXT): fdlibm/.patched
	cd fdlibm && $(call cc,*.c)
	cd fdlibm && $(call ar,fdlibm,*.$(OBJEXT))

include/fdlibm.h: fdlibm/.patched
	mkdir -p include && cp -f 'fdlibm/fdlibm.h' '$@'

lib/fdlibm.$(LIBEXT): fdlibm/fdlibm.$(LIBEXT)
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
	    -e 's/mtherr( *\("\|\(s\|st\|fname\) *,\)[^)]*)/return(NAN)/;1i\' \
	    -e 'extern double NAN;' cephes/*.c
	mv -f cephes/const.c_ cephes/const.c
	if test -e cephes/simq.c; then \
	    sed -i -e 's/printf([^)]*) *;//' cephes/simq.c; fi
	rm -f cephes/{floor,log,atan,pow,exp,lmdif,mtherr,mod2pi,dtestvec}.c
    ifeq "$(PLATFORM)" "w32"
	rm -f cephes/setprec.c
    endif
	rm -f cephes/setpmsvc.c
	touch '$@'

cephes/cephes.$(LIBEXT): cephes/.patched
	cd cephes && $(call cc,*.c)
	cd cephes && $(call ar,cephes,*.$(OBJEXT))

lib/cephes.$(LIBEXT): cephes/cephes.$(LIBEXT)
	mkdir -p lib && cp -f '$<' '$@'


# Install q glue
ifeq "$(TOOLCHAIN)" "mssdk"
    include/k.h: download/k.h
	mkdir -p include && cp -f '$<' '$@' && echo >>'$@'

    lib/q.lib: download/$(PLATFORM)/q.lib
	mkdir -p lib && cp -f '$<' '$@'
else
    include/k.h: download/k.h
	mkdir -p qlib && sed -e '/#define  *isnan/d;$$a\' -e '' '$<' >'$@'

    qlib/q.orig.lib: download/$(PLATFORM)/q.lib
	mkdir -p qlib && cp -f '$<' '$@'

    qlib/q.def: qlib/q.orig.lib
	$(NM) '$<' | awk "/^[[:alnum:].]+:/&&!name[1]\
	    {split(\$$1,name,\":\");print\"NAME\",name[1];print\"EXPORTS\"}\
	    /[[:space:]]T[[:space:]]+_?[[:alnum:]]/{sub(\"^_\",\"\",\$$3);\
	    print \$$3;c=1}END{exit !c}" >'$@.tmp'
	mv -f '$@.tmp' '$@'

    qlib/q.$(LIBEXT): qlib/q.def
	$(DLLTOOL) -d '$<' -l '$@'

    lib/q.$(LIBEXT): qlib/q.$(LIBEXT)
	mkdir -p lib && cp -f '$<' '$@'
endif


# Build CLAPACK libraries
ifeq "$(TOOLCHAIN)" "mssdk"
    include/clapack.h: download/clapack.h
	mkdir -p include && cp -f '$<' '$@'

    lib/f2c.$(LIBEXT): download/$(PLATFORM)/f2c.lib
	mkdir -p lib && cp -f '$<' '$@'

    lib/blas.$(LIBEXT): download/$(PLATFORM)/blas.lib
	mkdir -p lib && cp -f '$<' '$@'

    lib/clapack.$(LIBEXT): download/$(PLATFORM)/clapack.lib
	mkdir -p lib && cp -f '$<' '$@'
else
    clapack/.extracted: download/clapack.tgz
	mkdir -p clapack
	tar xf '$<' -C clapack --strip-components 1 --wildcards CLAPACK-*
	touch '$@'

    clapack/.patched: clapack/.extracted
	{ cat clapack/make.inc.example;  echo 'PLAT=_$(PLATFORM)'; \
	    echo 'CC=$(CC) -DNO_BLAS_WRAP'; echo 'CFLAGS=$(CFLAGS)'; \
	    echo 'NOOPT=$$(CFLAGS) -O0'; echo 'LOADOPTS=$$(CFLAGS)'; \
	    echo 'LOADER=$$(CC)'; echo 'LD=$(LD)'; echo 'ARCH=$(AR)'; \
	    echo 'RANLIB=$(RANLIB)'; } >clapack/make.inc
      ifeq "$(patsubst w%,,$(PLATFORM))" ""
	sed -i -e 's/\ba\.out\b/a.exe/g' clapack/F2CLIBS/libf2c/Makefile
      endif
	sed -i -e 's/^\(TIME[[:space:]]*=\).*/\1/' \
	    -e 's/\bcc\b/$$(CC)/;s/\bld\b/$$(LD)/;s/\bar\b/$$(ARCH)/' \
            -e 's/\branlib\b/$$(RANLIB)/' clapack/F2CLIBS/libf2c/Makefile
	sed -i -e '/^#define MSpc\b/d' clapack/F2CLIBS/libf2c/uninit.c
	{ echo '#include "f2c.h"'; echo 'doublereal second_(){return 0;}'; } \
	    >clapack/INSTALL/second.c
	sed -e 's/second/dsecnd/g' clapack/INSTALL/second.c \
	    >clapack/INSTALL/dsecnd.c
	find clapack -name '*.h' -type l -exec cp '{}' clapack/.tmp.h ';' \
	    -exec rm -- '{}' ';' -exec mv clapack/.tmp.h {} ';'
	touch '$@'

    clapack/.built: clapack/.patched
	cd clapack && make f2clib
	rm -f clapack/{BLAS/,}SRC/{f2c,blaswrap}.h
	cp -f clapack/INCLUDE/f2c.h clapack/BLAS/SRC
	cp -f clapack/INCLUDE/f2c.h clapack/SRC
	cp -f clapack/INCLUDE/blaswrap.h clapack/BLAS/SRC
	cp -f clapack/INCLUDE/blaswrap.h clapack/SRC
	cd clapack && make blaslib lapacklib
	touch '$@'

    include/clapack.h: clapack/.built
	mkdir -p include && cp -f clapack/INCLUDE/clapack.h '$@'

    lib/f2c.$(LIBEXT): clapack/.built
	mkdir -p lib && cp -f clapack/F2CLIBS/libf2c.$(LIBEXT) '$@'

    lib/blas.$(LIBEXT): clapack/.built
	mkdir -p lib && cp -f clapack/blas_$(PLATFORM).$(LIBEXT) '$@'

    lib/clapack.$(LIBEXT): clapack/.built
	mkdir -p lib && cp -f clapack/lapack_$(PLATFORM).$(LIBEXT) '$@'
endif


# Build cpoly
ifeq "$(BUILD)" "gpl"
    cpoly/cpoly.c: download/cpoly.c
	mkdir -p cpoly && cp -f '$<' '$@'
	sed -i -e 's/Rboolean *\*fail/int *fail, double* work/' \
	    -e 's/Rboolean/int/;s/FALSE/0/;s/TRUE/1/;s/R_PosInf/INFINITY/' \
	    -e 's/<R_ext\/Arith\.h>/<fdlibm.h>/;/<R_ext\/Memory\.h>/d' \
	    -e 's/R_alloc[^;]*/work/;/<Rmath\.h>/d;s/R_pow_di/pow/' \
	    -e 's/\bM_SQRT2\b/1.4142135623730950488/' \
	    -e 's/\bM_SQRT1_2\b/0.7071067811865475244/' \
	    -e '/<R_ext\/Applic\.h>/d;1i\' -e 'extern double INFINITY;' \
	    cpoly/cpoly.c

    cpoly/cpoly.$(LIBEXT): cpoly/cpoly.c include/fdlibm.h lib/cephes.$(LIBEXT)
	cd cpoly && $(call cc,cpoly.c,-I../include -DHAVE_HYPOT)
	cd cpoly && $(call ar,cpoly,cpoly.$(OBJEXT) ../lib/cephes.$(LIBEXT))

    lib/cpoly.$(LIBEXT): cpoly/cpoly.$(LIBEXT)
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

    cpoly/cpoly.$(LIBEXT): cpoly/cpoly.c cpoly/rpoly.c cpoly/f2c.h
	cd cpoly && $(call cc,cpoly.c rpoly.c,-DMSDOS)
	cd cpoly && $(call ar,cpoly,cpoly.$(OBJEXT) rpoly.$(OBJEXT))

    lib/cpoly.$(LIBEXT): cpoly/cpoly.$(LIBEXT)
	mkdir -p lib && cp -f '$<' '$@'
endif


# Build QML
VERSION=0.1.2
CONFIG=-DVERSION=$(VERSION)
ifeq "$(BUILD)" "gpl"
    CONFIG+= -DUSE_R_POLY
else ifeq "$(BUILD)" "bsd"
    CONFIG+= -DUSE_LAPACK_POLY
endif

SOURCES=qml.c include/fdlibm.h include/k.h include/clapack.h
LIBS=
ifeq "$(patsubst w%,,$(PLATFORM))" ""
    CONFIG+= -DDLLEXPORT
    LIBS+= lib/q.$(LIBEXT)
else
    CONFIG+= -DSOEXPORT
endif
ifneq "$(BUILD)" "bsd"
    LIBS+= lib/cpoly.$(LIBEXT)
endif
LIBS+=lib/cephes.$(LIBEXT) lib/fdlibm.$(LIBEXT)
LIBS+=lib/clapack.$(LIBEXT) lib/blas.$(LIBEXT) lib/f2c.$(LIBEXT)

ifeq "$(TOOLCHAIN)" "mssdk"
ifeq "$(PLATFORM)" "w64"
    lib/compat.c:
	mkdir -p lib && echo "void __fastcall __GSHandlerCheck() {}\
	    void __fastcall __security_check_cookie(unsigned* p) {}\
	    unsigned* __security_cookie;" >"$@"
    SOURCES+= lib/compat.c
endif
endif

build/qml.$(DLLEXT): $(SOURCES) $(LIBS)
	mkdir -p build
	cd build && $(call ccdll,qml,\
	    $(patsubst %,../%,$(filter-out %.h,$(SOURCES))),\
	    $(patsubst %,../%,$(LIBS)),\
	    -I../include $(CONFIG))

$(PLATFORM)/qml.$(DLLEXT): build/qml.$(DLLEXT)
	mkdir -p '$(PLATFORM)' && cp -pf '$<' '$@'


# Create distributable archive
DIST=LICENSE.txt LICENSE_CLAPACK.txt LICENSE_LIBF2C.txt LICENSE_Q.txt \
    README.txt Makefile qml.c qml.q test.q w32/qml.dll w64/qml.dll l32/qml.so

dist: $(DIST)
	rm -f 'qml-$(VERSION).zip'
	7z a -tzip -mx=9 'qml-$(VERSION).zip' $(DIST)


# Clean up
cleanpart:
	rm -rf cephes fdlibm cpoly qlib include lib build

clean: cleanpart
	rm -rf clapack

distclean: clean
	rm -rf download '$(PLATFORM)'
	rm -f 'qml-$(VERSION).zip'
