# Make configuration options:
#   make PLATFORM=l32 (default is to autodetect)
#   make BLAS=atlas   (default is blas)

ifeq "$(findstring $(PLATFORM),w32 w64 l32 l64 m32 m64 v32 v64 s32 s64)" ""
  ifeq "$(patsubst CYGWIN%,,$(shell uname -s).)" ""
    ifeq "$(patsubst %WOW64,,.$(shell uname -s))" ""
      PLATFORM=w64
    else
      PLATFORM=w32
    endif
  else
    ifeq "$(shell uname -s)" "Linux"
      ifeq "$(shell uname -m)" "x86_64"
        $(warning l64 platform is not tested)$(shell sleep 3)
        PLATFORM=l64
      else
        PLATFORM=l32
      endif
    else
      ifeq "$(shell uname -s)" "Darwin"
        ifeq "$(shell uname -m)" "x86_64"
          $(warning m64 platform is not tested)$(shell sleep 5)
          PLATFORM=m64
        else
          PLATFORM=m32
        endif
      else
        ifeq "$(shell uname -s)" "SunOS"
          ifneq "$(findstring i386,$(shell isainfo))" ""
            ifeq "$(shell isainfo -b)" "64"
              $(warning v64 platform is not tested)$(shell sleep 5)
              PLATFORM=v64
            else
              $(warning v32 platform is not tested)$(shell sleep 5)
              PLATFORM=v32
            endif
          else
            ifeq "$(shell isainfo -b)" "64"
              $(warning s64 platform is not tested)$(shell sleep 5)
              PLATFORM=s64
            else
              $(warning s32 platform is not tested)$(shell sleep 5)
              PLATFORM=s32
            endif
          endif
        else
          $(error couldn't determine platform, please set PLATFORM variable)
        endif
      endif
    endif
  endif
endif
$(info Platform:   $(PLATFORM))

PLATFORMCLASS=$(patsubst v,s,$(patsubst %64,%,$(patsubst %32,%,$(PLATFORM))))
PLATFORMBITS=$(patsubst %64,64,$(patsubst %32,32,$(PLATFORM)))

ifneq "$(BLAS)" "atlas"
    BLAS=blas
endif
$(info BLAS:       $(BLAS))


# Compiler and flags
DEFINES=-D_IEEE_LIBM -D__LITTLE_ENDIAN -D_REENTRANT
CFLAGS=

TOOLPREPEND=
ifneq "$(TOOLPREFIX)" ""
    TOOLPREPEND=$(TOOLPREFIX)-
endif
GCCAPPEND=
ifneq "$(GCCSUFFIX)" ""
    GCCAPPEND=-$(GCCSUFFIX)
endif
CC=$(TOOLPREPEND)gcc$(GCCAPPEND)
FC=$(TOOLPREPEND)gfortran$(GCCAPPEND)
LD=$(TOOLPREPEND)ld
AR=$(TOOLPREPEND)ar
RANLIB=$(TOOLPREPEND)ranlib
DLLTOOL=$(TOOLPREPEND)dlltool
NM=$(TOOLPREPEND)nm
AS=$(TOOLPREPEND)as
OBJEXT=o
CFLAGS+= -std=gnu99 -O2 -fno-strict-aliasing -ffloat-store -pipe
# -fstrict-aliasing breaks fdlibm
# -fno-float-store seems to break incbi() on Linux
CFLAGS+= -fno-builtin-sin -fno-builtin-cos -fno-builtin-sqrt
CFLAGS+= -Wall -Wno-parentheses -Wno-uninitialized
ENVFLAGS=
ifeq "$(PLATFORMCLASS)" "w"
    DLLEXT=dll
    EXEEXT=.exe
    CFLAGS+= -mno-cygwin
else
    DLLEXT=so
    EXEEXT=
    CFLAGS+= -fPIC
    # explicit binary type
    ifeq "$(patsubst %64,,$(PLATFORM))" ""
        CFLAGS+= -m64
    else
        CFLAGS+= -m32
    endif
endif
ifeq "$(PLATFORMCLASS)" "m"
    # Darwin flags
    SOFLAGS=-bundle -undefined dynamic_lookup -nodefaultlibs -Wl,-x
    ENVFLAGS=MACOSX_DEPLOYMENT_TARGET=10.4
else
    ifeq "$(PLATFORMCLASS)" "s"
        # OpenSolaris flags
        SOFLAGS=-G -Wl,-s,-Bsymbolic
    else
        # default flags
        SOFLAGS=-shared -Wl,-s,-Bsymbolic
    endif
endif
cc=$(CC) $(DEFINES) $(2) $(CFLAGS) -c $(1)
ar=$(AR) r $(1).a $(2)
ccdll=$(ENVFLAGS) $(CC) $(DEFINES) $(4) $(CFLAGS) $(SOFLAGS)
ifeq "$(PLATFORMCLASS)" "l"
    ccdll+= -Wl,--version-script,$(1).mapfile
endif
ifeq "$(PLATFORMCLASS)" "m"
    ccdll+= -exported_symbols_list $(1).symlist
endif
ifeq "$(PLATFORMCLASS)" "s"
ccdll+= -Wl,-M,$(1).mapfile
endif
ccdll+= -o $(1).$(DLLEXT) $(patsubst %,'%',$(2)) $(patsubst %,'%',$(3))
ccexe=$(ENVFLAGS) $(CC) $(DEFINES) $(4) $(CFLAGS) -o $(1)$(EXEEXT) \
    $(patsubst %,'%',$(2)) $(patsubst %,'%',$(3))


# Dependecy locations
SOURCEFORGE_DIR=http://prdownloads.sourceforge.net
NETLIB_DIR=http://www.netlib.org
#   FDLIBM
FDLIBM_SRC=$(SOURCEFORGE_DIR)/gnuwin32/fdlibm-5.2-2-src.zip
#   Cephes
CPROB_SRC=$(NETLIB_DIR)/cephes/cprob.tgz
#   CLAPACK
ATLAS_SRC=$(SOURCEFORGE_DIR)/math-atlas/Stable/3.8.4/atlas3.8.4.tar.bz2
CLAPACK_SRC=$(NETLIB_DIR)/clapack/clapack.tgz
#   f2c program
F2C_BIN_GZ=$(NETLIB_DIR)/f2c/mswin/f2c.exe.gz
F2C_SRC=ftp://ftp.netlib.org/f2c/src.tar.gz
F2C_H=$(NETLIB_DIR)/f2c/f2c.h
#   conmax
CONMAX_F=$(NETLIB_DIR)/opt/conmax.f
#   q stuff
Q_DIR=http://kx.com/q
K_H=$(Q_DIR)/c/c/k.h
Q_LIB=$(Q_DIR)/$(PLATFORM)/q.lib


# Main build target
.PHONY: all test cleanpart clean distclean dist install
all: $(PLATFORM)/qml.$(DLLEXT)


# Download targets
download/fdlibm-src.zip:
	mkdir -p download && wget -O '$@' '$(FDLIBM_SRC)'

download/cprob.tgz:
	mkdir -p download && wget -O '$@' '$(CPROB_SRC)'

download/k.h:
	mkdir -p download && wget -O '$@' '$(K_H)'

download/$(PLATFORM)/q.lib:
	mkdir -p 'download/$(PLATFORM)' && wget -O '$@' '$(Q_LIB)'

download/atlas.tar.bz2:
	mkdir -p download && wget -O '$@' '$(ATLAS_SRC)'

download/clapack.tgz:
	mkdir -p download && wget -O '$@' '$(CLAPACK_SRC)'

download/f2c.exe.gz:
	mkdir -p download && wget -O '$@' '$(F2C_BIN_GZ)'

download/f2c.tar.gz:
	mkdir -p download && wget -O '$@' '$(F2C_SRC)'

download/f2c.h:
	mkdir -p download && wget -O '$@' '$(F2C_H)'

download/conmax.f:
	mkdir -p download && wget -O '$@' '$(CONMAX_F)'


# Build FDLIBM
fdlibm/.extracted: download/fdlibm-src.zip
	mkdir -p fdlibm && unzip -o -j '$<' -d fdlibm
	touch '$@'

fdlibm/.patched: fdlibm/.extracted
	sed -i.tmp $(foreach const,\
	    DOMAIN SING OVERFLOW UNDERFLOW TLOSS PLOSS __P,\
	    -e '/#define[[:space:]]*$(const)/{h;s/.*/#undef $(const)/p;g;}') \
	    fdlibm/fdlibm.h
	rm -f fdlibm/libm-dllversion.c fdlibm/w_gamma.c fdlibm/*.o
	touch '$@'

fdlibm/fdlibm.a: fdlibm/.patched
	cd fdlibm && $(call cc,*.c)
	cd fdlibm && $(call ar,fdlibm,*.$(OBJEXT))

include/fdlibm.h: fdlibm/.patched
	mkdir -p include && cp -f 'fdlibm/fdlibm.h' '$@'

lib/fdlibm.a: fdlibm/fdlibm.a
	mkdir -p lib && cp -f '$<' '$@'
	$(RANLIB) '$@'


# Patch and build Cephes library
cephes/.extracted: download/cprob.tgz
	mkdir -p cephes
	tar xzf download/cprob.tgz -C cephes
	touch '$@'

cephes/.patched: cephes/.extracted
	sed -i.tmp -e '/#[[:space:]]*define[[:space:]]\{1,\}UNK/d' \
	    -e '$${p;s/.*/#define IBMPC 1/;}' cephes/mconf.h
	sed -i.tmp -e 's/erf/c_erf/g' \
	    -e '1{h;s/.*/extern double c_erf(double), c_erfc(double);/p;g;}' \
	    cephes/ndtr.c
	sed -i.tmp -e 's/y01/y0/g;/double  *igami(/{:a' \
	    -e 's/{/&if(y0<=0||y0>=1)return(NAN);/;t b' -e 'N;b a' -e ':b' \
	    -e '};' cephes/igami.c
	sed -i.tmp -e 's/\(for *( *i *= *0; *i *< *\)1\(00 *;\)/\15\2/' \
	    cephes/incbi.c
	for file in cephes/*.h; do \
	    sed -i.tmp -e 's/true_gamma/gamma/g' "$$file"; done
	for file in $(filter-out %/const.c,$(wildcard cephes/*.c)); do \
	    sed -i.tmp -e 's/true_gamma/gamma/g' \
	    -e 's/char st\{0,1\}\[\]="[[:alnum:]][[:alnum:]]*";//g' \
	    -e 's/mtherr *( *["sf][^,]*, *PLOSS *)//g' \
	    -e 's/mtherr *( *["sf][^,]*,[^)]*)/return(NAN)/g' \
	    -e '1{h;s/.*/extern double NAN;/p;g;}' "$$file"; done
	sed -i.tmp -n -e 'H;$${x;' \
	    -e 's/( *!isfinite(x) *)/(x==INFINITY||x==-INFINITY)/g' \
	    -e 's/\(extern  *\)\{0,1\}int  *sgngam[^;]*;//g' \
	    -e 's/double  *gamma([^)]*)[^;][^{]*{/&int sgngam;/g' \
	    -e 's/double  *lgam([^)]*)[^;][^{]*{/&int sgngam;/g' \
	    -e 'p;}' cephes/gamma.c
	sed -i.tmp -e '/if *( *n *<= *0[^0-9]/{N' \
	    -e 's/\(return *(\)-1\.[^)]*/\1NAN/;};/^ *kolmogorov *(/{:a' \
	    -e 's/{/&if(y<.12)return 1;/;t b' -e 'N;b a' -e ':b' -e '}' \
	    cephes/kolmogorov.c
	rm -f cephes/mtherr.c
	touch '$@'

cephes/cephes.a: cephes/.patched
	cd cephes && $(call cc,*.c)
	cd cephes && $(call ar,cephes,*.$(OBJEXT))

lib/cephes.a: cephes/cephes.a
	mkdir -p lib && cp -f '$<' '$@'
	$(RANLIB) '$@'


# Install q glue
include/k.h: download/k.h
	mkdir -p include && sed -e '/#define  *isnan/d;$${p;s/.*//;}' '$<' >'$@'

qlib/q.orig.lib: download/$(PLATFORM)/q.lib
	mkdir -p qlib && cp -f '$<' '$@'

qlib/q.def: qlib/q.orig.lib
	$(NM) '$<' | awk "/^[[:alnum:].]+:/&&!name[1]\
	    {split(\$$1,name,\":\");print\"NAME\",name[1];print\"EXPORTS\"}\
	    /[[:space:]]T[[:space:]]+_?[[:alnum:]]/{sub(\"^_\",\"\",\$$3);\
	    print \$$3;c=1}END{exit !c}" >'$@.tmp'
	mv -f '$@.tmp' '$@'

qlib/q.a: qlib/q.def
	$(DLLTOOL) -d '$<' -l '$@'

lib/q.a: qlib/q.a
	mkdir -p lib && cp -f '$<' '$@'


# Build CLAPACK libraries
atlas/.extracted: download/atlas.tar.bz2
	mkdir -p atlas && tar xjf '$<' -C atlas
	touch '$@'

atlas/.built: atlas/.extracted
	mkdir -p atlas/build
	cd atlas/build && ../ATLAS/configure -b $(PLATFORMBITS) \
	    -C ic '$(CC)' -F ic '$(CFLAGS)' \
	    -C if '$(FC)' -F if '$(CFLAGS)'
	cd atlas/build && make build
	touch '$@'

lib/atlas.a: atlas/.built
	mkdir -p lib && cp -f atlas/build/lib/libatlas.a '$@'

lib/cblas.a: atlas/.built
	mkdir -p lib && cp -f atlas/build/lib/libcblas.a '$@'

clapack/.extracted: download/clapack.tgz
	mkdir -p clapack && tar xzf '$<' -C clapack
	mv -f clapack/CLAPACK-*/* clapack/
	touch '$@'

clapack/.patched: clapack/.extracted
	# thread-safety fix:
	sed -i.tmp '/equiv_/s/static//g' \
	    clapack/SRC/slaln2.c clapack/SRC/dlaln2.c
	
	{ cat clapack/make.inc.example; \
	    echo 'PLAT=_$(PLATFORM)'; echo 'CC=$(CC)'; \
	    echo 'CFLAGS=$(CFLAGS) -I$$(TOPDIR)/INCLUDE'; \
	    echo 'NOOPT=$$(CFLAGS) -O0'; echo 'LOADOPTS=$$(CFLAGS)'; \
	    echo 'LOADER=$$(CC)'; echo 'LD=$(LD)'; echo 'ARCH=$(AR)'; \
	    echo 'RANLIB=$(RANLIB)'; } >clapack/make.inc
	sed -i.tmp -e 's/^\(TIME[[:space:]]*=\).*/\1/' \
	    -e 's/\([^_[:alnum:]]\)cc\([^_[:alnum:]]\)/\1$$(CC)\2/' \
	    -e 's/\([^_[:alnum:]]\)ld\([^_[:alnum:]]\)/\1$$(LD)\2/' \
	    -e 's/\([^_[:alnum:]]\)ar\([^_[:alnum:]]\)/\1$$(ARCH)\2/' \
	    -e 's/\([^_[:alnum:]]\)ranlib\([^_[:alnum:]]\)/\1$$(RANLIB)\2/' \
	    clapack/F2CLIBS/libf2c/Makefile
      ifeq "$(PLATFORMCLASS)" "w"
	sed -i.tmp -e 's/\([^_[:alnum:]]\)a\.out/\1a.exe/g' \
	    clapack/F2CLIBS/libf2c/Makefile
      endif
      ifneq "$(findstring $(PLATFORMCLASS),m s)" ""
	sed -i.tmp -e '/\$$(LD) .*\.xxx/d;/mv .*\.xxx/d' \
	    clapack/F2CLIBS/libf2c/Makefile
	sed -i.tmp 's,\./test[a-z]*,:,g' clapack/Makefile
	sed -i.tmp '/-o  *test[a-z]* /d' clapack/INSTALL/Makefile
      endif
	sed -i.tmp -e '/^#define MSpc\([^_[:alnum:]]\|$$\)/d' \
	    clapack/F2CLIBS/libf2c/uninit.c
	{ echo '#include "f2c.h"'; echo 'doublereal second_(){return 0;}'; } \
	    >clapack/INSTALL/second.c
	sed -e 's/second/dsecnd/g' clapack/INSTALL/second.c \
	    >clapack/INSTALL/dsecnd.c
      ifeq "$(PLATFORMCLASS)" "w"
	find clapack -name '*.h' -type l -exec cp '{}' clapack/.tmp.h ';' \
	    -exec rm -- '{}' ';' -exec mv clapack/.tmp.h '{}' ';'
      endif
	touch '$@'

clapack/.built: clapack/.patched
	cd clapack/F2CLIBS/libf2c && make hadd
	cd clapack && make f2clib
      ifeq "$(BLAS)" "blas"
	cd clapack && make blaslib
      else
	cd clapack && make cblaswrap
      endif
	cd clapack && make lapacklib
	touch '$@'

include/f2c.h: clapack/.built
	mkdir -p include && cp -f clapack/INCLUDE/f2c.h '$@'

include/clapack.h: clapack/.built
	mkdir -p include && cp -f clapack/INCLUDE/clapack.h '$@'

lib/f2c.a: clapack/.built
	mkdir -p lib && cp -f clapack/F2CLIBS/libf2c.a '$@'
	$(RANLIB) '$@'

lib/blas.a: clapack/.built
	mkdir -p lib && cp -f clapack/blas_$(PLATFORM).a '$@'
	$(RANLIB) '$@'

lib/cblaswr.a: clapack/.built
	mkdir -p lib && cp -f clapack/libcblaswr.a '$@'

lib/clapack.a: clapack/.built
	mkdir -p lib && cp -f clapack/lapack_$(PLATFORM).a '$@'
	$(RANLIB) '$@'


# Build f2c
ifeq "$(PLATFORMCLASS)" "w"
    f2c/f2c$(EXEEXT): download/f2c.exe.gz
	mkdir -p f2c && gunzip -c '$<' >'$@' && chmod +x '$@'
else
    f2c/.extracted: download/f2c.tar.gz
	mkdir -p f2c && tar xzf '$<' -C f2c
	mv -f f2c/src/* f2c/
	touch '$@'

    f2c/.patched: f2c/.extracted
	{ echo 'CC=$(CC)'; echo 'CFLAGS=$(DEFINES) $(CFLAGS)'; \
	    sed -e '/^ *CC *=/d;/^ *CFLAGS *=/d;/^f2c:/{n;' \
	        -e 's/\$$(CC) /&$$(CFLAGS) /;}' f2c/makefile.u; } >f2c/Makefile
	sed -i.tmp -e '/Cextern int \(unlink\|fork\)/d' f2c/sysdep.c
	touch '$@'

    f2c/f2c$(EXEEXT): f2c/.patched
	cd f2c && make f2c
endif


# Build conmax
conmax/conmax.orig.f: download/conmax.f
	mkdir -p conmax && cp -f '$<' '$@'

conmax/conmax.in.f: conmax/conmax.orig.f
	sed -e '/^C/d' \
	    -e 's/PTTBL(IPTB,INDM)/PTTBL(0:0,0:0)/g;s/IPTB,//g;s/INDM,//g' \
	    -e '/NE MULLER(/{h;s/(.*/(LIMMUL,NSRCH,/p;g;s/.*(/     */;}' \
	    -e '/LL MULLER(/{h;s/(.*/(5,     NSRCH,/p;g;s/.*(/     */;}' \
	    -e '/NE SEARSL(/{h;s/(.*/(INITLM,NADD,LIMS1,/p;g;s/.*(/     */;}' \
	    -e '/LL SEARSL(/{h;s/(.*/(6,     4,   LIMS1,/p;g;s/.*(/     */;}' \
	    -e '/^ *LIMMUL=5 *$$/d;/^ *INITLM=6 *$$/d;/^ *NADD=4 *$$/d' \
	    '$<' >'$@'

conmax/conmax.in.P conmax/conmax.in.c: conmax/conmax.in.f f2c/f2c$(EXEEXT)
	cd conmax && ../f2c/f2c$(EXEEXT) -P -a conmax.in.f

conmax/conmax.h: %.h: %.in.P
	sed -e 's/doublereal *\* *pttbl/void *pttbl/g' '$<' >'$@'

conmax/conmax.c: %.c: %.in.c
	sed -e '/int MAIN__(/,/} \/\* MAIN__ /d' \
	    -e '/^ *\(\/\*[^*]*\*\/\)* *int fnset_(/,/} \/\* fnset_ /d' \
	    -e '/extern *\(\/\*[^*]*\*\/\)* *int /,/;/d' \
	    -e 's/doublereal *\* *pttbl/void *pttbl/g' \
	    -e '/#include "f2c\.h"/{p;s/f2c/conmax/;}' '$<' >'$@'

conmax/conmax.a: conmax/conmax.c conmax/conmax.h include/f2c.h
	cd conmax && $(call cc,conmax.c,-I../include -DMSDOS)
	cd conmax && $(call ar,conmax,conmax.$(OBJEXT))

include/conmax.h: conmax/conmax.h
	mkdir -p include && cp -f '$<' '$@'

lib/conmax.a: conmax/conmax.a
	mkdir -p lib && cp -f '$<' '$@'


# Build QML
VERSION=0.3.1
CONFIG=QML_VERSION=$(VERSION)

SOURCES=qml.c include/fdlibm.h include/k.h include/f2c.h include/clapack.h
SOURCES+= include/conmax.h include/config.h
LIBS=
ifeq "$(PLATFORMCLASS)" "w"
    CONFIG+= QML_DLLEXPORT
    LIBS+= lib/q.a
endif
LIBS+=lib/conmax.a
LIBS+=lib/cephes.a lib/fdlibm.a
LIBS+=lib/clapack.a
ifeq "$(BLAS)" "atlas"
    LIBS+=lib/cblaswr.a lib/cblas.a lib/atlas.a
else
    LIBS+=lib/blas.a
endif
LIBS+=lib/f2c.a lib/fdlibm.a

build/config.c:
	mkdir -p build
	{ echo '#include <stdio.h>'; echo '#include <k.h>'; \
	    echo '#include <f2c.h>'; echo 'int main(){\
	        if(sizeof(I)>=sizeof(integer)) {\
		    puts("#define QML_kLONG(x) ((integer*)kI(x))"); \
		    puts("#define QML_KLONG KI");\
		} else if(sizeof(J)>=sizeof(integer)) {\
		    puts("#define QML_kLONG(x) ((integer*)kJ(x))"); \
		    puts("#define QML_KLONG KJ");\
		} else return 1;return 0;}'; } >'$@'

build/config$(EXEEXT): build/config.c include/k.h include/f2c.h
	cd build && $(call ccexe,config,config.c,,-I../include)

build/config.h: build/config$(EXEEXT)
	'$<' >'$@'
	$(foreach define,$(CONFIG),\
	    echo '#define $(subst =, ,$(define))' >>'$@' &&) :

include/config.h: build/config.h
	mkdir -p include && cp -f '$<' '$@'

ifeq "$(PLATFORMCLASS)" "s"
    build/compat.c:
	mkdir -p build && echo "int MAIN__() { return 0; }" >'$@'
    SOURCES+= build/compat.c
endif

build/qml.symlist: qml.c
	mkdir -p build
	sed -n -e '/#define/d;/WRAP/{' \
	    -e 's/WRAP[[:alnum:]]*(\([^,)]*\)[^)]*)/ _qml_\1 /g;H;}' \
	    -e '/qml_[_[:alnum:]]\{1,\}(/{s/^.*\(qml_[_[:alnum:]]*\).*/_\1/' \
	    -e 'H;};$${s/.*//;x;s/[[:space:]]\{1,\}/ /g;s/^ //;s/ $$//;x;G;:a' \
	    -e 's/^\(.\)\([^ ]*\) /\1\2\1/;t a' -e 's/^.//p;}' '$<' >'$@'

build/qml.mapfile: build/qml.symlist
	{ echo '{ global:'; sed -e 's/^_/    /;s/$$/;/' '$<'; \
	    echo '  local: *; };'; } >'$@'

build/qml.$(DLLEXT): $(SOURCES) $(LIBS) build/qml.mapfile build/qml.symlist
	mkdir -p build
	cd build && $(call ccdll,qml,\
	    $(patsubst %,../%,$(filter-out %.h,$(SOURCES))),\
	    $(patsubst %,../%,$(LIBS)),\
	    -I../include)

$(PLATFORM)/qml.$(DLLEXT): build/qml.$(DLLEXT)
	mkdir -p '$(PLATFORM)' && cp -pf '$<' '$@'

test: all
	q test.q -s 16


# Create distributable archive
DIST=LICENSE.txt LICENSE_CLAPACK.txt LICENSE_LIBF2C.txt LICENSE_Q.txt \
    README.txt CHANGES.txt \
    Makefile qml.c qml.q test.q

dist:
	rm -f 'qml-$(VERSION).zip'
	7z a -tzip -mx=9 'qml-$(VERSION).zip' $(DIST)


# Clean up
cleanpart:
	rm -rf cephes fdlibm f2c conmax qlib include lib build

clean: cleanpart
	rm -rf atlas clapack

distclean: clean
	rm -rf download '$(PLATFORM)'
	rm -f 'qml-$(VERSION).zip'
