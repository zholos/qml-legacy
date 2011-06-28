# Make configuration options:
#   make PLATFORM=l32 (default is to autodetect)
#   make BLAS=atlas   (default is blas)
#

ifneq "$(findstring $(PLATFORM),w32 w64 l32 l64 m32 m64 v32 v64 s32 s64)" ""
else ifdef PLATFORM
    $(error PLATFORM=$(PLATFORM) is not valid)
else ifeq "$(patsubst CYGWIN%,,$(shell uname -s).)" ""
    ifeq "$(patsubst %WOW64,,.$(shell uname -s))" ""
        PLATFORM=w64
    else
        PLATFORM=w32
    endif
else ifeq "$(shell uname -s)" "Linux"
    ifeq "$(shell uname -m)" "x86_64"
        PLATFORM=l64
    else
        PLATFORM=l32
    endif
else ifeq "$(shell uname -s)" "Darwin"
    ifeq "$(shell uname -m)" "x86_64"
        PLATFORM=m64
    else
        PLATFORM=m32
    endif
else ifeq "$(shell uname -s)" "SunOS"
    ifneq "$(findstring i386,$(shell isainfo))" ""
        ifeq "$(shell isainfo -b)" "64"
            PLATFORM=v64
        else
            PLATFORM=v32
        endif
    else
        ifeq "$(shell isainfo -b)" "64"
            PLATFORM=s64
        else
            PLATFORM=s32
        endif
    endif
else
    $(error couldn't determine platform, please set PLATFORM variable)
endif
$(info Platform:   $(PLATFORM))

PLATFORMCLASS=$(patsubst v,s,$(patsubst %64,%,$(patsubst %32,%,$(PLATFORM))))
PLATFORMBITS=$(patsubst %64,64,$(patsubst %32,32,$(PLATFORM)))

ifeq "$(BLAS)" "atlas"
else ifdef BLAS
    $(error BLAS="$(BLAS)" is not valid)
else
    BLAS=blas
endif
ifeq "$(ATLAS)" "blas"
    $(error use BLAS=atlas)
endif
$(info BLAS:       $(BLAS))

$(shell sleep 1)


#
# Compiler and flags
#

# CDEFS and COPTS are flags that control what code is generated.
# CFLAGS are flags necessary to make objects in the correct format.

CDEFS=-D_IEEE_LIBM -D__LITTLE_ENDIAN -D_REENTRANT
COPTS=-std=gnu99 -O2 -fno-strict-aliasing -ffloat-store \
    -Wall -Wno-parentheses -Wno-uninitialized
# -fstrict-aliasing breaks fdlibm
# -fno-float-store seems to break incbi() on Linux
CFLAGS=-pipe -fno-builtin-sin -fno-builtin-cos -fno-builtin-sqrt
SOFLAGS=
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
    SOFLAGS+= -bundle -undefined dynamic_lookup -nodefaultlibs -Wl,-x
    ENVFLAGS+= MACOSX_DEPLOYMENT_TARGET=10.4
else
    ifeq "$(PLATFORMCLASS)" "s"
        # OpenSolaris flags
        SOFLAGS+= -G -Wl,-s,-Bsymbolic
    else
        # default flags
        SOFLAGS+= -shared -Wl,-s,-Bsymbolic
    endif
endif

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
SEDI=sed -i.tmp

cc=$(CC) $(CDEFS) $(2) $(COPTS) $(CFLAGS) -c $(1)
ar=$(AR) r $(1).a $(2)
ccdll=$(ENVFLAGS) $(CC) $(CDEFS) $(4) $(COPTS) $(CFLAGS) $(SOFLAGS)
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
ccexe=$(ENVFLAGS) $(CC) $(CDEFS) $(4) $(COPTS) $(CFLAGS) -o $(1)$(EXEEXT) \
    $(patsubst %,'%',$(2)) $(patsubst %,'%',$(3))


#
# Main build target
#

.PHONY: all
all: $(PLATFORM)/qml.$(DLLEXT)


#
# Download external libraries
#

SOURCEFORGE_DIR=http://prdownloads.sourceforge.net
NETLIB_DIR=http://www.netlib.org

FDLIBM_SRC=$(SOURCEFORGE_DIR)/gnuwin32/fdlibm-5.2-2-src.zip

CPROB_SRC=$(NETLIB_DIR)/cephes/cprob.tgz

ATLAS_SRC=$(SOURCEFORGE_DIR)/math-atlas/Stable/3.8.4/atlas3.8.4.tar.bz2
CLAPACK_SRC=$(NETLIB_DIR)/clapack/clapack.tgz

F2C_BIN_GZ=$(NETLIB_DIR)/f2c/mswin/f2c.exe.gz
F2C_SRC=ftp://ftp.netlib.org/f2c/src.tar.gz
F2C_H=$(NETLIB_DIR)/f2c/f2c.h

CONMAX_F=$(NETLIB_DIR)/opt/conmax.f

Q_DIR=http://kx.com/q
K_H=$(Q_DIR)/c/c/k.h
Q_LIB=$(Q_DIR)/$(PLATFORM)/q.lib

download:
	mkdir $@

download/$(PLATFORM): | download
	mkdir $@

download/fdlibm-src.zip: | download
	wget -O $@ '$(FDLIBM_SRC)'

download/cprob.tgz: | download
	wget -O $@ '$(CPROB_SRC)'

download/k.h: | download
	wget -O $@ '$(K_H)'

download/$(PLATFORM)/q.lib: | download/$(PLATFORM)
	wget -O $@ '$(Q_LIB)'

download/atlas.tar.bz2: | download
	wget -O $@ '$(ATLAS_SRC)'

download/clapack.tgz: | download
	wget -O $@ '$(CLAPACK_SRC)'

download/f2c.exe.gz: | download
	wget -O $@ '$(F2C_BIN_GZ)'

download/f2c.tar.gz: | download
	wget -O $@ '$(F2C_SRC)'

download/f2c.h: | download
	wget -O $@ '$(F2C_H)'

download/conmax.f: | download
	wget -O $@ '$(CONMAX_F)'


#
# Build FDLIBM
#

fdlibm/.extracted: download/fdlibm-src.zip
	mkdir -p fdlibm
	unzip -o -j $< -d fdlibm
	touch $@

fdlibm/.patched: fdlibm/.extracted
	$(SEDI) $(foreach const,\
	    DOMAIN SING OVERFLOW UNDERFLOW TLOSS PLOSS __P,\
	    -e '/#define[[:space:]]*$(const)/{h;s/.*/#undef $(const)/p;g;}') \
	    fdlibm/fdlibm.h
	rm -f fdlibm/libm-dllversion.c fdlibm/w_gamma.c fdlibm/*.o
	touch $@

fdlibm/.built: fdlibm/.patched
	cd fdlibm && $(call cc,*.c)
	cd fdlibm && $(call ar,fdlibm,*.$(OBJEXT))
	$(RANLIB) fdlibm/fdlibm.a
	touch $@

include:
	mkdir $@

include/fdlibm.h: fdlibm/.patched | include
	cp fdlibm/fdlibm.h $@

lib:
	mkdir $@

lib/fdlibm.a: fdlibm/.built | lib
	cp fdlibm/fdlibm.a $@


#
# Build Cephes library
#

cephes/.extracted: download/cprob.tgz
	mkdir -p cephes
	tar xzf download/cprob.tgz -C cephes
	touch $@

cephes/.patched: cephes/.extracted
	$(SEDI) -e '/#[[:space:]]*define[[:space:]]\{1,\}UNK/d' \
	    -e '$${p;s/.*/#define IBMPC 1/;}' cephes/mconf.h
	$(SEDI) -e 's/erf/c_erf/g' \
	    -e '1{h;s/.*/extern double c_erf(double), c_erfc(double);/p;g;}' \
	    cephes/ndtr.c
	$(SEDI) -e 's/y01/y0/g;/double  *igami(/{:a' \
	    -e 's/{/&if(y0<=0||y0>=1)return(NAN);/;t b' -e 'N;b a' -e ':b' \
	    -e '};' cephes/igami.c
	$(SEDI) -e 's/\(for *( *i *= *0; *i *< *\)1\(00 *;\)/\15\2/' \
	    cephes/incbi.c
	for file in cephes/*.h; do \
	    $(SEDI) -e 's/true_gamma/gamma/g' "$$file"; done
	for file in $(filter-out %/const.c,$(wildcard cephes/*.c)); do \
	    $(SEDI) -e 's/true_gamma/gamma/g' \
	        -e 's/char st\{0,1\}\[\]="[[:alnum:]][[:alnum:]]*";//g' \
	        -e 's/mtherr *( *["sf][^,]*, *PLOSS *)//g' \
	        -e 's/mtherr *( *["sf][^,]*,[^)]*)/return(NAN)/g' \
	        -e '1{h;s/.*/extern double NAN;/p;g;}' "$$file"; done
	$(SEDI) -n -e 'H;$${x;' \
	    -e 's/( *!isfinite(x) *)/(x==INFINITY||x==-INFINITY)/g' \
	    -e 's/\(extern  *\)\{0,1\}int  *sgngam[^;]*;//g' \
	    -e 's/double  *gamma([^)]*)[^;][^{]*{/&int sgngam;/g' \
	    -e 's/double  *lgam([^)]*)[^;][^{]*{/&int sgngam;/g' \
	    -e 'p;}' cephes/gamma.c
	$(SEDI) -e '/if *( *n *<= *0[^0-9]/{N' \
	    -e 's/\(return *(\)-1\.[^)]*/\1NAN/;};/^ *kolmogorov *(/{:a' \
	    -e 's/{/&if(y<.12)return 1;/;t b' -e 'N;b a' -e ':b' -e '}' \
	    cephes/kolmogorov.c
	rm -f cephes/mtherr.c
	touch $@

cephes/.built: cephes/.patched
	cd cephes && $(call cc,*.c)
	cd cephes && $(call ar,cephes,*.$(OBJEXT))
	$(RANLIB) cephes/cephes.a
	touch $@

lib/cephes.a: cephes/.built | lib
	cp cephes/cephes.a $@


#
# Install q glue
#

include/k.h: download/k.h | include
	sed -e '/#define  *isnan/d;$${p;s/.*//;}' $< >$@

qlib/q.def: download/$(PLATFORM)/q.lib
	mkdir -p qlib
	$(NM) $< | awk "/^[[:alnum:].]+:/&&!name[1]\
	    {split(\$$1,name,\":\");print\"NAME\",name[1];print\"EXPORTS\"}\
	    /[[:space:]]T[[:space:]]+_?[[:alnum:]]/{sub(\"^_\",\"\",\$$3);\
	    print \$$3;c=1}END{exit !c}" >$@.tmp
	mv $@.tmp $@

qlib/q.a: qlib/q.def
	$(DLLTOOL) -d $< -l $@

lib/q.a: qlib/q.a | lib
	cp $< $@


#
# Build LAPACK libraries
#

# ATLAS is an optimized BLAS plus some LAPACK routines
atlas/.extracted: download/atlas.tar.bz2
	mkdir -p atlas
	tar xjf $< -C atlas
	touch $@

atlas/.patched: atlas/.extracted
	$(SEDI) -f patch/atlas-ilaenv.sed \
	    atlas/ATLAS/interfaces/lapack/F77/src/ilaenv.f
	$(SEDI) 's/- \(cd .* \/F77\)/\1/' atlas/ATLAS/makes/Make.bin
	touch $@

atlas/.configured: atlas/.patched
	mkdir -p atlas/build
	cd atlas/build && ../ATLAS/configure -t 0 -b $(PLATFORMBITS) \
	    -C ic $(CC) -C if $(FC) -Fa al '$(CFLAGS)' --prefix=../install
	$(SEDI) '/FLAGS *=/s/-O\($$\|[1 ]\)/-O2 /' atlas/build/Make.inc
	touch $@

atlas/.built: atlas/.configured
	make -C atlas/build build install
	touch $@

lib/atlas: | lib
	mkdir $@

lib/atlas/atlas.a: atlas/.built | lib/atlas
	cp atlas/install/lib/libatlas.a $@

lib/atlas/cblas.a: atlas/.built | lib/atlas
	cp atlas/install/lib/libcblas.a $@

lib/atlas/lapack.a: atlas/.built | lib/atlas
	cp atlas/install/lib/liblapack.a $@

# CLAPACK includes libf2c, reference BLAS and LAPACK translated to C
clapack/.extracted: download/clapack.tgz
	mkdir -p clapack
	tar xzf $< -C clapack
	mv clapack/CLAPACK-*/* clapack/
	touch $@

clapack/.patched: clapack/.extracted
	# thread-safety fix:
	$(SEDI) '/equiv_/s/static//g' clapack/SRC/slaln2.c clapack/SRC/dlaln2.c
	
	$(SEDI) -e 's/^\(TIME[[:space:]]*=\).*/\1/' \
	    -e 's/\([^_[:alnum:]]\)cc\([^_[:alnum:]]\)/\1$$(CC)\2/' \
	    -e 's/\([^_[:alnum:]]\)ld\([^_[:alnum:]]\)/\1$$(LD)\2/' \
	    -e 's/\([^_[:alnum:]]\)ar\([^_[:alnum:]]\)/\1$$(ARCH)\2/' \
	    -e 's/\([^_[:alnum:]]\)ranlib\([^_[:alnum:]]\)/\1$$(RANLIB)\2/' \
	    clapack/F2CLIBS/libf2c/Makefile
      ifeq "$(PLATFORMCLASS)" "w"
	$(SEDI) -e 's/\([^_[:alnum:]]\)a\.out/\1a.exe/g' \
	    clapack/F2CLIBS/libf2c/Makefile
      endif
      ifneq "$(findstring $(PLATFORMCLASS),m s)" ""
	$(SEDI) -e '/\$$(LD) .*\.xxx/d;/mv .*\.xxx/d' \
	    clapack/F2CLIBS/libf2c/Makefile
	$(SEDI) 's,\./test[a-z]*,:,g' clapack/Makefile
	$(SEDI) '/-o  *test[a-z]* /d' clapack/INSTALL/Makefile
      endif
	$(SEDI) -e '/^#define MSpc\([^_[:alnum:]]\|$$\)/d' \
	    clapack/F2CLIBS/libf2c/uninit.c
	{ echo '#include "f2c.h"'; echo 'doublereal second_(){return 0;}'; } \
	    >clapack/INSTALL/second.c
	sed -e 's/second/dsecnd/g' clapack/INSTALL/second.c \
	    >clapack/INSTALL/dsecnd.c
      ifeq "$(PLATFORMCLASS)" "w"
	find clapack -name '*.h' -type l -exec cp '{}' clapack/.tmp.h ';' \
	    -exec rm -- '{}' ';' -exec mv clapack/.tmp.h '{}' ';'
      endif
	touch $@

ifeq "$(BLAS)" "atlas"
    clapack/make.inc.cflags: atlas/.configured
	sed -n '/ICCFLAGS *=/s//CFLAGS=/p' atlas/build/Make.inc > $@.tmp
	mv $@.tmp $@
else
    clapack/make.inc.cflags:
	echo 'CFLAGS=$$(CDEFS) $(COPTS) $(CFLAGS)' > $@
endif

clapack/.configured: clapack/.patched clapack/make.inc.cflags
	cp clapack/make.inc.example clapack/make.inc
	cat clapack/make.inc.cflags >>clapack/make.inc
	{ echo 'PLAT=_$(PLATFORM)'; echo 'CC=$(CC)'; \
	    echo 'CDEFS=-I$$(TOPDIR)/INCLUDE'; \
	    echo 'NOOPT=$$(CFLAGS) -O0'; echo 'LOADOPTS=$$(CFLAGS)'; \
	    echo 'LOADER=$$(CC)'; echo 'LD=$(LD)'; echo 'ARCH=$(AR)'; \
	    echo 'RANLIB=$(RANLIB)'; } >>clapack/make.inc
	touch $@

clapack/.built: clapack/.configured
	make -C clapack/F2CLIBS/libf2c hadd
	make -C clapack f2clib
      ifeq "$(BLAS)" "blas"
	make -C clapack blaslib
      else
	make -C clapack cblaswrap
      endif
	make -C clapack lapacklib
	touch $@

include/f2c.h: clapack/.built | include
	cp clapack/INCLUDE/f2c.h $@

lib/f2c.a: clapack/.built | lib
	cp clapack/F2CLIBS/libf2c.a $@

include/clapack.h: clapack/.built
	cp clapack/INCLUDE/clapack.h $@

lib/lapack: | lib
	mkdir $@

lib/lapack/blas.a: clapack/.built | lib/lapack
	cp clapack/blas_$(PLATFORM).a $@

lib/lapack/cblaswr.a: clapack/.built | lib/lapack
	cp clapack/libcblaswr.a $@

lib/lapack/lapack.a: clapack/.built | lib/lapack
	cp clapack/lapack_$(PLATFORM).a $@

# Combine lapack.a from CLAPACK with optimized LAPACK routines from ATLAS
alapack/.built: lib/atlas/lapack.a lib/lapack/lapack.a
	mkdir -p alapack
	cp lib/lapack/lapack.a alapack/alapack.a
	cd alapack && $(AR) x ../lib/atlas/lapack.a
	cd alapack && $(AR) r alapack.a *.o
	$(RANLIB) alapack/alapack.a
	touch $@

lib/alapack.a: alapack/.built | lib
	cp alapack/alapack.a $@


# Build f2c program
ifeq "$(PLATFORMCLASS)" "w"
    f2c/.built: download/f2c.exe.gz
	mkdir -p f2c
	gunzip -c $< >f2c/f2c$(EXEEXT)
	chmod +x f2c/f2c$(EXEEXT)
	touch $@
else
    f2c/.extracted: download/f2c.tar.gz
	mkdir -p f2c
	tar xzf $< -C f2c
	mv f2c/src/* f2c/
	touch $@

    f2c/.patched: f2c/.extracted
	{ echo 'CC=$(CC)'; echo 'CFLAGS=$(CDEFS) $(COPTS) $(CFLAGS)'; \
	    sed -e '/^ *CC *=/d;/^ *CFLAGS *=/d;/^f2c:/{n;' \
	        -e 's/\$$(CC) /&$$(CFLAGS) /;}' f2c/makefile.u; } >f2c/Makefile
	$(SEDI) -e '/Cextern int \(unlink\|fork\)/d' f2c/sysdep.c
	touch $@

    f2c/.built: f2c/.patched
	make -C f2c f2c
	touch $@
endif

bin:
	mkdir $@

bin/f2c$(EXEEXT): f2c/.built | bin
	cp f2c/f2c$(EXEEXT) $@


# Build conmax
conmax/conmax.in.f: download/conmax.f
	mkdir -p conmax
	sed -e '/^C/d' \
	    -e 's/PTTBL(IPTB,INDM)/PTTBL(0:0,0:0)/g;s/IPTB,//g;s/INDM,//g' \
	    -e '/NE MULLER(/{h;s/(.*/(LIMMUL,NSRCH,/p;g;s/.*(/     */;}' \
	    -e '/LL MULLER(/{h;s/(.*/(5,     NSRCH,/p;g;s/.*(/     */;}' \
	    -e '/NE SEARSL(/{h;s/(.*/(INITLM,NADD,LIMS1,/p;g;s/.*(/     */;}' \
	    -e '/LL SEARSL(/{h;s/(.*/(6,     4,   LIMS1,/p;g;s/.*(/     */;}' \
	    -e '/^ *LIMMUL=5 *$$/d;/^ *INITLM=6 *$$/d;/^ *NADD=4 *$$/d' \
	    $< >$@

conmax/conmax.in.P conmax/conmax.in.c: conmax/conmax.in.f bin/f2c$(EXEEXT)
	cd conmax && ../bin/f2c$(EXEEXT) -P -a conmax.in.f

conmax/conmax.h: %.h: %.in.P
	sed -e 's/doublereal *\* *pttbl/void *pttbl/g' $< >$@

conmax/conmax.c: %.c: %.in.c
	sed -e '/int MAIN__(/,/} \/\* MAIN__ /d' \
	    -e '/^ *\(\/\*[^*]*\*\/\)* *int fnset_(/,/} \/\* fnset_ /d' \
	    -e '/extern *\(\/\*[^*]*\*\/\)* *int /,/;/d' \
	    -e 's/doublereal *\* *pttbl/void *pttbl/g' \
	    -e '/#include "f2c\.h"/{p;s/f2c/conmax/;}' $< >$@

conmax/conmax.a: conmax/conmax.c conmax/conmax.h include/f2c.h
	cd conmax && $(call cc,conmax.c,-I../include -DMSDOS)
	cd conmax && $(call ar,conmax,conmax.$(OBJEXT))

include/conmax.h: conmax/conmax.h | include
	cp $< $@

lib/conmax.a: conmax/conmax.a | lib
	cp $< $@


#
# Build QML
#

VERSION=0.3.2
CONFIG=QML_VERSION=$(VERSION)$(patsubst -,,-$(patsubst blas,,$(BLAS)))

SOURCES=qml.c include/fdlibm.h include/k.h include/f2c.h include/clapack.h \
    include/conmax.h include/config.h
LIBS=
ifeq "$(PLATFORMCLASS)" "w"
    CONFIG+= QML_DLLEXPORT
    LIBS+= lib/q.a
endif
LIBS+= lib/conmax.a
LIBS+= lib/cephes.a lib/fdlibm.a
ifeq "$(BLAS)" "atlas"
    LIBS+= lib/alapack.a
    LIBS+= lib/lapack/cblaswr.a lib/atlas/cblas.a lib/atlas/atlas.a
else
    LIBS+= lib/lapack/lapack.a lib/lapack/blas.a
endif
LIBS+= lib/f2c.a lib/fdlibm.a

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
		} else return 1;return 0;}'; } >$@

build/config$(EXEEXT): build/config.c include/k.h include/f2c.h
	cd build && $(call ccexe,config,config.c,,-I../include)

build/config.h: build/config$(EXEEXT)
	$< >$@
	$(foreach define,$(CONFIG),\
	    echo '#define $(subst =, ,$(define))' >>$@ &&) :

include/config.h: build/config.h | include
	cp $< $@

ifeq "$(PLATFORMCLASS)" "s"
    build/compat.c:
	mkdir -p build
	echo "int MAIN__() { return 0; }" >$@
    SOURCES+= build/compat.c
endif

build/qml.symlist: qml.c
	mkdir -p build
	sed -n -e '/#define/d;/WRAP/{' \
	    -e 's/WRAP[[:alnum:]]*(\([^,)]*\)[^)]*)/ _qml_\1 /g;H;}' \
	    -e '/qml_[_[:alnum:]]\{1,\}(/{s/^.*\(qml_[_[:alnum:]]*\).*/_\1/' \
	    -e 'H;};$${s/.*//;x;s/[[:space:]]\{1,\}/ /g;s/^ //;s/ $$//;x;G;:a' \
	    -e 's/^\(.\)\([^ ]*\) /\1\2\1/;t a' -e 's/^.//p;}' $< >$@

build/qml.mapfile: build/qml.symlist
	{ echo '{ global:'; sed -e 's/^_/    /;s/$$/;/' $<; \
	    echo '  local: *; };'; } >$@

build/qml.$(DLLEXT): $(SOURCES) $(LIBS) build/qml.mapfile build/qml.symlist
	mkdir -p build
	cd build && $(call ccdll,qml,\
	    $(patsubst %,../%,$(filter-out %.h,$(SOURCES))),\
	    $(patsubst %,../%,$(LIBS)),\
	    -I../include)

$(PLATFORM)/qml.$(DLLEXT): build/qml.$(DLLEXT)
	mkdir -p $(PLATFORM)
	cp -p $< $@

.PHONY: test
test: all
	q test.q -s 16


#
# Create distributable archive
#

DIST=LICENSE.txt LICENSE_CLAPACK.txt LICENSE_LIBF2C.txt LICENSE_Q.txt \
    README.txt CHANGES.txt \
    Makefile patch qml.c qml.q test.q

.PHONY: dist
dist:
	rm -f 'qml-$(VERSION).zip'
	7z a -tzip -mx=9 'qml-$(VERSION).zip' $(DIST)


#
# Clean up
#

.PHONY: cleanpart clean distclean
cleanpart:
	rm -rf cephes fdlibm alapack f2c conmax qlib include lib bin build

# These two take much longer to build, sometimes we want to keep them.
clean: cleanpart
	rm -rf atlas clapack

distclean: clean
	rm -rf download '$(PLATFORM)'
	rm -f 'qml-$(VERSION).zip'
