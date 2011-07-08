# In ATLAS ilaenv.f, replace string operations with character operations, so
# that calls to _gfortran_string_compare, memcpy and memmove aren't emitted.
# Beware the 72-character line limit.

# NAME:    123456
# SUBNAM = ******
# C1     = *
# C2     =  **
# C3     =    ***
# C4     =     **

/^  *CHARACTER.* C[1-4]/d
/^  *CHARACTER.* SUBNAM/d
/^  *I[CZ] *=/d
/^  *C[1-4] *=/d
/^  *SUBNAM *=/d
/^  *IF *( *IZ/,/^      END IF/d

s/^/     $   @/;h;s/.*//;G
s/C1 *.EQ. *'\([^']\)'/@(LSAME(NAME(1:1),'\1'))@/g
s/C2 *.EQ. *'\([^']\)\([^']\)'/@(LSAME(NAME(2:2),'\1').AND.LSAME(NAME(3:3),'\2'))@/g
s/C3 *.EQ. *'\([^']\)\([^']\)\([^']\)'/@(LSAME(NAME(4:4),'\1').AND.LSAME(NAME(5:5),'\2')@.AND.LSAME(NAME(6:6),'\3'))@/g
s/C4 *.EQ. *'\([^']\)\([^']\)'/@(LSAME(NAME(5:5),'\1').AND.LSAME(NAME(6:6),'\2'))@/g
s/LSAME( *C1 *,/LSAME(NAME(1:1),/g
s/C3 *( *1 *: *1 *) *.EQ. *'\([^']\)'/@(LSAME(NAME(4:4),'\1'))@/g
:a
s/^\([^@]*\)@\([^@]*\)@/\1@\2\1/;t a
s/^[^@]*@//
