# In relevant LAPACK routines, replace particular string concatenation with
# array operations, so that calls to _gfortran_concat_string aren't emitted.

/^  *CHARACTER *SIDE, *TRANS/s/$/, SIDETRANS(2)/
/^\(      .*\)SIDE *\/\/ *TRANS/{s//\1SIDETRANS/;h;
    s/[^ ].*/SIDETRANS(1) = SIDE/p
    s/[^ ].*/SIDETRANS(2) = TRANS/p
g;}

/^  *CHARACTER *JOBU, *JOBVT/s/$/, JOBUJOBVT(2)/
/^\(      .*\)JOBU *\/\/ *JOBVT/{s//\1JOBUJOBVT/;h;
    s/[^ ].*/JOBUJOBVT(1) = JOBU/p
    s/[^ ].*/JOBUJOBVT(2) = JOBVT/p
g;}

/^  *CHARACTER *COMPZ, *JOB/s/$/, JOBCOMPZ(2)/
/^\(      .*\)JOB( *: *1 *) *\/\/ *COMPZ( *: *1 *)/{s//\1JOBCOMPZ/;h;
    s/[^ ].*/JOBCOMPZ(1) = JOB/p
    s/[^ ].*/JOBCOMPZ(2) = COMPZ/p
g;}


# Similarly, replace 2**X by equivalent bit shifts, so that calls to
# _gfortran_pow_i4_i4 aren't emitted.

s/ 2 *\*\* *\([[:alnum:]]*\|([^)]*)\) *$/ ISHFT(1,\1)/
