\d .qml
{{y set{x 0^first where{count inv `$":",string[x],".dll"}each x}
    [`qml,`$"/"sv string .z.o,`qml]
    2:(`$"qml_",string y;x)}[x]'[y]}'[1 2 3;(
    `tan`asin`acos`atan`sinh`cosh`tanh`asinh`acosh`atanh`exp`expm1`log`log10,
    `logb`log1p`cbrt`floor`ceil`fabs`erf`erfc`lgamma`gamma`j0`j1`y0`y1,
    `ncdf`nicdf`mdet`minv`mevu`mchol`msvd`poly`const;
    `atan2`pow`hypot`fmod`beta`pgammar`pgammarc`ipgammarc,
    `c2cdf`c2icdf`stcdf`sticdf;
    `pbetar`ipbetar`fcdf`ficdf`gcdf`gicdf)];
version:const 0;
pi:const 1;
e:const 2;
eps:const 3;
pgamma:{gamma[x]*pgammar[x;y]};
pgammac:{gamma[x]*pgammarc[x;y]};
pbeta:{beta[x;y]*pbetar[x;y;z]};
mdiag:{{@[x#abs[type y]$0;z;:;y]}[count x]'[x;til count x]};
mpinv:{flip x[0]mmu{?[(y<x)|y=0;y;1%y]}[x[1;0;0]*eps*count[x 0]|count x 2]'[x 1]
    mmu flip(x:msvd x)2};
mev:{x@\:idesc{$[0>type x;x*x;sum x*x]}'[(x:mevu x)0]};

\d .
