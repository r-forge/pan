/* dist_dna.c       2008-12-22 */

/* Copyright 2005-2008 Emmanuel Paradis

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <R_ext/Lapack.h>

/* from R: print(log(4), d = 22) */
#define LN4 1.386294361119890572454

/* returns 8 if the base is known surely, 0 otherwise */
#define KnownBase(a) a & 8

/* returns 1 if the base is adenine surely, 0 otherwise */
#define IsAdenine(a) a == 136

/* returns 1 if the base is guanine surely, 0 otherwise */
#define IsGuanine(a) a == 72

/* returns 1 if the base is cytosine surely, 0 otherwise */
#define IsCytosine(a) a == 40

/* returns 1 if the base is thymine surely, 0 otherwise */
#define IsThymine(a) a == 24

/* returns 1 if the base is a purine surely, 0 otherwise */
#define IsPurine(a) a > 63

/* returns 1 if the base is a pyrimidine surely, 0 otherwise */
#define IsPyrimidine(a) a < 64

/* returns 1 if both bases are different surely, 0 otherwise */
#define DifferentBase(a, b) (a & b) < 16

/* returns 1 if both bases are the same surely, 0 otherwise */
#define SameBase(a, b) KnownBase(a) && a == b

#define CHECK_PAIRWISE_DELETION\
    if (KnownBase(x[s1]) && KnownBase(x[s2])) L++;\
    else continue;

#define COUNT_TS_TV\
    if (SameBase(x[s1], x[s2])) continue;\
    Nd++;\
    if (IsPurine(x[s1]) && IsPurine(x[s2])) {\
        Ns++;\
        continue;\
    }\
    if (IsPyrimidine(x[s1]) && IsPyrimidine(x[s2])) Ns++;



void Ts(unsigned char *x, int *n, int *s, double *d, int *scaled)
{
        int i1, i2, s1, s2, target, Nd, Ns;

        target = 0;
        for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
            Nd = Ns = 0;
            for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n)
              COUNT_TS_TV
                      if (*scaled) d[target] = ((double) Ns / *s);
                      else d[target] = ((double) Ns);
            target++;
        }
    }
}

void Ts_pairdel(unsigned char *x, int *n, int *s, double *d, int *scaled)
{
        int i1, i2, s1, s2, target, Nd, Ns, L;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
            Nd = Ns = L = 0;
            for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
                CHECK_PAIRWISE_DELETION
                COUNT_TS_TV
            }
            if (*scaled) d[target] = ((double) Ns/L);
            else d[target] = ((double) Ns);
            target++;
        }
    }
}

void nb_ts(unsigned char *x, int *n, int *s, double *d,
              int *pairdel, int *scaled)
{
        if(pairdel) Ts_pairdel(x, n, s, d, scaled);
        else Ts(x, n, s, d, scaled);
}
