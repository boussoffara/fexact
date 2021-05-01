//
// Created by Alaa on 01/05/2021.
//

#ifndef FEXACT_FEXACT_H
#define FEXACT_FEXACT_H
void
fexact(int *nrow, int *ncol, int *table, int *ldtabl,
       double *expect, double *percnt, double *emin, double *prt,
       double *pre, /* new in C : */ int *workspace,
        /* new arg, was const = 30*/int *mult);
#endif //FEXACT_FEXACT_H
