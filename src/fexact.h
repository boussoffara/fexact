//
// Created by Alaa on 01/05/2021.
//

#ifndef FEXACT_FEXACT_H
#define FEXACT_FEXACT_H
void
fexact(int nrow, int ncol, const int table[], int ldtabl,
       double expect, double percnt, double emin,
       double *prt, double *pre,
       int workspace, int mult);
#endif //FEXACT_FEXACT_H
