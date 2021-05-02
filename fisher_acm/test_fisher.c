#include <stdio.h>
#include "fexact.h"

int main() {
    printf("Hello, World!\n");
    int nrow =3;
    int ncol=2;
    int table[6]={5,0,0,5,5,5};
    int ldtabl=3;
    double expect=-1;
    double percnt=100;
    double emin=0;
    double prt;
    double pre;
    int workspace=100000;
    int mult=30;
    fexact(&nrow, &ncol, (int*)&table, &ldtabl,
    &expect, &percnt, &emin, &prt,
    &pre, /* new in C : */ &workspace,
    /* new arg, was const = 30*/&mult);
    printf("%f\n" , prt);
    printf("%f\n",pre);
    return 0;
}
