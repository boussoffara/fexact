//
// Created by Alaa on 02/05/2021.
//

#include "fexact.h"
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdbool.h>

#if PY_MAJOR_VERSION >= 3
#define PY3K
#endif


static PyObject *
trace(PyObject *self, PyObject *args){
  PyArrayObject *array;
//  double sum;
  int i,j;
  if (!PyArg_ParseTuple(args, "O!",&PyArray_Type, &array)){
    PyErr_SetString(PyExc_ValueError,"Error while parsing the trajectory coordinates in get_spinangle_traj");
    return NULL;
  }
  if (array->nd != 2 || array->descr->type_num != PyArray_INT) {
    PyErr_SetString(PyExc_ValueError,
    "array must be two-dimensional and of type integer");
    return NULL;
  }

    int nrow =array->dimensions[0];
    int ldtabl=array->dimensions[0];
    int ncol=array->dimensions[1];

    int * table=malloc(nrow*ncol*sizeof(int));

    for (j = 0; j < ncol; j++){
      for (i = 0; i < nrow; i++){
      table[i + j * ldtabl]=*(int *)(array->data + i*array->strides[0] + j*array->strides[1]);
      }
    }


    double expect=-1;
    double percnt=100;
    double emin=0;
    double prt;
    double pre;
    int workspace=200000;
    int mult=30;
    fexact(&nrow, &ncol, table, &ldtabl,
    &expect, &percnt, &emin, &prt,
    &pre, &workspace,
    &mult);

    return PyFloat_FromDouble(pre);
}



//fexact(contingency,workspace,)


static PyMethodDef fexact_methods[] = {
	{"fexact", trace,	METH_VARARGS,
	  "Doc string."},
	{NULL,		NULL}		/* sentinel */
};


#ifdef PY3K

static struct PyModuleDef fexactmodule = {
    PyModuleDef_HEAD_INIT,
    "fexact",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,
    fexact_methods
};


PyMODINIT_FUNC PyInit_fexact(void)
{   import_array();
    return PyModule_Create(&fexactmodule);
}

#else

PyMODINIT_FUNC initfexact() {
    Py_InitModule3("fexact", fexact_methods, "mod doc");
}
#endif
