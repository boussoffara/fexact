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
trace(PyObject *self, PyObject *args, PyObject *kwargs){
  PyArrayObject *array;
  int workspace=200000;
  int i,j;
  int mult=30;
  int hybrid=0;
  double expect=-1;
  double percnt=100;
  double emin=0;
  double prt=0;
  double pre=0;

  static char *kwlist[] = {"array", "workspace","hybrid", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!|ip", kwlist,&PyArray_Type, &array,&workspace,&hybrid))
    {
        return NULL;
    }
  if (hybrid){
    expect=5;
    percnt=80;
    emin= 1;
  }

  if (array->nd != 2 || array->descr->type_num != NPY_INT64) {
    PyErr_SetString(PyExc_ValueError,
    "array must be two-dimensional and cast of type np.int64, (np.array(x,dtype=np.int64)");
    return NULL;
  }

  // copy table and prepare corresponding vars for fexact
  int nrow = (int) array->dimensions[0];
  int ldtabl= (int) array->dimensions[0];
  int ncol= (int) array->dimensions[1];
  int * table=malloc(nrow*ncol*sizeof(int));
  for (j = 0; j < ncol; j++){
    for (i = 0; i < nrow; i++){
    table[i + j * ldtabl]=*(int *)(array->data + i*array->strides[0] + j*array->strides[1]);
    }
  }


  // call fexact
  fexact(nrow, ncol, table, ldtabl,
  expect, percnt, emin, &prt,
  &pre, workspace,
  mult);
  free(table);
  if (PyErr_Occurred()){
    return NULL;
  }

  return PyFloat_FromDouble(pre);
}



//fexact(contingency,workspace,)


static PyMethodDef fexact_methods[] = {
	{"fexact", trace,	METH_VARARGS | METH_KEYWORDS,
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
