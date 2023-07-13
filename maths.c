#include <Python.h>
#include "header.h"
#define PY_SSIZE_T_CLEAN
#include <numpy/arrayobject.h>


static PyObject *py_square(PyObject *self, PyObject *args) {
  int n_num, result;
  if (!PyArg_ParseTuple(args, "i", &n_num)) {
    return NULL;
  }
  result = square(n_num);

  return Py_BuildValue("i", result);
}
static PyObject *py_square_vec(PyObject *self, PyObject *args) {
  PyArrayObject* array;

  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &array)){
        return NULL;
  }

  

    int n1 = array->dimensions[0];
    printf("n1 = %d\n", n1);

    double * output = (double *) malloc(sizeof(double) * n1);
    for (int i = 0; i < n1; ++i)
    {
        output[i]=0;
    }

    for (int i = 0; i < n1; i++){
        output[i] = *((double *) array -> data + i) ;
        output[i] = output[i]*output[i];
    }

    
    return PyArray_SimpleNewFromData(1, PyArray_DIMS(array), PyArray_TYPE(array), output);
    //return Py_None;
}   

static PyObject *py_square_vec_noreturn(PyObject *self, PyObject *args) {
  PyArrayObject* array;

  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &array)){
        return NULL;
  }

  if (array -> nd != 1 || array->descr->type_num != PyArray_DOUBLE ) {
        PyErr_SetString(PyExc_ValueError, "arrays must be one-dimensional and of type float");
        return NULL;
    }

    int n1 = array->dimensions[0];

    double * arr = (double *) array -> data;

    square_vec(arr, n1);



    return Py_None;
} 

static PyMethodDef test4Methods[] = {
  {"square", py_square, METH_VARARGS, "Function for calculating square in C"},
  {"square_vec", py_square_vec, METH_VARARGS, "Function for calculating square in C"},
  {"square_vec_noreturn", py_square_vec_noreturn, METH_VARARGS, "Function for calculating square in C"},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef test4 = {
  PyModuleDef_HEAD_INIT,
  "test4",
  "Custom test4 module",
  -1,
  test4Methods
};

PyMODINIT_FUNC PyInit_test4(void)
{
    import_array();
    return PyModule_Create(&test4);
}
