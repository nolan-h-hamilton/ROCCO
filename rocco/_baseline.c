#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

#include "native/baseline_backend.h"

/*
 * Python wrapper for the Whittaker baseline code.
 *
 * The numerical core is adapted from the Consenrich baseline implementation:
 *   https://github.com/nolan-h-hamilton/Consenrich
 *   https://github.com/nolan-h-hamilton/Consenrich/blob/main/src/consenrich/cconsenrich.pyx
 */

static PyObject *crossfit_whittaker_baseline(
    PyObject *self,
    PyObject *args,
    PyObject *kwargs)
{
    static char *kwlist[] = {
        "values",
        "penalty_lambda",
        NULL,
    };
    PyObject *values_obj = NULL;
    PyArrayObject *values_arr = NULL;
    PyArrayObject *baseline_arr = NULL;
    double penalty_lambda = 0.0;
    int ndim = 0;
    int status = 0;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "Od",
            kwlist,
            &values_obj,
            &penalty_lambda))
    {
        return NULL;
    }

    values_arr = (PyArrayObject *)PyArray_FROM_OTF(
        values_obj,
        NPY_FLOAT64,
        NPY_ARRAY_IN_ARRAY);
    if (values_arr == NULL)
    {
        return NULL;
    }

    ndim = PyArray_NDIM(values_arr);
    if (ndim != 1 && ndim != 2)
    {
        PyErr_SetString(
            PyExc_ValueError,
            "`values` must be one-dimensional or two-dimensional");
        Py_DECREF(values_arr);
        return NULL;
    }

    baseline_arr = (PyArrayObject *)PyArray_ZEROS(
        ndim,
        PyArray_DIMS(values_arr),
        NPY_FLOAT64,
        0);
    if (baseline_arr == NULL)
    {
        Py_DECREF(values_arr);
        return NULL;
    }

    Py_BEGIN_ALLOW_THREADS
    if (ndim == 1)
    {
        status = rocco_crossfit_whittaker_baseline_f64(
            (const double *)PyArray_DATA(values_arr),
            (size_t)PyArray_DIM(values_arr, 0),
            penalty_lambda,
            (double *)PyArray_DATA(baseline_arr));
    }
    else
    {
        status = rocco_crossfit_whittaker_baseline_matrix_f64(
            (const double *)PyArray_DATA(values_arr),
            (size_t)PyArray_DIM(values_arr, 0),
            (size_t)PyArray_DIM(values_arr, 1),
            penalty_lambda,
            (double *)PyArray_DATA(baseline_arr));
    }
    Py_END_ALLOW_THREADS

    Py_DECREF(values_arr);

    if (status != 0)
    {
        Py_DECREF(baseline_arr);
        PyErr_NoMemory();
        return NULL;
    }

    return (PyObject *)baseline_arr;
}

static PyMethodDef baseline_methods[] = {
    {
        "crossfit_whittaker_baseline",
        (PyCFunction)crossfit_whittaker_baseline,
        METH_VARARGS | METH_KEYWORDS,
        "Estimate a Consenrich-style cross-fit Whittaker baseline."
    },
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef baseline_module = {
    PyModuleDef_HEAD_INIT,
    "_baseline",
    "Local baseline estimation for ROCCO.",
    -1,
    baseline_methods,
};

PyMODINIT_FUNC PyInit__baseline(void)
{
    PyObject *module = PyModule_Create(&baseline_module);
    if (module == NULL)
    {
        return NULL;
    }
    import_array();
    return module;
}
