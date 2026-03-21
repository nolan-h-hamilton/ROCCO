#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

#include "native/wls_backend.h"

/*
 * Python wrapper for centered-WLS scoring.
 *
 * This computes :math:`\hat\mu_j`, :math:`\hat{\mathrm{se}}_j`, and the
 * lower-bound score used during budget initialization.
 */

static PyObject *score_centered_wls(
    PyObject *self,
    PyObject *args,
    PyObject *kwargs)
{
    static char *kwlist[] = {
        "centered_matrix",
        "lower_bound_z",
        "prior_df",
        NULL,
    };
    PyObject *matrix_obj = NULL;
    PyArrayObject *matrix_arr = NULL;
    PyArrayObject *weights_arr = NULL;
    PyArrayObject *mean_arr = NULL;
    PyArrayObject *se_arr = NULL;
    PyArrayObject *z_arr = NULL;
    PyArrayObject *scores_arr = NULL;
    double lower_bound_z = 1.0;
    double prior_df = 5.0;
    int status = 0;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "O|dd",
            kwlist,
            &matrix_obj,
            &lower_bound_z,
            &prior_df))
    {
        return NULL;
    }

    matrix_arr = (PyArrayObject *)PyArray_FROM_OTF(
        matrix_obj,
        NPY_FLOAT64,
        NPY_ARRAY_IN_ARRAY);
    if (matrix_arr == NULL)
    {
        return NULL;
    }
    if (PyArray_NDIM(matrix_arr) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "`centered_matrix` must be two-dimensional");
        Py_DECREF(matrix_arr);
        return NULL;
    }

    {
        npy_intp sample_count = PyArray_DIM(matrix_arr, 0);
        npy_intp locus_count = PyArray_DIM(matrix_arr, 1);
        npy_intp sample_dims[1] = {sample_count};
        npy_intp locus_dims[1] = {locus_count};

        weights_arr = (PyArrayObject *)PyArray_ZEROS(1, sample_dims, NPY_FLOAT64, 0);
        mean_arr = (PyArrayObject *)PyArray_ZEROS(1, locus_dims, NPY_FLOAT64, 0);
        se_arr = (PyArrayObject *)PyArray_ZEROS(1, locus_dims, NPY_FLOAT64, 0);
        z_arr = (PyArrayObject *)PyArray_ZEROS(1, locus_dims, NPY_FLOAT64, 0);
        scores_arr = (PyArrayObject *)PyArray_ZEROS(1, locus_dims, NPY_FLOAT64, 0);
        if (weights_arr == NULL || mean_arr == NULL || se_arr == NULL ||
            z_arr == NULL || scores_arr == NULL)
        {
            Py_DECREF(matrix_arr);
            Py_XDECREF(weights_arr);
            Py_XDECREF(mean_arr);
            Py_XDECREF(se_arr);
            Py_XDECREF(z_arr);
            Py_XDECREF(scores_arr);
            return NULL;
        }

        Py_BEGIN_ALLOW_THREADS
        status = rocco_score_centered_wls_f64(
            (const double *)PyArray_DATA(matrix_arr),
            (size_t)sample_count,
            (size_t)locus_count,
            lower_bound_z,
            prior_df,
            (double *)PyArray_DATA(weights_arr),
            (double *)PyArray_DATA(mean_arr),
            (double *)PyArray_DATA(se_arr),
            (double *)PyArray_DATA(z_arr),
            (double *)PyArray_DATA(scores_arr));
        Py_END_ALLOW_THREADS
    }

    Py_DECREF(matrix_arr);

    if (status != 0)
    {
        Py_DECREF(weights_arr);
        Py_DECREF(mean_arr);
        Py_DECREF(se_arr);
        Py_DECREF(z_arr);
        Py_DECREF(scores_arr);
        if (status == -1)
        {
            PyErr_NoMemory();
        }
        else
        {
            PyErr_SetString(PyExc_ValueError, "Invalid centered-WLS inputs");
        }
        return NULL;
    }

    return Py_BuildValue(
        "NNNNN",
        scores_arr,
        weights_arr,
        mean_arr,
        se_arr,
        z_arr);
}

static PyMethodDef wls_methods[] = {
    {
        "score_centered_wls",
        (PyCFunction)score_centered_wls,
        METH_VARARGS | METH_KEYWORDS,
        "Score a centered sample-by-locus matrix.",
    },
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef wls_module = {
    PyModuleDef_HEAD_INIT,
    "_wls",
    "Centered-WLS scoring for ROCCO.",
    -1,
    wls_methods,
};

PyMODINIT_FUNC PyInit__wls(void)
{
    PyObject *module = PyModule_Create(&wls_module);
    if (module == NULL)
    {
        return NULL;
    }
    import_array();
    return module;
}
