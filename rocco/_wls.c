#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

#include "native/wls_backend.h"

static PyObject *score_centered_wls(
    PyObject *self,
    PyObject *args,
    PyObject *kwargs)
{
    static char *kwlist[] = {
        "centered_matrix",
        "lower_bound_z",
        "prior_df",
        "min_effect",
        "spatial_window",
        "precision_floor_ratio",
        NULL,
    };
    PyObject *matrix_obj = NULL;
    PyArrayObject *matrix_arr = NULL;
    PyArrayObject *mean_arr = NULL;
    PyArrayObject *raw_variance_arr = NULL;
    PyArrayObject *prior_variance_arr = NULL;
    PyArrayObject *moderated_variance_arr = NULL;
    PyArrayObject *se_arr = NULL;
    PyArrayObject *scores_arr = NULL;
    double lower_bound_z = 1.0;
    double prior_df = 5.0;
    PyObject *min_effect_obj = Py_None;
    double min_effect = 0.0;
    int use_min_effect = 0;
    int spatial_window = 31;
    double precision_floor_ratio = 0.01;
    double total_df = 0.0;
    int resolved_window = 0;
    int status = 0;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "O|ddOid",
            kwlist,
            &matrix_obj,
            &lower_bound_z,
            &prior_df,
            &min_effect_obj,
            &spatial_window,
            &precision_floor_ratio))
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
    if (min_effect_obj != NULL && min_effect_obj != Py_None)
    {
        min_effect = PyFloat_AsDouble(min_effect_obj);
        if (PyErr_Occurred() != NULL)
        {
            Py_DECREF(matrix_arr);
            return NULL;
        }
        if (min_effect < 0.0)
        {
            min_effect = 0.0;
        }
        use_min_effect = 1;
    }

    {
        npy_intp sample_count = PyArray_DIM(matrix_arr, 0);
        npy_intp locus_count = PyArray_DIM(matrix_arr, 1);
        npy_intp locus_dims[1] = {locus_count};

        mean_arr = (PyArrayObject *)PyArray_ZEROS(1, locus_dims, NPY_FLOAT64, 0);
        raw_variance_arr = (PyArrayObject *)PyArray_ZEROS(1, locus_dims, NPY_FLOAT64, 0);
        prior_variance_arr = (PyArrayObject *)PyArray_ZEROS(1, locus_dims, NPY_FLOAT64, 0);
        moderated_variance_arr = (PyArrayObject *)PyArray_ZEROS(1, locus_dims, NPY_FLOAT64, 0);
        se_arr = (PyArrayObject *)PyArray_ZEROS(1, locus_dims, NPY_FLOAT64, 0);
        scores_arr = (PyArrayObject *)PyArray_ZEROS(1, locus_dims, NPY_FLOAT64, 0);
        if (mean_arr == NULL || raw_variance_arr == NULL || prior_variance_arr == NULL ||
            moderated_variance_arr == NULL || se_arr == NULL || scores_arr == NULL)
        {
            Py_DECREF(matrix_arr);
            Py_XDECREF(mean_arr);
            Py_XDECREF(raw_variance_arr);
            Py_XDECREF(prior_variance_arr);
            Py_XDECREF(moderated_variance_arr);
            Py_XDECREF(se_arr);
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
                min_effect,
                use_min_effect,
                spatial_window,
                precision_floor_ratio,
                (double *)PyArray_DATA(mean_arr),
                (double *)PyArray_DATA(raw_variance_arr),
                (double *)PyArray_DATA(prior_variance_arr),
                (double *)PyArray_DATA(moderated_variance_arr),
                (double *)PyArray_DATA(se_arr),
                (double *)PyArray_DATA(scores_arr),
                &total_df,
                &resolved_window);
        Py_END_ALLOW_THREADS
    }

    Py_DECREF(matrix_arr);

    if (status != 0)
    {
        Py_DECREF(mean_arr);
        Py_DECREF(raw_variance_arr);
        Py_DECREF(prior_variance_arr);
        Py_DECREF(moderated_variance_arr);
        Py_DECREF(se_arr);
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
        "NNNNNNdi",
        scores_arr,
        mean_arr,
        raw_variance_arr,
        prior_variance_arr,
        moderated_variance_arr,
        se_arr,
        total_df,
        resolved_window);
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
