#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

#include <stdint.h>
#include <stdlib.h>

static PyObject *solve_penalized_chain(PyObject *self, PyObject *args)
{
    PyObject *scores_obj = NULL;
    PyObject *switch_costs_obj = NULL;
    PyArrayObject *scores_arr = NULL;
    PyArrayObject *switch_costs_arr = NULL;
    PyArrayObject *solution_arr = NULL;
    double *scores = NULL;
    double *switch_costs = NULL;
    npy_uint8 *solution = NULL;
    uint8_t *bt0 = NULL;
    uint8_t *bt1 = NULL;
    double selection_penalty = 0.0;
    npy_intp n = 0;
    npy_intp i = 0;
    double prev0_val = 0.0;
    double prev1_val = 0.0;
    long prev0_count = 0;
    long prev1_count = 0;
    double best_val = 0.0;
    long best_count = 0;
    int state = 0;

    if (!PyArg_ParseTuple(
            args,
            "OOd",
            &scores_obj,
            &switch_costs_obj,
            &selection_penalty))
    {
        return NULL;
    }

    scores_arr = (PyArrayObject *)PyArray_FROM_OTF(
        scores_obj,
        NPY_FLOAT64,
        NPY_ARRAY_IN_ARRAY);
    if (scores_arr == NULL)
    {
        return NULL;
    }
    switch_costs_arr = (PyArrayObject *)PyArray_FROM_OTF(
        switch_costs_obj,
        NPY_FLOAT64,
        NPY_ARRAY_IN_ARRAY);
    if (switch_costs_arr == NULL)
    {
        Py_DECREF(scores_arr);
        return NULL;
    }

    if (PyArray_NDIM(scores_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "`scores` must be one-dimensional");
        goto fail;
    }
    if (PyArray_NDIM(switch_costs_arr) != 1)
    {
        PyErr_SetString(
            PyExc_ValueError,
            "`switch_costs` must be one-dimensional");
        goto fail;
    }

    n = PyArray_DIM(scores_arr, 0);
    if (n <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "`scores` cannot be empty");
        goto fail;
    }
    if (n > 1 && PyArray_DIM(switch_costs_arr, 0) != (n - 1))
    {
        PyErr_SetString(
            PyExc_ValueError,
            "`switch_costs` must have length len(scores) - 1");
        goto fail;
    }

    scores = (double *)PyArray_DATA(scores_arr);
    switch_costs = (double *)PyArray_DATA(switch_costs_arr);

    bt0 = (uint8_t *)calloc((size_t)n, sizeof(uint8_t));
    bt1 = (uint8_t *)calloc((size_t)n, sizeof(uint8_t));
    if (bt0 == NULL || bt1 == NULL)
    {
        PyErr_NoMemory();
        goto fail;
    }

    solution_arr = (PyArrayObject *)PyArray_ZEROS(
        1,
        &n,
        NPY_UINT8,
        0);
    if (solution_arr == NULL)
    {
        goto fail;
    }
    solution = (npy_uint8 *)PyArray_DATA(solution_arr);

    prev0_val = 0.0;
    prev0_count = 0;
    prev1_val = scores[0] - selection_penalty;
    prev1_count = 1;

    Py_BEGIN_ALLOW_THREADS
    for (i = 1; i < n; ++i)
    {
        const double switch_cost = switch_costs[i - 1];
        const double stay0_val = prev0_val;
        const long stay0_count = prev0_count;
        const double switch0_val = prev1_val - switch_cost;
        const long switch0_count = prev1_count;
        double new0_val = 0.0;
        long new0_count = 0;

        const double stay1_val = prev1_val + scores[i] - selection_penalty;
        const long stay1_count = prev1_count + 1;
        const double switch1_val =
            prev0_val - switch_cost + scores[i] - selection_penalty;
        const long switch1_count = prev0_count + 1;
        double new1_val = 0.0;
        long new1_count = 0;

        if (switch0_val > stay0_val ||
            (switch0_val == stay0_val && switch0_count < stay0_count))
        {
            new0_val = switch0_val;
            new0_count = switch0_count;
            bt0[i] = 1U;
        }
        else
        {
            new0_val = stay0_val;
            new0_count = stay0_count;
            bt0[i] = 0U;
        }

        if (switch1_val > stay1_val ||
            (switch1_val == stay1_val && switch1_count < stay1_count))
        {
            new1_val = switch1_val;
            new1_count = switch1_count;
            bt1[i] = 0U;
        }
        else
        {
            new1_val = stay1_val;
            new1_count = stay1_count;
            bt1[i] = 1U;
        }

        prev0_val = new0_val;
        prev0_count = new0_count;
        prev1_val = new1_val;
        prev1_count = new1_count;
    }

    if (prev1_val > prev0_val ||
        (prev1_val == prev0_val && prev1_count < prev0_count))
    {
        best_val = prev1_val;
        best_count = prev1_count;
        state = 1;
    }
    else
    {
        best_val = prev0_val;
        best_count = prev0_count;
        state = 0;
    }

    solution[n - 1] = (npy_uint8)state;
    for (i = n - 1; i > 0; --i)
    {
        state = (state == 0) ? (int)bt0[i] : (int)bt1[i];
        solution[i - 1] = (npy_uint8)state;
    }
    Py_END_ALLOW_THREADS

    free(bt0);
    free(bt1);
    Py_DECREF(scores_arr);
    Py_DECREF(switch_costs_arr);

    return Py_BuildValue(
        "Ndl",
        solution_arr,
        best_val,
        best_count);

fail:
    if (bt0 != NULL)
    {
        free(bt0);
    }
    if (bt1 != NULL)
    {
        free(bt1);
    }
    Py_XDECREF(solution_arr);
    Py_XDECREF(scores_arr);
    Py_XDECREF(switch_costs_arr);
    return NULL;
}

static PyMethodDef chain_dp_methods[] = {
    {
        "solve_penalized_chain",
        solve_penalized_chain,
        METH_VARARGS,
        "Solve the penalized binary chain problem."
    },
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef chain_dp_module = {
    PyModuleDef_HEAD_INIT,
    "_chain_dp",
    "Chain solver for ROCCO.",
    -1,
    chain_dp_methods,
};

PyMODINIT_FUNC PyInit__chain_dp(void)
{
    PyObject *module = PyModule_Create(&chain_dp_module);
    if (module == NULL)
    {
        return NULL;
    }
    import_array();
    return module;
}
