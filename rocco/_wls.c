#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/random/bitgen.h>

#include <math.h>
#include <stdint.h>

#include "native/wls_backend.h"

static uint64_t stationaryBootstrapBoundedIndex(
    bitgen_t *bitGenerator,
    uint64_t bound)
{
    uint64_t threshold = (UINT64_C(0) - bound) % bound;
    uint64_t value = 0U;

    do
    {
        value = bitGenerator->next_uint64(bitGenerator->state);
    } while (value < threshold);
    return value % bound;
}

static npy_intp stationaryBootstrapLocalIndex(
    bitgen_t *bitGenerator,
    npy_intp outputIndex,
    npy_intp size,
    npy_intp radius)
{
    npy_intp leftRadius = outputIndex;
    npy_intp rightRadius = size - outputIndex - 1;
    npy_intp lower = 0;
    npy_intp upper = 0;

    if (leftRadius > radius)
    {
        leftRadius = radius;
    }
    if (rightRadius > radius)
    {
        rightRadius = radius;
    }
    lower = outputIndex - leftRadius;
    upper = outputIndex + rightRadius + 1;
    return lower + (npy_intp)stationaryBootstrapBoundedIndex(
                       bitGenerator,
                       (uint64_t)(upper - lower));
}

static void stationaryNullBootstrapKernel(
    const double *templateData,
    double *drawData,
    npy_intp size,
    double meanBlockLength,
    npy_intp maxLocalRadiusIntervals,
    bitgen_t *bitGenerator)
{
    const int useLocalRadius = maxLocalRadiusIntervals >= 0;
    const double restartProbability = 1.0 / meanBlockLength;
    npy_intp localRadius = 0;
    npy_intp templateIndex = 0;
    npy_intp outputIndex = 0;
    double drawSum = 0.0;
    double drawMean = 0.0;

    if (useLocalRadius && size > 1)
    {
        if (meanBlockLength >= (double)size)
        {
            localRadius = size - 1;
        }
        else
        {
            double radiusValue = ceil(sqrt(meanBlockLength * (double)size));
            localRadius = radiusValue >= (double)(size - 1)
                              ? size - 1
                              : (npy_intp)radiusValue;
        }
        if (localRadius > maxLocalRadiusIntervals)
        {
            localRadius = maxLocalRadiusIntervals;
        }
    }

    templateIndex = useLocalRadius
                        ? stationaryBootstrapLocalIndex(
                              bitGenerator,
                              0,
                              size,
                              localRadius)
                        : (npy_intp)stationaryBootstrapBoundedIndex(
                              bitGenerator,
                              (uint64_t)size);
    drawData[0] = templateData[templateIndex];
    drawSum = drawData[0];

    for (outputIndex = 1; outputIndex < size; ++outputIndex)
    {
        if (bitGenerator->next_double(bitGenerator->state) < restartProbability)
        {
            templateIndex = useLocalRadius
                                ? stationaryBootstrapLocalIndex(
                                      bitGenerator,
                                      outputIndex,
                                      size,
                                      localRadius)
                                : (npy_intp)stationaryBootstrapBoundedIndex(
                                      bitGenerator,
                                      (uint64_t)size);
        }
        else
        {
            ++templateIndex;
            if (templateIndex == size)
            {
                templateIndex = 0;
            }
        }
        drawData[outputIndex] = templateData[templateIndex];
        drawSum += drawData[outputIndex];
    }

    drawMean = drawSum / (double)size;
    for (outputIndex = 0; outputIndex < size; ++outputIndex)
    {
        drawData[outputIndex] -= drawMean;
    }
}

static PyObject *stationaryNullBootstrapDraw(
    PyObject *self,
    PyObject *args,
    PyObject *kwargs)
{
    static char *kwlist[] = {
        "template",
        "meanBlockLength",
        "rng",
        "maxLocalRadiusIntervals",
        NULL,
    };
    PyObject *templateObject = NULL;
    PyObject *rngObject = NULL;
    PyObject *randomModule = NULL;
    PyObject *generatorType = NULL;
    PyObject *bitGeneratorObject = NULL;
    PyObject *capsule = NULL;
    PyObject *lockObject = NULL;
    PyObject *lockResult = NULL;
    PyArrayObject *templateArray = NULL;
    PyArrayObject *drawArray = NULL;
    bitgen_t *bitGenerator = NULL;
    const double *templateData = NULL;
    double meanBlockLength = 0.0;
    Py_ssize_t maxLocalRadiusIntervals = -1;
    npy_intp size = 0;
    npy_intp outputDimensions[1] = {0};
    npy_intp index = 0;
    int isGenerator = 0;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OdO|n",
            kwlist,
            &templateObject,
            &meanBlockLength,
            &rngObject,
            &maxLocalRadiusIntervals))
    {
        return NULL;
    }
    if (!isfinite(meanBlockLength) || meanBlockLength < 1.0)
    {
        PyErr_SetString(
            PyExc_ValueError,
            "`meanBlockLength` must be finite and at least 1");
        return NULL;
    }
    if (maxLocalRadiusIntervals < -1)
    {
        PyErr_SetString(
            PyExc_ValueError,
            "`maxLocalRadiusIntervals` must be -1 or non-negative");
        return NULL;
    }

    randomModule = PyImport_ImportModule("numpy.random");
    if (randomModule == NULL)
    {
        goto error;
    }
    generatorType = PyObject_GetAttrString(randomModule, "Generator");
    if (generatorType == NULL)
    {
        goto error;
    }
    isGenerator = PyObject_IsInstance(rngObject, generatorType);
    if (isGenerator < 0)
    {
        goto error;
    }
    if (isGenerator == 0)
    {
        PyErr_SetString(PyExc_TypeError, "`rng` must be a numpy.random.Generator");
        goto error;
    }
    Py_CLEAR(generatorType);
    Py_CLEAR(randomModule);

    templateArray = (PyArrayObject *)PyArray_FROM_OTF(
        templateObject,
        NPY_FLOAT64,
        NPY_ARRAY_IN_ARRAY);
    if (templateArray == NULL)
    {
        goto error;
    }
    if (PyArray_NDIM(templateArray) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "`template` must be one-dimensional");
        goto error;
    }
    size = PyArray_DIM(templateArray, 0);
    if (size <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "`template` must be non-empty");
        goto error;
    }
    templateData = (const double *)PyArray_DATA(templateArray);
    for (index = 0; index < size; ++index)
    {
        if (!isfinite(templateData[index]))
        {
            PyErr_SetString(PyExc_ValueError, "`template` values must be finite");
            goto error;
        }
    }

    bitGeneratorObject = PyObject_GetAttrString(rngObject, "bit_generator");
    if (bitGeneratorObject == NULL)
    {
        goto error;
    }
    capsule = PyObject_GetAttrString(bitGeneratorObject, "capsule");
    if (capsule == NULL)
    {
        goto error;
    }
    if (!PyCapsule_IsValid(capsule, "BitGenerator"))
    {
        PyErr_SetString(PyExc_ValueError, "`rng` has an invalid BitGenerator capsule");
        goto error;
    }
    bitGenerator = (bitgen_t *)PyCapsule_GetPointer(capsule, "BitGenerator");
    if (bitGenerator == NULL)
    {
        goto error;
    }
    lockObject = PyObject_GetAttrString(bitGeneratorObject, "lock");
    if (lockObject == NULL)
    {
        goto error;
    }

    outputDimensions[0] = size;
    drawArray = (PyArrayObject *)PyArray_EMPTY(1, outputDimensions, NPY_FLOAT64, 0);
    if (drawArray == NULL)
    {
        goto error;
    }
    lockResult = PyObject_CallMethod(lockObject, "acquire", NULL);
    if (lockResult == NULL)
    {
        goto error;
    }
    if (lockResult == Py_False)
    {
        PyErr_SetString(PyExc_RuntimeError, "Failed to acquire the BitGenerator lock");
        goto error;
    }
    Py_CLEAR(lockResult);

    Py_BEGIN_ALLOW_THREADS
        stationaryNullBootstrapKernel(
            templateData,
            (double *)PyArray_DATA(drawArray),
            size,
            meanBlockLength,
            (npy_intp)maxLocalRadiusIntervals,
            bitGenerator);
    Py_END_ALLOW_THREADS

    lockResult = PyObject_CallMethod(lockObject, "release", NULL);
    if (lockResult == NULL)
    {
        goto error;
    }
    Py_CLEAR(lockResult);
    Py_CLEAR(lockObject);
    Py_CLEAR(capsule);
    Py_CLEAR(bitGeneratorObject);
    Py_CLEAR(templateArray);
    return (PyObject *)drawArray;

error:
    Py_XDECREF(lockResult);
    Py_XDECREF(drawArray);
    Py_XDECREF(lockObject);
    Py_XDECREF(capsule);
    Py_XDECREF(bitGeneratorObject);
    Py_XDECREF(templateArray);
    Py_XDECREF(generatorType);
    Py_XDECREF(randomModule);
    return NULL;
}

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
        "stationaryNullBootstrapDraw",
        (PyCFunction)stationaryNullBootstrapDraw,
        METH_VARARGS | METH_KEYWORDS,
        NULL,
    },
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
