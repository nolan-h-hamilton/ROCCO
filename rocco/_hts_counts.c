#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "native/ccounts_backend.h"

static ccounts_sourceConfig make_source_config(const char *alignment_path)
{
    ccounts_sourceConfig source_config;
    source_config.path = alignment_path;
    source_config.sourceKind = ccounts_sourceKindBAM;
    source_config.barcodeTag = NULL;
    source_config.barcodeAllowListFile = NULL;
    source_config.barcodeGroupMapFile = NULL;
    return source_config;
}

static void raise_native_error(ccounts_result result)
{
    if (result.errorCode == 0)
    {
        return;
    }
    if (result.errorMessage != NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, result.errorMessage);
        return;
    }
    PyErr_SetString(PyExc_RuntimeError, "HTS counting failed");
}

static int parse_count_mode(const char *count_mode, uint8_t *count_mode_out)
{
    if (count_mode == NULL || count_mode_out == NULL)
    {
        PyErr_SetString(PyExc_ValueError, "count mode is invalid");
        return 0;
    }
    if (
        strcmp(count_mode, "coverage") == 0 ||
        strcmp(count_mode, "cov") == 0 ||
        strcmp(count_mode, "0") == 0)
    {
        *count_mode_out = (uint8_t)ccounts_countModeCoverage;
        return 1;
    }
    if (
        strcmp(count_mode, "cutsite") == 0 ||
        strcmp(count_mode, "cut") == 0 ||
        strcmp(count_mode, "cutsites") == 0 ||
        strcmp(count_mode, "1") == 0)
    {
        *count_mode_out = (uint8_t)ccounts_countModeCutSite;
        return 1;
    }
    if (
        strcmp(count_mode, "fiveprime") == 0 ||
        strcmp(count_mode, "five_prime") == 0 ||
        strcmp(count_mode, "5p") == 0 ||
        strcmp(count_mode, "2") == 0)
    {
        *count_mode_out = (uint8_t)ccounts_countModeFivePrime;
        return 1;
    }
    if (
        strcmp(count_mode, "center") == 0 ||
        strcmp(count_mode, "centre") == 0 ||
        strcmp(count_mode, "midpoint") == 0 ||
        strcmp(count_mode, "3") == 0)
    {
        *count_mode_out = (uint8_t)ccounts_countModeCenter;
        return 1;
    }
    PyErr_Format(PyExc_ValueError, "unsupported count mode `%s`", count_mode);
    return 0;
}

static PyObject *is_alignment_paired_end(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {
        "alignment_path",
        "max_reads",
        "thread_count",
        NULL,
    };
    const char *alignment_path = NULL;
    int max_reads = 1000;
    int thread_count = 0;
    int is_paired_end_out = 0;
    ccounts_result result;
    ccounts_sourceConfig source_config;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "s|ii",
            kwlist,
            &alignment_path,
            &max_reads,
            &thread_count))
    {
        return NULL;
    }

    source_config = make_source_config(alignment_path);
    Py_BEGIN_ALLOW_THREADS
    result = ccounts_isPairedEnd(
        &source_config,
        thread_count,
        max_reads,
        &is_paired_end_out);
    Py_END_ALLOW_THREADS
    if (result.errorCode != 0)
    {
        raise_native_error(result);
        return NULL;
    }
    if (is_paired_end_out)
    {
        Py_RETURN_TRUE;
    }
    Py_RETURN_FALSE;
}

static PyObject *get_alignment_read_length(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {
        "alignment_path",
        "min_reads",
        "thread_count",
        "max_iterations",
        "flag_exclude",
        NULL,
    };
    const char *alignment_path = NULL;
    int min_reads = 32;
    int thread_count = 0;
    int max_iterations = 4096;
    int flag_exclude = 0;
    uint32_t read_length_out = 0;
    ccounts_result result;
    ccounts_sourceConfig source_config;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "s|iiii",
            kwlist,
            &alignment_path,
            &min_reads,
            &thread_count,
            &max_iterations,
            &flag_exclude))
    {
        return NULL;
    }

    source_config = make_source_config(alignment_path);
    Py_BEGIN_ALLOW_THREADS
    result = ccounts_getReadLength(
        &source_config,
        thread_count,
        min_reads,
        max_iterations,
        flag_exclude,
        &read_length_out);
    Py_END_ALLOW_THREADS
    if (result.errorCode != 0)
    {
        raise_native_error(result);
        return NULL;
    }
    return PyLong_FromUnsignedLong((unsigned long)read_length_out);
}

static PyObject *get_alignment_fragment_length(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {
        "alignment_path",
        "thread_count",
        "flag_exclude",
        "max_iterations",
        "max_insert_size",
        "block_size",
        "rolling_chunk_size",
        "lag_step",
        "early_exit",
        "fallback",
        NULL,
    };
    const char *alignment_path = NULL;
    int thread_count = 0;
    int flag_exclude = 0;
    int max_iterations = 1000;
    int max_insert_size = 1000;
    int block_size = 5000;
    int rolling_chunk_size = 250;
    int lag_step = 5;
    int early_exit = 250;
    int fallback = 0;
    uint32_t fragment_length_out = 0;
    ccounts_result result;
    ccounts_sourceConfig source_config;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "s|iiiiiiiii",
            kwlist,
            &alignment_path,
            &thread_count,
            &flag_exclude,
            &max_iterations,
            &max_insert_size,
            &block_size,
            &rolling_chunk_size,
            &lag_step,
            &early_exit,
            &fallback))
    {
        return NULL;
    }

    source_config = make_source_config(alignment_path);
    Py_BEGIN_ALLOW_THREADS
    result = ccounts_getFragmentLength(
        &source_config,
        thread_count,
        flag_exclude,
        max_iterations,
        max_insert_size,
        block_size,
        rolling_chunk_size,
        lag_step,
        early_exit,
        fallback,
        &fragment_length_out);
    Py_END_ALLOW_THREADS
    if (result.errorCode != 0)
    {
        raise_native_error(result);
        return NULL;
    }
    return PyLong_FromUnsignedLong((unsigned long)fragment_length_out);
}

static PyObject *get_alignment_chrom_range(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {
        "alignment_path",
        "chromosome",
        "chrom_length",
        "thread_count",
        "flag_exclude",
        NULL,
    };
    const char *alignment_path = NULL;
    const char *chromosome = NULL;
    unsigned long long chrom_length = 0ULL;
    int thread_count = 0;
    int flag_exclude = 0;
    uint64_t start_out = 0ULL;
    uint64_t end_out = 0ULL;
    ccounts_result result;
    ccounts_sourceConfig source_config;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "ssK|ii",
            kwlist,
            &alignment_path,
            &chromosome,
            &chrom_length,
            &thread_count,
            &flag_exclude))
    {
        return NULL;
    }

    source_config = make_source_config(alignment_path);
    Py_BEGIN_ALLOW_THREADS
    result = ccounts_getChromRange(
        &source_config,
        chromosome,
        (uint64_t)chrom_length,
        thread_count,
        flag_exclude,
        &start_out,
        &end_out);
    Py_END_ALLOW_THREADS
    if (result.errorCode != 0)
    {
        raise_native_error(result);
        return NULL;
    }

    return Py_BuildValue("KK", (unsigned long long)start_out, (unsigned long long)end_out);
}

static PyObject *get_alignment_mapped_read_count(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {
        "alignment_path",
        "exclude_chromosomes",
        "thread_count",
        "count_mode",
        "one_read_per_bin",
        NULL,
    };
    const char *alignment_path = NULL;
    PyObject *exclude_obj = Py_None;
    PyObject *exclude_seq = NULL;
    const char *count_mode = "coverage";
    const char **exclude_ptrs = NULL;
    Py_ssize_t exclude_count = 0;
    Py_ssize_t idx = 0;
    int thread_count = 0;
    int one_read_per_bin = 0;
    uint8_t count_mode_code = (uint8_t)ccounts_countModeCoverage;
    uint64_t mapped_out = 0ULL;
    uint64_t unmapped_out = 0ULL;
    ccounts_result result;
    ccounts_sourceConfig source_config;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "s|Oisi",
            kwlist,
            &alignment_path,
            &exclude_obj,
            &thread_count,
            &count_mode,
            &one_read_per_bin))
    {
        return NULL;
    }
    if (!parse_count_mode(count_mode, &count_mode_code))
    {
        return NULL;
    }

    if (exclude_obj != Py_None)
    {
        exclude_seq = PySequence_Fast(
            exclude_obj,
            "`exclude_chromosomes` must be a sequence of strings");
        if (exclude_seq == NULL)
        {
            return NULL;
        }
        exclude_count = PySequence_Fast_GET_SIZE(exclude_seq);
        if (exclude_count > 0)
        {
            exclude_ptrs = (const char **)calloc(
                (size_t)exclude_count,
                sizeof(const char *));
            if (exclude_ptrs == NULL)
            {
                Py_DECREF(exclude_seq);
                return PyErr_NoMemory();
            }
            for (idx = 0; idx < exclude_count; ++idx)
            {
                PyObject *item = PySequence_Fast_GET_ITEM(exclude_seq, idx);
                if (!PyUnicode_Check(item))
                {
                    PyErr_SetString(
                        PyExc_TypeError,
                        "`exclude_chromosomes` entries must be strings");
                    free(exclude_ptrs);
                    Py_DECREF(exclude_seq);
                    return NULL;
                }
                exclude_ptrs[idx] = PyUnicode_AsUTF8(item);
                if (exclude_ptrs[idx] == NULL)
                {
                    free(exclude_ptrs);
                    Py_DECREF(exclude_seq);
                    return NULL;
                }
            }
        }
    }

    source_config = make_source_config(alignment_path);
    Py_BEGIN_ALLOW_THREADS
    result = ccounts_getMappedReadCount(
        &source_config,
        thread_count,
        exclude_ptrs,
        (int)exclude_count,
        count_mode_code,
        (uint8_t)(one_read_per_bin != 0),
        &mapped_out,
        &unmapped_out);
    Py_END_ALLOW_THREADS

    free(exclude_ptrs);
    Py_XDECREF(exclude_seq);

    if (result.errorCode != 0)
    {
        raise_native_error(result);
        return NULL;
    }

    return Py_BuildValue(
        "KK",
        (unsigned long long)mapped_out,
        (unsigned long long)unmapped_out);
}

static PyObject *count_alignment_region(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {
        "alignment_path",
        "chromosome",
        "start",
        "end",
        "interval_size_bp",
        "read_length",
        "one_read_per_bin",
        "thread_count",
        "flag_include",
        "flag_exclude",
        "shift_forward_strand53",
        "shift_reverse_strand53",
        "extend_bp",
        "max_insert_size",
        "paired_end_mode",
        "infer_fragment_length",
        "min_mapping_quality",
        "min_template_length",
        "count_mode",
        NULL,
    };
    const char *alignment_path = NULL;
    const char *chromosome = NULL;
    const char *count_mode = "coverage";
    int start = 0;
    int end = 0;
    int interval_size_bp = 0;
    int read_length = 0;
    int one_read_per_bin = 0;
    int thread_count = 0;
    int flag_include = 0;
    int flag_exclude = 0;
    int shift_forward_strand53 = 0;
    int shift_reverse_strand53 = 0;
    int extend_bp = 0;
    int max_insert_size = 1000;
    int paired_end_mode = 0;
    int infer_fragment_length = 0;
    int min_mapping_quality = 0;
    int min_template_length = -1;
    uint8_t count_mode_code = (uint8_t)ccounts_countModeCoverage;
    ccounts_sourceConfig source_config;
    ccounts_sourceHandle *source_handle = NULL;
    ccounts_region region;
    ccounts_countOptions count_options;
    ccounts_result result;
    PyArrayObject *counts_arr = NULL;
    npy_intp num_intervals = 0;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "ssiiii|iiiiiiiiiiiis",
            kwlist,
            &alignment_path,
            &chromosome,
            &start,
            &end,
            &interval_size_bp,
            &read_length,
            &one_read_per_bin,
            &thread_count,
            &flag_include,
            &flag_exclude,
            &shift_forward_strand53,
            &shift_reverse_strand53,
            &extend_bp,
            &max_insert_size,
            &paired_end_mode,
            &infer_fragment_length,
            &min_mapping_quality,
            &min_template_length,
            &count_mode))
    {
        return NULL;
    }

    if (interval_size_bp <= 0 || end <= start)
    {
        PyErr_SetString(
            PyExc_ValueError,
            "invalid interval size or genomic segment");
        return NULL;
    }
    if (!parse_count_mode(count_mode, &count_mode_code))
    {
        return NULL;
    }

    num_intervals = (npy_intp)(((end - start - 1) / interval_size_bp) + 1);
    counts_arr = (PyArrayObject *)PyArray_ZEROS(
        1,
        &num_intervals,
        NPY_FLOAT32,
        0);
    if (counts_arr == NULL)
    {
        return NULL;
    }

    source_config = make_source_config(alignment_path);
    region.chromosome = chromosome;
    region.start = (uint32_t)start;
    region.end = (uint32_t)end;
    region.intervalSizeBP = (uint32_t)interval_size_bp;
    memset(&count_options, 0, sizeof(ccounts_countOptions));
    count_options.threadCount = (uint16_t)((thread_count > 0) ? thread_count : 0);
    count_options.flagInclude = (uint16_t)((flag_include > 0) ? flag_include : 0);
    count_options.flagExclude = (uint16_t)((flag_exclude > 0) ? flag_exclude : 0);
    count_options.countMode = count_mode_code;
    count_options.oneReadPerBin = (uint8_t)(one_read_per_bin != 0);
    count_options.shiftForwardStrand53 = (int64_t)shift_forward_strand53;
    count_options.shiftReverseStrand53 = (int64_t)shift_reverse_strand53;
    count_options.readLength = (int64_t)read_length;
    count_options.extendBP = (int64_t)extend_bp;
    count_options.minMappingQuality = (int64_t)min_mapping_quality;
    count_options.minTemplateLength = (int64_t)min_template_length;
    count_options.maxInsertSize = (int64_t)max_insert_size;
    count_options.pairedEndMode = (int64_t)paired_end_mode;
    count_options.inferFragmentLength = (int64_t)infer_fragment_length;

    Py_BEGIN_ALLOW_THREADS
    result = ccounts_openSource(&source_config, &source_handle);
    if (result.errorCode == 0)
    {
        result = ccounts_countRegion(
            source_handle,
            &region,
            &count_options,
            (float *)PyArray_DATA(counts_arr),
            (size_t)num_intervals);
    }
    if (source_handle != NULL)
    {
        ccounts_closeSource(source_handle);
    }
    Py_END_ALLOW_THREADS

    if (result.errorCode != 0)
    {
        Py_DECREF(counts_arr);
        raise_native_error(result);
        return NULL;
    }

    return (PyObject *)counts_arr;
}

static PyObject *count_alignment_intervals(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {
        "alignment_path",
        "chromosomes",
        "starts",
        "ends",
        "read_length",
        "one_read_per_bin",
        "thread_count",
        "flag_include",
        "flag_exclude",
        "shift_forward_strand53",
        "shift_reverse_strand53",
        "extend_bp",
        "max_insert_size",
        "paired_end_mode",
        "infer_fragment_length",
        "min_mapping_quality",
        "min_template_length",
        "count_mode",
        NULL,
    };
    const char *alignment_path = NULL;
    const char *count_mode = "coverage";
    PyObject *chromosomes_obj = NULL;
    PyObject *starts_obj = NULL;
    PyObject *ends_obj = NULL;
    PyObject *chrom_seq = NULL;
    PyObject *start_seq = NULL;
    PyObject *end_seq = NULL;
    PyArrayObject *counts_arr = NULL;
    float *count_data = NULL;
    npy_intp count_dims[1];
    Py_ssize_t interval_count = 0;
    Py_ssize_t interval_index = 0;
    int read_length = 0;
    int one_read_per_bin = 0;
    int thread_count = 0;
    int flag_include = 0;
    int flag_exclude = 0;
    int shift_forward_strand53 = 0;
    int shift_reverse_strand53 = 0;
    int extend_bp = 0;
    int max_insert_size = 1000;
    int paired_end_mode = 0;
    int infer_fragment_length = 0;
    int min_mapping_quality = 0;
    int min_template_length = -1;
    uint8_t count_mode_code = (uint8_t)ccounts_countModeCoverage;
    ccounts_sourceConfig source_config;
    ccounts_sourceHandle *source_handle = NULL;
    ccounts_region region;
    ccounts_countOptions count_options;
    ccounts_result result = {0, NULL};

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "sOOO|iiiiiiiiiiiiis",
            kwlist,
            &alignment_path,
            &chromosomes_obj,
            &starts_obj,
            &ends_obj,
            &read_length,
            &one_read_per_bin,
            &thread_count,
            &flag_include,
            &flag_exclude,
            &shift_forward_strand53,
            &shift_reverse_strand53,
            &extend_bp,
            &max_insert_size,
            &paired_end_mode,
            &infer_fragment_length,
            &min_mapping_quality,
            &min_template_length,
            &count_mode))
    {
        return NULL;
    }

    if (!parse_count_mode(count_mode, &count_mode_code))
    {
        return NULL;
    }

    chrom_seq = PySequence_Fast(
        chromosomes_obj,
        "`chromosomes` must be a sequence of strings");
    start_seq = PySequence_Fast(
        starts_obj,
        "`starts` must be a sequence of integers");
    end_seq = PySequence_Fast(
        ends_obj,
        "`ends` must be a sequence of integers");
    if (chrom_seq == NULL || start_seq == NULL || end_seq == NULL)
    {
        Py_XDECREF(chrom_seq);
        Py_XDECREF(start_seq);
        Py_XDECREF(end_seq);
        return NULL;
    }

    interval_count = PySequence_Fast_GET_SIZE(chrom_seq);
    if (PySequence_Fast_GET_SIZE(start_seq) != interval_count ||
        PySequence_Fast_GET_SIZE(end_seq) != interval_count)
    {
        Py_DECREF(chrom_seq);
        Py_DECREF(start_seq);
        Py_DECREF(end_seq);
        PyErr_SetString(
            PyExc_ValueError,
            "`chromosomes`, `starts`, and `ends` must have the same length");
        return NULL;
    }

    count_dims[0] = (npy_intp)interval_count;
    counts_arr = (PyArrayObject *)PyArray_ZEROS(
        1,
        count_dims,
        NPY_FLOAT32,
        0);
    if (counts_arr == NULL)
    {
        Py_DECREF(chrom_seq);
        Py_DECREF(start_seq);
        Py_DECREF(end_seq);
        return NULL;
    }

    source_config = make_source_config(alignment_path);
    result = ccounts_openSource(&source_config, &source_handle);
    if (result.errorCode != 0)
    {
        Py_DECREF(counts_arr);
        Py_DECREF(chrom_seq);
        Py_DECREF(start_seq);
        Py_DECREF(end_seq);
        raise_native_error(result);
        return NULL;
    }

    memset(&count_options, 0, sizeof(ccounts_countOptions));
    count_options.threadCount = (uint16_t)((thread_count > 0) ? thread_count : 0);
    count_options.flagInclude = (uint16_t)((flag_include > 0) ? flag_include : 0);
    count_options.flagExclude = (uint16_t)((flag_exclude > 0) ? flag_exclude : 0);
    count_options.countMode = count_mode_code;
    count_options.oneReadPerBin = (uint8_t)(one_read_per_bin != 0);
    count_options.shiftForwardStrand53 = (int64_t)shift_forward_strand53;
    count_options.shiftReverseStrand53 = (int64_t)shift_reverse_strand53;
    count_options.readLength = (int64_t)read_length;
    count_options.extendBP = (int64_t)extend_bp;
    count_options.minMappingQuality = (int64_t)min_mapping_quality;
    count_options.minTemplateLength = (int64_t)min_template_length;
    count_options.maxInsertSize = (int64_t)max_insert_size;
    count_options.pairedEndMode = (int64_t)paired_end_mode;
    count_options.inferFragmentLength = (int64_t)infer_fragment_length;
    count_data = (float *)PyArray_DATA(counts_arr);

    for (interval_index = 0; interval_index < interval_count; ++interval_index)
    {
        PyObject *chrom_obj = PySequence_Fast_GET_ITEM(chrom_seq, interval_index);
        PyObject *start_obj = PySequence_Fast_GET_ITEM(start_seq, interval_index);
        PyObject *end_obj = PySequence_Fast_GET_ITEM(end_seq, interval_index);
        const char *chromosome = NULL;
        unsigned long long start_val = 0ULL;
        unsigned long long end_val = 0ULL;

        chromosome = PyUnicode_AsUTF8(chrom_obj);
        if (chromosome == NULL)
        {
            ccounts_closeSource(source_handle);
            Py_DECREF(counts_arr);
            Py_DECREF(chrom_seq);
            Py_DECREF(start_seq);
            Py_DECREF(end_seq);
            return NULL;
        }
        start_val = PyLong_AsUnsignedLongLong(start_obj);
        if (PyErr_Occurred())
        {
            ccounts_closeSource(source_handle);
            Py_DECREF(counts_arr);
            Py_DECREF(chrom_seq);
            Py_DECREF(start_seq);
            Py_DECREF(end_seq);
            return NULL;
        }
        end_val = PyLong_AsUnsignedLongLong(end_obj);
        if (PyErr_Occurred())
        {
            ccounts_closeSource(source_handle);
            Py_DECREF(counts_arr);
            Py_DECREF(chrom_seq);
            Py_DECREF(start_seq);
            Py_DECREF(end_seq);
            return NULL;
        }
        if (start_val > (unsigned long long)UINT32_MAX)
        {
            ccounts_closeSource(source_handle);
            Py_DECREF(counts_arr);
            Py_DECREF(chrom_seq);
            Py_DECREF(start_seq);
            Py_DECREF(end_seq);
            PyErr_SetString(
                PyExc_ValueError,
                "interval start exceeds uint32 range");
            return NULL;
        }
        if (end_val <= start_val)
        {
            ccounts_closeSource(source_handle);
            Py_DECREF(counts_arr);
            Py_DECREF(chrom_seq);
            Py_DECREF(start_seq);
            Py_DECREF(end_seq);
            PyErr_SetString(
                PyExc_ValueError,
                "each interval must satisfy end > start");
            return NULL;
        }
        if (end_val > (unsigned long long)UINT32_MAX)
        {
            ccounts_closeSource(source_handle);
            Py_DECREF(counts_arr);
            Py_DECREF(chrom_seq);
            Py_DECREF(start_seq);
            Py_DECREF(end_seq);
            PyErr_SetString(
                PyExc_ValueError,
                "interval end exceeds uint32 range");
            return NULL;
        }

        region.chromosome = chromosome;
        region.start = (uint32_t)start_val;
        region.end = (uint32_t)end_val;
        region.intervalSizeBP = (uint32_t)(end_val - start_val);
        count_data[interval_index] = 0.0f;
        result = ccounts_countRegion(
            source_handle,
            &region,
            &count_options,
            &count_data[interval_index],
            1U);
        if (result.errorCode != 0)
        {
            ccounts_closeSource(source_handle);
            Py_DECREF(counts_arr);
            Py_DECREF(chrom_seq);
            Py_DECREF(start_seq);
            Py_DECREF(end_seq);
            raise_native_error(result);
            return NULL;
        }
    }

    ccounts_closeSource(source_handle);
    Py_DECREF(chrom_seq);
    Py_DECREF(start_seq);
    Py_DECREF(end_seq);
    return (PyObject *)counts_arr;
}

static PyMethodDef hts_counts_methods[] = {
    {
        "is_alignment_paired_end",
        (PyCFunction)is_alignment_paired_end,
        METH_VARARGS | METH_KEYWORDS,
        "Return whether an alignment file appears to be paired-end."
    },
    {
        "get_alignment_read_length",
        (PyCFunction)get_alignment_read_length,
        METH_VARARGS | METH_KEYWORDS,
        "Estimate a representative read length from an alignment file."
    },
    {
        "get_alignment_fragment_length",
        (PyCFunction)get_alignment_fragment_length,
        METH_VARARGS | METH_KEYWORDS,
        "Estimate a representative fragment length from an alignment file."
    },
    {
        "get_alignment_chrom_range",
        (PyCFunction)get_alignment_chrom_range,
        METH_VARARGS | METH_KEYWORDS,
        "Return the covered start and end coordinates for one chromosome."
    },
    {
        "get_alignment_mapped_read_count",
        (PyCFunction)get_alignment_mapped_read_count,
        METH_VARARGS | METH_KEYWORDS,
        "Return mapped and unmapped read counts from an alignment file."
    },
    {
        "count_alignment_region",
        (PyCFunction)count_alignment_region,
        METH_VARARGS | METH_KEYWORDS,
        "Count an evenly spaced region from an alignment file."
    },
    {
        "count_alignment_intervals",
        (PyCFunction)count_alignment_intervals,
        METH_VARARGS | METH_KEYWORDS,
        "Count a list of genomic intervals from an alignment file."
    },
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef hts_counts_module = {
    PyModuleDef_HEAD_INIT,
    "_hts_counts",
    "HTS counting for ROCCO.",
    -1,
    hts_counts_methods,
};

PyMODINIT_FUNC PyInit__hts_counts(void)
{
    PyObject *module = PyModule_Create(&hts_counts_module);
    if (module == NULL)
    {
        return NULL;
    }
    import_array();
    return module;
}
