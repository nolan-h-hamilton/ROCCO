#include "ccounts_backend.h"

#include <htslib/hts.h>
#include <htslib/khash.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>

#include <ctype.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief native handle for one open source
 *
 * keeps the htslib state for either alignments
 * or tabix fragments, plus any loaded allowlist
 */
struct ccounts_sourceHandle
{
    ccounts_sourceKind sourceKind;
    samFile *fileHandle;
    sam_hdr_t *header;
    hts_idx_t *indexHandle;
    htsFile *fragmentsHandle;
    tbx_t *tbxHandle;
    char **barcodeAllowList;
    int barcodeAllowCount;
};

typedef struct ccounts_scoredIndex
{
    double value;
    int index;
} ccounts_scoredIndex;

KHASH_SET_INIT_STR(ccounts_barcodeSet)

static ccounts_result ccounts_makeResult(int32_t errorCode, const char *errorMessage)
{
    ccounts_result result;
    result.errorCode = errorCode;
    result.errorMessage = errorMessage;
    return result;
}

static ccounts_result ccounts_makeOk(void)
{
    return ccounts_makeResult(0, NULL);
}

static int ccounts_hasValue(const char *value)
{
    return value != NULL && value[0] != '\0';
}

static int ccounts_isAlignmentKind(ccounts_sourceKind sourceKind)
{
    return sourceKind == ccounts_sourceKindBAM;
}

static int ccounts_compareUint32(const void *left, const void *right)
{
    uint32_t leftValue = *((const uint32_t *)left);
    uint32_t rightValue = *((const uint32_t *)right);
    if (leftValue < rightValue)
    {
        return -1;
    }
    if (leftValue > rightValue)
    {
        return 1;
    }
    return 0;
}

static int ccounts_compareCStringPtr(const void *left, const void *right)
{
    const char *leftValue = *((const char *const *)left);
    const char *rightValue = *((const char *const *)right);
    if (leftValue == NULL && rightValue == NULL)
    {
        return 0;
    }
    if (leftValue == NULL)
    {
        return -1;
    }
    if (rightValue == NULL)
    {
        return 1;
    }
    return strcmp(leftValue, rightValue);
}

static int ccounts_compareScoredIndexDesc(const void *left, const void *right)
{
    const ccounts_scoredIndex *leftValue = (const ccounts_scoredIndex *)left;
    const ccounts_scoredIndex *rightValue = (const ccounts_scoredIndex *)right;
    if (leftValue->value > rightValue->value)
    {
        return -1;
    }
    if (leftValue->value < rightValue->value)
    {
        return 1;
    }
    if (leftValue->index < rightValue->index)
    {
        return -1;
    }
    if (leftValue->index > rightValue->index)
    {
        return 1;
    }
    return 0;
}

static int ccounts_stringInList(
    const char *value,
    const char *const *stringList,
    int stringCount)
{
    int stringIndex = 0;
    if (value == NULL || stringList == NULL || stringCount <= 0)
    {
        return 0;
    }
    for (stringIndex = 0; stringIndex < stringCount; ++stringIndex)
    {
        if (stringList[stringIndex] != NULL && strcmp(value, stringList[stringIndex]) == 0)
        {
            return 1;
        }
    }
    return 0;
}

static int ccounts_stringFieldInList(
    const char *value,
    size_t valueLength,
    const char *const *stringList,
    int stringCount)
{
    int stringIndex = 0;
    size_t targetLength = 0U;
    if (value == NULL || stringList == NULL || stringCount <= 0)
    {
        return 0;
    }
    for (stringIndex = 0; stringIndex < stringCount; ++stringIndex)
    {
        if (stringList[stringIndex] == NULL)
        {
            continue;
        }
        targetLength = strlen(stringList[stringIndex]);
        if (targetLength == valueLength && strncmp(value, stringList[stringIndex], valueLength) == 0)
        {
            return 1;
        }
    }
    return 0;
}

static ccounts_result ccounts_applySourceConfig(
    samFile *fileHandle,
    const ccounts_sourceConfig *sourceConfig,
    int threadCount)
{
    if (fileHandle == NULL || sourceConfig == NULL)
    {
        return ccounts_makeResult(-1, "source config is invalid");
    }
    if (threadCount > 1 && hts_set_threads((htsFile *)fileHandle, threadCount) != 0)
    {
        return ccounts_makeResult(-1, "failed to set htslib thread count");
    }
    return ccounts_makeOk();
}

static void ccounts_freeBarcodeAllowList(
    char **barcodeAllowList,
    int barcodeAllowCount)
{
    int barcodeIndex = 0;
    if (barcodeAllowList == NULL)
    {
        return;
    }
    for (barcodeIndex = 0; barcodeIndex < barcodeAllowCount; ++barcodeIndex)
    {
        if (barcodeAllowList[barcodeIndex] != NULL)
        {
            free(barcodeAllowList[barcodeIndex]);
        }
    }
    free(barcodeAllowList);
}

static ccounts_result ccounts_loadBarcodeAllowList(
    const char *allowListPath,
    char ***barcodeAllowListOut,
    int *barcodeAllowCountOut)
{
    FILE *allowListFile = NULL;
    char **barcodeAllowList = NULL;
    int barcodeAllowCount = 0;
    int barcodeAllowCapacity = 0;
    char lineBuffer[4096];
    char *lineStart = NULL;
    char *lineEnd = NULL;
    size_t barcodeLength = 0;
    char *barcodeCopy = NULL;

    if (barcodeAllowListOut != NULL)
    {
        *barcodeAllowListOut = NULL;
    }
    if (barcodeAllowCountOut != NULL)
    {
        *barcodeAllowCountOut = 0;
    }
    if (!ccounts_hasValue(allowListPath))
    {
        return ccounts_makeOk();
    }

    allowListFile = fopen(allowListPath, "r");
    if (allowListFile == NULL)
    {
        return ccounts_makeResult(-1, "failed to open barcode allowlist");
    }

    while (fgets(lineBuffer, (int)sizeof(lineBuffer), allowListFile) != NULL)
    {
        lineStart = lineBuffer;
        while (*lineStart != '\0' && isspace((unsigned char)(*lineStart)))
        {
            ++lineStart;
        }
        if (*lineStart == '\0' || *lineStart == '#')
        {
            continue;
        }

        lineEnd = lineStart;
        while (*lineEnd != '\0' && *lineEnd != '\n' && *lineEnd != '\r' &&
               *lineEnd != '\t' && !isspace((unsigned char)(*lineEnd)))
        {
            ++lineEnd;
        }
        barcodeLength = (size_t)(lineEnd - lineStart);
        if (barcodeLength == 0U)
        {
            continue;
        }

        if (barcodeAllowCount >= barcodeAllowCapacity)
        {
            int nextCapacity = barcodeAllowCapacity > 0 ? (2 * barcodeAllowCapacity) : 256;
            char **nextAllowList = (char **)realloc(
                barcodeAllowList,
                (size_t)nextCapacity * sizeof(char *));
            if (nextAllowList == NULL)
            {
                fclose(allowListFile);
                ccounts_freeBarcodeAllowList(barcodeAllowList, barcodeAllowCount);
                return ccounts_makeResult(-1, "failed to expand barcode allowlist");
            }
            barcodeAllowList = nextAllowList;
            barcodeAllowCapacity = nextCapacity;
        }

        barcodeCopy = (char *)malloc(barcodeLength + 1U);
        if (barcodeCopy == NULL)
        {
            fclose(allowListFile);
            ccounts_freeBarcodeAllowList(barcodeAllowList, barcodeAllowCount);
            return ccounts_makeResult(-1, "failed to store barcode allowlist entry");
        }
        memcpy(barcodeCopy, lineStart, barcodeLength);
        barcodeCopy[barcodeLength] = '\0';
        barcodeAllowList[barcodeAllowCount] = barcodeCopy;
        ++barcodeAllowCount;
    }

    fclose(allowListFile);
    if (barcodeAllowCount > 1)
    {
        // keep allowlist sorted so membership checks stay cheap
        qsort(
            barcodeAllowList,
            (size_t)barcodeAllowCount,
            sizeof(char *),
            ccounts_compareCStringPtr);
    }
    if (barcodeAllowListOut != NULL)
    {
        *barcodeAllowListOut = barcodeAllowList;
    }
    if (barcodeAllowCountOut != NULL)
    {
        *barcodeAllowCountOut = barcodeAllowCount;
    }
    return ccounts_makeOk();
}

static int ccounts_barcodeAllowed(
    char **barcodeAllowList,
    int barcodeAllowCount,
    const char *barcodeStart,
    size_t barcodeLength)
{
    int low = 0;
    int high = barcodeAllowCount - 1;
    int mid = 0;
    int cmp = 0;
    const char *target = NULL;
    size_t targetLength = 0U;

    if (barcodeAllowCount <= 0 || barcodeAllowList == NULL)
    {
        return 1;
    }
    // allowlist stays sorted so this can stay a binary search
    while (low <= high)
    {
        mid = low + ((high - low) / 2);
        target = barcodeAllowList[mid];
        if (target == NULL)
        {
            return 0;
        }
        targetLength = strlen(target);
        cmp = strncmp(target, barcodeStart, barcodeLength);
        if (cmp == 0)
        {
            if (targetLength == barcodeLength)
            {
                return 1;
            }
            cmp = targetLength < barcodeLength ? -1 : 1;
        }
        if (cmp < 0)
        {
            low = mid + 1;
        }
        else
        {
            high = mid - 1;
        }
    }
    return 0;
}

static int ccounts_parseInt64Field(
    const char *fieldStart,
    size_t fieldLength,
    int64_t *valueOut)
{
    size_t fieldIndex = 0U;
    int sign = 1;
    int64_t value = 0;

    if (valueOut != NULL)
    {
        *valueOut = 0;
    }
    if (fieldStart == NULL || fieldLength == 0U)
    {
        return 0;
    }
    if (fieldStart[0] == '-')
    {
        sign = -1;
        fieldIndex = 1U;
    }
    for (; fieldIndex < fieldLength; ++fieldIndex)
    {
        if (fieldStart[fieldIndex] < '0' || fieldStart[fieldIndex] > '9')
        {
            return 0;
        }
        value = (10 * value) + (int64_t)(fieldStart[fieldIndex] - '0');
    }
    if (valueOut != NULL)
    {
        *valueOut = sign > 0 ? value : -value;
    }
    return 1;
}

static char *ccounts_copyStringField(
    const char *fieldStart,
    size_t fieldLength)
{
    char *fieldCopy = NULL;
    if (fieldStart == NULL || fieldLength == 0U)
    {
        return NULL;
    }
    fieldCopy = (char *)malloc(fieldLength + 1U);
    if (fieldCopy == NULL)
    {
        return NULL;
    }
    memcpy(fieldCopy, fieldStart, fieldLength);
    fieldCopy[fieldLength] = '\0';
    return fieldCopy;
}

static ccounts_result ccounts_openFragmentsSource(
    const ccounts_sourceConfig *sourceConfig,
    ccounts_sourceHandle *sourceHandle)
{
    ccounts_result result;

    if (sourceConfig == NULL || sourceHandle == NULL || !ccounts_hasValue(sourceConfig->path))
    {
        return ccounts_makeResult(-1, "fragments source config is invalid");
    }

    sourceHandle->fragmentsHandle = hts_open(sourceConfig->path, "r");
    if (sourceHandle->fragmentsHandle == NULL)
    {
        return ccounts_makeResult(-1, "failed to open fragments source");
    }

    // note, fragments support expects a bgzip file plus tabix index
    sourceHandle->tbxHandle = tbx_index_load(sourceConfig->path);
    if (sourceHandle->tbxHandle == NULL)
    {
        hts_close(sourceHandle->fragmentsHandle);
        sourceHandle->fragmentsHandle = NULL;
        return ccounts_makeResult(-1, "failed to load tabix index for fragments source");
    }

    result = ccounts_loadBarcodeAllowList(
        sourceConfig->barcodeAllowListFile,
        &sourceHandle->barcodeAllowList,
        &sourceHandle->barcodeAllowCount);
    if (result.errorCode != 0)
    {
        tbx_destroy(sourceHandle->tbxHandle);
        sourceHandle->tbxHandle = NULL;
        hts_close(sourceHandle->fragmentsHandle);
        sourceHandle->fragmentsHandle = NULL;
        return result;
    }
    return ccounts_makeOk();
}

static ccounts_result ccounts_openAlignmentFile(
    const ccounts_sourceConfig *sourceConfig,
    int threadCount,
    samFile **fileHandleOut,
    sam_hdr_t **headerOut)
{
    samFile *fileHandle = NULL;
    sam_hdr_t *header = NULL;
    ccounts_result result;

    if (fileHandleOut != NULL)
    {
        *fileHandleOut = NULL;
    }
    if (headerOut != NULL)
    {
        *headerOut = NULL;
    }

    if (sourceConfig == NULL || !ccounts_hasValue(sourceConfig->path))
    {
        return ccounts_makeResult(-1, "alignment path is empty");
    }
    if (!ccounts_isAlignmentKind(sourceConfig->sourceKind))
    {
        return ccounts_makeResult(-1, "source kind is not BAM");
    }

    fileHandle = sam_open(sourceConfig->path, "r");
    if (fileHandle == NULL)
    {
        return ccounts_makeResult(-1, "failed to open alignment source");
    }

    result = ccounts_applySourceConfig(fileHandle, sourceConfig, threadCount);
    if (result.errorCode != 0)
    {
        sam_close(fileHandle);
        return result;
    }

    header = sam_hdr_read(fileHandle);
    if (header == NULL)
    {
        sam_close(fileHandle);
        return ccounts_makeResult(-1, "failed to read alignment header");
    }

    if (fileHandleOut != NULL)
    {
        *fileHandleOut = fileHandle;
    }
    if (headerOut != NULL)
    {
        *headerOut = header;
    }
    return ccounts_makeOk();
}

static void ccounts_closeAlignmentFile(samFile *fileHandle, sam_hdr_t *header)
{
    if (header != NULL)
    {
        sam_hdr_destroy(header);
    }
    if (fileHandle != NULL)
    {
        sam_close(fileHandle);
    }
}

ccounts_result ccounts_checkAlignmentFile(
    const ccounts_sourceConfig *sourceConfig,
    int buildIndex,
    int threadCount,
    int *hasIndexOut)
{
    samFile *fileHandle = NULL;
    sam_hdr_t *header = NULL;
    hts_idx_t *indexHandle = NULL;
    htsFile *fragmentsHandle = NULL;
    tbx_t *tbxHandle = NULL;
    ccounts_result result;

    if (hasIndexOut != NULL)
    {
        *hasIndexOut = 0;
    }

    if (ccounts_isAlignmentKind(sourceConfig->sourceKind))
    {
        result = ccounts_openAlignmentFile(sourceConfig, threadCount, &fileHandle, &header);
        if (result.errorCode != 0)
        {
            return result;
        }

        indexHandle = sam_index_load((htsFile *)fileHandle, sourceConfig->path);
        if (indexHandle == NULL && buildIndex)
        {
            ccounts_closeAlignmentFile(fileHandle, header);
            if (sam_index_build3(sourceConfig->path, NULL, 0, threadCount > 0 ? threadCount : 1) < 0)
            {
                return ccounts_makeResult(-1, "failed to build alignment index");
            }
            result = ccounts_openAlignmentFile(sourceConfig, threadCount, &fileHandle, &header);
            if (result.errorCode != 0)
            {
                return result;
            }
            indexHandle = sam_index_load((htsFile *)fileHandle, sourceConfig->path);
        }
        if (hasIndexOut != NULL)
        {
            *hasIndexOut = indexHandle != NULL ? 1 : 0;
        }
        if (indexHandle != NULL)
        {
            hts_idx_destroy(indexHandle);
        }
        ccounts_closeAlignmentFile(fileHandle, header);
        return ccounts_makeOk();
    }

    fragmentsHandle = hts_open(sourceConfig->path, "r");
    if (fragmentsHandle == NULL)
    {
        return ccounts_makeResult(-1, "failed to open fragments source");
    }
    tbxHandle = tbx_index_load(sourceConfig->path);
    if (hasIndexOut != NULL)
    {
        *hasIndexOut = tbxHandle != NULL ? 1 : 0;
    }
    if (tbxHandle != NULL)
    {
        tbx_destroy(tbxHandle);
    }
    hts_close(fragmentsHandle);
    return ccounts_makeOk();
}

ccounts_result ccounts_isPairedEnd(
    const ccounts_sourceConfig *sourceConfig,
    int threadCount,
    int maxReads,
    int *isPairedEndOut)
{
    samFile *fileHandle = NULL;
    sam_hdr_t *header = NULL;
    bam1_t *record = NULL;
    ccounts_result result;
    int seen = 0;

    if (isPairedEndOut != NULL)
    {
        *isPairedEndOut = 0;
    }
    if (!ccounts_isAlignmentKind(sourceConfig->sourceKind))
    {
        return ccounts_makeOk();
    }

    result = ccounts_openAlignmentFile(sourceConfig, threadCount, &fileHandle, &header);
    if (result.errorCode != 0)
    {
        return result;
    }

    record = bam_init1();
    if (record == NULL)
    {
        ccounts_closeAlignmentFile(fileHandle, header);
        return ccounts_makeResult(-1, "failed to allocate bam record");
    }

    while (sam_read1(fileHandle, header, record) >= 0)
    {
        if ((record->core.flag & BAM_FPAIRED) != 0)
        {
            if (isPairedEndOut != NULL)
            {
                *isPairedEndOut = 1;
            }
            break;
        }
        ++seen;
        if (maxReads > 0 && seen >= maxReads)
        {
            break;
        }
    }

    bam_destroy1(record);
    ccounts_closeAlignmentFile(fileHandle, header);
    return ccounts_makeOk();
}

ccounts_result ccounts_getReadLength(
    const ccounts_sourceConfig *sourceConfig,
    int threadCount,
    int minReads,
    int maxIterations,
    int flagExclude,
    uint32_t *readLengthOut)
{
    samFile *fileHandle = NULL;
    sam_hdr_t *header = NULL;
    bam1_t *record = NULL;
    ccounts_result result;
    uint32_t *readLengths = NULL;
    int readCapacity = 0;
    int readCount = 0;
    int iterCount = 0;
    int midIndex = 0;
    hts_pos_t queryLength = 0;

    if (readLengthOut != NULL)
    {
        *readLengthOut = 0;
    }

    if (minReads < 1)
    {
        minReads = 1;
    }
    if (maxIterations < minReads)
    {
        maxIterations = minReads;
    }

    if (!ccounts_isAlignmentKind(sourceConfig->sourceKind))
    {
        // fragments -- FFR: check if necessary
        kstring_t lineBuffer = {0, 0, NULL};
        htsFile *fragmentsHandle = NULL;
        int iteratorCode = 0;
        int fragmentCount = 0;
        int64_t fragStart = 0;
        int64_t fragEnd = 0;
        const char *cursor = NULL;
        const char *fieldStart = NULL;
        size_t fieldLength = 0U;
        uint32_t *fragmentLengths = NULL;

        fragmentsHandle = hts_open(sourceConfig->path, "r");
        if (fragmentsHandle == NULL)
        {
            return ccounts_makeResult(-1, "failed to open fragments source");
        }
        fragmentLengths = (uint32_t *)malloc((size_t)maxIterations * sizeof(uint32_t));
        if (fragmentLengths == NULL)
        {
            hts_close(fragmentsHandle);
            return ccounts_makeResult(-1, "failed to allocate fragments read length buffer");
        }
        while (fragmentCount < maxIterations &&
               (iteratorCode = hts_getline(fragmentsHandle, '\n', &lineBuffer)) >= 0)
        {
            cursor = lineBuffer.s;
            int fieldIndex = 0;
            while (fieldIndex < 3 && cursor != NULL)
            {
                fieldStart = cursor;
                while (*cursor != '\0' && *cursor != '\t' && *cursor != '\n' && *cursor != '\r')
                {
                    ++cursor;
                }
                fieldLength = (size_t)(cursor - fieldStart);
                if (fieldIndex == 1 && !ccounts_parseInt64Field(fieldStart, fieldLength, &fragStart))
                {
                    fieldLength = 0U;
                    break;
                }
                if (fieldIndex == 2 && !ccounts_parseInt64Field(fieldStart, fieldLength, &fragEnd))
                {
                    fieldLength = 0U;
                    break;
                }
                if (*cursor == '\t')
                {
                    ++cursor;
                }
                ++fieldIndex;
            }
            if (fieldLength == 0U || fragEnd <= fragStart)
            {
                continue;
            }
            fragmentLengths[fragmentCount] = (uint32_t)(fragEnd - fragStart);
            ++fragmentCount;
            if (fragmentCount >= minReads)
            {
                break;
            }
        }
        if (lineBuffer.s != NULL)
        {
            free(lineBuffer.s);
        }
        if (fragmentCount > 0)
        {
            int midIndex = 0;
            qsort(fragmentLengths, (size_t)fragmentCount, sizeof(uint32_t), ccounts_compareUint32);
            midIndex = fragmentCount / 2;
            if ((fragmentCount & 1) == 0)
            {
                if (readLengthOut != NULL)
                {
                    *readLengthOut = (uint32_t)((fragmentLengths[midIndex - 1] + fragmentLengths[midIndex]) / 2U);
                }
            }
            else if (readLengthOut != NULL)
            {
                *readLengthOut = fragmentLengths[midIndex];
            }
        }
        free(fragmentLengths);
        hts_close(fragmentsHandle);
        if (fragmentCount == 0)
        {
            return ccounts_makeResult(-1, "failed to estimate fragment length from fragments source");
        }
        return ccounts_makeOk();
    }

    result = ccounts_openAlignmentFile(sourceConfig, threadCount, &fileHandle, &header);
    if (result.errorCode != 0)
    {
        return result;
    }

    record = bam_init1();
    if (record == NULL)
    {
        ccounts_closeAlignmentFile(fileHandle, header);
        return ccounts_makeResult(-1, "failed to allocate bam record");
    }

    readCapacity = maxIterations;
    readLengths = (uint32_t *)malloc((size_t)readCapacity * sizeof(uint32_t));
    if (readLengths == NULL)
    {
        bam_destroy1(record);
        ccounts_closeAlignmentFile(fileHandle, header);
        return ccounts_makeResult(-1, "failed to allocate read length buffer");
    }

    while (iterCount < maxIterations && sam_read1(fileHandle, header, record) >= 0)
    {
        ++iterCount;
        if ((record->core.flag & flagExclude) != 0)
        {
            continue;
        }

        queryLength = record->core.l_qseq;
        if (queryLength <= 0 && record->core.n_cigar > 0)
        {
            queryLength = bam_cigar2qlen(record->core.n_cigar, bam_get_cigar(record));
        }
        if (queryLength <= 0)
        {
            continue;
        }

        readLengths[readCount] = (uint32_t)queryLength;
        ++readCount;
        if (readCount >= minReads)
        {
            break;
        }
    }

    if (readCount > 0)
    {
        qsort(readLengths, (size_t)readCount, sizeof(uint32_t), ccounts_compareUint32);
        midIndex = readCount / 2;
        if ((readCount & 1) == 0)
        {
            if (readLengthOut != NULL)
            {
                *readLengthOut = (uint32_t)((readLengths[midIndex - 1] + readLengths[midIndex]) / 2U);
            }
        }
        else if (readLengthOut != NULL)
        {
            *readLengthOut = readLengths[midIndex];
        }
    }

    free(readLengths);
    bam_destroy1(record);
    ccounts_closeAlignmentFile(fileHandle, header);

    if (readCount == 0)
    {
        return ccounts_makeResult(-1, "failed to estimate read length");
    }
    return ccounts_makeOk();
}

/**
 * @brief estimate the fragment length (PE) or pseudo-fragment length (SE) via strand xcorr
 */
ccounts_result ccounts_getFragmentLength(
    const ccounts_sourceConfig *sourceConfig,
    int threadCount,
    int flagExclude,
    int maxIterations,
    int maxInsertSize,
    int blockSize,
    int rollingChunkSize,
    int lagStep,
    int earlyExit,
    int fallbackLength,
    uint32_t *fragmentLengthOut)
{
    samFile *fileHandle = NULL;
    sam_hdr_t *header = NULL;
    hts_idx_t *indexHandle = NULL;
    hts_itr_t *iteratorHandle = NULL;
    bam1_t *record = NULL;
    ccounts_result result;
    int topTid[3] = {-1, -1, -1};
    uint64_t topLen[3] = {0ULL, 0ULL, 0ULL};
    int tid = 0;
    int topIndex = 0;
    int iterCount = 0;
    int readSampleCount = 0;
    int pairedEnd = 0;
    double readLengthSum = 0.0;
    hts_pos_t queryLength = 0;
    int minInsertSize = 1;
    int requiredSamplesPE = 0;
    uint32_t *templateLengths = NULL;
    int templateCount = 0;
    int midIndex = 0;
    int contigIndex = 0;
    int64_t templateLength = 0;
    int64_t absoluteTemplateLength = 0;
    int numChunks = 0;
    double *rawCounts = NULL;
    double *prefixCounts = NULL;
    double *localDensity = NULL;
    ccounts_scoredIndex *rankedCenters = NULL;
    unsigned char *seen = NULL;
    int *blockCenters = NULL;
    int blockCenterCount = 0;
    int winSize = 0;
    int winHalf = 0;
    int takeK = 0;
    int sortedIndex = 0;
    int startIndex = 0;
    int endIndex = 0;
    int centerIndex = 0;
    int blockHalf = 0;
    int64_t contigLength = 0;
    int64_t blockStartBP = 0;
    int64_t blockEndBP = 0;
    double *fwd = NULL;
    double *rev = NULL;
    double fwdSum = 0.0;
    double revSum = 0.0;
    double fwdMean = 0.0;
    double revMean = 0.0;
    int blockOffset = 0;
    int64_t readStart = 0;
    int64_t readEnd = 0;
    int64_t fivePrimeIndex = 0;
    int maxValidLag = 0;
    int localMinLag = 0;
    int localMaxLag = 0;
    int bestLag = -1;
    double bestScore = 0.0;
    double score = 0.0;
    int lag = 0;
    int signalIndex = 0;
    int blockLength = 0;
    uint32_t *bestLags = NULL;
    int bestLagCapacity = 0;
    int bestLagCount = 0;
    uint32_t medianLag = 0U;

    if (fragmentLengthOut != NULL)
    {
        *fragmentLengthOut = (uint32_t)(fallbackLength > 0 ? fallbackLength : 0);
    }
    if (maxIterations < 1)
    {
        maxIterations = 1;
    }
    if (maxInsertSize < 1)
    {
        maxInsertSize = 1;
    }
    if (blockSize < 64)
    {
        blockSize = 64;
    }
    if (rollingChunkSize < 1)
    {
        rollingChunkSize = 1;
    }
    if (lagStep < 1)
    {
        lagStep = 1;
    }
    if (earlyExit < 1)
    {
        earlyExit = maxIterations;
    }
    if (!ccounts_isAlignmentKind(sourceConfig->sourceKind))
    {
        return ccounts_makeResult(-1, "source kind is not BAM");
    }

    result = ccounts_openAlignmentFile(sourceConfig, threadCount, &fileHandle, &header);
    if (result.errorCode != 0)
    {
        return result;
    }

    indexHandle = sam_index_load((htsFile *)fileHandle, sourceConfig->path);
    if (indexHandle == NULL)
    {
        ccounts_closeAlignmentFile(fileHandle, header);
        return ccounts_makeResult(-1, "alignment index is required for fragment-length estimation");
    }

    record = bam_init1();
    if (record == NULL)
    {
        hts_idx_destroy(indexHandle);
        ccounts_closeAlignmentFile(fileHandle, header);
        return ccounts_makeResult(-1, "failed to allocate bam record");
    }

    /* get 3 largest contigs -- we estimate fragment length using these */
    for (tid = 0; tid < header->n_targets; ++tid)
    {
        uint64_t targetLength = (uint64_t)header->target_len[tid];
        for (topIndex = 0; topIndex < 3; ++topIndex)
        {
            if (targetLength > topLen[topIndex])
            {
                int shiftIndex = 0;
                for (shiftIndex = 2; shiftIndex > topIndex; --shiftIndex)
                {
                    topLen[shiftIndex] = topLen[shiftIndex - 1];
                    topTid[shiftIndex] = topTid[shiftIndex - 1];
                }
                topLen[topIndex] = targetLength;
                topTid[topIndex] = tid;
                break;
            }
        }
    }

    for (contigIndex = 0; contigIndex < 3 && readSampleCount < maxIterations; ++contigIndex)
    {
        if (topTid[contigIndex] < 0 || topLen[contigIndex] == 0ULL)
        {
            continue;
        }
        iteratorHandle = sam_itr_queryi(
            indexHandle,
            topTid[contigIndex],
            0,
            (hts_pos_t)topLen[contigIndex]);
        if (iteratorHandle == NULL)
        {
            continue;
        }
        while (readSampleCount < maxIterations &&
               sam_itr_next((htsFile *)fileHandle, iteratorHandle, record) >= 0)
        {
            if ((record->core.flag & flagExclude) != 0)
            {
                continue;
            }
            if ((record->core.flag & BAM_FUNMAP) != 0)
            {
                continue;
            }
            if (!pairedEnd && (record->core.flag & BAM_FPAIRED) != 0)
            {
                pairedEnd = 1;
            }

            queryLength = record->core.l_qseq;
            if (queryLength <= 0 && record->core.n_cigar > 0)
            {
                queryLength = bam_cigar2qlen(record->core.n_cigar, bam_get_cigar(record));
            }
            if (queryLength <= 0)
            {
                continue;
            }
            readLengthSum += (double)queryLength;
            ++readSampleCount;
        }
        hts_itr_destroy(iteratorHandle);
        iteratorHandle = NULL;
    }

    if (readSampleCount <= 0)
    {
        if (fragmentLengthOut != NULL)
        {
            *fragmentLengthOut = (uint32_t)(fallbackLength > 0 ? fallbackLength : 0);
        }
        bam_destroy1(record);
        hts_idx_destroy(indexHandle);
        ccounts_closeAlignmentFile(fileHandle, header);
        return ccounts_makeOk();
    }

    minInsertSize = (int)(readLengthSum / (double)readSampleCount);
    if (minInsertSize < 1)
    {
        minInsertSize = 1;
    }
    if (minInsertSize > maxInsertSize)
    {
        minInsertSize = maxInsertSize;
    }

    if (pairedEnd)
    {
        requiredSamplesPE = maxIterations;
        if (requiredSamplesPE < 2000)
        {
            requiredSamplesPE = 2000;
        }
        templateLengths = (uint32_t *)malloc((size_t)requiredSamplesPE * sizeof(uint32_t));
        if (templateLengths == NULL)
        {
            bam_destroy1(record);
            hts_idx_destroy(indexHandle);
            ccounts_closeAlignmentFile(fileHandle, header);
            return ccounts_makeResult(-1, "failed to allocate fragment length buffer");
        }

        for (contigIndex = 0; contigIndex < 3 && templateCount < requiredSamplesPE; ++contigIndex)
        {
            if (topTid[contigIndex] < 0 || topLen[contigIndex] == 0ULL)
            {
                continue;
            }
            iteratorHandle = sam_itr_queryi(
                indexHandle,
                topTid[contigIndex],
                0,
                (hts_pos_t)topLen[contigIndex]);
            if (iteratorHandle == NULL)
            {
                continue;
            }
            while (templateCount < requiredSamplesPE &&
                   sam_itr_next((htsFile *)fileHandle, iteratorHandle, record) >= 0)
            {
                if ((record->core.flag & flagExclude) != 0)
                {
                    continue;
                }
                if ((record->core.flag & BAM_FPROPER_PAIR) == 0)
                {
                    continue;
                }
                if ((record->core.flag & BAM_FREAD2) != 0)
                {
                    continue;
                }
                if ((record->core.flag & BAM_FMUNMAP) != 0 || record->core.mtid != record->core.tid)
                {
                    continue;
                }

                templateLength = (int64_t)record->core.isize;
                absoluteTemplateLength = templateLength >= 0 ? templateLength : -templateLength;
                if (absoluteTemplateLength < (int64_t)minInsertSize ||
                    absoluteTemplateLength > (int64_t)maxInsertSize)
                {
                    continue;
                }
                templateLengths[templateCount] = (uint32_t)absoluteTemplateLength;
                ++templateCount;
            }
            hts_itr_destroy(iteratorHandle);
            iteratorHandle = NULL;
        }

        if (templateCount > 0)
        {
            qsort(templateLengths, (size_t)templateCount, sizeof(uint32_t), ccounts_compareUint32);
            midIndex = templateCount / 2;
            if ((templateCount & 1) == 0)
            {
                medianLag = (uint32_t)((templateLengths[midIndex - 1] + templateLengths[midIndex]) / 2U);
            }
            else
            {
                medianLag = templateLengths[midIndex];
            }
            if (medianLag < (uint32_t)minInsertSize)
            {
                medianLag = (uint32_t)minInsertSize;
            }
            if (medianLag > (uint32_t)maxInsertSize)
            {
                medianLag = (uint32_t)maxInsertSize;
            }
            if (fragmentLengthOut != NULL)
            {
                *fragmentLengthOut = medianLag;
            }
        }

        free(templateLengths);
        bam_destroy1(record);
        hts_idx_destroy(indexHandle);
        ccounts_closeAlignmentFile(fileHandle, header);
        return ccounts_makeOk();
    }

    bestLagCapacity = earlyExit;
    if (bestLagCapacity < 3)
    {
        bestLagCapacity = 3;
    }
    bestLags = (uint32_t *)malloc((size_t)bestLagCapacity * sizeof(uint32_t));
    fwd = (double *)malloc((size_t)blockSize * sizeof(double));
    rev = (double *)malloc((size_t)blockSize * sizeof(double));
    if (bestLags == NULL || fwd == NULL || rev == NULL)
    {
        free(bestLags);
        free(fwd);
        free(rev);
        bam_destroy1(record);
        hts_idx_destroy(indexHandle);
        ccounts_closeAlignmentFile(fileHandle, header);
        return ccounts_makeResult(-1, "failed to allocate single-end fragment length buffers");
    }

    blockHalf = blockSize / 2;

    /* now we run a rolling window across the largest contigs, looking for regions with high read density,
     *  and compute strand xcorr to find good candidate lags */
    for (contigIndex = 0; contigIndex < 3 && bestLagCount < earlyExit; ++contigIndex)
    {
        if (topTid[contigIndex] < 0 || topLen[contigIndex] == 0ULL)
        {
            continue;
        }
        contigLength = (int64_t)topLen[contigIndex];
        if (contigLength < (int64_t)blockSize)
        {
            continue;
        }

        numChunks = (int)((contigLength + (int64_t)rollingChunkSize - 1) / (int64_t)rollingChunkSize);
        if (numChunks < 1)
        {
            continue;
        }

        rawCounts = (double *)calloc((size_t)numChunks, sizeof(double));
        prefixCounts = (double *)calloc((size_t)(numChunks + 1), sizeof(double));
        localDensity = (double *)calloc((size_t)numChunks, sizeof(double));
        rankedCenters = (ccounts_scoredIndex *)malloc((size_t)numChunks * sizeof(ccounts_scoredIndex));
        seen = (unsigned char *)calloc((size_t)numChunks, sizeof(unsigned char));
        blockCenters = (int *)malloc((size_t)maxIterations * sizeof(int));
        if (rawCounts == NULL || prefixCounts == NULL || localDensity == NULL ||
            rankedCenters == NULL || seen == NULL || blockCenters == NULL)
        {
            free(rawCounts);
            free(prefixCounts);
            free(localDensity);
            free(rankedCenters);
            free(seen);
            free(blockCenters);
            free(bestLags);
            free(fwd);
            free(rev);
            bam_destroy1(record);
            hts_idx_destroy(indexHandle);
            ccounts_closeAlignmentFile(fileHandle, header);
            return ccounts_makeResult(-1, "failed to allocate single-end rolling windows");
        }

        iteratorHandle = sam_itr_queryi(
            indexHandle,
            topTid[contigIndex],
            0,
            (hts_pos_t)contigLength);
        if (iteratorHandle != NULL)
        {
            while (sam_itr_next((htsFile *)fileHandle, iteratorHandle, record) >= 0)
            {
                if ((record->core.flag & flagExclude) != 0)
                {
                    continue;
                }
                if ((record->core.flag & BAM_FUNMAP) != 0)
                {
                    continue;
                }
                blockOffset = (int)((int64_t)record->core.pos / (int64_t)rollingChunkSize);
                if (blockOffset >= 0 && blockOffset < numChunks)
                {
                    rawCounts[blockOffset] += 1.0;
                }
            }
            hts_itr_destroy(iteratorHandle);
            iteratorHandle = NULL;
        }

        winSize = blockSize / rollingChunkSize;
        if (winSize < 1)
        {
            winSize = 1;
        }
        if ((winSize & 1) == 0)
        {
            winSize += 1;
        }
        winHalf = winSize / 2;

        prefixCounts[0] = 0.0;
        for (sortedIndex = 0; sortedIndex < numChunks; ++sortedIndex)
        {
            prefixCounts[sortedIndex + 1] = prefixCounts[sortedIndex] + rawCounts[sortedIndex];
        }
        for (sortedIndex = 0; sortedIndex < numChunks; ++sortedIndex)
        {
            startIndex = sortedIndex - winHalf;
            endIndex = startIndex + winSize;
            if (startIndex < 0)
            {
                startIndex = 0;
                endIndex = winSize < numChunks ? winSize : numChunks;
            }
            if (endIndex > numChunks)
            {
                endIndex = numChunks;
                startIndex = endIndex - winSize;
                if (startIndex < 0)
                {
                    startIndex = 0;
                }
            }
            localDensity[sortedIndex] = prefixCounts[endIndex] - prefixCounts[startIndex];
            rankedCenters[sortedIndex].value = localDensity[sortedIndex];
            rankedCenters[sortedIndex].index = sortedIndex;
        }

        qsort(rankedCenters, (size_t)numChunks, sizeof(ccounts_scoredIndex), ccounts_compareScoredIndexDesc);
        takeK = maxIterations < numChunks ? maxIterations : numChunks;
        blockCenterCount = 0;
        for (sortedIndex = 0; sortedIndex < numChunks && blockCenterCount < takeK; ++sortedIndex)
        {
            centerIndex = rankedCenters[sortedIndex].index;
            if (rankedCenters[sortedIndex].value <= 0.0 || seen[centerIndex] != 0U)
            {
                continue;
            }
            blockCenters[blockCenterCount] = centerIndex;
            ++blockCenterCount;
            startIndex = centerIndex - winHalf;
            endIndex = startIndex + winSize;
            if (startIndex < 0)
            {
                startIndex = 0;
            }
            if (endIndex > numChunks)
            {
                endIndex = numChunks;
            }
            for (iterCount = startIndex; iterCount < endIndex; ++iterCount)
            {
                seen[iterCount] = 1U;
            }
        }

        for (sortedIndex = 0; sortedIndex < blockCenterCount && bestLagCount < earlyExit; ++sortedIndex)
        {
            centerIndex = blockCenters[sortedIndex];
            blockStartBP = (int64_t)centerIndex * (int64_t)rollingChunkSize +
                           ((int64_t)rollingChunkSize / 2) - (int64_t)blockHalf;
            if (blockStartBP < 0)
            {
                blockStartBP = 0;
            }
            blockEndBP = blockStartBP + (int64_t)blockSize;
            if (blockEndBP > contigLength)
            {
                blockEndBP = contigLength;
                blockStartBP = blockEndBP - (int64_t)blockSize;
                if (blockStartBP < 0)
                {
                    continue;
                }
            }

            memset(fwd, 0, (size_t)blockSize * sizeof(double));
            memset(rev, 0, (size_t)blockSize * sizeof(double));

            iteratorHandle = sam_itr_queryi(
                indexHandle,
                topTid[contigIndex],
                (hts_pos_t)blockStartBP,
                (hts_pos_t)blockEndBP);
            if (iteratorHandle == NULL)
            {
                continue;
            }
            while (sam_itr_next((htsFile *)fileHandle, iteratorHandle, record) >= 0)
            {
                if ((record->core.flag & flagExclude) != 0)
                {
                    continue;
                }
                if ((record->core.flag & BAM_FUNMAP) != 0)
                {
                    continue;
                }

                readStart = (int64_t)record->core.pos;
                readEnd = (int64_t)bam_endpos(record);
                if (readEnd <= readStart)
                {
                    continue;
                }
                if (readStart < blockStartBP || readEnd > blockEndBP)
                {
                    continue;
                }

                if ((record->core.flag & BAM_FREVERSE) == 0)
                {
                    blockOffset = (int)(readStart - blockStartBP);
                    if (blockOffset >= 0 && blockOffset < blockSize)
                    {
                        fwd[blockOffset] += 1.0;
                    }
                }
                else
                {
                    fivePrimeIndex = (readEnd - 1) - blockStartBP;
                    if (fivePrimeIndex >= 0 && fivePrimeIndex < (int64_t)blockSize)
                    {
                        rev[(int)fivePrimeIndex] += 1.0;
                    }
                }
            }
            hts_itr_destroy(iteratorHandle);
            iteratorHandle = NULL;

            fwdSum = 0.0;
            revSum = 0.0;
            for (blockOffset = 0; blockOffset < blockSize; ++blockOffset)
            {
                fwdSum += fwd[blockOffset];
                revSum += rev[blockOffset];
            }
            if (fwdSum < 10.0 || revSum < 10.0)
            {
                continue;
            }

            fwdMean = fwdSum / (double)blockSize;
            revMean = revSum / (double)blockSize;
            for (blockOffset = 0; blockOffset < blockSize; ++blockOffset)
            {
                fwd[blockOffset] -= fwdMean;
                rev[blockOffset] -= revMean;
            }

            maxValidLag = maxInsertSize < (blockSize - 1) ? maxInsertSize : (blockSize - 1);
            localMinLag = minInsertSize;
            localMaxLag = maxValidLag;
            if (localMaxLag < localMinLag)
            {
                continue;
            }

            bestLag = -1;
            bestScore = 0.0;
            for (lag = localMinLag; lag <= localMaxLag; lag += lagStep)
            {
                blockLength = blockSize - lag;
                if (blockLength <= 0)
                {
                    continue;
                }
                score = 0.0;
                for (signalIndex = 0; signalIndex < blockLength; ++signalIndex)
                {
                    score += fwd[signalIndex] * rev[signalIndex + lag];
                }
                if (bestLag < 0 || score > bestScore)
                {
                    bestScore = score;
                    bestLag = lag;
                }
            }

            if (bestLag > 0 && bestScore != 0.0)
            {
                bestLags[bestLagCount] = (uint32_t)(bestLag + 1);
                ++bestLagCount;
            }
        }

        free(rawCounts);
        free(prefixCounts);
        free(localDensity);
        free(rankedCenters);
        free(seen);
        free(blockCenters);
        rawCounts = NULL;
        prefixCounts = NULL;
        localDensity = NULL;
        rankedCenters = NULL;
        seen = NULL;
        blockCenters = NULL;
    }

    if (bestLagCount > 0)
    {
        qsort(bestLags, (size_t)bestLagCount, sizeof(uint32_t), ccounts_compareUint32);
        midIndex = bestLagCount / 2;
        if ((bestLagCount & 1) == 0)
        {
            medianLag = (uint32_t)((bestLags[midIndex - 1] + bestLags[midIndex]) / 2U);
        }
        else
        {
            medianLag = bestLags[midIndex];
        }
        if (medianLag < (uint32_t)minInsertSize)
        {
            medianLag = (uint32_t)minInsertSize;
        }
        if (medianLag > (uint32_t)maxInsertSize)
        {
            medianLag = (uint32_t)maxInsertSize;
        }
        if (fragmentLengthOut != NULL)
        {
            *fragmentLengthOut = medianLag;
        }
    }

    free(rawCounts);
    free(prefixCounts);
    free(localDensity);
    free(rankedCenters);
    free(seen);
    free(blockCenters);
    free(bestLags);
    free(fwd);
    free(rev);
    bam_destroy1(record);
    hts_idx_destroy(indexHandle);
    ccounts_closeAlignmentFile(fileHandle, header);
    return ccounts_makeOk();
}

/**
 * @brief get the start and end coordinates of the region covered by alignments/fragments on a contig
 */
ccounts_result ccounts_getChromRange(
    const ccounts_sourceConfig *sourceConfig,
    const char *chromosome,
    uint64_t chromLength,
    int threadCount,
    int flagExclude,
    uint64_t *startOut,
    uint64_t *endOut)
{
    ccounts_sourceHandle *sourceHandle = NULL;
    ccounts_result result;
    bam1_t *record = NULL;
    hts_itr_t *iteratorHandle = NULL;
    int32_t tid = -1;
    uint64_t tailStart = 0;

    if (startOut != NULL)
    {
        *startOut = 0;
    }
    if (endOut != NULL)
    {
        *endOut = 0;
    }

    result = ccounts_openSource(sourceConfig, &sourceHandle);
    if (result.errorCode != 0)
    {
        return result;
    }
    if (sourceHandle == NULL)
    {
        return ccounts_makeResult(-1, "failed to initialize source handle");
    }
    if (sourceHandle->sourceKind == ccounts_sourceKindFragments)
    {
        kstring_t lineBuffer = {0, 0, NULL};
        char regionBuffer[1024];
        hts_itr_t *iteratorHandle = NULL;
        int64_t fragStart = 0;
        int64_t fragEnd = 0;
        int iteratorCode = 0;
        int seenFragment = 0;

        if (sourceHandle->tbxHandle == NULL || sourceHandle->fragmentsHandle == NULL)
        {
            ccounts_closeSource(sourceHandle);
            return ccounts_makeResult(-1, "fragments index is required for chromosome range queries");
        }

        snprintf(regionBuffer, sizeof(regionBuffer), "%s:%d-%llu", chromosome, 1, (unsigned long long)chromLength);
        iteratorHandle = tbx_itr_querys(sourceHandle->tbxHandle, regionBuffer);
        if (iteratorHandle != NULL)
        {
            while ((iteratorCode = tbx_itr_next(sourceHandle->fragmentsHandle, sourceHandle->tbxHandle, iteratorHandle, &lineBuffer)) >= 0)
            {
                const char *cursor = lineBuffer.s;
                const char *fieldStart = NULL;
                const char *fieldEnd = NULL;
                size_t fieldLength = 0U;
                int fieldIndex = 0;
                while (fieldIndex < 3)
                {
                    fieldStart = cursor;
                    while (*cursor != '\0' && *cursor != '\t' && *cursor != '\n' && *cursor != '\r')
                    {
                        ++cursor;
                    }
                    fieldEnd = cursor;
                    fieldLength = (size_t)(fieldEnd - fieldStart);
                    if (fieldIndex == 1 && !ccounts_parseInt64Field(fieldStart, fieldLength, &fragStart))
                    {
                        fieldLength = 0U;
                        break;
                    }
                    if (fieldIndex == 2 && !ccounts_parseInt64Field(fieldStart, fieldLength, &fragEnd))
                    {
                        fieldLength = 0U;
                        break;
                    }
                    if (*cursor == '\t')
                    {
                        ++cursor;
                    }
                    ++fieldIndex;
                }
                if (fieldLength == 0U || fragEnd <= fragStart)
                {
                    continue;
                }
                if (!seenFragment)
                {
                    if (startOut != NULL)
                    {
                        *startOut = (uint64_t)fragStart;
                    }
                    seenFragment = 1;
                }
                if (endOut != NULL)
                {
                    *endOut = (uint64_t)fragEnd;
                }
            }
            hts_itr_destroy(iteratorHandle);
        }
        if (lineBuffer.s != NULL)
        {
            free(lineBuffer.s);
        }
        ccounts_closeSource(sourceHandle);
        return ccounts_makeOk();
    }
    if (sourceHandle->indexHandle == NULL)
    {
        ccounts_closeSource(sourceHandle);
        return ccounts_makeResult(-1, "alignment index is required for chromosome range queries");
    }

    if (threadCount > 1)
    {
        hts_set_threads((htsFile *)sourceHandle->fileHandle, threadCount);
    }

    tid = sam_hdr_name2tid(sourceHandle->header, chromosome);
    if (tid < 0)
    {
        ccounts_closeSource(sourceHandle);
        return ccounts_makeResult(-1, "chromosome not found in alignment header");
    }

    record = bam_init1();
    if (record == NULL)
    {
        ccounts_closeSource(sourceHandle);
        return ccounts_makeResult(-1, "failed to allocate bam record");
    }

    iteratorHandle = sam_itr_queryi(sourceHandle->indexHandle, tid, 0, (hts_pos_t)chromLength);
    if (iteratorHandle != NULL)
    {
        while (sam_itr_next((htsFile *)sourceHandle->fileHandle, iteratorHandle, record) >= 0)
        {
            if ((record->core.flag & flagExclude) != 0)
            {
                continue;
            }
            if (startOut != NULL)
            {
                *startOut = (uint64_t)record->core.pos;
            }
            break;
        }
        hts_itr_destroy(iteratorHandle);
    }
    /* start looking for end of effective region 2mb from contig bound*/
    const uint64_t tailCushion = 2000000ULL;
    tailStart = chromLength > tailCushion ? chromLength - tailCushion : 0;
    iteratorHandle = sam_itr_queryi(
        sourceHandle->indexHandle,
        tid,
        (hts_pos_t)tailStart,
        (hts_pos_t)chromLength);
    if (iteratorHandle != NULL)
    {
        while (sam_itr_next((htsFile *)sourceHandle->fileHandle, iteratorHandle, record) >= 0)
        {
            if ((record->core.flag & flagExclude) != 0)
            {
                continue;
            }
            if (endOut != NULL)
            {
                *endOut = (uint64_t)bam_endpos(record);
            }
        }
        hts_itr_destroy(iteratorHandle);
    }

    bam_destroy1(record);
    ccounts_closeSource(sourceHandle);
    return ccounts_makeOk();
}

ccounts_result ccounts_getMappedReadCount(
    const ccounts_sourceConfig *sourceConfig,
    int threadCount,
    const char *const *excludeChromosomes,
    int excludeChromosomeCount,
    uint8_t countMode,
    uint8_t oneReadPerBin,
    uint64_t *mappedReadCountOut,
    uint64_t *unmappedReadCountOut)
{
    ccounts_sourceHandle *sourceHandle = NULL;
    ccounts_result result;
    int targetCount = 0;
    int targetIndex = 0;
    uint64_t mappedCount = 0;
    uint64_t unmappedCount = 0;
    uint64_t mappedLocal = 0;
    uint64_t unmappedLocal = 0;
    const char *targetName = NULL;

    if (mappedReadCountOut != NULL)
    {
        *mappedReadCountOut = 0;
    }
    if (unmappedReadCountOut != NULL)
    {
        *unmappedReadCountOut = 0;
    }

    result = ccounts_openSource(sourceConfig, &sourceHandle);
    if (result.errorCode != 0)
    {
        return result;
    }
    if (sourceHandle == NULL)
    {
        return ccounts_makeResult(-1, "failed to initialize source handle");
    }
    /* we don't use index metadata for fragment/barcode sources */
    if (sourceHandle->sourceKind == ccounts_sourceKindFragments)
    {
        kstring_t lineBuffer = {0, 0, NULL};
        int iteratorCode = 0;
        uint64_t mappedCount = 0;
        int64_t fragmentCount = 1;
        uint64_t emittedCount = 0;
        const char *cursor = NULL;
        const char *fieldStart = NULL;
        const char *fieldEnd = NULL;
        size_t fieldLength = 0U;
        int fieldIndex = 0;
        const char *barcodeStart = NULL;
        size_t barcodeLength = 0U;
        if (sourceHandle->fragmentsHandle == NULL)
        {
            ccounts_closeSource(sourceHandle);
            return ccounts_makeResult(-1, "failed to initialize fragments handle");
        }
        while ((iteratorCode = hts_getline(sourceHandle->fragmentsHandle, '\n', &lineBuffer)) >= 0)
        {
            cursor = lineBuffer.s;
            fieldIndex = 0;
            barcodeStart = NULL;
            barcodeLength = 0U;
            fragmentCount = 1;
            while (fieldIndex < 5 && cursor != NULL)
            {
                fieldStart = cursor;
                while (*cursor != '\0' && *cursor != '\t' && *cursor != '\n' && *cursor != '\r')
                {
                    ++cursor;
                }
                fieldEnd = cursor;
                fieldLength = (size_t)(fieldEnd - fieldStart);
                if (fieldIndex == 0)
                {
                    if (ccounts_stringFieldInList(fieldStart, fieldLength, excludeChromosomes, excludeChromosomeCount))
                    {
                        fieldLength = 0U;
                        break;
                    }
                }
                else if (fieldIndex == 4)
                {
                    if (!ccounts_parseInt64Field(fieldStart, fieldLength, &fragmentCount))
                    {
                        fragmentCount = 1;
                    }
                }
                else if (fieldIndex == 3)
                {
                    barcodeStart = fieldStart;
                    barcodeLength = fieldLength;
                }
                if (*cursor == '\t')
                {
                    ++cursor;
                }
                ++fieldIndex;
            }
            if (fieldLength == 0U)
            {
                continue;
            }
            if (barcodeLength > 0U &&
                !ccounts_barcodeAllowed(
                    sourceHandle->barcodeAllowList,
                    sourceHandle->barcodeAllowCount,
                    barcodeStart,
                    barcodeLength))
            {
                continue;
            }
            // if the fragment count is missing or malformed, assume it's a single read
            // if its present, use it to scale the read count
            emittedCount = (uint64_t)(fragmentCount > 0 ? fragmentCount : 1);
            if (!oneReadPerBin &&
                ((ccounts_countMode)countMode == ccounts_countModeCutSite ||
                 (ccounts_countMode)countMode == ccounts_countModeFivePrime))
            {
                emittedCount *= 2U;
            }
            mappedCount += emittedCount;
        }
        if (mappedReadCountOut != NULL)
        {
            *mappedReadCountOut = mappedCount;
        }
        if (lineBuffer.s != NULL)
        {
            free(lineBuffer.s);
        }
        ccounts_closeSource(sourceHandle);
        return ccounts_makeOk();
    }
    if (sourceHandle->indexHandle == NULL)
    {
        ccounts_closeSource(sourceHandle);
        return ccounts_makeResult(-1, "alignment index is required for mapped read counts");
    }

    if (threadCount > 1)
    {
        hts_set_threads((htsFile *)sourceHandle->fileHandle, threadCount);
    }

    targetCount = sam_hdr_nref(sourceHandle->header);
    for (targetIndex = 0; targetIndex < targetCount; ++targetIndex)
    {
        mappedLocal = 0;
        unmappedLocal = 0;
        if (hts_idx_get_stat(sourceHandle->indexHandle, targetIndex, &mappedLocal, &unmappedLocal) != 0)
        {
            continue;
        }
        targetName = sam_hdr_tid2name(sourceHandle->header, targetIndex);
        if (ccounts_stringInList(targetName, excludeChromosomes, excludeChromosomeCount))
        {
            continue;
        }
        mappedCount += mappedLocal;
        unmappedCount += unmappedLocal;
    }
    unmappedCount += hts_idx_get_n_no_coor(sourceHandle->indexHandle);

    if (mappedReadCountOut != NULL)
    {
        *mappedReadCountOut = mappedCount;
    }
    if (unmappedReadCountOut != NULL)
    {
        *unmappedReadCountOut = unmappedCount;
    }

    ccounts_closeSource(sourceHandle);
    return ccounts_makeOk();
}

/**
 * @brief get the number of unique cell barcodes observed in a tabix/fragments file
 */
ccounts_result ccounts_getCellCount(
    const ccounts_sourceConfig *sourceConfig,
    uint64_t *cellCountOut)
{
    ccounts_sourceHandle *sourceHandle = NULL;
    ccounts_result result;
    kstring_t lineBuffer = {0, 0, NULL};
    khash_t(ccounts_barcodeSet) *barcodeSet = NULL;
    int iteratorCode = 0;
    const char *cursor = NULL;
    const char *fieldStart = NULL;
    size_t fieldLength = 0U;
    int fieldIndex = 0;
    const char *barcodeStart = NULL;
    size_t barcodeLength = 0U;

    if (cellCountOut != NULL)
    {
        *cellCountOut = 0;
    }

    result = ccounts_openSource(sourceConfig, &sourceHandle);
    if (result.errorCode != 0)
    {
        return result;
    }
    if (sourceHandle == NULL)
    {
        return ccounts_makeResult(-1, "failed to initialize source handle");
    }
    if (sourceHandle->sourceKind != ccounts_sourceKindFragments ||
        sourceHandle->fragmentsHandle == NULL)
    {
        ccounts_closeSource(sourceHandle);
        return ccounts_makeResult(-1, "cell count is only supported for fragments sources");
    }

    barcodeSet = kh_init(ccounts_barcodeSet);
    if (barcodeSet == NULL)
    {
        ccounts_closeSource(sourceHandle);
        return ccounts_makeResult(-1, "failed to initialize barcode set");
    }
    // iterate through the fragments file, collect specified barcodes
    while ((iteratorCode = hts_getline(sourceHandle->fragmentsHandle, '\n', &lineBuffer)) >= 0)
    {
        cursor = lineBuffer.s;
        fieldIndex = 0;
        barcodeStart = NULL;
        barcodeLength = 0U;
        while (*cursor != '\0' && *cursor != '\n' && *cursor != '\r')
        {
            fieldStart = cursor;
            while (*cursor != '\0' && *cursor != '\t' && *cursor != '\n' && *cursor != '\r')
            {
                ++cursor;
            }
            fieldLength = (size_t)(cursor - fieldStart);
            if (fieldIndex == 3)
            {
                barcodeStart = fieldStart;
                barcodeLength = fieldLength;
                break;
            }
            // move after hitting barcode
            if (*cursor == '\t')
            {
                ++cursor;
            }
            ++fieldIndex;
        }
        if (barcodeLength == 0U)
        {
            continue;
        }
        if (!ccounts_barcodeAllowed(
                sourceHandle->barcodeAllowList,
                sourceHandle->barcodeAllowCount,
                barcodeStart,
                barcodeLength))
        {
            continue;
        }
        {
            char *barcodeCopy = ccounts_copyStringField(barcodeStart, barcodeLength);
            /* hash for O(1) lookup */
            khiter_t barcodeIndex;
            int absent = 0;
            if (barcodeCopy == NULL)
            {
                if (lineBuffer.s != NULL)
                {
                    free(lineBuffer.s);
                }
                for (barcodeIndex = kh_begin(barcodeSet); barcodeIndex != kh_end(barcodeSet); ++barcodeIndex)
                {
                    if (kh_exist(barcodeSet, barcodeIndex))
                    {
                        free((char *)kh_key(barcodeSet, barcodeIndex));
                    }
                }
                kh_destroy(ccounts_barcodeSet, barcodeSet);
                ccounts_closeSource(sourceHandle);
                return ccounts_makeResult(-1, "failed to store barcode value");
            }
            // keep only distinct allowed barcodes
            barcodeIndex = kh_put(ccounts_barcodeSet, barcodeSet, barcodeCopy, &absent);
            if (absent < 0)
            {
                free(barcodeCopy);
                if (lineBuffer.s != NULL)
                {
                    free(lineBuffer.s);
                }
                for (barcodeIndex = kh_begin(barcodeSet); barcodeIndex != kh_end(barcodeSet); ++barcodeIndex)
                {
                    if (kh_exist(barcodeSet, barcodeIndex))
                    {
                        free((char *)kh_key(barcodeSet, barcodeIndex));
                    }
                }
                kh_destroy(ccounts_barcodeSet, barcodeSet);
                ccounts_closeSource(sourceHandle);
                return ccounts_makeResult(-1, "failed to add barcode to set");
            }
            if (absent == 0)
            {
                free(barcodeCopy);
            }
        }
    }

    if (cellCountOut != NULL)
    {
        *cellCountOut = (uint64_t)kh_size(barcodeSet);
    }
    if (lineBuffer.s != NULL)
    {
        free(lineBuffer.s);
    }
    if (barcodeSet != NULL)
    {
        khiter_t barcodeIndex;
        for (barcodeIndex = kh_begin(barcodeSet); barcodeIndex != kh_end(barcodeSet); ++barcodeIndex)
        {
            if (kh_exist(barcodeSet, barcodeIndex))
            {
                free((char *)kh_key(barcodeSet, barcodeIndex));
            }
        }
        kh_destroy(ccounts_barcodeSet, barcodeSet);
    }
    ccounts_closeSource(sourceHandle);
    return ccounts_makeOk();
}

ccounts_result ccounts_openSource(
    const ccounts_sourceConfig *sourceConfig,
    ccounts_sourceHandle **sourceHandleOut)
{
    ccounts_sourceHandle *sourceHandle = NULL;
    ccounts_result result;

    if (sourceHandleOut != NULL)
    {
        *sourceHandleOut = NULL;
    }

    sourceHandle = (ccounts_sourceHandle *)calloc(1, sizeof(ccounts_sourceHandle));
    if (sourceHandle == NULL)
    {
        return ccounts_makeResult(-1, "failed to allocate source handle");
    }
    sourceHandle->sourceKind = sourceConfig->sourceKind;

    // open once here so repeated region work can reuse the same handle
    if (ccounts_isAlignmentKind(sourceConfig->sourceKind))
    {
        result = ccounts_openAlignmentFile(sourceConfig, 0, &sourceHandle->fileHandle, &sourceHandle->header);
        if (result.errorCode != 0)
        {
            free(sourceHandle);
            return result;
        }
        sourceHandle->indexHandle = sam_index_load((htsFile *)sourceHandle->fileHandle, sourceConfig->path);
    }
    else
    {
        result = ccounts_openFragmentsSource(sourceConfig, sourceHandle);
        if (result.errorCode != 0)
        {
            free(sourceHandle);
            return result;
        }
    }

    if (sourceHandleOut != NULL)
    {
        *sourceHandleOut = sourceHandle;
    }
    return ccounts_makeOk();
}

void ccounts_closeSource(ccounts_sourceHandle *sourceHandle)
{
    if (sourceHandle == NULL)
    {
        return;
    }
    if (sourceHandle->indexHandle != NULL)
    {
        hts_idx_destroy(sourceHandle->indexHandle);
    }
    if (sourceHandle->tbxHandle != NULL)
    {
        tbx_destroy(sourceHandle->tbxHandle);
    }
    if (sourceHandle->fragmentsHandle != NULL)
    {
        hts_close(sourceHandle->fragmentsHandle);
    }
    ccounts_freeBarcodeAllowList(
        sourceHandle->barcodeAllowList,
        sourceHandle->barcodeAllowCount);
    ccounts_closeAlignmentFile(sourceHandle->fileHandle, sourceHandle->header);
    free(sourceHandle);
}

/**
 * @ brief count coverage for a specified region and options, writing to the pre-allocated buffer
 */
ccounts_result ccounts_countRegion(
    ccounts_sourceHandle *sourceHandle,
    const ccounts_region *region,
    const ccounts_countOptions *countOptions,
    float *countBuffer,
    size_t countBufferLength)
{
    bam1_t *record = NULL;
    hts_itr_t *iteratorHandle = NULL;
    float *deltaBuffer = NULL;
    int32_t tid = -1;
    size_t intervalIndex = 0;
    float deltaValue = 0.0f;
    int64_t start64 = 0;
    int64_t end64 = 0;
    int64_t step64 = 0;
    int64_t readStart = 0;
    int64_t readEnd = 0;
    int64_t adjStart = 0;
    int64_t adjEnd = 0;
    int64_t fivePrime = 0;
    int64_t minTemplateLength = 0;
    int64_t templateLength = 0;
    int64_t absoluteTemplateLength = 0;
    int64_t midPoint = 0;
    size_t index0 = 0;
    size_t index1 = 0;

    if (sourceHandle == NULL || region == NULL || countOptions == NULL || countBuffer == NULL)
    {
        return ccounts_makeResult(-1, "count request is invalid");
    }
    if (sourceHandle->sourceKind == ccounts_sourceKindFragments)
    {
        hts_itr_t *fragmentsIterator = NULL;
        kstring_t lineBuffer = {0, 0, NULL};
        char *regionString = NULL;
        size_t regionStringLength = 0U;
        int iteratorCode = 0;
        const char *cursor = NULL;
        const char *fieldStart = NULL;
        size_t fieldLength = 0U;
        int fieldIndex = 0;
        int64_t fragStart = 0;
        int64_t fragEnd = 0;
        int64_t fragCount = 1;
        const char *barcodeStart = NULL;
        size_t barcodeLength = 0U;
        int64_t cutPosition = 0;
        float incrementValue = 1.0f;

        if (sourceHandle->tbxHandle == NULL || sourceHandle->fragmentsHandle == NULL)
        {
            return ccounts_makeResult(-1, "fragments source handle is not initialized");
        }

        start64 = (int64_t)region->start;
        end64 = (int64_t)region->end;
        step64 = (int64_t)region->intervalSizeBP;

        regionStringLength = strlen(region->chromosome) + 64U;
        regionString = (char *)malloc(regionStringLength);
        if (regionString == NULL)
        {
            return ccounts_makeResult(-1, "failed to allocate fragments region string");
        }
        // tabix queries use 1-based closed coordinates
        snprintf(
            regionString,
            regionStringLength,
            "%s:%llu-%llu",
            region->chromosome,
            (unsigned long long)(region->start + 1U),
            (unsigned long long)region->end);
        fragmentsIterator = tbx_itr_querys(sourceHandle->tbxHandle, regionString);
        free(regionString);
        if (fragmentsIterator == NULL)
        {
            return ccounts_makeResult(-1, "failed to open fragments region iterator");
        }

        deltaBuffer = (float *)calloc(countBufferLength + 1U, sizeof(float));
        if (deltaBuffer == NULL)
        {
            hts_itr_destroy(fragmentsIterator);
            return ccounts_makeResult(-1, "failed to allocate fragments delta buffer");
        }

        while ((iteratorCode = tbx_itr_next(sourceHandle->fragmentsHandle, sourceHandle->tbxHandle, fragmentsIterator, &lineBuffer)) >= 0)
        {
            cursor = lineBuffer.s;
            fieldIndex = 0;
            fragStart = 0;
            fragEnd = 0;
            fragCount = 1;
            barcodeStart = NULL;
            barcodeLength = 0U;

            /* FFR: write inline helper for parsing fields */
            while (*cursor != '\0' && *cursor != '\n' && *cursor != '\r')
            {
                fieldStart = cursor;
                while (*cursor != '\0' && *cursor != '\t' && *cursor != '\n' && *cursor != '\r')
                {
                    ++cursor;
                }
                fieldLength = (size_t)(cursor - fieldStart);
                if (fieldIndex == 1 && !ccounts_parseInt64Field(fieldStart, fieldLength, &fragStart))
                {
                    fieldLength = 0U;
                    break;
                }
                if (fieldIndex == 2 && !ccounts_parseInt64Field(fieldStart, fieldLength, &fragEnd))
                {
                    fieldLength = 0U;
                    break;
                }
                if (fieldIndex == 3)
                {
                    barcodeStart = fieldStart;
                    barcodeLength = fieldLength;
                }
                if (fieldIndex == 4 && !ccounts_parseInt64Field(fieldStart, fieldLength, &fragCount))
                {
                    fragCount = 1;
                }
                if (*cursor == '\t')
                {
                    ++cursor;
                }
                ++fieldIndex;
            }

            if (fieldLength == 0U || fragEnd <= fragStart)
            {
                continue;
            }
            if (barcodeLength > 0U &&
                !ccounts_barcodeAllowed(
                    sourceHandle->barcodeAllowList,
                    sourceHandle->barcodeAllowCount,
                    barcodeStart,
                    barcodeLength))
            {
                continue;
            }

            incrementValue = (float)(fragCount > 0 ? fragCount : 1);
            if ((ccounts_countMode)countOptions->countMode == ccounts_countModeCenter ||
                countOptions->oneReadPerBin)
            {
                // oneReadPerBin --> fragments marked by their center
                midPoint = (fragStart + fragEnd) / 2;
                if (midPoint < start64 || midPoint >= end64)
                {
                    continue;
                }
                intervalIndex = (size_t)((midPoint - start64) / step64);
                if (intervalIndex < countBufferLength)
                {
                    countBuffer[intervalIndex] += incrementValue;
                }
                continue;
            }

            if ((ccounts_countMode)countOptions->countMode == ccounts_countModeCutSite ||
                (ccounts_countMode)countOptions->countMode == ccounts_countModeFivePrime)
            {
                // fragments are assumed to already represent insertion endpoints
                cutPosition = fragStart;
                if (cutPosition >= start64 && cutPosition < end64)
                {
                    intervalIndex = (size_t)((cutPosition - start64) / step64);
                    if (intervalIndex < countBufferLength)
                    {
                        countBuffer[intervalIndex] += incrementValue;
                    }
                }
                cutPosition = fragEnd - 1;
                if (cutPosition >= start64 && cutPosition < end64)
                {
                    intervalIndex = (size_t)((cutPosition - start64) / step64);
                    if (intervalIndex < countBufferLength)
                    {
                        countBuffer[intervalIndex] += incrementValue;
                    }
                }
                continue;
            }

            adjStart = fragStart;
            adjEnd = fragEnd;
            if (adjEnd <= start64 || adjStart >= end64)
            {
                continue;
            }
            if (adjStart < start64)
            {
                adjStart = start64;
            }
            if (adjEnd > end64)
            {
                adjEnd = end64;
            }
            index0 = (size_t)((adjStart - start64) / step64);
            index1 = (size_t)(((adjEnd - 1) - start64) / step64);
            if (index0 >= countBufferLength)
            {
                continue;
            }
            if (index1 >= countBufferLength)
            {
                index1 = countBufferLength - 1U;
            }
            if (index0 > index1)
            {
                continue;
            }
            /* deltaBuffer stores increments at the start of intervals and decrements at the end
             * st we can single-pass after processing all fragments */
            deltaBuffer[index0] += incrementValue;
            deltaBuffer[index1 + 1U] -= incrementValue;
        }

        deltaValue = 0.0f;
        for (intervalIndex = 0; intervalIndex < countBufferLength; ++intervalIndex)
        {
            deltaValue += deltaBuffer[intervalIndex];
            countBuffer[intervalIndex] += deltaValue;
        }

        if (lineBuffer.s != NULL)
        {
            free(lineBuffer.s);
        }
        hts_itr_destroy(fragmentsIterator);
        free(deltaBuffer);
        return ccounts_makeOk();
    }
    if (sourceHandle->fileHandle == NULL || sourceHandle->header == NULL)
    {
        return ccounts_makeResult(-1, "source handle is not initialized");
    }
    if (sourceHandle->indexHandle == NULL)
    {
        return ccounts_makeResult(-1, "alignment index is required for region counting");
    }
    if (countOptions->inferFragmentLength > 0 && countOptions->extendBP <= 0)
    {
        return ccounts_makeResult(-1, "native fragment length inference is not wired yet");
    }

    if (countOptions->threadCount > 1)
    {
        hts_set_threads((htsFile *)sourceHandle->fileHandle, countOptions->threadCount);
    }

    tid = sam_hdr_name2tid(sourceHandle->header, region->chromosome);
    if (tid < 0)
    {
        return ccounts_makeResult(-1, "chromosome not found in alignment header");
    }

    deltaBuffer = (float *)calloc(countBufferLength + 1U, sizeof(float));
    if (deltaBuffer == NULL)
    {
        return ccounts_makeResult(-1, "failed to allocate delta buffer");
    }

    record = bam_init1();
    if (record == NULL)
    {
        free(deltaBuffer);
        return ccounts_makeResult(-1, "failed to allocate bam record");
    }

    iteratorHandle = sam_itr_queryi(
        sourceHandle->indexHandle,
        tid,
        (hts_pos_t)region->start,
        (hts_pos_t)region->end);
    if (iteratorHandle == NULL)
    {
        bam_destroy1(record);
        free(deltaBuffer);
        return ccounts_makeResult(-1, "failed to open region iterator");
    }

    start64 = (int64_t)region->start;
    end64 = (int64_t)region->end;
    step64 = (int64_t)region->intervalSizeBP;
    // single-end paths fall back to readLength when minTemplateLength is unset
    minTemplateLength = countOptions->minTemplateLength >= 0
                            ? countOptions->minTemplateLength
                            : countOptions->readLength;

    while (sam_itr_next((htsFile *)sourceHandle->fileHandle, iteratorHandle, record) >= 0)
    {
        if (countOptions->flagInclude > 0 &&
            (record->core.flag & countOptions->flagInclude) != countOptions->flagInclude)
        {
            continue;
        }
        if ((record->core.flag & countOptions->flagExclude) != 0)
        {
            continue;
        }
        if (record->core.qual < countOptions->minMappingQuality)
        {
            continue;
        }

        readStart = (int64_t)record->core.pos;
        readEnd = (int64_t)bam_endpos(record);

        if (countOptions->pairedEndMode > 0)
        {
            if ((record->core.flag & BAM_FPROPER_PAIR) == 0)
            {
                continue;
            }
            if ((record->core.flag & BAM_FREAD2) != 0)
            {
                continue;
            }
            if ((record->core.flag & BAM_FMUNMAP) != 0 || record->core.mtid != record->core.tid)
            {
                continue;
            }

            templateLength = (int64_t)record->core.isize;
            absoluteTemplateLength = templateLength >= 0 ? templateLength : -templateLength;
            if (absoluteTemplateLength == 0 || absoluteTemplateLength < minTemplateLength)
            {
                continue;
            }
            if (countOptions->maxInsertSize > 0 && absoluteTemplateLength > countOptions->maxInsertSize)
            {
                continue;
            }

            // adjStart and adjEnd track the inferred fragment after shifting
            if (templateLength >= 0)
            {
                adjStart = readStart;
                adjEnd = readStart + absoluteTemplateLength;
            }
            else
            {
                adjEnd = readEnd;
                adjStart = adjEnd - absoluteTemplateLength;
            }

            if ((record->core.flag & BAM_FREVERSE) == 0)
            {
                adjStart += countOptions->shiftForwardStrand53;
                adjEnd += countOptions->shiftForwardStrand53;
            }
            else
            {
                adjStart -= countOptions->shiftReverseStrand53;
                adjEnd -= countOptions->shiftReverseStrand53;
            }
        }
        else
        {
            if ((record->core.flag & BAM_FREVERSE) == 0)
            {
                fivePrime = readStart + countOptions->shiftForwardStrand53;
                if (countOptions->extendBP > 0)
                {
                    adjStart = fivePrime;
                    adjEnd = fivePrime + countOptions->extendBP;
                }
                else
                {
                    adjStart = readStart + countOptions->shiftForwardStrand53;
                    adjEnd = readEnd + countOptions->shiftForwardStrand53;
                }
            }
            else
            {
                fivePrime = (readEnd - 1) - countOptions->shiftReverseStrand53;
                if (countOptions->extendBP > 0)
                {
                    adjEnd = fivePrime + 1;
                    adjStart = adjEnd - countOptions->extendBP;
                }
                else
                {
                    adjStart = readStart - countOptions->shiftReverseStrand53;
                    adjEnd = readEnd - countOptions->shiftReverseStrand53;
                }
            }
        }

        if (adjEnd <= start64 || adjStart >= end64)
        {
            continue;
        }
        if (adjStart < start64)
        {
            adjStart = start64;
        }
        if (adjEnd > end64)
        {
            adjEnd = end64;
        }

        if (countOptions->oneReadPerBin)
        {
            midPoint = (adjStart + adjEnd) / 2;
            intervalIndex = (size_t)((midPoint - start64) / step64);
            if (intervalIndex < countBufferLength)
            {
                countBuffer[intervalIndex] += 1.0f;
            }
            continue;
        }

        index0 = (size_t)((adjStart - start64) / step64);
        index1 = (size_t)(((adjEnd - 1) - start64) / step64);
        if (index0 >= countBufferLength)
        {
            continue;
        }
        if (index1 >= countBufferLength)
        {
            index1 = countBufferLength - 1U;
        }
        if (index0 > index1)
        {
            continue;
        }

        deltaBuffer[index0] += 1.0f;
        deltaBuffer[index1 + 1U] -= 1.0f;
    }

    deltaValue = 0.0f;
    for (intervalIndex = 0; intervalIndex < countBufferLength; ++intervalIndex)
    {
        deltaValue += deltaBuffer[intervalIndex];
        countBuffer[intervalIndex] += deltaValue;
    }

    hts_itr_destroy(iteratorHandle);
    bam_destroy1(record);
    free(deltaBuffer);
    return ccounts_makeOk();
}
