#ifndef CCOUNTS_BACKEND_H
#define CCOUNTS_BACKEND_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

    typedef enum ccounts_sourceKind
    {
        ccounts_sourceKindBAM = 0,
        ccounts_sourceKindFragments = 1
    } ccounts_sourceKind;

    typedef enum ccounts_countMode
    {
        ccounts_countModeCoverage = 0,
        ccounts_countModeCutSite = 1,
        ccounts_countModeFivePrime = 2,
        ccounts_countModeCenter = 3
    } ccounts_countMode;

    /**
     * @brief input source-specific configuration
     *
     * describes one BAM or fragments input plus barcode filter state
     */
    typedef struct ccounts_sourceConfig
    {
        const char *path;
        ccounts_sourceKind sourceKind;
        const char *barcodeTag;
        const char *barcodeAllowListFile;
        const char *barcodeGroupMapFile;
    } ccounts_sourceConfig;

    /**
     * @brief evenly spaced genomic interval request
     *
     * defines the chromosome span and bin width for one counting call
     */
    typedef struct ccounts_region
    {
        const char *chromosome;
        uint32_t start;
        uint32_t end;
        uint32_t intervalSizeBP;
    } ccounts_region;

    /**
     * @brief counting and filtering options for one region query
     *
     * controls read filtering, fragment interpretation, shifting,
     * extension, and how signal is distributed across bins (e.g. coverage vs cut site)
     */
    typedef struct ccounts_countOptions
    {
        uint16_t threadCount;
        uint16_t flagInclude;
        uint16_t flagExclude;
        uint8_t countMode;
        uint8_t oneReadPerBin;
        int64_t shiftForwardStrand53;
        int64_t shiftReverseStrand53;
        int64_t readLength;
        int64_t extendBP;
        int64_t minMappingQuality;
        int64_t minTemplateLength;
        int64_t maxInsertSize;
        int64_t pairedEndMode;
        int64_t inferFragmentLength;
    } ccounts_countOptions;

    typedef struct ccounts_result
    {
        int32_t errorCode;
        const char *errorMessage;
    } ccounts_result;

    typedef struct ccounts_sourceHandle ccounts_sourceHandle;

    ccounts_result ccounts_checkAlignmentFile(
        const ccounts_sourceConfig *sourceConfig,
        int buildIndex,
        int threadCount,
        int *hasIndexOut);

    ccounts_result ccounts_isPairedEnd(
        const ccounts_sourceConfig *sourceConfig,
        int threadCount,
        int maxReads,
        int *isPairedEndOut);

    ccounts_result ccounts_getReadLength(
        const ccounts_sourceConfig *sourceConfig,
        int threadCount,
        int minReads,
        int maxIterations,
        int flagExclude,
        uint32_t *readLengthOut);

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
        uint32_t *fragmentLengthOut);

    ccounts_result ccounts_getChromRange(
        const ccounts_sourceConfig *sourceConfig,
        const char *chromosome,
        uint64_t chromLength,
        int threadCount,
        int flagExclude,
        uint64_t *startOut,
        uint64_t *endOut);

    ccounts_result ccounts_getMappedReadCount(
        const ccounts_sourceConfig *sourceConfig,
        int threadCount,
        const char *const *excludeChromosomes,
        int excludeChromosomeCount,
        uint8_t countMode,
        uint8_t oneReadPerBin,
        uint64_t *mappedReadCountOut,
        uint64_t *unmappedReadCountOut);

    ccounts_result ccounts_getCellCount(
        const ccounts_sourceConfig *sourceConfig,
        uint64_t *cellCountOut);

    ccounts_result ccounts_openSource(
        const ccounts_sourceConfig *sourceConfig,
        ccounts_sourceHandle **sourceHandleOut);

    void ccounts_closeSource(ccounts_sourceHandle *sourceHandle);

    ccounts_result ccounts_countRegion(
        ccounts_sourceHandle *sourceHandle,
        const ccounts_region *region,
        const ccounts_countOptions *countOptions,
        float *countBuffer,
        size_t countBufferLength);

#ifdef __cplusplus
}
#endif

#endif
