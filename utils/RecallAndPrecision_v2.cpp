#define DLL_EXPORT_SYM

#include <mex.h>
#include <intrin.h> // __popcnt()头文件
#include <immintrin.h>
#include <algorithm>

#define IN_BincodeXtrain prhs[0]
#define IN_BincodeXtest prhs[1]
#define IN_IntSel prhs[2]
#define IN_Gt prhs[3]
#define IN_K prhs[4]
#define OUT_Recall plhs[0]
#define OUT_Precision plhs[1]
#define OUT_Recallk plhs[2]
#define OUT_Precisionk plhs[3]
#define OUT_MAP plhs[4]

typedef unsigned char byte;

struct HammingDistIdx
{
    int nDist;
    int nIdx;
};
typedef struct HammingDistIdx HammingDistIdx;

bool compare_dist(HammingDistIdx item_a, HammingDistIdx item_b){ return (item_a.nDist < item_b.nDist); }

void hamming_search_8bits(int* pNumSuccess, double* pNumSuccessk, double* pPrecisionk, double* pAp, int* pGt, const byte* pDataBin, const byte* pQueryBin,
        const int nNumData, const int nNumQueries, const int nSel, const int nNumNeighbor)
{
    typedef unsigned char EntryType;
    const int nBytes = 1;
    EntryType* pData;
    EntryType* pQuery;
    
    HammingDistIdx* pHammingDist = (HammingDistIdx*)_aligned_malloc(nNumData * sizeof(HammingDistIdx), 16);
    pQuery = (EntryType*)pQueryBin;
    for (int idx_query = 0; idx_query < nNumQueries; idx_query++)
    {
        pData = (EntryType*)pDataBin;
        
        for(int idx_data = 0; idx_data < nNumData; idx_data++)
        {
            int query = *pQuery;
            int data = *pData;
            pHammingDist[idx_data].nDist = (int)__popcnt(query ^ data);
            pData++;
            pHammingDist[idx_data].nIdx = idx_data;
        }
        
        std::partial_sort(pHammingDist, pHammingDist + nSel, pHammingDist + nNumData, compare_dist);
        
        double ap = 0;
        int trueNum = 0;
        for (int i = 0; i < nSel; i++) {
            int j = 0;
            for (j = 0; j < nNumNeighbor; j++) {
                if (pHammingDist[i].nIdx == (*(pGt + j))) {
                    break;
                }
            }
            if (j < nNumNeighbor) {
                pNumSuccess[i] ++;
                trueNum++;
                pNumSuccessk[(int)(pHammingDist[i].nDist)] += 1;
                ap += (double)trueNum / (i + 1);
            }
            pPrecisionk[(int)(pHammingDist[i].nDist)] += 1;
        }

        if (trueNum > 0)
            pAp[idx_query] = ap / trueNum;
        pGt += nNumNeighbor;
        pQuery = pQuery + 1;
    }
    
    _aligned_free(pHammingDist);
    // free(pHammingDist);
}

void hamming_search_16bits(int* pNumSuccess, double* pNumSuccessk, double* pPrecisionk, double* pAp, int* pGt, const byte* pDataBin, const byte* pQueryBin,
        const int nNumData, const int nNumQueries, const int nSel, const int nNumNeighbor)
{
    typedef unsigned char EntryType;
    const int nBytes = 2;
    EntryType* pData;
    EntryType* pQuery;
    
    HammingDistIdx* pHammingDist = (HammingDistIdx*)_aligned_malloc(nNumData * sizeof(HammingDistIdx), 16);
    pQuery = (EntryType*)pQueryBin;
    for (int idx_query = 0; idx_query < nNumQueries; idx_query++)
    {
        pData = (EntryType*)pDataBin;
        int query = *pQuery;
        query = (query << 8) + *(pQuery + 1);
        for(int idx_data = 0; idx_data < nNumData; idx_data++)
        {
            int data = *pData;
            data = (data << 8) + *(pData + 1);
            pHammingDist[idx_data].nDist = (int)__popcnt(query ^ data);
            pData = pData + 2;
            pHammingDist[idx_data].nIdx = idx_data;
        }
        
        std::partial_sort(pHammingDist, pHammingDist + nSel, pHammingDist + nNumData, compare_dist);
        
        double ap = 0;
        int trueNum = 0;
        for (int i = 0; i < nSel; i++) {
            int j = 0;
            for (j = 0; j < nNumNeighbor; j++) {
                if (pHammingDist[i].nIdx == (*(pGt + j))) {
                    break;
                }
            }
            if (j < nNumNeighbor) {
                pNumSuccess[i] ++;
                trueNum++;
                pNumSuccessk[(int)(pHammingDist[i].nDist)] += 1;
                ap += (double)trueNum / (i + 1);
            }
            pPrecisionk[(int)(pHammingDist[i].nDist)] += 1;
        }

        if (trueNum > 0)
            pAp[idx_query] = ap / trueNum;
        pGt += nNumNeighbor;
        pQuery = pQuery + 2;
    }
    
    _aligned_free(pHammingDist);
    // free(pHammingDist);
}

void hamming_search_32bits(int* pNumSuccess, double* pNumSuccessk, double* pPrecisionk, double* pAp, int* pGt, const byte* pDataBin, const byte* pQueryBin,
        const int nNumData, const int nNumQueries, const int nSel, const int nNumNeighbor)
{
    typedef unsigned int EntryType;
    const int nBytes = 4;
    EntryType* pData;
    EntryType* pQuery;
    
    HammingDistIdx* pHammingDist = (HammingDistIdx*)_aligned_malloc(nNumData * sizeof(HammingDistIdx), 16);
    // HammingDistIdx* pHammingDist = (HammingDistIdx*)malloc(nNumData * sizeof(HammingDistIdx));
    pQuery = (EntryType*)pQueryBin;
    int* trueNumk = new int[nBytes * 8 + 1];
    for (int idx_query = 0; idx_query < nNumQueries; idx_query++)
    {
        pData = (EntryType*)pDataBin;
        
        for(int idx_data = 0; idx_data < nNumData; idx_data++)
        {
            pHammingDist[idx_data].nDist = (int)__popcnt(*pData ^ *pQuery);
            pData++;
            pHammingDist[idx_data].nIdx = idx_data;
        }
        
        std::partial_sort(pHammingDist, pHammingDist + nSel, pHammingDist + nNumData, compare_dist);
        
        double ap = 0;
        int trueNum = 0;
        for (int i = 0; i < nSel; i++) {
            int j = 0;
            for (j = 0; j < nNumNeighbor; j++) {
                if (pHammingDist[i].nIdx == (*(pGt + j))) {
                    break;
                }
            }
            if (j < nNumNeighbor) {
                pNumSuccess[i] ++;
                trueNum++;
                pNumSuccessk[(int)(pHammingDist[i].nDist)] += 1;
                ap += (double)trueNum / (i + 1);
            }
            pPrecisionk[(int)(pHammingDist[i].nDist)] += 1;
        }

        if (trueNum > 0)
            pAp[idx_query] = ap / trueNum;
        pQuery++;
        pGt += nNumNeighbor;
    }
    
    _aligned_free(pHammingDist);
    // free(pHammingDist);
}

void hamming_search_64bits(int* pNumSuccess, double* pNumSuccessk, double* pPrecisionk, double* pAp, int* pGt, const byte* pDataBin, const byte* pQueryBin,
        const int nNumData, const int nNumQueries, const int nSel, const int nNumNeighbor)
{
    typedef unsigned __int64 EntryType;
    const int nBytes = 8;
    EntryType* pData;
    EntryType* pQuery;
    
    HammingDistIdx* pHammingDist = (HammingDistIdx*)_aligned_malloc(nNumData * sizeof(HammingDistIdx), 16);
    pQuery = (EntryType*)pQueryBin;
    for (int idx_query = 0; idx_query < nNumQueries; idx_query++)
    {
        pData = (EntryType*)pDataBin;
        
        for (int idx_data = 0; idx_data < nNumData; idx_data++)
        {
            pHammingDist[idx_data].nDist = (int) __popcnt64(*pQuery ^ *pData++);
            pHammingDist[idx_data].nIdx = idx_data;
        }
        
        std::partial_sort(pHammingDist, pHammingDist + nSel, pHammingDist + nNumData, compare_dist);
        
        double ap = 0;
        int trueNum = 0;
        for (int i = 0; i < nSel; i++) {
            int j = 0;
            for (j = 0; j < nNumNeighbor; j++) {
                if (pHammingDist[i].nIdx == (*(pGt + j))) {
                    break;
                }
            }
            if (j < nNumNeighbor) {
                pNumSuccess[i] ++;
                trueNum++;
                pNumSuccessk[(int)(pHammingDist[i].nDist)] += 1;
                ap += (double)trueNum / (i + 1);
            }
            pPrecisionk[(int)(pHammingDist[i].nDist)] += 1;
        }

        if (trueNum > 0)
            pAp[idx_query] = ap / trueNum;
        pQuery++;
        pGt += nNumNeighbor;
    }
    
    _aligned_free(pHammingDist);
}

void hamming_search_128bits(int* pNumSuccess, double* pNumSuccessk, double* pPrecisionk, double* pAp, int* pGt, const byte* pDataBin, const byte* pQueryBin,
        const int nNumData, const int nNumQueries, const int nSel, const int nNumNeighbor)
{
    typedef __m128i EntryType; // 128bits
    const int nBytes = 16;
    EntryType* pData;
    EntryType* pQuery;
    
    HammingDistIdx* pHammingDist = (HammingDistIdx*)_aligned_malloc(nNumData * sizeof(HammingDistIdx), 16);
    pQuery = (EntryType*)pQueryBin;
    int* trueNumk = new int[nBytes * 8 + 1];
    for (int idx_query = 0; idx_query < nNumQueries; idx_query++)
    {
        pData = (EntryType*)pDataBin;
        unsigned __int64 *pQuerySec0 = (unsigned __int64*)pQuery;
        unsigned __int64 *pQuerySec1 = pQuerySec0 + 1;
        for (int idx_data = 0; idx_data < nNumData; idx_data++)
        {
            unsigned __int64 *pDataSec0 = (unsigned __int64*)pData;
            unsigned __int64 *pDataSec1 = pDataSec0 + 1;
            pHammingDist[idx_data].nDist = (int)__popcnt64(*pDataSec0 ^ *pQuerySec0) + (int)__popcnt64(*pDataSec1 ^ *pQuerySec1); //popcnt计算64比特中1的个数——Hamming重量
            pHammingDist[idx_data].nIdx = idx_data;
            pData++;
        }
        std::partial_sort(pHammingDist, pHammingDist + nSel, pHammingDist + nNumData, compare_dist);
        
        double ap = 0;
        int trueNum = 0;
        for (int i = 0; i < nSel; i++) {
            int j = 0;
            for (j = 0; j < nNumNeighbor; j++) {
                if (pHammingDist[i].nIdx == (*(pGt + j))) {
                    break;
                }
            }
            if (j < nNumNeighbor) {
                pNumSuccess[i] ++;
                trueNum++;
                pNumSuccessk[(int)(pHammingDist[i].nDist)] += 1;
                ap += (double)trueNum / (i + 1);
            }
            pPrecisionk[(int)(pHammingDist[i].nDist)] += 1;
        }

        if (trueNum > 0)
            pAp[idx_query] = ap / trueNum;
        pGt += nNumNeighbor;
        pQuery = pQuery + 1;
    }
    
    _aligned_free(pHammingDist);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    byte* pDataBin = (byte*)mxGetPr(IN_BincodeXtrain);
    const int nNumData = (int)mxGetN(IN_BincodeXtrain);
    const int nBytes = (int)mxGetM(IN_BincodeXtrain);
    byte* pQueryBin = (byte*)mxGetPr(IN_BincodeXtest);
    const int nNumQueries = (int)mxGetN(IN_BincodeXtest);
    const int nNumSel = (int)mxGetScalar(IN_IntSel);
    int* pGt = (int*)mxGetPr(IN_Gt);
    const int nNumNeighbor = (int)mxGetM(IN_Gt);
    const int nNumN = (int)mxGetN(IN_Gt);
    const int nK = (int)mxGetScalar(IN_K);
    
    printf("nNumData: %d\n", nNumData);
    printf("nNumQueries: %d\n", nNumQueries);
    printf("nBytes: %d\n", nBytes);
    printf("nBytes: %d\n", (int)mxGetM(IN_BincodeXtest));
    printf("nSelectivity: %d\n", nNumSel);
    printf("nNumN: %d\n", nNumN);
    printf("nK: %d\n", nK);
    
    // Hamming distance based ranking
    const int nNumCores = 100;
    const int nSubNumQueries = nNumQueries / nNumCores;
    int** pNumSuccess = new int*[nNumCores];
    for (int core = 0; core < nNumCores; core++)
        pNumSuccess[core] = new int[nNumSel];
    double** pNumSuccessk = new double*[nNumCores];
    for (int core = 0; core < nNumCores; core++)
        pNumSuccessk[core] = new double[nBytes * 8 + 1];
    double** pPrecisionk = new double*[nNumCores];
    for (int core = 0; core < nNumCores; core++)
        pPrecisionk[core] = new double[nBytes * 8 + 1];
    double* pAp = new double[nNumQueries];
    memset(pAp, 0, sizeof(double) * nNumQueries);
    
    printf("success!\n");
    
    if (nBytes == 1)
    {
#pragma omp parallel for
        for (int core = 0; core < nNumCores; core++)
        {
            memset(pNumSuccess[core], 0, sizeof(int) * nNumSel);
            memset(pNumSuccessk[core], 0, sizeof(double) * (nBytes * 8 + 1));
            memset(pPrecisionk[core], 0, sizeof(double) * (nBytes * 8 + 1));
            hamming_search_8bits(pNumSuccess[core], pNumSuccessk[core], pPrecisionk[core], pAp + core * nSubNumQueries, pGt + core * nSubNumQueries * nNumNeighbor, pDataBin,
                    pQueryBin + core * nSubNumQueries * nBytes, nNumData, nSubNumQueries, nNumSel, nNumNeighbor);
            printf("success! %d\n", core);
        }
    }
    else if (nBytes == 2)
    {
#pragma omp parallel for
        for (int core = 0; core < nNumCores; core++)
        {
            memset(pNumSuccess[core], 0, sizeof(int) * nNumSel);
            memset(pNumSuccessk[core], 0, sizeof(double) * (nBytes * 8 + 1));
            memset(pPrecisionk[core], 0, sizeof(double) * (nBytes * 8 + 1));
            hamming_search_16bits(pNumSuccess[core], pNumSuccessk[core], pPrecisionk[core], pAp + core * nSubNumQueries, pGt + core * nSubNumQueries * nNumNeighbor, pDataBin,
                    pQueryBin + core * nSubNumQueries * nBytes, nNumData, nSubNumQueries, nNumSel, nNumNeighbor);
        }
    }
    else if (nBytes == 4)
    {
#pragma omp parallel for
        for (int core = 0; core < nNumCores; core++)
        {
            memset(pNumSuccess[core], 0, sizeof(int) * nNumSel);
            memset(pNumSuccessk[core], 0, sizeof(double) * (nBytes * 8 + 1));
            memset(pPrecisionk[core], 0, sizeof(double) * (nBytes * 8 + 1));
            hamming_search_32bits(pNumSuccess[core], pNumSuccessk[core], pPrecisionk[core], pAp + core * nSubNumQueries, pGt + core * nSubNumQueries * nNumNeighbor, pDataBin,
                    pQueryBin + core * nSubNumQueries * nBytes, nNumData, nSubNumQueries, nNumSel, nNumNeighbor);
        }
    }
    else if (nBytes == 8)
    {
#pragma omp parallel for
        for (int core = 0; core < nNumCores; core++)
        {
            memset(pNumSuccess[core], 0, sizeof(int) * nNumSel);
            memset(pNumSuccessk[core], 0, sizeof(double) * (nBytes * 8 + 1));
            memset(pPrecisionk[core], 0, sizeof(double) * (nBytes * 8 + 1));
            hamming_search_64bits(pNumSuccess[core],pNumSuccessk[core],  pPrecisionk[core], pAp + core * nSubNumQueries, pGt + core * nSubNumQueries * nNumNeighbor, pDataBin,
                    pQueryBin + core * nSubNumQueries * nBytes, nNumData, nSubNumQueries, nNumSel, nNumNeighbor);
        }
    }
    else if (nBytes == 16)
    {
#pragma omp parallel for
        for (int core = 0; core < nNumCores; core++)
        {
            memset(pNumSuccess[core], 0, sizeof(int) * nNumSel);
            memset(pNumSuccessk[core], 0, sizeof(double) * (nBytes * 8 + 1));
            memset(pPrecisionk[core], 0, sizeof(double) * (nBytes * 8 + 1));
            hamming_search_128bits(pNumSuccess[core], pNumSuccessk[core], pPrecisionk[core], pAp + core * nSubNumQueries, pGt + core * nSubNumQueries * nNumNeighbor, pDataBin,
                    pQueryBin + core * nSubNumQueries * nBytes, nNumData, nSubNumQueries, nNumSel, nNumNeighbor);
        }
    }
    else
    {
        printf("False bit number.\n");
        return;
    }
    
    OUT_Recall = mxCreateDoubleMatrix(nNumSel, 1, mxREAL);
    double* recall = (double*)mxGetPr(OUT_Recall);
    OUT_Precision = mxCreateDoubleMatrix(nNumSel, 1, mxREAL);
    double* precision = (double*)mxGetPr(OUT_Precision);
    OUT_Recallk = mxCreateDoubleMatrix(nBytes * 8 + 1, 1, mxREAL);
    double* recallk = (double*) mxGetPr(OUT_Recallk);
    OUT_Precisionk = mxCreateDoubleMatrix(nBytes * 8 + 1, 1, mxREAL);
    double* precisionk = (double*) mxGetPr(OUT_Precisionk);
    OUT_MAP = mxCreateDoubleMatrix(1, 1, mxREAL);
    double* map = (double*) mxGetPr(OUT_MAP);
    
    map[0] = 0;
    for (int i = 0; i < nNumQueries; i++) {
        map[0] += pAp[i] / nNumQueries;
    }
    
    for (int j = 0; j < nNumSel; j++)
    {
        recall[j] = 0;
        precision[j] = 0;
        for (int i = 0; i < nNumCores; i++)
        {
            if (j > 0)
                pNumSuccess[i][j] += pNumSuccess[i][j - 1];
            recall[j] += (double)pNumSuccess[i][j] / double(nNumQueries * nK);
            precision[j] += (double)pNumSuccess[i][j] / double(nNumQueries * (j + 1));
        }
    }
    
    
    for (int j = 0; j < nBytes * 8 + 1; j++)
    {
        precisionk[j] = 0;
        recallk[j] = 0;
        double tmp1 = 0, tmp2 = 0;
        for (int i = 0; i < nNumCores; i++)
        {
            if (j > 0) {
                pNumSuccessk[i][j] += pNumSuccessk[i][j - 1];
                pPrecisionk[i][j] += pPrecisionk[i][j - 1];
            }
            
            tmp1 += pNumSuccessk[i][j];
            tmp2 += pPrecisionk[i][j];
        }
        recallk[j] += tmp1 / (nNumQueries * nK);
        precisionk[j] += tmp1 / tmp2;
    }
    
    for (int i = 0; i < nNumCores; i++) {
        delete[] pNumSuccess[i];
        delete[] pNumSuccessk[i];
        delete[] pPrecisionk[i];
    }
    delete[] pNumSuccess;
    delete[] pNumSuccessk;
    delete[] pPrecisionk;
    delete[] pAp;
    return;
}