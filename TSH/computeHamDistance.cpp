#define DLL_EXPORT_SYM

#include <mex.h>
#include <math.h>

#define IN_TrainData prhs[0]
#define IN_TrainBData prhs[1]
#define IN_NeighborNum prhs[2]
#define IN_Index prhs[3]
#define OUT_HammingDistance plhs[0]
//#define OUT_EuclideanDistance plhs[1]

typedef unsigned char byte;

void ComputeHammingDistance(byte* pBData, int* pIndex, double* hammingDis, int nNumData, int nBytes, int nNeighborNum)
{
    for (int i = 0; i < nNumData; i++)
    {
        for (int j = 0; j < nNeighborNum; j++)
        {
            hammingDis[i * nNeighborNum + j] = 0;
            for (int k = 0; k < nBytes; k++)
            {
                byte tmp = (pBData[i * nBytes + k] ^ pBData[pIndex[i * nNeighborNum + j] * nBytes + k]);
                while (tmp)
                {
                    hammingDis[i * nNeighborNum + j]++;
                    tmp &= tmp - 1;
                }
            }
        }
    }
}

void ComputeEuclideanDistance(double* pData, int* pIndex, double* euclideanDis, int nNumData, int nDim, int nNeighborNum)
{
    for (int i = 0; i < nNumData; i++)
    {
        for (int j = 0; j < nNeighborNum; j++)
        {
            euclideanDis[i * nNeighborNum + j] = 0;
            for (int k = 0; k < nDim; k++)
                euclideanDis[i * nNeighborNum + j] += (pData[i * nDim + k] - pData[pIndex[i * nNeighborNum + j] * nDim + k]) *
                        (pData[i * nDim + k] - pData[pIndex[i * nNeighborNum + j] * nDim + k]);
            euclideanDis[i * nNeighborNum + j] = sqrt(euclideanDis[i * nNeighborNum + j]);
        }
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    byte* pBData = (byte*)mxGetPr(IN_TrainBData);
    const int nNumData = (int)mxGetN(IN_TrainBData);
    const int nBytes = (int)mxGetM(IN_TrainBData);
    double* pData = (double*)mxGetPr(IN_TrainData);
    const int nDim = (int)mxGetM(IN_TrainData);
    const int nNeighborNum = (int)mxGetScalar(IN_NeighborNum);
    int* pIndex = (int*)mxGetPr(IN_Index);
    
    printf("nNumData: %d,\tnBits: %d, \tnDim: %d", nNumData, nBytes, nDim);
    
    OUT_HammingDistance = mxCreateDoubleMatrix(nNumData, nNeighborNum, mxREAL);
    double* hammingDis = (double*)mxGetPr(OUT_HammingDistance);
 //   OUT_EuclideanDistance = mxCreateDoubleMatrix(nNumData, nNeighborNum, mxREAL);
 //   double* euclideanDis = (double*)mxGetPr(OUT_EuclideanDistance);
    
    // ¼ÆËãHamming¾àÀë
    ComputeHammingDistance(pBData, pIndex, hammingDis, nNumData, nBytes, nNeighborNum);
            
    // ¼ÆËãEuclidean¾àÀë
//   ComputeEuclideanDistance(pData, pIndex, euclideanDis, nNumData, nDim, nNeighborNum);
    
    return;
}