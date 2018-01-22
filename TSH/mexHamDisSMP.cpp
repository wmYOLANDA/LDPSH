#define DLL_EXPORT_SYM

#include <mex.h>
#include <math.h>

#define IN_TrainData prhs[0]
#define IN_TrainBData prhs[1]
#define IN_PairsNum prhs[2]
#define IN_Index prhs[3]
#define IN_Flag prhs[4]
#define OUT_HammingDistance plhs[0]
#define OUT_EuclideanDistance plhs[1]

#define PI 3.1415926

typedef unsigned char byte;

void ComputeHammingDistance(byte* pBData, int* pIndex, double* hammingDis, int nNumData, int nBytes, int nPairsNum)
{
    for (int i = 0; i < nPairsNum; i++)
    {
        hammingDis[i] = 0;
        for (int k = 0; k < nBytes; k++)
        {
            byte tmp = (pBData[pIndex[i * 2] * nBytes + k] ^ pBData[pIndex[i * 2 + 1] * nBytes + k]);
            while (tmp)
            {
                hammingDis[i]++;
                tmp &= tmp - 1;
            }
        }
    }
}

void ComputeEuclideanDistance(double* pData, int* pIndex, double* euclideanDis, int nNumData, int nDim, int nPairsNum)
{
    for (int i = 0; i < nPairsNum; i++)
    {
        euclideanDis[i] = 0;
        for (int k = 0; k < nDim; k++)
            euclideanDis[i] += (pData[pIndex[i * 2] * nDim + k] - pData[pIndex[i * 2 + 1] * nDim + k]) *
                    (pData[pIndex[i * 2] * nDim + k] - pData[pIndex[i * 2 + 1] * nDim + k]);
        euclideanDis[i] = sqrt(euclideanDis[i]);
    }
}

void ComputeNormalizedAngleDistance(double* pData, int* pIndex, double* euclideanDis, int nNumData, int nDim, int nPairsNum)
{
    for (int i = 0; i < nPairsNum; i++)
    {
        double inProd = 0;
		double uNorm = 0;
		double vNorm = 0;
        for (int k = 0; k < nDim; k++) 
        {
			inProd += pData[pIndex[i * 2] * nDim + k] * pData[pIndex[i * 2 + 1] * nDim + k];
			uNorm += pData[pIndex[i * 2] * nDim + k] * pData[pIndex[i * 2] * nDim + k];
			vNorm += pData[pIndex[i * 2 + 1] * nDim + k] * pData[pIndex[i * 2 + 1] * nDim + k];
		}
		euclideanDis[i] = acos(inProd / sqrt(uNorm * vNorm)) / PI;
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    byte* pBData = (byte*)mxGetPr(IN_TrainBData);
    const int nNumData = (int)mxGetN(IN_TrainBData);
    const int nBytes = (int)mxGetM(IN_TrainBData);
    double* pData = (double*)mxGetPr(IN_TrainData);
    const int nDim = (int)mxGetM(IN_TrainData);
    const int nPairsNum = (int)mxGetScalar(IN_PairsNum);
    int* pIndex = (int*)mxGetPr(IN_Index);
    const int flag = (int)mxGetScalar(IN_Flag);
    
    printf("nNumData: %d,\tnBits: %d, \tnDim: %d\n", nNumData, nBytes, nDim);
    
    OUT_HammingDistance = mxCreateDoubleMatrix(nPairsNum, 1, mxREAL);
    double* hammingDis = (double*)mxGetPr(OUT_HammingDistance);
    OUT_EuclideanDistance = mxCreateDoubleMatrix(nPairsNum, 1, mxREAL);
    double* euclideanDis = (double*)mxGetPr(OUT_EuclideanDistance);
    
    // ¼ÆËãHamming¾àÀë
    ComputeHammingDistance(pBData, pIndex, hammingDis, nNumData, nBytes, nPairsNum);
   
    // ¼ÆËãEuclidean¾àÀë
    if (flag == 1)
        ComputeEuclideanDistance(pData, pIndex, euclideanDis, nNumData, nDim, nPairsNum);
    if (flag == 2)
        ComputeNormalizedAngleDistance(pData, pIndex, euclideanDis, nNumData, nDim, nPairsNum);
    return;
}