#ifndef DEALIMAGE_H
#define DEALIMAGE_H

#include "string"
#include "vector"
#include "slic.h"
#include "image.h"

using namespace std;

void ndwiOTSU(float *ndwi, int num, int &iFlag, float &ndwiFlag);

//��Ӱ�����������
void kmeans(float *ndwi, int num, int kNum, int*& ndwiClass);

void otsuClass(float* ndwiOrigin, int num, int *&ndwiClass);

//���������㷨
void areaGrow(int width, int height, int *&classFlag);

//���������㷨
void areaGrowHole(int width, int height, int *&classFlag);

void expendArea(int width, int height, int *&classFlag);

//NDWI����½�غ���ˮ�����Ӱ��
template <typename T>
void ndwiLandClass(string imgPath, T* pBuf, int width, int height, int *klabels, int numlables, int *&classFlag)
{
	classFlag = new int[width*height]; //Ӱ��½�ء�ǳˮ����ˮ�������ǩ
	memset(classFlag, 0, sizeof(int)*width*height);

	//�ó����ض�����Ӱ��ndwi����ƽ��
	float *ndwiSum = new float[numlables];
	int *ndwiNum = new int[numlables];
	memset(ndwiSum, 0.0, sizeof(float)*numlables);
	memset(ndwiNum, 0, sizeof(int)*numlables);
	//����ÿ�����ص�ndwiֵ
	float *ndwiValue = new float[width*height];
	//vector<int> kFlag; //������¼�쳣�����صı��
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int temp = i*width + j;
			float G = (float)pBuf[4 * temp + 1];
			float RNear = (float)pBuf[4 * temp + 3];
			if (G + RNear != 0)
			{
				ndwiValue[temp] = (G - RNear) / (G + RNear);
				ndwiSum[klabels[temp]] += ndwiValue[temp];
				ndwiNum[klabels[temp]]++;
			}

			//if (G + RNear == 0)
			//{
			//	if (kFlag.size() == 0)
			//	{
			//		kFlag.push_back(klabels[temp]);
			//	}
			//	else
			//	{
			//		int flag = 0; //��ʾû�иó����ر��
			//		for (int k = 0; k < kFlag.size(); k++)
			//		{
			//			if (kFlag[k] == klabels[temp])
			//			{
			//				flag = 1; //��ʾ�Ѿ����ڸó����ر��
			//			}
			//		}
			//		if (flag == 0)
			//		{
			//			kFlag.push_back(klabels[temp]);
			//		}
			//	}	
			//}
		} //end for j
	}// end for i

	//���쳣������ֵ���е���
	/*for (int i = 0; i < kFlag.size(); i++)
	{
		float w1 = 0.0;
		float w2 = 0.0;
		if (kFlag[i] - 1 >= 0 && kFlag[i] - 1 < numlables)
		{
			w1 = ndwiSum[kFlag[i] - 1];
		}
		if (kFlag[i] + 1 >= 0 && kFlag[i] + 1 < numlables)
		{
			w2 = ndwiSum[kFlag[i] - 1];
		}
		ndwiSum[kFlag[i]] = (w1 + w2) / 2;
	}*/
	
	string maskTxtPath = imgPath + "_SlicNdwi.txt";
	FILE *fp;
	fp = fopen(maskTxtPath.c_str(), "w");
	for (int k = 0; k < numlables; k++)
	{
		if (ndwiNum[k] != 0)
		{
			ndwiSum[k] /= ndwiNum[k];
		}	
		fprintf(fp, "%d\t%f\n",k, ndwiSum[k]);
	}
	fclose(fp);

	//�����ɵĳ����ط�Ϊ����(½�ء�ǳˮ����ˮ)
	int * ndwiClass = new int[numlables];
	//kmeans(ndwiSum, numlables, 3, ndwiClass); //��һ�ֻ���kmeans����

	//�ڶ��ֻ����趨����ֵ����
	/*for (int i = 0; i < numlables; i++)
	{
		if (ndwiSum[i]>=3.4)
		{
			ndwiClass[i] = 2;
		}
		else if (ndwiSum[i]>-1 && ndwiSum[i]<3.4)
		{
			ndwiClass[i] = 1;
		}
		else 
		{
			ndwiClass[i] = 0;
		}
	}*/

	//�����ֻ�������Ӧ��ֵ����
	otsuClass(ndwiSum, numlables, ndwiClass);
	
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int temp = i*width + j;
			classFlag[temp] = ndwiClass[klabels[temp]];
		}
	}
	
	delete[] ndwiSum; ndwiSum = NULL;
	delete[] ndwiNum; ndwiNum = NULL;	
	delete[] ndwiValue; ndwiValue = NULL;
	delete[] ndwiClass; ndwiClass = NULL;
	//vector<int>().swap(kFlag);
}

void performSeaLandSplit(string imgPath);

void performSeaLandSplitRecycle(string imgPath);


//���ڵ������ص�NDWI����½�غ���ˮ�����Ӱ��
template <typename T>
void ndwi(string imgPath, T* pBuf, int width, int height, int *&classFlag)
{
	classFlag = new int[width*height]; //Ӱ��½�ء�ǳˮ����ˮ�������ǩ
	memset(classFlag, 0, sizeof(int)*width*height);

	string maskTxtPath = imgPath + "_ndvi.txt";
	FILE *fp;
	fp = fopen(maskTxtPath.c_str(), "w");
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int temp = i*width + j;
			float G = (float)pBuf[4 * temp + 1];
			float RNear = (float)pBuf[4 * temp + 3];
			float ndwiValue = (G - RNear) / (G + RNear);
			if (ndwiValue >= 0.16)
			{
				classFlag[temp] = 2; //ǳˮ
			}
			if (ndwiValue >-0.2 && ndwiValue < 0.16)
			{
				classFlag[temp] = 1; //��ˮ
			}
			if (ndwiValue <= -0.2)
			{
				classFlag[temp] = 0; //½��
			}

			fprintf(fp, "%f\n", ndwiValue);
		}
	}
	fclose(fp);
}

#endif