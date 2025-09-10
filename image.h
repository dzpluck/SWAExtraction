#ifndef IMAGE_H
#define IMAGE_H

#include<gdal.h>
#include<ogrsf_frmts.h>
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "iostream"
#include "string"
#include "vector"

using namespace std;


//GDAL�ֿ��ȡӰ��
template <typename T>
int readImage(string filePath, int& width, int& height, int &nBands, GDALDataType& gdt, T*& pBuf, double*& trans, OGRSpatialReference &projSRS)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	GDALDataset *pdataset = NULL;
	pdataset = (GDALDataset *)GDALOpen(filePath.c_str(), GA_ReadOnly);
	if (pdataset == NULL)
	{
		return 0;
	}

	gdt = pdataset->GetRasterBand(1)->GetRasterDataType();//��������
	width = pdataset->GetRasterXSize();//���ο��
	height = pdataset->GetRasterYSize();//���θ߶�
	nBands = pdataset->GetRasterCount();//������
	if (nBands < 4)
	{
		cout << "Ӱ�񲨶���С��4������" << endl;
		return 0;
	}

	int i, j;
	pBuf = new T[height*width * 4];
	memset(pBuf, 0, sizeof(T)*height*width * 4);
	//pl = new double[height*width * 2]; //�����洢ͶӰ�Ժ������

	//��ȡ����任ϵ��
	trans = new double[6];
	CPLErr error = pdataset->GetGeoTransform(trans); //��������任����
	const char* PrjRef = pdataset->GetProjectionRef();//��ȡӰ��ͶӰ��Ϣ;
	projSRS.importFromWkt((char**)&PrjRef);  //�����ж�Ӱ���ǵ������껹��ͶӰ����

	int Block = 2000;
	int xNum = (width - 1) / Block + 1; //�зֳ�xNum����
	int yNum = (height - 1) / Block + 1; //�зֳ�yNum����
	for (i = 0; i<yNum; i++)
	{
		for (j = 0; j<xNum; j++)
		{
			int BlockPosX;	//��ǰ�������ԭͼ����ʼx����
			int BlockPosY;	//��ǰ�������ԭͼ����ʼy����
			BlockPosX = Block * j;
			BlockPosY = Block * i;

			int curBlockXSize = Block;		//��ǰ��� width
			int curBlockYSize = Block;		//��ǰ��� height

			if (j == xNum - 1)   //��β
			{
				curBlockXSize = width - Block * j;
			}
			if (i == yNum - 1) //��β
			{
				curBlockYSize = height - Block *  i;
			}

			T *pBuf1 = new T[curBlockXSize*curBlockYSize];
			T *pBuf2 = new T[curBlockXSize*curBlockYSize];
			T *pBuf3 = new T[curBlockXSize*curBlockYSize];
			T *pBuf4 = new T[curBlockXSize*curBlockYSize];

			pdataset->GetRasterBand(1)->RasterIO(GF_Read, BlockPosX, BlockPosY, curBlockXSize,
				curBlockYSize, pBuf1, curBlockXSize, curBlockYSize, gdt, 0, 0);
			pdataset->GetRasterBand(2)->RasterIO(GF_Read, BlockPosX, BlockPosY, curBlockXSize,
				curBlockYSize, pBuf2, curBlockXSize, curBlockYSize, gdt, 0, 0);
			pdataset->GetRasterBand(3)->RasterIO(GF_Read, BlockPosX, BlockPosY, curBlockXSize,
				curBlockYSize, pBuf3, curBlockXSize, curBlockYSize, gdt, 0, 0);
			pdataset->GetRasterBand(4)->RasterIO(GF_Read, BlockPosX, BlockPosY, curBlockXSize,
				curBlockYSize, pBuf4, curBlockXSize, curBlockYSize, gdt, 0, 0);

			for (int k = 0; k<curBlockYSize; k++)
			{
				for (int w = 0; w<curBlockXSize; w++)
				{
					int temp = k*curBlockXSize + w;
					int x0 = BlockPosX;
					int y0 = BlockPosY;
					int x1 = x0 + w;
					int y1 = y0 + k;
					int temp1 = y1*width + x1;
					pBuf[4 * temp1] = pBuf1[temp];  //�洢�����̡��졢�������ĸ�����
					pBuf[4 * temp1 + 1] = pBuf2[temp];
					pBuf[4 * temp1 + 2] = pBuf3[temp];
					pBuf[4 * temp1 + 3] = pBuf4[temp];
         
				}
			}

			delete[]pBuf1; pBuf1 = NULL;
			delete[]pBuf2; pBuf2 = NULL;
			delete[]pBuf3; pBuf3 = NULL;
			delete[]pBuf4; pBuf4 = NULL;
		}//end for
	}//end for

	if (pdataset != NULL)
	{
		GDALClose(pdataset);
		pdataset = NULL;
	}

	return 1;
}


//GDAL�ֿ鱣�溣½����Ӱ����Ĥ���
void saveImageMask(string filepath, int width, int height, int nchannel, GDALDataType gdt, double * trans, OGRSpatialReference oSRS, unsigned char *wd);

//GDAL�ֿ鱣�泬�������ɽ��
void saveImage2(int width, int height, int nchannel, GDALDataType gdt, unsigned char *wd);

#endif