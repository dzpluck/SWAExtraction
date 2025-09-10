#include "stdafx.h"
#include "image.h"

//GDAL�ֿ鱣�溣½����Ӱ����Ĥ���
void saveImageMask(string filepath, int width, int height, int nchannel, GDALDataType gdt, double * trans, OGRSpatialReference oSRS, unsigned char *wd)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	GDALDriver *poDriver = NULL;
	GDALDataset *pImgOut = NULL;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTIFF");

	if (poDriver != NULL)
	{
		pImgOut = poDriver->Create(filepath.c_str(), width, height, nchannel, gdt, NULL);

		pImgOut->SetGeoTransform(trans);
		char *pszSRS_WKT = NULL;
		oSRS.exportToWkt(&pszSRS_WKT);
		pImgOut->SetProjection(pszSRS_WKT);
		CPLFree(pszSRS_WKT);//ʹ������ͷ�

		if (pImgOut != NULL)
		{
			int i, j;
			int Block = 2000;
			int xNum = (width - 1) / Block + 1; //�зֳ�xNum����
			int yNum = (height - 1) / Block + 1; //�зֳ�yNum����

			for (i = 0; i<yNum; i++)
			{
				for (j = 0; j<xNum; j++)
				{
					int BlockPosX = 0;	//��ǰ�������ԭͼ����ʼx����
					int BlockPosY = 0;	//��ǰ�������ԭͼ����ʼy����
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

					int widthx = curBlockXSize;
					int heighty = curBlockYSize;

					int m = widthx * heighty;
					unsigned char* ubuf1 = new unsigned char[m];
					unsigned char* ubuf2 = new unsigned char[m];
					unsigned char* ubuf3 = new unsigned char[m];

					for (int k = 0; k < heighty; k++)
					{
						for (int w = 0; w < widthx; w++)
						{
							int temp = k*widthx + w;
							int x0 = j*Block;//ÿ����ʼ�������
							int y0 = i*Block;
							int x1 = x0 + w;
							int y1 = y0 + k;
							int temp1 = y1*width + x1;

							ubuf1[temp] = wd[3 * temp1];
							ubuf2[temp] = wd[3 * temp1 + 1];
							ubuf3[temp] = wd[3 * temp1 + 2];
						}
					}

					//ԭʼӰ�񵽷ֿ�Ӱ���ת����ϵ
					pImgOut->GetRasterBand(1)->RasterIO(GF_Write, BlockPosX, BlockPosY, curBlockXSize,
						curBlockYSize, ubuf1, widthx, heighty, gdt, 0, 0);
					pImgOut->GetRasterBand(2)->RasterIO(GF_Write, BlockPosX, BlockPosY, curBlockXSize,
						curBlockYSize, ubuf2, widthx, heighty, gdt, 0, 0);
					pImgOut->GetRasterBand(3)->RasterIO(GF_Write, BlockPosX, BlockPosY, curBlockXSize,
						curBlockYSize, ubuf3, widthx, heighty, gdt, 0, 0);

					delete[]ubuf1; ubuf1 = NULL;
					delete[]ubuf2; ubuf2 = NULL;
					delete[]ubuf3; ubuf3 = NULL;
				}
			}
		}
	}

	if (pImgOut != NULL)
	{
		GDALClose(pImgOut);
		pImgOut = NULL;
	}
}


//GDAL�ֿ鱣�泬�������ɽ��
void saveImage2(int width, int height, int nchannel, GDALDataType gdt, unsigned char *wd)
{
	string filepath = "D://����һ��//����//ICESat-2����У��//����//Ӱ��//mask.tiff";


	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	GDALDriver *poDriver = NULL;
	GDALDataset *pImgOut = NULL;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTIFF");

	if (poDriver != NULL)
	{
		pImgOut = poDriver->Create(filepath.c_str(), width, height, nchannel, gdt, NULL);


		if (pImgOut != NULL)
		{
			int i, j;
			int Block = 2000;
			int xNum = (width - 1) / Block + 1; //�зֳ�xNum����
			int yNum = (height - 1) / Block + 1; //�зֳ�yNum����

			for (i = 0; i<yNum; i++)
			{
				for (j = 0; j<xNum; j++)
				{
					int BlockPosX = 0;	//��ǰ�������ԭͼ����ʼx����
					int BlockPosY = 0;	//��ǰ�������ԭͼ����ʼy����
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

					int widthx = curBlockXSize;
					int heighty = curBlockYSize;

					int m = widthx * heighty;
					unsigned char* ubuf1 = new unsigned char[m];
					unsigned char* ubuf2 = new unsigned char[m];
					unsigned char* ubuf3 = new unsigned char[m];

					for (int k = 0; k < heighty; k++)
					{
						for (int w = 0; w < widthx; w++)
						{
							int temp = k*widthx + w;
							int x0 = j*Block;//ÿ����ʼ�������
							int y0 = i*Block;
							int x1 = x0 + w;
							int y1 = y0 + k;
							int temp1 = y1*width + x1;

							ubuf1[temp] = wd[3 * temp1];
							ubuf2[temp] = wd[3 * temp1 + 1];
							ubuf3[temp] = wd[3 * temp1 + 2];
						}
					}

					//ԭʼӰ�񵽷ֿ�Ӱ���ת����ϵ
					pImgOut->GetRasterBand(1)->RasterIO(GF_Write, BlockPosX, BlockPosY, curBlockXSize,
						curBlockYSize, ubuf1, widthx, heighty, gdt, 0, 0);
					pImgOut->GetRasterBand(2)->RasterIO(GF_Write, BlockPosX, BlockPosY, curBlockXSize,
						curBlockYSize, ubuf2, widthx, heighty, gdt, 0, 0);
					pImgOut->GetRasterBand(3)->RasterIO(GF_Write, BlockPosX, BlockPosY, curBlockXSize,
						curBlockYSize, ubuf3, widthx, heighty, gdt, 0, 0);

					delete[]ubuf1; ubuf1 = NULL;
					delete[]ubuf2; ubuf2 = NULL;
					delete[]ubuf3; ubuf3 = NULL;
				}
			}
		}
	}

	if (pImgOut != NULL)
	{
		GDALClose(pImgOut);
		pImgOut = NULL;
	}
}