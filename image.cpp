#include "stdafx.h"
#include "image.h"

//GDAL分块保存海陆分离影像掩膜结果
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
		CPLFree(pszSRS_WKT);//使用完后释放

		if (pImgOut != NULL)
		{
			int i, j;
			int Block = 2000;
			int xNum = (width - 1) / Block + 1; //行分成xNum个块
			int yNum = (height - 1) / Block + 1; //列分成yNum个块

			for (i = 0; i<yNum; i++)
			{
				for (j = 0; j<xNum; j++)
				{
					int BlockPosX = 0;	//当前块相对于原图的起始x坐标
					int BlockPosY = 0;	//当前块相对于原图的起始y坐标
					BlockPosX = Block * j;
					BlockPosY = Block * i;

					int curBlockXSize = Block;		//当前块的 width
					int curBlockYSize = Block;		//当前块的 height

					if (j == xNum - 1)   //列尾
					{
						curBlockXSize = width - Block * j;
					}
					if (i == yNum - 1) //行尾
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
							int x0 = j*Block;//每块起始点的坐标
							int y0 = i*Block;
							int x1 = x0 + w;
							int y1 = y0 + k;
							int temp1 = y1*width + x1;

							ubuf1[temp] = wd[3 * temp1];
							ubuf2[temp] = wd[3 * temp1 + 1];
							ubuf3[temp] = wd[3 * temp1 + 2];
						}
					}

					//原始影像到分块影像的转换关系
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


//GDAL分块保存超像素生成结果
void saveImage2(int width, int height, int nchannel, GDALDataType gdt, unsigned char *wd)
{
	string filepath = "D://海洋一所//论文//ICESat-2折射校正//东岛//影像//mask.tiff";


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
			int xNum = (width - 1) / Block + 1; //行分成xNum个块
			int yNum = (height - 1) / Block + 1; //列分成yNum个块

			for (i = 0; i<yNum; i++)
			{
				for (j = 0; j<xNum; j++)
				{
					int BlockPosX = 0;	//当前块相对于原图的起始x坐标
					int BlockPosY = 0;	//当前块相对于原图的起始y坐标
					BlockPosX = Block * j;
					BlockPosY = Block * i;

					int curBlockXSize = Block;		//当前块的 width
					int curBlockYSize = Block;		//当前块的 height

					if (j == xNum - 1)   //列尾
					{
						curBlockXSize = width - Block * j;
					}
					if (i == yNum - 1) //行尾
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
							int x0 = j*Block;//每块起始点的坐标
							int y0 = i*Block;
							int x1 = x0 + w;
							int y1 = y0 + k;
							int temp1 = y1*width + x1;

							ubuf1[temp] = wd[3 * temp1];
							ubuf2[temp] = wd[3 * temp1 + 1];
							ubuf3[temp] = wd[3 * temp1 + 2];
						}
					}

					//原始影像到分块影像的转换关系
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