#include "stdafx.h"
#include "dealImge.h"
#include "iostream"

using namespace std;

void bubbleSortNdwiSeed(float *ndwiOrigin, int num, int kNum, float *&ndwiSeed)
{
	float * ndwi = new float[num];
	for (int i = 0; i < num; i++)
	{
		ndwi[i] = ndwiOrigin[i];
	}

	for (int i = 0; i < num - 1; i++)
	{
		int flag = 1; // 以 flag 作为标志，如果遍历完没有改变说明有序不再遍历排序
		for (int j = 0; j < num - 1 - i; j++)
		{
			if (ndwi[j]>ndwi[j + 1])
			{
				float temp = ndwi[j];
				ndwi[j] = ndwi[j + 1];
				ndwi[j + 1] = temp;
				flag = 0;
			}
		}
		if (flag)
		{
			break;
		}
	}

	int dis = num / kNum;//种子点距离间隔
	for (int k = 0; k < kNum; k++)
	{
		/*int w = dis*k;
		if (w < 0)
		{
			w = 0;
		}
		if (w > num - 1)
		{
			w = num - 1;
		}*/
		int w = dis*(0.5 + k);
		ndwiSeed[k] = ndwi[w];
	}

	delete[] ndwi; ndwi = NULL;
}

//计算方差
float sd(float*ndwi, int num)
{
	float sum = 0.0;
	for (int i = 0; i < num; i++)
	{
		sum += ndwi[i];
	}
	sum /= num;

	float ndwiSd = 0.0;
	for (int i = 0; i < num; i++)
	{
		ndwiSd += pow(ndwi[i] - sum, 2);
	}
	ndwiSd /= (num - 1);

	return ndwiSd;
}

void ndwiSD(float *ndwi, int num, int &iFlag, float &ndwiFlag)
{
	iFlag = 1;
	float sdFlag = 1000000.0;
	ndwiFlag = 0.0;
	for (int i = 1; i < num - 2; i++)
	{
		float *ndwi11 = new float[i + 1];
		float *ndwi12 = new float[num - i - 1];
		for (int j = 0; j < i + 1; j++)
		{
			ndwi11[j] = ndwi[j];
		}
		for (int j = 0; j < num - i - 1; j++)
		{
			ndwi12[j] = ndwi[i + 1 + j];
		}

		float sd11 = sd(ndwi11, i + 1);
		float sd12 = sd(ndwi12, num - i - 1);
		if (sdFlag > sd11 + sd12)
		{
			sdFlag = sd11 + sd12;
			iFlag = i;
			ndwiFlag = ndwi[i];
		}

		delete[]ndwi11; ndwi11 = NULL;
		delete[]ndwi12; ndwi12 = NULL;
	}
}

void otsuClass(float* ndwiOrigin, int num, int *&ndwiClass)
{
	//ndwiClass = new int[num]; //每个超像素的分类类别
	memset(ndwiClass, 0, sizeof(int)*num);

	//对ndwi进行升序排序
	float *ndwi = new float[num];
	for (int i = 0; i < num; i++)
	{
		ndwi[i] = ndwiOrigin[i];
	}
	for (int i = 0; i < num - 1; i++)
	{
		int flag = 1; // 以 flag 作为标志，如果遍历完没有改变说明有序不再遍历排序
		for (int j = 0; j < num - 1 - i; j++)
		{
			if (ndwi[j]>ndwi[j + 1])
			{
				float temp = ndwi[j];
				ndwi[j] = ndwi[j + 1];
				ndwi[j + 1] = temp;
				flag = 0;
			}
		}
		if (flag)
		{
			break;
		}
	}

	//前半段数组自适应处理
	int num1 = num / 2;
	float* ndwi1 = new float[num1];
	for (int i = 0; i < num1; i++)
	{
		ndwi1[i] = ndwi[i];
	}
	int iFlag1;
	float ndwiFlag1;
	//ndwiSD(ndwi1, num1, iFlag1, ndwiFlag1);
	ndwiOTSU(ndwi1, num1, iFlag1, ndwiFlag1);
	cout << iFlag1 << "  " << ndwiFlag1 << endl;

	//后半段进行自适应处理
	int num2 = num - iFlag1 - 1;
	float* ndwi2 = new float[num2];
	for (int i = 0; i < num2; i++)
	{
		ndwi2[i] = ndwi[iFlag1 + 1 + i];
	}
	int iFlag2;
	float ndwiFlag2;
	//ndwiSD(ndwi2, num2, iFlag2, ndwiFlag2);	
	ndwiOTSU(ndwi2, num2, iFlag2, ndwiFlag2);
	cout << iFlag2 + iFlag1 + 1 << "  " << ndwiFlag2 << endl;

	//对异常分块情况进行调整
	float wp1 = (float)(iFlag1 + 1) / (float)num;
	float wp2 = (float)(num - (iFlag2 + iFlag1 + 1)) / (float)num;

	for (int i = 0; i < num; i++)
	{
		if (ndwiOrigin[i]>ndwiFlag2)
		{
			ndwiClass[i] = 2;  //浅水
		}
		else if (ndwiOrigin[i]>ndwiFlag1 && ndwiOrigin[i] <= ndwiFlag2)
		{
			ndwiClass[i] = 1; //深水
		}
		else
		{
			ndwiClass[i] = 0; //陆地
		}

		//对异常情况进行调整
		if (wp1 < 0.01 && ndwiClass[i] == 1)
		{
			ndwiClass[i] = 0;//深水区改为陆地
		}
		if (wp2 < 0.01 && ndwiClass[i] == 1)
		{
			ndwiClass[i] = 2;//深水改为浅水区
		}
	}

	delete[] ndwi1; ndwi1 = NULL;
	delete[] ndwi2; ndwi2 = NULL;
	delete[] ndwi; ndwi = NULL;

	/*int iFlag = 1;
	float sdFlag = 1000000.0;
	float ndwiFlag = 0.0;
	for (int i = 1; i < num1 - 2; i++)
	{
		float *ndwi11 = new float[i + 1];
		float *ndwi12 = new float[num1 - i - 1];
		for (int j = 0; j < i + 1; j++)
		{
			ndwi11[j] = ndwi1[j];
		}
		for (int j = 0; j < num1 - i - 1; j++)
		{
			ndwi12[j] = ndwi1[i + 1 + j];
		}

		float sd11 = sd(ndwi11, i + 1);
		float sd12 = sd(ndwi12, num1 - i - 1);
		if (sdFlag > sd11 + sd12)
		{
			sdFlag = sd11 + sd12;
			iFlag = i;
			ndwiFlag = ndwi1[i];
		}

		delete[]ndwi11; ndwi11 = NULL;
		delete[]ndwi12; ndwi12 = NULL;
	}*/
}

void ndwiOTSU(float *ndwi, int num, int &iFlag, float &ndwiFlag)
{
	float varValue = 0; //类间方差中间值保存
	float w0 = 0; //前景像素点数所占比例
	float w1 = 0; //背景像素点数所占比例
	float u0 = 0; //前景平均灰度
	float u1 = 0; //背景平均灰度

	for (int i = 0; i<num; i++)
	{
		//每次遍历之前初始化各变量
		w1 = 0;		u1 = 0;		w0 = 0;		u0 = 0;
		//***********背景各分量值计算**************************
		for (int j = 0; j <= i; j++) //背景部分各值计算
		{
			w1 += 1;  //背景部分像素点总数
			u1 += ndwi[j]; //背景部分像素总灰度和
		}
		//if (w1 == 0) //背景部分像素点数为0时退出
		//{
		//	break;
		//}
		u1 = u1 / w1; //背景像素平均灰度
		w1 = w1 / num; // 背景部分像素点数所占比例
		//***********背景各分量值计算**************************

		//***********前景各分量值计算**************************
		for (int k = i + 1; k<num; k++)
		{
			w0 += 1;  //前景部分像素点总数
			u0 += ndwi[k]; //前景部分像素总灰度和
		}
		//if (w0 == 0) //前景部分像素点数为0时退出
		//{
		//	break;
		//}
		u0 = u0 / w0; //前景像素平均灰度
		w0 = w0 / num; // 前景部分像素点数所占比例
		//***********前景各分量值计算**************************

		//***********类间方差计算******************************
		float varValueI = w0*w1*(u1 - u0)*(u1 - u0); //当前类间方差计算
		if (varValue<varValueI && w1 != 0 && w0 != 0)
		{
			varValue = varValueI;
			iFlag = i;
			ndwiFlag = ndwi[i];
		}
	}

}
 

//对影像进行三分类
void kmeans(float *ndwi, int num, int kNum, int*& ndwiClass)
{
	//ndwiClass = new int[num]; //每个超像素的分类类别
	memset(ndwiClass, 0, sizeof(int)*num);

	//初始种子点
	float *ndwiSeed = new float[kNum];
	float *ndwiSeedNew = new float[kNum];
	int *ndwiSeedNum = new int[kNum];
	bubbleSortNdwiSeed(ndwi, num, kNum, ndwiSeed);

	int circle = 0; //循环判断次数
	while (circle < 100000)
	{
		//对每个像素进行循环分类判断
		for (int i = 0; i < num; i++)
		{
			float temp = 100000.0;
			for (int j = 0; j < kNum; j++)
			{
				if (abs(ndwi[i] - ndwiSeed[j]) < temp)
				{
					ndwiClass[i] = j;
					temp = abs(ndwi[i] - ndwiSeed[j]);
				}
			}
		}

		//更新种子点
		for (int k = 0; k < kNum; k++)
		{
			ndwiSeedNew[k] = 0.0;
			ndwiSeedNum[k] = 0;
		}
		for (int i = 0; i < num; i++)
		{
			ndwiSeedNew[ndwiClass[i]] += ndwi[i];
			ndwiSeedNum[ndwiClass[i]]++;
		}
		for (int k = 0; k < kNum; k++)
		{
			ndwiSeedNew[k] /= ndwiSeedNum[k]; //计算新的种子点
		}

		//判断新旧种子点是否有变化
		int flag = 0;
		float *distance = new float[kNum];
		for (int k = 0; k < kNum; k++)
		{
			distance[k] = abs(ndwiSeed[k] - ndwiSeedNew[k]);
			if (distance[k]>0.00001)
			{
				flag = 1; //表示种子点没有稳定
			}
		}
		if (flag == 1)
		{
			for (int k = 0; k < kNum; k++)
			{
				ndwiSeed[k] = ndwiSeedNew[k]; //更新种子点
			}
		}
		else
		{
			break;
		}

	}//end while

	delete[] ndwiSeed; ndwiSeed = NULL;
	delete[] ndwiSeedNew; ndwiSeedNew = NULL;
	delete[] ndwiSeedNum; ndwiSeedNum = NULL;
}

//区域增长算法
void areaGrow(int width, int height, int *&classFlag)
{
	vector<int> num;  num.clear();//记录每块浅水区域像素数	
	int *flag = new int[width*height]; //记录每块浅水区域像素类别标签
	memset(flag, 0, sizeof(int)*width*height);

	int sz = width*height;
	static int nDx[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };//定义8领域数组
	static int nDy[8] = { 0, -1, -1, -1, 0, 1, 1, 1 };

	// 初始化,用来记录当前像素点属于哪个类别;
	// 用来标志当前像素点有没有被处理
	unsigned char *pUnRegion = new unsigned char[sz];
	memset(pUnRegion, 0, sz*sizeof(unsigned char));

	int nStart;
	int nEnd;
	int nCurrX;
	int nCurrY;
	// 图象的横纵坐标,用来对当前象素的8邻域进行遍历
	int xx;
	int yy;
	int areaFlag = 0;
	int areaNum = 0;
	for (int i = 0; i<height; i++)
	{
		for (int j = 0; j<width; j++)
		{
			int temp = i*width + j;

			// 定义堆栈，存储坐标
			vector<int> pnGrowQueX; pnGrowQueX.clear();
			vector<int> pnGrowQueY; pnGrowQueY.clear();
			areaNum = 0;
			if (pUnRegion[temp] == 0 && classFlag[temp] == 2) //如果当前点没有被处理
			{
				areaFlag++; //种子点标签加1
				areaNum++; //该区域像素数据加1
				flag[temp] = areaFlag;

				pUnRegion[temp] = 1;
				//初始化
				nStart = 0;
				nEnd = 0;

				// 把种子点的坐标压入栈
				pnGrowQueX.push_back(j);
				pnGrowQueY.push_back(i);

				while (nStart <= nEnd)
				{
					// 当前种子点的坐标
					nCurrX = pnGrowQueX[nStart];
					nCurrY = pnGrowQueY[nStart];

					// 对当前点的4邻域进行遍历
					for (int k = 0; k<8; k++)
					{
						// 8邻域象素的坐标
						xx = nCurrX + nDx[k];
						yy = nCurrY + nDy[k];
						int temp2 = yy*width + xx;

						// pUnRegion[yy*nWidth+xx]==0 表示还没有处理 
						// 生长条件：判断象素(xx，yy)和当前象素(nCurrX,nCurrY) 象素值差的绝对值
						if ((xx < width) && (xx >= 0) && (yy < height) && (yy >= 0) && (pUnRegion[temp2] == 0) && (classFlag[temp2] == 2))
						{
							// 堆栈的尾部指针后移一位
							nEnd++;

							// 像素(xx，yy) 压入栈
							pnGrowQueX.push_back(xx);
							pnGrowQueY.push_back(yy);

							// 同时也表明该象素处理过
							pUnRegion[temp2] = 1;
							
							flag[temp2] = areaFlag;
							areaNum++;
						}
					}//end for k

					nStart++;
				}//end while

				num.push_back(areaNum);
			}//end if

			vector <int>().swap(pnGrowQueX);
			vector <int>().swap(pnGrowQueY);

		}//end for j
	}//end for i

	int numThreshold = width*height*0.002; //浅水区域块像素阈值数
	//int numThreshold = 125;
	vector<int> liuFlag; liuFlag.clear();
	for (int i = 0; i < num.size(); i++)
	{
		if (num[i] >= numThreshold)
		{
			liuFlag.push_back(i + 1); //保留主要的浅水区域标签
		}
	}

	for (int i = 0; i < sz; i++)
	{
		int count = 0;
		for (int j = 0; j < liuFlag.size(); j++)
		{
			if (flag[i] == liuFlag[j])
			{
				count = 1;
				break;
			}
		}
		if (count != 1 && classFlag[i]==2)
		{
			classFlag[i] = 1;  //噪声浅水区域变为深水区域
		}
	}

	delete[]pUnRegion; pUnRegion = NULL;
	delete[] flag; flag = NULL;
	vector<int>().swap(num);
	vector<int>().swap(liuFlag);
}

//区域增长算法
void areaGrowHole(int width, int height, int *&classFlag)
{
	vector<int> num;  num.clear();//记录每块浅水区域像素数	
	int *flag = new int[width*height]; //记录每块浅水区域像素类别标签
	memset(flag, 0, sizeof(int)*width*height);

	int sz = width*height;
	static int nDx[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };//定义8领域数组
	static int nDy[8] = { 0, -1, -1, -1, 0, 1, 1, 1 };
	//static int nDx[4] = { -1,  0,  1,  0 };//定义4领域数组
	//static int nDy[4] = { 0,  -1,  0, 1 };

	// 初始化,用来记录当前像素点属于哪个类别;
	// 用来标志当前像素点有没有被处理
	unsigned char *pUnRegion = new unsigned char[sz];
	memset(pUnRegion, 0, sz*sizeof(unsigned char));

	int nStart;
	int nEnd;
	int nCurrX;
	int nCurrY;
	// 图象的横纵坐标,用来对当前象素的8邻域进行遍历
	int xx;
	int yy;
	int areaFlag = 0;
	int areaNum = 0;
	for (int i = 0; i<height; i++)
	{
		for (int j = 0; j<width; j++)
		{
			int temp = i*width + j;

			// 定义堆栈，存储坐标
			vector<int> pnGrowQueX; pnGrowQueX.clear();
			vector<int> pnGrowQueY; pnGrowQueY.clear();
			areaNum = 0;
			if (pUnRegion[temp] == 0 && classFlag[temp] != 2) //如果当前点没有被处理
			{
				areaFlag++; //种子点标签加1
				areaNum++; //该区域像素数据加1
				flag[temp] = areaFlag;

				pUnRegion[temp] = 1;
				//初始化
				nStart = 0;
				nEnd = 0;

				// 把种子点的坐标压入栈
				pnGrowQueX.push_back(j);
				pnGrowQueY.push_back(i);

				while (nStart <= nEnd)
				{
					// 当前种子点的坐标
					nCurrX = pnGrowQueX[nStart];
					nCurrY = pnGrowQueY[nStart];

					// 对当前点的4邻域进行遍历
					for (int k = 0; k<8; k++)
					{
						// 8邻域象素的坐标
						xx = nCurrX + nDx[k];
						yy = nCurrY + nDy[k];
						int temp2 = yy*width + xx;

						// pUnRegion[yy*nWidth+xx]==0 表示还没有处理 
						// 生长条件：判断象素(xx，yy)和当前象素(nCurrX,nCurrY) 象素值差的绝对值
						if ((xx < width) && (xx >= 0) && (yy < height) && (yy >= 0) && (pUnRegion[temp2] == 0) && (classFlag[temp2] != 2))
						{
							// 堆栈的尾部指针后移一位
							nEnd++;

							// 像素(xx，yy) 压入栈
							pnGrowQueX.push_back(xx);
							pnGrowQueY.push_back(yy);

							// 同时也表明该象素处理过
							pUnRegion[temp2] = 1;

							flag[temp2] = areaFlag;
							areaNum++;
						}
					}//end for k

					nStart++;
				}//end while

				num.push_back(areaNum);
			}//end if

			vector <int>().swap(pnGrowQueX);
			vector <int>().swap(pnGrowQueY);

		}//end for j
	}//end for i

	int numThreshold = width*height*0.0005; //浅水区域块像素阈值数
	//int numThreshold = 500;
	if (numThreshold < 250)
	{
		numThreshold = 250;
	}
	cout << "The number of Hole is :" << numThreshold << endl;
	vector<int> liuFlag; liuFlag.clear();
	for (int i = 0; i < num.size(); i++)
	{
		if (num[i] < numThreshold)
		{
			liuFlag.push_back(i + 1); //保留主要的浅水区域标签
		}
	}

	for (int i = 0; i < sz; i++)
	{
		int count = 0;
		for (int j = 0; j < liuFlag.size(); j++)
		{
			if (flag[i] == liuFlag[j])
			{
				count = 1;
				break;
			}
		}
		if (count == 1 && classFlag[i] != 2)
		{
			classFlag[i] = 2;  //浅水空洞区域填补为浅水
		}
	}

	delete[]pUnRegion; pUnRegion = NULL;
	delete[] flag; flag = NULL;
	vector<int>().swap(num);
	vector<int>().swap(liuFlag);
}

void expendArea(int width, int height, int *&classFlag)
{
	static int nDx[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };//定义8领域数组
	static int nDy[8] = { 0, -1, -1, -1, 0, 1, 1, 1 };

	int *newClassFlag = new int[width*height];
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int temp = i*width + j;
			newClassFlag[temp] = classFlag[temp];
		}
	}

	int num = 3; //循环膨胀次数
	int n = 0;
	while (n < num)
	{
		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				int temp = i*width + j;
				//if (classFlag[temp] == 2)
				//{
				//	for (int k = 0; k < 8; k++)
				//	{
				//		int xx = j + nDx[k];
				//		int yy = i + nDy[k];
				//		if (xx >= 0 && xx < width && yy >= 0 && yy < height)
				//		{
				//			int temp1 = yy*width + xx;
				//			newClassFlag[temp1] = 2;
				//		}
				//	}//end for k	
				//}

				//周围像素有大于等于5个为前景像素，则该像素也为前景像素
				if (classFlag[temp] != 2)
				{
					int count = 0;
					for (int k = 0; k < 8; k++)
					{
						int xx = j + nDx[k];
						int yy = i + nDy[k];
						if (xx >= 0 && xx < width && yy >= 0 && yy < height)
						{
							int temp1 = yy*width + xx;
							if (classFlag[temp1] == 2)
							{
								count++;
							}
						}
					}
					if (count >= 5)
					{
						newClassFlag[temp] = 2;
					}
				}
			}
		}

		//进行下一轮循环迭代
		n++;
		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				int temp = i*width + j;
				classFlag[temp] = newClassFlag[temp];
			}
		}
	}

	delete[] newClassFlag; newClassFlag = NULL;
}



void performSeaLandSplit(string imgPath)
{
	//*************************
	//影像读取与处理
	//*************************
	int width, height, bands;
	GDALDataType gdt; //数据类型
	GUInt16 * ubuff;
	double* trans;   //仿射变换参数
	OGRSpatialReference projSRS; //投影坐标信息
	readImage(imgPath, width, height, bands, gdt, ubuff, trans, projSRS);

	unsigned char* ubuffTemp = new unsigned char[3 * width*height]; //用于生成超像素进行平滑处理
	int *max = new int[3];
	int *min = new int[3];
	max[0] = ubuff[0];
	max[1] = ubuff[1];
	max[2] = ubuff[2];
	min[0] = ubuff[0];
	min[1] = ubuff[1];
	min[2] = ubuff[2];
	float* pBuf = new float[width*height*bands]; //转换到辐射率
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int temp = i*width + j;
			pBuf[4 * temp] = ((float)ubuff[4 * temp] * 0.0001) - 0.1;
			pBuf[4 * temp + 1] = ((float)ubuff[4 * temp + 1] * 0.0001) - 0.1;
			pBuf[4 * temp + 2] = ((float)ubuff[4 * temp + 2] * 0.0001) - 0.1;
			pBuf[4 * temp + 3] = ((float)ubuff[4 * temp + 3] * 0.0001) - 0.1;

			for (int k = 0; k < 3; k++)
			{
				if (max[k] < ubuff[4 * temp + k])
				{
					max[k] = ubuff[4 * temp + k];
				}
				if (min[k] > ubuff[4 * temp + k])
				{
					min[k] = ubuff[4 * temp + k];
				}
			}
		}
	}


	//*******************************
	//对算法消耗时间进行统计
	clock_t start, end;
	start = clock(); //开始计时
	//*******************************


	int flag = 0;//0表示对比方法， 1表示超像素自适应方法
	int *classFlag;
	//*******************************
	//基于单像素的NDWI对影像进行海陆分类（对比方法）
	//*******************************
	//ndwi(imgPath, pBuf, width, height, classFlag);

	
	//*************************
	//基于超像素对影像进行自适应阈值三分类
	//*************************
	flag = 1;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int temp = i*width + j;
			for (int k = 0; k < 3; k++)
			{
				float tempValue = ((float)(ubuff[4 * temp + k] - min[k])) / ((float)(max[k] - min[k]));
				ubuffTemp[3 * temp + (2 - k)] = (unsigned char)(tempValue * 255 > 255 ? 255 : tempValue * 255);
			}
		}
	}
	int *klabels;
	int numlabels;
	performSLIC(ubuffTemp, width, height, 1000, klabels, numlabels);//生成超像素
	ndwiLandClass(imgPath, pBuf, width, height, klabels, numlabels, classFlag);//基于超像素对影像进行自适应阈值三分类
	delete[] klabels; klabels = NULL;

	//*************************
	//对提取后的图形进行后续处理
	//*************************
	areaGrow(width, height, classFlag);
	areaGrowHole(width, height, classFlag);
	expendArea(width, height, classFlag);


	//*******************************
	end = clock();//结束计时
	float duration = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "The consuming time of shallow-sea detection is : " << duration << endl;
	//*******************************


	unsigned char *color = new unsigned char[width*height * 3];
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int temp = i*width + j;
			if (classFlag[temp] == 0) 
			{
				color[3 * temp] = 0;
				color[3 * temp + 1] = 0;
				color[3 * temp + 2] = 0;
			}
			if (classFlag[temp] == 1)
			{
				color[3 * temp] = 0;
				color[3 * temp + 1] = 0;
				color[3 * temp + 2] = 0;
			}
			if (classFlag[temp] == 2)
			{
				color[3 * temp] = 0;
				color[3 * temp + 1] = 255;
				color[3 * temp + 2] = 0;
			}
		}
	}

	string wdPath;
	if (flag == 0)
	{
		wdPath = imgPath + "_mask2.tiff";
	}
	else
	{
		wdPath = imgPath + "_mask.tiff";
		//wdPath = imgPath + "_mask_BL.tiff"; //基准检测结果
		//wdPath = imgPath + "_13.tiff";
	}

	saveImageMask(wdPath, width, height, 3, GDT_Byte, trans, projSRS, color);

	delete[] ubuff; ubuff = NULL;
	delete[] ubuffTemp; ubuffTemp = NULL;
	delete[] max; max = NULL;
	delete[] min; min = NULL;
	delete[] pBuf; pBuf = NULL;
	delete[] trans; trans = NULL;
	delete[] classFlag; classFlag = NULL;
	delete[] color; color = NULL;
}

void performSeaLandSplitRecycle(string imgPath)
{
	//*************************
	//影像读取与处理
	//*************************
	int width, height, bands;
	GDALDataType gdt; //数据类型
	GUInt16 * ubuff;
	double* trans;   //仿射变换参数
	OGRSpatialReference projSRS; //投影坐标信息
	readImage(imgPath, width, height, bands, gdt, ubuff, trans, projSRS);

	unsigned char* ubuffTemp = new unsigned char[3 * width*height]; //用于生成超像素进行平滑处理
	int *max = new int[3];
	int *min = new int[3];
	max[0] = ubuff[0];
	max[1] = ubuff[1];
	max[2] = ubuff[2];
	min[0] = ubuff[0];
	min[1] = ubuff[1];
	min[2] = ubuff[2];
	float* pBuf = new float[width*height*bands]; //转换到辐射率
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int temp = i*width + j;
			pBuf[4 * temp] = ((float)ubuff[4 * temp] * 0.0001) - 0.1;
			pBuf[4 * temp + 1] = ((float)ubuff[4 * temp + 1] * 0.0001) - 0.1;
			pBuf[4 * temp + 2] = ((float)ubuff[4 * temp + 2] * 0.0001) - 0.1;
			pBuf[4 * temp + 3] = ((float)ubuff[4 * temp + 3] * 0.0001) - 0.1;

			for (int k = 0; k < 3; k++)
			{
				if (max[k] < ubuff[4 * temp + k])
				{
					max[k] = ubuff[4 * temp + k];
				}
				if (min[k] > ubuff[4 * temp + k])
				{
					min[k] = ubuff[4 * temp + k];
				}
			}
		}
	}


	//*******************************
	//对算法消耗时间进行统计
	clock_t start, end;
	start = clock(); //开始计时
	//*******************************


	int flag = 0;//0表示对比方法， 1表示超像素自适应方法
	//*******************************
	//基于单像素的NDWI对影像进行海陆分类（对比方法）
	//*******************************
	//int *classFlag;
	//ndwi(imgPath, pBuf, width, height, classFlag);


	//*************************
	//基于超像素对影像进行自适应阈值三分类
	//*************************
	flag = 1;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int temp = i*width + j;
			for (int k = 0; k < 3; k++)
			{
				float tempValue = ((float)(ubuff[4 * temp + k] - min[k])) / ((float)(max[k] - min[k]));
				ubuffTemp[3 * temp + (2 - k)] = (unsigned char)(tempValue * 255 > 255 ? 255 : tempValue * 255);
			}
		}
	}
	int kNum[7] = {  16, 25, 36, 49, 64, 81, 100 };
	for (int k = 0; k < 7; k++)
	{
		int *classFlag;
		int *klabels;
		int numlabels;
		performSLIC(ubuffTemp, width, height, kNum[k], klabels, numlabels);//生成超像素
		ndwiLandClass(imgPath, pBuf, width, height, klabels, numlabels, classFlag);//基于超像素对影像进行自适应阈值三分类
		delete[] klabels; klabels = NULL;

		//*************************
		//对提取后的图形进行后续处理
		//*************************
		areaGrow(width, height, classFlag);
		areaGrowHole(width, height, classFlag);
		expendArea(width, height, classFlag);


		//*******************************
		end = clock();//结束计时
		float duration = (double)(end - start) / CLOCKS_PER_SEC;
		cout << "The consuming time of shallow-sea detection is : " << duration << endl;
		//*******************************


		unsigned char *color = new unsigned char[width*height * 3];
		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				int temp = i*width + j;
				if (classFlag[temp] == 0)
				{
					color[3 * temp] = 0;
					color[3 * temp + 1] = 0;
					color[3 * temp + 2] = 0;
				}
				if (classFlag[temp] == 1)
				{
					color[3 * temp] = 0;
					color[3 * temp + 1] = 0;
					color[3 * temp + 2] = 0;
				}
				if (classFlag[temp] == 2)
				{
					color[3 * temp] = 0;
					color[3 * temp + 1] = 255;
					color[3 * temp + 2] = 0;
				}
			}
		}

		string wdPath;
		if (flag == 0)
		{
			wdPath = imgPath + "_mask2.tiff";
		}
		else
		{
			char buffer[64]; // 分配足够的空间来存储转换后的字符串  
			_itoa(kNum[k], buffer, 10);
			wdPath = imgPath + "_"+ buffer + ".tiff";
		}

		saveImageMask(wdPath, width, height, 3, GDT_Byte, trans, projSRS, color);

		delete[] color; color = NULL;
		delete[] classFlag; classFlag = NULL;
	}
	
	delete[] ubuff; ubuff = NULL;
	delete[] ubuffTemp; ubuffTemp = NULL;
	delete[] max; max = NULL;
	delete[] min; min = NULL;
	delete[] pBuf; pBuf = NULL;
	delete[] trans; trans = NULL;
	//delete[] classFlag; classFlag = NULL;
}