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
		int flag = 1; // �� flag ��Ϊ��־�����������û�иı�˵�������ٱ�������
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

	int dis = num / kNum;//���ӵ������
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

//���㷽��
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
	//ndwiClass = new int[num]; //ÿ�������صķ������
	memset(ndwiClass, 0, sizeof(int)*num);

	//��ndwi������������
	float *ndwi = new float[num];
	for (int i = 0; i < num; i++)
	{
		ndwi[i] = ndwiOrigin[i];
	}
	for (int i = 0; i < num - 1; i++)
	{
		int flag = 1; // �� flag ��Ϊ��־�����������û�иı�˵�������ٱ�������
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

	//ǰ�����������Ӧ����
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

	//���ν�������Ӧ����
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

	//���쳣�ֿ�������е���
	float wp1 = (float)(iFlag1 + 1) / (float)num;
	float wp2 = (float)(num - (iFlag2 + iFlag1 + 1)) / (float)num;

	for (int i = 0; i < num; i++)
	{
		if (ndwiOrigin[i]>ndwiFlag2)
		{
			ndwiClass[i] = 2;  //ǳˮ
		}
		else if (ndwiOrigin[i]>ndwiFlag1 && ndwiOrigin[i] <= ndwiFlag2)
		{
			ndwiClass[i] = 1; //��ˮ
		}
		else
		{
			ndwiClass[i] = 0; //½��
		}

		//���쳣������е���
		if (wp1 < 0.01 && ndwiClass[i] == 1)
		{
			ndwiClass[i] = 0;//��ˮ����Ϊ½��
		}
		if (wp2 < 0.01 && ndwiClass[i] == 1)
		{
			ndwiClass[i] = 2;//��ˮ��Ϊǳˮ��
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
	float varValue = 0; //��䷽���м�ֵ����
	float w0 = 0; //ǰ�����ص�����ռ����
	float w1 = 0; //�������ص�����ռ����
	float u0 = 0; //ǰ��ƽ���Ҷ�
	float u1 = 0; //����ƽ���Ҷ�

	for (int i = 0; i<num; i++)
	{
		//ÿ�α���֮ǰ��ʼ��������
		w1 = 0;		u1 = 0;		w0 = 0;		u0 = 0;
		//***********����������ֵ����**************************
		for (int j = 0; j <= i; j++) //�������ָ�ֵ����
		{
			w1 += 1;  //�����������ص�����
			u1 += ndwi[j]; //�������������ܻҶȺ�
		}
		//if (w1 == 0) //�����������ص���Ϊ0ʱ�˳�
		//{
		//	break;
		//}
		u1 = u1 / w1; //��������ƽ���Ҷ�
		w1 = w1 / num; // �����������ص�����ռ����
		//***********����������ֵ����**************************

		//***********ǰ��������ֵ����**************************
		for (int k = i + 1; k<num; k++)
		{
			w0 += 1;  //ǰ���������ص�����
			u0 += ndwi[k]; //ǰ�����������ܻҶȺ�
		}
		//if (w0 == 0) //ǰ���������ص���Ϊ0ʱ�˳�
		//{
		//	break;
		//}
		u0 = u0 / w0; //ǰ������ƽ���Ҷ�
		w0 = w0 / num; // ǰ���������ص�����ռ����
		//***********ǰ��������ֵ����**************************

		//***********��䷽�����******************************
		float varValueI = w0*w1*(u1 - u0)*(u1 - u0); //��ǰ��䷽�����
		if (varValue<varValueI && w1 != 0 && w0 != 0)
		{
			varValue = varValueI;
			iFlag = i;
			ndwiFlag = ndwi[i];
		}
	}

}
 

//��Ӱ�����������
void kmeans(float *ndwi, int num, int kNum, int*& ndwiClass)
{
	//ndwiClass = new int[num]; //ÿ�������صķ������
	memset(ndwiClass, 0, sizeof(int)*num);

	//��ʼ���ӵ�
	float *ndwiSeed = new float[kNum];
	float *ndwiSeedNew = new float[kNum];
	int *ndwiSeedNum = new int[kNum];
	bubbleSortNdwiSeed(ndwi, num, kNum, ndwiSeed);

	int circle = 0; //ѭ���жϴ���
	while (circle < 100000)
	{
		//��ÿ�����ؽ���ѭ�������ж�
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

		//�������ӵ�
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
			ndwiSeedNew[k] /= ndwiSeedNum[k]; //�����µ����ӵ�
		}

		//�ж��¾����ӵ��Ƿ��б仯
		int flag = 0;
		float *distance = new float[kNum];
		for (int k = 0; k < kNum; k++)
		{
			distance[k] = abs(ndwiSeed[k] - ndwiSeedNew[k]);
			if (distance[k]>0.00001)
			{
				flag = 1; //��ʾ���ӵ�û���ȶ�
			}
		}
		if (flag == 1)
		{
			for (int k = 0; k < kNum; k++)
			{
				ndwiSeed[k] = ndwiSeedNew[k]; //�������ӵ�
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

//���������㷨
void areaGrow(int width, int height, int *&classFlag)
{
	vector<int> num;  num.clear();//��¼ÿ��ǳˮ����������	
	int *flag = new int[width*height]; //��¼ÿ��ǳˮ������������ǩ
	memset(flag, 0, sizeof(int)*width*height);

	int sz = width*height;
	static int nDx[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };//����8��������
	static int nDy[8] = { 0, -1, -1, -1, 0, 1, 1, 1 };

	// ��ʼ��,������¼��ǰ���ص������ĸ����;
	// ������־��ǰ���ص���û�б�����
	unsigned char *pUnRegion = new unsigned char[sz];
	memset(pUnRegion, 0, sz*sizeof(unsigned char));

	int nStart;
	int nEnd;
	int nCurrX;
	int nCurrY;
	// ͼ��ĺ�������,�����Ե�ǰ���ص�8������б���
	int xx;
	int yy;
	int areaFlag = 0;
	int areaNum = 0;
	for (int i = 0; i<height; i++)
	{
		for (int j = 0; j<width; j++)
		{
			int temp = i*width + j;

			// �����ջ���洢����
			vector<int> pnGrowQueX; pnGrowQueX.clear();
			vector<int> pnGrowQueY; pnGrowQueY.clear();
			areaNum = 0;
			if (pUnRegion[temp] == 0 && classFlag[temp] == 2) //�����ǰ��û�б�����
			{
				areaFlag++; //���ӵ��ǩ��1
				areaNum++; //�������������ݼ�1
				flag[temp] = areaFlag;

				pUnRegion[temp] = 1;
				//��ʼ��
				nStart = 0;
				nEnd = 0;

				// �����ӵ������ѹ��ջ
				pnGrowQueX.push_back(j);
				pnGrowQueY.push_back(i);

				while (nStart <= nEnd)
				{
					// ��ǰ���ӵ������
					nCurrX = pnGrowQueX[nStart];
					nCurrY = pnGrowQueY[nStart];

					// �Ե�ǰ���4������б���
					for (int k = 0; k<8; k++)
					{
						// 8�������ص�����
						xx = nCurrX + nDx[k];
						yy = nCurrY + nDy[k];
						int temp2 = yy*width + xx;

						// pUnRegion[yy*nWidth+xx]==0 ��ʾ��û�д��� 
						// �����������ж�����(xx��yy)�͵�ǰ����(nCurrX,nCurrY) ����ֵ��ľ���ֵ
						if ((xx < width) && (xx >= 0) && (yy < height) && (yy >= 0) && (pUnRegion[temp2] == 0) && (classFlag[temp2] == 2))
						{
							// ��ջ��β��ָ�����һλ
							nEnd++;

							// ����(xx��yy) ѹ��ջ
							pnGrowQueX.push_back(xx);
							pnGrowQueY.push_back(yy);

							// ͬʱҲ���������ش����
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

	int numThreshold = width*height*0.002; //ǳˮ�����������ֵ��
	//int numThreshold = 125;
	vector<int> liuFlag; liuFlag.clear();
	for (int i = 0; i < num.size(); i++)
	{
		if (num[i] >= numThreshold)
		{
			liuFlag.push_back(i + 1); //������Ҫ��ǳˮ�����ǩ
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
			classFlag[i] = 1;  //����ǳˮ�����Ϊ��ˮ����
		}
	}

	delete[]pUnRegion; pUnRegion = NULL;
	delete[] flag; flag = NULL;
	vector<int>().swap(num);
	vector<int>().swap(liuFlag);
}

//���������㷨
void areaGrowHole(int width, int height, int *&classFlag)
{
	vector<int> num;  num.clear();//��¼ÿ��ǳˮ����������	
	int *flag = new int[width*height]; //��¼ÿ��ǳˮ������������ǩ
	memset(flag, 0, sizeof(int)*width*height);

	int sz = width*height;
	static int nDx[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };//����8��������
	static int nDy[8] = { 0, -1, -1, -1, 0, 1, 1, 1 };
	//static int nDx[4] = { -1,  0,  1,  0 };//����4��������
	//static int nDy[4] = { 0,  -1,  0, 1 };

	// ��ʼ��,������¼��ǰ���ص������ĸ����;
	// ������־��ǰ���ص���û�б�����
	unsigned char *pUnRegion = new unsigned char[sz];
	memset(pUnRegion, 0, sz*sizeof(unsigned char));

	int nStart;
	int nEnd;
	int nCurrX;
	int nCurrY;
	// ͼ��ĺ�������,�����Ե�ǰ���ص�8������б���
	int xx;
	int yy;
	int areaFlag = 0;
	int areaNum = 0;
	for (int i = 0; i<height; i++)
	{
		for (int j = 0; j<width; j++)
		{
			int temp = i*width + j;

			// �����ջ���洢����
			vector<int> pnGrowQueX; pnGrowQueX.clear();
			vector<int> pnGrowQueY; pnGrowQueY.clear();
			areaNum = 0;
			if (pUnRegion[temp] == 0 && classFlag[temp] != 2) //�����ǰ��û�б�����
			{
				areaFlag++; //���ӵ��ǩ��1
				areaNum++; //�������������ݼ�1
				flag[temp] = areaFlag;

				pUnRegion[temp] = 1;
				//��ʼ��
				nStart = 0;
				nEnd = 0;

				// �����ӵ������ѹ��ջ
				pnGrowQueX.push_back(j);
				pnGrowQueY.push_back(i);

				while (nStart <= nEnd)
				{
					// ��ǰ���ӵ������
					nCurrX = pnGrowQueX[nStart];
					nCurrY = pnGrowQueY[nStart];

					// �Ե�ǰ���4������б���
					for (int k = 0; k<8; k++)
					{
						// 8�������ص�����
						xx = nCurrX + nDx[k];
						yy = nCurrY + nDy[k];
						int temp2 = yy*width + xx;

						// pUnRegion[yy*nWidth+xx]==0 ��ʾ��û�д��� 
						// �����������ж�����(xx��yy)�͵�ǰ����(nCurrX,nCurrY) ����ֵ��ľ���ֵ
						if ((xx < width) && (xx >= 0) && (yy < height) && (yy >= 0) && (pUnRegion[temp2] == 0) && (classFlag[temp2] != 2))
						{
							// ��ջ��β��ָ�����һλ
							nEnd++;

							// ����(xx��yy) ѹ��ջ
							pnGrowQueX.push_back(xx);
							pnGrowQueY.push_back(yy);

							// ͬʱҲ���������ش����
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

	int numThreshold = width*height*0.0005; //ǳˮ�����������ֵ��
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
			liuFlag.push_back(i + 1); //������Ҫ��ǳˮ�����ǩ
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
			classFlag[i] = 2;  //ǳˮ�ն������Ϊǳˮ
		}
	}

	delete[]pUnRegion; pUnRegion = NULL;
	delete[] flag; flag = NULL;
	vector<int>().swap(num);
	vector<int>().swap(liuFlag);
}

void expendArea(int width, int height, int *&classFlag)
{
	static int nDx[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };//����8��������
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

	int num = 3; //ѭ�����ʹ���
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

				//��Χ�����д��ڵ���5��Ϊǰ�����أ��������ҲΪǰ������
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

		//������һ��ѭ������
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
	//Ӱ���ȡ�봦��
	//*************************
	int width, height, bands;
	GDALDataType gdt; //��������
	GUInt16 * ubuff;
	double* trans;   //����任����
	OGRSpatialReference projSRS; //ͶӰ������Ϣ
	readImage(imgPath, width, height, bands, gdt, ubuff, trans, projSRS);

	unsigned char* ubuffTemp = new unsigned char[3 * width*height]; //�������ɳ����ؽ���ƽ������
	int *max = new int[3];
	int *min = new int[3];
	max[0] = ubuff[0];
	max[1] = ubuff[1];
	max[2] = ubuff[2];
	min[0] = ubuff[0];
	min[1] = ubuff[1];
	min[2] = ubuff[2];
	float* pBuf = new float[width*height*bands]; //ת����������
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
	//���㷨����ʱ�����ͳ��
	clock_t start, end;
	start = clock(); //��ʼ��ʱ
	//*******************************


	int flag = 0;//0��ʾ�Աȷ����� 1��ʾ����������Ӧ����
	int *classFlag;
	//*******************************
	//���ڵ����ص�NDWI��Ӱ����к�½���ࣨ�Աȷ�����
	//*******************************
	//ndwi(imgPath, pBuf, width, height, classFlag);

	
	//*************************
	//���ڳ����ض�Ӱ���������Ӧ��ֵ������
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
	performSLIC(ubuffTemp, width, height, 1000, klabels, numlabels);//���ɳ�����
	ndwiLandClass(imgPath, pBuf, width, height, klabels, numlabels, classFlag);//���ڳ����ض�Ӱ���������Ӧ��ֵ������
	delete[] klabels; klabels = NULL;

	//*************************
	//����ȡ���ͼ�ν��к�������
	//*************************
	areaGrow(width, height, classFlag);
	areaGrowHole(width, height, classFlag);
	expendArea(width, height, classFlag);


	//*******************************
	end = clock();//������ʱ
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
		//wdPath = imgPath + "_mask_BL.tiff"; //��׼�����
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
	//Ӱ���ȡ�봦��
	//*************************
	int width, height, bands;
	GDALDataType gdt; //��������
	GUInt16 * ubuff;
	double* trans;   //����任����
	OGRSpatialReference projSRS; //ͶӰ������Ϣ
	readImage(imgPath, width, height, bands, gdt, ubuff, trans, projSRS);

	unsigned char* ubuffTemp = new unsigned char[3 * width*height]; //�������ɳ����ؽ���ƽ������
	int *max = new int[3];
	int *min = new int[3];
	max[0] = ubuff[0];
	max[1] = ubuff[1];
	max[2] = ubuff[2];
	min[0] = ubuff[0];
	min[1] = ubuff[1];
	min[2] = ubuff[2];
	float* pBuf = new float[width*height*bands]; //ת����������
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
	//���㷨����ʱ�����ͳ��
	clock_t start, end;
	start = clock(); //��ʼ��ʱ
	//*******************************


	int flag = 0;//0��ʾ�Աȷ����� 1��ʾ����������Ӧ����
	//*******************************
	//���ڵ����ص�NDWI��Ӱ����к�½���ࣨ�Աȷ�����
	//*******************************
	//int *classFlag;
	//ndwi(imgPath, pBuf, width, height, classFlag);


	//*************************
	//���ڳ����ض�Ӱ���������Ӧ��ֵ������
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
		performSLIC(ubuffTemp, width, height, kNum[k], klabels, numlabels);//���ɳ�����
		ndwiLandClass(imgPath, pBuf, width, height, klabels, numlabels, classFlag);//���ڳ����ض�Ӱ���������Ӧ��ֵ������
		delete[] klabels; klabels = NULL;

		//*************************
		//����ȡ���ͼ�ν��к�������
		//*************************
		areaGrow(width, height, classFlag);
		areaGrowHole(width, height, classFlag);
		expendArea(width, height, classFlag);


		//*******************************
		end = clock();//������ʱ
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
			char buffer[64]; // �����㹻�Ŀռ����洢ת������ַ���  
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