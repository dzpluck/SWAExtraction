// seaLandSplit.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include "dealImge.h"


int main(int argc, _TCHAR* argv[])
{
	while (1)
	{
		cout << "���뷴�ݵ�ң��Ӱ�����ݣ�";
		string imgPath;
		cin >> imgPath;

		performSeaLandSplit(imgPath);
		//performSeaLandSplitRecycle(imgPath);
	}
	
	return 0;
}

