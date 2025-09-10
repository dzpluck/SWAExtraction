// seaLandSplit.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "dealImge.h"


int main(int argc, _TCHAR* argv[])
{
	while (1)
	{
		cout << "输入反演的遥感影像数据：";
		string imgPath;
		cin >> imgPath;

		performSeaLandSplit(imgPath);
		//performSeaLandSplitRecycle(imgPath);
	}
	
	return 0;
}

