#ifndef SLIC_H
#define SLIC_H


#include <vector>
using namespace std;


#define DBL_MAX         1.7976931348623158e+308 /* max value */
///////////////////////////生成SLIC超像素/////////////////////////////
///////////////////////////////////////////////////////////////////////////
//进行图像转换RGB-->LAB

void  RGB2XYZ(const int& sR,
	const int& sG,
	const int& sB,
	double&    X,
	double&    Y,
	double&    Z);
void  RGB2LAB(const int& sR,
	const int& sG,
	const int& sB,
	double&   lval,
	double&   aval,
	double&   bval);
void DoRGBtoLABConversion(
	const unsigned int*&		ubuff,
	double*&					lvec,
	double*&					avec,
	double*&					bvec,
	int                         width,
	int                         height);
//种子点均匀洒在图像上
void GetLABXYSeeds_ForGivenStepSize(
	vector<double>&				kseedsl,
	vector<double>&				kseedsa,
	vector<double>&				kseedsb,
	vector<double>&				kseedsx,
	vector<double>&				kseedsy,
	int                         m_width,
	int                         m_height,
	double*&                     m_lvec,
	double*&                     m_avec,
	double*&                     m_bvec,
	const int&					STEP,
	const bool&					perturbseeds,
	const vector<double>&       edgemag);
//在3*3窗口内调整种子点
void DetectLabEdges(
	const double*				lvec,
	const double*				avec,
	const double*				bvec,
	const int&					width,
	const int&					height,
	vector<double>&				edges);
void PerturbSeeds(
	vector<double>&				kseedsl,
	vector<double>&				kseedsa,
	vector<double>&				kseedsb,
	vector<double>&				kseedsx,
	vector<double>&			    kseedsy,
	double*&                     m_lvec,
	double*&                     m_avec,
	double*&                     m_bvec,
	int                         m_width,
	int                         m_height,
	const vector<double>&                   edges);
//在2S窗口内进行聚类
void PerformSuperpixelSLIC(
	vector<double>&				kseedsl,
	vector<double>&				kseedsa,
	vector<double>&				kseedsb,
	vector<double>&				kseedsx,
	vector<double>&				kseedsy,
	double*&                    m_lvec,
	double*&                    m_avec,
	double*&                    m_bvec,
	int                         m_width,
	int                         m_height,
	int*&						klabels,
	const int&					STEP,
	const vector<double>&		edgemag,
	const double&				m = 10.0);
//后续处理
void EnforceLabelConnectivity(
	const int*					labels,//input labels that need to be corrected to remove stray labels
	const int					width,
	const int					height,
	int*&						nlabels,//new labels
	int&						numlabels,//the number of labels changes in the end if segments are removed
	const int&					K);//the number of superpixels desired by the user
//画出边界
void DrawContoursAroundSegments(
	unsigned char*&			ubuff,
	int*&					labels,
	const int&              lLineBytes,
	const int&				width,
	const int&				height,
	const unsigned int&				color);


//求对生成的超像素灰度均值及方差
void  RGBtoGray(unsigned char *ubuff,//图像数组
	int width,//图像宽度
	int height,//图像高度
	int* klabels,//像素点的标签
	int  seednum,//超像素的分割数
	unsigned char *&gray,//超像素的灰度数组
	float *&SD//超像素的灰度方差
	);

//生成超像素标签
void performSLIC(unsigned char *ubuff, //判断数组
	int width,//图像宽度
	int height,//图像高度
	int k,//超像素的个数
	int *&klabels,//每个像素的超像素的标记
	int &numlabels//超像素的个数
	//unsigned char *&gray,//超像素数组均值
	//float*& SD//超像素的方差
	);


//求每个超像素的邻接超像素
void SuperpixelAdjacent(int *klabels,//每个像素的超像素标签
	int height,//图像高度
	int width,//图像宽度
	int seedNum,//超像素的个数
	vector< vector<int> >& adjacent//超像素的邻接矩阵容器
	);

//邻接超像素进行区域增长
void seedGrow(unsigned char *ubuff,
	int height,
	int width,
	int *klabels,
	int numlabels,
	unsigned char* gray,
	unsigned char cloudThreshold,//区域增长的云的阈值
	vector<int> seedFlag);//经过滤保留下来的超像素的标签



#endif