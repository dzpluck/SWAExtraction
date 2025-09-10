#ifndef SLIC_H
#define SLIC_H


#include <vector>
using namespace std;


#define DBL_MAX         1.7976931348623158e+308 /* max value */
///////////////////////////����SLIC������/////////////////////////////
///////////////////////////////////////////////////////////////////////////
//����ͼ��ת��RGB-->LAB

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
//���ӵ��������ͼ����
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
//��3*3�����ڵ������ӵ�
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
//��2S�����ڽ��о���
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
//��������
void EnforceLabelConnectivity(
	const int*					labels,//input labels that need to be corrected to remove stray labels
	const int					width,
	const int					height,
	int*&						nlabels,//new labels
	int&						numlabels,//the number of labels changes in the end if segments are removed
	const int&					K);//the number of superpixels desired by the user
//�����߽�
void DrawContoursAroundSegments(
	unsigned char*&			ubuff,
	int*&					labels,
	const int&              lLineBytes,
	const int&				width,
	const int&				height,
	const unsigned int&				color);


//������ɵĳ����ػҶȾ�ֵ������
void  RGBtoGray(unsigned char *ubuff,//ͼ������
	int width,//ͼ����
	int height,//ͼ��߶�
	int* klabels,//���ص�ı�ǩ
	int  seednum,//�����صķָ���
	unsigned char *&gray,//�����صĻҶ�����
	float *&SD//�����صĻҶȷ���
	);

//���ɳ����ر�ǩ
void performSLIC(unsigned char *ubuff, //�ж�����
	int width,//ͼ����
	int height,//ͼ��߶�
	int k,//�����صĸ���
	int *&klabels,//ÿ�����صĳ����صı��
	int &numlabels//�����صĸ���
	//unsigned char *&gray,//�����������ֵ
	//float*& SD//�����صķ���
	);


//��ÿ�������ص��ڽӳ�����
void SuperpixelAdjacent(int *klabels,//ÿ�����صĳ����ر�ǩ
	int height,//ͼ��߶�
	int width,//ͼ����
	int seedNum,//�����صĸ���
	vector< vector<int> >& adjacent//�����ص��ڽӾ�������
	);

//�ڽӳ����ؽ�����������
void seedGrow(unsigned char *ubuff,
	int height,
	int width,
	int *klabels,
	int numlabels,
	unsigned char* gray,
	unsigned char cloudThreshold,//�����������Ƶ���ֵ
	vector<int> seedFlag);//�����˱��������ĳ����صı�ǩ



#endif