#include "stdafx.h"
#include "slic.h"
#include <math.h>
#include "iostream"
#include "image.h"

using namespace std;

double max(double a, double b)
{
	if (a >= b)
		return a;
	else
		return b;
}
double min(double a, double b)
{
	if (a <= b)
		return a;
	else
		return b;
}

//==============================================================================
///	RGB2XYZ
///
/// sRGB (D65 illuninant assumption) to XYZ conversion
//==============================================================================
void  RGB2XYZ(
	const int&		sR,
	const int&		sG,
	const int&		sB,
	double&			X,
	double&			Y,
	double&			Z)
{
	double R = sR / 255.0;
	double G = sG / 255.0;
	double B = sB / 255.0;

	double r, g, b;

	if (R <= 0.04045)	r = R / 12.92;
	else				r = pow((R + 0.055) / 1.055, 2.4);
	if (G <= 0.04045)	g = G / 12.92;
	else				g = pow((G + 0.055) / 1.055, 2.4);
	if (B <= 0.04045)	b = B / 12.92;
	else				b = pow((B + 0.055) / 1.055, 2.4);

	X = r*0.4124564 + g*0.3575761 + b*0.1804375;
	Y = r*0.2126729 + g*0.7151522 + b*0.0721750;
	Z = r*0.0193339 + g*0.1191920 + b*0.9503041;
}

//===========================================================================
///	RGB2LAB
//===========================================================================
void RGB2LAB(const int& sR, const int& sG, const int& sB, double& lval, double& aval, double& bval)
{
	//------------------------
	// sRGB to XYZ conversion
	//------------------------
	double X, Y, Z;
	::RGB2XYZ(sR, sG, sB, X, Y, Z);

	//------------------------
	// XYZ to LAB conversion
	//------------------------
	double epsilon = 0.008856;	//actual CIE standard
	double kappa = 903.3;		//actual CIE standard

	double Xr = 0.950456;	//reference white
	double Yr = 1.0;		//reference white
	double Zr = 1.088754;	//reference white

	double xr = X / Xr;
	double yr = Y / Yr;
	double zr = Z / Zr;

	double fx, fy, fz;
	if (xr > epsilon)	fx = pow(xr, 1.0 / 3.0);
	else				fx = (kappa*xr + 16.0) / 116.0;
	if (yr > epsilon)	fy = pow(yr, 1.0 / 3.0);
	else				fy = (kappa*yr + 16.0) / 116.0;
	if (zr > epsilon)	fz = pow(zr, 1.0 / 3.0);
	else				fz = (kappa*zr + 16.0) / 116.0;

	lval = 116.0*fy - 16.0;
	aval = 500.0*(fx - fy);
	bval = 200.0*(fy - fz);
}

//===========================================================================
///	DoRGBtoLABConversion
///
///	For whole image: overlaoded floating point version
//===========================================================================
void DoRGBtoLABConversion(
	const unsigned int*&		ubuff,
	double*&					lvec,
	double*&					avec,
	double*&					bvec,
	int                         width,
	int                         height)
{
	int sz = width*height;
	lvec = new double[sz];
	avec = new double[sz];
	bvec = new double[sz];
	int k = 0;
	int i, j;
	for (i = 0; i<height; i++)
	{
		for (j = 0; j<width; j++)
		{
			int temp = i*width + j;
			int r = ubuff[3 * temp];
			int g = ubuff[3 * temp + 1];
			int b = ubuff[3 * temp + 2];
			::RGB2LAB(r, g, b, lvec[k], avec[k], bvec[k]);
			k++;
		}
	}
}
//===========================================================================
///	GetLABXYSeeds_ForGivenStepSize
///
/// The k seed values are taken as uniform spatial pixel samples.
//===========================================================================
void GetLABXYSeeds_ForGivenStepSize(
	vector<double>&				kseedsl,
	vector<double>&				kseedsa,
	vector<double>&				kseedsb,
	vector<double>&				kseedsx,
	vector<double>&				kseedsy,
	int                         m_width,
	int                         m_height,
	double*&                    m_lvec,
	double*&                    m_avec,
	double*&                    m_bvec,
	const int&					STEP,
	const bool&					perturbseeds,
	const vector<double>&       edgemag)
{
	const bool hexgrid = false;
	int numseeds(0);
	int n(0);

	//int xstrips = m_width/STEP;
	//int ystrips = m_height/STEP;
	int xstrips = (0.5 + double(m_width) / double(STEP));
	int ystrips = (0.5 + double(m_height) / double(STEP));

	int xerr = m_width - STEP*xstrips; if (xerr < 0){ xstrips--; xerr = m_width - STEP*xstrips; }
	int yerr = m_height - STEP*ystrips; if (yerr < 0){ ystrips--; yerr = m_height - STEP*ystrips; }

	double xerrperstrip = double(xerr) / double(xstrips);
	double yerrperstrip = double(yerr) / double(ystrips);

	int xoff = STEP / 2;
	int yoff = STEP / 2;
	//-------------------------
	numseeds = xstrips*ystrips;
	//-------------------------
	kseedsl.resize(numseeds);
	kseedsa.resize(numseeds);
	kseedsb.resize(numseeds);
	kseedsx.resize(numseeds);
	kseedsy.resize(numseeds);

	for (int y = 0; y < ystrips; y++)
	{
		int ye = y*yerrperstrip;
		for (int x = 0; x < xstrips; x++)
		{
			int xe = x*xerrperstrip;
			int seedx = (x*STEP + xoff + xe);
			if (hexgrid){ seedx = x*STEP + (xoff << (y & 0x1)) + xe; seedx = min(m_width - 1, seedx); }//for hex grid sampling
			int seedy = (y*STEP + yoff + ye);
			int i = seedy*m_width + seedx;

			kseedsl[n] = m_lvec[i];
			kseedsa[n] = m_avec[i];
			kseedsb[n] = m_bvec[i];
			kseedsx[n] = seedx;
			kseedsy[n] = seedy;
			n++;
		}
	}


	if (perturbseeds)
	{
		::PerturbSeeds(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, m_lvec, m_avec, m_bvec, m_width, m_height, edgemag);
	}
}
void DetectLabEdges(
	const double*				lvec,
	const double*				avec,
	const double*				bvec,
	const int&					width,
	const int&					height,
	vector<double>&				edges)
{
	int sz = width*height;

	edges.resize(sz, 0);
	for (int j = 1; j < height - 1; j++)
	{
		for (int k = 1; k < width - 1; k++)
		{
			int i = j*width + k;

			double dx = (lvec[i - 1] - lvec[i + 1])*(lvec[i - 1] - lvec[i + 1]) +
				(avec[i - 1] - avec[i + 1])*(avec[i - 1] - avec[i + 1]) +
				(bvec[i - 1] - bvec[i + 1])*(bvec[i - 1] - bvec[i + 1]);

			double dy = (lvec[i - width] - lvec[i + width])*(lvec[i - width] - lvec[i + width]) +
				(avec[i - width] - avec[i + width])*(avec[i - width] - avec[i + width]) +
				(bvec[i - width] - bvec[i + width])*(bvec[i - width] - bvec[i + width]);

			//edges[i] = fabs(dx) + fabs(dy);
			edges[i] = dx*dx + dy*dy;
		}
	}
}

//===========================================================================
///	PerturbSeeds
//===========================================================================
void  PerturbSeeds(
	vector<double>&				kseedsl,
	vector<double>&				kseedsa,
	vector<double>&				kseedsb,
	vector<double>&				kseedsx,
	vector<double>&				kseedsy,
	double*&             m_lvec,
	double*&             m_avec,
	double*&             m_bvec,
	int                         m_width,
	int                         m_height,
	const vector<double>&                   edges)
{
	const int dx8[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1, 1 };

	int numseeds = kseedsl.size();

	for (int n = 0; n < numseeds; n++)
	{
		int ox = kseedsx[n];//original x
		int oy = kseedsy[n];//original y
		int oind = oy*m_width + ox;

		int storeind = oind;
		for (int i = 0; i < 8; i++)
		{
			int nx = ox + dx8[i];//new x
			int ny = oy + dy8[i];//new y

			if (nx >= 0 && nx < m_width && ny >= 0 && ny < m_height)
			{
				int nind = ny*m_width + nx;
				if (edges[nind] < edges[storeind])
				{
					storeind = nind;
				}
			}
		}
		if (storeind != oind)
		{
			kseedsx[n] = storeind%m_width;
			kseedsy[n] = storeind / m_width;
			kseedsl[n] = m_lvec[storeind];
			kseedsa[n] = m_avec[storeind];
			kseedsb[n] = m_bvec[storeind];
		}
	}
}

void  PerformSuperpixelSLIC(
	vector<double>&				kseedsl,
	vector<double>&				kseedsa,
	vector<double>&				kseedsb,
	vector<double>&				kseedsx,
	vector<double>&				kseedsy,
	double*&             m_lvec,
	double*&             m_avec,
	double*&             m_bvec,
	int                         m_width,
	int                         m_height,
	int*&					klabels,
	const int&				STEP,
	const vector<double>&                   edgemag,
	const double&				M)
{
	int sz = m_width*m_height;
	const int numk = kseedsl.size();
	//----------------
	int offset = STEP;
	//if(STEP < 8) offset = STEP*1.5;//to prevent a crash due to a very small step size
	//----------------

	vector<double> clustersize(numk, 0);
	vector<double> inv(numk, 0);//to store 1/clustersize[k] values

	vector<double> sigmal(numk, 0);
	vector<double> sigmaa(numk, 0);
	vector<double> sigmab(numk, 0);
	vector<double> sigmax(numk, 0);
	vector<double> sigmay(numk, 0);
	vector<double> distvec(sz, DBL_MAX);

	double invwt = 1.0 / ((STEP / M)*(STEP / M));

	int x1, y1, x2, y2;
	double l, a, b;
	double dist;
	double distxy;
	for (int itr = 0; itr < 10; itr++)
	{
		distvec.assign(sz, DBL_MAX);
		for (int n = 0; n < numk; n++)
		{
			y1 = max(0.0, kseedsy[n] - offset);
			y2 = min((double)m_height, kseedsy[n] + offset);
			x1 = max(0.0, kseedsx[n] - offset);
			x2 = min((double)m_width, kseedsx[n] + offset);


			for (int y = y1; y < y2; y++)
			{
				for (int x = x1; x < x2; x++)
				{
					int i = y*m_width + x;

					l = m_lvec[i];
					a = m_avec[i];
					b = m_bvec[i];

					dist = (l - kseedsl[n])*(l - kseedsl[n]) +
						(a - kseedsa[n])*(a - kseedsa[n]) +
						(b - kseedsb[n])*(b - kseedsb[n]);

					distxy = (x - kseedsx[n])*(x - kseedsx[n]) +
						(y - kseedsy[n])*(y - kseedsy[n]);

					//------------------------------------------------------------------------
					dist += distxy*invwt;//dist = sqrt(dist) + sqrt(distxy*invwt);//this is more exact
					//------------------------------------------------------------------------
					if (dist < distvec[i])
					{
						distvec[i] = dist;
						klabels[i] = n;
					}
				}
			}
		}
		//-----------------------------------------------------------------
		// Recalculate the centroid and store in the seed values
		//-----------------------------------------------------------------
		//instead of reassigning memory on each iteration, just reset.

		sigmal.assign(numk, 0);
		sigmaa.assign(numk, 0);
		sigmab.assign(numk, 0);
		sigmax.assign(numk, 0);
		sigmay.assign(numk, 0);
		clustersize.assign(numk, 0);
		//------------------------------------
		//edgesum.assign(numk, 0);
		//------------------------------------

		{int ind(0);
		for (int r = 0; r < m_height; r++)
		{
			for (int c = 0; c < m_width; c++)
			{
				sigmal[klabels[ind]] += m_lvec[ind];
				sigmaa[klabels[ind]] += m_avec[ind];
				sigmab[klabels[ind]] += m_bvec[ind];
				sigmax[klabels[ind]] += c;
				sigmay[klabels[ind]] += r;
				//------------------------------------
				//edgesum[klabels[ind]] += edgemag[ind];
				//------------------------------------
				clustersize[klabels[ind]] += 1.0;
				ind++;
			}
		}}

		{for (int k = 0; k < numk; k++)
		{
			if (clustersize[k] <= 0) clustersize[k] = 1;
			inv[k] = 1.0 / clustersize[k];//computing inverse now to multiply, than divide later
		}}

		{for (int k = 0; k < numk; k++)
		{
			kseedsl[k] = sigmal[k] * inv[k];
			kseedsa[k] = sigmaa[k] * inv[k];
			kseedsb[k] = sigmab[k] * inv[k];
			kseedsx[k] = sigmax[k] * inv[k];
			kseedsy[k] = sigmay[k] * inv[k];
			//------------------------------------
			//edgesum[k] *= inv[k];
			//------------------------------------
		}}
	}
}
void  EnforceLabelConnectivity(
	const int*					labels,//input labels that need to be corrected to remove stray labels
	const int					width,
	const int					height,
	int*&						nlabels,//new labels
	int&						numlabels,//the number of labels changes in the end if segments are removed
	const int&					K) //the number of superpixels desired by the user
{
	//	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	//	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

	const int dx4[4] = { -1, 0, 1, 0 };
	const int dy4[4] = { 0, -1, 0, 1 };

	const int sz = width*height;
	const int SUPSZ = sz / K;
	//nlabels.resize(sz, -1);
	for (int i = 0; i < sz; i++) nlabels[i] = -1;
	int label(0);
	int* xvec = new int[sz];
	int* yvec = new int[sz];
	int oindex(0);
	int adjlabel(0);//adjacent label
	for (int j = 0; j < height; j++)
	{
		for (int k = 0; k < width; k++)
		{
			if (0 > nlabels[oindex])
			{
				nlabels[oindex] = label;
				//--------------------
				// Start a new segment
				//--------------------
				xvec[0] = k;
				yvec[0] = j;
				//-------------------------------------------------------
				// Quickly find an adjacent label for use later if needed
				//-------------------------------------------------------
				{for (int n = 0; n < 4; n++)
				{
					int x = xvec[0] + dx4[n];
					int y = yvec[0] + dy4[n];
					if ((x >= 0 && x < width) && (y >= 0 && y < height))
					{
						int nindex = y*width + x;
						if (nlabels[nindex] >= 0) adjlabel = nlabels[nindex];
					}
				}}

				int count(1);
				for (int c = 0; c < count; c++)
				{
					for (int n = 0; n < 4; n++)
					{
						int x = xvec[c] + dx4[n];
						int y = yvec[c] + dy4[n];

						if ((x >= 0 && x < width) && (y >= 0 && y < height))
						{
							int nindex = y*width + x;

							if (0 > nlabels[nindex] && labels[oindex] == labels[nindex])
							{
								xvec[count] = x;
								yvec[count] = y;
								nlabels[nindex] = label;
								count++;
							}
						}

					}
				}
				//-------------------------------------------------------
				// If segment size is less then a limit, assign an
				// adjacent label found before, and decrement label count.
				//-------------------------------------------------------
				if (count <= SUPSZ >> 2)
				{
					for (int c = 0; c < count; c++)
					{
						int ind = yvec[c] * width + xvec[c];
						nlabels[ind] = adjlabel;
					}
					label--;
				}
				label++;
			}
			oindex++;
		}
	}
	numlabels = label;

	if (xvec) delete[] xvec;
	if (yvec) delete[] yvec;
}
/////////////////////////////////////////////////////
//画出边界点
void  DrawContoursAroundSegments(
	unsigned char*&			ubuff,
	int*&					labels,
	const int&              lLineBytes,
	const int&				width,
	const int&				height,
	const unsigned int&		color)
{
	const int dx8[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1, 1 };
	int sz = width*height;
	vector<bool> istaken(sz, false);
	vector<int> contourx(sz); vector<int> contoury(sz);
	int mainindex(0); int cind(0);
	int j;
	for (j = 0; j < height; j++)
	{
		for (int k = 0; k < width; k++)
		{
			int np(0);
			for (int i = 0; i < 8; i++)
			{
				int x = k + dx8[i];
				int y = j + dy8[i];

				if ((x >= 0 && x < width) && (y >= 0 && y < height))
				{
					int index = y*width + x;

					//if( false == istaken[index] )//comment this to obtain internal contours
					{
						if (labels[mainindex] != labels[index]) np++;
					}
				}
			}
			if (np > 1)
			{
				contourx[cind] = k;
				contoury[cind] = j;
				istaken[mainindex] = true;
				//img[mainindex] = color;
				cind++;
			}
			mainindex++;
		}
	}

	int numboundpix = cind;//int(contourx.size());
	for (j = 0; j < numboundpix; j++)
	{
		//画出边界点
		int ii = contoury[j] * width + contourx[j];
		ubuff[ii * 3] = 255;
		ubuff[ii * 3 + 1] = 0;
		ubuff[ii * 3 + 2] = 0;


		//为了更加显示边界上的点，对边界点的八邻域点扩展
		/*for (int n = 0; n < 8; n++)
		{
			int x = contourx[j] + dx8[n];
			int y = contoury[j] + dy8[n];
			if ((x >= 0 && x < width) && (y >= 0 && y < height))
			{
				int ind = y*width + x;
				if (!istaken[ind])
				{
					ubuff[ind * 3] = 255;
					ubuff[ind * 3 + 1] = 0;
					ubuff[ind * 3 + 2] = 0;
				}
			}
		}*/
	}
}


//生成超像素标签
void performSLIC(unsigned char *ubuff, //判断数组
	int width,//图像宽度
	int height,//图像高度
	int k,//超像素的个数
	int *&klabels,//每个像素的超像素的标记
	int &numlabels//超像素的个数
	//unsigned char *&gray,//超像素数组均值
	//float*& SD//超像素的方差
	)
{
	//k = width*height / k + 1;  //循环判断时

	k = width*height / 25 + 1;
	int compactness = 10;//紧凑系数
	const int superpixelsize = 0.5 + double(width*height) / double(k);
	const int STEP = sqrt(double(superpixelsize)) + 0.5;
	vector<double>kseedsl(0);
	vector<double>kseedsa(0);
	vector<double>kseedsb(0);
	vector<double>kseedsx(0);
	vector<double>kseedsy(0);
	int sz = width*height;

	//像素的超像素的标签
	klabels = new int[sz];
	int i, j;
	for (i = 0; i<sz; i++)
		klabels[i] = -1;

	double* m_lvec;
	double* m_bvec;
	double* m_avec;

	unsigned int* ubuff1 = new unsigned int[3 * width*height];
	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			int temp = i*width + j;
			ubuff1[3 * temp] = (unsigned int)ubuff[3 * temp];
			ubuff1[3 * temp + 1] = (unsigned int)ubuff[3 * temp + 1];
			ubuff1[3 * temp + 2] = (unsigned int)ubuff[3 * temp + 2];
		}
	}
	DoRGBtoLABConversion((const unsigned int*&)ubuff1, m_lvec, m_avec, m_bvec, width, height);
	delete[] ubuff1; ubuff1 = NULL;

	bool perturbseeds(true);
	vector<double>edgemag(0);
	if (perturbseeds)
	{
		DetectLabEdges(m_lvec, m_avec, m_bvec, width, height, edgemag);
	}
	GetLABXYSeeds_ForGivenStepSize(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, width, height, m_lvec, m_avec,
		m_bvec, STEP, perturbseeds, edgemag);
	PerformSuperpixelSLIC(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, m_lvec, m_avec, m_bvec, width,
		height, klabels, STEP, edgemag, compactness);
	numlabels = kseedsl.size();//超像素的个数
	int* nlabels = new int[sz];
	EnforceLabelConnectivity(klabels, width, height, nlabels, numlabels, double(sz) / double(STEP*STEP));
	for (i = 0; i<sz; i++)
	{
		klabels[i] = nlabels[i];
	}
	if (nlabels)
	{
		delete[]nlabels;
	}

	//****************************
	//画出超像素的分割边界
	//****************************
	/*unsigned char* ubuff2 = new unsigned char[width*height * 3];
	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			int temp = i*width + j;
			ubuff2[3 * temp] = ubuff[3 * temp];
			ubuff2[3 * temp + 1] = ubuff[3 * temp + 1];
			ubuff2[3 * temp + 2] = ubuff[3 * temp + 2];
		}
	}
	DrawContoursAroundSegments(ubuff2, klabels, 1, width, height, 1);

	saveImage2(width, height, 3, GDT_Byte, ubuff2);
	delete[] ubuff2; ubuff2 = NULL;*/



	//经超像素光滑后，每个超像素中像素点的灰度数组
	//RGBtoGray(ubuff, width, height, klabels, numlabels, gray, SD);

	delete[]m_avec; m_avec = NULL;
	delete[]m_bvec; m_bvec = NULL;
	delete[]m_lvec; m_lvec = NULL;
	vector <double>().swap(kseedsl);
	vector <double>().swap(kseedsa);
	vector <double>().swap(kseedsb);
	vector <double>().swap(kseedsx);
	vector <double>().swap(kseedsy);
	vector <double>().swap(edgemag);
}


//求对生成的超像素灰度均值及方差
void  RGBtoGray(unsigned char *ubuff,//图像数组
	int width,//图像宽度
	int height,//图像高度
	int* klabels,//像素点的标签
	int  seednum,//超像素的分割数
	unsigned char *&gray,//超像素的灰度数组
	float *&SD//超像素的灰度方差
	)
{
	int i, j;

	gray = new unsigned char[seednum];
	float *R = new float[seednum];
	int *num = new int[seednum];
	SD = new float[seednum];

	for (i = 0; i < seednum; i++)
	{
		gray[i] = 0;
		R[i] = 0.0;
		num[i] = 0;
		SD[i] = 0.0;
	}

	for (i = 0; i<height; i++)
	{
		for (j = 0; j < width; j++)
		{
			int temp = i*width + j;
			R[klabels[temp]] += (float)ubuff[temp];
			num[klabels[temp]]++;
		}
	}
	for (i = 0; i<seednum; i++)
	{
		if (num[i] != 0)
		{
			R[i] /= num[i];		//平均值	
		}
	}

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			int temp = i*width + j;
			SD[klabels[temp]] += pow((ubuff[temp] - R[klabels[temp]]), 2);
		}
	}


	for (i = 0; i < seednum; i++)
	{
		//求方差
		SD[i] /= num[i];

		//求灰度均值
		if (R[i]>255)
		{
			R[i] = 255;
		}
		gray[i] = (unsigned char)R[i];
	}




	delete[]R; R = NULL;
	delete[]num; num = NULL;
}


//求每个超像素的邻接超像素
void SuperpixelAdjacent(int *klabels,//每个像素的超像素标签
	int height,//图像高度
	int width,//图像宽度
	int seedNum,//超像素的个数
	vector< vector<int> >& adjacent//超像素的邻接矩阵容器
	)
{
	int i, j;

	adjacent.resize(seedNum);
	for (i = 0; i < seedNum; i++)
	{
		adjacent[i].clear();
	}

	const int dx8[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1, 1 };

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			int temp = i*width + j;
			for (int k = 0; k < 8; k++)
			{
				int x = j + dx8[k];
				int y = i + dy8[k];
				if (x >= 0 && x < width && y >= 0 && y < height)
				{
					int temp1 = y*width + x;
					if (klabels[temp] != klabels[temp1])
					{
						//对temp标签超像素进行邻域增长
						if (adjacent[klabels[temp]].size() == 0)
						{
							adjacent[klabels[temp]].push_back(klabels[temp1]);
						}
						else
						{
							int flag1 = 0;
							for (int w1 = 0; w1 < adjacent[klabels[temp]].size(); w1++)
							{
								if (klabels[temp1] == adjacent[klabels[temp]][w1])
								{
									flag1 = 1;
									break;
								}
							}
							if (flag1 == 0)
							{
								adjacent[klabels[temp]].push_back(klabels[temp1]);
							}
						}

						//对temp1标签超像素进行邻域增长
						if (adjacent[klabels[temp1]].size() == 0)
						{
							adjacent[klabels[temp1]].push_back(klabels[temp]);
						}
						else
						{
							int flag2 = 0;
							for (int w2 = 0; w2 < adjacent[klabels[temp1]].size(); w2++)
							{
								if (klabels[temp] == adjacent[klabels[temp1]][w2])
								{
									flag2 = 1;
									break;
								}
							}
							if (flag2 == 0)
							{
								adjacent[klabels[temp1]].push_back(klabels[temp]);
							}
						}


					}
				}
			}//end for k 8邻域判断

		}//end for j
	}//end for i

}


//邻接超像素进行区域增长
void seedGrow(unsigned char *ubuff,
	int height,
	int width,
	int *klabels,
	int numlabels,
	unsigned char* gray,
	unsigned char cloudThreshold,//区域增长的云的阈值
	vector<int> seedFlag)//经过滤保留下来的超像素的标签
{
	int i, j;

	//多超像素进行区域增长
	vector< vector<int> > adj;//超像素间的邻接超像素矩阵
	SuperpixelAdjacent(klabels, height, width, numlabels, adj);

	int *seedGrowFlag = new int[numlabels];//对超像素是否需要区域增长进行标记
	for (i = 0; i < numlabels; i++)
	{
		seedGrowFlag[i] = 0;
	}
	for (i = 0; i < seedFlag.size(); i++)
	{
		for (j = 0; j < adj[seedFlag[i]].size(); j++)
		{
			if (gray[adj[seedFlag[i]][j]] >= cloudThreshold)
			{
				seedGrowFlag[adj[seedFlag[i]][j]] = 1;
			}
		}
	}

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			int temp = i*width + j;
			if (seedGrowFlag[klabels[temp]] == 1)
			{
				ubuff[temp] = 255;
			}
		}
	}

	delete[]seedGrowFlag; seedGrowFlag = NULL;
	for (i = 0; i < numlabels; i++)
	{
		vector<int>().swap(adj[i]);
	}
}