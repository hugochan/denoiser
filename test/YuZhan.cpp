#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <fstream>

// add
#include <assert.h>
#include <math.h>
#include <time.h>

using namespace std;

int process_count = 0;

int test = 1;
#define DIM 3
#define BIG  1000000000
#define HISTOGNUM  30
#define ASC_LINE_LEN 14
#define ERROR 0
#define TRUE  1
#define FALSE 0
#define M_PI  3.1415926
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l] = h + s*(g - h*tau);
static int iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
	(iminarg1) : (iminarg2))
#define IMAX(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
	(iminarg2) : (iminarg1))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define ABS(a) ((a) >= 0.0 ? (a) : -(a))

#define  NR_END  1
#define  FREE_ARG  char*
#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;}
typedef float REAL;

typedef float REAL_DIM[3];
typedef float REAL;
typedef int INT32;
typedef short INT16;
typedef INT16 INT16_DIM[DIM];
typedef REAL REAL_DIM[DIM];
typedef REAL REAL_33[3][3];

struct VoxelDataStruc
{
	int count;
	int connectivityFlag_backup;
	int connectivityFlag;   //  false -- isolated,  true -- connected
	// -1 ---  untouched  0 -- unprocessed,  1 -- processed   2 -- true data  3 -- outlier 
	float maxError;
	// local coord directions
	REAL_DIM n;  // z direction  the third eigen vector

	char hlevel;  // -1 -- default,  0 -- 100 histogram bins

	REAL centroid[DIM];
};



static REAL maxarg1;
static REAL maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))
#define FMIN(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg2) : (maxarg1))

// declare
REAL angle_vect(REAL_DIM a, REAL_DIM b);
void propagateVoxel(int i0, int j0, int k0, VoxelDataStruc ***vds, int cubelen, int **** list);
void calcCentroid(int num_nlist, INT32 *nlist, REAL_DIM* vert_v, REAL centroid[]);
void tellTrueFalse(int num_nlist, INT32 *nlist, REAL_DIM* vert_v, bool* istatus);
void calcEigenVector(int num_nlist, INT32 *nlist, REAL_DIM vtmp2, REAL_DIM* vert_v);
void calcPointNoise(int num_nlist, INT32 *nlist, REAL_DIM vtmp2, REAL noise[], REAL_DIM* vert_v);
void jacobi(REAL **a, int n, REAL d[], REAL **v, int *nrot);
void eigsrt(REAL d[], REAL **v, int n);
REAL *dvector(int nl, int nh);
void free_dvector(REAL *v, int nl, int nh);
void dvector_cross_product(REAL *A, REAL *B, REAL *C, INT16 n, INT16 m);
void dvector_cross_product(REAL *A, REAL *B, REAL *C, INT16 n, INT16 m);
void dvector_assign(REAL_DIM source, REAL_DIM target, int dim);
void calc_Cartesian_frame_to_frame_transformation_matrix(REAL_DIM U, REAL_DIM V, REAL_DIM W,
	REAL_DIM O, REAL E[][4]);
REAL **dmatrix(int nrl, int nrh, int ncl, int nch);
REAL dvector_dot_product(REAL *a, REAL *b, int n, int m);
void calc_frame_to_Cartesian_frame_transformation_matrix(REAL_DIM U, REAL_DIM V, REAL_DIM W,
	REAL_DIM O, REAL E[][4]);
void svdcmp(REAL **a, int m, int n, REAL w[], REAL **v);
void svbksb(REAL **u, REAL w[], REAL **v, int m, int n, REAL b[], REAL x[]);
void free_dmatrix(REAL **m, int nrl, int nrh, int ncl, int nch);
REAL pythag(REAL a, REAL b);
void dvector_negation(REAL_DIM x, int n, int m);
REAL NormalizeEigenValue(REAL_DIM vtmp2);
void calcEigenFactor(int num_list, INT32* nlist, REAL_DIM* vert_v, REAL* EigenFactor);



//typedef int INTN[NUM_DENOISING_NEIGHBOR];
// global variables
REAL *D = dvector(1, 3);
REAL **V2 = dmatrix(1, 6, 1, 6);
REAL **V3 = dmatrix(1, 3, 1, 3);
REAL **A = dmatrix(1, 3, 1, 3);
REAL *W = dvector(1, 6);
REAL *X = dvector(1, 6);
REAL E[4][4], S[4], F[4][4];
REAL_DIM V[3];
int curvert;
REAL min_cubesize = BIG;

struct CoordDataStruc
{
	char i;
	char j;
	char k;
};



typedef int INT32;
//




int main()
{
	time_t t_start, t_end;
	t_start = time(NULL);
	ifstream scr0;
	int begin, end, size;
	REAL_DIM *vert_v;

	//add
	//char input_filename[50] = "NoiseModel1.asc";
	char input_filename[50] = "NoiseModel2.asc";

	//

	scr0.open(input_filename, ios::in);   // .asc file
	if (!scr0)
	{
		cout << "Open file error with " << input_filename << endl;
		return -1;
	}

	begin = scr0.tellg();
	scr0.seekg(0, ios::end);
	end = scr0.tellg();

	// add
	int VERSION = 1;
	//

	if (end > 4000000000)
	{
		scr0.close(); scr0.clear();
		cout << "Error: Version " << VERSION << " cannot handle data models with more than 50 million points !"
			<< endl;
		return -1;
	}


	size = end - begin;

	if (size <= 0)
	{
		scr0.close(); scr0.clear();
		return 0;
	}


	char *str, *str_sav;
	int mem_loop = 0;

	// add 
	char* memblock;
	//

	memblock = new char[size];

	scr0.seekg(0, ios::beg);
	scr0.read(memblock, size);
	scr0.close(); scr0.clear();

	// find the number of vertices
	//char *str, *str_sav;
	str = memblock;

	// add
	int k;
	//

	k = 0;
	while (str = strstr(str, "\n"))
	{
		str++;
		if (str[0] == '\r')
			str++;
		k++;
	}

	vert_v = new REAL_DIM[k];

	// add
	int i, j;
	//

	i = 0;

	str = memblock;
	str_sav = str;
	// add
	int len0;
	char* token;
	char* current_line = new char[size];
	char * seps = " ";
	//

	while (str = strstr(str, "\n"))
	{

		len0 = str - str_sav;

		//if(str_sav[0] != '#' && len0 >= 24)
		if (str_sav[0] != '#' && len0 >= ASC_LINE_LEN)
		{
			for (j = 0; j < len0; j++)
				current_line[j] = str_sav[j];
			current_line[j] = '\0';

			token = strtok(current_line, seps); // Ignore the card identifier

			//x[i]._v[0] = (REAL) atof(token);
			vert_v[i][0] = (REAL)atof(token);

			token = strtok(NULL, seps); // Ignore the card identifier

			//x[i]._v[1] = (REAL) atof(token);
			vert_v[i][1] = (REAL)atof(token);

			token = strtok(NULL, seps); // Ignore the card identifier

			//x[i]._v[2] = (REAL) atof(token);
			vert_v[i][2] = (REAL)atof(token);

			i++;
		} // end of if(str_sav[0] != '#')

		str++;
		if (str[0] == '\r')
			str++;
		str_sav = str;
	} // end of while(str=strstr(str, "\n"))
	delete[] current_line;


	curvert = i;
	delete[] memblock;


	// create rectangular grids
	float xmin[DIM], xmax[DIM];
	for (j = 0; j < 3; j++)
	{
		xmin[j] = BIG;
		xmax[j] = -BIG;
	}

	for (i = 0; i < curvert; i++)
	for (j = 0; j<3; j++)
	{
		if (xmin[j]>vert_v[i][j]) xmin[j] = vert_v[i][j];
		if (xmax[j] < vert_v[i][j]) xmax[j] = vert_v[i][j];
	}

	int m;

	//m = (int) pow(10, (1/3.0)*log10(double(curvert)));

	//if(m > 200)  
	//m = 100;
	m = 50;

	for (j = 0; j < DIM; j++)
	{
		if (min_cubesize > (xmax[j] - xmin[j]) / m)
			min_cubesize = (xmax[j] - xmin[j]) / m;
	}


	int ***count = new int**[m];
	int ****list = new int***[m];
	VoxelDataStruc ***vds = new VoxelDataStruc**[m];  // Create volxel data structure
	for (i = 0; i < m; i++)
	{
		try
		{
			count[i] = new int*[m];
			list[i] = new int**[m];
			vds[i] = new VoxelDataStruc*[m];
		}
		catch (std::bad_alloc&)
		{

		}
		for (j = 0; j < m; j++)
		{
			try
			{
				count[i][j] = new int[m];
				list[i][j] = new int*[m];
				vds[i][j] = new VoxelDataStruc[m];
			}
			catch (std::bad_alloc&)
			{

			}
			for (k = 0; k < m; k++)
			{
				try
				{
					count[i][j][k] = 0;
					list[i][j][k] = NULL;
					vds[i][j][k].connectivityFlag = false;
				}
				catch (std::bad_alloc&)
				{

				}
			}
		}
	}

	int index[DIM];
	for (i = 0; i < curvert; i++)
	{
		for (j = 0; j < 3; j++)
		{
			index[j] = floor((vert_v[i][j] - xmin[j]) / (xmax[j] - xmin[j])*(m - 1) + 0.5);
			if (index[j] >= m)
				assert(index[j] < m);
		}
		count[index[0]][index[1]][index[2]]++;
	}

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			for (k = 0; k < m; k++)
			{
				list[i][j][k] = new int[count[i][j][k]];
				count[i][j][k] = 0;
			}
		}
	}

	for (j = 0; j<3; j++)
		assert((xmax[j] - xmin[j])> 0.00001);

	for (i = 0; i < curvert; i++)
	{
		for (j = 0; j < 3; j++)
		{
			index[j] = floor((vert_v[i][j] - xmin[j]) / (xmax[j] - xmin[j])*(m - 1) + 0.5);
		}
		list[index[0]][index[1]][index[2]][count[index[0]][index[1]][index[2]]] = i;
		count[index[0]][index[1]][index[2]]++;
	}


	ofstream tmpoup1, tmpoup2, tmpoup;
	// debug output
	tmpoup1.open("tmp1.asc", ios::out);
	tmpoup2.open("tmp2.asc", ios::out);


	// add
	float* vtmp = new float[DIM];
	float* GmaxerrorArray = new float[curvert];
	float* EigenFactor = new float;
	//
	// denoising
	int itmp = 0, ithreshold = 1, l, id, local_sum, ihistogram = 15;
	int *icount = new int[curvert];
	bool *istatus = new bool[curvert];
	REAL noise[2000];
	int choice = 8;  //  1 -- eigen value    2 --  least-square fitting error  3 -- density 4 -- eigen value plus
	//  6 -- non-connectivity factor (outlier should be greater) 7 -- conn + telltruefalse
	// 8 -- connectivity + Bi-means Cluster + surface propogation
	int ic = 0;
	int inum = 0;
	int istep = 3;
	if (choice == 6 || choice == 7 || choice == 8)
	{
		ithreshold = 2;
		// setup voxel information
		for (i = 0; i < m; i++)
		{
			for (j = 0; j < m; j++)
			{
				for (k = 0; k<m; k++)
				{
					vds[i][j][k].count = count[i][j][k]; //  count
					if (vds[i][j][k].connectivityFlag == false)
					{

						if (count[i][j][k] > ithreshold)
						{

							vds[i][j][k].connectivityFlag = true;
							if (i + 1 < m)  vds[i + 1][j][k].connectivityFlag = true;
							if (j + 1 < m)  vds[i][j + 1][k].connectivityFlag = true;
							if (k + 1 < m)  vds[i][j][k + 1].connectivityFlag = true;
							if (i - 1 >= 0)  vds[i - 1][j][k].connectivityFlag = true;
							if (j - 1 >= 0)  vds[i][j - 1][k].connectivityFlag = true;
							if (k - 1 >= 0)  vds[i][j][k - 1].connectivityFlag = true;
						}
						else
						{

							if ((i + 1 < m) && count[i + 1][j][k] > ithreshold)
							{
								vds[i][j][k].connectivityFlag = true; continue;
							}
							if ((j + 1 < m) && count[i][j + 1][k] > ithreshold)
							{
								vds[i][j][k].connectivityFlag = true; continue;
							}

							if ((k + 1 < m) && count[i][j][k + 1] > ithreshold)
							{
								vds[i][j][k].connectivityFlag = true; continue;
							}
							if ((i - 1 >= 0) && count[i - 1][j][k] > ithreshold)
							{
								vds[i][j][k].connectivityFlag = true; continue;
							}
							if ((j - 1 >= 0) && count[i][j - 1][k] > ithreshold)
							{
								vds[i][j][k].connectivityFlag = true; continue;
							}
							if ((k - 1 >= 0) && count[i][j][k - 1] > ithreshold)
							{
								vds[i][j][k].connectivityFlag = true; continue;
							}
						}
					}
				}
			}
		}
	}

	/*
	for (i = 0; i < curvert; i++)
		icount[i] = 0;
	*/
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			for (k = 0; k<m; k++)
			{
				itmp = count[i][j][k];
				if (vds[i][j][k].connectivityFlag == true
					&& itmp > 0)// connectivity check is really useful to skip very smooth portions, and catch noisy part.

				{
					int* localIdx = new int[itmp];
					for (l = 0; l < itmp; l++)
						localIdx[l] = list[i][j][k][l];

					switch (choice)
					{
					case 1:
						//calcEigenVector(itmp, localIdx[0], vtmp);
						calcEigenVector(itmp, localIdx, vtmp, vert_v);

						break;
					case 2:
						calcPointNoise(itmp, localIdx, vtmp, noise, vert_v);
						break;
					case 4:
						calcEigenVector(itmp, localIdx, vtmp, vert_v);
						break;
					case 5:
						calcEigenFactor(itmp, localIdx, vert_v, EigenFactor);
						break;
					case 6:
						calcEigenVector(itmp, localIdx, vtmp, vert_v);
						break;
					case 8:
						calcCentroid(itmp, localIdx, vert_v, vds[i][j][k].centroid);
						calcEigenVector(itmp, localIdx, vtmp, vert_v);
						vds[i][j][k].n[0] = V3[3][1];
						vds[i][j][k].n[1] = V3[3][2];
						vds[i][j][k].n[2] = V3[3][3];
						break;
					case 7:

						if ((i + 1 < m) && count[i + 1][j][k] >= ithreshold)
						{
							inum += 1;
						}
						if ((j + 1 < m) && count[i][j + 1][k] >= ithreshold)
						{
							inum += 1;
						}

						if ((k + 1 < m) && count[i][j][k + 1] >= ithreshold)
						{
							inum += 1;
						}
						if ((i - 1 >= 0) && count[i - 1][j][k] >= ithreshold)
						{
							inum += 1;
						}
						if ((j - 1 >= 0) && count[i][j - 1][k] >= ithreshold)
						{							
							inum += 1;
						}
						if ((k - 1 >= 0) && count[i][j][k - 1] >= ithreshold)
						{
							inum += 1;
						}

						if (inum <= 4)
						{
							for (l = 0; l < itmp; l++)
							{
								id = localIdx[l];
								istatus[id] = true;
							}
						}
						else
						{

							tellTrueFalse(itmp, localIdx, vert_v, istatus);
						}
						/*
						if (i == 0 || i == m - 1 || j == 0 || j == m - 1 || k == 0 || k == m - 1) // margin points
						{
							tellTrueFalse(itmp, localIdx, vert_v, istatus);
						}
						else
						{
							for (l = 0; l < itmp; l++)
							{
								id = localIdx[l];
								istatus[id] = true;
							}
						}
						*/
						break;
					}
					/*
					if (choice == 6)
					{
					ic = 0;
					if (i<m - istep && count[(i + istep) % m][j][k] == 0)   ic++;
					if (j<m - istep && count[i][(j + istep) % m][k] == 0)   ic++;
					if (k<m - istep && count[i][j][(k + istep) % m] == 0)   ic++;
					if (i>2 && count[(i - istep) % m][j][k] == 0)   ic++;
					if (j>2 && count[i][(j + m - istep) % m][k] == 0)   ic++;
					if (k>2 && count[i][j][(k + m - istep) % m] == 0)   ic++;
					}
					*/
					if (choice != 7)
					{
						for (l = 0; l < itmp; l++)
						{
							id = localIdx[l];

							switch (choice)
							{
							case 1:
								GmaxerrorArray[id] = vtmp[2];
								break;
							case 2:
								GmaxerrorArray[id] = noise[l];
								break;
							case 3:
								GmaxerrorArray[id] = count[i][j][k];
								break;
							case 4:
								GmaxerrorArray[id] = NormalizeEigenValue(vtmp);
							case 5:
								GmaxerrorArray[id] = *EigenFactor;
							case 6: // non-connectivity factor
								GmaxerrorArray[id] = vtmp[2];
								break;
							case 8:
								GmaxerrorArray[id] = vtmp[2];
							}

							istatus[id] = true;
						}
					}
					delete[]localIdx;
				}
				else  // small point cloud
				{
					for (l = 0; l < itmp; l++)
					{
						id = list[i][j][k][l];
						istatus[id] = false;
						GmaxerrorArray[id] = 0.0;
					}
				}
			}
		}
	}

	if (choice != 7 && choice != 8)
	{
		// Bi-means clustering
		double vmax = -BIG, vmin = BIG, va = 0.0;
		int num_point = 0;
		for (i = 0; i<curvert; i++)
		{
			if (istatus[i] == true)
			{
				if (GmaxerrorArray[i]> vmax) vmax = GmaxerrorArray[i];
				if (GmaxerrorArray[i] < vmin) vmin = GmaxerrorArray[i];
				num_point++;
			}
		}

		ihistogram = 15;
		double vsize = vmax - vmin;
		for (i = 0; i < curvert; i++)
		{
			if (istatus[i] == true)
			{
				j = int((GmaxerrorArray[i] - vmin)*float(HISTOGNUM - 1) / vsize);

				if (j >= ihistogram)
				{
					istatus[i] = false;
					tmpoup2 << vert_v[i][0] << " " << vert_v[i][1] << " " << vert_v[i][2] << endl;
				}
				else
				{
					tmpoup1 << vert_v[i][0] << " " << vert_v[i][1] << " " << vert_v[i][2] << endl;
				}
			}
			else
			{
				tmpoup2 << vert_v[i][0] << " " << vert_v[i][1] << " " << vert_v[i][2] << endl;

			}
		}

	}
	else if (choice == 7)
	{
		for (i = 0; i < curvert; i++)
		{
			if (istatus[i] == true)
			{

				tmpoup1 << vert_v[i][0] << " " << vert_v[i][1] << " " << vert_v[i][2] << endl;
			}
			else
			{
				tmpoup2 << vert_v[i][0] << " " << vert_v[i][1] << " " << vert_v[i][2] << endl;
			}
		}
	}
	else if (choice == 8)
	{
		// Bi-means clustering
		double vmax = -BIG, vmin = BIG;
//		int num_point = 0;
		for (i = 0; i<curvert; i++)
		{
			if (istatus[i] == true)
			{
				if (GmaxerrorArray[i]> vmax) vmax = GmaxerrorArray[i];
				if (GmaxerrorArray[i] < vmin) vmin = GmaxerrorArray[i];
//				num_point++;
			}
		}

		double vsize = vmax - vmin;

		// back up connectivityFlag
		for (i = 0; i < m; i++)
		{
			for (j = 0; j < m; j++)
			{
				for (k = 0; k < m; k++)
				{
					vds[i][j][k].connectivityFlag_backup = vds[i][j][k].connectivityFlag;
				}
			}
		}


		// *** Voxel propagation  ***
		// Initialization
		int i2;
		for (i = 0; i<m; i++)
		{
			for (j = 0; j<m; j++)
			{
				for (k = 0; k<m; k++)
				{
					itmp = count[i][j][k];

					if (itmp == 0 || vds[i][j][k].connectivityFlag == false)
					{
						vds[i][j][k].hlevel = HISTOGNUM;  // for propagation purpose
					}
					else
					{
						for (l = 0; l<itmp; l++)
						{
							id = list[i][j][k][l];
							if (istatus[id] == true)
							{
								i2 = int((GmaxerrorArray[id] - vmin)*float(HISTOGNUM - 1) / vsize);
								if (i2 >= HISTOGNUM) i2 = HISTOGNUM - 1;
								vds[i][j][k].hlevel = i2;
							}
						}
					}
					vds[i][j][k].connectivityFlag = -1; // for propagation purpose
				}
			}
		}


		// Propagation
		int count_untouched = 0;

		i2 = 2;
		//i2 = HISTOGNUM-1;
		int break_flag = false;

		for (l = 1; l<i2; l++)
		{
			for (i = 0; i<m; i++)
			{
				for (j = 0; j<m; j++)
				{
					for (k = 0; k<m; k++)
					{
						if (vds[i][j][k].hlevel == l)  // cubes with no or very few points are set to have a level = HISTOGNUM
						{
							if (vds[i][j][k].connectivityFlag == -1) // untouched
							{
								propagateVoxel(i, j, k, vds, m, list);
								itmp = count[i][j][k];
								for (l = 0; l < itmp; l++)
								{
									id = list[i][j][k][l];
								}
								break_flag = true; // single initial point
							}
						}
						if (break_flag) break;
					}
					if (break_flag) break;
				}
				if (break_flag) break;
			}
			if (break_flag) break;
		}
		cout << "propogation process count: ";
		cout << process_count << endl;
		
		for (i = 0; i<m; i++)
		{
			for (j = 0; j<m; j++)
			{
				for (k = 0; k<m; k++)
				{
					itmp = count[i][j][k];
					if (vds[i][j][k].connectivityFlag == 3 \
						|| vds[i][j][k].connectivityFlag == -1) // outlier
					{
						for (l = 0; l<itmp; l++)
						{
							id = list[i][j][k][l];
							tmpoup2 << vert_v[id][0] << " " << vert_v[id][1] << " " << vert_v[id][2] << endl;
						}
					}

					else  // true data + untouched
					{
						for (l = 0; l<itmp; l++)
						{
							id = list[i][j][k][l];
							tmpoup1 << vert_v[id][0] << " " << vert_v[id][1] << " " << vert_v[id][2] << endl;
						}
					}
				}
			}
		}
		cout << "count_untouched: ";
		cout << count_untouched << endl;

	}
	tmpoup1.close();
	tmpoup2.close();
	t_end = time(NULL);
	cout << "run time:";
	cout << t_end-t_start << endl;

	// deallocation
	delete []vtmp;
	delete [] GmaxerrorArray;
	delete EigenFactor;
	delete [] istatus;
//	delete [] icount;


	for (i = 0; i<m; i++)
	{

		for (j = 0; j < m; j++)
		{
			
			for (k = 0; k<m; k++)
			{
				delete[] list[i][j][k];
			}
			delete [] count[i][j];
			delete [] vds[i][j];
			delete [] list[i][j];
		}
		delete [] count[i];
		delete [] vds[i];
		delete[] list[i];
	}
	delete [] count;
	delete [] vds;
	delete [] list;
	system("pause");
	return 0;
}


void calcEigenFactor(int num_nlist, INT32 *nlist, REAL_DIM* vert_v, REAL* EigenFactor)
{
	int j, k, l, id2, nrot;
	REAL_DIM centroid, vtmp;
	float* vtmp2 = new float[3];

	// determine local coordinate of neighboring points
	REAL sum;

	// determine centroid
	for (k = 0; k<DIM; k++) centroid[k] = 0.0;
	for (j = 0; j<num_nlist; j++)
	for (k = 0; k<DIM; k++)
		centroid[k] += vert_v[nlist[j]][k];
	for (k = 0; k<DIM; k++) centroid[k] /= float(num_nlist);

	//begin build covariance matrix
	//first diagonal element
	sum = 0;
	for (j = 0; j<num_nlist; j++)
	{
		//sum+=(vert[nlist[j]]._v[0]-centroid[0])*(vert[nlist[j]]._v[0]-centroid[0]);
		sum += (vert_v[nlist[j]][0] - centroid[0])*(vert_v[nlist[j]][0] - centroid[0]);
	}


	A[1][1] = sum / num_nlist;

	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[0]-centroid[0])*(vert[nlist[j]]._v[1]-centroid[1]);
		sum += (vert_v[nlist[j]][0] - centroid[0])*(vert_v[nlist[j]][1] - centroid[1]);

	A[1][2] = sum / num_nlist;
	A[2][1] = A[1][2];


	//second diagonal element
	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[1]-centroid[1])*(vert[nlist[j]]._v[1]-centroid[1]);
		sum += (vert_v[nlist[j]][1] - centroid[1])*(vert_v[nlist[j]][1] - centroid[1]);

	A[2][2] = sum / num_nlist;

	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[0]-centroid[0])*(vert[nlist[j]]._v[2]-centroid[2]);
		sum += (vert_v[nlist[j]][0] - centroid[0])*(vert_v[nlist[j]][2] - centroid[2]);

	A[1][3] = sum / num_nlist;
	A[3][1] = A[1][3];

	//third diagonal element
	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[2]-centroid[2])*(vert[nlist[j]]._v[2]-centroid[2]);
		sum += (vert_v[nlist[j]][2] - centroid[2])*(vert_v[nlist[j]][2] - centroid[2]);

	A[3][3] = sum / num_nlist;

	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[1]-centroid[1])*(vert[nlist[j]]._v[2]-centroid[2]);
		sum += (vert_v[nlist[j]][1] - centroid[1])*(vert_v[nlist[j]][2] - centroid[2]);

	A[2][3] = sum / num_nlist;
	A[3][2] = A[2][3];

	/* d[1..3] returns the eigenvalues of a. v[1..3][1..3] is
	a matrix whose columns contain, on output, the normalized
	eigenvectors of a.
	*/


	jacobi(A, 3, D, V3, &nrot);

	// abs
	for (j = 0; j < DIM; j++)
	{
		//D[j + 1] = fabs(D[j+1]);
	}

	/* output of eigenvalues in descending order */
	eigsrt(D, V3, 3);

	for (j = 0; j<DIM; j++)
	{
		vtmp2[j] = D[j + 1];
	}

	// test
	REAL *vert_vector = dvector(1, 3);

	REAL error = 0;
	REAL max_error = 0;
	REAL tmp_error = 0;

	for (j = 0; j < num_nlist; j++)
	{
		for (k = 0; k < DIM; k++)
		{
			vert_vector[k + 1] = vert_v[nlist[j]][k] - centroid[k];
		}
		tmp_error = fabs(dvector_dot_product(V3[3], vert_vector, 3, 1) / sqrt(V3[3][1] * V3[3][1] \
			+ V3[3][2] * V3[3][2] + V3[3][3] * V3[3][3]));
		error += tmp_error;
		if (tmp_error > max_error) max_error = tmp_error;
	}
	error /= num_nlist;
	float ratio = 0.5;
	*EigenFactor = ratio*fabs(vtmp2[2] / vtmp2[0]) + (1-ratio)*(error/max_error);
	delete [] vtmp2;
}

void tellTrueFalse(int num_nlist, INT32 *nlist, REAL_DIM* vert_v, bool* istatus)
{
	int j, k, l, id2, nrot;
	REAL_DIM centroid, vtmp;
	float* vtmp2 = new float[3];

	// determine local coordinate of neighboring points
	REAL sum;

	// determine centroid
	for (k = 0; k<DIM; k++) centroid[k] = 0.0;
	for (j = 0; j<num_nlist; j++)
	for (k = 0; k<DIM; k++)
		centroid[k] += vert_v[nlist[j]][k];
	for (k = 0; k<DIM; k++) centroid[k] /= float(num_nlist);

	//begin build covariance matrix
	//first diagonal element
	sum = 0;
	for (j = 0; j<num_nlist; j++)
	{
		//sum+=(vert[nlist[j]]._v[0]-centroid[0])*(vert[nlist[j]]._v[0]-centroid[0]);
		sum += (vert_v[nlist[j]][0] - centroid[0])*(vert_v[nlist[j]][0] - centroid[0]);
	}


	A[1][1] = sum / num_nlist;

	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[0]-centroid[0])*(vert[nlist[j]]._v[1]-centroid[1]);
		sum += (vert_v[nlist[j]][0] - centroid[0])*(vert_v[nlist[j]][1] - centroid[1]);

	A[1][2] = sum / num_nlist;
	A[2][1] = A[1][2];


	//second diagonal element
	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[1]-centroid[1])*(vert[nlist[j]]._v[1]-centroid[1]);
		sum += (vert_v[nlist[j]][1] - centroid[1])*(vert_v[nlist[j]][1] - centroid[1]);

	A[2][2] = sum / num_nlist;

	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[0]-centroid[0])*(vert[nlist[j]]._v[2]-centroid[2]);
		sum += (vert_v[nlist[j]][0] - centroid[0])*(vert_v[nlist[j]][2] - centroid[2]);

	A[1][3] = sum / num_nlist;
	A[3][1] = A[1][3];

	//third diagonal element
	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[2]-centroid[2])*(vert[nlist[j]]._v[2]-centroid[2]);
		sum += (vert_v[nlist[j]][2] - centroid[2])*(vert_v[nlist[j]][2] - centroid[2]);

	A[3][3] = sum / num_nlist;

	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[1]-centroid[1])*(vert[nlist[j]]._v[2]-centroid[2]);
		sum += (vert_v[nlist[j]][1] - centroid[1])*(vert_v[nlist[j]][2] - centroid[2]);

	A[2][3] = sum / num_nlist;
	A[3][2] = A[2][3];

	/* d[1..3] returns the eigenvalues of a. v[1..3][1..3] is
	a matrix whose columns contain, on output, the normalized
	eigenvectors of a.
	*/


	jacobi(A, 3, D, V3, &nrot);



	/* output of eigenvalues in descending order */
	eigsrt(D, V3, 3);

	for (j = 0; j<DIM; j++)
	{
		vtmp2[j] = D[j + 1];
	}

	// tell true or false
	REAL* distance_array = new REAL[num_nlist];
	REAL DD = -(V3[3][1] * centroid[0] + V3[3][2] * centroid[1] + V3[3][3] * centroid[2]);
	REAL rootsquare_ABC = sqrt(V3[3][1] * V3[3][1] + V3[3][2] * V3[3][2] + V3[3][3] * V3[3][3]);
	double vmax = -BIG, vmin = BIG, va = 0.0;
	int ihistogram = 9999;
	int numOfHis = 10000;
	int i_bin;

	for (j = 0; j < num_nlist; j++)
	{
		distance_array[j] = abs(V3[3][1] * vert_v[nlist[j]][0] + V3[3][2] * vert_v[nlist[j]][1] + V3[3][3] * vert_v[nlist[j]][2] + DD)/rootsquare_ABC;
		if (distance_array[j]> vmax) vmax = distance_array[j];
		if (distance_array[j]< vmin) vmin = distance_array[j];
	}

	double VRatio = (vmax - vmin == 0) ? 0 : float(numOfHis - 1) / (vmax - vmin);
	for (j = 0; j < num_nlist; j++)
	{
		i_bin = int((distance_array[j] - vmin)*VRatio);

		if (i_bin >= ihistogram)
		{
			istatus[nlist[j]] = false;
		}
		else
		{
			istatus[nlist[j]] = true;
		}
	}
	delete [] vtmp2;
	delete [] distance_array;
}

void calcEigenVector(int num_nlist, INT32 *nlist, REAL_DIM vtmp2, REAL_DIM* vert_v)
{
	int j, k, l, id2, nrot;
	REAL_DIM centroid, vtmp;

	// determine local coordinate of neighboring points
	REAL sum;

	// determine centroid
	for (k = 0; k<DIM; k++) centroid[k] = 0.0;
	for (j = 0; j<num_nlist; j++)
	for (k = 0; k<DIM; k++)
		centroid[k] += vert_v[nlist[j]][k];
	for (k = 0; k<DIM; k++) centroid[k] /= float(num_nlist);

	//begin build covariance matrix
	//first diagonal element
	sum = 0;
	for (j = 0; j<num_nlist; j++)
	{
		//sum+=(vert[nlist[j]]._v[0]-centroid[0])*(vert[nlist[j]]._v[0]-centroid[0]);
		sum += (vert_v[nlist[j]][0] - centroid[0])*(vert_v[nlist[j]][0] - centroid[0]);
	}


	A[1][1] = sum / num_nlist;

	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[0]-centroid[0])*(vert[nlist[j]]._v[1]-centroid[1]);
		sum += (vert_v[nlist[j]][0] - centroid[0])*(vert_v[nlist[j]][1] - centroid[1]);

	A[1][2] = sum / num_nlist;
	A[2][1] = A[1][2];


	//second diagonal element
	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[1]-centroid[1])*(vert[nlist[j]]._v[1]-centroid[1]);
		sum += (vert_v[nlist[j]][1] - centroid[1])*(vert_v[nlist[j]][1] - centroid[1]);

	A[2][2] = sum / num_nlist;

	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[0]-centroid[0])*(vert[nlist[j]]._v[2]-centroid[2]);
		sum += (vert_v[nlist[j]][0] - centroid[0])*(vert_v[nlist[j]][2] - centroid[2]);

	A[1][3] = sum / num_nlist;
	A[3][1] = A[1][3];

	//third diagonal element
	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[2]-centroid[2])*(vert[nlist[j]]._v[2]-centroid[2]);
		sum += (vert_v[nlist[j]][2] - centroid[2])*(vert_v[nlist[j]][2] - centroid[2]);

	A[3][3] = sum / num_nlist;

	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[1]-centroid[1])*(vert[nlist[j]]._v[2]-centroid[2]);
		sum += (vert_v[nlist[j]][1] - centroid[1])*(vert_v[nlist[j]][2] - centroid[2]);

	A[2][3] = sum / num_nlist;
	A[3][2] = A[2][3];

	/* d[1..3] returns the eigenvalues of a. v[1..3][1..3] is
	a matrix whose columns contain, on output, the normalized
	eigenvectors of a.
	*/


	jacobi(A, 3, D, V3, &nrot);

	// abs
	for (j = 0; j < DIM; j++)
	{
		//D[j + 1] = fabs(D[j+1]);
	}

	/* output of eigenvalues in descending order */
	eigsrt(D, V3, 3);

	for (j = 0; j<DIM; j++)
	{
		vtmp2[j] = D[j + 1];
	}
}

void jacobi(REAL **a, int n, REAL d[], REAL **v, int *nrot)
{
	/* Computer all eigenvalues and eigenvectors of a real symmetric matrix
	a[1..n][1..n]. On output, elements of a above the diagonal are
	destroyed. d[1..n] returns the eigenvalues of a. v[1..n][1..n] is
	a matrix whose columns contain, on output, the normalized eigenvectors
	of a. nrot returns the number of Jacobi rotations that were requried.
	*/

	int j, iq, ip, i;
	REAL tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

	b = dvector(1, n);
	z = dvector(1, n);
	for (ip = 1; ip <= n; ip++)
	{
		for (iq = 1; iq <= n; iq++) v[ip][iq] = 0.0;
		v[ip][ip] = 1.0;
	}
	for (ip = 1; ip <= n; ip++) {
		b[ip] = d[ip] = a[ip][ip];
		z[ip] = 0.0;
	}

	*nrot = 0;

	for (i = 1; i <= 50; i++) {
		sm = 0.0;
		for (ip = 1; ip <= n - 1; ip++) {
			for (iq = ip + 1; iq <= n; iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free_dvector(z, 1, n);
			free_dvector(b, 1, n);
			return;
		}
		if (i<4)
			tresh = 0.2*sm / (n*n);
		else
			tresh = 0.0;
		for (ip = 1; ip <= n - 1; ip++) {
			for (iq = ip + 1; iq <= n; iq++) {
				g = 100.0*fabs(a[ip][iq]);
				if (i>4 && (REAL)(fabs(d[ip]) + g) == (REAL)fabs(d[ip])
					&& (REAL)(fabs(d[iq]) + g) == (REAL)fabs(d[iq]))
					a[ip][iq] = 0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h = d[iq] - d[ip];
					if ((REAL)(fabs(h) + g) == (REAL)fabs(h))
						t = (a[ip][iq]) / h;
					else {
						theta = 0.5*h / (a[ip][iq]);
						t = 1.0 / (fabs(theta) + sqrt(1.0 + theta*theta));
						if (theta<0.0) t = -t;
					}
					c = 1.0 / sqrt(1 + t*t);
					s = t*c;
					tau = s / (1.0 + c);
					h = t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq] = 0.0;
					for (j = 1; j <= ip - 1; j++) {
						ROTATE(a, j, ip, j, iq);
					}
					for (j = ip + 1; j <= iq - 1; j++) {
						ROTATE(a, ip, j, j, iq);
					}
					for (j = iq + 1; j <= n; j++) {
						ROTATE(a, ip, j, iq, j);
					}
					for (j = 1; j <= n; j++) {
						ROTATE(v, j, ip, j, iq);
					}
					++(*nrot);
				}
			}
		}
		for (ip = 1; ip <= n; ip++) {
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
	//cout << "Too many iterations in routine jacobi" << endl;

}

void eigsrt(REAL d[], REAL **v, int n)
{
	/* Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] as
	output from jacobi, this routine sorts the eigenvalues into descending
	order, and rearranges the columns of v correspondingly. The method
	is straight insertion.
	*/
	int k, j, i;
	REAL p;

	for (i = 1; i<n; i++)
	{
		p = d[k = i];
		for (j = i + 1; j <= n; j++)
		if (d[j] >= p) p = d[k = j];
		if (k != i) {
			d[k] = d[i];
			d[i] = p;
			for (j = 1; j <= n; j++) {
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	}
}


void calcPointNoise(int num_nlist, INT32 *nlist, REAL_DIM vtmp2, REAL noise[], REAL_DIM* vert_v)
{

	REAL **U = dmatrix(1, num_nlist, 1, 6);
	REAL *B = dvector(1, num_nlist);

	int j, k, l, id2, nrot;
	REAL_DIM centroid, vtmp;

	// determine local coordinate of neighboring points
	REAL sum;
	REAL_DIM *c = new REAL_DIM[num_nlist];



	// determine centroid
	for (k = 0; k<DIM; k++) centroid[k] = 0.0;
	for (j = 0; j<num_nlist; j++)
	for (k = 0; k<DIM; k++)
		//centroid[k] += vert[nlist[j]]._v[k];
		centroid[k] += vert_v[nlist[j]][k];
	for (k = 0; k<DIM; k++) centroid[k] /= float(num_nlist);

	//begin build covariance matrix
	//first diagonal element
	REAL *w = new REAL[num_nlist];

	sum = 0;
	for (j = 0; j<num_nlist; j++)
	{
		//sum+=(vert[nlist[j]]._v[0]-centroid[0])*(vert[nlist[j]]._v[0]-centroid[0]);
		sum += (vert_v[nlist[j]][0] - centroid[0])*(vert_v[nlist[j]][0] - centroid[0]);
	}

	A[1][1] = sum / num_nlist;

	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[0]-centroid[0])*(vert[nlist[j]]._v[1]-centroid[1]);
		sum += (vert_v[nlist[j]][0] - centroid[0])*(vert_v[nlist[j]][1] - centroid[1]);

	A[1][2] = sum / num_nlist;
	A[2][1] = A[1][2];


	//second diagonal element
	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[1]-centroid[1])*(vert[nlist[j]]._v[1]-centroid[1]);
		sum += (vert_v[nlist[j]][1] - centroid[1])*(vert_v[nlist[j]][1] - centroid[1]);

	A[2][2] = sum / num_nlist;

	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[0]-centroid[0])*(vert[nlist[j]]._v[2]-centroid[2]);
		sum += (vert_v[nlist[j]][0] - centroid[0])*(vert_v[nlist[j]][2] - centroid[2]);

	A[1][3] = sum / num_nlist;
	A[3][1] = A[1][3];

	//third diagonal element
	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[2]-centroid[2])*(vert[nlist[j]]._v[2]-centroid[2]);
		sum += (vert_v[nlist[j]][2] - centroid[2])*(vert_v[nlist[j]][2] - centroid[2]);

	A[3][3] = sum / num_nlist;

	sum = 0;
	for (j = 0; j<num_nlist; j++)
		//sum+=(vert[nlist[j]]._v[1]-centroid[1])*(vert[nlist[j]]._v[2]-centroid[2]);
		sum += (vert_v[nlist[j]][1] - centroid[1])*(vert_v[nlist[j]][2] - centroid[2]);

	A[2][3] = sum / num_nlist;
	A[3][2] = A[2][3];

	/* d[1..3] returns the eigenvalues of a. v[1..3][1..3] is
	a matrix whose columns contain, on output, the normalized
	eigenvectors of a.
	*/
	
	jacobi(A, 3, D, V3, &nrot);

	/* output of eigenvalues in descending order */
	eigsrt(D, V3, 3);

	for (j = 0; j<DIM; j++)
	{
		vtmp2[j] = D[j + 1];
	}

	for (k = 0; k<DIM; k++)
	for (l = 0; l<DIM; l++)
		V[k][l] = V3[l + 1][k + 1];

	dvector_cross_product(V[0], V[1], vtmp, 3, 0);
	if (dvector_dot_product(vtmp, V[2], 3, 0) < 0.0)
	{
		dvector_assign(V[0], vtmp, 3);
		dvector_assign(V[1], V[0], 3);
		dvector_assign(vtmp, V[1], 3);
	}

	// ***determine local coordinate of neighboring points
	calc_Cartesian_frame_to_frame_transformation_matrix(V[0], V[1], V[2],
		centroid, E);

	__int16 num_n = 0;
	for (j = 0; j<num_nlist; j++)
	{
		id2 = nlist[j];			// global vertex id

		for (k = 0; k<DIM; k++)
		{
			c[num_n][k] = 0.0;
			for (l = 0; l<DIM; l++)
				S[l] = vert_v[id2][l];
			S[3] = 1.0;
			for (l = 0; l<4; l++)
				c[num_n][k] += E[l][k] * S[l];
		}
		num_n++;
	}

	// local coordinate for centroid
	for (k = 0; k<DIM; k++)
	{
		vtmp[k] = 0.0;
		for (l = 0; l<DIM; l++)
			S[l] = centroid[l];
		S[3] = 1.0;
		for (l = 0; l<4; l++)
			vtmp[k] += E[l][k] * S[l];
	}


	// ***transform back to the Cartesian coordinate
	calc_frame_to_Cartesian_frame_transformation_matrix(V[0], V[1], V[2],
		centroid, F);

	__int16 m, m2;
	INT32 id;
	REAL x, y, z, wmax, wmin;
	// determine curvatures
	REAL a[7];


		//Error *w3 = new Error[num_nlist];

	m = num_nlist;
	for (j = 0; j<m; j++)
		w[j] = 1.0;


	// ***determine local quadric coordinate patch
	for (j = 0; j<m; j++)
	{
		x = c[j][0], y = c[j][1], z = c[j][2];

		U[j + 1][1] = w[j] * x*x, U[j + 1][2] = w[j] * x*y, U[j + 1][3] = w[j] * y*y;
		U[j + 1][4] = w[j] * x, U[j + 1][5] = w[j] * y, U[j + 1][6] = w[j] * 1.0;
		B[j + 1] = w[j] * z;
	}

	// SVD the square matrix U
	svdcmp(U, m, 6, W, V2);
	wmax = 0.0;              /// Will be the maximum sigular value obtained
	for (j = 1; j <= 6; j++)
	if (W[j] > wmax)   wmax = W[j];
	// This is where we set the threshold for singular values allowed to be 
	// nonzero. The constant is typical, but not universal. You have to 
	// experiment with your own application
	wmin = wmax*1.0e-6;
	for (j = 1; j <= 6; j++)
	if (W[j] < wmin)  W[j] = 0.0;
	// Backsubstitute
	svbksb(U, W, V2, m, 6, B, X);

	// ***determine curvatures
	// quadratic surface: z = a1*x^2 + a2*xy + a3*y^2 + a4 x + a5 y + a6
	a[1] = X[1], a[2] = X[2], a[3] = X[3], a[4] = X[4], a[5] = X[5], a[6] = X[6];

	REAL tmp, tmp2;

	for (j = 0; j<m; j++)
	{
		x = c[j][0], y = c[j][1];
		tmp = a[1] * x*x + a[2] * x*y + a[3] * y*y + a[4] * x + a[5] * y + a[6];

		tmp2 = tmp - c[j][2];
		//noise += fabs(tmp2);
		noise[j] = fabs(tmp2);
	}

	//noise /= float(m);

	free_dmatrix(U, 1, num_nlist, 1, 6);
	free_dvector(B, 1, num_nlist);

	delete[]c;
	delete[]w;


}

void svbksb(REAL **u, REAL w[], REAL **v, int m, int n, REAL b[], REAL x[])
{
	// Solves A X = B for a vector X, where A is specified by the arrays 
	// u[1..m][1..n], w[1..n], v[1..n][1..n] as returned by svdcmp. m and n
	// are the dimensions of a, and will be equal for square matrices. 
	// b[1..m] is the input right-hand side. x[1..n] is the outptu solution
	// vector. No input quantities are destroyed, so the routine may be called
	// seqeuntially with different b's. 

	int jj, j, i;
	REAL s, *tmp;

	tmp = dvector(1, n);
	// Calculate U^T B
	for (j = 1; j <= n; j++)
	{
		s = 0.0;
		if (w[j])
		{
			for (i = 1; i <= m; i++) s += u[i][j] * b[i];
			s /= w[j];
		}
		tmp[j] = s;
	}

	// Matrix multiply by V to get answer
	for (j = 1; j <= n; j++)
	{
		s = 0.0;
		for (jj = 1; jj <= n; jj++) s += v[j][jj] * tmp[jj];
		x[j] = s;
	}
	free_dvector(tmp, 1, n);
}
void svdcmp(REAL **a, int m, int n, REAL w[], REAL **v)
{
	// Given a matrix a[1..m][1..n], this routine computes its singular value decomposition,
	// A = U W V^T. The matrix U replaces a on output. The diagonal matrix of singular vaules
	// W is output as a vector w[1..n]. The matrix V (not the transpose V^T) is output as 
	// v[1..n][1..n].

	int flag, i, its, j, jj, k, l, nm;
	REAL anorm, c, f, g, h, s, scale, x, y, z, *rv1;

	rv1 = dvector(1, n);
	g = scale = anorm = 0.0;
	// Householder reduction to bidiagonal form
	for (i = 1; i <= n; i++)
	{
		l = i + 1;
		rv1[i] = scale*g;
		g = s = scale = 0.0;
		if (i <= m)
		{
			for (k = i; k <= m; k++) scale += fabs(a[k][i]);
			if (scale)
			{
				for (k = i; k <= m; k++)
				{
					a[k][i] /= scale;
					s += a[k][i] * a[k][i];
				}
				f = a[i][i];
				g = -SIGN(sqrt(s), f);
				h = f*g - s;
				a[i][i] = f - g;
				for (j = l; j <= n; j++)
				{
					for (s = 0.0, k = i; k <= m; k++) s += a[k][i] * a[k][j];
					f = s / h;
					for (k = i; k <= m; k++) a[k][j] += f*a[k][i];
				}
				for (k = i; k <= m; k++) a[k][i] *= scale;
			}
		}
		w[i] = scale*g;
		g = s = scale = 0.0;
		if (i <= m && i != n)
		{
			for (k = l; k <= n; k++) scale += fabs(a[i][k]);
			if (scale)
			{
				for (k = l; k <= n; k++)
				{
					a[i][k] /= scale;
					s += a[i][k] * a[i][k];
				}
				f = a[i][l];
				g = -SIGN(sqrt(s), f);
				h = f*g - s;
				a[i][l] = f - g;
				for (k = l; k <= n; k++) rv1[k] = a[i][k] / h;
				for (j = l; j <= m; j++)
				{
					for (s = 0.0, k = l; k <= n; k++) s += a[j][k] * a[i][k];
					for (k = l; k <= n; k++) a[j][k] += s*rv1[k];
				}
				for (k = l; k <= n; k++) a[i][k] *= scale;
			}
		}
		anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
	}

	// Accumulation of right hand transformations
	for (i = n; i >= 1; i--)
	{
		if (i<n)
		{
			if (g)
			{
				for (j = l; j <= n; j++)    // Double division to avoid possible underflow
					v[j][i] = (a[i][j] / a[i][l]) / g;
				for (j = l; j <= n; j++)
				{
					for (s = 0.0, k = l; k <= n; k++) s += a[i][k] * v[k][j];
					for (k = l; k <= n; k++)  v[k][j] += s*v[k][i];
				}
			}
			for (j = l; j <= n; j++) v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}

	// Accumulation of left-hand transformations
	for (i = IMIN(m, n); i >= 1; i--)
	{
		l = i + 1;
		g = w[i];
		for (j = l; j <= n; j++)  a[i][j] = 0.0;
		if (g)
		{
			g = 1.0 / g;
			for (j = l; j <= n; j++)
			{
				for (s = 0.0, k = l; k <= m; k++) s += a[k][i] * a[k][j];
				f = (s / a[i][i])*g;
				for (k = i; k <= m; k++) a[k][j] += f*a[k][i];
			}
			for (j = i; j <= m; j++) a[j][i] *= g;
		}
		else
		for (j = i; j <= m; j++) a[j][i] = 0.0;
		++a[i][i];
	}

	// Diagonalization of the bidiagonal form: Loop over singular values, and
	// over allowed iterations
	for (k = n; k >= 1; k--)
	{
		for (its = 1; its <= 30; its++)
		{
			flag = 1;
			// test for splitting
			for (l = k; l >= 1; l--)
			{
				nm = l - 1;                       // Note that rv1[1] is always zero
				if ((REAL)(fabs(rv1[l]) + anorm) == anorm)
				{
					flag = 0;
					break;
				}
				if ((REAL)(fabs(w[nm]) + anorm) == anorm) break;
			}
			if (flag)
			{
				c = 0.0;                     // Cancellation of rv1[l], if l>1
				s = 1.0;
				for (i = l; i <= k; i++)
				{
					f = s*rv1[i];
					rv1[i] = c*rv1[i];
					if ((REAL)(fabs(f) + anorm) == anorm) break;
					g = w[i];
					h = pythag(f, g);
					w[i] = h;
					h = 1.0 / h;
					c = g*h;
					s = -f*h;
					for (j = 1; j <= m; j++)
					{
						y = a[j][nm];
						z = a[j][i];
						a[j][nm] = y*c + z*s;
						a[j][i] = z*c - y*s;
					}
				}
			}
			z = w[k];
			if (l == k)                 // Convergence
			{
				if (z<0.0)            // Singular value is made nonnegative
				{
					w[k] = -z;
					for (j = 1; j <= n; j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30)
				cout << "No convergence in 30 svdcmp iterations." << endl;
			x = w[l];
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z)*(y + z) + (g - h)*(g + h)) / (2.0*h*y);
			g = pythag(f, 1.0);
			f = ((x - z)*(x + z) + h*((y / (f + SIGN(g, f))) - h)) / x;
			c = s = 1.0;
			// Next QR transformation
			for (j = l; j <= nm; j++)
			{
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x*c + g*s;
				g = g*c - x*s;
				h = y*s;
				y *= c;
				for (jj = 1; jj <= n; jj++)
				{
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x*c + z*s;
					v[jj][i] = z*c - x*s;
				}
				z = pythag(f, h);
				w[j] = z;                    // Rotation can be arbitrary if z = 0
				if (z)
				{
					z = 1.0 / z;
					c = f*z;
					s = h*z;
				}
				f = c*g + s*y;
				x = c*y - s*g;
				for (jj = 1; jj <= m; jj++)
				{
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y*c + z*s;
					a[jj][i] = z*c - y*s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
	free_dvector(rv1, 1, n);
}

REAL pythag(REAL a, REAL b)
{
	// Computes (a^2+b^2)^0.5 without destructive underflow or overflow
	REAL absa, absb;
	absa = fabs(a);
	absb = fabs(b);
	if (absa>absb)
		return absa*sqrt(1.0 + SQR(absb / absa));
	else
		return (absb == 0.0 ? 0.0 : absb*sqrt(1.0 + SQR(absa / absb)));
}

void calc_Cartesian_frame_to_frame_transformation_matrix(REAL_DIM U, REAL_DIM V, REAL_DIM W,
	REAL_DIM O, REAL E[][4])
{
	// Input:  U, V, W are three direction vector of the general frame
	//         O - is the origin of the general frame
	REAL d1, d2, d3;
	REAL_DIM nO;

	// case 1:  t = [1 0 0]
	d1 = U[0], d2 = V[0], d3 = W[0];
	E[0][0] = d1, E[0][1] = d2, E[0][2] = d3, E[0][3] = 0.0;

	// case 2: t = [0 1 0]
	d1 = U[1], d2 = V[1], d3 = W[1];
	E[1][0] = d1, E[1][1] = d2, E[1][2] = d3, E[1][3] = 0.0;

	// case 3: t = [0 0 1]
	d1 = U[2], d2 = V[2], d3 = W[2];
	E[2][0] = d1, E[2][1] = d2, E[2][2] = d3, E[2][3] = 0.0;

	// case 4: t = [0 0 0] - O
	dvector_assign(O, nO, 3);
	dvector_negation(nO, 3, 0);
	d1 = dvector_dot_product(U, nO, 3, 0);
	d2 = dvector_dot_product(V, nO, 3, 0);
	d3 = dvector_dot_product(W, nO, 3, 0);

	E[3][0] = d1, E[3][1] = d2, E[3][2] = d3, E[3][3] = 1.0;

}
void calc_frame_to_Cartesian_frame_transformation_matrix(REAL_DIM U, REAL_DIM V, REAL_DIM W,
	REAL_DIM O, REAL E[][4])
{
	// Input:  U, V, W are three direction vector of the general frame
	//         O - is the origin of the general frame

	// vector u
	E[0][0] = U[0], E[0][1] = U[1], E[0][2] = U[2], E[0][3] = 0.0;

	// vector v
	E[1][0] = V[0], E[1][1] = V[1], E[1][2] = V[2], E[1][3] = 0.0;

	// vector w
	E[2][0] = W[0], E[2][1] = W[1], E[2][2] = W[2], E[2][3] = 0.0;

	// vector O
	E[3][0] = O[0], E[3][1] = O[1], E[3][2] = O[2], E[3][3] = 1.0;

}


REAL dvector_dot_product(char a0[], REAL* b, int n, int m)
{
	REAL_DIM a;
	for (int i = 0; i<DIM; i++)
	{
		a[i] = a0[i];
	}
	return dvector_dot_product(a, b, n, m);
}

REAL dvector_dot_product(char a0[], char b0[], int n, int m)
{
	REAL_DIM a, b;
	for (int i = 0; i<DIM; i++)
	{
		a[i] = a0[i], b[i] = b0[i];
	}
	return dvector_dot_product(a, b, n, m);
}

REAL dvector_dot_product(REAL *a, REAL *b, int n, int m)
{
	// Input: n -- dimension
	//        m -- starting index 0 or 1

	int i;
	REAL result = 0.0;

	if (m == 1)
	for (i = 1; i <= n; i++)
		result += a[i] * b[i];
	else if (m == 0)
	for (i = 0; i< n; i++)
		result += a[i] * b[i];

	return(result);
}

void dvector_assign(REAL_DIM source, REAL_DIM target, int dim)
{
	for (int i = 0; i< dim; i++)
		target[i] = source[i];
}

void dvector_negation(REAL_DIM x, int n, int m)
{
	// Input: n -- dimension
	//        m -- starting index 0 or 1

	int i;
	REAL result = 0.0;

	if (m == 1)
	for (i = 1; i <= n; i++)
		x[i] = -x[i];
	else if (m == 0)
	for (i = 0; i< n; i++)
		x[i] = -x[i];
}



REAL **dmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a REAL matrix with subscript range m[nrl..nrh][ncl..nch] */
/* **m <==> *m[row]  */
{
	int i, nrow, ncol;
	REAL **m;

	nrow = nrh - nrl + 1;
	ncol = nch - ncl + 1;

	/* allocate pointers to rows */
	m = (REAL **)malloc((size_t)((nrow + NR_END)*sizeof(REAL*)));
	if (!m) cout << "allocation failure 1 in dmatrix()\n";
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl] = (REAL *)malloc((size_t)((nrow*ncol + NR_END)*sizeof(REAL)));
	if (!m[nrl]) cout << "allocation failure 2 in dmatrix()";
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

REAL *dvector(int nl, int nh)
/* allocate a REAL vector with subscript range v[nl..nh] */
{
	REAL *v;

	v = (REAL *)malloc((size_t)((nh - nl + 1 + NR_END)*sizeof(REAL)));
	if (!v) cout << "allocation failure in dvector() \n";
	return v - nl + NR_END;
}


void free_dvector(REAL *v, int nl, int nh)
/* free a REAL vector allocated by dvector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}

void free_dmatrix(REAL **m, int nrl, int nrh, int ncl, int nch)
/* free a REAL matrix allocated by dmatrix() */
{
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
}

void add_M(REAL_33 M, REAL_DIM p)
{
	REAL x = p[0], y = p[1], z = p[2], t;
	M[0][0] += x*x;
	t = x*y;
	M[0][1] += t;
	M[1][0] += t;
	t = x*z;
	M[0][2] += t;
	M[2][0] += t;
	M[1][1] += y*y;
	t = y*z;
	M[1][2] += t;
	M[2][1] += t;
	M[2][2] += z*z;
}

void dvector_cross_product(REAL *A, REAL *B, REAL *C, INT16 n, INT16 m)
{
	// Input: n -- dimension 
	//        m -- starting index (0, 1) 



	if (n == 3)  // dimension
	{
		if (m == 1)
		{
			C[1] = A[2] * B[3] - B[2] * A[3];
			C[2] = B[1] * A[3] - A[1] * B[3];
			C[3] = A[1] * B[2] - B[1] * A[2];
		}
		else if (m == 0)
		{
			C[0] = A[1] * B[2] - B[1] * A[2];
			C[1] = B[0] * A[2] - A[0] * B[2];
			C[2] = A[0] * B[1] - B[0] * A[1];
		}
	}
	else
	{
		cerr << "Error: wrong n in dvector_cross_product()." << endl;
	}
}

int dvector_normalize(REAL *V_in, REAL *V_out, INT16 n, INT16 m)
{
	// Input: n -- dimension 
	//        m -- starting index (0, 1) 

	REAL d = 0.0, Tol = 1.0E-36;
	int i;

	if (m == 1)
	{
		for (i = 1; i <= n; i++)
		{
			d += V_in[i] * V_in[i];
		}

		d = sqrt(d);
		if (d > Tol)
		{
			// normalize X_new		
			for (i = 1; i <= 3; i++)
				V_out[i] = V_in[i] / d;
		}
		else
		{
#ifdef _DEBUG
			cerr << "Error: d is zero in dvector_normalize() d= " <<
				d << endl;
#endif
			return(ERROR);
		}
	}
	else if (m == 0)
	{
		for (i = 0; i< n; i++)
		{
			d += V_in[i] * V_in[i];
		}

		d = sqrt(d);
		if (d > Tol)
		{
			// normalize X_new		
			for (i = 0; i< 3; i++)
				V_out[i] = V_in[i] / d;
		}
		else
		{
#ifdef _DEBUG
			cerr << "Error: d is zero in dvector_normalize(). d= " <<
				d << endl;
#endif
			return(ERROR);
		}
	}

	return(TRUE);
}


REAL dist_pt(const REAL_DIM a, const REAL_DIM b)
{
	return sqrt(SQR(a[0] - b[0]) + SQR(a[1] - b[1]) + SQR(a[2] - b[2]));
}

REAL NormalizeEigenValue(REAL_DIM vtmp2)
{
	//return 3*vtmp2[2] / (vtmp2[0]+vtmp2[1]+vtmp2[2]);
	return vtmp2[2] / vtmp2[0];

	//return (vtmp2[2]) / (vtmp2[0] * vtmp2[1]);

}

// circular propagation for noised data
void propagateVoxel(int i0, int j0, int k0, VoxelDataStruc ***vds, int cubelen, int **** list)
{
	INT16 num_nlist = 0, num_nlist2 = 0;
	int i, j, nr = 30, id2, m;   // for noised model: nr=30; for non-noised model: nr=15
	int nR = 4 * nr, nr2, nR2;
	int num_s1, num_s2, s1[100], s2[100], id, k;
	//int dim=FRONT_SIZE;
	int dim = curvert / 4;
	CoordDataStruc *wave_front = new CoordDataStruc[dim];
	CoordDataStruc *wave_front2 = new CoordDataStruc[dim];
	int num_wave_front = 0, num_wave_front2 = 0;
	//bool *flag=new bool[curvert];
	REAL size = 1.0;
	bool flag2 = false;
	CoordDataStruc nid;
	nid.i = i0, nid.j = j0, nid.k = k0;
	int i1, j1, k1, i2, j2, k2, l;
	float A, B, C, D, d;
	float scale = 1*min_cubesize; // best for this case
	float tmp, angle_t = 0; //0.5*M_PI; // 9 degree, best for this case
	//REAL angle_t2 = 0.75*M_PI;
	float dist_angl_ratio = 0.44;
	float dist_angl;
	//float dist_angl_threshold = 0.3;
	float dist_angl_threshold = 0;

	num_wave_front = 1;
	wave_front[0] = nid;
	REAL_DIM v1, v2;

	do
	{
		num_wave_front2 = 0;
		for (m = 0; m<num_wave_front; m++)
		{

			i1 = wave_front[m].i, j1 = wave_front[m].j, k1 = wave_front[m].k; // x0, y0, z0
			//if( (Gstatus_c[id]==UNPROCESSED && Gindex_c[id]<=false) )
			if (vds[i1][j1][k1].connectivityFlag <= 0) // untouched or unprocessed
			{
				vds[i1][j1][k1].connectivityFlag = 2; // true data

				// construct the plane equation
				// centroid is shifted to the lower left corner
				
				// distance-based
				A = vds[i1][j1][k1].n[0], B = vds[i1][j1][k1].n[1], C = vds[i1][j1][k1].n[2];
				D = -(A*vds[i1][j1][k1].centroid[0] + B*vds[i1][j1][k1].centroid[1] + C*vds[i1][j1][k1].centroid[2]);
				for (l = 0; l<DIM; l++)
					v1[l] = vds[i1][j1][k1].n[l];

				for (i = -1; i <= 1; i++)
				for (j = -1; j <= 1; j++)
				for (k = -1; k <= 1; k++)
				{

					i2 = i1 + i, j2 = j1 + j, k2 = k1 + k;  // x1, y1, z1

					if (i == 0 && j == 0 && k == 0)
						continue;
					if (i2<0 || j2 <0 || k2<0)
						continue;
					if (i2>cubelen - 1 || j2>cubelen - 1 || k2>cubelen - 1)
						continue;
//					if (abs(i) + abs(j) + abs(k) != 1)
//						continue;
					
					// no difference
					int kk;
					int aa, bb, cc;
					for (kk = 0; kk < vds[i2][j2][k2].count; kk++)
					{
						if (list[i2][j2][k2][kk] == 14270)
						{
							aa = i2;
							bb = j2;
							cc = k2;
						}
					}




					if (vds[i2][j2][k2].count == 0)
					{
						vds[i2][j2][k2].connectivityFlag = 3; // outlier
						continue;
					}
					

					// distance-based
					// distance from (x1, y1, z1) to the plane (A,B,C,D)
					d = fabs(A*vds[i2][j2][k2].centroid[0] + B*vds[i2][j2][k2].centroid[1] + C*vds[i2][j2][k2].centroid[2] + D) / sqrt(A*A + B*B + C*C);

//					for (int l = 0; l<DIM; l++)
//						v2[l] = vds[i2][j2][k2].centroid[l] - vds[i1][j1][k1].centroid[l];
					for (l = 0; l<DIM; l++)
						v2[l] = vds[i2][j2][k2].n[l];


					// anlge-based
					tmp = angle_vect(v1, v2);
					/*
					if (tmp > M_PI / 2)
					{
						cout << tmp << endl;
						cout<< list[i2][j2][k2][0] << endl;
					}
					*/

					// linear combination of distance and angle
					//dist_angl = (1 - dist_angl_ratio)*(tmp - angle_t) / (M_PI / 2) + \
						//dist_angl_ratio*(d - scale) / min_cubesize;
					dist_angl = (1 - dist_angl_ratio)*(tmp - angle_t) / (M_PI/2) + \
						dist_angl_ratio*(d - scale) / min_cubesize;

					//if (d < scale) // in the propagation plane
					//if (tmp < angle_t)
					//if ((tmp < angle_t) || d < scale)
					if (dist_angl < dist_angl_threshold)
					// The purpose of the second condition is to get rid of cubes with no or very few points
					{
						if (i2 == aa && j2 == bb && k2 == cc)
						{
							aa = aa;
						}
						if (vds[i2][j2][k2].connectivityFlag == -1)  // untouched
						{
							nid.i = i2, nid.j = j2, nid.k = k2;
							wave_front2[num_wave_front2++] = nid;
							vds[i2][j2][k2].connectivityFlag = 0;  // unprocessed
						}
						else if (vds[i2][j2][k2].connectivityFlag == 0) // set in former loop
						{
							// do nothing
						}
						else if (vds[i2][j2][k2].connectivityFlag == 2)
						{
							// do nothing
						}
						else // treated as noise before
						{
							vds[i2][j2][k2].connectivityFlag = 2; // true data

						}

					}
					else // outside the propagation plane
					{
						if (vds[i2][j2][k2].connectivityFlag == -1)
						{
							vds[i2][j2][k2].connectivityFlag = 3; // outlier
						}
					}
				}

			}
		}  // end of for(m=0; m<num_wave_front; m++)


		if (num_wave_front2 == 0)
			break;

		num_wave_front = num_wave_front2;
		for (i = 0; i<num_wave_front; i++)
			wave_front[i] = wave_front2[i];

	} while (num_wave_front > 0);   // end of do


	delete[]wave_front;
	delete[]wave_front2;

}

// determine centroid
void calcCentroid(int num_nlist, INT32 *nlist, REAL_DIM* vert_v, REAL centroid[])
{
	int j, k;
	for (k = 0; k<DIM; k++) centroid[k] = 0.0;
	for (j = 0; j<num_nlist; j++)
	for (k = 0; k<DIM; k++)
		centroid[k] += vert_v[nlist[j]][k];
	for (k = 0; k<DIM; k++) centroid[k] /= float(num_nlist);
}

REAL angle_vect(REAL_DIM a, REAL_DIM b)
{
	int i;
	REAL cp = 0;
	REAL d1 = 0;
	REAL d2 = 0;
	REAL tol = 1.0e-6;

	for (i = 0; i < DIM; i++)
	{
		cp += a[i] * b[i];
		d1 += a[i] * a[i];
		d2 += b[i] * b[i];
	}

	//if (d1 == 0 || d2 == 0) 
	if (d1 <= tol || d2 <= tol)
	{
#ifdef _DEBUG
		cerr << "d1= " << d1 << " d2= " << d2 << endl;
		cerr << "Identical nodes !" << endl;
#endif
		cp = 0;
	}
	else
	{
		cp /= sqrt(d1*d2);
		cp = fabs(cp); // added line
		if (fabs(cp - 1.0) < tol)
			cp = 0.0;
		else if (fabs(cp + 1.0) < tol)
			cp = M_PI;
		else
			cp = acos((cp));
	}

	return cp;
}