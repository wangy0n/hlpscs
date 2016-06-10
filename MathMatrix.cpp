//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	MathMatrix.cpp
//OBJECTIVE:	Mathematical matrix operations
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#include "MathMatrix.h"

//---------------------------------------------------------------------------
//构造函数

//构造函数（M*N）初始化零
MathMatrix::MathMatrix(const int& M, const int& N)
{
	vector<double> temp(N,0.0);
	element.assign(M,temp);
}
MathMatrix::MathMatrix(const vector<double>& vec2D)
{
	int N = int(vec2D.size());
	vector<double> temp(N);
	element.assign(1,temp);
	for (int i=0; i<N; i++)
	{
		element[0][i] = vec2D[i];
	}
}
MathMatrix::MathMatrix(const vector<vector<double> >& vec2D)
{
	element = vec2D;
}
//构造函数（数组构造1*N）
MathMatrix::MathMatrix(const double *vec, const int& N)
{
	if (N<=0)
	{
		cout << "申请矩阵有错！N<=0 " << endl;
		exit (0);
	}
	vector<double> temp(N);
	element.assign(1,temp);
	for (int i=0; i<N; i++)
	{
		element[0][i] = vec[i];
	}
}
//构造函数（数组构造M*N）
MathMatrix::MathMatrix(const double *vec, const int& M, const int& N)
{
	if (M<=0 || N<=0)
	{
		cout << "申请矩阵有错！M<=0或N<=0 " << endl;
		exit (0);
	}
	vector<double> temp(N);
	element.assign(M,temp);
	if (N!=1)
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				element[i][j] = vec[M*i+j];
			}
		}
	}
	else 
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				element[i][j] = vec[i];
			}
		}
	}
}
//---------------------------------------------------------------------------
//一维数组赋给矩阵 
void MathMatrix::Evalu(const double *vec, const int& N)
{
	if (N<=0)
	{
		cout << "申请矩阵有错！N<=0 " << endl;
		exit (0);
	}
	vector<double> temp(N);
	element.assign(1,temp);
	for (int i=0; i<N; i++)
	{
		element[0][i] = vec[i];
	}
}
//---------------------------------------------------------------------------
//二维数组赋给矩阵 
void MathMatrix::Evalu(const double *vec, const int& M, const int& N)
{
	if (M<=0 || N<=0)
	{
		cout << "申请矩阵有错！M<=0或N<=0 " << endl;
		exit (0);
	}
	vector<double> temp(N);
	element.assign(M,temp);
	if (N!=1)
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				element[i][j] = vec[M*i+j];
			}
		}
	}
	else 
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				element[i][j] = vec[i];
			}
		}
	}
}
//---------------------------------------------------------------------------
//重载赋值操作符（一维向量赋给矩阵） 
void MathMatrix::operator=(const vector<double>& vec)
{
	element.assign(1,vec);
}
//重载赋值操作符（二维向量赋给矩阵）
void MathMatrix::operator=(const vector<vector<double> >& vec)
{
	element = vec;
}
//---------------------------------------------------------------------------
//重载乘法操作符（矩阵乘矩阵）
MathMatrix MathMatrix::operator*(const MathMatrix& matrix)
{
	int L = 0;
	if (element[0].size() != matrix.element.size())
	{
		cout << "矩阵不能相乘！列号不等于行号！" << endl;
		exit(0);
	}
	else
	{
		L = int(matrix.element.size());
	}
	//待求矩阵：M行号与N列号
	int M = int(element.size());
	int N = int(matrix.element[0].size());
	//定义待求矩阵Tem_mat
	MathMatrix Tem_mat(M,N);
	//累加器
	double conter;

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			conter = 0.0;
			for(int k=0; k<L; k++)
			{
				conter = conter + element[i][k] * matrix.element[k][j];
			}
			Tem_mat.element[i][j] = conter;
		}
	}

	return Tem_mat;
}
//---------------------------------------------------------------------------
//重载乘法操作符（矩阵乘实数）
MathMatrix MathMatrix::operator*(const double& R)
{
	//待求矩阵：M行号与N列号
	int M = int(element.size());
	int N = int(element[0].size());
	//定义待求矩阵Tem_mat
	MathMatrix Tem_mat(M,N);

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			Tem_mat.element[i][j] = element[i][j] * R;
		}
	}

	return Tem_mat;
}
//---------------------------------------------------------------------------
//重载加法操作符（矩阵加矩阵）
MathMatrix MathMatrix::operator+(const MathMatrix& matrix)
{
	if ((element.size() != matrix.element.size())||(element[0].size() != matrix.element[0].size()))
	{
		cout << "矩阵不能相减！列号或行号不相等！" << endl;
		exit(0);
	}
	//待求矩阵：M行号与N列号
	int M = int(element.size());
	int N = int(element[0].size());
	//定义待求矩阵Tem_mat
	MathMatrix Tem_mat(M,N);

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			Tem_mat.element[i][j] = element[i][j] + matrix.element[i][j];
		}
	}

	return Tem_mat;
}
//---------------------------------------------------------------------------
//重载加法操作符（矩阵加实数）
double MathMatrix::operator+(const double& R)
{
	//待求矩阵：M行号与N列号
	int M = int(element.size());
	int N = int(element[0].size());
	if (M!=N||M==1)
	{
		cout << "矩阵不是1*1的，不能和一个实数相加！" << endl;
		exit(0);
	}
	//定义待求矩阵
	double R2 = element[0][0] + R;

	return R2;
}
//---------------------------------------------------------------------------
//重载减法操作符（矩阵减矩阵）
MathMatrix MathMatrix::operator-(const MathMatrix& matrix)
{
	if ((element.size() != matrix.element.size())||(element[0].size() != matrix.element[0].size()))
	{
		cout << "矩阵不能相减！列号或行号不相等！" << endl;
		exit(0);
	}
	//待求矩阵：M行号与N列号
	int M = int(element.size());
	int N = int(element[0].size());
	//定义待求矩阵Tem_mat
	MathMatrix Tem_mat(M,N);

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			Tem_mat.element[i][j] = element[i][j] - matrix.element[i][j];
		}
	}

	return Tem_mat;
}
//---------------------------------------------------------------------------
//判定对称矩阵 输出-1：不是方阵，输出0：不对称，输出1：对称
int MathMatrix::Symmetry(void)
{
	int M = int(element.size());
	int N = int(element[0].size());

	if (M==N)
	{
		for (int i=0; i<M; i++)
		{
			for (int j=i+1; j<M; j++)
			{
				if (fabs(element[i][j]-element[j][i])>Zero)
				{
					return 0;
				}
			}
		}
	}
	else
	{
		return -1;
	}
	return 1;
}
//---------------------------------------------------------------------------
//矩阵求逆
MathMatrix MathMatrix::Inverse(void)
{
	vector<vector<double> > vec = element;
	//--------------------------------------------------
	//判定对称正定矩阵

	//--------------------------------------------------
	//求逆矩阵
	int stRank = int(vec.size());
	vector<double> b(stRank,0.0);
	for(int k=0; k<stRank; k++)
	{
		double w= vec[0][0];
		int m = stRank - k -1;
		for(int i=1; i<stRank; i++)
		{
			double g = vec[i][0];
			b[i] = g / w;
			if (i<=m)
			{
				b[i] = -b[i];
			}
			for(int j=1; j<=i; j++)
			{
				vec[i-1][j-1] = vec[i][j] + g * b[j];
			}
		}
		vec[stRank-1][stRank-1] = 1.0 / w;
		for(int i= 1; i<stRank; i++)
		{
			vec[stRank-1][i-1] =  b[i];
		}
	}
	for(int i=0; i<stRank-1; i++)
	{
		for(int j = i+1; j<stRank; j++)
		{
			vec[i][j] = vec[j][i];
		}
	}

	MathMatrix matrix = vec;
	return matrix;
}
//---------------------------------------------------------------------------
//矩阵转置
MathMatrix MathMatrix::Transpose(void)
{
	int M = int(element.size());
	int N = int(element[0].size());

	vector<vector<double> > vec;
	vector<double> temp(M);
	vec.assign(N,temp);

	for(int i=0; i<N; i++)
	{
		for(int j=0; j<M; j++)
		{
			vec[i][j] = element[j][i];
		}
	}

	MathMatrix matrix = vec;
	return matrix;
}
//---------------------------------------------------------------------------
//求出矩阵一列，将其转为矩阵
MathMatrix MathMatrix::GetCol(int cn){
	MathMatrix matrix(int(element.size()),1);
	for(int i=0;i<int(element.size());i++)
		matrix.element[i][0]=element[i][cn];
	return matrix;
}                  
//---------------------------------------------------------------------------
//求矩阵行数
int MathMatrix::RowN(void)
{
	return int(element.size());
}
//---------------------------------------------------------------------------
//求矩阵行数
int MathMatrix::CalN(void)
{
	return int(element[0].size());
}
//============================================================================
