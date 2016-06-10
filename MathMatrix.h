//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	MathMatrix.h
//OBJECTIVE:	Mathematical matrix operations
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#ifndef MATHMATRIX_H
#define MATHMATRIX_H

#include<cmath>
#include<iomanip>
#include<iostream>
#include<vector>
#include<cstdlib>
using namespace std;

const double Zero = 1.0E-8;

//---------------------------------------------------------------------------
//矩阵类
class MathMatrix
{
public:
	//-----------------------------------------------
	//数据成员
	vector<vector<double> > element;

	//-----------------------------------------------
	//成员函数
	//构造函数
	MathMatrix(){};
	MathMatrix(const int&, const int&);							//构造大小为M*N，值为零的矩阵
	MathMatrix(const vector<double>&);						//一维向量构造
	MathMatrix(const vector<vector<double> >&);			//二维向量构造					
	MathMatrix(const double *, const int&);					//一维数组构造
	MathMatrix(const double *, const int&, const int&);	//二维数组构造，实参为（&vec[0][0],m,n）
	//-----------------------------------------------
	void Evalu(const double *, const int&);						//一维数组赋给矩阵 
	void Evalu(const double *, const int& , const int&);	//二维数组赋给矩阵，实参为（&vec[0][0],m,n） 
	void operator=(const vector<double>&);					//重载赋值操作符（一维向量赋给矩阵）
	void operator=(const vector<vector<double> >&);	//重载赋值操作符（二维向量赋给矩阵）
	MathMatrix operator*(const MathMatrix&);				//重载乘法操作符（矩阵乘矩阵）
	MathMatrix operator*(const double&);						//重载乘法操作符（矩阵乘实数）
	MathMatrix operator+(const MathMatrix&);				//重载加法操作符（矩阵加矩阵）			
	double operator+(const double&);								//重载加法操作符（1*1矩阵加实数）			
	MathMatrix operator-(const MathMatrix&);				//重载减法操作符（矩阵减矩阵）
	MathMatrix Inverse(void);					//矩阵求逆
	MathMatrix Transpose(void);			//矩阵转置
	int RowN(void);								//求矩阵行数
	int CalN(void);									//求矩阵列数
	int Symmetry(void);							//判定矩阵对称
	//-----------------------------------------------------------------------------
	MathMatrix GetCol(int cn);                  //求出矩阵一列，将其转为矩阵

	//重载输出流操作符
	friend ostream& operator<<(ostream& o, const MathMatrix& matrix);
};
//重载输出流操作符
inline ostream& operator<<(ostream& o,const MathMatrix& matrix)
{
	for (int i=0; i<int(matrix.element.size()); i++)
	{
		for(int j=0; j<int(matrix.element[i].size()); j++)
		{
//			if (fabs(matrix.element[i][j])<1.0E-8)
//			{
//				o <<setw(8) << setfill(' ') << 0 << " ";
//			}
//			else
//			{
				o <<setw(16) << setfill(' ') << matrix.element[i][j] << " ";
//			}
		}
		o << endl;
	}
	return o;
}
//---------------------------------------------------------------------------
#endif
//===========================================================================
