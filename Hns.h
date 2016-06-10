//===========================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Hns.h (Fei Han's namespace definition)
//OBJECTIVE:	Output stream definition for outputing into a file
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa;  harrisonhanfei@gmail.com
//===========================================================================

#ifndef HNS_H
#define HNS_H

#include<fstream>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<string>
#include<vector>
#include<cmath>
#include "MathMatrix.h"
using namespace std;

//---------------------------------------------------------------------------
class SetWP
{
public:
	SetWP& operator()(int w)
	{
		width = w;
		precision = 0;
		return *this;
	}
	SetWP& operator()(int w,int p)
	{
		width = w;
		precision = p;
		return *this;
	}
	int width;
	int precision;
};
//---------------------------------------------------------------------------
class foutstream 
{
public:
	//构造函数
	foutstream();

	//操作符定义
	foutstream& operator << (short n);
	foutstream& operator << (unsigned short n);
	foutstream& operator << (int n);
	foutstream& operator << (unsigned int n);
//	foutstream& operator << (size_t n);
	foutstream& operator << (long n);
	foutstream& operator << (unsigned long n);
	foutstream& operator << (long long n);
	foutstream& operator << (unsigned long long n);
	foutstream& operator << (float f);
	foutstream& operator << (double f);
	foutstream& operator << (long double f);
	foutstream& operator << (char c);
	foutstream& operator << (bool b);
	foutstream& operator << (string s);
	foutstream& operator << (foutstream& (*os)(foutstream& ys));
	foutstream& operator << (SetWP& sw);
	foutstream& operator << (const char* c);
	foutstream& operator << (const MathMatrix& matrix);

	int width(int n);
	int precision(int n);
	int setf(int n);
	int setf(int n1,int n2); 

//	void open_deffo_stream(char* file,int mod=ios::app);	//app表示续写
	void open_deffo_stream(char* file,int mod=ios::out);	//out是表示重新写
	void close_deffo_stream( );

	string output_file;
private:
	ostream *out_stream; 
	ofstream of_stream;
	bool is_open_file;
};
//---------------------------------------------------------------------------
namespace hns
{
	extern foutstream hout;
//	void open_deffo_stream(char* file,int mod=ios::app);	//app表示续写
	void open_deffo_stream(char* file,int mod=ios::out);	//out是表示重新写

	void close_deffo_stream( );
	foutstream& endl(foutstream& os);
	extern SetWP setwp;
}

//---------------------------------------------------------------------------
#endif
//===========================================================================
