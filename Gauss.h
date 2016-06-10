//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Gauss.h
//OBJECTIVE:	Compute the gaussian points
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#ifndef GAUSS_H
#define GAUSS_H

#include<iostream>
#include<vector>
#include<cmath>
#include "Hns.h"
#include "Fem_3D.h"
using namespace hns;

const double PI = 3.1415926535897932384626433832795;

//----------------------------------------------------------

class Gauss
{
	public:
		//数据变量
		vector<Node> gauss; 
		vector<double> weight;

		//构造函数；
		Gauss(){};
	
		//成员函数
		//生成三维高斯节点
		int Generate_gauss(const int &ipre);
		//生成二维高斯节点
		int Generate_gauss_2D(const int &ipre);

	private:

		//数据变量；
		int precision;

		//用于生成高斯点序列
		int Generate_gauss_array(vector<double> &gauss, vector<double> &weight)const;
		//读入信息一行，跳过注释行（以%开头）
		string Get_Line(ifstream &infile)const;
};//-----------------------------------------------------
#endif
//===========================================================================
