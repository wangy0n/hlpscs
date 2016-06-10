//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	SolveEqu.h
//OBJECTIVE:	Solving linear equations
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================
#ifndef SOVLEEQU_H
#define SOVLEEQU_H

#include<assert.h>
#include<iostream>
#include<cmath>
#include<vector>
#include "Geometry_3D.h"
#include "Geometry_2D.h"
#include "Fem_3D.h"
#include "Hns.h"
#include "Input_Reader.h"
using namespace hns;
#include <omp.h>

#define SCHUNKSIZE 1		//defines the chunk size as 1 contiguous iteration
//---------------------------------------------------------------------------
//定义解方程类
class SolveEqu
{
public:
	//紧缩存储刚度矩阵(紧缩存储节点对应关系)
	int izig(const vector<Node> &nodes, vector<int> &Iz, vector<int> &Ig)const;

	//Set the fixed displacement constraints
	int Fixed_displacement_constraints(const struct Displace &displace, const double &rmp, const vector<Node> &nodes, int &bnod_num, vector<int> &ip, vector<double> &vp)const;
	//Set the fixed displacement constraints for 2D
	int Fixed_displacement_constraints_2D(const struct Displace &displace, const double &rmp, const vector<Node> &nodes, int &bnod_num, vector<int> &ip, vector<double> &vp)const;

	//处理三维固定位移约束点的零位移值
	void Deal_with_displacement_zero_value(const int &bnod_num, const int &nod_size, const vector<int> &Iz, const vector<int> &Ig, 
																			 const vector<int> &ip, const vector<double> &vp, vector <double> &F, vector <double> &AK)const;
	
	//求解线性方程组
	void Solve_linear_equations(const int &bnod_num, const int &N, const vector<int> &Iz, const vector<int> &Ig, const vector<int> &ip, const vector<double> &vp, 
													  const vector<double> &A, const vector<double> &B, vector<double> &X)const;
	void Solve_linear_equations_omp(const int &bnod_num, const int &N, const vector<int> &Iz, const vector<int> &Ig, const vector<int> &ip, const vector<double> &vp, 
																const vector<double> &A, const vector<double> &B, vector<double> &X)const;

	//用完全的刚度矩阵及右端项检测
	void Complete_matrix_equright_testing(const vector<double> &A, const vector <double> &F, const vector<Node> &nodes, const vector<int> &Iz, const vector<int> &Ig)const;

private:
	//处理右端向量
	void deal_with_F(const double E[][3], const int num[], const double dist[3], const int &Nsize, const vector<int> &Iz, const vector<int> &Ig, const vector<int> &Iw, const vector<int>* Ik,
								   vector<double> &F, const vector<double> &AK)const;

	//处理三维固定位移约束点的非零位移值
	void displacement_nonzero_value(const int &Kw, const int &Kg, const int &bnod_num, const vector<int> &Ip, const vector<double> &vp, vector<double> &W, vector<double> &G)const;

	//矩阵AK乘向量U
	void mabvm(const int &N, const int &N1, const vector<int> &Iz, const vector<int> &Ig, const vector<double> &AK, const vector<double> &U, vector<double> &V)const;
	//(OpenMP)矩阵AK乘向量U
	void mabvm_omp(const int &N, const int &N1, const vector<int> &Iz, const vector<int> &Ig, const vector<double> &AK, const vector<double> &U, vector<double> &V)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
