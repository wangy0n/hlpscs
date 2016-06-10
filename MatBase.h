//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	MatBase.h
//OBJECTIVE:	Generate material database
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================
#ifndef MATBASE_H
#define MATBASE_H

#include "Input_Reader.h"
#include "MatPro.h"
#include "Gauss.h"
#include "Mesher.h"
#include "Hns.h"
using namespace hns;

//---------------------------------------------------------------------------
class MatBase
{
	public:

		//数据成员
		vector<MatPro> mats_vec;

		//构造函数
		MatBase(){};

		//成员函数
		//Generate the material base
		int Generate_matbase(const struct Stif_nonloc &stif_nonloc, const struct Grid_size &nonloc_gsize, const int &gaupoi_num, struct Peri_para &peri_para);
		//Generate the material base with damage information 
		int Generate_matbase(const struct Stif_nonloc &stif_nonloc, const struct Grid_size &nonloc_gsize, const int &gaupoi_num, const struct Damage &damages, struct Peri_para &peri_para);
		//Generate the material base with damage information for 2D
		int Generate_matbase_2D(const struct Stif_nonloc &stif_nonloc, const struct Grid_size &nonloc_gsize, const int &gaupoi_num, const struct Damage &damages, struct Peri_para &peri_para);

	private:
		//成员函数
		//计算长程力等效刚度值在三维球坐标系下
		double Stiffness_long_range_interaction_Spherical(const int &ni, const int &nj, const vector<Node> &gauss, const vector<double> &weight, const double &decayR, const double &inst_len);
		//计算长程力等效刚度值在三维笛卡尔坐标系下
		double Stiffness_long_range_interaction_Cartesian(const int &ni, const int &nj, const vector<Node> &gauss, const vector<double> &weight, const double &decayR, const double &inst_len);
		//计算球面坐标向量
		double Spherial_coordinates(const int &num, const double &r, const double &sita, const double &phi);
		//计算笛卡尔坐标向量
		double Cartesian_coordinates(const int &num, const double &x, const double &y, const double &z);
		//输出所有材料的刚度矩阵
		void Print_stiffness_matrix(string print_name, const double acoe[])const;
		void Print_stiffness_matrix(string print_name, const double acoe[], const int &key)const;
		void Print_stiffness_matrix(string print_name, MatPro &mat, const int &key)const;
		//预估局部模型等效刚度阵
		int Estimate_stiffness_long_range_interaction(const vector<Element> &elements, const vector<Node> &nodes, const double &grid_vol, const vector<Node> &gauss, 
													const vector<double> &weight, const double a[], MatPro &mat, const double &decayR, const double &inst_len, const int &output_key)const;
		int Estimate_stiffness_long_range_interaction_2D(const vector<Element> &elements, const vector<Node> &nodes, const double &grid_vol, const vector<Node> &gauss, 
													const vector<double> &weight, const double a[], MatPro &mat, const double &decayR, const double &inst_len, const int &output_key)const;
};
//-------------------------------------------------------
#endif
//===========================================================================
