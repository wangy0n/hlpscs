//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	MatPro.h
//OBJECTIVE:	Material properties
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================

#ifndef MATPRO_H
#define MATPRO_H

#include "Input_Reader.h"
#include "Gauss.h"
#include "MathMatrix.h"
#include "Hns.h"

using namespace hns;

//---------------------------------------------------------------------------
class MatPro
{
	public:	
	
		//数据成员
		string type_val;   //三维材料性质类型，目前定义三种("Isotropic"：各向同性; "Transverse"：横观各向同性; "Orthotropic"：正交各向异性)
		double E11, E22, E33, Nu12, Nu23, Nu13, G12, G23, G13;	//材料的弹性模量，泊松比，剪切模量
		double elas_matrix[6][6];		//Dij 弹性矩阵

		//构造函数；
		MatPro(string itype="Isotropic"){ type_val = itype; }
		MatPro(const struct Stif_nonloc &stif_nonloc);

		//---------------------------------------------
		//成员函数；
		//设置各向同性材料参数的弹性模量、泊松比和剪切模量
		void set_ela_para(const double &E);
		//设置横观各向同性材料的弹性模量、泊松比和剪切模量
		void set_ela_para(const double &E1, const double &E2, const double &Nu2);
		//设置正交各向异性材料的弹性模量、泊松比和剪切模量
		void set_ela_para(const double &iE1, const double &iE2, const double &iE3, const double &iNu12, const double &iNu23, const double &iNu13);
		//生成Dij弹性矩阵，此弹性矩阵对应的应变向量为（e11,e22,e33,2*e12,2*e23,2*e31）转置
		int Generate_elas_matrix();
		int Generate_elas_matrix_2D();
		//根据Dij弹性矩阵，反求材料参数
		int Get_ele_para_by_ela_matrix();
		//根据理论公式反求系数，并验证一致性（由于要除以Horizon半径，所以Horizon半径应该大于1，避免误差过大）
		void Compare_coef_by_analysis_formula(const int mat_type, const double &decayR, const double acoe[])const;
		//根据理论公式反求刚度矩阵，并验证一致性
		void Compare_matrix_by_analysis_formula(const double &decayR, const double decay_acoe[])const;

		void print()const;
};
//------------------------------------------------
#endif
//===========================================================================
