//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Postprocessor.h
//OBJECTIVE:	Implement the postprocessor (the visual output)
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================

#ifndef POSTPROCESSOR_H
#define POSTPROCESSOR_H

#include<iomanip>
#include<iostream>
#include<fstream>
#include<sstream>
#include"Fem_3D.h"
#include "MatBase.h"
#include "Geometry_3D.h"
#include "Global_Stiff_Matrix.h"
#include "time.h"
#include "Hns.h"
using namespace hns;

//---------------------------------------------------------------------------
class Postprocessor 
{
	public:
		//构造函数
		Postprocessor(){};
		Postprocessor(const vector<double> &post_u, const vector<Node> &post_nodes, const vector<Element> &post_elements);
		Postprocessor(const vector<double> &post_u, const vector<double> &post_dissi, const vector<Node> &post_nodes, const vector<Element> &post_elements);

		//成员函数
		int Treatment(const vector<vector<Node> > &nodes_tot, const vector<vector<Element> > &elements_tot, const vector<MatPro> &mats,
							    const vector<vector<double> > &U_tot, const string &wr_mod, const vector<vector<double> > &Dissipation);
		int Treatment(Input *Init, const vector<vector<Node> > &nodes_tot, const vector<vector<Element> > &elements_tot, const vector<MatPro> &mats,
							    const vector<vector<double> > &U_tot, const string &wr_mod, const vector<vector<double> > &Dissipation);
		int Treatment_2D(Input *Init, const vector<vector<Node> > &nodes_tot, const vector<vector<Element> > &elements_tot, const vector<MatPro> &mats,
										const vector<vector<double> > &U_tot, const string &wr_mod);
		//Immediately write result data
		int Immedi_write_data(const int &Unum, const string &wr_mod);
		//Immediately write result data with weight functions
		int Immedi_write_data(const int &Unum, const string &wr_mod, const struct Weight_func &weight_func);
		//Export deformed mesh for testing
		int Testing_export_deformed_mesh(const string &output_file_name)const;

	private:

		//数据变量
		vector<double> u;
		vector<Node> nodes;
		vector<Element> elements;
		vector<double> dissi;
		vector<double> total_ener;
		vector<vector<double> > ele_strain;
		vector<vector<double> > nod_strain;
		vector<double> nod_local_dam;			//the local damage value of every node
		vector<double> nod_nonlocal_dam;		//the nonlocal damage value of every node
		vector<double> nod_damage;				//damaged energy of every node
		vector<double> nod_dissip;					//bond broken energy of every node
		vector<double> nod_dissip_ratio;			//the ration of bond broken energy of every node (%)
		vector<vector<double> > ele_stress;
		vector<vector<double> > nod_stress;
		vector<double> ele_energy;
		vector<double> nod_energy;
		vector<double> total_break_energy;

		Weight_func postp_weight_func;

		//成员函数
		//输出或读取文件数和加权函数
		int wr_files_num_weight_funcs(const string &output_file_name, const string &wr_mod, const struct Weight_func &weight_func, int &Unum);
		//输出或读取位移解
		int write_or_read_data(const string &output_file_name, const string &wr_mod);
		//计算单元形心处应变向量
		int Calculate_Elements_Strain(const vector<Node> &nodes, const vector<Element> &elements);
		int Calculate_Elements_Strain_2D(const vector<Node> &nodes, const vector<Element> &elements);
		//计算节点处应变向量
		int Calculate_Nodes_Strain(vector<Node> &nodes, const vector<Element> &elements);
		//计算节点处应变向量、应力向量和应变能密度
		int Calculate_Nodes_Strain_Stress_Energy(string mod, double d_crit, vector<Node> &nodes, const vector<Element> &elements);
		//输出Tecplot可视化位移场云图
		int Export_Displacement_Contour(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const;
		int Export_Displacement_Contour_2D(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const;
		//输出Tecplot可视化边界面位移场云图
		int Export_Boundary_Displacement_Contour(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements, const int key)const;
		//输出网格变形图
		int Export_Deformed_Mesh(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const;
		//输出应变场云图
		int Export_Strain_Mesh(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const;
		//输出应变场云图 for 2D problem
		int Export_Strain_Mesh_2D(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const;
		//输出应变场、应力场和应变能密度云图 for 2D problem
		int Export_Strain_Stress_Energy_Mesh_2D(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const;
		//坐标变换并输出应变场云图
		int Export_Strain_After_Coordinates_Transformation(const vector<Node> &nodes, const vector<Element> &elements)const;
		//Export total dissipatve energy in tecplot
		int Export_Total_Dissipative_Energy(const string &output_file_name, const vector<vector<double> > &Dissipation, const string &wr_mod);
		//Export 2D dissipatve energy in tecplot
		int Export_Dissipative_Energy_2D(const string &output_file_name, const vector<double> &damage_ener, const vector<double> &break_ener)const;
		//向量叉积
		Point_3D Vector_Products(const Point_3D &p1, const Point_3D &p2)const;
		//法向向量旋转至Z轴变换矩阵
		void Coordinates_Transformation_Matrix(const Point_3D &poi, MathMatrix &transf)const;
		//Estimate the average stress and strain at the boundary layer elements
		void Estimate_boundary_stress_strain(const vector<Node> &nodes, const vector<Element> &elements, const vector<MatPro> &mats, const int &inum, const struct Geom_RVE &geom_rve,
																   const struct Iterative &iter, const struct Displace &displace, vector<double> &ave_strain_xx, vector<double> &ave_stress_xx)const;
		int Estimate_boundary_stress_strain(const vector<Node> &nodes, const vector<Element> &elements, const vector<MatPro> &mats, const int &inum, const struct Geom_RVE &geom_rve,
															   vector<double> &ave_strain_xx, vector<double> &ave_stress_xx)const;
		//Estimate the average force at the top center of beam under three point bending
		int Estimate_boundary_force_2D(const vector<Node> &nodes, const vector<Element> &elements, const vector<MatPro> &mats, const int &inum, const struct Force_Disp_TPB &force_disp, vector<double> &TPB_force_yy)const;
		//Export stress strain curve
		int Export_Stress_Strain_Curve(const string &output_file_name, const vector<double> &ave_strain_xx, const vector<double> &ave_stress_xx)const;
		//Export boundary force displacement curve for 2D three point bending
		int Export_Force_Displacement_Curve(const string &output_file_name, const vector<double> &TPB_force_yy, double &delta_disp)const;
		//Find the same nodes seperated by DG elements 
		int Find_same_DG_nodes(const int &nnum, const vector<Node> &nodes, vector<int> &same_nods)const;
		//计算单元形心处应力向量和应变能密度for 2D problem (Plane Strain Assumption)
		int Calculate_Elements_Stress_Energy_2D(const vector<MatPro> &mats, const vector<Element> &elements);

};
//---------------------------------------------------------------------------
#endif
//===========================================================================
