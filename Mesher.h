//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Mesher.h
//OBJECTIVE:	Generate the mesh of RVE, including cracks.
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================
#ifndef MESHER_H
#define MESHER_H

#include "Input_Reader.h"
#include "Geometry_3D.h"
#include"Gauss.h"
#include "Fem_3D.h"
#include "WeightFunc.h"
#include "time.h"
#include "Hns.h"
using namespace hns;

//---------------------------------------------------------------------------
class Mesher 
{
	public:

		//数据成员
		vector<Node> nodes;				//节点类向量
		vector<Element> elements;		//六面体单元类向量
		
		//Constructor
		Mesher(){};
		//Generate mesh with cracks
		int Generate_mesh(const struct Grid_size &grid_size, const struct Ele_prop &ele_prop, const struct Crack &Cracks, const struct Geom_RVE &geom_rve, const string disc, const struct Weight_func &weight_func);
		//Generate meshes for damage and fracture propagation
		int Generate_mesh(const struct Grid_size &grid_size, const struct Ele_prop &ele_prop, const struct Geom_RVE &geom_rve, const string disc, const vector<bool> &full_dam_eles);
		//Generate meshes with cracks for damage and fracture propagation
		int Generate_mesh(const struct Grid_size &grid_size, const struct Ele_prop &ele_prop, const struct Crack &Cracks, const struct Geom_RVE &geom_rve, const string disc, const vector<bool> &full_dam_eles);
		//Import mesh data from a file and reconfiguration mesh data(DGFEM)
		int Import_mesh_reconfiguration(string output_file_name, const string disc, const vector<bool> &full_dam_eles);
		//Import mesh data from a file and reconfiguration mesh data(DGFEM) for 2D
		int Import_mesh_reconfiguration_2D(string output_file_name, const string disc, const vector<bool> &full_dam_eles, const struct Weight_func &weight_func);
		//Import mesh data from a file and reconfiguration mesh (DGFEM) by weighting function for 2D
		int Mesh_reconfiguration_weight_func_2D(string output_file_name, const string disc, const struct Weight_func &weight_func);
		//Import mesh data from a file and reconfiguration mesh data(DGFEM) by weigted element for 2D
		int Mesh_reconfiguration_weighted_eles_2D(string output_file_name, const string disc, const vector<bool> &weighted_eles);
		//Import mesh data from a file and reconfiguration mesh (DGFEM) by weighting function and delta legnth for 2D
		int Mesh_reconfiguration_weight_delta_2D(string output_file_name, const string disc, const struct Weight_func &weight_func, const struct Peri_para &peri_para);
		//Import mesh data from a file and reconfiguration mesh data(DGFEM) for pure nonlocal fracture propagation
		int Import_mesh_reconfiguration(string output_file_name, const string disc, const struct Weight_func &weight_func);


		//生成格子用于预估局部模型等效刚度阵
		int Generate_grids_for_effective_stiffness(const struct Grid_size &nonloc_gsize, const double &decayR, double &grid_vol);
		int Generate_grids_for_effective_stiffness_2D(const struct Grid_size &nonloc_gsize, const double &decayR, double &grid_vol);
		//Find the neighbour nodes of every node amd the neighbour elements of every element (New Version with weight_func determination)
		int Deter_relative_nodes_elements(vector<Node> &nodes, vector<Element> &elements, const double &dist, const string &mod, const struct Crack &cracks, const struct Weight_func &weight_func)const;
		//Find the neighbour elements of every element
		int Deter_relative_elements(vector<Node> &nodes, vector<Element> &elements, const double &dist, const string &mod, const struct Crack &cracks, const struct Weight_func &weight_func)const;
		//Find the neighbour nodes of every node amd the neighbour elements of every element
		int Deter_relative_nodes_elements(vector<Node> &nodes, vector<Element> &elements, const double &dist, const string &mod, const struct Crack &cracks)const;
		//Estimate the volume of brick element
		double Calculate_brick_volume(const Node elenod[])const;
		double Calculate_brick_volume(const double (*elenod)[3])const;
		//Estimate the area of quadrilateral element
		double Calculate_quadri_area(const Node elenod[])const;

	private:
		double dx, dy, dz; //各个方向的剖分细度
		struct Crack_Paras //裂纹的参数
		{	
			int num, n0, n1, n2, n3;
			char ch0, ch1, ch2;
			double d10, d20, d21, d30, d31;
		};		

		//成员函数
		//Generate brick background grids
		int Brick_background_mesh(const struct Grid_size &grid_size, const struct Geom_RVE &geom_rve);
		//Generate brick grids for static fracture problems
		int Brick_background_mesh(const struct Grid_size &grid_size, const struct Geom_RVE &geom_rve, const string &disc, const struct Weight_func &weight_func);
		//Generate brick grids for damage and fracture propagation
		int Brick_background_mesh(const struct Grid_size &grid_size, const struct Geom_RVE &geom_rve, const string &disc, const vector<bool> &full_dam_eles);
		//Read 3D mesh data from the file
		int Read_mesh(const string &output_file_name);
		//Read 2D mesh data from the file
		int Read_mesh_2DTet(const string &output_file_name);
		//Read 2D mesh data from 3D mesh data  file
		int Read_mesh_3To2D(const string &output_file_name, const double &Const_Z);
		//赋予网格单元材料属性
		int Element_material_property(const struct Ele_prop &ele_prop);
		//输出Tecplot可视化网格数据
		int Export_mesh_data_tecplot(string output_file_name)const;
		//Export mesh data file for 2D grids shown in Tecplot
		int Export_mesh_data_tecplot_2D(string output_file_name)const;
		//根据给定的i j k以及最大i_max j_max k_max, 决定给定节点的位置（角点、边界线、边界面、内部）
		int Deter_node_type(const int i, const int j, const int k, const int i_max, const int j_max, const int k_max)const;
		//生成网格裂纹数据
		int Add_crack_to_element(const struct Crack &Cracks, const struct Geom_RVE &geom_rve);
		//the parameters of single crack
		int Single_crack_parameters(const int &nc, const struct Crack &Cracks, const struct Geom_RVE &geom_rve, struct Crack_Paras &crapa);
		//Find the neigbour elements of every node
		int Deter_nodes_relative_eles();
		//根据裂纹调整网格
		int Adjust_mesh_with_crack(const struct Crack_Paras &crapa);
		//Evaluate the relative elements under a crack condition
		void Evaluate_relative_element_cracks(const double &vi1, const double &vi2, const double &vj1, const double &vj2, const double &d10, const double &d20, const double &d21, Element &ele, const int &nj)const;
};
//-------------------------------------------------------
#endif
//=====================================================
