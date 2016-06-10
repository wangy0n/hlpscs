//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Fem_3D.h
//OBJECTIVE:	Basic finite element computation classes, including notes and elements, etc.
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#ifndef FEM_3D_H
#define FEM_3D_H

#include<iostream>
#include<cmath>
#include<vector>
using namespace std;


//----------------------------------------------------------------------------
// class containing the individual elementary matrices
class Matrices
{
public:

	double Matelem[24][24];
	vector<double> Relative_mat;
};

//---------------------------------------------------------------------------
//定义节点类
class Node
{	
	public:
		double x, y, z; //节点的坐标值

		int type;	//0内点，1边界面点，2边界线点，3角点
		bool DG_node;						//the node for Discontinuous Galerkin computation
		vector<int> relative_nods;	//相关节点向量（节点编号）
		vector<int> relative_eles;	//相关单元向量（单元编号）

		//构造函数
		Node(const double ix=0, const double iy=0, const double iz=0);
		Node(const double ix, const double iy, const double iz, const int itype);
		
		//成员函数；
		double	distance_to(const Node& n); //距离
		double	distance_to(Node& n); //距离

};
//---------------------------------------------------------------------------
//单元类头文件；
class Element
{
	public:	
		//Data Member
		int type;	//表示单元的形状类型(三个数字xyz： x 维度；y单元的节点个数；z单元的形函数幂次)；
						//例如，121: 一维两节点(线段)线性形函数；231: 二维三节点(三角形)线性形函数；241: 二维四节点(四边形)线性形函数；
						//341: 三维四节点(四面体)线性形函数；361: 三维六节点(三棱柱)线性形函数；381：三维八节点(六面体)线性形函数
		int mat;	//表示单元的材料物性；
		int flag;   //用于做单元变化的标记；
		bool DG_elem;						//The element for Discontinuous Galerkin computation
		double Damage_energy;		//The damaged energy of this element for continuum damage mechanics
		double Dissipative_energy;	//The dissipative energy of this element for bond broken in peridynamics
		double dam;							//The average damaged value in the center of an element
		vector<int> nodes_id;
		vector<int> relative_eles;	//相关单元向量（单元编号）
		
		//Constructor
		Element(){};
};
//---------------------------------------------------------------------------
//定义六面体单元
class Hexahedron
{
	public:
        int nodesId[8];	 //eight nodes;
        int materialId;
        int facesId[6];

		//构造函数；
		Hexahedron(){};
        Hexahedron(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);
		
		//member function
		//Calculate the volume of a Hexahedron
		double cal_hex_volume(const vector<Node> &nodes)const;
	private:
		//计算3个节点坐标组成的3×3行列式的值
		double cal_det(const Node &n1, const Node &n2, const Node &n3)const ;
		//计算4节点四面体单元的体积(节点顺序满足右手定则）
		double cal_tet_volume(const Node &node1, const Node &node2, const Node &node3, const Node &node4)const;
};
//---------------------------------------------------------------------------
//边线类头文件；
class Side
{
public:	
	int type;	//表示边线的类型；0：内线段；1：力边线；2：自由边线；3：固定边线
	vector<int> nodes_id;
};
#endif
//===========================================================================
