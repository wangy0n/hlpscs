//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Fem_3D.cpp
//OBJECTIVE:	Basic finite element computation classes, including notes and elements, etc.
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#include "Fem_3D.h"

//---------------------------------------------------------------------------
//Node类构造函数；
Node::Node(const double ix, const double iy, const double iz)
{
	x=ix;
	y=iy;
	z=iz;
}//---------------------------------------
Node::Node(const double ix, const double iy, const double iz, const int itype)
{
	x=ix;
	y=iy;
	z=iz;
	type=itype;
}//---------------------------------------
//---------------------------------------------------------------------------
//距离；
double Node::distance_to(const Node& n)
{
	return(sqrt( (x-n.x)*(x-n.x)+(y-n.y)*(y-n.y)+(z-n.z)*(z-n.z)) );
}
double Node::distance_to(Node& n)
{
	return(sqrt( (x-n.x)*(x-n.x)+(y-n.y)*(y-n.y)+(z-n.z)*(z-n.z)) );
}//----------------------------------------

//=================================================================
//Element 构造/成员函数；

//=================================================================
//Hexahedron 构造/成员函数；
//---------------------------------------------------------------------------
//六面体单元构造函数
Hexahedron::Hexahedron(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{
	nodesId[0] = n1; nodesId[1] = n2;
	nodesId[2] = n3; nodesId[3] = n4;
	nodesId[4] = n5; nodesId[5] = n6;
	nodesId[6] = n7; nodesId[7] = n8;
}
//---------------------------------------------------------------------------
//Calculate the volume of a Hexahedron
double Hexahedron::cal_hex_volume(const vector<Node> &nodes)const
{
	double volume = cal_tet_volume(nodes[nodesId[0]], nodes[nodesId[2]], nodes[nodesId[1]], nodes[nodesId[4]])
							  + cal_tet_volume(nodes[nodesId[1]], nodes[nodesId[4]], nodes[nodesId[2]], nodes[nodesId[5]])
							  + cal_tet_volume(nodes[nodesId[1]], nodes[nodesId[2]], nodes[nodesId[3]], nodes[nodesId[5]])
							  + cal_tet_volume(nodes[nodesId[2]], nodes[nodesId[4]], nodes[nodesId[6]], nodes[nodesId[7]])
							  + cal_tet_volume(nodes[nodesId[2]], nodes[nodesId[3]], nodes[nodesId[5]], nodes[nodesId[7]])
							  + cal_tet_volume(nodes[nodesId[2]], nodes[nodesId[5]], nodes[nodesId[4]], nodes[nodesId[7]]);
	return volume;
}
//---------------------------------------------------------------------------
//计算3个节点坐标组成的3×3行列式的值
double Hexahedron::cal_det(const Node &n1, const Node &n2, const Node &n3)const
{
	double det	=	n1.x*(n2.y*n3.z-n2.z*n3.y)+
							n1.y*(n2.z*n3.x-n2.x*n3.z)+
							n1.z*(n2.x*n3.y-n2.y*n3.x);
	return det ;
}
//---------------------------------------------------------------------------
//计算4节点四面体单元的体积(节点顺序满足右手定则）
double Hexahedron::cal_tet_volume(const Node &node1, const Node &node2, const Node &node3, const Node &node4)const 
{
	double volume;
	double aai = cal_det( node2, node3, node4 );
	double aaj = cal_det( node3, node4, node1 );
	double aak = cal_det( node4, node1, node2 );
	double aal = cal_det( node1, node2, node3 );

	volume = 1.0/6.0*(aai-aaj+aak-aal) ;
	return volume ;
}
//=================================================================
