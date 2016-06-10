//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Geometry_3D.cpp
//OBJECTIVE:	The definitions of point, line and plane in 3D
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#include "Geometry_3D.h"

//---------------------------------------------------------------------------
//三维点类
Point_3D::Point_3D( double px, double py, double pz )
{
	x = px;
	y = py;
	z = pz;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator+( Point_3D &pt )
{
	Point_3D rp( x + pt.x, y + pt.y, z + pt.z );
	return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator+( const Point_3D &pt )
{
	Point_3D rp( x + pt.x, y + pt.y, z + pt.z );
	return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator+( double d )
{
	Point_3D rp( x + d, y + d, z + d );
	return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator-( Point_3D &pt )
{
	Point_3D rp( x - pt.x, y - pt.y, z - pt.z );
	return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator-( const Point_3D &pt )
{
	Point_3D rp( x - pt.x, y - pt.y, z - pt.z );
	return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator-( double d )
{
	Point_3D rp( x - d, y - d, z - d );
	return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator*( double d )
{
	Point_3D rp( x*d, y*d, z*d );
	return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator/( double d )
{
	Point_3D rp( x/d, y/d, z/d );
	return rp;
}
//---------------------------------------------------------------------------
bool Point_3D::operator==( Point_3D &pt )
{
	return (x==pt.x&&y==pt.y&&z==pt.z);
}
//---------------------------------------------------------------------------
bool Point_3D::operator!=( Point_3D &pt )
{
	return (x!=pt.x||y!=pt.y||z!=pt.z);
}
//---------------------------------------------------------------------------
double Point_3D::distance_to(const Point_3D &pt )const
{
	double rv2 = (x-pt.x)*(x-pt.x)+(y-pt.y)*(y-pt.y)+(z-pt.z)*(z-pt.z);
	return sqrt(rv2);
}
//---------------------------------------------------------------------------
double Point_3D::distance_to(const double &px, const double &py, const double &pz)const
{
	double rv2 = (x-px)*(x-px)+(y-py)*(y-py)+(z-pz)*(z-pz);
	return sqrt(rv2);
}
//===========================================================================

//三维线类
//---------------------------------------------------------------------------
//构造函数
Line_3D::Line_3D(Point_3D p0, Point_3D p1)
{
	point[0] = p0;
	point[1] = p1;
	xm = p1.x - p0.x;
	yn = p1.y - p0.y;
	zl = p1.z - p0.z;
	len = length();
	if(len==0) virtual_line = false;
	else virtual_line = true;
}
//---------------------------------------------------------------------------   
//线段的长度
double Line_3D::length() 
{
	double dx = point[1].x-point[0].x;
	double dy = point[1].y-point[0].y;
	double dz = point[1].z-point[0].z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}
//---------------------------------------------------------------------------
//点到直线的距离
double Line_3D::distance_point_to_line(const Point_3D *point_temp)const
{
	if(xm==0&&yn==0&&zl==0) hout << "注意！该线段退化成一个点！" <<endl;
	double X = yn*(point_temp->y-point[0].y)-zl*(point_temp->z-point[0].z);
	double Y = zl*(point_temp->z-point[0].z)-xm*(point_temp->x-point[0].x);
	double Z = xm*(point_temp->x-point[0].x)-yn*(point_temp->y-point[0].y);
	return sqrt(X*X+Y*Y+Z*Z)/sqrt(xm*xm+yn*yn+zl*zl);
}
double Line_3D::distance_point_to_line(const Point_3D &point_temp)const
{
	if(xm==0&&yn==0&&zl==0) hout << "注意！该线段退化成一个点！" <<endl;
	double X = yn*(point_temp.y-point[0].y)-zl*(point_temp.z-point[0].z);
	double Y = zl*(point_temp.z-point[0].z)-xm*(point_temp.x-point[0].x);
	double Z = xm*(point_temp.x-point[0].x)-yn*(point_temp.y-point[0].y);
	return sqrt(X*X+Y*Y+Z*Z)/sqrt(xm*xm+yn*yn+zl*zl);
}
double Line_3D::distance_point_to_line(const double dx, const double dy, const double dz)const
{
	if(xm==0&&yn==0&&zl==0) hout << "注意！该线段退化成一个点！" <<endl;
	double X = yn*(dy-point[0].y)-zl*(dz-point[0].z);
	double Y = zl*(dz-point[0].z)-xm*(dx-point[0].x);
	double Z = xm*(dx-point[0].x)-yn*(dy-point[0].y);
	return sqrt(X*X+Y*Y+Z*Z)/sqrt(xm*xm+yn*yn+zl*zl);
}
//---------------------------------------------------------------------------
//判断线段包含一个点
int Line_3D::contain(const Point_3D &point_temp)const
{
	//点到线段两个端点的距离大于两个端点之间的距离
	if( fabs(point_temp.distance_to(point[0])+point_temp.distance_to(point[1])-point[0].distance_to(point[1]))>Zero ) return 0; 
	return 1;
}
//===========================================================================

//三维面类
//---------------------------------------------------------------------------
//构造函数
Plane_3D::Plane_3D(double para[])
{
	for(int i=0; i<4; i++)	coef[i] = para[i];
	if(coef[0]==0.0&&coef[1]==0.0&&coef[2]==0.0) virtual_plane = false;
	else virtual_plane = true;
}
//---------------------------------------------------------------------------
Plane_3D::Plane_3D(vector<double> &para)
{
	for(int i=0; i<4; i++)	coef[i] = para[i];
	if(coef[0]==0.0&&coef[1]==0.0&&coef[2]==0.0) virtual_plane = false;
	else virtual_plane = true;
}
//---------------------------------------------------------------------------
Plane_3D::Plane_3D(const vector<double> &para)
{
	for(int i=0; i<4; i++)	coef[i] = para[i];
	if(coef[0]==0.0&&coef[1]==0.0&&coef[2]==0.0) virtual_plane = false;
	else virtual_plane = true;
}
//---------------------------------------------------------------------------
//判断空间平面包含一个点
int Plane_3D::contain(const Point_3D &point_temp)const
{
	if( coef[0]*point_temp.x+coef[1]*point_temp.y+coef[2]*point_temp.z+coef[3]==0 ) return 1;	//在平面上
	return 0;	//在平面外
}
//---------------------------------------------------------------------------
//判断空间平面包含一个点
int Plane_3D::contain(const double dx, const double dy, const double dz)const
{
	if( coef[0]*dx+coef[1]*dy+coef[2]*dz+coef[3]==0 ) return 1; //在平面上
	return 0;  //在平面外
}
//===========================================================================
