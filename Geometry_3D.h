//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Geometry_3D.h
//OBJECTIVE:	The definitions of point, line and plane in 3D
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#ifndef GEOMETRY_3D_H
#define GEOMETRY_3D_H

#include<cmath>
#include<stdlib.h>
#include<vector>
#include "Hns.h"
using namespace hns;

#include "MathMatrix.h"
//---------------------------------------------------------------------------
//定义三维空间点类
class Point_3D
{
	public:

		//数据成员
		double x, y, z;
		int flag;

		//构造函数
		Point_3D(){};
		Point_3D( double px, double py, double pz );

		//成员函数
        Point_3D operator+( Point_3D &pt );
        Point_3D operator+( const Point_3D &pt );
        Point_3D operator+( double d );
        Point_3D operator-( double d );
        Point_3D operator-( Point_3D &pt );
        Point_3D operator-( const Point_3D &pt );
        Point_3D operator*( double d );
        Point_3D operator/( double d );
		bool operator==( Point_3D &pt );
		bool operator!=( Point_3D &pt );
        double distance_to(const Point_3D &pt)const;
		double distance_to(const double &px, const double &py,  const double &pz)const;
};

//---------------------------------------------------------------------------
//定义三维线(段)类
class Line_3D
{
	public:

		//数据成员
		 Point_3D point[2];	//表示线段的两个端点坐标
		 double len;					//表示线段的长度
		 bool virtual_line;		//表示线段是否被认定为虚线段(false:虚线段（收缩为一点）; true:非虚线段)


		//构造函数
		Line_3D(Point_3D p0, Point_3D p1);
		
		//成员函数
		double length();	//线段的长度
        double distance_point_to_line(const Point_3D *point_temp)const;    //点到直线的距离
        double distance_point_to_line(const Point_3D &point_temp)const;    //点到直线的距离
        double distance_point_to_line(const double dx, const double dy, const double dz)const;  //点到直线的距离
		int contain(const Point_3D &point_temp)const;    //判断线段包含一个点

	private:
				
		//数据成员
        double xm, yn, zl ;       //直线方程的比例系数 (x-x0)/xm=(y-y0)/yn=(z-z0)/zl
};
//---------------------------------------------------------------------------
//定义三维面类
class Plane_3D
{
	public:

		//数据成员
		double coef[4];			//表示空间平面方程的四个系数ax+by+cz+d=0
		bool virtual_plane;		//表示平面是否被认定为虚平面(false:虚平面（法向向量为(0,0,0)）; true:非虚平面)

		//构造函数
		Plane_3D(){};
		Plane_3D(double para[]);
		Plane_3D(vector<double> &para);
		Plane_3D(const vector<double> &para);
		
		//成员函数
        double distance_point_to_plane(const Point_3D *point_temp)const;    //点到空间平面的距离
        double distance_point_to_plane(const Point_3D &point_temp)const;    //点到空间平面的距离
        double distance_point_to_line(const double dx, const double dy, const double dz)const;  //点到空间平面的距离
		int contain(const Point_3D &point_temp)const;    //判断空间平面包含一个点
		int contain(const double dx, const double dy, const double dz)const;    //判断空间平面包含一个点
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
