//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Geometry_2D.h
//OBJECTIVE:	The definitions of point, line and shape in 2D
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#ifndef GEOMETRY_2D_H
#define GEOMETRY_2D_H

#include<cmath>
#include<stdlib.h>
#include<vector>
#include "Hns.h"
using namespace hns;

#include "MathMatrix.h"

//---------------------------------------------------------------------------
//定义二维空间点类
class Point_2D
{
	public:

		//数据成员
		double x, y;
		int flag;

		//构造函数
		Point_2D(){};
		Point_2D( double px, double py );

		//成员函数
        Point_2D operator+( Point_2D &pt );
        Point_2D operator+( const Point_2D &pt );
        Point_2D operator+( double d );
        Point_2D operator-( double d );
        Point_2D operator-( Point_2D &pt );
        Point_2D operator-( const Point_2D &pt );
        Point_2D operator*( double d );
        Point_2D operator/( double d );
		bool operator==( Point_2D &pt );
		bool operator!=( Point_2D &pt );
        double distance_to(const Point_2D &pt)const;
		double distance_to(const double &px, const double &py)const;
};

//---------------------------------------------------------------------------
//定义二维线(段)类
class Line_2D
{
	public:

		//数据成员
		 Point_2D point[2];	//表示线段的两个端点坐标
		 double len;					//表示线段的长度
		 bool virtual_line;			//表示线段是否被认定为虚线段(fase:虚线段（收缩为一点）; true:非虚线段)


		//构造函数
		Line_2D(){};
		Line_2D(Point_2D p0, Point_2D p1);
		
		//成员函数
		double length();	//线段的长度
        double distance_point_to_line(const Point_2D *point_temp)const;    //点到直线的距离
        double distance_point_to_line(const Point_2D &point_temp)const;    //点到直线的距离
        double distance_point_to_line(double dx, double dy)const;    //点到直线的距离
        double path_point_to_line(const Point_2D &point_temp)const;    //点到直线的路径，判断点是在直线的一侧还是另一侧
		int contain(const Point_2D &point_temp)const;    //判断线段包含一个点
		int contain(const double &px, const double &py)const;    //判断线段包含一个点

	private:
				
		//数据成员
        double A, B, D ;       //直线方程的系数 Ax+By+D=0
};

//---------------------------------------------------------------------------
//定义矩形类
class Rectangle
{
	public:

		//数据成员
		Point_2D point[4];	//定义矩阵的四个角点，从左下角开始顺时针转
		double length;			//定义矩阵的长
		double width;				//定义矩阵的宽
		double area;				//表示矩形的面积
		int virtual_rect;			//表示矩形是否被认定为虚矩形(0:虚矩形; 1:非虚矩形)

		//构造函数
		Rectangle(){};
		Rectangle (Point_2D p0, Point_2D p1, Point_2D p2, Point_2D p3);
		Rectangle (Point_2D p0, double len, double wid);
		Rectangle (Point_2D p0, Point_2D p2);

		//成员函数
		double calculate_area();			//计算矩形的面积
		int contain(const Point_2D &poi)const;		//判断矩形区域是否包含了一个点（返回值：0没包含，1包含）
		int contain_in(const Point_2D &poi)const;//判断一个点是否被包含在矩形区域内部(边界除外)（返回值：0没包含，1包含）
		int contain(const Line_2D &line)const;			//判断矩形区域是否包含了一段线段（返回值：0没包含，1包含）
		int contain(const Rectangle &rect)const;		//判断矩形区域是否包含了另一矩形（返回值：0没包含，1包含）
		int contain_in(const Line_2D &line)const;	//判断一段线段是否被包含在矩形内部，在举行边界上不算（返回值：0没包含，1包含） 
		int overlap(const Line_2D &line)const;			//判断线段是否和矩形相交，相交包括线段的两个端点都在矩形的边界上，或至少一个端点在矩形内
																					//但如果仅有一个端点在矩形边界上，另一端点在矩形外部则不算相交)(返回值：0不相交，1相交)
		void output_parameters();		//输出矩形的全部信息（包括四个顶点坐标，长，宽，面积，虚矩形标记）
		int rectangles_overlapping(const Rectangle &Rect1, const Rectangle &Rect2);	 //计算矩形Rect1和Rect2的重叠矩形区
		int make_nine_grids_by(const Rectangle rect, vector<Rectangle> &grids);	//将矩形以内部小矩形为中心做九宫格分割
};
//---------------------------------------------------------------------------
//记录耦合区单元与区域1或2网格的相交单元和相交小矩形
struct Relative_Ele_Rect
{	
	vector<int> relative_ele_num[2];
	vector<Rectangle> relative_rect[2];
};
#endif
//===========================================================================
