//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Geometry_2D.cpp
//OBJECTIVE:	The definitions of point, line and shape in 2D
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#include "Geometry_2D.h"

//---------------------------------------------------------------------------
//二维点类
Point_2D::Point_2D( double px, double py )
{
	x = px;
	y = py;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator+( Point_2D &pt )
{
	Point_2D rp( x + pt.x, y + pt.y);
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator+( const Point_2D &pt )
{
	Point_2D rp( x + pt.x, y + pt.y);
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator+( double d )
{
	Point_2D rp( x + d, y + d);
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator-( Point_2D &pt )
{
	Point_2D rp( x - pt.x, y - pt.y );
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator-( const Point_2D &pt )
{
	Point_2D rp( x - pt.x, y - pt.y );
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator-( double d )
{
	Point_2D rp( x - d, y - d );
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator*( double d )
{
	Point_2D rp( x*d, y*d );
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator/( double d )
{
	Point_2D rp( x/d, y/d );
	return rp;
}
//---------------------------------------------------------------------------
bool Point_2D::operator==( Point_2D &pt )
{
	return (x==pt.x&&y==pt.y);
}
//---------------------------------------------------------------------------
bool Point_2D::operator!=( Point_2D &pt )
{
	return (x!=pt.x||y!=pt.y);
}
//---------------------------------------------------------------------------
double Point_2D::distance_to(const Point_2D &pt )const
{
	double rv2 = (x-pt.x)*(x-pt.x)+(y-pt.y)*(y-pt.y);
	return sqrt(rv2);
}
//---------------------------------------------------------------------------
double Point_2D::distance_to(const double &px, const double &py)const
{
	double rv2 = (x-px)*(x-px)+(y-py)*(y-py);
	return sqrt(rv2);
}
//===========================================================================

//二维线类
//---------------------------------------------------------------------------
//构造函数
Line_2D::Line_2D(Point_2D p0, Point_2D p1)
{
	point[0] = p0;
	point[1] = p1;
	A = p1.y - p0.y ;
	B = -(p1.x - p0.x);
	D = -(A*p0.x+B*p0.y);
	len = length();
	if(len==0) virtual_line = false;
	else virtual_line = true;
}
//---------------------------------------------------------------------------   
//线段的长度
double Line_2D::length() 
{
	double dx = point[1].x-point[0].x;
	double dy = point[1].y-point[0].y;
	return sqrt(dx*dx+dy*dy);
}
//---------------------------------------------------------------------------
//点到直线的距离
double Line_2D::distance_point_to_line(const Point_2D *point_temp)const
{
	if(A==0&&B==0) hout << "注意！该线段退化成一个点！求点到直线的距离时，分母为零！" <<endl;
	return fabs(A*point_temp->x+B*point_temp->y+D)/sqrt(A*A+B*B);
}
double Line_2D::distance_point_to_line(const Point_2D &point_temp)const
{
	if(A==0&&B==0) hout << "注意！该线段退化成一个点！求点到直线的距离时，分母为零！" <<endl;
	return fabs(A*point_temp.x+B*point_temp.y+D)/sqrt(A*A+B*B);
}
double Line_2D::distance_point_to_line(double dx, double dy)const
{
	if(A==0&&B==0) hout << "注意！该线段退化成一个点！求点到直线的距离时，分母为零！" <<endl;
	return fabs(A*dx+B*dy+D)/sqrt(A*A+B*B);
}
//---------------------------------------------------------------------------
//点到直线的路径，判断点是在直线的一侧还是另一侧
double Line_2D::path_point_to_line(const Point_2D &point_temp)const
{
	if(A==0&&B==0) hout << "注意！该线段退化成一个点！求点到直线的距离时，分母为零！" <<endl;
	return A*point_temp.x+B*point_temp.y+D;
}
//---------------------------------------------------------------------------
//判断线段包含一个点
int Line_2D::contain(const Point_2D &point_temp)const
{
	//点到线段两个端点的距离大于两个端点之间的距离
	if( fabs(point_temp.distance_to(point[0])+point_temp.distance_to(point[1])-point[0].distance_to(point[1]))>Zero ) return 0; 
	return 1;
}
int Line_2D::contain(const double &px, const double &py)const
{
	Point_2D point_temp(px, py);
	//点到线段两个端点的距离大于两个端点之间的距离
	if( fabs(point_temp.distance_to(point[0])+point_temp.distance_to(point[1])-point[0].distance_to(point[1]))>Zero ) return 0; 
	return 1;
}
//===========================================================================

//矩形类
//---------------------------------------------------------------------------
//构造函数
Rectangle::Rectangle(Point_2D p0, Point_2D p1, Point_2D p2, Point_2D p3)
{
	point[0] = p0;
	point[1] = p1;
	point[2] = p2;
	point[3] = p3;
	length = p2.x-p0.x;
	width = p2.y-p0.y;
	area = length*width;
	if(length>0&&width>0) virtual_rect = 1;
	else virtual_rect = 0;
}
Rectangle::Rectangle(Point_2D p0, double len, double wid)
{
	point[0].x = p0.x;		point[0].y = p0.y;
	point[1].x = p0.x;		point[1].y = p0.y+wid;
	point[2].x = p0.x+len;		point[2].y = p0.y+wid;
	point[3].x = p0.x+len;		point[3].y = p0.y;
	length = len;
	width = wid;
	area = length*width;
	if(length>0&&width>0) virtual_rect = 1;
	else virtual_rect = 0;
}
Rectangle::Rectangle(Point_2D p0, Point_2D p2)
{
	point[0] = p0;
	point[1].x = p0.x;
	point[1].y = p2.y;
	point[2] = p2;
	point[3].x = p2.x;
	point[3].y = p0.y;
	length = p2.x-p0.x;
	width = p2.y-p0.y;
	area = length*width;
	if(length>0&&width>0) virtual_rect = 1;
	else virtual_rect = 0;
}
//---------------------------------------------------------------------------
//计算矩形的面积
double Rectangle::calculate_area()
{
	return length*width;
}
//---------------------------------------------------------------------------
//判断矩形区域是否包含了一个点（返回值：0没包含，1包含）
int Rectangle::contain(const Point_2D &poi)const
{
	if(  poi.x>=point[0].x && poi.x<=point[2].x &&
		  poi.y>=point[0].y && poi.y<=point[2].y ) return 1;
	else return 0;
}
//---------------------------------------------------------------------------
//判断一个点是否被包含在矩形区域内部(边界除外)（返回值：0没包含，1包含）
int Rectangle::contain_in(const Point_2D &poi)const
{
	if(  poi.x>point[0].x && poi.x<point[2].x &&
		  poi.y>point[0].y && poi.y<point[2].y ) return 1;
	else return 0;
}
//---------------------------------------------------------------------------
//判断矩形区域是否包含了一段线段（返回值：0没包含，1包含）
int Rectangle::contain(const Line_2D &line)const
{
	for(int i=0; i<2; i++)
	if( contain(line.point[i])==0 ) return 0;
	return 1;
}
//---------------------------------------------------------------------------
//判断矩形区域是否包含了另一矩形（返回值：0没包含，1包含）
int Rectangle::contain(const Rectangle &rect)const
{
	for(int i=0; i<4; i++)
	if( contain(rect.point[i])==0 ) return 0;
	return 1;
}
//---------------------------------------------------------------------------
//判断一段线段是否被包含在矩形内部，在矩形边界上不算（返回值：0没包含，1包含） 
int Rectangle::contain_in(const Line_2D &line)const
{
	if( contain(line.point[0])==1&&contain(line.point[1])==1 ) 
	{
		if( line.point[0].x == line.point[1].x&&(line.point[0].x==point[0].x||line.point[0].x==point[2].x) ) return 0;
		if( line.point[0].y == line.point[1].y&&(line.point[0].y==point[0].y||line.point[0].y==point[2].y) ) return 0;
		return 1;
	}
	return 0;
}
//---------------------------------------------------------------------------
//判断线段是否和矩形相交，相交包括线段的两个端点都在矩形的边界上，或至少一个端点在矩形内
//但如果仅有一个端点在矩形边界上，另一端点在矩形外部则不算相交)(返回值：0不相交，1相交)
int Rectangle::overlap(const Line_2D &line)const
{
	if( contain(line.point[0])==1||contain(line.point[1])==1 )
	{
		if( line.point[0].x == line.point[1].x) 
		{
			if( line.point[0].y > line.point[1].y&&(line.point[0].y <= point[0].y||line.point[1].y>=point[2].y) ) return 0;
			if( line.point[0].y < line.point[1].y&&(line.point[1].y <= point[0].y||line.point[0].y>=point[2].y) ) return 0;
		}
		if( line.point[0].y == line.point[1].y) 
		{
			if( line.point[0].x > line.point[1].x&&(line.point[0].x <= point[0].x||line.point[1].x>=point[2].x) ) return 0;
			if( line.point[0].x < line.point[1].x&&(line.point[1].x <= point[0].x||line.point[0].x>=point[2].x) ) return 0;
		}
		return 1;
	}
	return 0;
}
//---------------------------------------------------------------------------
//计算矩形Rect1和Rect2的重叠矩形区
int Rectangle::rectangles_overlapping(const Rectangle &Rect1, const Rectangle &Rect2)
{
	//两个矩形的公共部分，即两个矩形的最低点与最高点坐标依次排序，选中间的两个值
	//两个矩形的最低点中取坐标最大值
	if(Rect1.point[0].x<=Rect2.point[0].x)	point[0].x = Rect2.point[0].x;
	else point[0].x = Rect1.point[0].x;
	if(Rect1.point[0].y<=Rect2.point[0].y)	point[0].y = Rect2.point[0].y;
	else point[0].y = Rect1.point[0].y;
	//两个矩形的最高点中取坐标最小值
	if(Rect1.point[2].x<=Rect2.point[2].x)	point[2].x = Rect1.point[2].x;
	else point[2].x = Rect2.point[2].x;
	if(Rect1.point[2].y<=Rect2.point[2].y)	point[2].y = Rect1.point[2].y;
	else point[2].y = Rect2.point[2].y;
	
	point[1].x = point[0].x;
	point[1].y = point[2].y;
	point[3].x = point[2].x;
	point[3].y = point[0].y;
	length = point[2].x-point[0].x;
	width = point[2].y-point[0].y;
	area = length*width;
	if(length>0&&width>0) virtual_rect = 1;
	else virtual_rect = 0;  //没有重叠矩形，或者仅仅重叠一点或一条线段

	return virtual_rect;
}
//---------------------------------------------------------------------------
//将矩形以内部小矩形为中心做九宫格分割
int Rectangle::make_nine_grids_by(const Rectangle rect, vector<Rectangle> &grids)
{
	if( contain(rect)==0 ) 
	{
		hout << "该矩形不在将要被分割的大矩形区域内！ 请检查！" << endl;
		return 0;
	}	
	//记录基本的分割点
	Point_2D point_base[4];
	point_base[0] = point[0];
	point_base[1] = rect.point[0];
	point_base[2] = rect.point[2];
	point_base[3] = point[2];
	//计算九宫格上所有的节点（从左向右，从下向上）
	Point_2D point_ex[16];
	for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
		{
			point_ex[i*4+j].x = point_base[j].x;
			point_ex[i*4+j].y = point_base[i].y;
		}
	//记录区域的各个子矩形
	grids.clear();  //情况向量
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
		{
			Rectangle rect_temp( point_ex[i*4+j], point_ex[(i+1)*4+j], point_ex[(i+1)*4+j+1], point_ex[i*4+j+1] );
			grids.push_back(rect_temp);
		}

	return 1;
}

//---------------------------------------------------------------------------
//输出矩形的全部信息（包括四个顶点坐标，长，宽，面积，虚矩形标记）
void Rectangle::output_parameters()
{
	hout << "The coordinates of four corner points (starting from lower-left corner point in clockwise direction):" << endl;
	for(int i=0; i<4; i++) hout << "(" << point[i].x << ", " << point[i].y << ")" << "   ";
	hout << endl;
	hout << "The length, width, area and virtual_rect tag:" << endl;
	hout << "Length " << length << ", Width " << width << ", area " << area << ", virtual_rect " << virtual_rect << endl;
}
//===========================================================================
