//===========================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Hns.cpp (Fei Han's namespace definition)
//OBJECTIVE:	Output stream definition for outputing into a file
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa;  harrisonhanfei@gmail.com
//===========================================================================

#include "Hns.h"
//---------------------------------------------------------------------------
//名空间成员定义
foutstream hns::hout;
SetWP hns::setwp;
//---------------------------------------------------------------------------
void hns::open_deffo_stream(char* file,int mod)
{
	hout.open_deffo_stream(file,mod);
}
//---------------------------------------------------------------------------
void hns::close_deffo_stream( )
{
	hout.close_deffo_stream( );
}                                                                       
//---------------------------------------------------------------------------
foutstream& hns::endl(foutstream& os)
{
	return (os << "\n");
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//文件输出流(foutstream)类成员定义
foutstream::foutstream()
{          
	is_open_file = false;
	out_stream   = &cout;
}
//---------------------------------------------------------------------------
void foutstream::open_deffo_stream(char* file,int mod)
{
	if( is_open_file ) close_deffo_stream();
	of_stream.open(file);
	if( !of_stream )
	{
		cout << "Failed to open file!" << file << std::endl;
		return ;
	}
	out_stream = &of_stream;
	output_file = string(file);  //foutstream public part
	is_open_file = true;
}
//---------------------------------------------------------------------------
void foutstream::close_deffo_stream( )
{
	of_stream.close();
	out_stream   = &cout;
	is_open_file = false;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( short n )
{
	(*out_stream) << n ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( unsigned short n )
{
	(*out_stream) << n ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( int n )
{
	(*out_stream) << n ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( unsigned int n )
{
	(*out_stream) << n ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
//foutstream& foutstream:: operator <<( size_t n )
//{
//	(*out_stream) << n ;
//	out_stream->flush();
//	return *this;
//}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( long n )
{
	(*out_stream) << n ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( unsigned long n )
{
	(*out_stream) << n ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( long long n )
{
	(*out_stream) << n ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( unsigned long long n )
{
	(*out_stream) << n ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( float f )
{
	(*out_stream) << f ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( double f )
{
	(*out_stream) << f ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( long double f )
{
	(*out_stream) << f ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( char c )
{
	(*out_stream) << c ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( bool b )
{
	(*out_stream) << b ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( string s )
{
	(*out_stream) << s ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( foutstream& (*os)(foutstream& yt) )
{
	(*os)(*this);
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<(const char* c )
{
	(*out_stream) << c ;
	out_stream->flush();
	return *this;
}
foutstream& foutstream:: operator << (const MathMatrix& matrix)
{
	for (int i=0; i<int(matrix.element.size()); i++)
	{
		for(int j=0; j<int(matrix.element[i].size()); j++)
		{
//			if (fabs(matrix.element[i][j])<1.0E-8)
//			{
//				(*out_stream) <<setw(8) << setfill(' ') << 0 << " ";
//			}
//			else
//			{
				(*out_stream) <<setw(16) << setfill(' ') << matrix.element[i][j] << " ";
//			}
		}
		(*out_stream) << endl;
	}
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator << (SetWP& sw)
{
	this->width(sw.width);
	if( sw.precision != 0 )
	{
		this->precision(sw.precision);
	}
	return *this;
}
//---------------------------------------------------------------------------
int foutstream:: width ( int n )
{
	return int(out_stream->width(n));
}
//---------------------------------------------------------------------------
int foutstream:: precision( int n )
{
	return int(out_stream->precision(n));
}
//---------------------------------------------------------------------------
int foutstream:: setf( int n )
{
//	return out_stream->setf(n);
	return 1;
}
//---------------------------------------------------------------------------
int foutstream:: setf( int n1, int n2 )
{
//	return out_stream->setf(n1,n2);
	return 1;
}
//===========================================================================
