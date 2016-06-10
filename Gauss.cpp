//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Gauss.cpp
//OBJECTIVE:	Compute the gaussian points
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#include "Gauss.h"

//---------------------------------------------------------------------------
//生成三维高斯节点
int Gauss::Generate_gauss(const int &ipre)
{
	//Define the number of gaussian points
	precision = ipre;
	//---------------------------------------------------------------------------
	//生成高斯节点
	vector<double> temp_gauss;
	vector<double> temp_wight;
	if(Generate_gauss_array(temp_gauss, temp_wight)==0) return 0;

	for(int i=0; i<(int)temp_gauss.size(); i++)
	    for(int j=0; j<(int)temp_gauss.size(); j++)
			for(int k=0; k<(int)temp_gauss.size(); k++)
			{
				weight.push_back(temp_wight[i]*temp_wight[j]*temp_wight[k]);
				gauss.push_back(Node(temp_gauss[i], temp_gauss[j], temp_gauss[k]));
			}

	return 1;
}
//---------------------------------------------------------------------------
//生成二维高斯节点
int Gauss::Generate_gauss_2D(const int &ipre)
{
	//Define the number of gaussian points
	precision = ipre;
	//---------------------------------------------------------------------------
	//生成高斯节点
	vector<double> temp_gauss;
	vector<double> temp_wight;
	if(Generate_gauss_array(temp_gauss, temp_wight)==0) return 0;

	for(int i=0; i<(int)temp_gauss.size(); i++)
	    for(int j=0; j<(int)temp_gauss.size(); j++)
		{
			weight.push_back(temp_wight[i]*temp_wight[j]);
			gauss.push_back(Node(temp_gauss[i], temp_gauss[j], 0));
		}

	return 1;
}
//---------------------------------------------------------------------------
//用于生成高斯点序列
int Gauss::Generate_gauss_array(vector<double> &gauss, vector<double> &weight)const
{
	//开始计算
	int n = precision;
	int iter = 10;   //迭代次数
	int m=int((n+1)/2); 
	int e1=n*(n+1);
	int mm=4*m-1;
	vector<long double> t;
	for(int i=3; i<=mm; i=i+4)
	{
		t.push_back((PI/(4*n+2))*i);
	}
	int nn=(1-(1-1/n)/(8*n*n));
	vector<long double> xo;
	for(int i=0; i<m; i++)
	{
		xo.push_back(nn*cos(t[i]));
	}
	vector<long double> pk;
	vector<long double> den(m);
	vector<long double> d1(m);
	vector<long double> dpn(m);
	vector<long double> d2pn(m);
	vector<long double> d3pn(m);
	vector<long double> d4pn(m);
	vector<long double> u(m);
	vector<long double> v(m);
	vector<long double> h(m);
	vector<long double> p(m);
	vector<long double> dp(m);

	for(int j=1; j<=iter; j++)
	{
		vector<long double> pkm1(m, 1.0);
		pk=xo;
		for(int k=2; k<=n; k++)
		{
			vector<long double> t1(m);
			vector<long double> pkp1(m);
			for(int i=0; i<m; i++)
			{
				t1[i]=xo[i]*pk[i];
				pkp1[i]=t1[i]-pkm1[i]-(t1[i]-pkm1[i])/k+t1[i];
			}
			pkm1=pk;
			pk=pkp1;
		}
		for(int i=0; i<m; i++)
		{
			den[i]=1.0-xo[i]*xo[i];
			d1[i]=n*(pkm1[i]-xo[i]*pk[i]);
			dpn[i]=d1[i]/den[i];
			d2pn[i]=(2.0*xo[i]*dpn[i]-e1*pk[i])/den[i];
			d3pn[i]=(4*xo[i]*d2pn[i]+(2-e1)*dpn[i])/den[i];
			d4pn[i]=(6*xo[i]*d3pn[i]+(6-e1)*d2pn[i])/den[i];
			u[i]=pk[i]/dpn[i]; 
			v[i]=d2pn[i]/dpn[i];
			h[i]=-u[i]*(1+(0.5*u[i])*(v[i]+u[i]*(v[i]*v[i]-u[i]*d3pn[i]/(3*dpn[i]))));
			p[i]=pk[i]+h[i]*(dpn[i]+(0.5*h[i])*(d2pn[i]+(h[i]/3)*(d3pn[i]+0.25*h[i]*d4pn[i])));
			dp[i]=dpn[i]+h[i]*(d2pn[i]+(0.5*h[i])*(d3pn[i]+h[i]*d4pn[i]/3));
			h[i]=h[i]-p[i]/dp[i]; 
			xo[i]=xo[i]+h[i];
		}
	}
	vector<long double> bp(m);
	vector<long double> fx(m);
	vector<long double> wf(m);
	for(int i=0; i<m; i++)
	{
		bp[i]=-xo[i]-h[i];
		fx[i]=d1[i]-h[i]*e1*(pk[i]+(h[i]/2)*(dpn[i]+(h[i]/3)*(d2pn[i]+(h[i]/4)*(d3pn[i]+(0.2*h[i])*d4pn[i]))));
		wf[i]=2*(1-bp[i]*bp[i])/(fx[i]*fx[i]);
	}

	if((m+m)>n) bp[m-1]=0;
	if((m+m) != n) m=m-1;
	for(int i=m-1; i>=0; i--)
	{
		bp.push_back(-bp[i]);
		wf.push_back(wf[i]);
	}
 	for(int i=0; i<n; i++)
	{
		gauss.push_back(bp[i]);
		weight.push_back(wf[i]);
	}

	//cout << "bp=" << endl;
	//for(int i=0; i<n; i++)
	//{
	//	cout << setprecision(20) << bp[i] << endl;
	//}
	//cout << "wf=" << endl;
	//long double sum=0.0;
	//for(int i=0; i<n; i++)
	//{
	//	sum += wf[i];
	//	cout << setprecision(20) << wf[i] << endl;
	//}
	//cout << setprecision(20) << "sum=" << sum << endl;

	return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//读入一行信息，并跳过注释行（以"%"开头）；
string Gauss::Get_Line(ifstream &infile)const
{
	string s;
	//读入信息一行
	getline(infile,s);
	//跳过注释行     
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}
//===========================================================================
