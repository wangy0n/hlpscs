//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	SolveEqu.cpp
//OBJECTIVE:	Solving linear equations
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#include "SolveEqu.h"

//---------------------------------------------------------------------------
//紧缩存储刚度矩阵
int SolveEqu::izig(const vector<Node> &nodes, vector<int> &Iz, vector<int> &Ig)const
{
	//iz[i-1], iz[i]中记录与第i个点相关的数据在ig中存储的开始位置和结束位置
	//第一个点的编号最小(没有比它编号还小的节点)，所以iz数据从第二个点(i==1)开始
	Iz.push_back(0);

	//从第二个节点开始循环所有的节点
	for(int i=1; i<(int)nodes.size(); i++)
	{
		for(int j=0; j<(int)nodes[i].relative_nods.size(); j++)	//注意这里节点的相关节点是从小到大排列的，并且不包括节点自身的编号
		{
			if(nodes[i].relative_nods[j]>=i) break;						//一旦一个编号大了（等于的情况其实不存在），说明后面的编号都大了
			else Ig.push_back(nodes[i].relative_nods[j]);
		}
		Iz.push_back((int)Ig.size());
	}

	//输出用于检查
	//hout << "Iz:" << endl;
	//for(int i=0; i<(int)Iz.size(); i++)
	//	hout << i << " " << Iz[i] << endl;
	//hout << endl << "Ig:" << endl;
	//for(int i=0; i<(int)Ig.size(); i++)
	//	hout << i << " " << Ig[i] << endl;

	return 1;
}
//-----------------------------------------------------------
//处理右端向量F
void SolveEqu::deal_with_F(const double E[][3], const int num[], const double dist[3], const int &Nsize, const vector<int> &Iz, const vector<int> &Ig, const vector<int> &Iw, const vector<int>* Ik,
											   vector<double> &F, const vector<double> &AK)const
{
	//-----------------------------------------------------------
	//计算向量Q (xj=alpha*xi - Q; alpha =1《有限元方法及其应用》), dist是xi到xj的向量, 所以udiffi相当于Q
	double udiff[3] = {0, 0, 0};
	for(int j=0; j<3; j++)
		for(int k=0; k<3; k++)
			udiff[j] += E[j][k]*dist[k];

	//-----------------------------------------------------------
	//节点编号k<节点编号j
	int Ih=0;
	if(num[0]!=0)
	{
		for(int j=Iz[num[0]-1]; j<Iz[num[0]]; j++)
		{
			Ih = 6*num[0] + 9*j;
			for(int k=0; k<=2; k++)
				for(int m=0; m<=2; m++)
					F[3*Ig[j]+k] += AK[Ih+k+3*m]*udiff[m];  //注意是乘矩阵的转置
		}
	}
	//-----------------------------------------------------------
	//节点编号k>节点编号j
	if(num[0]!=Nsize)
	{
		for(int j=Iw[num[0]]; j<Iw[num[0]+1]; j++)
		{
			Ih = 6*Ik[0][j] + 9*Ik[1][j];
			for(int k=0; k<=2; k++)
				for(int m=0; m<=2; m++)
					F[3*Ik[0][j]+k] += AK[Ih+3*k+m]*udiff[m];  //注意是乘矩阵本身
		}
	}

	//-----------------------------------------------------------
	//节点编号i
	
	//二分法在i的相关节点中查找j的位置，如果不相关矩阵就是零矩阵	
	bool mark = false;
	int left = Iz[num[1]-1];
	int middle = 0;
	int right = Iz[num[1]]-1;
	while(right>=left)
	{
		middle = (left + right)/2;
		if(Ig[middle] == num[0]) { mark = true; break; }				//找到相关位置
		else if(Ig[middle] > num[0]) right = middle - 1;
		else left = middle + 1;
	}

	//计算j点的对角阵在AK中的起始位置
	Ih = 6*num[0] + 9*Iz[num[0]];
	if(mark) //i与j相关
	{
		int Ij = 6*num[1] + 9*middle;
		for(int j=0; j<=2; j++)
		{
			F[3*num[1]+j] += F[3*num[0]+j];
			for(int k=0; k<=2; k++)
			{
				if(j==2||k==2) 
					F[3*num[1]+j] += (AK[Ij+3*j+k]+AK[Ih+j+k+1])*udiff[k];
				else
					F[3*num[1]+j] += (AK[Ij+3*j+k]+AK[Ih+j+k])*udiff[k];
			}
		}
	}
	else //i与j不相关, AK[ij]=0
	{
		for(int j=0; j<=2; j++) 
		{
			F[3*num[1]+j] += F[3*num[0]+j];
			for(int k=0; k<=2; k++)
			{
				if(j==2||k==2) 
					F[3*num[1]+j] += AK[Ih+j+k+1]*udiff[k];
				else
					F[3*num[1]+j] += AK[Ih+j+k]*udiff[k];
			}
		}
	}

	//-----------------------------------------------------------
	//节点编号j
	for(int j=0; j<=2; j++) F[3*num[0]+j] = 0;

}
//---------------------------------------------------------------------------
//Set the fixed displacement constraints
int SolveEqu::Fixed_displacement_constraints(const struct Displace &displace, const double &rmp, const vector<Node> &nodes, int &bnod_num, vector<int> &ip, vector<double> &vp)const
{
	for(int i=0; i<displace.num; i++)
	{
		//-----------------------------------------------------------
		//Define variables
		Point_3D fpoi;		//Point
		Plane_3D fpla;		//Surface
		struct
		{
			double xmin;
			double xmax;
			double ymin;
			double ymax;
			double zmin;
			double zmax;
		}zone;						//Zone
		zone.xmin = 0.0;
		zone.xmax = 0.0;
		zone.ymin = 0.0;
		zone.ymax = 0.0;
		zone.zmin = 0.0;
		zone.zmax = 0.0;

		//-----------------------------------------------------------
		//Define the displacement domain
		if(displace.domain_type[i]=="Point") 
		{
			fpoi.x = displace.coef[i][0];
			fpoi.y = displace.coef[i][1];
			fpoi.z = displace.coef[i][2];
		}
		else if(displace.domain_type[i]=="Surface") 
		{
			for(int j=0; j<4; j++) fpla.coef[j]	=	displace.coef[i][j];
		}
		else if(displace.domain_type[i]=="Zone")
		{
			zone.xmin = displace.coef[i][0];
			zone.xmax = displace.coef[i][1];
			zone.ymin = displace.coef[i][2];
			zone.ymax = displace.coef[i][3];
			zone.zmin = displace.coef[i][4];
			zone.zmax = displace.coef[i][5];
		}
		else if(displace.domain_type[i]!="All_surfaces")
		{
			cout << "Error: Failed to find the type of displacement domain!" << endl; 
			hout << "Error: Failed to find the type of displacement domain!" << endl;
			return 0; 
		}

		//-----------------------------------------------------------
		//Define the sign of displacement constraint
		int sign = 0;
		for(int j=0; j<(int)displace.disp_type[i].size(); j++) 
		{
			if(displace.disp_type[i][j]=="Fixed_displacement_x") sign += 1;
			else if(displace.disp_type[i][j]=="Fixed_displacement_y") sign += 2;
			else if(displace.disp_type[i][j]=="Fixed_displacement_z") sign += 4;
			else if(displace.disp_type[i][j]!="Pure_shear")
			{
				cout << "Error: Failed to find the type of displacement constraint!" << endl; 
				hout << "Error: Failed to find the type of displacement constraint!" << endl;
				return 0; 
			}
		}

		//-----------------------------------------------------------
		//Define the value of displacement constraint by the ratio
		vector<double> value;
		for(int j=0; j<(int)displace.value[i].size(); j++)  value.push_back(displace.value[i][j]*rmp);
		
		//-----------------------------------------------------------
		int count = 0;
		for(int j=0; j<(int)nodes.size(); j++)
		{
			if((displace.domain_type[i]=="Point"&&fpoi.distance_to(nodes[j].x, nodes[j].y, nodes[j].z)<=Zero)||
				(displace.domain_type[i]=="Surface"&&fpla.contain(nodes[j].x, nodes[j].y, nodes[j].z))||
				(displace.domain_type[i]=="Zone"&&nodes[j].x-zone.xmin>=Zero&&zone.xmax-nodes[j].x>=Zero&&
				  nodes[j].y-zone.ymin>=Zero&&zone.ymax-nodes[j].y>=Zero&&nodes[j].z-zone.zmin>=Zero&&zone.zmax-nodes[j].z>=Zero))
			{
				count ++;	//Find a point
				ip.push_back(j);
				ip.push_back(sign);
				//Initialize the displacement constraints
				//--------------------------------------------------------------------------------------------------------
				//Type	u		v		w		空表示这一维是自由边界(在此处统一赋零值，但是不会被调用)  
				//   7		0		0		0		其它0处表示是零值或非零值 但必须有值
				//	 6		空	0		0		
				//   5		0		空	0    
				//	 4		空	空	0
				//   3		0		0		空
				//   2		空	0		空
				//   1		0		空	空
				//--------------------------------------------------------------------------------------------------------		
				switch(sign)
				{
					case 1:	vp.push_back(value[0]);	vp.push_back(0);				vp.push_back(0);	break;
					case 2:	vp.push_back(0);				vp.push_back(value[0]);	vp.push_back(0);	break;
					case 3:	vp.push_back(value[0]);	vp.push_back(value[1]);	vp.push_back(0);	break;
					case 4:	vp.push_back(0);				vp.push_back(0);				vp.push_back(value[0]);	break;
					case 5:	vp.push_back(value[0]);	vp.push_back(0);				vp.push_back(value[1]);	break;
					case 6:	vp.push_back(0);				vp.push_back(value[0]);	vp.push_back(value[1]);	break;
					case 7:	vp.push_back(value[0]);	vp.push_back(value[1]);	vp.push_back(value[2]);	break;
					default:	
					{
						hout << "Error: the signs of nodes conformed to displacement constraint!" << endl;						
						hout << "Error: the signs of nodes conformed to displacement constraint!" << endl;	
						return 0;
					}
				}
				bnod_num++;  //the number of nodes on the boundary
			}
			else if(displace.domain_type[i]=="All_surfaces"&&displace.disp_type[i][0]=="Pure_shear"&&nodes[j].type!=0) //Points on the sufaces
			{
				count ++;	//Find a point
				
				const double ux = (nodes[j].y+nodes[j].z)*value[0];		//相当于	[	0,		val,	val	]		[x]
				const double uy = (nodes[j].x+nodes[j].z)*value[0];		//				[	val,	0,		val	]	*	[y]		
				const double uz = (nodes[j].x+nodes[j].y)*value[0];		//				[	val,	val,	0		]		[z]

				ip.push_back(j);		ip.push_back(7);
				vp.push_back(ux);		vp.push_back(uy);		vp.push_back(uz); 

				bnod_num++;  //the number of nodes on the boundary
			}
		}
		if(count==0) 
		{ 
			cout << "Error: Failed to find any nodes conformed to the displacement constraint!" << endl; 
			hout << "Error: Failed to find any nodes conformed to the displacement constraint!" << endl; 		
			return 0; 
		}
	}

	return 1;
}
//---------------------------------------------------------------------------
//处理三维固定位移约束点的零位移值
void SolveEqu::Deal_with_displacement_zero_value(const int &bnod_num, const int &nod_size, const vector<int> &Iz, const vector<int> &Ig, 
																						   const vector<int> &ip, const vector<double> &vp, vector <double> &F, vector <double> &AK)const
//--------------------------------------------------------------------------------------------------------
//Type	u		v		w		空表示这一维是自由边界(在此处统一赋零值，但是不会被调用)  
//   7		0		0		0		其它0处表示是零值或非零值 但必须有值
//	 6		空	0		0		
//   5		0		空	0    
//	 4		空	空	0
//   3		0		0		空
//   2		空	0		空
//   1		0		空	空
//--------------------------------------------------------------------------------------------------------		
{
	for(int i=0; i<bnod_num; i++)
	{
		int key[3] = {1, 2, 4}; 
		int Ia = 2*i;
		int Ih = 3*i;
		int Mt = ip[Ia];
		int Lkt = ip[Ia+1];
		for(int j=2; j>=0; j--)	 //以下一小段循环体是要在key[0]~key[2]中分别记录u,v,w哪项有值(key==1)而哪项为空(key==0)
		{
			int item = Lkt%key[j];	//取余数
			key[j] = Lkt/key[j];		//取商
			Lkt = item;
		}

		int Kmt = 0;		//记录初始位置
		//---------------------------------------------------------------------
		for(int j=0; j<=2; j++)
			if(key[j]==1&&fabs(vp[Ih+j])<Zero)  //小于一个极小量就是等于0
			{
				//---------------------------------------------------------------------
				F[3*Mt+j]=0.0;		//右端项
				//---------------------------------------------------------------------
				if(Mt!=0)				//Mt的相关节点中比Mt编号小的
				{
					for(int k=Iz[Mt-1]; k<Iz[Mt]; k++)
						for(int m=0; m<=2; m++)
						{
							Kmt = 6*Mt+9*k+3*j+m;
							AK[Kmt]=0.0;
						}
				}
				//---------------------------------------------------------------------
				for(int k=0; k<=2; k++)  //计算对角矩阵
				{
					if(j==2||k==2) Kmt=6*Mt+9*Iz[Mt]+j+k+1;				//相当于	[0,1,	3]		当j==0时, 消 0,	1,	3; 
					else Kmt=6*Mt+9*Iz[Mt]+j+k;										//				[1,2,	4]		当j==1时, 消	1,	2,	4;
 					if(k==j) AK[Kmt]=1.0;														//				[3,4,	5]		当 j==2时 消 3,	4,	5;
					else AK[Kmt]=0.0;
				}
			}
		//---------------------------------------------------------------------
		for(int k=Mt+1; k<nod_size; k++)  //Mt的相关节点中比Mt编号大的
		{
			//二分法在k的相关节点中查找Mt的位置
			bool mark = false;
			int left = Iz[k-1];
			int middle = 0;
			int right = Iz[k]-1;
			while(right>=left)
			{
				middle = (left + right)/2;
				if(Ig[middle] == Mt) { mark = true; break; }				//找到相关位置
				else if(Ig[middle] > Mt) right = middle - 1;
				else left = middle + 1;
			}		
			if(mark)
			{
				for(int j=0; j<=2; j++)
					if(key[j]==1&&fabs(vp[Ih+j])<Zero)
					{
						for(int m=0; m<=6; m+=3) //二分法
						{
							Kmt = 6*k+9*middle+j+m;
							AK[Kmt]=0.0;
						}
					}
			}
		}
	}
}
//---------------------------------------------------------------------------
//求解线性方程组函数(共轭梯度CONJUGATE GRADIENT METHOD)
void SolveEqu::Solve_linear_equations(const int &bnod_num, const int &N, const vector<int> &Iz, const vector<int> &Ig, const vector<int> &ip, 
																    const vector<double> &vp, const vector<double> &A, const vector<double> &B, vector<double> &X)const
{
	//---------------------------------------------------------------------
	vector<double> P(3*N, 0);
	vector<double> R(3*N, 0);
	vector<double> S(3*N, 0);
	vector<double> V(3*N, 0);
	double Rk,R0,RR0,APP,AK,RRk,BK;
	int N1=3*N;  //节点总自由度                             
	int Nup=N1; //最大迭代步数
	int K=0;
	//---------------------------------------------------------------------
	if(bnod_num!=0) displacement_nonzero_value(1,1,bnod_num,ip,vp,X,V);  //处理位移非零值约束条件

	//---------------------------------------------------------------------
	mabvm(N,N1,Iz,Ig,A,X,V);		//CALCULATE PRODUCT A*X0=> AP  

	//---------------------------------------------------------------------
	for(int i=0; i<N1; i++)		//CALCULATE R0=P0=B-A*X0
	{
		R[i]=B[i]-V[i];
		P[i]=R[i];
	}

	//---------------------------------------------------------------------
	if(bnod_num!= 0) displacement_nonzero_value(0,1,bnod_num,ip,vp,P,R);		//处理位移非零值约束条件

	//---------------------------------------------------------------------
	RR0=0.0;
	for(int i=0; i<N1; i++)			//COMPUTE R(K)*R(K)
	{
		RR0=RR0+R[i]*R[i];
	}
	R0=sqrt(RR0);
	//--------------------------------------------------------------------- 
	//ENTER TO ITERATE                                                
	//---------------------------------------------------------------------
	while(K<Nup)
	{
		K=K+1;
		//---------------------------------------------------------------------
		mabvm(N,N1,Iz,Ig,A,P,S);		//CALCULATE AP=A*P=>S

   		//---------------------------------------------------------------------
		if(bnod_num!=0) displacement_nonzero_value(0,0,bnod_num,ip,vp,S,P);	//处理位移非零值约束条件

		APP=0.0;
		for(int i=0; i<N1; i++)	//CALCULATE INNER PRODUCT (AP,P)                          
		{
			APP=APP+S[i]*P[i];
		}                                    
		AK=-RR0/APP;                                     
		//---------------------------------------------------------------------
		//CALCULATE X(K)=X(K-1)-AK*P(K)                  
		//CALCULATE R(K)=R(K-1)+AK*AP(K)                
		RRk=0.0;
		for(int i=0; i<N1; i++)
		{
			X[i]=X[i]-AK*P[i];
			R[i]=R[i]+AK*S[i];
			RRk=RRk+R[i]*R[i];
		}
		Rk=sqrt(RRk);
		if(fabs(Rk)>=Zero*fabs(R0))
		{
			BK=RRk/RR0;
			RR0=RRk;
			for(int i=0; i<N1; i++)
			{
				P[i]=R[i]+BK*P[i];
			}
		}
		else	break;
	}
	if(K==Nup) hout << "注意！共轭梯度法解线性方程组迭代" <<K<< "次仍未收敛，请检查！" << endl;
}
//---------------------------------------------------------------------------
//(OpenMP)求解线性方程组函数(共轭梯度CONJUGATE GRADIENT METHOD)
void SolveEqu::Solve_linear_equations_omp(const int &bnod_num, const int &N, const vector<int> &Iz, const vector<int> &Ig, const vector<int> &ip, 
																				const vector<double> &vp, const vector<double> &A, const vector<double> &B, vector<double> &X)const
{
	//---------------------------------------------------------------------
	vector<double> P(3*N, 0);
	vector<double> R(3*N, 0);
	vector<double> S(3*N, 0);
	vector<double> V(3*N, 0);
	double Rk,R0,RR0,APP,AK,RRk,BK;
	int N1=3*N;  //节点总自由度                             
	int Nup=N1; //最大迭代步数
	int K=0;
	//---------------------------------------------------------------------
	if(bnod_num!=0) displacement_nonzero_value(1,1,bnod_num,ip,vp,X,V);  //处理位移非零值约束条件

	//---------------------------------------------------------------------
	mabvm_omp(N,N1,Iz,Ig,A,X,V);		//CALCULATE PRODUCT A*X0=> AP

	//---------------------------------------------------------------------
	//CALCULATE R0=P0=B-A*X0
	#pragma omp parallel
	{
		vector<int> inum;
		vector<double> tempR;

		#pragma omp for schedule(dynamic, SCHUNKSIZE)
		for(int i=0; i<N1; i++)
		{
			inum.push_back(i);
			tempR.push_back(B[i]-V[i]);
		}

		#pragma omp critical
		for(int i=0; i<(int)inum.size(); i++)
		{
			R[inum[i]]=tempR[i];
			P[inum[i]]=tempR[i];
		}
	}

	//---------------------------------------------------------------------
	if(bnod_num!= 0) displacement_nonzero_value(0,1,bnod_num,ip,vp,P,R);		//处理位移非零值约束条件

	//---------------------------------------------------------------------
	//COMPUTE R(K)*R(K)
	RR0=0.0;
	#pragma omp parallel
	{
		double temp_RR0=0.0;

		#pragma omp for schedule(dynamic, SCHUNKSIZE)
		for(int i=0; i<N1; i++)			
		{
			temp_RR0+=R[i]*R[i];
		}

		#pragma omp critical
		{
			RR0+=temp_RR0;
		}
	}
	R0=sqrt(RR0);
	//--------------------------------------------------------------------- 
	//ENTER TO ITERATE                                                
	//---------------------------------------------------------------------
	while(K<Nup)
	{
		K=K+1;
		//---------------------------------------------------------------------
		mabvm_omp(N,N1,Iz,Ig,A,P,S);		//CALCULATE AP=A*P=>S

   		//---------------------------------------------------------------------
		if(bnod_num!=0) displacement_nonzero_value(0,0,bnod_num,ip,vp,S,P);	//处理位移非零值约束条件

		//---------------------------------------------------------------------
		//CALCULATE INNER PRODUCT (AP,P)
		APP=0.0;
		#pragma omp parallel
		{
			double temp_APP=0.0;

			#pragma omp for schedule(dynamic, SCHUNKSIZE)
			for(int i=0; i<N1; i++)
			{
				temp_APP+=S[i]*P[i];
			}

			#pragma omp critical
			{
				APP+=temp_APP;
			}
		}
		AK=-RR0/APP;

		//---------------------------------------------------------------------
		//CALCULATE X(K)=X(K-1)-AK*P(K)                  
		//CALCULATE R(K)=R(K-1)+AK*AP(K)                
		RRk=0.0;
		#pragma omp parallel
		{
			vector<int> inum;
			vector<double> tempX;
			vector<double> tempR;
			double temp_RRk = 0.0;

			#pragma omp for schedule(dynamic, SCHUNKSIZE)
			for(int i=0; i<N1; i++)
			{
				inum.push_back(i);
				tempX.push_back(X[i]-AK*P[i]);
				tempR.push_back(R[i]+AK*S[i]);
				double Ri = tempR.back();
				temp_RRk += Ri*Ri;
			}

			#pragma omp critical
			{
				for(int i=0; i<(int)inum.size(); i++)
				{
					X[inum[i]] = tempX[i];
					R[inum[i]] = tempR[i];
				}
				RRk+=temp_RRk;
			}
		}
		Rk=sqrt(RRk);

		//---------------------------------------------------------------------
		if(fabs(Rk)>=Zero*fabs(R0))
		{
			BK=RRk/RR0;
			RR0=RRk;

			#pragma omp parallel
			{
				vector<int> inum;
				vector<double> tempP;

				#pragma omp for schedule(dynamic, SCHUNKSIZE)
				for(int i=0; i<N1; i++)
				{
					inum.push_back(i);
					tempP.push_back(R[i]+BK*P[i]);
				}

				#pragma omp critical
				for(int i=0; i<(int)inum.size(); i++)
				{
					P[inum[i]]=tempP[i];
				}
			}
		}
		else	break;
	}
	if(K==Nup) hout << "注意！共轭梯度法解线性方程组迭代" <<K<< "次仍未收敛，请检查！" << endl;
}
//---------------------------------------------------------------------------
//处理三维固定位移约束点的非零位移值
void SolveEqu::displacement_nonzero_value(const int &Kw, const int &Kg, const int &bnod_num, const vector<int> &ip, const vector<double> &vp, vector<double> &W, vector<double> &G)const
//Kw	当Kw!=0, 将固定值赋予W, 无论Kg的值, 否则将零值赋予W;
//Kg		当Kg!=0, 将零值赋予G.
//--------------------------------------------------------------------------------------------------------
//Type	u		v		w		空表示这一维是自由边界(在此处统一赋零值，但是不会被调用)  
//   7		0		0		0		其它0处表示是零值或非零值 但必须有值
//	 6		空	0		0		
//   5		0		空	0    
//	 4		空	空	0
//   3		0		0		空
//   2		空	0		空
//   1		0		空	空
//--------------------------------------------------------------------------------------------------------		
{
	for(int i=0; i<bnod_num; i++)
	{
		int key[3] = {1, 2, 4};
		int Ia = 2*i;
		int Ib = 3*i;
		int Ih = 3*ip[Ia];
		int Lkt = ip[Ia+1];
		for(int j=2; j>=0; j--)	//以下一小段循环体是要在key[0]~key[2]中分别记录u,v,w哪项有值(key==1)而哪项为空(key==0)
		{
			int item = Lkt%key[j];	//取余数
			key[j] = Lkt/key[j];		//取商
			Lkt = item;
		}
		for(int j=0; j<=2; j++)
		{
			if(key[j]==1&&fabs(vp[Ib+j])>=Zero)  //大于一个极小量就是不等于0
			{
				if(Kw!=0) W[Ih+j]=vp[Ib+j];
				else
				{
					if(Kg!=0) G[Ih+j]=0.0;
					W[Ih+j]=0.0;
				}
			}
		}
	}
}
//---------------------------------------------------------------------------
//矩阵AK乘向量U
void SolveEqu::mabvm(const int &N, const int &N1, const vector<int> &Iz, const vector<int> &Ig, const vector<double> &AK, const vector<double> &U, vector<double> &V)const
//AK	刚度矩阵
//U		输入向量
//V		结果向量
{
	for(int i=0; i< N1; i++) V[i]=0.0;		//赋初值

	for(int i=0; i<N; i++)
	{
		int Ik=6*i+9*Iz[i];                                            
		int Kh=3*i;
		for(int j=0; j<=2; j++)
		{
			for(int k=0; k<=2; k++)
			{
				if(j==2||k==2) V[Kh+j]=V[Kh+j]+AK[Ik+j+k+1]*U[Kh+k];
				else	V[Kh+j]=V[Kh+j]+AK[Ik+j+k]*U[Kh+k];                              
			}
		}
		if(i==0) continue;
		for(int m=Iz[i-1]; m<Iz[i]; m++)
		{
			Ik=6*i+9*m;
			int Ia=3*Ig[m];
			for(int j=0; j<=2; j++)
			{
				for(int k=0; k<=2; k++)
				{
					V[Ia+j]=V[Ia+j]+AK[Ik+3*k+j]*U[Kh+k];                         
					V[Kh+j]=V[Kh+j]+AK[Ik+3*j+k]*U[Ia+k];                        
				}
			}
		}
	}
}
//---------------------------------------------------------------------------
//(OpenMP)矩阵AK乘向量U
void SolveEqu::mabvm_omp(const int &N, const int &N1, const vector<int> &Iz, const vector<int> &Ig, const vector<double> &AK, const vector<double> &U, vector<double> &V)const
//AK	刚度矩阵
//U		输入向量
//V		结果向量
{
	for(int i=0; i< N1; i++) V[i]=0.0;		//赋初值

	#pragma omp parallel
	{
		vector<int> inum;
		vector<double> tempV;

		#pragma omp for schedule(dynamic, SCHUNKSIZE)
		for(int i=0; i<N; i++)
		{
			int Ik=6*i+9*Iz[i];                                            
			int Kh=3*i;
			for(int j=0; j<=2; j++)
			{
				double VKh = 0.0;
				for(int k=0; k<=2; k++)
				{
					if(j==2||k==2) VKh += AK[Ik+j+k+1]*U[Kh+k];
					else	VKh +=  AK[Ik+j+k]*U[Kh+k];
				}
				inum.push_back(Kh+j);
				tempV.push_back(VKh);
			}

			if(i==0) continue;

			for(int m=Iz[i-1]; m<Iz[i]; m++)
			{
				Ik=6*i+9*m;
				int Ia=3*Ig[m];
				for(int j=0; j<=2; j++)
				{
					double VIa = 0.0;
					double VKh = 0.0;
					for(int k=0; k<=2; k++)
					{
						VIa+=AK[Ik+3*k+j]*U[Kh+k];                         
						VKh+=AK[Ik+3*j+k]*U[Ia+k];                        
					}

					inum.push_back(Ia+j);
					tempV.push_back(VIa);

					inum.push_back(Kh+j);
					tempV.push_back(VKh);
				}
			}
		}

		#pragma omp critical
		for(int i=0; i<(int)inum.size(); i++)
		{
			V[inum[i]]+=tempV[i];
		}
	}
}
//---------------------------------------------------------------------------
//用完全的刚度矩阵及右端项检测
void SolveEqu::Complete_matrix_equright_testing(const vector<double> &A, const vector <double> &F, const vector<Node> &nodes, const vector<int> &Iz, const vector<int> &Ig)const
{
	//---------------------------------------------------------------------------
	//计算完全的刚度矩阵和右端项
	vector<double> temp_vec(3*(int)nodes.size(),0.0);
	vector<vector<double> > full_matrix(3*(int)nodes.size(), temp_vec);
	for(int i=0; i<(int)Iz.size(); i++)
	{
		//对角矩阵
		int Mt =  6*i + 9*Iz[i];
		for(int j=0; j<=2; j++)
			for(int k=0; k<=2; k++)
			{
				if(j==2||k==2) full_matrix[3*i+j][3*i+k] = A[Mt+j+k+1];
				else full_matrix[3*i+j][3*i+k] = A[Mt+j+k];
			}
		if(i!=0)
		{
			for(int j=Iz[i-1]; j<Iz[i]; j++)
			{
				Mt =  6*i + 9*j;
				for(int k=0; k<=2; k++)
					for(int m=0; m<=2; m++)
					{
						full_matrix[3*i+k][3*Ig[j]+m] = A[Mt+3*k+m];
						full_matrix[3*Ig[j]+m][3*i+k] = A[Mt+3*k+m];
					}
			}
		}
	}
	vector<double> fullright = F;

	//---------------------------------------------------------------------------
	//输出完全的刚度矩阵查看
	ofstream odata;
	odata.open("Testing_complete_matrix_equright.dat");

	odata << "Complete_matrix:" << endl;
	odata << "//-------------------------------------------------------------------------------------------------" << endl;
	odata << setw(8) << setprecision(4)  << " ";
	for(int i=0; i<(int)nodes.size(); i++)
	{
		for(int j=0; j<3; j++)	odata << setw(8) << setprecision(4) << i;
		odata << " | ";
	}
	odata << endl;
	for(int i=0; i<(int)nodes.size(); i++)
	{
		for(int k=0; k<3; k++)
		{
			odata << setw(8) << setprecision(4) << i;
			for(int j=0; j<(int)nodes.size(); j++)
			{
				for(int m=0; m<3; m++)
				{
					if(fabs(full_matrix[3*i+k][3*j+m])<=Zero) odata << setw(8) << setprecision(4) << 0;
					else odata << setw(8) << setprecision(4) << full_matrix[3*i+k][3*j+m];
				}
				odata << " | ";
			}
			odata << endl;
		}
		odata << endl;
	}
	odata << endl << endl;

	//---------------------------------------------------------------------------
	//输出完全的右端项查看
	odata << "Complete_equright:" << endl;
	odata << "//-------------------------------------------------------------------------------------------------" << endl;
	for(int i=0; i<(int)nodes.size(); i++)
	{
		for(int k=0; k<3; k++)
		{
			odata << setw(10) << setprecision(4) << i;
			if(fabs(fullright[3*i+k])<=Zero) odata << setw(15) << setprecision(4) << 0;
			else odata << setw(15) << setprecision(4) << fullright[3*i+k];
			odata << endl;
		}
		odata << endl;
	}
	odata << endl << endl;

}
//---------------------------------------------------------------------------
//Set the fixed displacement constraints for 2D problem
int SolveEqu::Fixed_displacement_constraints_2D(const struct Displace &displace, const double &rmp, const vector<Node> &nodes, int &bnod_num, vector<int> &ip, vector<double> &vp)const
{
	for(int i=0; i<displace.num; i++)
	{
		//-----------------------------------------------------------
		//Define variables
		Point_2D fpoi;		//Point
		Line_2D	  flin;		//Line
		struct
		{
			double xmin;
			double xmax;
			double ymin;
			double ymax;
		}zone;						//Zone
		zone.xmin = 0.0;
		zone.xmax = 0.0;
		zone.ymin = 0.0;
		zone.ymax = 0.0;

		//-----------------------------------------------------------
		//Define the displacement domain
		if(displace.domain_type[i]=="Point") 
		{
			fpoi.x = displace.coef[i][0];
			fpoi.y = displace.coef[i][1];
		}
		else if(displace.domain_type[i]=="Line") 
		{
			flin.point[0].x	=	displace.coef[i][0];
			flin.point[0].y	=	displace.coef[i][1];
			flin.point[1].x	=	displace.coef[i][2];
			flin.point[1].y	=	displace.coef[i][3];
		}
		else if(displace.domain_type[i]=="Zone")
		{
			zone.xmin = displace.coef[i][0];
			zone.xmax = displace.coef[i][1];
			zone.ymin = displace.coef[i][2];
			zone.ymax = displace.coef[i][3];
		}
		else if(displace.domain_type[i]!="All_boundaries")
		{
			cout << "Error: Failed to find the type of displacement domain!" << endl; 
			hout << "Error: Failed to find the type of displacement domain!" << endl;
			return 0; 
		}

		//-----------------------------------------------------------
		//Define the sign of displacement constraint
		int sign = 4;   //the number start from 4
		for(int j=0; j<(int)displace.disp_type[i].size(); j++) 
		{
			if(displace.disp_type[i][j]=="Fixed_displacement_x") sign += 1;
			else if(displace.disp_type[i][j]=="Fixed_displacement_y") sign += 2;
			else if(displace.disp_type[i][j]!="Pure_shear")
			{
				cout << "Error: Failed to find the type of displacement constraint!" << endl; 
				hout << "Error: Failed to find the type of displacement constraint!" << endl;
				return 0; 
			}
		}

		//-----------------------------------------------------------
		//Define the value of displacement constraint by the ratio
		vector<double> value;
		for(int j=0; j<(int)displace.value[i].size(); j++)  value.push_back(displace.value[i][j]*rmp);
		
		//-----------------------------------------------------------
		int count = 0;
		for(int j=0; j<(int)nodes.size(); j++)
		{
			if((displace.domain_type[i]=="Point"&&fpoi.distance_to(nodes[j].x, nodes[j].y)<=Zero)||
				(displace.domain_type[i]=="Line"&&flin.contain(nodes[j].x, nodes[j].y))||
				(displace.domain_type[i]=="Zone"&&nodes[j].x-zone.xmin>=Zero&&zone.xmax-nodes[j].x>=Zero&&
				  nodes[j].y-zone.ymin>=Zero&&zone.ymax-nodes[j].y>=Zero))
			{
				count ++;	//Find a point
				ip.push_back(j);
				ip.push_back(sign);
				//Initialize the displacement constraints
				//--------------------------------------------------------------------------------------------------------
				//Type	u		v		w		空表示这一维是自由边界(在此处统一赋零值，但是不会被调用)  
				//   7		0		0		0		其它0处表示是零值或非零值 但必须有值
				//	 6		空	0		0		
				//   5		0		空	0    
				//	 4		空	空	0
				//   3		0		0		空
				//   2		空	0		空
				//   1		0		空	空
				//--------------------------------------------------------------------------------------------------------		
				switch(sign)
				{
					case 1:	vp.push_back(value[0]);	vp.push_back(0);				vp.push_back(0);	break;
					case 2:	vp.push_back(0);				vp.push_back(value[0]);	vp.push_back(0);	break;
					case 3:	vp.push_back(value[0]);	vp.push_back(value[1]);	vp.push_back(0);	break;
					case 4:	vp.push_back(0);				vp.push_back(0);				vp.push_back(0);	break;
					case 5:	vp.push_back(value[0]);	vp.push_back(0);				vp.push_back(0);	break;
					case 6:	vp.push_back(0);				vp.push_back(value[0]);	vp.push_back(0);	break;
					case 7:	vp.push_back(value[0]);	vp.push_back(value[1]);	vp.push_back(0);	break;
					default:	
					{
						hout << "Error: the signs of nodes conformed to displacement constraint!" << endl;						
						hout << "Error: the signs of nodes conformed to displacement constraint!" << endl;	
						return 0;
					}
				}
				bnod_num++;  //the number of nodes on the boundary
			}
			else if(displace.domain_type[i]=="All_boundaries"&&displace.disp_type[i][0]=="Pure_shear"&&nodes[j].type!=0)		//Points on the boundaries
			{
				count ++;	//Find a point
				
				const double ux = nodes[j].y*value[0];		//相当于	[	0,		val	]		[x]
				const double uy = nodes[j].x*value[0];		//				[	val,	0,		]	*	[y]		

				ip.push_back(j);		ip.push_back(7);
				vp.push_back(ux);		vp.push_back(uy);		vp.push_back(0); 

				bnod_num++;  //the number of nodes on the boundary
			}
		}
		if(count==0) 
		{ 
			cout << "Error: Failed to find any nodes conformed to the displacement constraint!" << endl; 
			hout << "Error: Failed to find any nodes conformed to the displacement constraint!" << endl; 		
			return 0; 
		}
	}

	return 1;
}
//===========================================================================
