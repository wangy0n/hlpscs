//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	MatPro.cpp
//OBJECTIVE:	Material properties
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================

#include"MatPro.h"

MatPro::MatPro(const struct Stif_nonloc &stif_nonloc)
{
	type_val = stif_nonloc.type;
	E11 = stif_nonloc.E11;
	E22 = stif_nonloc.E22;
	E33 = stif_nonloc.E33;
	Nu12 = stif_nonloc.Nu12;
	Nu23 = stif_nonloc.Nu23;
	Nu13 = stif_nonloc.Nu13;
	G12 = stif_nonloc.G12;
	G23 = stif_nonloc.G23;
	G13 = stif_nonloc.G13;
}
//------------------------------------------------------------------------------
//设置各向同性材料参数的弹性模量、泊松比和剪切模量
void MatPro::set_ela_para(const double &E)
{
    E11 = E; E22 = E; E33 = E;
	Nu12 = 0.25; Nu23 = 0.25; Nu13 = 0.25;
	type_val = "Isotropic";
}
//---------------------------------------------------------------------
//设置横观各向同性材料的弹性模量、泊松比和剪切模量
void MatPro::set_ela_para(const double &E1, const double &E2, const double &Nu2)
{
	E11 = E1; E22 = E1; E33 = E2;
	Nu23 = Nu2; Nu13 = Nu2;
	Nu12 = (E1-4*Nu2*Nu2*E2)/(3*E1); 
	type_val = 1;
}
//---------------------------------------------------------------------
//设置正交各向异性材料的弹性模量、泊松比和剪切模量
void MatPro::set_ela_para(const double &iE1, const double &iE2, const double &iE3, const double &iNu12, const double &iNu23, const double &iNu13)
{
	E11 = iE1; E22 = iE2; E33 = iE3;
	Nu12 = iNu12; Nu23 = iNu23; Nu13 = iNu13;
	type_val = 2;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//生成Dij弹性矩阵，此弹性矩阵对应的应变向量为（e11,e22,e33,2*e12,2*e23,2*e31）转置
int MatPro::Generate_elas_matrix()
{
	for(int i=0; i<6; i++)
		for(int j=0; j<6; j++)
			elas_matrix[i][j] = 0.0;	

	elas_matrix[0][0] = E11*E11*(E22-Nu23*Nu23*E33)/(E11*E22-Nu12*Nu12*E22*E22-Nu23*Nu23*E11*E33-Nu13*Nu13*E22*E33-2*Nu12*Nu23*Nu13*E22*E33);
	elas_matrix[1][1] = E22*E22*(E11-Nu13*Nu13*E33)/(E11*E22-Nu12*Nu12*E22*E22-Nu23*Nu23*E11*E33-Nu13*Nu13*E22*E33-2*Nu12*Nu23*Nu13*E22*E33);
	elas_matrix[2][2] = E22*E33*(E11-Nu12*Nu12*E22)/(E11*E22-Nu12*Nu12*E22*E22-Nu23*Nu23*E11*E33-Nu13*Nu13*E22*E33-2*Nu12*Nu23*Nu13*E22*E33);
	elas_matrix[0][1] = E11*E22*(Nu12*E22+Nu23*Nu13*E33)/(E11*E22-Nu12*Nu12*E22*E22-Nu23*Nu23*E11*E33-Nu13*Nu13*E22*E33-2*Nu12*Nu23*Nu13*E22*E33);
	elas_matrix[0][2] = E11*E22*E33*(Nu13+Nu12*Nu23)/(E11*E22-Nu12*Nu12*E22*E22-Nu23*Nu23*E11*E33-Nu13*Nu13*E22*E33-2*Nu12*Nu23*Nu13*E22*E33);
	elas_matrix[1][2] = E22*E33*(Nu23*E11+Nu12*Nu13*E22)/(E11*E22-Nu12*Nu12*E22*E22-Nu23*Nu23*E11*E33-Nu13*Nu13*E22*E33-2*Nu12*Nu23*Nu13*E22*E33);
	elas_matrix[1][0] = elas_matrix[0][1];
	elas_matrix[2][0] = elas_matrix[0][2];
	elas_matrix[2][1] = elas_matrix[1][2];
	elas_matrix[3][3] = elas_matrix[0][1];
	elas_matrix[4][4] = elas_matrix[1][2];
	elas_matrix[5][5] = elas_matrix[0][2];

	return 1;
}
//------------------------------------------------------------------------------
//生成Dij弹性矩阵，此弹性矩阵对应的应变向量为（e11,e22,e33,2*e12,2*e23,2*e31）转置
int MatPro::Generate_elas_matrix_2D()
{
	for(int i=0; i<6; i++)
		for(int j=0; j<6; j++)
			elas_matrix[i][j] = 0.0;	

	elas_matrix[0][0] = E11*E11*(E22-Nu23*Nu23*E33)/(E11*E22-Nu12*Nu12*E22*E22-Nu23*Nu23*E11*E33-Nu13*Nu13*E22*E33-2*Nu12*Nu23*Nu13*E22*E33);
	elas_matrix[1][1] = E22*E22*(E11-Nu13*Nu13*E33)/(E11*E22-Nu12*Nu12*E22*E22-Nu23*Nu23*E11*E33-Nu13*Nu13*E22*E33-2*Nu12*Nu23*Nu13*E22*E33);
	elas_matrix[2][2] = E22*E33*(E11-Nu12*Nu12*E22)/(E11*E22-Nu12*Nu12*E22*E22-Nu23*Nu23*E11*E33-Nu13*Nu13*E22*E33-2*Nu12*Nu23*Nu13*E22*E33);
	elas_matrix[0][1] = E11*E22*(Nu12*E22+Nu23*Nu13*E33)/(E11*E22-Nu12*Nu12*E22*E22-Nu23*Nu23*E11*E33-Nu13*Nu13*E22*E33-2*Nu12*Nu23*Nu13*E22*E33);
	elas_matrix[0][2] = E11*E22*E33*(Nu13+Nu12*Nu23)/(E11*E22-Nu12*Nu12*E22*E22-Nu23*Nu23*E11*E33-Nu13*Nu13*E22*E33-2*Nu12*Nu23*Nu13*E22*E33);
	elas_matrix[1][2] = E22*E33*(Nu23*E11+Nu12*Nu13*E22)/(E11*E22-Nu12*Nu12*E22*E22-Nu23*Nu23*E11*E33-Nu13*Nu13*E22*E33-2*Nu12*Nu23*Nu13*E22*E33);
	elas_matrix[1][0] = elas_matrix[0][1];
	elas_matrix[2][0] = elas_matrix[0][2];
	elas_matrix[2][1] = elas_matrix[1][2];
	elas_matrix[3][3] = elas_matrix[0][1];
	elas_matrix[4][4] = G23;
	elas_matrix[5][5] = G13;

	return 1;
}
//------------------------------------------------------------------------------
//根据Dij弹性矩阵，反求材料参数
int  MatPro::Get_ele_para_by_ela_matrix()
{
	//首先生成柔度矩阵S;
    MathMatrix D(6,6);
	for(int i=0; i<6; i++)
		for(int j=0; j<6; j++)
			D.element[i][j]=elas_matrix[i][j];
	
	//求其逆
    MathMatrix S(6,6);
	S=D.Inverse();

	//求材料参数
	E11=1.0/S.element[0][0];
	E22=1.0/S.element[1][1];
	E33=1.0/S.element[2][2];
	Nu12=-S.element[0][1]*E11;
	Nu13=-S.element[0][2]*E11;
	Nu23=-S.element[1][2]*E22;
	G12=1.0/S.element[3][3];
	G23=1.0/S.element[4][4];
	G13=1.0/S.element[5][5];

	return 1;
}
//------------------------------------------------------------------------------
//根据理论公式反求系数，并验证一致性（由于要除以Horizon半径，所以Horizon半径应该大于1，避免误差过大）
void MatPro::Compare_coef_by_analysis_formula(const int mat_type, const double &decayR, const double decay_acoe[])const
{
	double delt5 = pow(decayR,5.0);
	double acoe[6] = { 0 };
	
	if(mat_type==1)
	{
		acoe[0] = 5*(8*elas_matrix[0][0]+12*elas_matrix[0][2]+3*elas_matrix[2][2])/(6*delt5*PI);
		acoe[1] = -25*(4*elas_matrix[0][0]-3*elas_matrix[0][2]-3*elas_matrix[2][2])/(6*delt5*PI);
		acoe[3] = 45*(elas_matrix[0][0]-6*elas_matrix[0][2]+elas_matrix[2][2])/(2*delt5*PI);
	}
	else if(mat_type==2)
	{
		acoe[0] = 5*(elas_matrix[0][0]+2*elas_matrix[0][1]+2*elas_matrix[0][2]+elas_matrix[1][1]+2*elas_matrix[1][2]+elas_matrix[2][2])/(2*delt5*PI);
		acoe[1] = -25*(elas_matrix[0][0]+2*elas_matrix[0][1]-elas_matrix[0][2]+elas_matrix[1][1]-elas_matrix[1][2]-2*elas_matrix[2][2])/(4*delt5*PI);
		acoe[2] = 25*(elas_matrix[0][0]+elas_matrix[0][2]-elas_matrix[1][1]-elas_matrix[1][2])/(8*delt5*PI);
		acoe[3] = 45*(3*elas_matrix[0][0]+6*elas_matrix[0][1]-24*elas_matrix[0][2]+3*elas_matrix[1][1]-24*elas_matrix[1][2]+8*elas_matrix[2][2])/(16*delt5*PI);
		acoe[4] = -15*(elas_matrix[0][0]-6*elas_matrix[0][2]-elas_matrix[1][1]+6*elas_matrix[1][2])/(16*delt5*PI);
		acoe[5] = 15*(elas_matrix[0][0]-6*elas_matrix[0][1]+elas_matrix[1][1])/(128*delt5*PI);
	}

	hout << endl << "The values of analysis formula: " << endl;
	for(int i=0; i<6; i++) hout << acoe[i] << "  ";
	hout << endl;

	hout << "The values of numerical formula: " << endl;
	for(int i=0; i<6; i++) hout << decay_acoe[i] << "  ";
	hout << endl << endl;
}
//------------------------------------------------------------------------------
//根据理论公式反求刚度矩阵，并验证一致性
void MatPro::Compare_matrix_by_analysis_formula(const double &decayR, const double acoe[])const
{
	double delt5 = pow(decayR,5.0);
	double elasmat[6][6] = { {0}, {0}, {0}, {0}, {0}, {0} };

	elasmat[0][0] = 2*(21*acoe[0]-6*acoe[1]+36*acoe[2]+acoe[3]-20*acoe[4]+280*acoe[5])*delt5*PI/525;
	elasmat[1][1] = 2*(21*acoe[0]-6*acoe[1]-36*acoe[2]+acoe[3]+20*acoe[4]+280*acoe[5])*delt5*PI/525;
	elasmat[2][2] = 2*(63*acoe[0]+36*acoe[1]+8*acoe[3])*delt5*PI/1575;
	elasmat[0][1] = 2*(21*acoe[0]-6*acoe[1]+acoe[3]-840*acoe[5])*delt5*PI/1575;
	elasmat[1][2] = 2*(21*acoe[0]+3*acoe[1]-18*acoe[2]-4*acoe[3]-60*acoe[4])*delt5*PI/1575;
	elasmat[0][2] = 2*(21*acoe[0]+3*acoe[1]+18*acoe[2]-4*acoe[3]+60*acoe[4])*delt5*PI/1575;

	elasmat[1][0] = elasmat[0][1];
	elasmat[2][1] = elasmat[1][2];
	elasmat[2][0] = elasmat[0][2];
	elasmat[3][3] = elasmat[0][1];
	elasmat[4][4] = elasmat[1][2];
	elasmat[5][5] = elasmat[0][2];

	hout<<"The matrix by analysis forlmula:"<<endl;
	for(int i=0; i<6; i++){
		for(int j=0; j<6; j++)
			hout << elasmat[i][j] << "  ";
		hout<<endl;}
	hout<<endl;

	hout<<"The matrix by numerical forlmula:"<<endl;
	for(int i=0; i<6; i++){
		for(int j=0; j<6; j++)
			hout << elas_matrix[i][j] << "  ";
		hout<<endl;}
	hout<<endl;
}
//------------------------------------------------------------------------------
void MatPro::print()const
{
	//弹性矩阵输出测试
	cout<<"Stiffness_matrix_Dij:"<<endl;
	for(int i=0; i<6; i++){
		for(int j=0; j<6; j++)
			cout<<elas_matrix[i][j]<<"  ";
		cout<<endl;}
	cout<<endl;
}
//===========================================================================
