//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Global_Stiff_Matrix.cpp
//OBJECTIVE:	Compute and assemble the global stiffness matrix
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================

#include "Global_Stiff_Matrix.h"

//---------------------------------------------------------------------------
//Estimate element matrices
int Global_Stiff_Matrix::Gen_element_matrices(const int &gaussnum, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func,
																					const vector<MatPro> &mats, vector<Node> &nodes, const vector<Element> &elements, const vector<bool> &break_table, 
																					vector<vector<double> > &Gp_val, vector<MathMatrix> &ele_self_matrix, vector<vector<MathMatrix> > &ele_relative_matrix)const
{
	//---------------------------------------------------------------------------
	//Generating gaussian points
	Gauss gau;		//gaussian potins
	if(gau.Generate_gauss(gaussnum)==0) return 0;

	//------------------------------------------------------------------------
	//Initialization of variables
	cout << "    -_- Initializing variables..." << endl;
	hout << "    -_- Initializing variables..." << endl;

	const int GS = (int)gau.gauss.size();
	const int ES = (int)elements.size();
	double (*gauss_ns)[8] = new double [GS][8];				//记录单元高斯点的形函数
	double (*gauss_nw)[8] = new double [GS][8];				//记录单元高斯点的形函数(带权重系数)
	double (*gauss_dif)[3][8] = new double [GS][3][8];		//记录单元高斯点形函数的导数
	double Jacobi = 0;																//标准立方体单元的雅可比值
	double (*gauss_po)[3] = new double [GS][3];				//标准立方体单元的高斯点坐标
	double (*ele_cent)[4] = new double [ES][4];					//记录单元中心点位置向量(x,y,z)分别放于[0],[1]和[2]中（以标准立方体单元中心点为原点）, [3]中放单元顶点权重的平均值
	
	//Calculate data of every gaussian point of element
	Generate_element_gauss_data(elements, nodes, com_mod, weight_func, gau.gauss, gau.weight, gauss_ns, gauss_nw, gauss_dif, Jacobi, gauss_po, ele_cent);

	WeightFunc Wei_Fun;
	//Export the values of weighting function at every nodes for testing
	Wei_Fun.Export_weight_function_contour("Weight_function_contour.dat", nodes, elements, gauss_ns, gauss_dif, Jacobi, gauss_po, gau.weight, weight_func, ele_cent);
	//Calculate the values of weighting function of gaussian points in all elements
	Wei_Fun.Value_weightfunc_gausspois_elements(weight_func, nodes, elements, GS, gauss_ns, gauss_po, ele_cent, Gp_val);

	cout << "    ^_^ Variables initialized successfully!" << endl << endl;
	hout << "    ^_^ Variables initialized successfully!" << endl << endl;

	//------------------------------------------------------------------------
	//Loops all elements
	cout << "    -_- Calculating the element matrices..." << endl;
	hout << "    -_- Calculating the element matrices..." << endl;
	//执行openmp
	#pragma omp parallel
	{
		double element_stiff_matrix1[24][24];
		double element_stiff_matrix2[24][24];
		double (*Long_range_stiffness)[6][6] = new double [GS][6][6];
		double *alpha_key = new double [GS];

		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<ES; i++) 
		{
			//Output calculating progress
			#pragma omp critical
			if(i%(ES/10)==0)	cout << i*100/ES << "% done" << endl;
			
			//Extracting coordinates of nodes in the principal element 
			double elenodes[8][3];	
			for(int j=0; j<8; j++)
			{
				elenodes[j][0] = nodes[elements[i].nodes_id[j]].x;
				elenodes[j][1] = nodes[elements[i].nodes_id[j]].y;
				elenodes[j][2] = nodes[elements[i].nodes_id[j]].z;
			}
			//---------------------------------------------------------------------------
			//Calculating hybrid models
			if(com_mod=="Hybrid")
			{
				//Time markers
				//clock_t ctn1,ctn2;
				//ctn1 = clock();
				
				//------------------------------------------------------------------------
				//Initialazing long-range stiffness tensor
				for(int j=0; j<GS; j++)
				{
					for(int k=0; k<6; k++)
						for(int m=0; m<6; m++)
							Long_range_stiffness[j][k][m] = 0.0;
					alpha_key[j] = -1.0;
				}

				//Defining parameters of the center of the principal element 
				double elec_left[4] = { ele_cent[i][0], ele_cent[i][1], ele_cent[i][2], ele_cent[i][3] };

				//------------------------------------------------------------------------
				//Accumulate the number of relative elements (less than i) //attention: Openmp Parallel
				long int ref_tot = 0;
				for(int j=0; j<i; j++) ref_tot += elements[j].relative_eles.size();

				//------------------------------------------------------------------------
				//Loops element neighbours
				for(int j=0; j<(int)elements[i].relative_eles.size(); j++)
				{
					ref_tot++; //Accumlation
					//--------------------------------------------------
					//Defining parameters of the center of the associated element 
					const int ere = elements[i].relative_eles[j]; 
					double elec_right[4] = { ele_cent[ere][0], ele_cent[ere][1], ele_cent[ere][2], ele_cent[ere][3] };
					
					//--------------------------------------------------
					//A preliminary examination for skipping both elements completely in the local continuum area
					if(elec_left[3]<=Zero&&elec_right[3]<=Zero) 
					{
						for(int count=0; count<GS; count++)
						{
							if(alpha_key[count]==-1.0)  alpha_key[count] = 0.0;
							else if(alpha_key[count]!=0.0) alpha_key[count] = 0.5;
						}
						continue;
					}

					//--------------------------------------------------
					//Extracting coordinates of nodes in the associated element
					double elenodes_relative[8][3];
					for(int k=0; k<8; k++)
					{
						elenodes_relative[k][0] = nodes[elements[ere].nodes_id[k]].x;
						elenodes_relative[k][1] = nodes[elements[ere].nodes_id[k]].y;
						elenodes_relative[k][2] = nodes[elements[ere].nodes_id[k]].z;
					}

					//--------------------------------------------------
					//Calculate the position of break_table
					const long int pos_bt = (ref_tot-1)*GS*GS;
					//------------------------------------------------------------------------
					//Generate the element matrices based on non-local continuum model (long-range forces)
					Generate_Longforce_Elestiff_Brokeless(element_stiff_matrix1, element_stiff_matrix2, alpha_key, Long_range_stiffness, weight_func, elements[i].flag, elements[ere].flag, elenodes, 
																							elenodes_relative, gauss_ns, gauss_nw, gauss_dif, gau.weight, peri_para, Jacobi, gauss_po, elec_left, elec_right, break_table, pos_bt, Gp_val, i, ere);

					#pragma omp critical
					{
						//---------------------------------------------------------------------------
						//Accumulate the new element matrices
						for(int k=0; k<24; k++)
							for(int m=0; m<24; m++)
							{
								ele_self_matrix[i].element[k][m] += element_stiff_matrix1[k][m];
								ele_relative_matrix[i][j].element[k][m]=element_stiff_matrix2[k][m];
							}
					}
				}
				//Attention: a half of parameters
				for(int j=0; j<GS; j++)
					for(int k=0; k<6; k++)
						for(int m=0; m<6; m++)
							Long_range_stiffness[j][k][m] = Long_range_stiffness[j][k][m]/2.0;

				//------------------------------------------------------------------------
				//ctn2 = clock();
				//hout << "Total num of elements: " << (int)elements.size() << "; Element " << i << " took time: " << (double)(ctn2-ctn1)/CLOCKS_PER_SEC << "sec; " << endl;
			}
			else
			{
				//------------------------------------------------------------------------
				//Initialazing long-range stiffness tensor
				for(int j=0; j<GS; j++)
				{
					for(int k=0; k<6; k++)
						for(int m=0; m<6; m++)
							Long_range_stiffness[j][k][m] = 0.0;
					alpha_key[j] = 0.0;
				}
			}
			 
			//Generate the element matrix based on local continuum model (contact forces)
			Generate_Contacforce_Elestiff(element_stiff_matrix1, alpha_key, Long_range_stiffness, elenodes, mats, elements[i].flag, Jacobi, gauss_dif, gau.weight);
			#pragma omp critical
			{
				for(int k=0; k<24; k++)
					for(int m=0; m<24; m++)
						ele_self_matrix[i].element[k][m] += element_stiff_matrix1[k][m];
			}
		}
		//---------------------------------------------------------------------------
		//Delete pointers
		delete[] Long_range_stiffness;
		delete[] alpha_key;
	}

	//---------------------------------------------------------------------------
	//Delete pointers
	delete[] gauss_ns;
	delete[] gauss_nw;
	delete[] gauss_dif;
	delete[] gauss_po;
	delete[] ele_cent;

	cout << "^_^ The element matrices calculated successfully!" << endl;
	hout << "^_^ The element matrices calculated successfully!" << endl;

	return 1;
}
//---------------------------------------------------------------------------
//Estimate element matrices with local damage
int Global_Stiff_Matrix::Gen_element_matrices_damage(const int &gaussnum, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func, const struct Damage &damages,
																									const vector<MatPro> &mats, vector<Node> &nodes, const vector<Element> &elements, const vector<bool> &full_dam_eles, const vector<vector<double> > &damage_table,
																									const vector<bool> &break_table, vector<vector<double> > &Gp_val, vector<MathMatrix> &ele_self_matrix, vector<vector<MathMatrix> > &ele_relative_matrix)const
{
	//---------------------------------------------------------------------------
	//Generating gaussian points
	Gauss gau;		//gaussian potins
	if(gau.Generate_gauss(gaussnum)==0) return 0;

	//------------------------------------------------------------------------
	//Initialization of variables
	cout << "    -_- Initializing variables..." << endl;
	hout << "    -_- Initializing variables..." << endl;

	const int GS = (int)gau.gauss.size();
	const int ES = (int)elements.size();
	double (*gauss_ns)[8] = new double [GS][8];				//记录单元高斯点的形函数
	double (*gauss_nw)[8] = new double [GS][8];				//记录单元高斯点的形函数(带权重系数)
	double (*gauss_dif)[3][8] = new double [GS][3][8];		//记录单元高斯点形函数的导数
	double Jacobi = 0;																//标准立方体单元的雅可比值
	double (*gauss_po)[3] = new double [GS][3];				//标准立方体单元的高斯点坐标
	double (*ele_cent)[4] = new double [ES][4];					//记录单元中心点位置向量(x,y,z)分别放于[0],[1]和[2]中（以标准立方体单元中心点为原点）, [3]中放单元顶点权重的平均值
	
	//Calculate data of every gaussian point of element
	Generate_element_gauss_data(elements, nodes, com_mod, weight_func, gau.gauss, gau.weight, gauss_ns, gauss_nw, gauss_dif, Jacobi, gauss_po, ele_cent);

	WeightFunc Wei_Fun;
	//Export the values of weighting function at every nodes for testing
//	Wei_Fun.Export_weight_function_contour("Weight_function_contour.dat", nodes, elements, gauss_ns, gauss_dif, Jacobi, gauss_po, gau.weight, weight_func, ele_cent);
	
	//Calculate the values of weighting function of gaussian points in all elements
	Wei_Fun.Value_weightfunc_gausspois_elements(weight_func, nodes, elements, GS, gauss_ns, gauss_po, ele_cent, Gp_val);

	cout << "    ^_^ Variables initialized successfully!" << endl << endl;
	hout << "    ^_^ Variables initialized successfully!" << endl << endl;

	//------------------------------------------------------------------------
	//Loops all elements
	cout << "    -_- Calculating the element matrices..." << endl;
	hout << "    -_- Calculating the element matrices..." << endl;
	//执行openmp
	#pragma omp parallel
	{
		double element_stiff_matrix1[24][24];
		double element_stiff_matrix2[24][24];
		double (*Long_range_stiffness)[6][6] = new double [GS][6][6];
		double *alpha_key = new double [GS];

		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<ES; i++) 
		{
			//Output calculating progress
			#pragma omp critical
			if(i%(ES/10)==0)	cout << i*100/ES << "% done" << endl;
			
			//Extracting coordinates of nodes in the principal element 
			double elenodes[8][3];	
			for(int j=0; j<8; j++)
			{
				elenodes[j][0] = nodes[elements[i].nodes_id[j]].x;
				elenodes[j][1] = nodes[elements[i].nodes_id[j]].y;
				elenodes[j][2] = nodes[elements[i].nodes_id[j]].z;
			}
			//---------------------------------------------------------------------------
			//Calculating hybrid models
			if(com_mod=="Hybrid")
			{
				//Time markers
				//clock_t ctn1,ctn2;
				//ctn1 = clock();
				
				//------------------------------------------------------------------------
				//Initialazing long-range stiffness tensor
				for(int j=0; j<GS; j++)
				{
					for(int k=0; k<6; k++)
						for(int m=0; m<6; m++)
							Long_range_stiffness[j][k][m] = 0.0;
					alpha_key[j] = -1.0;
				}

				//Defining parameters of the center of the principal element 
				double elec_left[4] = { ele_cent[i][0], ele_cent[i][1], ele_cent[i][2], ele_cent[i][3] };

				//------------------------------------------------------------------------
				//Accumulate the number of relative elements (less than i) //attention: Openmp Parallel
				long int ref_tot = 0;
				for(int j=0; j<i; j++) ref_tot += elements[j].relative_eles.size();

				//------------------------------------------------------------------------
				//Loops element neighbours
				for(int j=0; j<(int)elements[i].relative_eles.size(); j++)
				{
					ref_tot++; //Accumlation
					//--------------------------------------------------
					//Defining parameters of the center of the associated element 
					const int ere = elements[i].relative_eles[j]; 
					double elec_right[4] = { ele_cent[ere][0], ele_cent[ere][1], ele_cent[ere][2], ele_cent[ere][3] };
					
					//--------------------------------------------------
					//A preliminary examination for skipping both elements completely in the local continuum area
					if(elec_left[3]<=Zero&&elec_right[3]<=Zero) 
					{
						for(int count=0; count<GS; count++)
						{
							if(alpha_key[count]==-1.0)  alpha_key[count] = 0.0;
							else if(alpha_key[count]!=0.0) alpha_key[count] = 0.5;
						}
						continue;
					}

					//--------------------------------------------------
					//Extracting coordinates of nodes in the associated element
					double elenodes_relative[8][3];
					for(int k=0; k<8; k++)
					{
						elenodes_relative[k][0] = nodes[elements[ere].nodes_id[k]].x;
						elenodes_relative[k][1] = nodes[elements[ere].nodes_id[k]].y;
						elenodes_relative[k][2] = nodes[elements[ere].nodes_id[k]].z;
					}

					//--------------------------------------------------
					//Calculate the position of break_table
					const long int pos_bt = (ref_tot-1)*GS*GS;
					//------------------------------------------------------------------------
					//Generate the element matrices based on non-local continuum model (long-range forces)
					Generate_Longforce_Elestiff_Brokeless(element_stiff_matrix1, element_stiff_matrix2, alpha_key, Long_range_stiffness, weight_func, elements[i].flag, elements[ere].flag, elenodes, 
																							elenodes_relative, gauss_ns, gauss_nw, gauss_dif, gau.weight, peri_para, Jacobi, gauss_po, elec_left, elec_right, break_table, pos_bt, Gp_val, i, ere);

					#pragma omp critical
					{
						//---------------------------------------------------------------------------
						//Accumulate the new element matrices
						for(int k=0; k<24; k++)
							for(int m=0; m<24; m++)
							{
								ele_self_matrix[i].element[k][m] += element_stiff_matrix1[k][m];
								ele_relative_matrix[i][j].element[k][m]=element_stiff_matrix2[k][m];
							}
					}
				}
				//Attention: a half of parameters
				for(int j=0; j<GS; j++)
					for(int k=0; k<6; k++)
						for(int m=0; m<6; m++)
							Long_range_stiffness[j][k][m] = Long_range_stiffness[j][k][m]/2.0;

				//------------------------------------------------------------------------
				//ctn2 = clock();
				//hout << "Total num of elements: " << (int)elements.size() << "; Element " << i << " took time: " << (double)(ctn2-ctn1)/CLOCKS_PER_SEC << "sec; " << endl;
			}
			else
			{
				//------------------------------------------------------------------------
				//Initialazing long-range stiffness tensor
				for(int j=0; j<GS; j++)
				{
					for(int k=0; k<6; k++)
						for(int m=0; m<6; m++)
							Long_range_stiffness[j][k][m] = 0.0;
					alpha_key[j] = 0.0;
				}
			}
			
			//Generate the element matrix based on local continuum model (coupled part + damaged part)
			vector<double> ele_dam_table(damage_table[i]);
			Generate_Contacforce_Elestiff_Coupled_Damaged(element_stiff_matrix1, alpha_key, Long_range_stiffness, damages, full_dam_eles[i], ele_dam_table, elenodes, mats, elements[i].flag, elements[i].mat, Jacobi, gauss_dif, gau.weight);
			#pragma omp critical
			{
				for(int k=0; k<24; k++)
					for(int m=0; m<24; m++)
						ele_self_matrix[i].element[k][m] += element_stiff_matrix1[k][m];
			}
		}
		//---------------------------------------------------------------------------
		//Delete pointers
		delete[] Long_range_stiffness;
		delete[] alpha_key;
	}

	//---------------------------------------------------------------------------
	//Delete pointers
	delete[] gauss_ns;
	delete[] gauss_nw;
	delete[] gauss_dif;
	delete[] gauss_po;
	delete[] ele_cent;

	cout << "^_^ The element matrices calculated successfully!" << endl;
	hout << "^_^ The element matrices calculated successfully!" << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------
//Calculate data of every gaussian point of element
void Global_Stiff_Matrix::Generate_element_gauss_data(const vector<Element> &elements, const vector<Node> &nodes, const string &com_mod, const struct Weight_func &weight_func, const vector<Node> &gauss, 
																								    const vector<double> &weight, double (*gauss_ns)[8], double (*gauss_nw)[8], double (*gauss_dif)[3][8], double &Jacobi,  double (*gauss_po)[3], double (*ele_cent)[4])const
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Element loop
	int inum = -1;
	for(int i=0; i<(int)elements.size(); i++)
	{
		if(elements[i].flag==0) //Standard cubic element
		{
			inum = i;
			break;
		}
	}

	for(int count=0; count<(int)gauss.size(); count++)
	{	
		double Nshape[8] = {0};
		if(com_mod=="Hybrid") 
		{
			//--------------------------------------------
			//计算高斯积分点整体坐标
			Nshape[0]=0.125*(1.0-gauss[count].x)*(1.0-gauss[count].y)*(1.0-gauss[count].z);
			Nshape[1]=0.125*(1.0+gauss[count].x)*(1.0-gauss[count].y)*(1.0-gauss[count].z);
			Nshape[2]=0.125*(1.0+gauss[count].x)*(1.0+gauss[count].y)*(1.0-gauss[count].z);
			Nshape[3]=0.125*(1.0-gauss[count].x)*(1.0+gauss[count].y)*(1.0-gauss[count].z);
			Nshape[4]=0.125*(1.0-gauss[count].x)*(1.0-gauss[count].y)*(1.0+gauss[count].z);
			Nshape[5]=0.125*(1.0+gauss[count].x)*(1.0-gauss[count].y)*(1.0+gauss[count].z);
			Nshape[6]=0.125*(1.0+gauss[count].x)*(1.0+gauss[count].y)*(1.0+gauss[count].z);
			Nshape[7]=0.125*(1.0-gauss[count].x)*(1.0+gauss[count].y)*(1.0+gauss[count].z);

			for(int j=0; j<8; j++) 
			{
				gauss_ns[count][j] = Nshape[j]; 
				gauss_nw[count][j] = Nshape[j]*weight[count];
			} 				
			//--------------------------------------------
			//计算高斯积分点坐标
			Point_3D gaupoi(0, 0, 0);		
			if(inum>=0)
			{
				for(int j=0; j<8; j++) 
				{
					gaupoi.x += Nshape[j]*nodes[elements[inum].nodes_id[j]].x;
					gaupoi.y += Nshape[j]*nodes[elements[inum].nodes_id[j]].y;
					gaupoi.z += Nshape[j]*nodes[elements[inum].nodes_id[j]].z;
				}
			}

			//记录高斯积分点坐标
			gauss_po[count][0] = gaupoi.x;
			gauss_po[count][1] = gaupoi.y;
			gauss_po[count][2] = gaupoi.z;
		}
		else
		{
			//初始化
			gauss_po[count][0] = 0.0;
			gauss_po[count][1] = 0.0;
			gauss_po[count][2] = 0.0;
		}

		//--------------------------------------------
		//计算Ｊ矩阵
		//--------------------------------------------
		//形函数N对gauss[count].x, gauss[count].y, gauss[count].z的偏导矩阵
		double diff[3][8];
		diff[0][0]=-0.125*(1.0-gauss[count].y)*(1.0-gauss[count].z);
		diff[0][1]=-diff[0][0];                         
		diff[0][2]=0.125*(1.0+gauss[count].y)*(1.0-gauss[count].z);
		diff[0][3]=-diff[0][2];
		diff[0][4]=-0.125*(1.0-gauss[count].y)*(1.0+gauss[count].z);
		diff[0][5]=-diff[0][4];
		diff[0][6]=0.125*(1.0+gauss[count].y)*(1.0+gauss[count].z);
		diff[0][7]=-diff[0][6];

		diff[1][0]=-0.125*(1.0-gauss[count].x)*(1.0-gauss[count].z);
		diff[1][1]=-0.125*(1.0+gauss[count].x)*(1.0-gauss[count].z);
		diff[1][2]=-diff[1][1];
		diff[1][3]=-diff[1][0];
		diff[1][4]=-0.125*(1.0-gauss[count].x)*(1.0+gauss[count].z);
		diff[1][5]=-0.125*(1.0+gauss[count].x)*(1.0+gauss[count].z);
		diff[1][6]=-diff[1][5];
		diff[1][7]=-diff[1][4];

		diff[2][0]=-0.125*(1.0-gauss[count].x)*(1.0-gauss[count].y);
		diff[2][1]=-0.125*(1.0+gauss[count].x)*(1.0-gauss[count].y);
		diff[2][2]=-0.125*(1.0+gauss[count].x)*(1.0+gauss[count].y);
		diff[2][3]=-0.125*(1.0-gauss[count].x)*(1.0+gauss[count].y);
		diff[2][4]=-diff[2][0];
		diff[2][5]=-diff[2][1];
		diff[2][6]=-diff[2][2];
		diff[2][7]=-diff[2][3];

		//记录每个高斯点形函数偏导数的值
		for(int j=0; j<3; j++)
			for(int k=0; k<8; k++)
				gauss_dif[count][j][k] = diff[j][k];

		if(inum>=0)
		{
			//--------------------------------------------------
			//单元节点坐标矩阵
			double elenode[8][3];
			for(int j=0; j<8; j++)
			{
				elenode[j][0]=nodes[elements[inum].nodes_id[j]].x;
				elenode[j][1]=nodes[elements[inum].nodes_id[j]].y;
				elenode[j][2]=nodes[elements[inum].nodes_id[j]].z;
			}
			//--------------------------------------------------
			//J矩阵
			double Jmatrix[3][3];
			//以上两个矩阵的积
			for(int j=0; j<3; j++)
				for(int k=0; k<3; k++)
				{
					Jmatrix[j][k]=0;
					for(int m=0; m<8; m++)
					Jmatrix[j][k] += diff[j][m]*elenode[m][k];
				}
			//--------------------------------------------------
			//求出J矩阵的行列式
			Jacobi = Jmatrix[0][0]*(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])
							-Jmatrix[0][1]*(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])
							+Jmatrix[0][2]*(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0]);
		}
		else 
		{
			Jacobi = 0.0;
		}
	}

	//The coordinates of center
	Point_3D scenter(0,0,0);	//Standard cubic element
	if(inum>=0)
	{
		for(int j=0; j<8; j++) 
		{
			scenter.x += nodes[elements[inum].nodes_id[j]].x;
			scenter.y += nodes[elements[inum].nodes_id[j]].y;
			scenter.z += nodes[elements[inum].nodes_id[j]].z;
		}
		scenter = scenter/8;
	}
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Element loop
	for(int i=0; i<(int)elements.size(); i++)
	{
		for(int j=0; j<4; j++) ele_cent[i][j] = 0.0;
		//计算该单元的中心点坐标
		for(int j=0; j<8; j++) 
		{
			ele_cent[i][0] += nodes[elements[i].nodes_id[j]].x;
			ele_cent[i][1] += nodes[elements[i].nodes_id[j]].y;
			ele_cent[i][2] += nodes[elements[i].nodes_id[j]].z;
			const Point_3D poi_nod(nodes[elements[i].nodes_id[j]].x, nodes[elements[i].nodes_id[j]].y, nodes[elements[i].nodes_id[j]].z);
			WeightFunc Wei_Fun;  //the class of weighting function
			ele_cent[i][3] += Wei_Fun.Value_weight_function(poi_nod, weight_func);
		}
		for(int j=0; j<4; j++) ele_cent[i][j] = ele_cent[i][j]/8;
		ele_cent[i][0] = ele_cent[i][0] - scenter.x;
		ele_cent[i][1] = ele_cent[i][1] - scenter.y;
		ele_cent[i][2] = ele_cent[i][2] - scenter.z;
	}
}
//-----------------------------------------------------------------------------------------------
//Generate the element matrix based on non-local continuum model (long-range forces) without broken bonds
void Global_Stiff_Matrix::Generate_Longforce_Elestiff_Brokeless(double (*element_stiff_matrix1)[24], double (*element_stiff_matrix2)[24], double *alpha_key, double (*Long_range_stiffness)[6][6], const struct Weight_func &weight_func,
																													const int &flag_left, const int &flag_right, const double (*elenodes_left)[3], const double (*elenodes_right)[3], const double (*gauss_ns)[8], const double (*gauss_nw)[8],
																													const double (*gauss_dif)[3][8], const vector<double> &weight, const struct Peri_para &peri_para, const double &Jacobi, const double (*gauss_po)[3], const double elec_left[],
																													const double elec_right[], const vector<bool> &break_table, const long int &pos_bt, const vector<vector<double> > &Gp_val, const int &pri_ele, const int &ass_ele)const
{
	long int cogau = pos_bt;
	//--------------------------------------------------
	//Clearing
	for(int i=0; i<24; i++) 
		for (int j=0; j<24; j++)
		{
			element_stiff_matrix1[i][j] = 0;
			element_stiff_matrix2[i][j] = 0;
		}

	//--------------------------------------------------	
	//Loop the gaussian points in the current element
	const int gau_num = (int)weight.size();
	for(int count1=0; count1<gau_num; count1++)
	{
		//--------------------------------------------------	
		double Gmat[6] = {0};
		double GNmat[6][8] = { {0}, {0}, {0}, {0}, {0}, {0} };
		//--------------------------------------------------	
		//左端高斯点坐标
		Point_3D gaupoi_left(0, 0, 0);
		if(flag_left==0)
		{
			gaupoi_left.x = gauss_po[count1][0] + elec_left[0];
			gaupoi_left.y = gauss_po[count1][1] + elec_left[1];
			gaupoi_left.z = gauss_po[count1][2] + elec_left[2];
		}
		else
		{
			for(int i=0; i<8; i++) 
			{
				gaupoi_left.x += gauss_ns[count1][i]*elenodes_left[i][0];
				gaupoi_left.y += gauss_ns[count1][i]*elenodes_left[i][1];
				gaupoi_left.z += gauss_ns[count1][i]*elenodes_left[i][2];
			}
		}

		//--------------------------------------------
		//计算左端高斯点的权重
		double weight_left;
		if(elec_left[3]<=Zero) weight_left = 0.0;
		else if(fabs(elec_left[3]-1.0)<=Zero) weight_left = 1.0;
		else weight_left = Gp_val[pri_ele][count1];

		//计算Ｊ矩阵；
		//--------------------------------------------------
		double Jleft_val;
		if(flag_left==0) Jleft_val = Jacobi;
		else
		{
			//J矩阵
			double Jmatrix_left[3][3];
			//以上两个矩阵的积
			for(int i=0; i<3; i++)
				for(int j=0; j<3; j++)
				{
					Jmatrix_left[i][j]=0;
					for(int k=0; k<8; k++)
						Jmatrix_left[i][j] += gauss_dif[count1][i][k]*elenodes_left[k][j];
				}

			//求出J矩阵的行列式
			Jleft_val = Jmatrix_left[0][0]*(Jmatrix_left[1][1]*Jmatrix_left[2][2]-Jmatrix_left[1][2]*Jmatrix_left[2][1])
								-Jmatrix_left[0][1]*(Jmatrix_left[1][0]*Jmatrix_left[2][2]-Jmatrix_left[1][2]*Jmatrix_left[2][0])
								+Jmatrix_left[0][2]*(Jmatrix_left[1][0]*Jmatrix_left[2][1]-Jmatrix_left[1][1]*Jmatrix_left[2][0]);	
		}
		
		//--------------------------------------------------
		//左端带权重形函数
		double gaus_nwct1[8];
		for(int i=0; i<8; i++) gaus_nwct1[i] = gauss_nw[count1][i]*Jleft_val;

		//------------------------------------------------------------------------------------------------------------------------
		//循环外单元高斯点计算积分
		for(int count2=0; count2<gau_num; count2++)
		{
			if(break_table[cogau++]) continue;  //First use and then ++, same to: { break_table[cogau] continue; cogau++; } 
			//--------------------------------------------
			//右端高斯点坐标
			Point_3D gaupoi_right(0, 0, 0);
			if(flag_right==0)
			{
				gaupoi_right.x = gauss_po[count2][0] + elec_right[0];
				gaupoi_right.y = gauss_po[count2][1] + elec_right[1];
				gaupoi_right.z = gauss_po[count2][2] + elec_right[2];
			}
			else
			{
				for(int i=0; i<8; i++) 
				{
					gaupoi_right.x +=  gauss_ns[count2][i]*elenodes_right[i][0];
					gaupoi_right.y +=  gauss_ns[count2][i]*elenodes_right[i][1];
					gaupoi_right.z +=  gauss_ns[count2][i]*elenodes_right[i][2];
				}
			}

			//--------------------------------------------
			//计算右端高斯点的权重
			double weight_right;
			if(elec_right[3]<=Zero) weight_right = 0.0;
			else if(fabs(elec_right[3]-1.0)<=Zero) weight_right = 1.0;
			else weight_right = Gp_val[ass_ele][count2];			

			//--------------------------------------------
			//计算权重函数值
			double alpha = 0.5*(weight_left+weight_right);

			if(alpha<1.0&&alpha>0.0) alpha_key[count1] = 0.5;
			else if(alpha_key[count1]==-1.0)  alpha_key[count1] = alpha;
			else if(alpha_key[count1]!=alpha) alpha_key[count1] = 0.5;

			if(alpha<=Zero) continue;  //等于零值，没有长程效果

			//--------------------------------------------
			//计算长程作用衰减函数值
			const double x = gaupoi_right.x-gaupoi_left.x;
			const double y = gaupoi_right.y-gaupoi_left.y;
			const double z = gaupoi_right.z-gaupoi_left.z;

			const double dis_squr = x*x + y*y + z*z;
			const double poi_dis = sqrt(dis_squr);		//两点之间的距离
			
			if(poi_dis>peri_para.horizon_R||poi_dis<Zero) continue;		//圆形积分域

			//计算Ｊ矩阵；
			double Jright_val;
			if(flag_right==0) Jright_val = Jacobi;
			else
			{
				//--------------------------------------------------
				//J矩阵
				double Jmatrix_right[3][3];
				//以上两个矩阵的积
				for(int i=0; i<3; i++)
					for(int j=0; j<3; j++)
					{
						Jmatrix_right[i][j]=0;
						for(int k=0; k<8; k++)
							Jmatrix_right[i][j]=Jmatrix_right[i][j] + gauss_dif[count2][i][k]*elenodes_right[k][j];
					}

				//求出J矩阵的行列式
				Jright_val = Jmatrix_right[0][0]*(Jmatrix_right[1][1]*Jmatrix_right[2][2]-Jmatrix_right[1][2]*Jmatrix_right[2][1])
												-Jmatrix_right[0][1]*(Jmatrix_right[1][0]*Jmatrix_right[2][2]-Jmatrix_right[1][2]*Jmatrix_right[2][0])
												+Jmatrix_right[0][2]*(Jmatrix_right[1][0]*Jmatrix_right[2][1]-Jmatrix_right[1][1]*Jmatrix_right[2][0]);
			}
	
			//--------------------------------------------
			//计算长程力矩阵
			const double cos2sita = z*z/dis_squr;		//cos(sita)^2
			double cos2pha;
			if(fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
			else cos2pha = x*x/(x*x + y*y);					//cos(pha)^2

			double sum = peri_para.acoe[0] + peri_para.acoe[1]*0.5*(3*cos2sita-1) + peri_para.acoe[2]*(2*cos2pha-1)*3*(1-cos2sita) + peri_para.acoe[3]*0.125*(35*cos2sita*cos2sita - 30*cos2sita + 3)
									+ peri_para.acoe[4]*(2*cos2pha-1)*7.5*(7*cos2sita-1)*(1-cos2sita) + peri_para.acoe[5]*(8*cos2pha*cos2pha-8*cos2pha+1)*105*(1-cos2sita)*(1-cos2sita);
			
			const double gv = exp(-poi_dis/peri_para.intrinsic_L)*sum*alpha*Jright_val*weight[count2]; //衰减以及权重函数值

			const double r_comp[6] = {x*x, y*y, z*z, x*y, y*z, z*x}; //等效长程力衰减计算
			const double temp_gmat[6] = {gv*r_comp[0], gv*r_comp[1], gv*r_comp[2], gv*r_comp[3], gv*r_comp[4], gv*r_comp[5]}; //衰减函数矩阵的对称项

			for(int i=0; i<6; i++) Gmat[i] += temp_gmat[i];  //关于count2叠加

			for(int i=0; i<6; i++)
				for(int j=0; j<8; j++) GNmat[i][j] +=  temp_gmat[i]*gauss_ns[count2][j]; //关于count2并与形函数乘积叠加

			//--------------------------------------------
			//计算长程力等效刚度矩阵
			for(int i=0; i<6; i++)
				for(int j=0; j<6; j++)
					Long_range_stiffness[count1][i][j] += temp_gmat[i]*r_comp[j]; //相当于gv*r_comp[i]*r_comp[j]
		}

		//------------------------------------------------------------------------------------------------------------------------
		int ni = 0, nj = 0;
		double temxy = 0, nxy[6] = {0};
		for(int i=0; i<8; i++)
		{
			nj = 0;
			for(int j=0; j<8; j++)
			{
				temxy = gaus_nwct1[i]*gauss_ns[count1][j];
				for(int k=0; k<6; k++)
				nxy[k] = temxy*Gmat[k];

				element_stiff_matrix1[ni][nj] += nxy[0];
				element_stiff_matrix1[ni+1][nj+1] += nxy[1];
				element_stiff_matrix1[ni+2][nj+2] += nxy[2];
				element_stiff_matrix1[ni][nj+1] += nxy[3];
				element_stiff_matrix1[ni+1][nj] += nxy[3];
				element_stiff_matrix1[ni+1][nj+2] += nxy[4];
				element_stiff_matrix1[ni+2][nj+1] += nxy[4];
				element_stiff_matrix1[ni][nj+2] += nxy[5];
				element_stiff_matrix1[ni+2][nj] += nxy[5];

				for(int k=0; k<6; k++)
				nxy[k] = gaus_nwct1[i]*GNmat[k][j];

				element_stiff_matrix2[ni][nj] -= nxy[0];
				element_stiff_matrix2[ni+1][nj+1] -= nxy[1];
				element_stiff_matrix2[ni+2][nj+2] -= nxy[2];
				element_stiff_matrix2[ni][nj+1] -= nxy[3];
				element_stiff_matrix2[ni+1][nj] -= nxy[3];
				element_stiff_matrix2[ni+1][nj+2] -= nxy[4];
				element_stiff_matrix2[ni+2][nj+1] -= nxy[4];
				element_stiff_matrix2[ni][nj+2] -= nxy[5];
				element_stiff_matrix2[ni+2][nj] -= nxy[5];

				nj += 3;
			}
			ni += 3;
		}
	}

	//输出单刚阵，用于检查
	//hout << "element stiffness matrix1: " << endl;
	//for(int i=0; i<24; i++){
	//	for(int j=0; j<24;j++)
	//		hout << element_stiff_matrix1[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl;

	//hout << "element stiffness matrix2: " << endl;
	//for(int i=0; i<24; i++){
	//	for(int j=0; j<24;j++)
	//		hout << element_stiff_matrix2[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl;
}
//-----------------------------------------------------------------------------------------------
//Generate the element matrix based on local continuum model (contact forces)
void Global_Stiff_Matrix::Generate_Contacforce_Elestiff(double (*element_stiff_matrix)[24], const double *alpha_key, const double (*Long_range_stiffness)[6][6], const double (*elenodes)[3], 
																									  const vector<MatPro> &mats, const int &flag, const double &Jacobi, const double (*gauss_dif)[3][8], const vector<double> &weight)const
{
	//--------------------------------------------------
	//Clearing
	for(int i=0; i<24; i++) 
		for (int j=0; j<24; j++) 
			element_stiff_matrix[i][j] = 0;
	//--------------------------------------------------
	//Loop gaussian points
	const int ws = (int)weight.size();
	for(int count=0; count<ws; count++)
	{
		//-----------------------------------------------------------------
		//计算此单元所对应的材料弹性矩阵
		double ele_elas[6][6];		
		if(alpha_key[count]==0.0)  //continuum区
		{
			for(int i=0; i<6; i++)
				for(int j=0; j<6; j++)
					ele_elas[i][j] = mats[0].elas_matrix[i][j];
		}
		else if(alpha_key[count]==1.0) continue; //nonlocal区
		else  //耦合区
		{
			for(int i=0; i<6; i++)
				for(int j=0; j<6; j++)
					ele_elas[i][j] = mats[0].elas_matrix[i][j] - Long_range_stiffness[count][i][j];
		}
		
		//--------------------------------------------------
		//形函数N对gauss[count].x, gauss[count].y, gauss[count].z的偏导矩阵
		double diff[3][8];
		for(int i=0; i<3; i++)
			for(int j=0; j<8; j++)
				diff[i][j] = gauss_dif[count][i][j];

		//--------------------------------------------------
		//J矩阵
		double Jmatrix[3][3];
		//以上两个矩阵的积
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
			{
				Jmatrix[i][j]=0;
				for(int k=0; k<8; k++)
				Jmatrix[i][j] += diff[i][k]*elenodes[k][j];
			}
		//--------------------------------------------------
		//求出J矩阵的行列式
		double J_val;
		if(flag==0) J_val = Jacobi;
		else
		{
			J_val = Jmatrix[0][0]*(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])
						 -Jmatrix[0][1]*(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])
						 +Jmatrix[0][2]*(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0]);
		}

		//----------------------------------------------------
		//求出J矩阵的逆矩阵
		double Jinverse[3][3];
			
		Jinverse[0][0]=(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])/J_val;
		Jinverse[1][1]=(Jmatrix[0][0]*Jmatrix[2][2]-Jmatrix[0][2]*Jmatrix[2][0])/J_val;
		Jinverse[2][2]=(Jmatrix[0][0]*Jmatrix[1][1]-Jmatrix[0][1]*Jmatrix[1][0])/J_val;

		Jinverse[0][1]=-(Jmatrix[0][1]*Jmatrix[2][2]-Jmatrix[0][2]*Jmatrix[2][1])/J_val;
		Jinverse[0][2]=(Jmatrix[0][1]*Jmatrix[1][2]-Jmatrix[0][2]*Jmatrix[1][1])/J_val;

		Jinverse[1][0]=-(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])/J_val;
		Jinverse[1][2]=-(Jmatrix[0][0]*Jmatrix[1][2]-Jmatrix[0][2]*Jmatrix[1][0])/J_val;

		Jinverse[2][0]=(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0])/J_val;
		Jinverse[2][1]=-(Jmatrix[0][0]*Jmatrix[2][1]-Jmatrix[0][1]*Jmatrix[2][0])/J_val;
		//-------------------------------------------------------
		//求出N对x,y,z的偏导
		double diffxy[3][8];
		for(int i=0; i<3; i++)
			for(int j=0; j<8; j++)
			{
				diffxy[i][j]=0;
				for(int k=0; k<3; k++)
					diffxy[i][j] += Jinverse[i][k]*diff[k][j];
			}
		//--------------------------------------------------------
		//求出B矩阵
		double B[6][24];
		for(int i=0; i<6; i++)
			for(int j=0; j<24; j++)
				B[i][j]=0;
		for(int i=0; i<8; i++)
		{
			B[0][i*3+0]=diffxy[0][i];
			B[1][i*3+1]=diffxy[1][i];
			B[2][i*3+2]=diffxy[2][i];
			B[3][i*3+0]=diffxy[1][i];
			B[3][i*3+1]=diffxy[0][i];
			B[4][i*3+1]=diffxy[2][i];
			B[4][i*3+2]=diffxy[1][i];
			B[5][i*3+0]=diffxy[2][i];
			B[5][i*3+2]=diffxy[0][i];
		}
		//--------------------------------------------------------
		//计算B_trans
		double B_trans[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
				B_trans[i][j] = B[j][i];
		//--------------------------------------------------------------------------------------------
		//求出B_trans矩阵与ele_elas矩阵的乘积array1
		double array1[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
			{
				array1[i][j]=0; 
				for(int k=0; k<6; k++)
					array1[i][j] += B_trans[i][k]*ele_elas[k][j];
			}
		//求出array1矩阵与B矩阵的乘积array2
		double array2[24][24];
		for(int i=0; i<24; i++)
			for(int j=0; j<24; j++)
			{
				array2[i][j]=0;
				for(int k=0; k<6; k++)
					array2[i][j] += array1[i][k]*B[k][j];
				element_stiff_matrix[i][j] += array2[i][j]*J_val*weight[count];
			}
	}

	//输出单刚阵，用于检查
	//hout << "element stiffness matrix: " << endl;
	//for(int i=0; i<24; i++){
	//	for(int j=0; j<24; j++)
	//		hout << element_stiff_matrix[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl << endl;
}
//-----------------------------------------------------------------------------------------------
//Generate the element matrix based on local continuum model (contact forces)
void Global_Stiff_Matrix::Generate_Contacforce_Elestiff_Coupled_Damaged(double (*element_stiff_matrix)[24], const double *alpha_key, const double (*Long_range_stiffness)[6][6], const Damage &damages, const bool full_dam_i, 
											  const vector<double> &ele_dam_table, const double (*elenodes)[3], const vector<MatPro> &mats, const int &flag, const int &elemat, const double &Jacobi, const double (*gauss_dif)[3][8], const vector<double> &weight)const
{
	//--------------------------------------------------
	//Clearing
	for(int i=0; i<24; i++) 
		for (int j=0; j<24; j++) 
			element_stiff_matrix[i][j] = 0;
	//--------------------------------------------------
	//Loop gaussian points
	const int ws = (int)weight.size();
	for(int count=0; count<ws; count++)
	{
		//-----------------------------------------------------------------
		//Initializing elastic matrix
		double ele_elas[6][6];
		for(int i=0; i<6; i++) 
			for (int j=0; j<6; j++)
				ele_elas[i][j] = 0.0;

		//-------------------------------------------------
		//Estimate elastic matrix for coupling part
		if(alpha_key[count]==0.0)  //continuum区
		{
			for(int i=0; i<6; i++)
				for(int j=0; j<6; j++)
					ele_elas[i][j] = mats[0].elas_matrix[i][j];
		}
		else if(alpha_key[count]!=1.0) //耦合区
		{
			for(int i=0; i<6; i++)
				for(int j=0; j<6; j++)
					ele_elas[i][j] = mats[0].elas_matrix[i][j] - Long_range_stiffness[count][i][j];
		}
		
		//-------------------------------------------------
		//Estimate elastic matrix for damage part
		if(!full_dam_i) //Include the nonlocal area: alpha_key[count]==1.0.
		{
			//Compute elastic matrix of element for damaged part
			double damage_differ = damages.d_crit - ele_dam_table[count];
			if(damage_differ>Zero)
			{
				for(int i=0; i<6; i++)
					for(int j=0; j<6; j++)
						ele_elas[i][j] += mats[elemat+1].elas_matrix[i][j]*damage_differ;
			}
		}

		//--------------------------------------------------
		//形函数N对gauss[count].x, gauss[count].y, gauss[count].z的偏导矩阵
		double diff[3][8];
		for(int i=0; i<3; i++)
			for(int j=0; j<8; j++)
				diff[i][j] = gauss_dif[count][i][j];

		//--------------------------------------------------
		//J矩阵
		double Jmatrix[3][3];
		//以上两个矩阵的积
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
			{
				Jmatrix[i][j]=0;
				for(int k=0; k<8; k++)
				Jmatrix[i][j] += diff[i][k]*elenodes[k][j];
			}
		//--------------------------------------------------
		//求出J矩阵的行列式
		double J_val;
		if(flag==0) J_val = Jacobi;
		else
		{
			J_val = Jmatrix[0][0]*(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])
						 -Jmatrix[0][1]*(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])
						 +Jmatrix[0][2]*(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0]);
		}

		//----------------------------------------------------
		//求出J矩阵的逆矩阵
		double Jinverse[3][3];
			
		Jinverse[0][0]=(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])/J_val;
		Jinverse[1][1]=(Jmatrix[0][0]*Jmatrix[2][2]-Jmatrix[0][2]*Jmatrix[2][0])/J_val;
		Jinverse[2][2]=(Jmatrix[0][0]*Jmatrix[1][1]-Jmatrix[0][1]*Jmatrix[1][0])/J_val;

		Jinverse[0][1]=-(Jmatrix[0][1]*Jmatrix[2][2]-Jmatrix[0][2]*Jmatrix[2][1])/J_val;
		Jinverse[0][2]=(Jmatrix[0][1]*Jmatrix[1][2]-Jmatrix[0][2]*Jmatrix[1][1])/J_val;

		Jinverse[1][0]=-(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])/J_val;
		Jinverse[1][2]=-(Jmatrix[0][0]*Jmatrix[1][2]-Jmatrix[0][2]*Jmatrix[1][0])/J_val;

		Jinverse[2][0]=(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0])/J_val;
		Jinverse[2][1]=-(Jmatrix[0][0]*Jmatrix[2][1]-Jmatrix[0][1]*Jmatrix[2][0])/J_val;
		//-------------------------------------------------------
		//求出N对x,y,z的偏导
		double diffxy[3][8];
		for(int i=0; i<3; i++)
			for(int j=0; j<8; j++)
			{
				diffxy[i][j]=0;
				for(int k=0; k<3; k++)
					diffxy[i][j] += Jinverse[i][k]*diff[k][j];
			}
		//--------------------------------------------------------
		//求出B矩阵
		double B[6][24];
		for(int i=0; i<6; i++)
			for(int j=0; j<24; j++)
				B[i][j]=0;
		for(int i=0; i<8; i++)
		{
			B[0][i*3+0]=diffxy[0][i];
			B[1][i*3+1]=diffxy[1][i];
			B[2][i*3+2]=diffxy[2][i];
			B[3][i*3+0]=diffxy[1][i];
			B[3][i*3+1]=diffxy[0][i];
			B[4][i*3+1]=diffxy[2][i];
			B[4][i*3+2]=diffxy[1][i];
			B[5][i*3+0]=diffxy[2][i];
			B[5][i*3+2]=diffxy[0][i];
		}
		//--------------------------------------------------------
		//计算B_trans
		double B_trans[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
				B_trans[i][j] = B[j][i];
		//--------------------------------------------------------------------------------------------
		//求出B_trans矩阵与ele_elas矩阵的乘积array1
		double array1[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
			{
				array1[i][j]=0; 
				for(int k=0; k<6; k++)
					array1[i][j] += B_trans[i][k]*ele_elas[k][j];
			}
		//求出array1矩阵与B矩阵的乘积array2
		double array2[24][24];
		for(int i=0; i<24; i++)
			for(int j=0; j<24; j++)
			{
				array2[i][j]=0;
				for(int k=0; k<6; k++)
					array2[i][j] += array1[i][k]*B[k][j];
				element_stiff_matrix[i][j] += array2[i][j]*J_val*weight[count];
			}
	}

	//输出单刚阵，用于检查
	//hout << "element stiffness matrix: " << endl;
	//for(int i=0; i<24; i++){
	//	for(int j=0; j<24; j++)
	//		hout << element_stiff_matrix[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl << endl;
}
//-----------------------------------------------------------------------------------------------
//Assemble the global stiffness matrix
int Global_Stiff_Matrix::Update_global_matrix(const vector<MathMatrix> &ele_self_matrix, const vector<vector<MathMatrix> > &ele_relative_matrix, const vector<Element> &elements, 
																	   const vector<int> &Iz, const vector<int> &Ig, vector<double> &total_matrix)const
{
	const int ES = (int)elements.size();
	for(int i=0; i<ES; i++) 
	{
		Add_to_gsmatrix(ele_self_matrix[i], Iz, Ig, total_matrix, elements[i]);	//由于在openmp中, 所以不能返回0值, 终止程序

		for(int j=0; j<(int)elements[i].relative_eles.size(); j++)
		{
			Add_to_gsmatrix(ele_relative_matrix[i][j], Iz, Ig, total_matrix, elements[i], elements[elements[i].relative_eles[j]]);	//由于在openmp中, 所以不能返回0值, 终止程序
		}
	}

	return 1;
}
//-----------------------------------------------------------------------------------------------
//Add element matrix to globla matrix(in one element)
void Global_Stiff_Matrix::Add_to_gsmatrix(const double (*element_stiff_matrix)[24], const vector<int> &Iz, const vector<int> &Ig, vector<double> &total_matrix, const Element &element)const
{
	for(int i=0; i<(int)element.nodes_id.size(); i++)
		for(int j=0; j<(int)element.nodes_id.size(); j++)
		{
			const int row_node = element.nodes_id[i];
			const int col_node = element.nodes_id[j];
			if(row_node < col_node) continue;
			if(row_node > col_node )
			{
				//二分法查找位置
				bool mark = false;
				int left = Iz[row_node-1];
				int middle = 0;
				int right = Iz[row_node];
				while(right>=left)
				{
					middle = (left + right)/2;
					if(Ig[middle] == col_node) { mark = true; break; } //找到对应节点
					else if(Ig[middle] > col_node) right = middle - 1;
					else left = middle + 1;
				}					
				if(!mark) { hout << "错误！在节点" << row_node << "的Iz, Ig中没有找到节点" << col_node << "请检查！" << endl; }
				
				const int Mt = 6*row_node + 9*middle;	//找到对应一维存储的起始位置

				for(int k=0; k<=2; k++)
					for(int m=0; m<=2; m++)
						total_matrix[Mt+3*k+m] += element_stiff_matrix[3*i+k][3*j+m];
			}
			else
			{
				const int Mt = 6*row_node + 9*Iz[row_node]; //找到对应一维存储的起始位置
				
				for(int k=0; k<=2; k++)
					for(int m=0; m<=k; m++)	//注意此处是小于等于k, 因为是下三角阵
					{
						if(k==2)
							total_matrix[Mt+k+m+1] += element_stiff_matrix[3*i+k][3*j+m];
						else
							total_matrix[Mt+k+m] += element_stiff_matrix[3*i+k][3*j+m];
					}
			}			             
		}
}
//-----------------------------------------------------------------------------------------------
//Add element matrix to globla matrix(in one element)
void Global_Stiff_Matrix::Add_to_gsmatrix(const MathMatrix &element_stiff_matrix, const vector<int> &Iz, const vector<int> &Ig, vector<double> &total_matrix, const Element &element)const
{
	for(int i=0; i<(int)element.nodes_id.size(); i++)
		for(int j=0; j<(int)element.nodes_id.size(); j++)
		{
			const int row_node = element.nodes_id[i];
			const int col_node = element.nodes_id[j];
			if(row_node < col_node) continue;
			if(row_node > col_node )
			{
				//二分法查找位置
				bool mark = false;
				int left = Iz[row_node-1];
				int middle = 0;
				int right = Iz[row_node];
				while(right>=left)
				{
					middle = (left + right)/2;
					if(Ig[middle] == col_node) { mark = true; break; } //找到对应节点
					else if(Ig[middle] > col_node) right = middle - 1;
					else left = middle + 1;
				}					
				if(!mark) { hout << "错误！在节点" << row_node << "的Iz, Ig中没有找到节点" << col_node << "请检查！" << endl; }
				
				const int Mt = 6*row_node + 9*middle;	//找到对应一维存储的起始位置

				for(int k=0; k<=2; k++)
					for(int m=0; m<=2; m++)
						total_matrix[Mt+3*k+m] += element_stiff_matrix.element[3*i+k][3*j+m];
			}
			else
			{
				const int Mt = 6*row_node + 9*Iz[row_node]; //找到对应一维存储的起始位置
				
				for(int k=0; k<=2; k++)
					for(int m=0; m<=k; m++)	//注意此处是小于等于k, 因为是下三角阵
					{
						if(k==2)
							total_matrix[Mt+k+m+1] += element_stiff_matrix.element[3*i+k][3*j+m];
						else
							total_matrix[Mt+k+m] += element_stiff_matrix.element[3*i+k][3*j+m];
					}
			}			             
		}
}
//-----------------------------------------------------------------------------------------------
//Add element matrix to globla matrix(between different elements)
void Global_Stiff_Matrix::Add_to_gsmatrix(const double (*element_stiff_matrix)[24], const vector<int> &Iz, const vector<int> &Ig, vector<double> &total_matrix, const Element &ele_row, const Element &ele_col)const
{
	for(int i=0; i<(int)ele_row.nodes_id.size(); i++)
		for(int j=0; j<(int)ele_col.nodes_id.size(); j++)
		{
			const int row_node = ele_row.nodes_id[i];
			const int col_node = ele_col.nodes_id[j];
			if(row_node < col_node) continue;
			if(row_node > col_node )
			{
				//二分法查找位置
				bool mark = false;
				int left = Iz[row_node-1];
				int middle = 0;
				int right = Iz[row_node];
				while(right>=left)
				{
					middle = (left + right)/2;
					if(Ig[middle] == col_node) { mark = true; break; } //找到对应节点
					else if(Ig[middle] > col_node) right = middle - 1;
					else left = middle + 1;
				}					
				if(!mark) { hout << "错误！在节点" << row_node << "的Iz, Ig中没有找到节点" << col_node << "请检查！" << endl; }
				
				const int Mt = 6*row_node + 9*middle;	//找到对应一维存储的起始位置

				for(int k=0; k<=2; k++)
					for(int m=0; m<=2; m++)
						total_matrix[Mt+3*k+m] += element_stiff_matrix[3*i+k][3*j+m];
			}
			else
			{
				const int Mt = 6*row_node + 9*Iz[row_node]; //找到对应一维存储的起始位置
				
				for(int k=0; k<=2; k++)
					for(int m=0; m<=k; m++)	//注意此处是小于等于k, 因为是下三角阵
					{
						if(k==2)
							total_matrix[Mt+k+m+1] += element_stiff_matrix[3*i+k][3*j+m];
						else
							total_matrix[Mt+k+m] += element_stiff_matrix[3*i+k][3*j+m];
					}
			}			             
		}
}
//-----------------------------------------------------------------------------------------------
//Add element matrix to globla matrix(between different elements)
void Global_Stiff_Matrix::Add_to_gsmatrix(const MathMatrix &element_stiff_matrix, const vector<int> &Iz, const vector<int> &Ig, vector<double> &total_matrix, const Element &ele_row, const Element &ele_col)const
{
	for(int i=0; i<(int)ele_row.nodes_id.size(); i++)
		for(int j=0; j<(int)ele_col.nodes_id.size(); j++)
		{
			const int row_node = ele_row.nodes_id[i];
			const int col_node = ele_col.nodes_id[j];
			if(row_node < col_node) continue;
			if(row_node > col_node )
			{
				//二分法查找位置
				bool mark = false;
				int left = Iz[row_node-1];
				int middle = 0;
				int right = Iz[row_node];
				while(right>=left)
				{
					middle = (left + right)/2;
					if(Ig[middle] == col_node) { mark = true; break; } //找到对应节点
					else if(Ig[middle] > col_node) right = middle - 1;
					else left = middle + 1;
				}					
				if(!mark) { hout << "错误！在节点" << row_node << "的Iz, Ig中没有找到节点" << col_node << "请检查！" << endl; }
				
				const int Mt = 6*row_node + 9*middle;	//找到对应一维存储的起始位置

				for(int k=0; k<=2; k++)
					for(int m=0; m<=2; m++)
						total_matrix[Mt+3*k+m] += element_stiff_matrix.element[3*i+k][3*j+m];
			}
			else
			{
				const int Mt = 6*row_node + 9*Iz[row_node]; //找到对应一维存储的起始位置
				
				for(int k=0; k<=2; k++)
					for(int m=0; m<=k; m++)	//注意此处是小于等于k, 因为是下三角阵
					{
						if(k==2)
							total_matrix[Mt+k+m+1] += element_stiff_matrix.element[3*i+k][3*j+m];
						else
							total_matrix[Mt+k+m] += element_stiff_matrix.element[3*i+k][3*j+m];
					}
			}			             
		}
}
//-----------------------------------------------------------------------------------------------
//Update damaged element matrices
int Global_Stiff_Matrix::Update_damaged_element_matrices(const int &gaussnum, const struct Damage damages, const vector<MatPro> &mats, const vector<Node> &nodes, const vector<Element> &elements, const vector<double> &U_s,
																										   const vector<bool> &full_dam_eles, vector<MathMatrix> &ele_self_matrix, vector<vector<double> > &damage_table, int &damaged_sum, vector<bool> &dam_iter)const
{
	//---------------------------------------------------------------------------
	//Generating gaussian points
	Gauss gau;		//gaussian potins
	if(gau.Generate_gauss(gaussnum)==0) return 0;

	//------------------------------------------------------------------------
	//Loops all elements
	cout << "-_- Updating the element matrices in damage part (local)..." << endl;
	hout << "-_- Updating the element matrices in damage part (local)..." << endl;

	//Implement openmp
	#pragma omp parallel
	{
		const int ES = (int)elements.size();
		double element_stiff_matrix[24][24];
		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<ES; i++) 
		{
			//Output calculating progress
			#pragma omp critical
			if(i%(ES/10)==0)	cout << i*100/ES <<"% done"<<endl;

			//---------------------------------------------------------------------------
			if(full_dam_eles[i]) continue;

			//---------------------------------------------------------------------------
			//Extracting coordinates of nodes in the element 
			double elenodes[8][3];
			for(int j=0; j<8; j++)
			{
				elenodes[j][0] = nodes[elements[i].nodes_id[j]].x;
				elenodes[j][1] = nodes[elements[i].nodes_id[j]].y;
				elenodes[j][2] = nodes[elements[i].nodes_id[j]].z;
			}
			//Extracting displacement solutions of nodes in the element
			double U_ele_nod[24] = { 0.0 };
			for(int j=0; j<8; j++)
			{
				U_ele_nod[3*j] = U_s[3*elements[i].nodes_id[j]];
				U_ele_nod[3*j+1] = U_s[3*elements[i].nodes_id[j]+1];
				U_ele_nod[3*j+2] = U_s[3*elements[i].nodes_id[j]+2];
			}

			//Define a record for updated broken bonds
			int damaged_num=0;
			//-------------------------------------------------------------------
			//updated stiffness matrix of an element
			if(elements[i].type==381&&(int)elements[i].nodes_id.size()==8)	Brick_line_update_damaged(element_stiff_matrix, damage_table[i], damaged_num, damages, elenodes, U_ele_nod, mats, elements[i].mat, gau.gauss, gau.weight); //For a damaged brick element
			else 
			{
				cout << "Error: the type of element or the number of nodes in this element is wrong, when updating stiffness matrix of an element!"  << endl;
				hout << "Error: the type of element or the number of nodes in this element is wrong, when updating stiffness matrix of an element!"  << endl;

			}

			if(damaged_num>0)
			{
				//Updating data
				#pragma omp critical
				{
					damaged_sum += damaged_num;

					dam_iter[i]=true;

					for(int k=0; k<24; k++)
						for(int m=0; m<24; m++)
						{
							ele_self_matrix[i].element[k][m] -= element_stiff_matrix[k][m];
						}
				}
			}
		}
	}

	cout << "^_^ The element matrices updated successfully!" << endl;
	hout << "^_^ The element matrices updated successfully!" << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------
//Update element matrices with broken bonds
int Global_Stiff_Matrix::Update_nonlocal_element_matrices(const int &gaussnum, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func, const vector<MatPro> &mats,
																										   const vector<Node> &nodes, const vector<vector<double> > &Gp_val, const vector<double> &U_s, vector<Element> &elements,
																											vector<MathMatrix> &ele_self_matrix, vector<vector<MathMatrix> > &ele_relative_matrix, vector<bool> &break_table, int &broken_sum, vector<bool> &break_iter)const
{
	//---------------------------------------------------------------------------
	//Generating gaussian points
	Gauss gau;		//gaussian potins
	if(gau.Generate_gauss(gaussnum)==0) return 0;

	//------------------------------------------------------------------------
	//Initialization of variables
	cout << "-_- Initializing variables..." << endl;
	hout << "-_- Initializing variables..." << endl;

	const int GS = (int)gau.gauss.size();
	const int ES = (int)elements.size();
	double (*gauss_ns)[8] = new double [GS][8];				//记录单元高斯点的形函数
	double (*gauss_nw)[8] = new double [GS][8];				//记录单元高斯点的形函数(带权重系数)
	double (*gauss_dif)[3][8] = new double [GS][3][8];		//记录单元高斯点形函数的导数
	double Jacobi = 0;																//标准立方体单元的雅可比值
	double (*gauss_po)[3] = new double [GS][3];				//标准立方体单元的高斯点坐标
	double (*ele_cent)[4] = new double [ES][4];					//记录单元中心点位置向量(x,y,z)分别放于[0],[1]和[2]中（以标准立方体单元中心点为原点）, [3]中放单元顶点权重的平均值
	
	//Calculate data of every gaussian point of element
	Generate_element_gauss_data(elements, nodes, com_mod, weight_func, gau.gauss, gau.weight, gauss_ns, gauss_nw, gauss_dif, Jacobi, gauss_po, ele_cent);
	
	cout << "^_^ Variables initialized successfully!" << endl << endl;
	hout << "^_^ Variables initialized successfully!" << endl << endl;

	//------------------------------------------------------------------------
	//Loops all elements
	cout << "-_- Updating the element matrices in nonlocal part..." << endl;
	hout << "-_- Updating the element matrices in nonlocal part..." << endl;
	//执行openmp
	#pragma omp parallel
	{
		double element_stiff_matrix1[24][24];
		double element_stiff_matrix2[24][24];
		double *alpha_key = new double [GS];
		
		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<ES; i++) 
		{
			//Output calculating progress
			#pragma omp critical
			if(i%(ES/10)==0)	cout << i*100/ES <<"% done"<<endl;
			
			//---------------------------------------------------------------------------
			//Extracting coordinates of nodes in the principal element 
			double elenodes[8][3];
			for(int j=0; j<8; j++)
			{
				elenodes[j][0] = nodes[elements[i].nodes_id[j]].x;
				elenodes[j][1] = nodes[elements[i].nodes_id[j]].y;
				elenodes[j][2] = nodes[elements[i].nodes_id[j]].z;
			}
			//Extracting coordinates of nodes in the updated principal element with displacement solutions
			double new_elenodes[8][3];
			for(int j=0; j<8; j++)
			{
				new_elenodes[j][0] = elenodes[j][0] + U_s[3*elements[i].nodes_id[j]];
				new_elenodes[j][1] = elenodes[j][1] + U_s[3*elements[i].nodes_id[j]+1];
				new_elenodes[j][2] = elenodes[j][2] + U_s[3*elements[i].nodes_id[j]+2];
			}

			//---------------------------------------------------------------------------
			//Calculating hybrid models
			if(com_mod=="Hybrid")
			{	
				//------------------------------------------------------------------------
				//Initialazing alpha weighting coefficients
				for(int j=0; j<GS; j++) alpha_key[j] = -1.0;

				//Defining parameters of the center of the principal element
				double elec_left[4] = { ele_cent[i][0], ele_cent[i][1], ele_cent[i][2], ele_cent[i][3] };

				//------------------------------------------------------------------------
				//Accumulate the number of relative elements (less than i) //attention: Openmp Parallel
				long int ref_tot = 0;
				for(int j=0; j<i; j++) ref_tot += elements[j].relative_eles.size();

				//------------------------------------------------------------------------
				//Loops element neighbours
				for(int j=0; j<(int)elements[i].relative_eles.size(); j++)
				{
					ref_tot++; //Accumlation
					//--------------------------------------------------
					//Defining parameters of the center of the associated element 
					const int ere = elements[i].relative_eles[j]; 
					double elec_right[4] = { ele_cent[ere][0], ele_cent[ere][1], ele_cent[ere][2], ele_cent[ere][3] };

					//--------------------------------------------------
					//A preliminary examination for skipping both elements completely in the local continuum area
					if(elec_left[3]<=Zero&&elec_right[3]<=Zero) 
					{
						for(int count=0; count<GS; count++)
						{
							if(alpha_key[count]==-1.0)  alpha_key[count] = 0.0;
							else if(alpha_key[count]!=0.0) alpha_key[count] = 0.5;
						}
						continue;
					}

					//--------------------------------------------------
					//Extracting coordinates of nodes in the associated element
					double elenodes_relative[8][3];
					for(int k=0; k<8; k++)
					{
						elenodes_relative[k][0] = nodes[elements[ere].nodes_id[k]].x;
						elenodes_relative[k][1] = nodes[elements[ere].nodes_id[k]].y;
						elenodes_relative[k][2] = nodes[elements[ere].nodes_id[k]].z;
					}
					//Extracting coordinates of nodes in the updated associated element with displacement solutions
					//extracting new relative element coordinates
					double new_elenodes_relative[8][3];
					for(int k=0; k<8; k++)
					{
						new_elenodes_relative[k][0] = elenodes_relative[k][0] + U_s[3*elements[ere].nodes_id[k]];
						new_elenodes_relative[k][1] = elenodes_relative[k][1] + U_s[3*elements[ere].nodes_id[k]+1];
						new_elenodes_relative[k][2] = elenodes_relative[k][2] + U_s[3*elements[ere].nodes_id[k]+2];
					}

					//--------------------------------------------------
					//Calculate the position of break_table
					const long int pos_bt = (ref_tot-1)*GS*GS;
					//Define a record for updated broken bonds
					int broken_num=0;
					//------------------------------------------------------------------------
					//Generate the element matrices based on non-local continuum model (long-range forces) judging broken bonds 
					Generate_Longforce_Elestiff_Breaks(element_stiff_matrix1, element_stiff_matrix2, alpha_key, break_table, broken_num, elements[i].Dissipative_energy,
																					   weight_func, elements[i].flag, elements[ere].flag, elenodes, new_elenodes, elenodes_relative, new_elenodes_relative, gauss_ns, gauss_nw,
																					   gauss_dif, gau.weight, peri_para, Jacobi, gauss_po, elec_left, elec_right, pos_bt, Gp_val, i, ere);
					
					if(broken_num>0)
					{
						//Updating data
						#pragma omp critical	
						{
							broken_sum += broken_num;

							break_iter[i]=true;
							break_iter[ere]=true;

							for(int k=0; k<24; k++)
								for(int m=0; m<24; m++)
								{
									ele_self_matrix[i].element[k][m] -= element_stiff_matrix1[k][m];
									ele_relative_matrix[i][j].element[k][m] -= element_stiff_matrix2[k][m];
								}
						}
					}
				}
			}		
		}
		//---------------------------------------------------------------------------
		//Delete pointer
		delete[] alpha_key;
	}

	//---------------------------------------------------------------------------
	//Delete pointers
	delete[] gauss_ns;
	delete[] gauss_nw;
	delete[] gauss_dif;
	delete[] gauss_po;
	delete[] ele_cent;

	cout << "^_^ The element matrices updated successfully!" << endl;
	hout << "^_^ The element matrices updated successfully!" << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------
//Update element matrices with broken bonds (after up to damage criterion)
int Global_Stiff_Matrix::Update_nonlocal_element_matrices(const int &gaussnum, const struct Damage damages, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func, const vector<MatPro> &mats,
																										   const vector<Node> &nodes, const vector<vector<double> > &Gp_val, const vector<double> &U_s, vector<Element> &elements, const vector<vector<double> > &damage_table,
																											vector<MathMatrix> &ele_self_matrix, vector<vector<MathMatrix> > &ele_relative_matrix, vector<bool> &break_table, int &broken_sum, vector<bool> &break_iter)const
{
	//---------------------------------------------------------------------------
	//Generating gaussian points
	Gauss gau;		//gaussian potins
	if(gau.Generate_gauss(gaussnum)==0) return 0;

	//------------------------------------------------------------------------
	//Initialization of variables
	cout << "-_- Initializing variables..." << endl;
	hout << "-_- Initializing variables..." << endl;

	const int GS = (int)gau.gauss.size();
	const int ES = (int)elements.size();
	double (*gauss_ns)[8] = new double [GS][8];				//记录单元高斯点的形函数
	double (*gauss_nw)[8] = new double [GS][8];				//记录单元高斯点的形函数(带权重系数)
	double (*gauss_dif)[3][8] = new double [GS][3][8];		//记录单元高斯点形函数的导数
	double Jacobi = 0;																//标准立方体单元的雅可比值
	double (*gauss_po)[3] = new double [GS][3];				//标准立方体单元的高斯点坐标
	double (*ele_cent)[4] = new double [ES][4];					//记录单元中心点位置向量(x,y,z)分别放于[0],[1]和[2]中（以标准立方体单元中心点为原点）, [3]中放单元顶点权重的平均值
	
	//Calculate data of every gaussian point of element
	Generate_element_gauss_data(elements, nodes, com_mod, weight_func, gau.gauss, gau.weight, gauss_ns, gauss_nw, gauss_dif, Jacobi, gauss_po, ele_cent);
	
	cout << "^_^ Variables initialized successfully!" << endl << endl;
	hout << "^_^ Variables initialized successfully!" << endl << endl;

	//------------------------------------------------------------------------
	//Loops all elements
	cout << "-_- Updating the element matrices in nonlocal part..." << endl;
	hout << "-_- Updating the element matrices in nonlocal part..." << endl;
	//执行openmp
	#pragma omp parallel
	{
		double element_stiff_matrix1[24][24];
		double element_stiff_matrix2[24][24];
		double *alpha_key = new double [GS];
		
		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<ES; i++) 
		{
			//Output calculating progress
			#pragma omp critical
			if(i%(ES/10)==0)	cout << i*100/ES <<"% done"<<endl;
			
			//---------------------------------------------------------------------------
			//Extracting coordinates of nodes in the principal element 
			double elenodes[8][3];
			for(int j=0; j<8; j++)
			{
				elenodes[j][0] = nodes[elements[i].nodes_id[j]].x;
				elenodes[j][1] = nodes[elements[i].nodes_id[j]].y;
				elenodes[j][2] = nodes[elements[i].nodes_id[j]].z;
			}
			//Extracting coordinates of nodes in the updated principal element with displacement solutions
			double new_elenodes[8][3];
			for(int j=0; j<8; j++)
			{
				new_elenodes[j][0] = elenodes[j][0] + U_s[3*elements[i].nodes_id[j]];
				new_elenodes[j][1] = elenodes[j][1] + U_s[3*elements[i].nodes_id[j]+1];
				new_elenodes[j][2] = elenodes[j][2] + U_s[3*elements[i].nodes_id[j]+2];
			}

			//---------------------------------------------------------------------------
			//Calculating hybrid models
			if(com_mod=="Hybrid")
			{	
				//------------------------------------------------------------------------
				//Initialazing alpha weighting coefficients
				for(int j=0; j<GS; j++) alpha_key[j] = -1.0;

				//Defining parameters of the center of the principal element
				double elec_left[4] = { ele_cent[i][0], ele_cent[i][1], ele_cent[i][2], ele_cent[i][3] };

				//------------------------------------------------------------------------
				//Accumulate the number of relative elements (less than i) //attention: Openmp Parallel
				long int ref_tot = 0;
				for(int j=0; j<i; j++) ref_tot += elements[j].relative_eles.size();

				//------------------------------------------------------------------------
				//Loops element neighbours
				for(int j=0; j<(int)elements[i].relative_eles.size(); j++)
				{
					ref_tot++; //Accumlation
					//--------------------------------------------------
					//Defining parameters of the center of the associated element 
					const int ere = elements[i].relative_eles[j]; 
					double elec_right[4] = { ele_cent[ere][0], ele_cent[ere][1], ele_cent[ere][2], ele_cent[ere][3] };

					//--------------------------------------------------
					//A preliminary examination for skipping both elements completely in the local continuum area
					if(elec_left[3]<=Zero&&elec_right[3]<=Zero) 
					{
						for(int count=0; count<GS; count++)
						{
							if(alpha_key[count]==-1.0)  alpha_key[count] = 0.0;
							else if(alpha_key[count]!=0.0) alpha_key[count] = 0.5;
						}
						continue;
					}

					//--------------------------------------------------
					//Extracting coordinates of nodes in the associated element
					double elenodes_relative[8][3];
					for(int k=0; k<8; k++)
					{
						elenodes_relative[k][0] = nodes[elements[ere].nodes_id[k]].x;
						elenodes_relative[k][1] = nodes[elements[ere].nodes_id[k]].y;
						elenodes_relative[k][2] = nodes[elements[ere].nodes_id[k]].z;
					}
					//Extracting coordinates of nodes in the updated associated element with displacement solutions
					//extracting new relative element coordinates
					double new_elenodes_relative[8][3];
					for(int k=0; k<8; k++)
					{
						new_elenodes_relative[k][0] = elenodes_relative[k][0] + U_s[3*elements[ere].nodes_id[k]];
						new_elenodes_relative[k][1] = elenodes_relative[k][1] + U_s[3*elements[ere].nodes_id[k]+1];
						new_elenodes_relative[k][2] = elenodes_relative[k][2] + U_s[3*elements[ere].nodes_id[k]+2];
					}

					//--------------------------------------------------
					//Calculate the position of break_table
					const long int pos_bt = (ref_tot-1)*GS*GS;
					//Define a record for updated broken bonds
					int broken_num=0;
					//------------------------------------------------------------------------
					//Generate the element matrices based on non-local continuum model (long-range forces) judging broken bonds 
					Generate_Longforce_Elestiff_Breaks(element_stiff_matrix1, element_stiff_matrix2, alpha_key, break_table, broken_num, damages, damage_table, elements[i].Dissipative_energy,
																					   weight_func, elements[i].flag, elements[ere].flag, elenodes, new_elenodes, elenodes_relative, new_elenodes_relative, gauss_ns, gauss_nw,
																					   gauss_dif, gau.weight, peri_para, Jacobi, gauss_po, elec_left, elec_right, pos_bt, Gp_val, i, ere);
					
					if(broken_num>0)
					{
						//Updating data
						#pragma omp critical	
						{
							broken_sum += broken_num;

							break_iter[i]=true;
							break_iter[ere]=true;

							for(int k=0; k<24; k++)
								for(int m=0; m<24; m++)
								{
									ele_self_matrix[i].element[k][m] -= element_stiff_matrix1[k][m];
									ele_relative_matrix[i][j].element[k][m] -= element_stiff_matrix2[k][m];
								}
						}
					}
				}
			}		
		}
		//---------------------------------------------------------------------------
		//Delete pointer
		delete[] alpha_key;
	}

	//---------------------------------------------------------------------------
	//Delete pointers
	delete[] gauss_ns;
	delete[] gauss_nw;
	delete[] gauss_dif;
	delete[] gauss_po;
	delete[] ele_cent;

	cout << "^_^ The element matrices updated successfully!" << endl;
	hout << "^_^ The element matrices updated successfully!" << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------
//Generate the element matrix based on non-local continuum model (long-range forces) with broken bonds
void Global_Stiff_Matrix::Generate_Longforce_Elestiff_Breaks(double (*element_stiff_matrix1)[24], double (*element_stiff_matrix2)[24], double *alpha_key, vector<bool> &break_table, int &broken_num, double &Dissip_ener,
																												const struct Weight_func &weight_func, const int &flag_left, const int &flag_right, const double (*elenodes_left)[3], const double(*n_ele_l)[3], const double (*elenodes_right)[3],
																												const double(*n_ele_r)[3], const double (*gauss_ns)[8], const double (*gauss_nw)[8], const double (*gauss_dif)[3][8], const vector<double> &weight, const struct Peri_para &peri_para,
																												const double &Jacobi, const double (*gauss_po)[3], const double elec_left[], const double elec_right[], const long int &pos_bt, const vector<vector<double> > &Gp_val, const int &pri_ele,
																												const int &ass_ele)const
{
	long int cogau = pos_bt;
	//--------------------------------------------------
	//Clearing
	for(int i=0; i<24; i++) 
		for (int j=0; j<24; j++)
		{
			element_stiff_matrix1[i][j] = 0;
			element_stiff_matrix2[i][j] = 0;
		}

	//--------------------------------------------------	
	//Loop the gaussian points in the current element
	const int gau_num = (int)weight.size();
	for(int count1=0; count1<gau_num; count1++)
	{
		//--------------------------------------------------
		int bk_count = 0;
		double Gmat[6] = {0};
		double GNmat[6][8] = { {0}, {0}, {0}, {0}, {0}, {0} };
		//--------------------------------------------------	
		//Coordinates of gaussian point (left)
		Point_3D gaupoi_left(0, 0, 0);
		if(flag_left==0)
		{
			gaupoi_left.x = gauss_po[count1][0] + elec_left[0];
			gaupoi_left.y = gauss_po[count1][1] + elec_left[1];
			gaupoi_left.z = gauss_po[count1][2] + elec_left[2];
		}
		else
		{
			for(int i=0; i<8; i++) 
			{
				gaupoi_left.x += gauss_ns[count1][i]*elenodes_left[i][0];
				gaupoi_left.y += gauss_ns[count1][i]*elenodes_left[i][1];
				gaupoi_left.z += gauss_ns[count1][i]*elenodes_left[i][2];
			}
		}
		//--------------------------------------------------	
		//New coordinates of gaussian point (left)
		Point_3D n_gaupoi_left(0, 0, 0);
		for(int i=0; i<8; i++) 
		{
			n_gaupoi_left.x += gauss_ns[count1][i]*n_ele_l[i][0];
			n_gaupoi_left.y += gauss_ns[count1][i]*n_ele_l[i][1];
			n_gaupoi_left.z += gauss_ns[count1][i]*n_ele_l[i][2];
		}
		//--------------------------------------------------
		//Evaluate values of weighting function at gaussian point (left)
		double weight_left;
		if(elec_left[3]<=Zero) weight_left = 0.0;
		else if(fabs(elec_left[3]-1.0)<=Zero) weight_left = 1.0;
		else weight_left = Gp_val[pri_ele][count1];

		//--------------------------------------------------
		//Evaluate Jacobi matix
		double Jleft_val;
		if(flag_left==0) Jleft_val = Jacobi;
		else
		{
			//Jacobi matix
			double Jmatrix_left[3][3];
			//以上两个矩阵的积
			for(int i=0; i<3; i++)
				for(int j=0; j<3; j++)
				{
					Jmatrix_left[i][j]=0;
					for(int k=0; k<8; k++)
						Jmatrix_left[i][j] += gauss_dif[count1][i][k]*elenodes_left[k][j];
				}

			//求出J矩阵的行列式
			Jleft_val = Jmatrix_left[0][0]*(Jmatrix_left[1][1]*Jmatrix_left[2][2]-Jmatrix_left[1][2]*Jmatrix_left[2][1])
								-Jmatrix_left[0][1]*(Jmatrix_left[1][0]*Jmatrix_left[2][2]-Jmatrix_left[1][2]*Jmatrix_left[2][0])
								+Jmatrix_left[0][2]*(Jmatrix_left[1][0]*Jmatrix_left[2][1]-Jmatrix_left[1][1]*Jmatrix_left[2][0]);	
		}
		
		//--------------------------------------------------
		//Evaluate shape-function values with weighting
		double gaus_nwct1[8];
		for(int i=0; i<8; i++) gaus_nwct1[i] = gauss_nw[count1][i]*Jleft_val;

		//------------------------------------------------------------------------------------------------------------------------
		//iteration on second ghost point
		for(int count2=0; count2<gau_num; count2++)
		{
			//--------------------------------------------------
			if(break_table[cogau++]) continue;  //First use and then ++, same to: { break_table[cogau] continue; cogau++; } 
			//--------------------------------------------
			//Coordinates of gaussian point (right)
			Point_3D gaupoi_right(0, 0, 0);
			if(flag_right==0)
			{
				gaupoi_right.x = gauss_po[count2][0] + elec_right[0];
				gaupoi_right.y = gauss_po[count2][1] + elec_right[1];
				gaupoi_right.z = gauss_po[count2][2] + elec_right[2];
			}
			else
			{
				for(int i=0; i<8; i++) 
				{
					gaupoi_right.x +=  gauss_ns[count2][i]*elenodes_right[i][0];
					gaupoi_right.y +=  gauss_ns[count2][i]*elenodes_right[i][1];
					gaupoi_right.z +=  gauss_ns[count2][i]*elenodes_right[i][2];
				}
			}
			//--------------------------------------------
			//New coordinates of gaussian point (right)
			Point_3D n_gaupoi_right(0, 0, 0);
			for(int i=0; i<8; i++) 
			{
				n_gaupoi_right.x +=  gauss_ns[count2][i]*n_ele_r[i][0];
				n_gaupoi_right.y +=  gauss_ns[count2][i]*n_ele_r[i][1];
				n_gaupoi_right.z +=  gauss_ns[count2][i]*n_ele_r[i][2];
			}

			//--------------------------------------------
			//Evaluate values of weighting function at gaussian point (right)
			double weight_right;
			if(elec_right[3]<=Zero) weight_right = 0.0;
			else if(fabs(elec_right[3]-1.0)<=Zero) weight_right = 1.0;
			else weight_right =Gp_val[ass_ele][count2];

			//--------------------------------------------
			//计算权重函数值
			double alpha = 0.5*(weight_left+weight_right);

			if(alpha<1.0&&alpha>0.0) alpha_key[count1] = 0.5;
			else if(alpha_key[count1]==-1.0)  alpha_key[count1] = alpha;
			else if(alpha_key[count1]!=alpha) alpha_key[count1] = 0.5;

			if(alpha<=Zero) continue;  //等于零值，没有长程效果

			//--------------------------------------------
			//计算长程作用衰减函数值
			const double x = gaupoi_right.x-gaupoi_left.x;
			const double y = gaupoi_right.y-gaupoi_left.y;
			const double z = gaupoi_right.z-gaupoi_left.z;

			const double dis_squr = x*x + y*y + z*z;
			const double poi_dis = sqrt(dis_squr);		//两点之间的距离

			if(poi_dis>peri_para.horizon_R||poi_dis<Zero) continue;		//圆形积分域

			//-------------------------------------------------
			//Evaluation: does the bond broke?
			const double n_x = n_gaupoi_right.x-n_gaupoi_left.x;
			const double n_y = n_gaupoi_right.y-n_gaupoi_left.y;
			const double n_z = n_gaupoi_right.z-n_gaupoi_left.z;

			const double n_dis_squr = n_x*n_x + n_y*n_y + n_z*n_z;
			const double n_poi_dis = sqrt(n_dis_squr);		//两点之间的距离

			if(n_poi_dis<peri_para.broken_factor*poi_dis) continue;

			//updating data
			long int btco = cogau - 1;
			break_table[btco]=true;
			bk_count++;

			//-------------------------------------------------
			//计算Ｊ矩阵；
			double Jright_val;
			if(flag_right==0) Jright_val = Jacobi;
			else
			{
				//--------------------------------------------------
				//J矩阵
				double Jmatrix_right[3][3];
				//以上两个矩阵的积
				for(int i=0; i<3; i++)
					for(int j=0; j<3; j++)
					{
						Jmatrix_right[i][j]=0;
						for(int k=0; k<8; k++)
							Jmatrix_right[i][j]=Jmatrix_right[i][j] + gauss_dif[count2][i][k]*elenodes_right[k][j];
					}

					//求出J矩阵的行列式
					Jright_val = Jmatrix_right[0][0]*(Jmatrix_right[1][1]*Jmatrix_right[2][2]-Jmatrix_right[1][2]*Jmatrix_right[2][1])
						-Jmatrix_right[0][1]*(Jmatrix_right[1][0]*Jmatrix_right[2][2]-Jmatrix_right[1][2]*Jmatrix_right[2][0])
						+Jmatrix_right[0][2]*(Jmatrix_right[1][0]*Jmatrix_right[2][1]-Jmatrix_right[1][1]*Jmatrix_right[2][0]);
			}

			//--------------------------------------------
			//Calculate the matrix contribution
			const double cos2sita = z*z/dis_squr;		//cos(sita)^2
			double cos2pha;
			if(fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
			else cos2pha = x*x/(x*x + y*y);					//cos(pha)^2

			double sum = peri_para.acoe[0] + peri_para.acoe[1]*0.5*(3*cos2sita-1) + peri_para.acoe[2]*(2*cos2pha-1)*3*(1-cos2sita) + peri_para.acoe[3]*0.125*(35*cos2sita*cos2sita - 30*cos2sita + 3)
									+ peri_para.acoe[4]*(2*cos2pha-1)*7.5*(7*cos2sita-1)*(1-cos2sita) + peri_para.acoe[5]*(8*cos2pha*cos2pha-8*cos2pha+1)*105*(1-cos2sita)*(1-cos2sita);

			const double gv = exp(-poi_dis/peri_para.intrinsic_L)*sum*alpha*Jright_val*weight[count2]; //衰减以及权重函数值
			Dissip_ener += exp(-poi_dis/peri_para.intrinsic_L)*sum*alpha*weight[count2]*(peri_para.broken_factor-1.0)*(peri_para.broken_factor-1.0)*0.25*weight[count1]*dis_squr*dis_squr*Jright_val*Jleft_val;

			const double r_comp[6] = {x*x, y*y, z*z, x*y, y*z, z*x}; //等效长程力衰减计算
			const double temp_gmat[6] = {gv*r_comp[0], gv*r_comp[1], gv*r_comp[2], gv*r_comp[3], gv*r_comp[4], gv*r_comp[5]}; //衰减函数矩阵的对称项

			for(int i=0; i<6; i++) Gmat[i] += temp_gmat[i];  //关于count2叠加

			for(int i=0; i<6; i++)
				for(int j=0; j<8; j++) GNmat[i][j] +=  temp_gmat[i]*gauss_ns[count2][j]; //关于count2并与形函数乘积叠加
		}

		if(bk_count>0)
		{
			broken_num += bk_count;
			//------------------------------------------------------------------------------------------------------------------------
			int ni = 0, nj = 0;
			double temxy = 0, nxy[6] = {0};
			for(int i=0; i<8; i++)
			{
				nj = 0;
				for(int j=0; j<8; j++)
				{
					temxy = gaus_nwct1[i]*gauss_ns[count1][j];
					for(int k=0; k<6; k++)
						nxy[k] = temxy*Gmat[k];

					element_stiff_matrix1[ni][nj] += nxy[0];
					element_stiff_matrix1[ni+1][nj+1] += nxy[1];
					element_stiff_matrix1[ni+2][nj+2] += nxy[2];
					element_stiff_matrix1[ni][nj+1] += nxy[3];
					element_stiff_matrix1[ni+1][nj] += nxy[3];
					element_stiff_matrix1[ni+1][nj+2] += nxy[4];
					element_stiff_matrix1[ni+2][nj+1] += nxy[4];
					element_stiff_matrix1[ni][nj+2] += nxy[5];
					element_stiff_matrix1[ni+2][nj] += nxy[5];

					for(int k=0; k<6; k++)
						nxy[k] = gaus_nwct1[i]*GNmat[k][j];

					element_stiff_matrix2[ni][nj] -= nxy[0];
					element_stiff_matrix2[ni+1][nj+1] -= nxy[1];
					element_stiff_matrix2[ni+2][nj+2] -= nxy[2];
					element_stiff_matrix2[ni][nj+1] -= nxy[3];
					element_stiff_matrix2[ni+1][nj] -= nxy[3];
					element_stiff_matrix2[ni+1][nj+2] -= nxy[4];
					element_stiff_matrix2[ni+2][nj+1] -= nxy[4];
					element_stiff_matrix2[ni][nj+2] -= nxy[5];
					element_stiff_matrix2[ni+2][nj] -= nxy[5];

					nj += 3;
				}
				ni += 3;
			}
		}
	}
}
//-----------------------------------------------------------------------------------------------
//Generate the element matrix based on non-local continuum model (long-range forces) with broken bonds (after up to damage criterion)
void Global_Stiff_Matrix::Generate_Longforce_Elestiff_Breaks(double (*element_stiff_matrix1)[24], double (*element_stiff_matrix2)[24], double *alpha_key, vector<bool> &break_table, int &broken_num, const struct Damage damages,  const vector<vector<double> > &damage_table,
																												double &Dissip_ener, const struct Weight_func &weight_func, const int &flag_left, const int &flag_right, const double (*elenodes_left)[3], const double(*n_ele_l)[3], const double (*elenodes_right)[3],
																												const double(*n_ele_r)[3], const double (*gauss_ns)[8], const double (*gauss_nw)[8], const double (*gauss_dif)[3][8], const vector<double> &weight, const struct Peri_para &peri_para,
																												const double &Jacobi, const double (*gauss_po)[3], const double elec_left[], const double elec_right[], const long int &pos_bt, const vector<vector<double> > &Gp_val, const int &pri_ele,
																												const int &ass_ele)const
{
	long int cogau = pos_bt;
	//--------------------------------------------------
	//Clearing
	for(int i=0; i<24; i++) 
		for (int j=0; j<24; j++)
		{
			element_stiff_matrix1[i][j] = 0;
			element_stiff_matrix2[i][j] = 0;
		}

	//--------------------------------------------------	
	//Loop the gaussian points in the current element
	const int gau_num = (int)weight.size();
	for(int count1=0; count1<gau_num; count1++)
	{
		//--------------------------------------------------
		//The bond can break until up to damage criterion
		if(damages.d_crit-damage_table[pri_ele][count1]>Zero) continue;

		//--------------------------------------------------
		int bk_count = 0;
		double Gmat[6] = {0};
		double GNmat[6][8] = { {0}, {0}, {0}, {0}, {0}, {0} };
		//--------------------------------------------------	
		//Coordinates of gaussian point (left)
		Point_3D gaupoi_left(0, 0, 0);
		if(flag_left==0)
		{
			gaupoi_left.x = gauss_po[count1][0] + elec_left[0];
			gaupoi_left.y = gauss_po[count1][1] + elec_left[1];
			gaupoi_left.z = gauss_po[count1][2] + elec_left[2];
		}
		else
		{
			for(int i=0; i<8; i++) 
			{
				gaupoi_left.x += gauss_ns[count1][i]*elenodes_left[i][0];
				gaupoi_left.y += gauss_ns[count1][i]*elenodes_left[i][1];
				gaupoi_left.z += gauss_ns[count1][i]*elenodes_left[i][2];
			}
		}
		//--------------------------------------------------	
		//New coordinates of gaussian point (left)
		Point_3D n_gaupoi_left(0, 0, 0);
		for(int i=0; i<8; i++) 
		{
			n_gaupoi_left.x += gauss_ns[count1][i]*n_ele_l[i][0];
			n_gaupoi_left.y += gauss_ns[count1][i]*n_ele_l[i][1];
			n_gaupoi_left.z += gauss_ns[count1][i]*n_ele_l[i][2];
		}
		//--------------------------------------------------
		//Evaluate values of weighting function at gaussian point (left)
		double weight_left;
		if(elec_left[3]<=Zero) weight_left = 0.0;
		else if(fabs(elec_left[3]-1.0)<=Zero) weight_left = 1.0;
		else weight_left = Gp_val[pri_ele][count1];

		//--------------------------------------------------
		//Evaluate Jacobi matix
		double Jleft_val;
		if(flag_left==0) Jleft_val = Jacobi;
		else
		{
			//Jacobi matix
			double Jmatrix_left[3][3];
			//以上两个矩阵的积
			for(int i=0; i<3; i++)
				for(int j=0; j<3; j++)
				{
					Jmatrix_left[i][j]=0;
					for(int k=0; k<8; k++)
						Jmatrix_left[i][j] += gauss_dif[count1][i][k]*elenodes_left[k][j];
				}

			//求出J矩阵的行列式
			Jleft_val = Jmatrix_left[0][0]*(Jmatrix_left[1][1]*Jmatrix_left[2][2]-Jmatrix_left[1][2]*Jmatrix_left[2][1])
								-Jmatrix_left[0][1]*(Jmatrix_left[1][0]*Jmatrix_left[2][2]-Jmatrix_left[1][2]*Jmatrix_left[2][0])
								+Jmatrix_left[0][2]*(Jmatrix_left[1][0]*Jmatrix_left[2][1]-Jmatrix_left[1][1]*Jmatrix_left[2][0]);	
		}
		
		//--------------------------------------------------
		//Evaluate shape-function values with weighting
		double gaus_nwct1[8];
		for(int i=0; i<8; i++) gaus_nwct1[i] = gauss_nw[count1][i]*Jleft_val;

		//------------------------------------------------------------------------------------------------------------------------
		//iteration on second ghost point
		for(int count2=0; count2<gau_num; count2++)
		{
			//--------------------------------------------------
			//The bond can break until up to damage criterion
			if(damages.d_crit-damage_table[ass_ele][count2]>Zero) continue;			

			//--------------------------------------------------
			if(break_table[cogau++]) continue;  //First use and then ++, same to: { break_table[cogau] continue; cogau++; } 
			//--------------------------------------------
			//Coordinates of gaussian point (right)
			Point_3D gaupoi_right(0, 0, 0);
			if(flag_right==0)
			{
				gaupoi_right.x = gauss_po[count2][0] + elec_right[0];
				gaupoi_right.y = gauss_po[count2][1] + elec_right[1];
				gaupoi_right.z = gauss_po[count2][2] + elec_right[2];
			}
			else
			{
				for(int i=0; i<8; i++) 
				{
					gaupoi_right.x +=  gauss_ns[count2][i]*elenodes_right[i][0];
					gaupoi_right.y +=  gauss_ns[count2][i]*elenodes_right[i][1];
					gaupoi_right.z +=  gauss_ns[count2][i]*elenodes_right[i][2];
				}
			}
			//--------------------------------------------
			//New coordinates of gaussian point (right)
			Point_3D n_gaupoi_right(0, 0, 0);
			for(int i=0; i<8; i++) 
			{
				n_gaupoi_right.x +=  gauss_ns[count2][i]*n_ele_r[i][0];
				n_gaupoi_right.y +=  gauss_ns[count2][i]*n_ele_r[i][1];
				n_gaupoi_right.z +=  gauss_ns[count2][i]*n_ele_r[i][2];
			}

			//--------------------------------------------
			//Evaluate values of weighting function at gaussian point (right)
			double weight_right;
			if(elec_right[3]<=Zero) weight_right = 0.0;
			else if(fabs(elec_right[3]-1.0)<=Zero) weight_right = 1.0;
			else weight_right = Gp_val[ass_ele][count2];

			//--------------------------------------------
			//计算权重函数值
			double alpha = 0.5*(weight_left+weight_right);

			if(alpha<1.0&&alpha>0.0) alpha_key[count1] = 0.5;
			else if(alpha_key[count1]==-1.0)  alpha_key[count1] = alpha;
			else if(alpha_key[count1]!=alpha) alpha_key[count1] = 0.5;

			if(alpha<=Zero) continue;  //等于零值，没有长程效果

			//--------------------------------------------
			//计算长程作用衰减函数值
			const double x = gaupoi_right.x-gaupoi_left.x;
			const double y = gaupoi_right.y-gaupoi_left.y;
			const double z = gaupoi_right.z-gaupoi_left.z;

			const double dis_squr = x*x + y*y + z*z;
			const double poi_dis = sqrt(dis_squr);		//两点之间的距离

			if(poi_dis>peri_para.horizon_R||poi_dis<Zero) continue;		//圆形积分域

			//-------------------------------------------------
			//Evaluation: does the bond broke?
			const double n_x = n_gaupoi_right.x-n_gaupoi_left.x;
			const double n_y = n_gaupoi_right.y-n_gaupoi_left.y;
			const double n_z = n_gaupoi_right.z-n_gaupoi_left.z;

			const double n_dis_squr = n_x*n_x + n_y*n_y + n_z*n_z;
			const double n_poi_dis = sqrt(n_dis_squr);		//两点之间的距离

			if(n_poi_dis<peri_para.broken_factor*poi_dis) continue;

			//updating data
			long int btco = cogau - 1;
			break_table[btco]=true;
			bk_count++;

			//-------------------------------------------------
			//计算Ｊ矩阵；
			double Jright_val;
			if(flag_right==0) Jright_val = Jacobi;
			else
			{
				//--------------------------------------------------
				//J矩阵
				double Jmatrix_right[3][3];
				//以上两个矩阵的积
				for(int i=0; i<3; i++)
					for(int j=0; j<3; j++)
					{
						Jmatrix_right[i][j]=0;
						for(int k=0; k<8; k++)
							Jmatrix_right[i][j]=Jmatrix_right[i][j] + gauss_dif[count2][i][k]*elenodes_right[k][j];
					}

					//求出J矩阵的行列式
					Jright_val = Jmatrix_right[0][0]*(Jmatrix_right[1][1]*Jmatrix_right[2][2]-Jmatrix_right[1][2]*Jmatrix_right[2][1])
						-Jmatrix_right[0][1]*(Jmatrix_right[1][0]*Jmatrix_right[2][2]-Jmatrix_right[1][2]*Jmatrix_right[2][0])
						+Jmatrix_right[0][2]*(Jmatrix_right[1][0]*Jmatrix_right[2][1]-Jmatrix_right[1][1]*Jmatrix_right[2][0]);
			}

			//--------------------------------------------
			//Calculate the matrix contribution
			const double cos2sita = z*z/dis_squr;		//cos(sita)^2
			double cos2pha;
			if(fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
			else cos2pha = x*x/(x*x + y*y);					//cos(pha)^2

			double sum = peri_para.acoe[0] + peri_para.acoe[1]*0.5*(3*cos2sita-1) + peri_para.acoe[2]*(2*cos2pha-1)*3*(1-cos2sita) + peri_para.acoe[3]*0.125*(35*cos2sita*cos2sita - 30*cos2sita + 3)
									+ peri_para.acoe[4]*(2*cos2pha-1)*7.5*(7*cos2sita-1)*(1-cos2sita) + peri_para.acoe[5]*(8*cos2pha*cos2pha-8*cos2pha+1)*105*(1-cos2sita)*(1-cos2sita);

			const double gv = exp(-poi_dis/peri_para.intrinsic_L)*sum*alpha*Jright_val*weight[count2]; //衰减以及权重函数值
			Dissip_ener += exp(-poi_dis/peri_para.intrinsic_L)*sum*alpha*weight[count2]*(peri_para.broken_factor-1.0)*(peri_para.broken_factor-1.0)*0.25*weight[count1]*dis_squr*dis_squr*Jright_val*Jleft_val;

			const double r_comp[6] = {x*x, y*y, z*z, x*y, y*z, z*x}; //等效长程力衰减计算
			const double temp_gmat[6] = {gv*r_comp[0], gv*r_comp[1], gv*r_comp[2], gv*r_comp[3], gv*r_comp[4], gv*r_comp[5]}; //衰减函数矩阵的对称项

			for(int i=0; i<6; i++) Gmat[i] += temp_gmat[i];  //关于count2叠加

			for(int i=0; i<6; i++)
				for(int j=0; j<8; j++) GNmat[i][j] +=  temp_gmat[i]*gauss_ns[count2][j]; //关于count2并与形函数乘积叠加
		}

		if(bk_count>0)
		{
			broken_num += bk_count;
			//------------------------------------------------------------------------------------------------------------------------
			int ni = 0, nj = 0;
			double temxy = 0, nxy[6] = {0};
			for(int i=0; i<8; i++)
			{
				nj = 0;
				for(int j=0; j<8; j++)
				{
					temxy = gaus_nwct1[i]*gauss_ns[count1][j];
					for(int k=0; k<6; k++)
						nxy[k] = temxy*Gmat[k];

					element_stiff_matrix1[ni][nj] += nxy[0];
					element_stiff_matrix1[ni+1][nj+1] += nxy[1];
					element_stiff_matrix1[ni+2][nj+2] += nxy[2];
					element_stiff_matrix1[ni][nj+1] += nxy[3];
					element_stiff_matrix1[ni+1][nj] += nxy[3];
					element_stiff_matrix1[ni+1][nj+2] += nxy[4];
					element_stiff_matrix1[ni+2][nj+1] += nxy[4];
					element_stiff_matrix1[ni][nj+2] += nxy[5];
					element_stiff_matrix1[ni+2][nj] += nxy[5];

					for(int k=0; k<6; k++)
						nxy[k] = gaus_nwct1[i]*GNmat[k][j];

					element_stiff_matrix2[ni][nj] -= nxy[0];
					element_stiff_matrix2[ni+1][nj+1] -= nxy[1];
					element_stiff_matrix2[ni+2][nj+2] -= nxy[2];
					element_stiff_matrix2[ni][nj+1] -= nxy[3];
					element_stiff_matrix2[ni+1][nj] -= nxy[3];
					element_stiff_matrix2[ni+1][nj+2] -= nxy[4];
					element_stiff_matrix2[ni+2][nj+1] -= nxy[4];
					element_stiff_matrix2[ni][nj+2] -= nxy[5];
					element_stiff_matrix2[ni+2][nj] -= nxy[5];

					nj += 3;
				}
				ni += 3;
			}
		}
	}
}
//-----------------------------------------------------------------------------------------------
//Update the damaged brick element
void Global_Stiff_Matrix::Brick_line_update_damaged(double (*element_stiff_matrix)[24], vector<double> &ele_dam_table, int &damaged_num, const struct Damage damages, const double (*elenodes)[3],
																								const double *U_ele_nod, const vector<MatPro> &mats, const int &elemat, const vector<Node> &gauss, const vector<double> &wight)const
{
	double Ymax = 0.0;
	//--------------------------------------------------
	//Clearing
	for(int i=0; i<24; i++) 
		for (int j=0; j<24; j++)
			element_stiff_matrix[i][j] = 0.0;

	//--------------------------------------------------
	//循环高斯点计算积分
	for(int count=0; count<(int)wight.size(); count++)
	{
		if(damages.d_crit-ele_dam_table[count]<=Zero) continue; //Jump over full damaged gaussian points

		//-----------------------------------------------------------------
		//计算此单元所对应的材料弹性矩阵
		double ele_elas[6][6];		
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
				ele_elas[i][j] = mats[elemat+1].elas_matrix[i][j];
		
		//计算Ｊ矩阵
		//--------------------------------------------
		//形函数N对gauss[count].x, gauss[count].y, gauss[count].z的偏导矩阵
		double diff[3][8];
		diff[0][0]=-0.125*(1.0-gauss[count].y)*(1.0-gauss[count].z);
		diff[0][1]=-diff[0][0];                         
		diff[0][2]=0.125*(1.0+gauss[count].y)*(1.0-gauss[count].z);
		diff[0][3]=-diff[0][2];
		diff[0][4]=-0.125*(1.0-gauss[count].y)*(1.0+gauss[count].z);
		diff[0][5]=-diff[0][4];
		diff[0][6]=0.125*(1.0+gauss[count].y)*(1.0+gauss[count].z);
		diff[0][7]=-diff[0][6];

		diff[1][0]=-0.125*(1.0-gauss[count].x)*(1.0-gauss[count].z);
		diff[1][1]=-0.125*(1.0+gauss[count].x)*(1.0-gauss[count].z);
		diff[1][2]=-diff[1][1];
		diff[1][3]=-diff[1][0];
		diff[1][4]=-0.125*(1.0-gauss[count].x)*(1.0+gauss[count].z);
		diff[1][5]=-0.125*(1.0+gauss[count].x)*(1.0+gauss[count].z);
		diff[1][6]=-diff[1][5];
		diff[1][7]=-diff[1][4];

		diff[2][0]=-0.125*(1.0-gauss[count].x)*(1.0-gauss[count].y);
		diff[2][1]=-0.125*(1.0+gauss[count].x)*(1.0-gauss[count].y);
		diff[2][2]=-0.125*(1.0+gauss[count].x)*(1.0+gauss[count].y);
		diff[2][3]=-0.125*(1.0-gauss[count].x)*(1.0+gauss[count].y);
		diff[2][4]=-diff[2][0];
		diff[2][5]=-diff[2][1];
		diff[2][6]=-diff[2][2];
		diff[2][7]=-diff[2][3];

		//--------------------------------------------------
		//J矩阵
		double Jmatrix[3][3];
		//以上两个矩阵的积
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
			{
				Jmatrix[i][j]=0;
				for(int k=0; k<8; k++)
				Jmatrix[i][j] += diff[i][k]*elenodes[k][j];
			}
		//--------------------------------------------------
		//求出J矩阵的行列式
		double J_val = Jmatrix[0][0]*(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])
								 -Jmatrix[0][1]*(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])
								 +Jmatrix[0][2]*(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0]);
		//----------------------------------------------------
		//求出J矩阵的逆矩阵
		double Jinverse[3][3];
			
		Jinverse[0][0]=(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])/J_val;
		Jinverse[1][1]=(Jmatrix[0][0]*Jmatrix[2][2]-Jmatrix[0][2]*Jmatrix[2][0])/J_val;
		Jinverse[2][2]=(Jmatrix[0][0]*Jmatrix[1][1]-Jmatrix[0][1]*Jmatrix[1][0])/J_val;

		Jinverse[0][1]=-(Jmatrix[0][1]*Jmatrix[2][2]-Jmatrix[0][2]*Jmatrix[2][1])/J_val;
		Jinverse[0][2]=(Jmatrix[0][1]*Jmatrix[1][2]-Jmatrix[0][2]*Jmatrix[1][1])/J_val;

		Jinverse[1][0]=-(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])/J_val;
		Jinverse[1][2]=-(Jmatrix[0][0]*Jmatrix[1][2]-Jmatrix[0][2]*Jmatrix[1][0])/J_val;

		Jinverse[2][0]=(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0])/J_val;
		Jinverse[2][1]=-(Jmatrix[0][0]*Jmatrix[2][1]-Jmatrix[0][1]*Jmatrix[2][0])/J_val;
		//-------------------------------------------------------
		//求出N对x,y,z的偏导
		double diffxy[3][8];
		for(int i=0; i<3; i++)
			for(int j=0; j<8; j++)
			{
				diffxy[i][j]=0;
				for(int k=0; k<3; k++)
					diffxy[i][j] += Jinverse[i][k]*diff[k][j];
			}
		//--------------------------------------------------------
		//求出B矩阵
		double B[6][24];
		for(int i=0; i<6; i++)
			for(int j=0; j<24; j++)
				B[i][j]=0;
		for(int i=0; i<8; i++)
		{
			B[0][i*3+0]=diffxy[0][i];
			B[1][i*3+1]=diffxy[1][i];
			B[2][i*3+2]=diffxy[2][i];
			B[3][i*3+0]=diffxy[1][i];
			B[3][i*3+1]=diffxy[0][i];
			B[4][i*3+1]=diffxy[2][i];
			B[4][i*3+2]=diffxy[1][i];
			B[5][i*3+0]=diffxy[2][i];
			B[5][i*3+2]=diffxy[0][i];
		}
		//--------------------------------------------------------
		//计算B_trans
		double B_trans[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
				B_trans[i][j] = B[j][i];
		//--------------------------------------------------------------------------------------------
		//求出B_trans矩阵与ele_elas矩阵的乘积array1
		double array1[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
			{
				array1[i][j]=0; 
				for(int k=0; k<6; k++)
					array1[i][j] += B_trans[i][k]*ele_elas[k][j];
			}
		//--------------------------------------------------------
		double Y_dam_force = 0.0;  //Damge force by energy decrease due to damage
		double tem[24][24];
		//求出array1矩阵与B矩阵的乘积array2
		double array2[24][24];
		for(int i=0; i<24; i++)
			for(int j=0; j<24; j++)
			{
				array2[i][j]=0;
				for(int k=0; k<6; k++)
					array2[i][j] += array1[i][k]*B[k][j];
				tem[i][j] = array2[i][j]*J_val*wight[count];
				Y_dam_force += 0.5*U_ele_nod[i]*array2[i][j]*U_ele_nod[j];
			}

if(Ymax<Y_dam_force) Ymax=Y_dam_force;

double K1=0.0;
double K0=0.0;
if(elemat==0) 
{
//	K0 = damages.k0*100;
	K0 = damages.k0;
	K1 = damages.k1;
}
else if(elemat==1)
{
	K0 = damages.k0*0.8;
	K1 = damages.k1*0.8;
}

		if((Y_dam_force - K1*ele_dam_table[count] - K0)>=Zero)
		{
			double new_dam = Y_dam_force/K1;
			if(damages.d_crit-new_dam<=Zero) new_dam = damages.d_crit;
			for(int i=0; i<24; i++)
				for(int j=0; j<24; j++)
					element_stiff_matrix[i][j] += tem[i][j]*(new_dam-ele_dam_table[count]);
			ele_dam_table[count] = new_dam;
			damaged_num ++;
		}
	}
//cout << "Ymax= " << Ymax << " EleMat= " << elemat << endl;
//hout << "Ymax= " << Ymax << " EleMat= " << elemat << endl;
	//输出单刚阵，用于检查
	//hout << "element stiffness matrix: " << endl;
	//for(int i=0; i<24; i++){
	//	for(int j=0; j<24; j++)
	//		hout << element_stiff_matrix[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl << endl;

}
//---------------------------------------------------------------------------
//Estimate element matrices with local damage for 2D problem
int Global_Stiff_Matrix::Gen_element_matrices_damage_2D(const int &gaussnum, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func, const struct Damage &damages,
											const vector<MatPro> &mats, vector<Node> &nodes, const vector<Element> &elements, const vector<bool> &full_dam_eles, const vector<vector<double> > &damage_table,
											const vector<bool> &break_table, vector<vector<double> > &Gp_val, vector<vector<double> > &Ak_val, vector<vector<double> > &Ele_local_stiff, vector<MathMatrix> &ele_self_matrix, 
											vector<vector<MathMatrix> > &ele_relative_matrix)const
{
	//---------------------------------------------------------------------------
	//Generating gaussian points
	Gauss gau;		//gaussian potins
	if(gau.Generate_gauss_2D(gaussnum)==0) return 0;

	//------------------------------------------------------------------------
	//Initialization of variables
	cout << "    -_- Initializing variables..." << endl;
	hout << "    -_- Initializing variables..." << endl;

	const int GS = (int)gau.gauss.size();
	const int ES = (int)elements.size();
	double (*gauss_ns)[4] = new double [GS][4];				//记录单元高斯点的形函数
	double (*gauss_nw)[4] = new double [GS][4];				//记录单元高斯点的形函数(带权重系数)
	double (*gauss_dif)[2][4] = new double [GS][2][4];		//记录单元高斯点形函数的导数
	double Jacobi = 0;																//标准正方形单元的雅可比值
	double (*gauss_po)[3] = new double [GS][3];				//标准正方形单元的高斯点坐标
	double (*ele_cent)[4] = new double [ES][4];					//记录单元中心点位置向量(x,y,z)分别放于[0],[1]和[2]中（以标准正方形单元中心点为原点）, [3]中放单元顶点权重的平均值

	//Calculate data of every gaussian point of element for 2D problem
	Generate_element_gauss_data_2D(elements, nodes, com_mod, weight_func, gau.gauss, gau.weight, gauss_ns, gauss_nw, gauss_dif, Jacobi, gauss_po, ele_cent);

	WeightFunc Wei_Fun;
	//Calculate the values of weighting function of gaussian points in all elements
	Wei_Fun.Value_weightfunc_gausspois_elements_2D(weight_func, nodes, elements, GS, gauss_ns, gauss_po, ele_cent, Gp_val);

	cout << "    ^_^ Variables initialized successfully!" << endl << endl;
	hout << "    ^_^ Variables initialized successfully!" << endl << endl;

	//------------------------------------------------------------------------
	//Loops all elements
	cout << "    -_- Calculating the element matrices..." << endl;
	hout << "    -_- Calculating the element matrices..." << endl;
	//执行openmp
	#pragma omp parallel
	{
		double element_stiff_matrix1[12][12];
		double element_stiff_matrix2[12][12];
		double (*Long_range_stiffness)[6][6] = new double [GS][6][6];
		double *alpha_key = new double [GS];
		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<ES; i++) 
		{
			//Output calculating progress
			#pragma omp critical
			if(i%(ES/10)==0)	cout << i*100/ES << "% done" << endl;
			
			//Extracting coordinates of nodes in the principal element 
			double elenodes[4][2];	
			for(int j=0; j<4; j++)
			{
				elenodes[j][0] = nodes[elements[i].nodes_id[j]].x;
				elenodes[j][1] = nodes[elements[i].nodes_id[j]].y;
			}
			//---------------------------------------------------------------------------
			//Calculating hybrid models
			if(com_mod=="Hybrid"&&!(weight_func.num==1&&weight_func.shape[0]=="Null"))
			{
				//Time markers
				//clock_t ctn1,ctn2;
				//ctn1 = clock();
				
				//------------------------------------------------------------------------
				//Initialazing long-range stiffness tensor
				for(int j=0; j<GS; j++)
				{
					for(int k=0; k<6; k++)
						for(int m=0; m<6; m++)
							Long_range_stiffness[j][k][m] = 0.0;
					alpha_key[j] = -1.0;
				}

				//------------------------------------------------------------------------
				//Defining parameters of the center of the principal element 
				double elec_left[4] = { ele_cent[i][0], ele_cent[i][1], ele_cent[i][2], ele_cent[i][3] };

				//------------------------------------------------------------------------
				//Accumulate the number of relative elements (less than i) //attention: Openmp Parallel
				long int ref_tot = 0;
				for(int j=0; j<i; j++) ref_tot += elements[j].relative_eles.size();

				//------------------------------------------------------------------------
				//Loops element neighbours
				for(int j=0; j<(int)elements[i].relative_eles.size(); j++)
				{
					ref_tot++; //Accumlation (注意不能移动到下面 因为以下有continue语句)
					//--------------------------------------------------
					//Defining parameters of the center of the associated element 
					const int ere = elements[i].relative_eles[j]; 
					double elec_right[4] = { ele_cent[ere][0], ele_cent[ere][1], ele_cent[ere][2], ele_cent[ere][3] };
					
					//--------------------------------------------------
					//A preliminary examination for skipping both elements completely in the local continuum area 
					//（注意elec_left[3]和elec_right[3]中分别是两个四边形顶点权重值得平均）
					if(elec_left[3]<=Zero&&elec_right[3]<=Zero) 
					{
						for(int count=0; count<GS; count++)
						{
							if(alpha_key[count]==-1.0)  alpha_key[count] = 0.0;
							else if(alpha_key[count]!=0.0) alpha_key[count] = 0.5;
						}
						continue;
					}

					//--------------------------------------------------
					//Extracting coordinates of nodes in the associated element
					double elenodes_relative[4][2];
					for(int k=0; k<4; k++)
					{
						elenodes_relative[k][0] = nodes[elements[ere].nodes_id[k]].x;
						elenodes_relative[k][1] = nodes[elements[ere].nodes_id[k]].y;
					}

					//--------------------------------------------------
					//Calculate the position of break_table
					const long int pos_bt = (ref_tot-1)*GS*GS;

					//------------------------------------------------------------------------
					//Generate the element matrices based on non-local continuum model (long-range forces)
					Generate_Longforce_Elestiff_Brokeless_2D(element_stiff_matrix1, element_stiff_matrix2, alpha_key, Long_range_stiffness, full_dam_eles[i], elements[i].flag, elements[ere].flag, elenodes,
																						elenodes_relative, gauss_ns, gauss_nw, gauss_dif, gau.weight, peri_para, Jacobi, gauss_po, elec_left, elec_right, break_table, pos_bt, Gp_val, i, ere);

					#pragma omp critical
					{
						//---------------------------------------------------------------------------
						//Accumulate the new element matrices
						for(int k=0; k<12; k++)
							for(int m=0; m<12; m++)
							{
								ele_self_matrix[i].element[k][m] += element_stiff_matrix1[k][m];
								ele_relative_matrix[i][j].element[k][m]=element_stiff_matrix2[k][m];
							}
					}
				}
				//Attention: a half of parameters
				for(int j=0; j<GS; j++)
					for(int k=0; k<6; k++)
						for(int m=0; m<6; m++)
							Long_range_stiffness[j][k][m] = Long_range_stiffness[j][k][m]/2.0;

				//------------------------------------------------------------------------
				//ctn2 = clock();
				//hout << "Total num of elements: " << (int)elements.size() << "; Element " << i << " took time: " << (double)(ctn2-ctn1)/CLOCKS_PER_SEC << "sec; " << endl;
			}
			else
			{
				//------------------------------------------------------------------------
				//Initialazing long-range stiffness tensor
				for(int j=0; j<GS; j++)
				{
					for(int k=0; k<6; k++)
						for(int m=0; m<6; m++)
							Long_range_stiffness[j][k][m] = 0.0;
					alpha_key[j] = 0.0;
				}
			}

			//Generate the element matrix based on local continuum model (coupled part + damaged part)
			vector<double> ele_dam_table(damage_table[i]);
			double(*temp_elas)[6][6] = new double[GS][6][6];
			Generate_Contacforce_Elestiff_Coupled_Damaged_2D(element_stiff_matrix1, alpha_key, Long_range_stiffness, temp_elas, damages, full_dam_eles[i], ele_dam_table, elenodes, mats, elements[i].flag, elements[i].mat, Jacobi, gauss_dif, gau.weight);

			#pragma omp critical
			{
				for(int k=0; k<12; k++)
					for(int m=0; m<12; m++)
						ele_self_matrix[i].element[k][m] += element_stiff_matrix1[k][m];
				
				for(int k=0; k<GS; k++) Ak_val[i][k] = alpha_key[k];
				
				int counter = 0;
				for(int k=0; k<GS; k++)
					for(int m=0; m<6; m++)
						for(int n=0; n<6; n++)
							Ele_local_stiff[i][counter++] = temp_elas[k][m][n];
			}
		}
		//---------------------------------------------------------------------------
		//Delete pointers
		delete[] Long_range_stiffness;
		delete[] alpha_key;
	}

	//---------------------------------------------------------------------------
	//Delete pointers
	delete[] gauss_ns;
	delete[] gauss_nw;
	delete[] gauss_dif;
	delete[] gauss_po;
	delete[] ele_cent;

	cout << "^_^ The element matrices calculated successfully!" << endl;
	hout << "^_^ The element matrices calculated successfully!" << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------
//Calculate data of every gaussian point of element for 2D problem
void Global_Stiff_Matrix::Generate_element_gauss_data_2D(const vector<Element> &elements, const vector<Node> &nodes, const string &com_mod, const struct Weight_func &weight_func, const vector<Node> &gauss, 
																						    const vector<double> &weight, double (*gauss_ns)[4], double (*gauss_nw)[4], double (*gauss_dif)[2][4], double &Jacobi,  double (*gauss_po)[3], double (*ele_cent)[4])const
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Element loop
	int inum = -1;
	for(int i=0; i<(int)elements.size(); i++)
	{
		if(elements[i].flag==0) //Standard cubic element
		{
			inum = i;
			break;
		}
	}

	for(int count=0; count<(int)gauss.size(); count++)
	{	
		double Nshape[4] = {0};
		if(com_mod=="Hybrid") 
		{
			//--------------------------------------------
			//计算高斯积分点整体坐标
			Nshape[0]=0.25*(1.0-gauss[count].x)*(1.0-gauss[count].y);
			Nshape[1]=0.25*(1.0+gauss[count].x)*(1.0-gauss[count].y);
			Nshape[2]=0.25*(1.0+gauss[count].x)*(1.0+gauss[count].y);
			Nshape[3]=0.25*(1.0-gauss[count].x)*(1.0+gauss[count].y);

			for(int j=0; j<4; j++) 
			{
				gauss_ns[count][j] = Nshape[j]; 
				gauss_nw[count][j] = Nshape[j]*weight[count];
			} 				
			//--------------------------------------------
			//计算高斯积分点坐标
			Point_3D gaupoi(0, 0, 0);		
			if(inum>=0)
			{
				for(int j=0; j<4; j++) 
				{
					gaupoi.x += Nshape[j]*nodes[elements[inum].nodes_id[j]].x;
					gaupoi.y += Nshape[j]*nodes[elements[inum].nodes_id[j]].y;
				}
			}

			//记录高斯积分点坐标
			gauss_po[count][0] = gaupoi.x;
			gauss_po[count][1] = gaupoi.y;
			gauss_po[count][2] = gaupoi.z;
		}
		else
		{
			//初始化
			gauss_po[count][0] = 0.0;
			gauss_po[count][1] = 0.0;
			gauss_po[count][2] = 0.0;
		}

		//--------------------------------------------
		//计算Ｊ矩阵
		//--------------------------------------------
		//形函数N对gauss[count].x, gauss[count].y的偏导矩阵
		double diff[2][4];
		diff[0][0]=-0.25*(1.0-gauss[count].y);
		diff[0][1]=-diff[0][0];
		diff[0][2]=0.25*(1.0+gauss[count].y);
		diff[0][3]=-diff[0][2];

		diff[1][0]=-0.25*(1.0-gauss[count].x);
		diff[1][1]=-0.25*(1.0+gauss[count].x);
		diff[1][2]=-diff[1][1];
		diff[1][3]=-diff[1][0];

		//记录每个高斯点形函数偏导数的值
		for(int j=0; j<2; j++)
			for(int k=0; k<4; k++)
				gauss_dif[count][j][k] = diff[j][k];

		if(inum>=0)
		{
			//--------------------------------------------------
			//单元节点坐标矩阵
			double elenode[4][2];
			for(int j=0; j<4; j++)
			{
				elenode[j][0]=nodes[elements[inum].nodes_id[j]].x;
				elenode[j][1]=nodes[elements[inum].nodes_id[j]].y;
			}
			//--------------------------------------------------
			//J矩阵
			double Jmatrix[2][2];
			//以上两个矩阵的积
			for(int j=0; j<2; j++)
				for(int k=0; k<2; k++)
				{
					Jmatrix[j][k]=0;
					for(int m=0; m<4; m++)
					Jmatrix[j][k] += diff[j][m]*elenode[m][k];
				}
			//--------------------------------------------------
			//求出J矩阵的行列式
			Jacobi = Jmatrix[0][0]*Jmatrix[1][1]-Jmatrix[0][1]*Jmatrix[1][0];
		}
		else 
		{
			Jacobi = 0.0;
		}
	}

	//The coordinates of center
	Point_3D scenter(0,0,0);	//Standard cubic element
	if(inum>=0)
	{
		for(int j=0; j<4; j++) 
		{
			scenter.x += nodes[elements[inum].nodes_id[j]].x;
			scenter.y += nodes[elements[inum].nodes_id[j]].y;
		}
		scenter = scenter/4;
	}
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Element loop
	for(int i=0; i<(int)elements.size(); i++)
	{
		for(int j=0; j<4; j++) ele_cent[i][j] = 0.0;
		//计算该单元的中心点坐标
		for(int j=0; j<4; j++) 
		{
			ele_cent[i][0] += nodes[elements[i].nodes_id[j]].x;
			ele_cent[i][1] += nodes[elements[i].nodes_id[j]].y;
			const Point_3D poi_nod(nodes[elements[i].nodes_id[j]].x, nodes[elements[i].nodes_id[j]].y, nodes[elements[i].nodes_id[j]].z);
			WeightFunc Wei_Fun;  //the class of weighting function
			ele_cent[i][3] += Wei_Fun.Value_weight_function(poi_nod, weight_func);  //注意是四个顶点的权重求和
		}
		for(int j=0; j<4; j++) ele_cent[i][j] = ele_cent[i][j]/4;
		ele_cent[i][0] = ele_cent[i][0] - scenter.x;
		ele_cent[i][1] = ele_cent[i][1] - scenter.y;
		ele_cent[i][2] = ele_cent[i][2] - scenter.z;
	}
}
//-----------------------------------------------------------------------------------------------
//Generate the element matrix based on non-local continuum model (long-range forces) without broken bonds for 2D problem
void Global_Stiff_Matrix::Generate_Longforce_Elestiff_Brokeless_2D(double(*element_stiff_matrix1)[12], double(*element_stiff_matrix2)[12], double *alpha_key, double(*Long_range_stiffness)[6][6], const bool full_dam_i,
												const int &flag_left, const int &flag_right, const double (*elenodes_left)[2], const double (*elenodes_right)[2], const double (*gauss_ns)[4], const double (*gauss_nw)[4],
												const double (*gauss_dif)[2][4], const vector<double> &weight, const struct Peri_para &peri_para, const double &Jacobi, const double (*gauss_po)[3], const double elec_left[],
												const double elec_right[], const vector<bool> &break_table, const long int &pos_bt, const vector<vector<double> > &Gp_val, const int &pri_ele, const int &ass_ele)const
{
	long int cogau = pos_bt;
	//--------------------------------------------------
	//Clearing
	for(int i=0; i<12; i++) 
		for (int j=0; j<12; j++)
		{
			element_stiff_matrix1[i][j] = 0;
			element_stiff_matrix2[i][j] = 0;
		}

	//--------------------------------------------------	
	//Loop the gaussian points in the current element
	const int gau_num = (int)weight.size();
	for(int count1=0; count1<gau_num; count1++)
	{
		//--------------------------------------------------	
		double Gmat[6] = {0};
		double GNmat[6][4] = { {0}, {0}, {0}, {0}, {0}, {0} };

		//--------------------------------------------------	
		//左端高斯点坐标
		Point_3D gaupoi_left(0, 0, 0);
		if(flag_left==0)
		{
			gaupoi_left.x = gauss_po[count1][0] + elec_left[0];
			gaupoi_left.y = gauss_po[count1][1] + elec_left[1];
		}
		else
		{
			for(int i=0; i<4; i++) 
			{
				gaupoi_left.x += gauss_ns[count1][i]*elenodes_left[i][0];
				gaupoi_left.y += gauss_ns[count1][i]*elenodes_left[i][1];
			}
		}

		//计算Ｊ矩阵；
		//--------------------------------------------------
		double Jleft_val;
		if(flag_left==0) Jleft_val = Jacobi;
		else
		{
			//J矩阵
			double Jmatrix_left[2][2];
			//以上两个矩阵的积
			for(int i=0; i<2; i++)
				for(int j=0; j<2; j++)
				{
					Jmatrix_left[i][j]=0;
					for(int k=0; k<4; k++)
						Jmatrix_left[i][j] += gauss_dif[count1][i][k]*elenodes_left[k][j];
				}

			//求出J矩阵的行列式
			Jleft_val = Jmatrix_left[0][0]*Jmatrix_left[1][1]-Jmatrix_left[0][1]*Jmatrix_left[1][0];	
		}
		
		//--------------------------------------------------
		//左端带权重形函数
		double gaus_nwct1[4];
		for(int i=0; i<4; i++) gaus_nwct1[i] = gauss_nw[count1][i]*Jleft_val;

		//------------------------------------------------------------------------------------------------------------------------
		//循环外单元高斯点计算积分
		for(int count2=0; count2<gau_num; count2++)
		{
			//--------------------------------------------
			//右端高斯点坐标
			Point_3D gaupoi_right(0, 0, 0);
			if(flag_right==0)
			{
				gaupoi_right.x = gauss_po[count2][0] + elec_right[0];
				gaupoi_right.y = gauss_po[count2][1] + elec_right[1];
			}
			else
			{
				for(int i=0; i<4; i++) 
				{
					gaupoi_right.x +=  gauss_ns[count2][i]*elenodes_right[i][0];
					gaupoi_right.y +=  gauss_ns[count2][i]*elenodes_right[i][1];
				}
			}

			//--------------------------------------------
			//计算长程作用衰减函数值
			const double x = gaupoi_right.x-gaupoi_left.x;
			const double y = gaupoi_right.y-gaupoi_left.y;

			const double dis_squr = x*x + y*y;
			const double poi_dis = sqrt(dis_squr);		//两点之间的距离
			
			if(poi_dis>peri_para.horizon_R||poi_dis<Zero) { cogau++; continue; }		//圆形积分域 //注意在continue之前 bond number 要++

			//--------------------------------------------
			//计算权重函数值(注意这段程序应该在计算完bond的长度之后, 因为长度太长活太短都不属于alpha考虑的范围)
			double alpha = 0.5*(Gp_val[pri_ele][count1] + Gp_val[ass_ele][count2]);
			if (fabs(alpha - 1.0) <= Zero) alpha = 1.0;  //为下面精确赋值
			else if (alpha <= Zero) alpha = 0.0;

			if (alpha<1.0&&alpha>0.0) alpha_key[count1] = 0.5;
			else if (alpha_key[count1] == -1.0)  alpha_key[count1] = alpha;
			else if (alpha_key[count1] != alpha) alpha_key[count1] = 0.5;

			if (alpha <= Zero) { cogau++; continue; }  //等于零值，没有长程效果 //注意在continue之前 bond number 要++

			//bond断裂后本来是不用计算了, 但是为了计算Long_range_stiffness, 还是要算的除非单元完全损伤可省去这部分
			if (full_dam_i&&break_table[cogau]){ cogau++; continue; }

			//--------------------------------------------
			//计算Ｊ矩阵；
			double Jright_val;
			if(flag_right==0) Jright_val = Jacobi;
			else
			{
				//--------------------------------------------------
				//J矩阵
				double Jmatrix_right[2][2];
				//以上两个矩阵的积
				for(int i=0; i<2; i++)
					for(int j=0; j<2; j++)
					{
						Jmatrix_right[i][j]=0;
						for(int k=0; k<4; k++)
							Jmatrix_right[i][j]=Jmatrix_right[i][j] + gauss_dif[count2][i][k]*elenodes_right[k][j];
					}

				//求出J矩阵的行列式
				Jright_val = Jmatrix_right[0][0]*Jmatrix_right[1][1]-Jmatrix_right[0][1]*Jmatrix_right[1][0];
			}
	
			//--------------------------------------------
			//计算长程力矩阵
			double cos2pha;
			if(fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
			else cos2pha = x*x/(x*x + y*y);					//cos(pha)^2

			double sum = peri_para.acoe[0] + peri_para.acoe[1]*(2*cos2pha-1) + peri_para.acoe[2]*(8*cos2pha*cos2pha-8*cos2pha+1);

			const double gv = exp(-poi_dis / peri_para.intrinsic_L)*sum*alpha*Jright_val*weight[count2]/dis_squr; //衰减以及权重函数值
//			const double gv = sum*alpha*Jright_val*weight[count2]/dis_squr; //constant micro-modulus

			const double r_comp[6] = {x*x, y*y, 0, x*y, y*0, 0*x}; //等效长程力衰减计算
			const double temp_gmat[6] = {gv*r_comp[0], gv*r_comp[1], gv*r_comp[2], gv*r_comp[3], gv*r_comp[4], gv*r_comp[5]}; //衰减函数矩阵的对称项

			//--------------------------------------------
			//计算长程力等效刚度矩阵
			for (int i = 0; i<6; i++)
				for (int j = 0; j<6; j++)
					Long_range_stiffness[count1][i][j] += temp_gmat[i] * r_comp[j]; //相当于gv*r_comp[i]*r_comp[j]

			//--------------------------------------------
			//如果bond已经断裂，则不需要累加如下结果
			if (break_table[cogau++]) continue; //First use and then ++, same to: { break_table[cogau] continue; cogau++; } 

			//--------------------------------------------
			for(int i=0; i<6; i++) Gmat[i] += temp_gmat[i];  //关于count2叠加

			for(int i=0; i<6; i++)
				for(int j=0; j<4; j++) GNmat[i][j] +=  temp_gmat[i]*gauss_ns[count2][j]; //关于count2并与形函数乘积叠加
		}

		//------------------------------------------------------------------------------------------------------------------------
		int ni = 0, nj = 0;
		double temxy = 0, nxy[6] = {0};
		for(int i=0; i<4; i++)
		{
			nj = 0;
			for(int j=0; j<4; j++)
			{
				temxy = gaus_nwct1[i]*gauss_ns[count1][j];
				for(int k=0; k<6; k++)
				nxy[k] = temxy*Gmat[k];

				element_stiff_matrix1[ni][nj] += nxy[0];
				element_stiff_matrix1[ni+1][nj+1] += nxy[1];
				element_stiff_matrix1[ni+2][nj+2] += nxy[2];
				element_stiff_matrix1[ni][nj+1] += nxy[3];
				element_stiff_matrix1[ni+1][nj] += nxy[3];
				element_stiff_matrix1[ni+1][nj+2] += nxy[4];
				element_stiff_matrix1[ni+2][nj+1] += nxy[4];
				element_stiff_matrix1[ni][nj+2] += nxy[5];
				element_stiff_matrix1[ni+2][nj] += nxy[5];

				for(int k=0; k<6; k++)
				nxy[k] = gaus_nwct1[i]*GNmat[k][j];

				element_stiff_matrix2[ni][nj] -= nxy[0];
				element_stiff_matrix2[ni+1][nj+1] -= nxy[1];
				element_stiff_matrix2[ni+2][nj+2] -= nxy[2];
				element_stiff_matrix2[ni][nj+1] -= nxy[3];
				element_stiff_matrix2[ni+1][nj] -= nxy[3];
				element_stiff_matrix2[ni+1][nj+2] -= nxy[4];
				element_stiff_matrix2[ni+2][nj+1] -= nxy[4];
				element_stiff_matrix2[ni][nj+2] -= nxy[5];
				element_stiff_matrix2[ni+2][nj] -= nxy[5];

				nj += 3;
			}
			ni += 3;
		}
	}

	//输出单刚阵，用于检查
	//hout << "element stiffness matrix1: " << endl;
	//for(int i=0; i<12; i++){
	//	for(int j=0; j<12;j++)
	//		hout << element_stiff_matrix1[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl;

	//hout << "element stiffness matrix2: " << endl;
	//for(int i=0; i<12; i++){
	//	for(int j=0; j<12;j++)
	//		hout << element_stiff_matrix2[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl;
}
//-----------------------------------------------------------------------------------------------
//Generate the element matrix based on local continuum model (contact forces) for 2D problem
void Global_Stiff_Matrix::Generate_Contacforce_Elestiff_Coupled_Damaged_2D(double (*element_stiff_matrix)[12], const double *alpha_key, const double (*Long_range_stiffness)[6][6], double (*temp_elas)[6][6], const Damage &damages, const bool full_dam_i, 
											  const vector<double> &ele_dam_table, const double (*elenodes)[2], const vector<MatPro> &mats, const int &flag, const int &elemat, const double &Jacobi, const double (*gauss_dif)[2][4], const vector<double> &weight)const
{
	//--------------------------------------------------
	//Initialization
	for(int i=0; i<12; i++) 
		for (int j=0; j<12; j++) 
			element_stiff_matrix[i][j] = 0.0;

	//--------------------------------------------------
	if(full_dam_i)  //Fully damaged element
	{
		const int ws = (int)weight.size();
		for(int count=0; count<ws; count++)
		{
			if (alpha_key[count] != 1.0) hout << "Attention: the alpha_key[" << count << "]!=1.0, in the full damaged element!" << endl;
			for(int i=0; i<6; i++)
				for(int j=0; j<6; j++)
					temp_elas[count][i][j] =  0.0;	//Initialization
		}
	}
	else
	{
		//--------------------------------------------------
		//Loop gaussian points
		const int ws = (int)weight.size();
		for(int count=0; count<ws; count++)
		{	
			//-----------------------------------------------------------------
			//Initializing elastic matrix
			double ele_elas[6][6];
			for(int i=0; i<6; i++) 
				for (int j=0; j<6; j++)
					ele_elas[i][j] = 0.0;

			//-------------------------------------------------
			//Compute elastic matrix of element without damaged part
			double damage_differ = 1.0 - ele_dam_table[count];
			if(damage_differ>Zero)
			{
				for(int i=0; i<6; i++)
					for(int j=0; j<6; j++)
						ele_elas[i][j] = mats[elemat+1].elas_matrix[i][j]*damage_differ;
			}

			//Remove the nonlocal part
			for(int i=0; i<6; i++)
				for(int j=0; j<6; j++)
				{
					ele_elas[i][j] -= Long_range_stiffness[count][i][j];
					temp_elas[count][i][j] =  mats[elemat+1].elas_matrix[i][j] - Long_range_stiffness[count][i][j];
				}	

			ele_elas[2][2] = mats[elemat+1].elas_matrix[2][2];
			ele_elas[4][4] = mats[elemat+1].elas_matrix[4][4];
			ele_elas[5][5] = mats[elemat+1].elas_matrix[5][5];

			//--------------------------------------------------
			//形函数N对gauss[count].x, gauss[count].y的偏导矩阵
			double diff[2][4];
			for(int i=0; i<2; i++)
				for(int j=0; j<4; j++)
					diff[i][j] = gauss_dif[count][i][j];

			//--------------------------------------------------
			//J矩阵
			double Jmatrix[2][2];
			//以上两个矩阵的积
			for(int i=0; i<2; i++)
				for(int j=0; j<2; j++)
				{
					Jmatrix[i][j]=0;
					for(int k=0; k<4; k++)
					Jmatrix[i][j] += diff[i][k]*elenodes[k][j];
				}
			//--------------------------------------------------
			//求出J矩阵的行列式
			double J_val;
			if(flag==0) J_val = Jacobi;
			else
			{
				J_val = Jmatrix[0][0]*Jmatrix[1][1]-Jmatrix[0][1]*Jmatrix[1][0];
			}

			//----------------------------------------------------
			//求出J矩阵的逆矩阵
			double Jinverse[2][2];
			Jinverse[0][0]=Jmatrix[1][1]/J_val;
			Jinverse[1][1]=Jmatrix[0][0]/J_val;
			Jinverse[0][1]=-Jmatrix[0][1]/J_val;
			Jinverse[1][0]=-Jmatrix[1][0]/J_val;

			//-------------------------------------------------------
			//求出N对x,y的偏导
			double diffxy[2][4];
			for(int i=0; i<2; i++)
				for(int j=0; j<4; j++)
				{
					diffxy[i][j]=0;
					for(int k=0; k<2; k++)
						diffxy[i][j] += Jinverse[i][k]*diff[k][j];
				}
			//--------------------------------------------------------
			//求出B矩阵
			double B[6][12];
			for(int i=0; i<6; i++)
				for(int j=0; j<12; j++)
					B[i][j]=0;
			for(int i=0; i<4; i++)
			{
				B[0][i*3+0]=diffxy[0][i];
				B[1][i*3+1]=diffxy[1][i];
				B[3][i*3+0]=diffxy[1][i];
				B[3][i*3+1]=diffxy[0][i];
				B[4][i*3+2]=diffxy[1][i];
				B[5][i*3+2]=diffxy[0][i];
			}
			//--------------------------------------------------------
			//计算B_trans
			double B_trans[12][6];
			for(int i=0; i<12; i++)
				for(int j=0; j<6; j++)
					B_trans[i][j] = B[j][i];
			//--------------------------------------------------------------------------------------------
			//求出B_trans矩阵与ele_elas矩阵的乘积array1
			double array1[12][6];
			for(int i=0; i<12; i++)
				for(int j=0; j<6; j++)
				{
					array1[i][j]=0; 
					for(int k=0; k<6; k++)
						array1[i][j] += B_trans[i][k]*ele_elas[k][j];
				}
			//求出array1矩阵与B矩阵的乘积array2
			double array2[12][12];
			for(int i=0; i<12; i++)
				for(int j=0; j<12; j++)
				{
					array2[i][j]=0;
					for(int k=0; k<6; k++)
						array2[i][j] += array1[i][k]*B[k][j];
					element_stiff_matrix[i][j] += array2[i][j]*J_val*weight[count];
				}
		}

		//输出单刚阵，用于检查
		//hout << "element stiffness matrix: " << endl;
		//for(int i=0; i<12; i++){
		//	for(int j=0; j<12; j++)
		//		hout << element_stiff_matrix[i][j] << " ";
		//		hout << endl;
		//}
		//hout << endl << endl << endl;
	}
}
//-----------------------------------------------------------------------------------------------
//(Constant stress assumption in a horizon) generate the element matrix based on non-local continuum model (long-range forces) without broken bonds for 2D problem
void Global_Stiff_Matrix::Generate_Longforce_Elestiff_Brokeless_2D_Stress_Assump(double(*element_stiff_matrix1)[12], double(*element_stiff_matrix2)[12], double *alpha_key, 
	double(*Long_range_stiffness)[6][6], const struct Weight_func &weight_func, const int &flag_left, const int &flag_right, const double(*elenodes_left)[2], const double(*elenodes_right)[2], 
	const double(*gauss_ns)[4], const double(*gauss_nw)[4], const double(*gauss_dif)[2][4], const vector<double> &weight, const struct Peri_para &peri_para, const double &Jacobi, 
	const double(*gauss_po)[3], const double elec_left[], const double elec_right[], const vector<vector<double> > &damage_table, const vector<bool> &break_table, const long int &pos_bt, 
	const vector<vector<double> > &Gp_val, const int &pri_ele, const int &ass_ele)const
{
	long int cogau = pos_bt;
	//--------------------------------------------------
	//Clearing
	for (int i = 0; i<12; i++)
	for (int j = 0; j<12; j++)
	{
		element_stiff_matrix1[i][j] = 0;
		element_stiff_matrix2[i][j] = 0;
	}

	//--------------------------------------------------	
	//Loop the gaussian points in the current element
	const int gau_num = (int)weight.size();
	for (int count1 = 0; count1<gau_num; count1++)
	{
		//--------------------------------------------------	
		double Gmat[6] = { 0 };
		double GNmat[6][4] = { { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 } };
		//--------------------------------------------------	
		//左端高斯点坐标
		Point_3D gaupoi_left(0, 0, 0);
		if (flag_left == 0)
		{
			gaupoi_left.x = gauss_po[count1][0] + elec_left[0];
			gaupoi_left.y = gauss_po[count1][1] + elec_left[1];
		}
		else
		{
			for (int i = 0; i<4; i++)
			{
				gaupoi_left.x += gauss_ns[count1][i] * elenodes_left[i][0];
				gaupoi_left.y += gauss_ns[count1][i] * elenodes_left[i][1];
			}
		}

		//--------------------------------------------
		//计算左端高斯点的权重
		double weight_left;
		if (elec_left[3] <= Zero) weight_left = 0.0;
		else if (fabs(elec_left[3] - 1.0) <= Zero) weight_left = 1.0;
		else weight_left = Gp_val[pri_ele][count1];

		//计算Ｊ矩阵；
		//--------------------------------------------------
		double Jleft_val;
		if (flag_left == 0) Jleft_val = Jacobi;
		else
		{
			//J矩阵
			double Jmatrix_left[2][2];
			//以上两个矩阵的积
			for (int i = 0; i<2; i++)
			for (int j = 0; j<2; j++)
			{
				Jmatrix_left[i][j] = 0;
				for (int k = 0; k<4; k++)
					Jmatrix_left[i][j] += gauss_dif[count1][i][k] * elenodes_left[k][j];
			}

			//求出J矩阵的行列式
			Jleft_val = Jmatrix_left[0][0] * Jmatrix_left[1][1] - Jmatrix_left[0][1] * Jmatrix_left[1][0];
		}

		//--------------------------------------------------
		//左端带权重形函数
		double gaus_nwct1[4];
		for (int i = 0; i<4; i++) gaus_nwct1[i] = gauss_nw[count1][i] * Jleft_val;

		//------------------------------------------------------------------------------------------------------------------------
		//循环外单元高斯点计算积分
		for (int count2 = 0; count2<gau_num; count2++)
		{
			//if(break_table[cogau++]) continue;  //First use and then ++, same to: { break_table[cogau] continue; cogau++; } 
			//--------------------------------------------
			//右端高斯点坐标
			Point_3D gaupoi_right(0, 0, 0);
			if (flag_right == 0)
			{
				gaupoi_right.x = gauss_po[count2][0] + elec_right[0];
				gaupoi_right.y = gauss_po[count2][1] + elec_right[1];
			}
			else
			{
				for (int i = 0; i<4; i++)
				{
					gaupoi_right.x += gauss_ns[count2][i] * elenodes_right[i][0];
					gaupoi_right.y += gauss_ns[count2][i] * elenodes_right[i][1];
				}
			}

			//--------------------------------------------
			//计算右端高斯点的权重
			double weight_right;
			if (elec_right[3] <= Zero) weight_right = 0.0;
			else if (fabs(elec_right[3] - 1.0) <= Zero) weight_right = 1.0;
			else weight_right = Gp_val[ass_ele][count2];

			//--------------------------------------------
			//计算权重函数值
			double alpha = 0.5*(weight_left + weight_right);

			if (alpha<1.0&&alpha>0.0) alpha_key[count1] = 0.5;
			else if (alpha_key[count1] == -1.0)  alpha_key[count1] = alpha;
			else if (alpha_key[count1] != alpha) alpha_key[count1] = 0.5;

			if (alpha <= Zero) continue;  //等于零值，没有长程效果

			//--------------------------------------------
			if (break_table[cogau++]) continue; //First use and then ++, same to: { break_table[cogau] continue; cogau++; } 

			//--------------------------------------------
			//计算长程作用衰减函数值
			const double x = gaupoi_right.x - gaupoi_left.x;
			const double y = gaupoi_right.y - gaupoi_left.y;

			const double dis_squr = x*x + y*y;
			const double poi_dis = sqrt(dis_squr);		//两点之间的距离

			if (poi_dis>peri_para.horizon_R || poi_dis<Zero) continue;		//圆形积分域

			//计算Ｊ矩阵；
			double Jright_val;
			if (flag_right == 0) Jright_val = Jacobi;
			else
			{
				//--------------------------------------------------
				//J矩阵
				double Jmatrix_right[2][2];
				//以上两个矩阵的积
				for (int i = 0; i<2; i++)
				for (int j = 0; j<2; j++)
				{
					Jmatrix_right[i][j] = 0;
					for (int k = 0; k<4; k++)
						Jmatrix_right[i][j] = Jmatrix_right[i][j] + gauss_dif[count2][i][k] * elenodes_right[k][j];
				}

				//求出J矩阵的行列式
				Jright_val = Jmatrix_right[0][0] * Jmatrix_right[1][1] - Jmatrix_right[0][1] * Jmatrix_right[1][0];
			}

			//--------------------------------------------
			//计算长程力矩阵
			double cos2pha;
			if (fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
			else cos2pha = x*x / (x*x + y*y);					//cos(pha)^2

			double sum = peri_para.acoe[0] + peri_para.acoe[1] * (2 * cos2pha - 1) + peri_para.acoe[2] * (8 * cos2pha*cos2pha - 8 * cos2pha + 1);

			double d_pri = 1.0 - damage_table[pri_ele][count1];
			double d_ass = 1.0 - damage_table[ass_ele][count2];
			const double gv = exp(-poi_dis / peri_para.intrinsic_L)*sum*alpha*Jright_val*weight[count2]*d_pri*d_pri/(d_ass*d_ass); //衰减以及权重函数值
			//const double gv = sum*alpha*Jright_val*weight[count2]/dis_squr; //constant micro-modulus (for stress)

			const double r_comp[6] = { x*x, y*y, 0, x*y, y * 0, 0 * x }; //等效长程力衰减计算
			const double temp_gmat[6] = { gv*r_comp[0], gv*r_comp[1], gv*r_comp[2], gv*r_comp[3], gv*r_comp[4], gv*r_comp[5] }; //衰减函数矩阵的对称项

			for (int i = 0; i<6; i++) Gmat[i] += temp_gmat[i];  //关于count2叠加

			for (int i = 0; i<6; i++)
			for (int j = 0; j<4; j++) GNmat[i][j] += temp_gmat[i] * gauss_ns[count2][j]; //关于count2并与形函数乘积叠加

			//--------------------------------------------
			//计算长程力等效刚度矩阵
			for (int i = 0; i<6; i++)
			for (int j = 0; j<6; j++)
				Long_range_stiffness[count1][i][j] += temp_gmat[i] * r_comp[j]; //相当于gv*r_comp[i]*r_comp[j]
		}

		//------------------------------------------------------------------------------------------------------------------------
		int ni = 0, nj = 0;
		double temxy = 0, nxy[6] = { 0 };
		for (int i = 0; i<4; i++)
		{
			nj = 0;
			for (int j = 0; j<4; j++)
			{
				temxy = gaus_nwct1[i] * gauss_ns[count1][j];
				for (int k = 0; k<6; k++)
					nxy[k] = temxy*Gmat[k];

				element_stiff_matrix1[ni][nj] += nxy[0];
				element_stiff_matrix1[ni + 1][nj + 1] += nxy[1];
				element_stiff_matrix1[ni + 2][nj + 2] += nxy[2];
				element_stiff_matrix1[ni][nj + 1] += nxy[3];
				element_stiff_matrix1[ni + 1][nj] += nxy[3];
				element_stiff_matrix1[ni + 1][nj + 2] += nxy[4];
				element_stiff_matrix1[ni + 2][nj + 1] += nxy[4];
				element_stiff_matrix1[ni][nj + 2] += nxy[5];
				element_stiff_matrix1[ni + 2][nj] += nxy[5];

				for (int k = 0; k<6; k++)
					nxy[k] = gaus_nwct1[i] * GNmat[k][j];

				element_stiff_matrix2[ni][nj] -= nxy[0];
				element_stiff_matrix2[ni + 1][nj + 1] -= nxy[1];
				element_stiff_matrix2[ni + 2][nj + 2] -= nxy[2];
				element_stiff_matrix2[ni][nj + 1] -= nxy[3];
				element_stiff_matrix2[ni + 1][nj] -= nxy[3];
				element_stiff_matrix2[ni + 1][nj + 2] -= nxy[4];
				element_stiff_matrix2[ni + 2][nj + 1] -= nxy[4];
				element_stiff_matrix2[ni][nj + 2] -= nxy[5];
				element_stiff_matrix2[ni + 2][nj] -= nxy[5];

				nj += 3;
			}
			ni += 3;
		}
	}

	//输出单刚阵，用于检查
	//hout << "element stiffness matrix1: " << endl;
	//for(int i=0; i<12; i++){
	//	for(int j=0; j<12;j++)
	//		hout << element_stiff_matrix1[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl;

	//hout << "element stiffness matrix2: " << endl;
	//for(int i=0; i<12; i++){
	//	for(int j=0; j<12;j++)
	//		hout << element_stiff_matrix2[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl;
}
//-----------------------------------------------------------------------------------------------
//Update damaged element matrices for 2D
int Global_Stiff_Matrix::Update_damaged_element_matrices_2D(const int &gaussnum, const struct Damage damages, const vector<MatPro> &mats, const vector<Node> &nodes, const vector<Element> &elements, const vector<double> &U_s,
																												   const vector<bool> &full_dam_eles, vector<MathMatrix> &ele_self_matrix, vector<vector<double> > &damage_table, int &damaged_sum, vector<bool> &dam_iter)const
{
	//---------------------------------------------------------------------------
	//Generating gaussian points
	Gauss gau;		//gaussian potins
	if(gau.Generate_gauss_2D(gaussnum)==0) return 0;

	//------------------------------------------------------------------------
	//Loops all elements
	cout << "-_- Updating the element matrices in damage part (local)..." << endl;
	hout << "-_- Updating the element matrices in damage part (local)..." << endl;

	//Implement openmp
	#pragma omp parallel
	{
		const int ES = (int)elements.size();
		double element_stiff_matrix[12][12];
		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<ES; i++) 
		{
			//Output calculating progress
			#pragma omp critical
			if(i%(ES/10)==0)	cout << i*100/ES <<"% done"<<endl;

			//---------------------------------------------------------------------------
			if(full_dam_eles[i]) continue;

			//---------------------------------------------------------------------------
			//Extracting coordinates of nodes in the element 
			double elenodes[4][2];
			for(int j=0; j<4; j++)
			{
				elenodes[j][0] = nodes[elements[i].nodes_id[j]].x;
				elenodes[j][1] = nodes[elements[i].nodes_id[j]].y;
			}
			//Extracting displacement solutions of nodes in the element
			double U_ele_nod[12] = { 0.0 };
			for(int j=0; j<4; j++)
			{
				U_ele_nod[3*j] = U_s[3*elements[i].nodes_id[j]];
				U_ele_nod[3*j+1] = U_s[3*elements[i].nodes_id[j]+1];
				U_ele_nod[3*j+2] = U_s[3*elements[i].nodes_id[j]+2];
			}

			//Define a record for updated broken bonds
			int damaged_num=0;
			//-------------------------------------------------------------------
			//updated stiffness matrix of an element
			if(elements[i].type==241&&(int)elements[i].nodes_id.size()==4)	Quadri_line_update_damaged(element_stiff_matrix, damage_table[i], damaged_num, damages, elenodes, U_ele_nod, mats, elements[i].mat, gau.gauss, gau.weight); //For a damaged brick element
			else 
			{
				cout << "Error: the type of element or the number of nodes in this element is wrong, when updating stiffness matrix of an element!"  << endl;
				hout << "Error: the type of element or the number of nodes in this element is wrong, when updating stiffness matrix of an element!"  << endl;

			}

			if(damaged_num>0)
			{
				//Updating data
				#pragma omp critical
				{
					damaged_sum += damaged_num;

					dam_iter[i]=true;

					for(int k=0; k<12; k++)
						for(int m=0; m<12; m++)
						{
							ele_self_matrix[i].element[k][m] -= element_stiff_matrix[k][m];
						}
				}
			}
		}
	}

	cout << "^_^ The element matrices updated successfully!" << endl;
	hout << "^_^ The element matrices updated successfully!" << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------
//Update the damaged brick element
void Global_Stiff_Matrix::Quadri_line_update_damaged(double (*element_stiff_matrix)[12], vector<double> &ele_dam_table, int &damaged_num, const struct Damage damages, const double (*elenodes)[2],
																								  const double *U_ele_nod, const vector<MatPro> &mats, const int &elemat, const vector<Node> &gauss, const vector<double> &wight)const
{
	double Ymax = 0.0;
	//--------------------------------------------------
	//Clearing
	for(int i=0; i<12; i++) 
		for (int j=0; j<12; j++)
			element_stiff_matrix[i][j] = 0.0;

	//--------------------------------------------------
	//循环高斯点计算积分
	for(int count=0; count<(int)wight.size(); count++)
	{
		if(damages.d_crit-ele_dam_table[count]<=Zero) continue; //Jump over full damaged gaussian points

		//-----------------------------------------------------------------
		//计算此单元所对应的材料弹性矩阵
		double ele_elas[6][6];		
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
				ele_elas[i][j] = mats[elemat+1].elas_matrix[i][j];
		
		//计算Ｊ矩阵
		//--------------------------------------------
		//形函数N对gauss[count].x, gauss[count].y的偏导矩阵
		double diff[2][4];
		diff[0][0]=-0.25*(1.0-gauss[count].y);
		diff[0][1]=-diff[0][0];
		diff[0][2]=0.25*(1.0+gauss[count].y);
		diff[0][3]=-diff[0][2];

		diff[1][0]=-0.25*(1.0-gauss[count].x);
		diff[1][1]=-0.25*(1.0+gauss[count].x);
		diff[1][2]=-diff[1][1];
		diff[1][3]=-diff[1][0];

		//--------------------------------------------------
		//J矩阵
		double Jmatrix[2][2];
		//以上两个矩阵的积
		for(int i=0; i<2; i++)
			for(int j=0; j<2; j++)
			{
				Jmatrix[i][j]=0;
				for(int k=0; k<4; k++)
				Jmatrix[i][j] += diff[i][k]*elenodes[k][j];
			}

		//--------------------------------------------------
		//求出J矩阵的行列式
		double J_val = Jmatrix[0][0]*Jmatrix[1][1]-Jmatrix[0][1]*Jmatrix[1][0];

		//----------------------------------------------------
		//求出J矩阵的逆矩阵
		double Jinverse[2][2];
		Jinverse[0][0]=Jmatrix[1][1]/J_val;
		Jinverse[1][1]=Jmatrix[0][0]/J_val;
		Jinverse[0][1]=-Jmatrix[0][1]/J_val;
		Jinverse[1][0]=-Jmatrix[1][0]/J_val;

		//-------------------------------------------------------
		//求出N对x,y的偏导
		double diffxy[2][4];
		for(int i=0; i<2; i++)
			for(int j=0; j<4; j++)
			{
				diffxy[i][j]=0;
				for(int k=0; k<2; k++)
					diffxy[i][j] += Jinverse[i][k]*diff[k][j];
			}

		//--------------------------------------------------------
		//求出B矩阵
		double B[6][12];
		for(int i=0; i<6; i++)
			for(int j=0; j<12; j++)
				B[i][j]=0;
		for(int i=0; i<4; i++)
		{
			B[0][i*3+0]=diffxy[0][i];
			B[1][i*3+1]=diffxy[1][i];
			B[3][i*3+0]=diffxy[1][i];
			B[3][i*3+1]=diffxy[0][i];
			B[4][i*3+2]=diffxy[1][i];
			B[5][i*3+2]=diffxy[0][i];
		}

		//--------------------------------------------------------
		//计算B_trans
		double B_trans[12][6];
		for(int i=0; i<12; i++)
			for(int j=0; j<6; j++)
				B_trans[i][j] = B[j][i];

		//--------------------------------------------------------------------------------------------
		//求出B_trans矩阵与ele_elas矩阵的乘积array1
		double array1[12][6];
		for(int i=0; i<12; i++)
			for(int j=0; j<6; j++)
			{
				array1[i][j]=0; 
				for(int k=0; k<6; k++)
					array1[i][j] += B_trans[i][k]*ele_elas[k][j];
			}

		//--------------------------------------------------------
		double Y_dam_force = 0.0;  //Damge force by energy decrease due to damage
		double tem[12][12];
		//求出array1矩阵与B矩阵的乘积array2
		double array2[12][12];
		for(int i=0; i<12; i++)
			for(int j=0; j<12; j++)
			{
				array2[i][j]=0;
				for(int k=0; k<6; k++)
					array2[i][j] += array1[i][k]*B[k][j];
				tem[i][j] = array2[i][j]*J_val*wight[count];
				Y_dam_force += 0.5*U_ele_nod[i]*array2[i][j]*U_ele_nod[j];
			}

			if(Ymax<Y_dam_force) Ymax=Y_dam_force;

			double K1=0.0;
			double K0=0.0;
			if(elemat==0) 
			{
				K0 = damages.k0;
				K1 = damages.k1;
			}
			else if(elemat==1)
			{
				K0 = damages.k0*0.8;
				K1 = damages.k1*0.8;
			}

		if((Y_dam_force - K1*ele_dam_table[count] - K0)>=Zero)
		{
			double new_dam = Y_dam_force/K1;
			if(damages.d_crit-new_dam<=Zero) new_dam = damages.d_crit;
			for(int i=0; i<12; i++)
				for(int j=0; j<12; j++)
					element_stiff_matrix[i][j] += tem[i][j]*(new_dam-ele_dam_table[count]);
			ele_dam_table[count] = new_dam;
			damaged_num ++;
		}
	}

//	cout << "Ymax= " << Ymax << " EleMat= " << elemat << endl;
//	hout << "Ymax= " << Ymax << " EleMat= " << elemat << endl;
	//输出单刚阵，用于检查
	//hout << "element stiffness matrix: " << endl;
	//for(int i=0; i<12; i++){
	//	for(int j=0; j<12; j++)
	//		hout << element_stiff_matrix[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl << endl;
}
//-----------------------------------------------------------------------------------------------
//Update damaged element matrices for 2D (New Algorithm)
int Global_Stiff_Matrix::Update_damaged_matrices_2D(const int &gaussnum, const struct Damage damages, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func, 
																					 const vector<MatPro> &mats, const vector<Node> &nodes, const vector<Element> &elements, const vector<double> &U_s, const vector<bool> &full_dam_eles,
																					 const vector<vector<double> > &Ak_val, const vector<vector<double> > &Ele_local_stiff, const vector<bool> *break_table, const vector<vector<double> > &Gp_val, 
																					 vector<MathMatrix> &ele_self_matrix, vector<vector<double> > *damage_table, vector<vector<double> > *Y_force, vector<double> &damage_en, int &damaged_sum, vector<bool> &dam_iter)const
{
	//---------------------------------------------------------------------------
	//Generating gaussian points
	Gauss gau;		//gaussian potins
	if(gau.Generate_gauss_2D(gaussnum)==0) return 0;

	//------------------------------------------------------------------------
	//Initialization of variables
	cout << "-_- Initializing variables..." << endl;
	hout << "-_- Initializing variables..." << endl;

	const int GS = (int)gau.gauss.size();
	const int ES = (int)elements.size();
	double (*gauss_ns)[4] = new double [GS][4];				//记录单元高斯点的形函数
	double (*gauss_nw)[4] = new double [GS][4];				//记录单元高斯点的形函数(带权重系数)
	double (*gauss_dif)[2][4] = new double [GS][2][4];		//记录单元高斯点形函数的导数
	double Jacobi = 0;																//标准正方形单元的雅可比值
	double (*gauss_po)[3] = new double [GS][3];				//标准正方形单元的高斯点坐标
	double (*ele_cent)[4] = new double [ES][4];					//记录单元中心点位置向量(x,y,z)分别放于[0],[1]和[2]中（以标准正方形单元中心点为原点）, [3]中放单元顶点权重的平均值

	if(com_mod=="Hybrid")
	{
		if(!(weight_func.num==1&&weight_func.shape[0]=="Null"))
		{
			//Calculate data of every gaussian point of element for 2D problem
			Generate_element_gauss_data_2D(elements, nodes, com_mod, weight_func, gau.gauss, gau.weight, gauss_ns, gauss_nw, gauss_dif, Jacobi, gauss_po, ele_cent);
		}
	}
	else if(com_mod!="Local") hout << "Attention: the computational mod is wrong out of Hybrid and Local." << endl;

	cout << "^_^ Variables initialized successfully!" << endl << endl;
	hout << "^_^ Variables initialized successfully!" << endl << endl;

	//------------------------------------------------------------------------
	//Loops all elements
	cout << "-_- Updating the element matrices in damage part (local)..." << endl;
	hout << "-_- Updating the element matrices in damage part (local)..." << endl;

	//------------------------------------------------------------------------
	double Y_eles_max = 0.0;
	//Implement openmp
	#pragma omp parallel
	{
		double Y_max = 0.0;
		const int ES = (int)elements.size();
		double element_stiff_matrix[12][12];
		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<ES; i++) 
		{
			//Output calculating progress
			#pragma omp critical
			if(i%(ES/10)==0)	cout << i*100/ES <<"% done"<<endl;

			//---------------------------------------------------------------------------
			//Fully damaged element
			if(full_dam_eles[i]) continue;

			//---------------------------------------------------------------------------
			//Extracting coordinates of nodes in the element 
			double elenodes[4][2];
			for(int j=0; j<4; j++)
			{
				elenodes[j][0] = nodes[elements[i].nodes_id[j]].x;
				elenodes[j][1] = nodes[elements[i].nodes_id[j]].y;
			}
			//Extracting displacement solutions of nodes in the element
			double U_ele_nod[12] = { 0.0 };
			for(int j=0; j<4; j++)
			{
				U_ele_nod[3*j] = U_s[3*elements[i].nodes_id[j]];
				U_ele_nod[3*j+1] = U_s[3*elements[i].nodes_id[j]+1];
				U_ele_nod[3*j+2] = U_s[3*elements[i].nodes_id[j]+2];
			}

			//------------------------------------------------------------------------		
			//Initialization of variable for nonlocal energy density at every Gaussian points
			vector<double> Y_nonlocal_energy(GS, 0.0);
			vector<double> Bond_dissip_en(GS, 0.0);

			//------------------------------------------------------------------------		
			if(com_mod=="Hybrid")
			{
				if(!(weight_func.num==1&&weight_func.shape[0]=="Null"))
				{
					//------------------------------------------------------------------------
					//Defining parameters of the center of the principal element
					double elec_left[4] = { ele_cent[i][0], ele_cent[i][1], ele_cent[i][2], ele_cent[i][3] };

					//------------------------------------------------------------------------
					//Accumulate the number of relative elements (less than i) //attention: Openmp Parallel
					long int ref_tot = 0;
					for(int j=0; j<i; j++) ref_tot += elements[j].relative_eles.size();

					//------------------------------------------------------------------------
					//Loops element neighbours
					for(int j=0; j<(int)elements[i].relative_eles.size(); j++)
					{
						ref_tot++; //Accumlation (注意不能移动到下面 因为以下有continue语句)
						//--------------------------------------------------
						//Defining parameters of the center of the associated element 
						const int ere = elements[i].relative_eles[j]; 
						double elec_right[4] = { ele_cent[ere][0], ele_cent[ere][1], ele_cent[ere][2], ele_cent[ere][3] };

						//--------------------------------------------------
						//A preliminary examination for skipping both elements completely in the local continuum area
						//（注意elec_left[3]和elec_right[3]中分别是两个四边形顶点权重值得平均）
						if(elec_left[3]<=Zero&&elec_right[3]<=Zero)  continue;

						//--------------------------------------------------
						//Extracting coordinates of nodes in the associated element
						double elenodes_relative[4][2];
						for(int k=0; k<4; k++)
						{
							elenodes_relative[k][0] = nodes[elements[ere].nodes_id[k]].x;
							elenodes_relative[k][1] = nodes[elements[ere].nodes_id[k]].y;
						}
						//Extracting displacement solutions of nodes in the element
						double U_elenod_relat[12] = { 0.0 };
						for(int k=0; k<4; k++)
						{
							U_elenod_relat[3*k] = U_s[3*elements[ere].nodes_id[k]];
							U_elenod_relat[3*k+1] = U_s[3*elements[ere].nodes_id[k]+1];
							U_elenod_relat[3*k+2] = U_s[3*elements[ere].nodes_id[k]+2];
						}

						//--------------------------------------------------
						//Calculate the position of break_table
						const long int pos_bt = (ref_tot-1)*GS*GS;

						//------------------------------------------------------------------------
						//Generate the energy density at every gaussian point based on non-local continuum model judging broken bonds for 2D problem
						Longforce_Energy_Density_Brokeless_2D(Y_nonlocal_energy, Bond_dissip_en, U_ele_nod, U_elenod_relat, elements[i].flag, elements[ere].flag, elenodes,
										elenodes_relative, gauss_ns, gauss_dif, gau.weight, peri_para, Jacobi, gauss_po, elec_left, elec_right, break_table, pos_bt, Gp_val, i, ere);
					}
				}
			}
			else if(com_mod!="Local") hout << "Attention: the computational mod is wrong out of Hybrid and Local." << endl;

			//Define a record for updated damaged elements
			int damaged_num=0;
			int increased_num = 0;
			double temp_Ymax = 0.0;
			double Damage_ener = 0.0;
			//-------------------------------------------------------------------
			//updated stiffness matrix of an element
			if(elements[i].type==241&&(int)elements[i].nodes_id.size()==4)
			{
//bool QDmark = false;
//for(int j=0; j<4; j++)
//{
//	if(nodes[elements[i].nodes_id[j]].y < 78) QDmark = true;
//	if(QDmark) break;
//}
				//For a damaged brick element
//if(QDmark)				
				Quadri_damage_update(element_stiff_matrix, Ele_local_stiff[i], Y_nonlocal_energy, Bond_dissip_en, damage_table[0][i], Y_force[0][i], damage_table[1][i], Y_force[1][i],
					                                 temp_Ymax, Damage_ener, increased_num, damaged_num, damages, elenodes, U_ele_nod, mats, elements[i].mat, gau.gauss, gau.weight);
			}
			else
			{
				cout << "Error: the type of element or the number of nodes in this element is wrong, when updating stiffness matrix of an element!"  << endl;
				hout << "Error: the type of element or the number of nodes in this element is wrong, when updating stiffness matrix of an element!"  << endl;
			}

			//Updating data
			if(Y_max<temp_Ymax) Y_max = temp_Ymax;
			if(damaged_num>0)
			{
				#pragma omp critical
				{
					if(increased_num>0)	damaged_sum += increased_num;
					damage_en[i] += Damage_ener;
					dam_iter[i]=true;

					for(int k=0; k<12; k++)
						for(int m=0; m<12; m++)
						{
							ele_self_matrix[i].element[k][m] -= element_stiff_matrix[k][m];
						}
				 }
			 }
		}

		#pragma omp critical
		if(Y_eles_max<Y_max) Y_eles_max = Y_max;
	}
	cout << "Y_eles_max = " << Y_eles_max << endl;
	hout << "Y_eles_max = " << Y_eles_max << endl;

	cout << "^_^ The element matrices updated successfully!" << endl;
	hout << "^_^ The element matrices updated successfully!" << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------
//Update the damaged brick element (New Algorithm)
void Global_Stiff_Matrix::Quadri_damage_update(double (*element_stiff_matrix)[12], const vector<double> &Ele_stiff, const vector<double> Y_nonlocal_energy, vector<double> &Bond_dissip_en, 
												const vector<double> &ele_damtab0, const vector<double> &Yfors0, vector<double> &ele_damtab1, vector<double> &Yfors1, double &Ymax, double &Damage_ener,
												int &increased_num, int &damaged_num, const struct Damage damages, const double (*elenodes)[2], const double *U_ele_nod, const vector<MatPro> &mats, 
												const int &elemat, const vector<Node> &gauss, const vector<double> &wight)const
{
	//--------------------------------------------------
	//Clearing
	for(int i=0; i<12; i++) 
		for (int j=0; j<12; j++)
			element_stiff_matrix[i][j] = 0.0;

	//--------------------------------------------------
	//循环高斯点计算积分
	for(int count=0; count<(int)wight.size(); count++)
	{
		if(damages.d_crit-ele_damtab0[count]<=Zero) continue; //Jump over full damaged gaussian points

		//-----------------------------------------------------------------
		//计算此单元所对应的材料弹性矩阵
		double ele_elas[6][6];
		int numb = 36*count;
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
				ele_elas[i][j] = Ele_stiff[numb++];  //注意无论Local和Hybrid这里都成立
		
		//计算Ｊ矩阵
		//--------------------------------------------
		//形函数N对gauss[count].x, gauss[count].y的偏导矩阵
		double diff[2][4];
		diff[0][0]=-0.25*(1.0-gauss[count].y);
		diff[0][1]=-diff[0][0];
		diff[0][2]=0.25*(1.0+gauss[count].y);
		diff[0][3]=-diff[0][2];

		diff[1][0]=-0.25*(1.0-gauss[count].x);
		diff[1][1]=-0.25*(1.0+gauss[count].x);
		diff[1][2]=-diff[1][1];
		diff[1][3]=-diff[1][0];

		//--------------------------------------------------
		//J矩阵
		double Jmatrix[2][2];
		//以上两个矩阵的积
		for(int i=0; i<2; i++)
			for(int j=0; j<2; j++)
			{
				Jmatrix[i][j]=0;
				for(int k=0; k<4; k++)
				Jmatrix[i][j] += diff[i][k]*elenodes[k][j];
			}

		//--------------------------------------------------
		//求出J矩阵的行列式
		double J_val = Jmatrix[0][0]*Jmatrix[1][1]-Jmatrix[0][1]*Jmatrix[1][0];

		//----------------------------------------------------
		//求出J矩阵的逆矩阵
		double Jinverse[2][2];
		Jinverse[0][0]=Jmatrix[1][1]/J_val;
		Jinverse[1][1]=Jmatrix[0][0]/J_val;
		Jinverse[0][1]=-Jmatrix[0][1]/J_val;
		Jinverse[1][0]=-Jmatrix[1][0]/J_val;

		//-------------------------------------------------------
		//求出N对x,y的偏导
		double diffxy[2][4];
		for(int i=0; i<2; i++)
			for(int j=0; j<4; j++)
			{
				diffxy[i][j]=0;
				for(int k=0; k<2; k++)
					diffxy[i][j] += Jinverse[i][k]*diff[k][j];
			}

		//--------------------------------------------------------
		//求出B矩阵
		double B[6][12];
		for(int i=0; i<6; i++)
			for(int j=0; j<12; j++)
				B[i][j]=0;
		for(int i=0; i<4; i++)
		{
			B[0][i*3+0]=diffxy[0][i];
			B[1][i*3+1]=diffxy[1][i];
			B[3][i*3+0]=diffxy[1][i];
			B[3][i*3+1]=diffxy[0][i];
			B[4][i*3+2]=diffxy[1][i];
			B[5][i*3+2]=diffxy[0][i];
		}

		//--------------------------------------------------------
		//计算B_trans
		double B_trans[12][6];
		for(int i=0; i<12; i++)
			for(int j=0; j<6; j++)
				B_trans[i][j] = B[j][i];

		//--------------------------------------------------------------------------------------------
		//求出B_trans矩阵与ele_elas矩阵的乘积array1
		double array1[12][6];
		for(int i=0; i<12; i++)
			for(int j=0; j<6; j++)
			{
				array1[i][j]=0; 
				for(int k=0; k<6; k++)
					array1[i][j] += B_trans[i][k]*ele_elas[k][j];
			}

		//--------------------------------------------------------
		double Y_local_energy = 0.0;  //Damge force by energy decrease due to damage
		//求出array1矩阵与B矩阵的乘积array2
		double array2[12][12];
		for(int i=0; i<12; i++)
			for(int j=0; j<12; j++)
			{
				array2[i][j]=0.0;
				for(int k=0; k<6; k++)	array2[i][j] += array1[i][k]*B[k][j];
				Y_local_energy += 0.5*U_ele_nod[i]*array2[i][j]*U_ele_nod[j];
			}

		//--------------------------------------------------------
		double Y_dam_force = Y_local_energy + Y_nonlocal_energy[count];

		if(Ymax<Y_dam_force) Ymax=Y_dam_force;

		if ((Y_dam_force - damages.k1*ele_damtab0[count] - damages.k0) >= Zero)
		{
			double delta_dam = (Y_dam_force - Yfors0[count]) / damages.k1;
			double temp_ele_dam = ele_damtab0[count] + delta_dam; 

			//----------------------------------------------------------------------------------------------------------------------------------------------------------
			//Explicit Algorithm
			if (temp_ele_dam>ele_damtab1[count]&&ele_damtab1[count]<damages.d_crit)
			{
				increased_num++;

				if (temp_ele_dam <= damages.d_crit)
				{
					ele_damtab1[count] = temp_ele_dam;
					Yfors1[count] = Y_dam_force;
				}
				else
				{
					ele_damtab1[count] = damages.d_crit;
					Yfors1[count] = (Y_dam_force - Yfors0[count])*(damages.d_crit - ele_damtab0[count]) / delta_dam + Yfors0[count];
				}
			}

			//----------------------------------------------------------------------------------------------------------------------------------------------------------
			//Implicit Algorithm
			//损伤值有变动，变动比上一次可大可小，但是都需要修改，
			//除了已经达到完全损伤了	, 再次计算后还是完全损伤的情况
			//注意这里有个计算精度的问题（精度越高，收敛越慢），这里的取值是100*Zero
			//if (fabs(temp_ele_dam - ele_damtab1[count])>1.0e-6&&
			//	(!(temp_ele_dam >= damages.d_crit && ele_damtab1[count] == damages.d_crit)))
			//{
			//	increased_num++;
			//	if (temp_ele_dam <= damages.d_crit)
			//	{
			//		ele_damtab1[count] = temp_ele_dam;
			//		Yfors1[count] = Y_dam_force;
			//	}
			//	else
			//	{
			//		ele_damtab1[count] = damages.d_crit;
			//		Yfors1[count] = (Y_dam_force - Yfors0[count])*(damages.d_crit - ele_damtab0[count]) / delta_dam + Yfors0[count];
			//	}
			//}
		}

		//注意这里是根据所有ele_damtab1（无论是否是这次更新的ele_damtab1）
		//来计算刚度的减少量element_stiff_matrix（这里是加号，因为外部的累计是减号）
		double diff_damtab = ele_damtab1[count] - ele_damtab0[count];
		if (diff_damtab>Zero)
		{
			for (int i=0; i<12; i++)
				for (int j=0; j<12; j++)
					element_stiff_matrix[i][j] += array2[i][j]*J_val*wight[count]*diff_damtab;

			Damage_ener += Yfors1[count]*diff_damtab*J_val*wight[count];

			damaged_num++;
		}
	}

//	cout << "Ymax= " << Ymax << " EleMat= " << elemat << endl;
//	hout << "Ymax= " << Ymax << " EleMat= " << elemat << endl;
	//输出单刚阵，用于检查
	//hout << "element stiffness matrix: " << endl;
	//for(int i=0; i<12; i++){
	//	for(int j=0; j<12; j++)
	//		hout << element_stiff_matrix[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl << endl;
}
//-----------------------------------------------------------------------------------------------
//Update the damaged brick element (Last Algorithm)
void Global_Stiff_Matrix::Quadri_damage_update(double (*element_stiff_matrix)[12], const vector<double> &Ele_stiff, const vector<double> Y_nonlocal_energy, const vector<double> &ele_damtab0, 
												const vector<double> &Yfors0, vector<double> &ele_damtab1, vector<double> &Yfors1, vector<double> &ele_damtab2, vector<double> &Yfors2, double &Ymax, 
												double &Damage_ener, int &increased_num, int &damaged_num, const struct Damage damages, const double (*elenodes)[2], const double *U_ele_nod, 
												const vector<MatPro> &mats, const int &elemat, const vector<Node> &gauss, const vector<double> &wight)const
{
	//--------------------------------------------------
	//Clearing
	for(int i=0; i<12; i++) 
		for (int j=0; j<12; j++)
			element_stiff_matrix[i][j] = 0.0;

	//--------------------------------------------------
	//循环高斯点计算积分
	for(int count=0; count<(int)wight.size(); count++)
	{
		if(damages.d_crit-ele_damtab0[count]<=Zero) continue; //Jump over full damaged gaussian points

		//-----------------------------------------------------------------
		//计算此单元所对应的材料弹性矩阵
		double ele_elas[6][6];
		int numb = 36*count;
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
				ele_elas[i][j] = Ele_stiff[numb++];  //注意无论Local和Hybrid这里都成立
		
		//计算Ｊ矩阵
		//--------------------------------------------
		//形函数N对gauss[count].x, gauss[count].y的偏导矩阵
		double diff[2][4];
		diff[0][0]=-0.25*(1.0-gauss[count].y);
		diff[0][1]=-diff[0][0];
		diff[0][2]=0.25*(1.0+gauss[count].y);
		diff[0][3]=-diff[0][2];

		diff[1][0]=-0.25*(1.0-gauss[count].x);
		diff[1][1]=-0.25*(1.0+gauss[count].x);
		diff[1][2]=-diff[1][1];
		diff[1][3]=-diff[1][0];

		//--------------------------------------------------
		//J矩阵
		double Jmatrix[2][2];
		//以上两个矩阵的积
		for(int i=0; i<2; i++)
			for(int j=0; j<2; j++)
			{
				Jmatrix[i][j]=0;
				for(int k=0; k<4; k++)
				Jmatrix[i][j] += diff[i][k]*elenodes[k][j];
			}

		//--------------------------------------------------
		//求出J矩阵的行列式
		double J_val = Jmatrix[0][0]*Jmatrix[1][1]-Jmatrix[0][1]*Jmatrix[1][0];

		//----------------------------------------------------
		//求出J矩阵的逆矩阵
		double Jinverse[2][2];
		Jinverse[0][0]=Jmatrix[1][1]/J_val;
		Jinverse[1][1]=Jmatrix[0][0]/J_val;
		Jinverse[0][1]=-Jmatrix[0][1]/J_val;
		Jinverse[1][0]=-Jmatrix[1][0]/J_val;

		//-------------------------------------------------------
		//求出N对x,y的偏导
		double diffxy[2][4];
		for(int i=0; i<2; i++)
			for(int j=0; j<4; j++)
			{
				diffxy[i][j]=0;
				for(int k=0; k<2; k++)
					diffxy[i][j] += Jinverse[i][k]*diff[k][j];
			}

		//--------------------------------------------------------
		//求出B矩阵
		double B[6][12];
		for(int i=0; i<6; i++)
			for(int j=0; j<12; j++)
				B[i][j]=0;
		for(int i=0; i<4; i++)
		{
			B[0][i*3+0]=diffxy[0][i];
			B[1][i*3+1]=diffxy[1][i];
			B[3][i*3+0]=diffxy[1][i];
			B[3][i*3+1]=diffxy[0][i];
			B[4][i*3+2]=diffxy[1][i];
			B[5][i*3+2]=diffxy[0][i];
		}

		//--------------------------------------------------------
		//计算B_trans
		double B_trans[12][6];
		for(int i=0; i<12; i++)
			for(int j=0; j<6; j++)
				B_trans[i][j] = B[j][i];

		//--------------------------------------------------------------------------------------------
		//求出B_trans矩阵与ele_elas矩阵的乘积array1
		double array1[12][6];
		for(int i=0; i<12; i++)
			for(int j=0; j<6; j++)
			{
				array1[i][j]=0; 
				for(int k=0; k<6; k++)
					array1[i][j] += B_trans[i][k]*ele_elas[k][j];
			}

		//--------------------------------------------------------
		double Y_local_energy = 0.0;  //Damge force by energy decrease due to damage
		//求出array1矩阵与B矩阵的乘积array2
		double array2[12][12];
		for(int i=0; i<12; i++)
			for(int j=0; j<12; j++)
			{
				array2[i][j]=0.0;
				for(int k=0; k<6; k++)	array2[i][j] += array1[i][k]*B[k][j];
				Y_local_energy += 0.5*U_ele_nod[i]*array2[i][j]*U_ele_nod[j];
			}

		//--------------------------------------------------------
		double Y_dam_force = Y_local_energy + Y_nonlocal_energy[count];

		if(Ymax<Y_dam_force) Ymax=Y_dam_force;

		if ((Y_dam_force - damages.k1*ele_damtab2[count] - damages.k0) >= Zero)
		{
			double delta_dam = (Y_dam_force - Yfors2[count]) / damages.k1;
			double temp_ele_dam = ele_damtab2[count] + delta_dam; 

			//----------------------------------------------------------------------------------------------------------------------------------------------------------
			//Explicit Algorithm
			if (temp_ele_dam>ele_damtab1[count]&&ele_damtab1[count]<damages.d_crit)
			{
				increased_num++;

				if (temp_ele_dam <= damages.d_crit)
				{
					ele_damtab1[count] = temp_ele_dam;
					Yfors1[count] = Y_dam_force;
				}
				else
				{
					ele_damtab1[count] = damages.d_crit;
					Yfors1[count] = (Y_dam_force - Yfors0[count])*(damages.d_crit - ele_damtab0[count]) / delta_dam + Yfors0[count];
				}
			}

			//----------------------------------------------------------------------------------------------------------------------------------------------------------
			//Implicit Algorithm
			//损伤值有变动，变动比上一次可大可小，但是都需要修改，
			//除了已经达到完全损伤了	, 再次计算后还是完全损伤的情况
			//注意这里有个计算精度的问题（精度越高，收敛越慢），这里的取值是100*Zero
			//if (fabs(temp_ele_dam - ele_damtab1[count])>1.0e-6&&
			//	(!(temp_ele_dam >= damages.d_crit && ele_damtab1[count] == damages.d_crit)))
			//{
			//	increased_num++;
			//	if (temp_ele_dam <= damages.d_crit)
			//	{
			//		ele_damtab1[count] = temp_ele_dam;
			//		Yfors1[count] = Y_dam_force;
			//	}
			//	else
			//	{
			//		ele_damtab1[count] = damages.d_crit;
			//		Yfors1[count] = (Y_dam_force - Yfors0[count])*(damages.d_crit - ele_damtab0[count]) / delta_dam + Yfors0[count];
			//	}
			//}
		}

		//注意这里是根据所有ele_damtab1（无论是否是这次更新的ele_damtab1）
		//来计算刚度的减少量element_stiff_matrix（这里是加号，因为外部的累计是减号）
		double diff_damtab = ele_damtab1[count] - ele_damtab0[count];
		if (diff_damtab>Zero)
		{
			for (int i=0; i<12; i++)
				for (int j=0; j<12; j++)
					element_stiff_matrix[i][j] += array2[i][j]*J_val*wight[count]*diff_damtab;

			Damage_ener += Yfors1[count]*diff_damtab*J_val*wight[count];

			damaged_num++;
		}
	}

//	cout << "Ymax= " << Ymax << " EleMat= " << elemat << endl;
//	hout << "Ymax= " << Ymax << " EleMat= " << elemat << endl;
	//输出单刚阵，用于检查
	//hout << "element stiffness matrix: " << endl;
	//for(int i=0; i<12; i++){
	//	for(int j=0; j<12; j++)
	//		hout << element_stiff_matrix[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl << endl;
}
//-----------------------------------------------------------------------------------------------
//Update element matrices with broken bonds (after up to damage criterion) for 2D problem
int Global_Stiff_Matrix::Update_nonlocal_element_matrices_2D(const int &gaussnum, const struct Damage damages, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func, const vector<MatPro> &mats,
																												   const vector<Node> &nodes, const vector<vector<double> > &Gp_val, const vector<double> &U_s, vector<Element> &elements, const vector<vector<double> > &damage_table,
																												   vector<bool> &full_dam_eles,  vector<MathMatrix> &ele_self_matrix, vector<vector<MathMatrix> > &ele_relative_matrix, vector<bool> &break_table, int &broken_sum, vector<bool> &break_iter)const
{
	//---------------------------------------------------------------------------
	//Generating gaussian points
	Gauss gau;		//gaussian potins
	if(gau.Generate_gauss_2D(gaussnum)==0) return 0;

	//------------------------------------------------------------------------
	//Initialization of variables
	cout << "-_- Initializing variables..." << endl;
	hout << "-_- Initializing variables..." << endl;

	const int GS = (int)gau.gauss.size();
	const int ES = (int)elements.size();
	double (*gauss_ns)[4] = new double [GS][4];				//记录单元高斯点的形函数
	double (*gauss_nw)[4] = new double [GS][4];				//记录单元高斯点的形函数(带权重系数)
	double (*gauss_dif)[2][4] = new double [GS][2][4];		//记录单元高斯点形函数的导数
	double Jacobi = 0;																//标准正方形单元的雅可比值
	double (*gauss_po)[3] = new double [GS][3];				//标准正方形单元的高斯点坐标
	double (*ele_cent)[4] = new double [ES][4];					//记录单元中心点位置向量(x,y,z)分别放于[0],[1]和[2]中（以标准正方形单元中心点为原点）, [3]中放单元顶点权重的平均值
	
	//Calculate data of every gaussian point of element for 2D problem
	Generate_element_gauss_data_2D(elements, nodes, com_mod, weight_func, gau.gauss, gau.weight, gauss_ns, gauss_nw, gauss_dif, Jacobi, gauss_po, ele_cent);
	
	cout << "^_^ Variables initialized successfully!" << endl << endl;
	hout << "^_^ Variables initialized successfully!" << endl << endl;

	//------------------------------------------------------------------------
	//Loops all elements
	cout << "-_- Updating the element matrices in nonlocal part..." << endl;
	hout << "-_- Updating the element matrices in nonlocal part..." << endl;
	//执行openmp
	#pragma omp parallel
	{
		double element_stiff_matrix1[12][12];
		double element_stiff_matrix2[12][12];
		double *alpha_key = new double [GS];
		
		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<ES; i++) 
		{
			//Output calculating progress
			#pragma omp critical
			if(i%(ES/10)==0)	cout << i*100/ES <<"% done"<<endl;
			
			//---------------------------------------------------------------------------
			//Extracting coordinates of nodes in the principal element 
			double elenodes[4][2];
			for(int j=0; j<4; j++)
			{
				elenodes[j][0] = nodes[elements[i].nodes_id[j]].x;
				elenodes[j][1] = nodes[elements[i].nodes_id[j]].y;
			}
			//Extracting coordinates of nodes in the updated principal element with displacement solutions
			double new_elenodes[4][2];
			for(int j=0; j<4; j++)
			{
				new_elenodes[j][0] = elenodes[j][0] + U_s[3*elements[i].nodes_id[j]];
				new_elenodes[j][1] = elenodes[j][1] + U_s[3*elements[i].nodes_id[j]+1];
			}

			//---------------------------------------------------------------------------
			//Calculating hybrid models
			if(com_mod=="Hybrid")
			{	
				//------------------------------------------------------------------------
				//Initialazing alpha weighting coefficients
				for(int j=0; j<GS; j++) alpha_key[j] = -1.0;

				//Defining parameters of the center of the principal element
				double elec_left[4] = { ele_cent[i][0], ele_cent[i][1], ele_cent[i][2], ele_cent[i][3] };

				//------------------------------------------------------------------------
				//Accumulate the number of relative elements (less than i) //attention: Openmp Parallel
				long int ref_tot = 0;
				for(int j=0; j<i; j++) ref_tot += elements[j].relative_eles.size();

				//------------------------------------------------------------------------
				//Loops element neighbours
				for(int j=0; j<(int)elements[i].relative_eles.size(); j++)
				{
					ref_tot++; //Accumlation
					//--------------------------------------------------
					//Defining parameters of the center of the associated element 
					const int ere = elements[i].relative_eles[j]; 
					double elec_right[4] = { ele_cent[ere][0], ele_cent[ere][1], ele_cent[ere][2], ele_cent[ere][3] };

if((!full_dam_eles[i])&&(!full_dam_eles[ere])) continue;

					//--------------------------------------------------
					//A preliminary examination for skipping both elements completely in the local continuum area
					if(elec_left[3]<=Zero&&elec_right[3]<=Zero) 
					{
						for(int count=0; count<GS; count++)
						{
							if(alpha_key[count]==-1.0)  alpha_key[count] = 0.0;
							else if(alpha_key[count]!=0.0) alpha_key[count] = 0.5;
						}
						continue;
					}

					//--------------------------------------------------
					//Extracting coordinates of nodes in the associated element
					double elenodes_relative[4][2];
					for(int k=0; k<4; k++)
					{
						elenodes_relative[k][0] = nodes[elements[ere].nodes_id[k]].x;
						elenodes_relative[k][1] = nodes[elements[ere].nodes_id[k]].y;
					}
					//Extracting coordinates of nodes in the updated associated element with displacement solutions
					//extracting new relative element coordinates
					double new_elenodes_relative[4][2];
					for(int k=0; k<4; k++)
					{
						new_elenodes_relative[k][0] = elenodes_relative[k][0] + U_s[3*elements[ere].nodes_id[k]];
						new_elenodes_relative[k][1] = elenodes_relative[k][1] + U_s[3*elements[ere].nodes_id[k]+1];
					}

					//--------------------------------------------------
					//Calculate the position of break_table
					const long int pos_bt = (ref_tot-1)*GS*GS;
					//Define a record for updated broken bonds
					int broken_num=0;
					//------------------------------------------------------------------------
					//Generate the element matrices based on non-local continuum model (long-range forces) judging broken bonds for 2D problem
					Generate_Longforce_Elestiff_Breaks_2D(element_stiff_matrix1, element_stiff_matrix2, alpha_key, break_table, broken_num, damages, damage_table, elements[i].Dissipative_energy,
																					   weight_func, elements[i].flag, elements[ere].flag, elenodes, new_elenodes, elenodes_relative, new_elenodes_relative, gauss_ns, gauss_nw,
																					   gauss_dif, gau.weight, peri_para, Jacobi, gauss_po, elec_left, elec_right, pos_bt, Gp_val, i, ere);
					
					if(broken_num>0)
					{
						//Updating data
						#pragma omp critical	
						{
							broken_sum += broken_num;

							break_iter[i]=true;
							break_iter[ere]=true;

							for(int k=0; k<12; k++)
								for(int m=0; m<12; m++)
								{
									ele_self_matrix[i].element[k][m] -= element_stiff_matrix1[k][m];
									ele_relative_matrix[i][j].element[k][m] -= element_stiff_matrix2[k][m];
								}
						}
					}
				}
			}		
		}
		//---------------------------------------------------------------------------
		//Delete pointer
		delete[] alpha_key;
	}

	//---------------------------------------------------------------------------
	//Delete pointers
	delete[] gauss_ns;
	delete[] gauss_nw;
	delete[] gauss_dif;
	delete[] gauss_po;
	delete[] ele_cent;

	cout << "^_^ The element matrices updated successfully!" << endl;
	hout << "^_^ The element matrices updated successfully!" << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------
//Generate the element matrix based on non-local continuum model (long-range forces) with broken bonds (after up to damage criterion) for 2D problem
void Global_Stiff_Matrix::Generate_Longforce_Elestiff_Breaks_2D(double (*element_stiff_matrix1)[12], double (*element_stiff_matrix2)[12], double *alpha_key, vector<bool> &break_table, int &broken_num, const struct Damage damages,  const vector<vector<double> > &damage_table,
																												double &Dissip_ener, const struct Weight_func &weight_func, const int &flag_left, const int &flag_right, const double (*elenodes_left)[2], const double(*n_ele_l)[2], const double (*elenodes_right)[2],
																												const double(*n_ele_r)[2], const double (*gauss_ns)[4], const double (*gauss_nw)[4], const double (*gauss_dif)[2][4], const vector<double> &weight, const struct Peri_para &peri_para,
																												const double &Jacobi, const double (*gauss_po)[3], const double elec_left[], const double elec_right[], const long int &pos_bt, const vector<vector<double> > &Gp_val, const int &pri_ele,
																												const int &ass_ele)const
{
	long int cogau = pos_bt;
	//--------------------------------------------------
	//Clearing
	for(int i=0; i<12; i++) 
		for (int j=0; j<12; j++)
		{
			element_stiff_matrix1[i][j] = 0;
			element_stiff_matrix2[i][j] = 0;
		}

	//--------------------------------------------------	
	//Loop the gaussian points in the current element
	const int gau_num = (int)weight.size();
	for(int count1=0; count1<gau_num; count1++)
	{
		//--------------------------------------------------
		//The bond can break until up to damage criterion
		if(damages.d_crit-damage_table[pri_ele][count1]>Zero) continue;

		//--------------------------------------------------
		int bk_count = 0;
		double Gmat[6] = {0};
		double GNmat[6][4] = { {0}, {0}, {0}, {0}, {0}, {0} };
		//--------------------------------------------------	
		//Coordinates of gaussian point (left)
		Point_3D gaupoi_left(0, 0, 0);
		if(flag_left==0)
		{
			gaupoi_left.x = gauss_po[count1][0] + elec_left[0];
			gaupoi_left.y = gauss_po[count1][1] + elec_left[1];
		}
		else
		{
			for(int i=0; i<4; i++) 
			{
				gaupoi_left.x += gauss_ns[count1][i]*elenodes_left[i][0];
				gaupoi_left.y += gauss_ns[count1][i]*elenodes_left[i][1];
			}
		}
		//--------------------------------------------------	
		//New coordinates of gaussian point (left)
		Point_3D n_gaupoi_left(0, 0, 0);
		for(int i=0; i<4; i++) 
		{
			n_gaupoi_left.x += gauss_ns[count1][i]*n_ele_l[i][0];
			n_gaupoi_left.y += gauss_ns[count1][i]*n_ele_l[i][1];
		}
		//--------------------------------------------------
		//Evaluate values of weighting function at gaussian point (left)
		double weight_left;
		if(elec_left[3]<=Zero) weight_left = 0.0;
		else if(fabs(elec_left[3]-1.0)<=Zero) weight_left = 1.0;
		else weight_left = Gp_val[pri_ele][count1];

		//--------------------------------------------------
		//Evaluate Jacobi matix
		double Jleft_val;
		if(flag_left==0) Jleft_val = Jacobi;
		else
		{
			//Jacobi matix
			double Jmatrix_left[2][2];
			//以上两个矩阵的积
			for(int i=0; i<2; i++)
				for(int j=0; j<2; j++)
				{
					Jmatrix_left[i][j]=0;
					for(int k=0; k<4; k++)
						Jmatrix_left[i][j] += gauss_dif[count1][i][k]*elenodes_left[k][j];
				}

			//求出J矩阵的行列式
			Jleft_val = Jmatrix_left[0][0]*Jmatrix_left[1][1]-Jmatrix_left[0][1]*Jmatrix_left[1][0];	
		}
		
		//--------------------------------------------------
		//Evaluate shape-function values with weighting
		double gaus_nwct1[4];
		for(int i=0; i<4; i++) gaus_nwct1[i] = gauss_nw[count1][i]*Jleft_val;

		//------------------------------------------------------------------------------------------------------------------------
		//iteration on second ghost point
		for(int count2=0; count2<gau_num; count2++)
		{
			//--------------------------------------------------
			//The bond can break until up to damage criterion
			if(damages.d_crit-damage_table[ass_ele][count2]>Zero) continue;			

			//--------------------------------------------------
			if(break_table[cogau++]) continue;  //First use and then ++, same to: { break_table[cogau] continue; cogau++; } 
			//--------------------------------------------
			//Coordinates of gaussian point (right)
			Point_3D gaupoi_right(0, 0, 0);
			if(flag_right==0)
			{
				gaupoi_right.x = gauss_po[count2][0] + elec_right[0];
				gaupoi_right.y = gauss_po[count2][1] + elec_right[1];
			}
			else
			{
				for(int i=0; i<4; i++) 
				{
					gaupoi_right.x +=  gauss_ns[count2][i]*elenodes_right[i][0];
					gaupoi_right.y +=  gauss_ns[count2][i]*elenodes_right[i][1];
				}
			}
			//--------------------------------------------
			//New coordinates of gaussian point (right)
			Point_3D n_gaupoi_right(0, 0, 0);
			for(int i=0; i<4; i++) 
			{
				n_gaupoi_right.x +=  gauss_ns[count2][i]*n_ele_r[i][0];
				n_gaupoi_right.y +=  gauss_ns[count2][i]*n_ele_r[i][1];
			}

			//--------------------------------------------
			//Evaluate values of weighting function at gaussian point (right)
			double weight_right;
			if(elec_right[3]<=Zero) weight_right = 0.0;
			else if(fabs(elec_right[3]-1.0)<=Zero) weight_right = 1.0;
			else weight_right = Gp_val[ass_ele][count2];

			//--------------------------------------------
			//计算权重函数值
			double alpha = 0.5*(weight_left+weight_right);

			if(alpha<1.0&&alpha>0.0) alpha_key[count1] = 0.5;
			else if(alpha_key[count1]==-1.0)  alpha_key[count1] = alpha;
			else if(alpha_key[count1]!=alpha) alpha_key[count1] = 0.5;

			if(alpha<=Zero) continue;  //等于零值，没有长程效果

			//--------------------------------------------
			//计算长程作用衰减函数值
			const double x = gaupoi_right.x-gaupoi_left.x;
			const double y = gaupoi_right.y-gaupoi_left.y;

			const double dis_squr = x*x + y*y ;
			const double poi_dis = sqrt(dis_squr);		//两点之间的距离

			if(poi_dis>peri_para.horizon_R||poi_dis<Zero) continue;		//圆形积分域

			//-------------------------------------------------
			//Evaluation: does the bond broke?
			const double n_x = n_gaupoi_right.x-n_gaupoi_left.x;
			const double n_y = n_gaupoi_right.y-n_gaupoi_left.y;

			const double n_dis_squr = n_x*n_x + n_y*n_y;
			const double n_poi_dis = sqrt(n_dis_squr);		//两点之间的距离

			if(n_poi_dis<peri_para.broken_factor*poi_dis) continue;

			//updating data
			long int btco = cogau - 1;
			break_table[btco]=true;
			bk_count++;

			//-------------------------------------------------
			//计算Ｊ矩阵；
			double Jright_val;
			if(flag_right==0) Jright_val = Jacobi;
			else
			{
				//--------------------------------------------------
				//J矩阵
				double Jmatrix_right[2][2];
				//以上两个矩阵的积
				for(int i=0; i<2; i++)
					for(int j=0; j<2; j++)
					{
						Jmatrix_right[i][j]=0;
						for(int k=0; k<4; k++)
							Jmatrix_right[i][j]=Jmatrix_right[i][j] + gauss_dif[count2][i][k]*elenodes_right[k][j];
					}

					//求出J矩阵的行列式
					Jright_val = Jmatrix_right[0][0]*Jmatrix_right[1][1]-Jmatrix_right[0][1]*Jmatrix_right[1][0];
			}

			//--------------------------------------------
			//Calculate the matrix contribution
			double cos2pha;
			if(fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
			else cos2pha = x*x/(x*x + y*y);					//cos(pha)^2

			double sum = peri_para.acoe[0] + peri_para.acoe[1]*(2*cos2pha-1) + peri_para.acoe[2]*(8*cos2pha*cos2pha-8*cos2pha+1);

			const double gv = exp(-poi_dis/peri_para.intrinsic_L)*sum*alpha*Jright_val*weight[count2]; //衰减以及权重函数值
			Dissip_ener += exp(-poi_dis/peri_para.intrinsic_L)*sum*alpha*weight[count2]*(peri_para.broken_factor-1.0)*(peri_para.broken_factor-1.0)*0.25*weight[count1]*dis_squr*dis_squr*Jright_val*Jleft_val;

			const double r_comp[6] = {x*x, y*y, 0, x*y, y*0, 0*x}; //等效长程力衰减计算
			const double temp_gmat[6] = {gv*r_comp[0], gv*r_comp[1], gv*r_comp[2], gv*r_comp[3], gv*r_comp[4], gv*r_comp[5]}; //衰减函数矩阵的对称项

			for(int i=0; i<6; i++) Gmat[i] += temp_gmat[i];  //关于count2叠加

			for(int i=0; i<6; i++)
				for(int j=0; j<4; j++) GNmat[i][j] +=  temp_gmat[i]*gauss_ns[count2][j]; //关于count2并与形函数乘积叠加
		}

		if(bk_count>0)
		{
			broken_num += bk_count;
			//------------------------------------------------------------------------------------------------------------------------
			int ni = 0, nj = 0;
			double temxy = 0, nxy[6] = {0};
			for(int i=0; i<4; i++)
			{
				nj = 0;
				for(int j=0; j<4; j++)
				{
					temxy = gaus_nwct1[i]*gauss_ns[count1][j];
					for(int k=0; k<6; k++)
						nxy[k] = temxy*Gmat[k];

					element_stiff_matrix1[ni][nj] += nxy[0];
					element_stiff_matrix1[ni+1][nj+1] += nxy[1];
					element_stiff_matrix1[ni+2][nj+2] += nxy[2];
					element_stiff_matrix1[ni][nj+1] += nxy[3];
					element_stiff_matrix1[ni+1][nj] += nxy[3];
					element_stiff_matrix1[ni+1][nj+2] += nxy[4];
					element_stiff_matrix1[ni+2][nj+1] += nxy[4];
					element_stiff_matrix1[ni][nj+2] += nxy[5];
					element_stiff_matrix1[ni+2][nj] += nxy[5];

					for(int k=0; k<6; k++)
						nxy[k] = gaus_nwct1[i]*GNmat[k][j];

					element_stiff_matrix2[ni][nj] -= nxy[0];
					element_stiff_matrix2[ni+1][nj+1] -= nxy[1];
					element_stiff_matrix2[ni+2][nj+2] -= nxy[2];
					element_stiff_matrix2[ni][nj+1] -= nxy[3];
					element_stiff_matrix2[ni+1][nj] -= nxy[3];
					element_stiff_matrix2[ni+1][nj+2] -= nxy[4];
					element_stiff_matrix2[ni+2][nj+1] -= nxy[4];
					element_stiff_matrix2[ni][nj+2] -= nxy[5];
					element_stiff_matrix2[ni+2][nj] -= nxy[5];

					nj += 3;
				}
				ni += 3;
			}
		}
	}
}
//-----------------------------------------------------------------------------------------------
//Update element matrices with broken bonds (after up to damage criterion) for 2D problem (New Algorithm)
int Global_Stiff_Matrix::Update_nonlocal_matrices_2D(const int &gaussnum, const struct Damage damages, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func,
																	                              const vector<Node> &nodes, const vector<vector<double> > &Gp_val, const vector<vector<double> > &Ak_val, const vector<double> &U_s, vector<Element> &elements,
																								  const vector<vector<double> > &damage_table, vector<MathMatrix> &ele_self_matrix, vector<vector<MathMatrix> > &ele_relative_matrix, 
																								  vector<bool> *break_table, vector<double> &dissip_en, int &broken_sum, vector<bool> &break_iter)const
{
	//---------------------------------------------------------------------------
	//Generating gaussian points
	Gauss gau;		//gaussian potins
	if(gau.Generate_gauss_2D(gaussnum)==0) return 0;

	//------------------------------------------------------------------------
	//Initialization of variables
	cout << "-_- Initializing variables..." << endl;
	hout << "-_- Initializing variables..." << endl;

	const int GS = (int)gau.gauss.size();
	const int ES = (int)elements.size();
	double (*gauss_ns)[4] = new double [GS][4];				//记录单元高斯点的形函数
	double (*gauss_nw)[4] = new double [GS][4];				//记录单元高斯点的形函数(带权重系数)
	double (*gauss_dif)[2][4] = new double [GS][2][4];		//记录单元高斯点形函数的导数
	double Jacobi = 0;																//标准正方形单元的雅可比值
	double (*gauss_po)[3] = new double [GS][3];				//标准正方形单元的高斯点坐标
	double (*ele_cent)[4] = new double [ES][4];					//记录单元中心点位置向量(x,y,z)分别放于[0],[1]和[2]中（以标准正方形单元中心点为原点）, [3]中放单元顶点权重的平均值
	
	//Calculate data of every gaussian point of element for 2D problem
	Generate_element_gauss_data_2D(elements, nodes, com_mod, weight_func, gau.gauss, gau.weight, gauss_ns, gauss_nw, gauss_dif, Jacobi, gauss_po, ele_cent);
	
	cout << "^_^ Variables initialized successfully!" << endl << endl;
	hout << "^_^ Variables initialized successfully!" << endl << endl;

	//------------------------------------------------------------------------
	//Loops all elements
	cout << "-_- Updating the element matrices in nonlocal part..." << endl;
	hout << "-_- Updating the element matrices in nonlocal part..." << endl;
	//执行openmp
	#pragma omp parallel
	{
		double element_stiff_matrix1[12][12];
		double element_stiff_matrix2[12][12];
		
		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<ES; i++) 
		{
			//Output calculating progress
			#pragma omp critical
			if(i%(ES/10)==0)	cout << i*100/ES <<"% done"<<endl;
			
			//---------------------------------------------------------------------------
			//Extracting coordinates of nodes in the principal element 
			double elenodes[4][2];
			for(int j=0; j<4; j++)
			{
				elenodes[j][0] = nodes[elements[i].nodes_id[j]].x;
				elenodes[j][1] = nodes[elements[i].nodes_id[j]].y;
			}
			//Extracting coordinates of nodes in the updated principal element with displacement solutions
			double new_elenodes[4][2];
			for(int j=0; j<4; j++)
			{
				new_elenodes[j][0] = elenodes[j][0] + U_s[3*elements[i].nodes_id[j]];
				new_elenodes[j][1] = elenodes[j][1] + U_s[3*elements[i].nodes_id[j]+1];
			}

			//Defining parameters of the center of the principal element
			double elec_left[4] = { ele_cent[i][0], ele_cent[i][1], ele_cent[i][2], ele_cent[i][3] };

			//------------------------------------------------------------------------
			//Accumulate the number of relative elements (less than i) //attention: Openmp Parallel
			long int ref_tot = 0;
			for(int j=0; j<i; j++) ref_tot += elements[j].relative_eles.size();

			//------------------------------------------------------------------------
			//Loops element neighbours
			for(int j=0; j<(int)elements[i].relative_eles.size(); j++)
			{
				ref_tot++; //Accumlation (注意不能移动到下面 因为以下有continue语句)
				//--------------------------------------------------
				//Defining parameters of the center of the associated element 
				const int ere = elements[i].relative_eles[j]; 
				double elec_right[4] = { ele_cent[ere][0], ele_cent[ere][1], ele_cent[ere][2], ele_cent[ere][3] };

				//--------------------------------------------------
				//The bond inside the element cannot be broken.
				if(i==ere) continue;

				//--------------------------------------------------
				//The bonds can be broken only between two elements at least one of which is a fully damaged element.
				if(!(elements[i].DG_elem||elements[ere].DG_elem)) continue;
//				if(!(elements[i].DG_elem&&elements[ere].DG_elem)) continue;

				//--------------------------------------------------
				//A preliminary examination for skipping both elements completely in the local continuum area
				//（注意elec_left[3]和elec_right[3]中分别是两个四边形顶点权重值得平均）
				if(elec_left[3]<=Zero&&elec_right[3]<=Zero) continue;

				//--------------------------------------------------
				//Extracting coordinates of nodes in the associated element
				double elenodes_relative[4][2];
				for(int k=0; k<4; k++)
				{
					elenodes_relative[k][0] = nodes[elements[ere].nodes_id[k]].x;
					elenodes_relative[k][1] = nodes[elements[ere].nodes_id[k]].y;
				}
				//Extracting coordinates of nodes in the updated associated element with displacement solutions
				//extracting new relative element coordinates
				double new_elenodes_relative[4][2];
				for(int k=0; k<4; k++)
				{
					new_elenodes_relative[k][0] = elenodes_relative[k][0] + U_s[3*elements[ere].nodes_id[k]];
					new_elenodes_relative[k][1] = elenodes_relative[k][1] + U_s[3*elements[ere].nodes_id[k]+1];
				}

				//--------------------------------------------------
				//Calculate the position of break_table
				const long int pos_bt = (ref_tot-1)*GS*GS;
				//Define a record for updated broken bonds
				int broken_num = 0;
				int increased_num = 0;
				double Dissip_ener = 0.0;
				//------------------------------------------------------------------------
				//Generate the element matrices based on non-local continuum model (long-range forces) judging broken bonds for 2D problem
				Longforce_Elestiff_Breaks_2D(element_stiff_matrix1, element_stiff_matrix2, break_table[0], break_table[1], increased_num, broken_num, damages, damage_table, Dissip_ener,
																weight_func, elements[i].flag, elements[ere].flag, elenodes, new_elenodes, elenodes_relative, new_elenodes_relative, gauss_ns, gauss_nw,
																gauss_dif, gau.weight, peri_para, Jacobi, gauss_po, elec_left, elec_right, pos_bt, Gp_val, Ak_val, i, ere);
					

				//Updating data
				#pragma omp critical	
				{
					if(increased_num>0)	broken_sum += increased_num;  //注意到最后可能会出现少量bond一步一步又恢复的不收敛情况
					if(broken_num>0)
					{
						dissip_en[i] += Dissip_ener;

						break_iter[i]=true;
						break_iter[ere]=true;

						for(int k=0; k<12; k++)
							for(int m=0; m<12; m++)
							{
								ele_self_matrix[i].element[k][m] -= element_stiff_matrix1[k][m];
								ele_relative_matrix[i][j].element[k][m] -= element_stiff_matrix2[k][m];
							}
					}
				}
			}
		}
	}

	//---------------------------------------------------------------------------
	//Delete pointers
	delete[] gauss_ns;
	delete[] gauss_nw;
	delete[] gauss_dif;
	delete[] gauss_po;
	delete[] ele_cent;

	cout << "^_^ The element matrices updated successfully!" << endl;
	hout << "^_^ The element matrices updated successfully!" << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------
//Generate the element matrix based on non-local continuum model (long-range forces) with broken bonds (after up to damage criterion) for 2D problem (New Algorithm)
void Global_Stiff_Matrix::Longforce_Elestiff_Breaks_2D(double (*element_stiff_matrix1)[12], double (*element_stiff_matrix2)[12], const vector<bool> &break_table0, vector<bool> &break_table1, int &increased_num, int &broken_num,
																						const struct Damage damages, const vector<vector<double> > &damage_table, double &Dissip_ener, const struct Weight_func &weight_func, const int &flag_left, const int &flag_right,
																						const double (*elenodes_left)[2], const double(*n_ele_l)[2], const double (*elenodes_right)[2], const double(*n_ele_r)[2], const double (*gauss_ns)[4], const double (*gauss_nw)[4], 
																						const double (*gauss_dif)[2][4], const vector<double> &weight, const struct Peri_para &peri_para,const double &Jacobi, const double (*gauss_po)[3], const double elec_left[], 
																						const double elec_right[], const long int &pos_bt, const vector<vector<double> > &Gp_val, const vector<vector<double> > &Ak_val, const int &pri_ele,const int &ass_ele)const
{
	long int cogau = pos_bt;
	//--------------------------------------------------
	//Clearing
	for(int i=0; i<12; i++) 
		for (int j=0; j<12; j++)
		{
			element_stiff_matrix1[i][j] = 0;
			element_stiff_matrix2[i][j] = 0;
		}

	//--------------------------------------------------	
	//Loop the gaussian points in the current element
	const int gau_num = (int)weight.size();
	for(int count1=0; count1<gau_num; count1++)
	{
		//--------------------------------------------------
		int bk_count = 0;
		double Gmat[6] = {0};
		double GNmat[6][4] = { {0}, {0}, {0}, {0}, {0}, {0} };
		//--------------------------------------------------	
		//Coordinates of gaussian point (left)
		Point_3D gaupoi_left(0, 0, 0);
		if(flag_left==0)
		{
			gaupoi_left.x = gauss_po[count1][0] + elec_left[0];
			gaupoi_left.y = gauss_po[count1][1] + elec_left[1];
		}
		else
		{
			for(int i=0; i<4; i++) 
			{
				gaupoi_left.x += gauss_ns[count1][i]*elenodes_left[i][0];
				gaupoi_left.y += gauss_ns[count1][i]*elenodes_left[i][1];
			}
		}
		//--------------------------------------------------	
		//New coordinates of gaussian point (left)
		Point_3D n_gaupoi_left(0, 0, 0);
		for(int i=0; i<4; i++) 
		{
			n_gaupoi_left.x += gauss_ns[count1][i]*n_ele_l[i][0];
			n_gaupoi_left.y += gauss_ns[count1][i]*n_ele_l[i][1];
		}

		//--------------------------------------------------
		//Evaluate Jacobi matix
		double Jleft_val;
		if(flag_left==0) Jleft_val = Jacobi;
		else
		{
			//Jacobi matix
			double Jmatrix_left[2][2];
			//以上两个矩阵的积
			for(int i=0; i<2; i++)
				for(int j=0; j<2; j++)
				{
					Jmatrix_left[i][j]=0;
					for(int k=0; k<4; k++)
						Jmatrix_left[i][j] += gauss_dif[count1][i][k]*elenodes_left[k][j];
				}

			//求出J矩阵的行列式
			Jleft_val = Jmatrix_left[0][0]*Jmatrix_left[1][1]-Jmatrix_left[0][1]*Jmatrix_left[1][0];	
		}
		
		//--------------------------------------------------
		//Evaluate shape-function values with weighting
		double gaus_nwct1[4];
		for(int i=0; i<4; i++) gaus_nwct1[i] = gauss_nw[count1][i]*Jleft_val;

		//------------------------------------------------------------------------------------------------------------------------
		//iteration on second ghost point
		for(int count2=0; count2<gau_num; count2++)
		{
			//--------------------------------------------------
			if(break_table0[cogau++]) continue;  //First use and then ++, same to: { break_table[cogau] continue; cogau++; } 
			
			//--------------------------------------------
			//Coordinates of gaussian point (right)
			Point_3D gaupoi_right(0, 0, 0);
			if(flag_right==0)
			{
				gaupoi_right.x = gauss_po[count2][0] + elec_right[0];
				gaupoi_right.y = gauss_po[count2][1] + elec_right[1];
			}
			else
			{
				for(int i=0; i<4; i++) 
				{
					gaupoi_right.x +=  gauss_ns[count2][i]*elenodes_right[i][0];
					gaupoi_right.y +=  gauss_ns[count2][i]*elenodes_right[i][1];
				}
			}

			//--------------------------------------------
			//New coordinates of gaussian point (right)
			Point_3D n_gaupoi_right(0, 0, 0);
			for(int i=0; i<4; i++) 
			{
				n_gaupoi_right.x +=  gauss_ns[count2][i]*n_ele_r[i][0];
				n_gaupoi_right.y +=  gauss_ns[count2][i]*n_ele_r[i][1];
			}

			//--------------------------------------------
			//计算权重函数值
			double alpha = 0.5*(Gp_val[pri_ele][count1]+Gp_val[ass_ele][count2]);
			if (fabs(alpha - 1.0) <= Zero) alpha = 1.0;  //为下面精确赋值
			else if (alpha <= Zero) alpha = 0.0;

			if(alpha<=Zero) continue;  //等于零值，没有长程效果

			//--------------------------------------------
			//计算长程作用衰减函数值
			const double x = gaupoi_right.x-gaupoi_left.x;
			const double y = gaupoi_right.y-gaupoi_left.y;

			const double dis_squr = x*x + y*y ;
			const double poi_dis = sqrt(dis_squr);		//两点之间的距离

			if(poi_dis>peri_para.horizon_R||poi_dis<Zero) continue;		//圆形积分域

			//-------------------------------------------------
			//Evaluation: does the bond broke?
			const double n_x = n_gaupoi_right.x-n_gaupoi_left.x;
			const double n_y = n_gaupoi_right.y-n_gaupoi_left.y;

			const double n_dis_squr = n_x*n_x + n_y*n_y;
			const double n_poi_dis = sqrt(n_dis_squr);		//两点之间的距离

			long int btco = cogau - 1;
			if(n_poi_dis<peri_para.broken_factor*poi_dis&&(!break_table1[btco])) continue;

			if(!break_table1[btco])
			{
				break_table1[btco]=true;  //记录从联结到断的bond
				increased_num++;				//记录状态变化的bond的个数
			}
			bk_count++;

			//-------------------------------------------------
			//计算Ｊ矩阵；
			double Jright_val;
			if(flag_right==0) Jright_val = Jacobi;
			else
			{
				//--------------------------------------------------
				//J矩阵
				double Jmatrix_right[2][2];
				//以上两个矩阵的积
				for(int i=0; i<2; i++)
					for(int j=0; j<2; j++)
					{
						Jmatrix_right[i][j]=0;
						for(int k=0; k<4; k++)
							Jmatrix_right[i][j]=Jmatrix_right[i][j] + gauss_dif[count2][i][k]*elenodes_right[k][j];
					}

					//求出J矩阵的行列式
					Jright_val = Jmatrix_right[0][0]*Jmatrix_right[1][1]-Jmatrix_right[0][1]*Jmatrix_right[1][0];
			}

			//--------------------------------------------
			//Calculate the matrix contribution
			double cos2pha;
			if(fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
			else cos2pha = x*x/(x*x + y*y);					//cos(pha)^2

			double sum = peri_para.acoe[0] + peri_para.acoe[1]*(2*cos2pha-1) + peri_para.acoe[2]*(8*cos2pha*cos2pha-8*cos2pha+1);

			const double gv = exp(-poi_dis / peri_para.intrinsic_L)*sum*alpha*Jright_val*weight[count2]/dis_squr; //衰减以及权重函数值
			double tem_ener = exp(-poi_dis / peri_para.intrinsic_L)*sum*alpha*weight[count2]*(peri_para.broken_factor - 1.0)*(peri_para.broken_factor - 1.0)*0.25*weight[count1]*dis_squr*Jright_val*Jleft_val;
//			const double gv = sum*alpha*Jright_val*weight[count2]/dis_squr; //constant micro-modulus
//			double tem_ener = sum*alpha*weight[count2]*(peri_para.broken_factor-1.0)*(peri_para.broken_factor-1.0)*0.25*weight[count1]*dis_squr*Jright_val*Jleft_val;
			Dissip_ener += tem_ener;

			const double r_comp[6] = {x*x, y*y, 0, x*y, y*0, 0*x}; //等效长程力衰减计算
			const double temp_gmat[6] = {gv*r_comp[0], gv*r_comp[1], gv*r_comp[2], gv*r_comp[3], gv*r_comp[4], gv*r_comp[5]}; //衰减函数矩阵的对称项

			for(int i=0; i<6; i++) Gmat[i] += temp_gmat[i];  //关于count2叠加

			for(int i=0; i<6; i++)
				for(int j=0; j<4; j++) GNmat[i][j] +=  temp_gmat[i]*gauss_ns[count2][j]; //关于count2并与形函数乘积叠加
		}

		if(bk_count>0)
		{
			broken_num += bk_count;
			//------------------------------------------------------------------------------------------------------------------------
			int ni = 0, nj = 0;
			double temxy = 0, nxy[6] = {0};
			for(int i=0; i<4; i++)
			{
				nj = 0;
				for(int j=0; j<4; j++)
				{
					temxy = gaus_nwct1[i]*gauss_ns[count1][j];
					for(int k=0; k<6; k++)
						nxy[k] = temxy*Gmat[k];

					element_stiff_matrix1[ni][nj] += nxy[0];
					element_stiff_matrix1[ni+1][nj+1] += nxy[1];
					element_stiff_matrix1[ni+2][nj+2] += nxy[2];
					element_stiff_matrix1[ni][nj+1] += nxy[3];
					element_stiff_matrix1[ni+1][nj] += nxy[3];
					element_stiff_matrix1[ni+1][nj+2] += nxy[4];
					element_stiff_matrix1[ni+2][nj+1] += nxy[4];
					element_stiff_matrix1[ni][nj+2] += nxy[5];
					element_stiff_matrix1[ni+2][nj] += nxy[5];

					for(int k=0; k<6; k++)
						nxy[k] = gaus_nwct1[i]*GNmat[k][j];

					element_stiff_matrix2[ni][nj] -= nxy[0];
					element_stiff_matrix2[ni+1][nj+1] -= nxy[1];
					element_stiff_matrix2[ni+2][nj+2] -= nxy[2];
					element_stiff_matrix2[ni][nj+1] -= nxy[3];
					element_stiff_matrix2[ni+1][nj] -= nxy[3];
					element_stiff_matrix2[ni+1][nj+2] -= nxy[4];
					element_stiff_matrix2[ni+2][nj+1] -= nxy[4];
					element_stiff_matrix2[ni][nj+2] -= nxy[5];
					element_stiff_matrix2[ni+2][nj] -= nxy[5];

					nj += 3;
				}
				ni += 3;
			}
		}
	}
}
//-----------------------------------------------------------------------------------------------
//Calculate the total energy available for break in every element in Postprocessor
int Global_Stiff_Matrix::Available_total_break_energy(const int &gaussnum, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func,
																				  const vector<Node> &nodes, const vector<Element> &elements, vector<double> &dissip_en)const
{
	//---------------------------------------------------------------------------
	//Generating gaussian points
	Gauss gau;		//gaussian potins
	if (gau.Generate_gauss_2D(gaussnum) == 0) return 0;

	//------------------------------------------------------------------------
	//Initialization of variables
	cout << "-_- Initializing variables..." << endl;
	hout << "-_- Initializing variables..." << endl;

	const int GS = (int)gau.gauss.size();
	const int ES = (int)elements.size();
	double(*gauss_ns)[4] = new double[GS][4];				//记录单元高斯点的形函数
	double(*gauss_nw)[4] = new double[GS][4];				//记录单元高斯点的形函数(带权重系数)
	double(*gauss_dif)[2][4] = new double[GS][2][4];		//记录单元高斯点形函数的导数
	double Jacobi = 0;																//标准正方形单元的雅可比值
	double(*gauss_po)[3] = new double[GS][3];				//标准正方形单元的高斯点坐标
	double(*ele_cent)[4] = new double[ES][4];					//记录单元中心点位置向量(x,y,z)分别放于[0],[1]和[2]中（以标准正方形单元中心点为原点）, [3]中放单元顶点权重的平均值

	//Calculate data of every gaussian point of element for 2D problem
	Generate_element_gauss_data_2D(elements, nodes, com_mod, weight_func, gau.gauss, gau.weight, gauss_ns, gauss_nw, gauss_dif, Jacobi, gauss_po, ele_cent);

	cout << "^_^ Variables initialized successfully!" << endl << endl;
	hout << "^_^ Variables initialized successfully!" << endl << endl;

	//------------------------------------------------------------------------
	//Loops all elements
	cout << "-_- Updating the element matrices in nonlocal part..." << endl;
	hout << "-_- Updating the element matrices in nonlocal part..." << endl;
	//执行openmp
	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for (int i = 0; i<ES; i++)
		{
			//---------------------------------------------------------------------------
			//Extracting coordinates of nodes in the principal element 
			double elenodes[4][2];
			for (int j = 0; j<4; j++)
			{
				elenodes[j][0] = nodes[elements[i].nodes_id[j]].x;
				elenodes[j][1] = nodes[elements[i].nodes_id[j]].y;
			}

			//---------------------------------------------------------------------------
			//Defining parameters of the center of the principal element
			double elec_left[4] = { ele_cent[i][0], ele_cent[i][1], ele_cent[i][2], ele_cent[i][3] };

			//------------------------------------------------------------------------
			//Accumulate the number of relative elements (less than i) //attention: Openmp Parallel
			long int ref_tot = 0;
			for (int j = 0; j<i; j++) ref_tot += elements[j].relative_eles.size();

			//------------------------------------------------------------------------
			//Loops element neighbours
			double Ava_dissip_ener = 0.0;
			for (int j = 0; j<(int)elements[i].relative_eles.size(); j++)
			{
				ref_tot++; //Accumlation
				//--------------------------------------------------
				//Defining parameters of the center of the associated element 
				const int ere = elements[i].relative_eles[j];
				double elec_right[4] = { ele_cent[ere][0], ele_cent[ere][1], ele_cent[ere][2], ele_cent[ere][3] };

				//--------------------------------------------------
				//Extracting coordinates of nodes in the associated element
				double elenodes_relative[4][2];
				for (int k = 0; k<4; k++)
				{
					elenodes_relative[k][0] = nodes[elements[ere].nodes_id[k]].x;
					elenodes_relative[k][1] = nodes[elements[ere].nodes_id[k]].y;
				}

				//------------------------------------------------------------------------
				//Generate the element matrices based on non-local continuum model (long-range forces) judging broken bonds for 2D problem
				Available_Breaks_Energy_2D(Ava_dissip_ener, elements[i].flag, elements[ere].flag, elenodes, elenodes_relative,
															  gauss_ns, gauss_nw, gauss_dif, gau.weight, peri_para, Jacobi, gauss_po, elec_left, elec_right);
			}
			#pragma omp critical	
			{
				//For those elements with very large size but far from the damage zone (to reduce the computational cost), 
				//the nonlocal dissipated energy is calculated to be zero because of any bond between two quadrature (gaussian) points is longer than the horizon.
				//To avoid the singularity when calculate the nonlocal damage at any element, so we enforce the dissp_en of this elelment equal to 1,
				//it doesn't change the damage results, because these elements are never damaged. 
				if(Ava_dissip_ener==0.0) dissip_en[i] = 1.0; 
				else dissip_en[i] = Ava_dissip_ener;
			}
		}
	}

	//---------------------------------------------------------------------------
	//Delete pointers
	delete[] gauss_ns;
	delete[] gauss_nw;
	delete[] gauss_dif;
	delete[] gauss_po;
	delete[] ele_cent;

	cout << "^_^ The available break energy in every element updated successfully!" << endl;
	hout << "^_^ The available break energy in every element updated successfully!" << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------
//Available break energy in every element
void Global_Stiff_Matrix::Available_Breaks_Energy_2D(double &Ava_dissip_ener, const int &flag_left, const int &flag_right, const double(*elenodes_left)[2],
																					 const double(*elenodes_right)[2], const double(*gauss_ns)[4], const double(*gauss_nw)[4], 
																					 const double(*gauss_dif)[2][4], const vector<double> &weight, const struct Peri_para &peri_para, const double &Jacobi,
																					 const double(*gauss_po)[3], const double elec_left[], const double elec_right[])const
{
	//--------------------------------------------------	
	//Loop the gaussian points in the current element
	const int gau_num = (int)weight.size();
	for (int count1 = 0; count1<gau_num; count1++)
	{
		//--------------------------------------------------	
		//Coordinates of gaussian point (left)
		Point_3D gaupoi_left(0, 0, 0);
		if (flag_left == 0)
		{
			gaupoi_left.x = gauss_po[count1][0] + elec_left[0];
			gaupoi_left.y = gauss_po[count1][1] + elec_left[1];
		}
		else
		{
			for (int i = 0; i<4; i++)
			{
				gaupoi_left.x += gauss_ns[count1][i] * elenodes_left[i][0];
				gaupoi_left.y += gauss_ns[count1][i] * elenodes_left[i][1];
			}
		}

		//--------------------------------------------------
		//Evaluate Jacobi matix
		double Jleft_val;
		if (flag_left == 0) Jleft_val = Jacobi;
		else
		{
			//Jacobi matix
			double Jmatrix_left[2][2];
			//以上两个矩阵的积
			for (int i = 0; i<2; i++)
			for (int j = 0; j<2; j++)
			{
				Jmatrix_left[i][j] = 0;
				for (int k = 0; k<4; k++)
					Jmatrix_left[i][j] += gauss_dif[count1][i][k] * elenodes_left[k][j];
			}

			//求出J矩阵的行列式
			Jleft_val = Jmatrix_left[0][0] * Jmatrix_left[1][1] - Jmatrix_left[0][1] * Jmatrix_left[1][0];
		}

		//------------------------------------------------------------------------------------------------------------------------
		//iteration on second ghost point
		for (int count2 = 0; count2<gau_num; count2++)
		{
			//--------------------------------------------
			//Coordinates of gaussian point (right)
			Point_3D gaupoi_right(0, 0, 0);
			if (flag_right == 0)
			{
				gaupoi_right.x = gauss_po[count2][0] + elec_right[0];
				gaupoi_right.y = gauss_po[count2][1] + elec_right[1];
			}
			else
			{
				for (int i = 0; i<4; i++)
				{
					gaupoi_right.x += gauss_ns[count2][i] * elenodes_right[i][0];
					gaupoi_right.y += gauss_ns[count2][i] * elenodes_right[i][1];
				}
			}

			//--------------------------------------------
			//计算长程作用衰减函数值
			const double x = gaupoi_right.x - gaupoi_left.x;
			const double y = gaupoi_right.y - gaupoi_left.y;

			const double dis_squr = x*x + y*y;
			const double poi_dis = sqrt(dis_squr);		//两点之间的距离

			if (poi_dis>peri_para.horizon_R||poi_dis<Zero) continue;		//圆形积分域

			//-------------------------------------------------
			//计算Ｊ矩阵；
			double Jright_val;
			if (flag_right == 0) Jright_val = Jacobi;
			else
			{
				//--------------------------------------------------
				//J矩阵
				double Jmatrix_right[2][2];
				//以上两个矩阵的积
				for (int i = 0; i<2; i++)
				for (int j = 0; j<2; j++)
				{
					Jmatrix_right[i][j] = 0;
					for (int k = 0; k<4; k++)
						Jmatrix_right[i][j] = Jmatrix_right[i][j] + gauss_dif[count2][i][k] * elenodes_right[k][j];
				}

				//求出J矩阵的行列式
				Jright_val = Jmatrix_right[0][0] * Jmatrix_right[1][1] - Jmatrix_right[0][1] * Jmatrix_right[1][0];
			}

			//--------------------------------------------
			//Calculate the matrix contribution
			double cos2pha;
			if (fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
			else cos2pha = x*x/(x*x+y*y);					//cos(pha)^2

			double sum = peri_para.acoe[0] + peri_para.acoe[1]*(2*cos2pha-1) + peri_para.acoe[2]*(8*cos2pha*cos2pha-8*cos2pha+1);
//			double  tem_ener = sum*weight[count2]*(peri_para.broken_factor-1.0)*(peri_para.broken_factor-1.0)*0.25*weight[count1]*dis_squr*Jright_val*Jleft_val;	//constant micro-modulus
			double  tem_ener = exp(-poi_dis/peri_para.intrinsic_L)*sum*weight[count2]*(peri_para.broken_factor-1.0)*(peri_para.broken_factor-1.0)*0.25*weight[count1]*dis_squr*Jright_val*Jleft_val;
			Ava_dissip_ener += tem_ener;
		}
	}
}
//-----------------------------------------------------------------------------------------------
//Generate the energy density at every gaussian point based on non-local continuum model judging broken bonds for 2D problem
void Global_Stiff_Matrix::Longforce_Energy_Density_Brokeless_2D(vector<double> &Y_nonlocal_energy, vector<double> &Bond_dissip_en, const double *U_ele_l, const double *U_ele_r, 
										const int &flag_left, const int &flag_right, const double (*elenodes_left)[2], const double (*elenodes_right)[2], const double (*gauss_ns)[4], const double (*gauss_dif)[2][4], 
										const vector<double> &weight, const struct Peri_para &peri_para, const double &Jacobi, const double (*gauss_po)[3], const double elec_left[], const double elec_right[], 
										const vector<bool> *break_table, const long int &pos_bt, const vector<vector<double> > &Gp_val, const int &pri_ele, const int &ass_ele)const
{
	long int cogau = pos_bt;

	//--------------------------------------------------	
	//Loop the gaussian points in the current element
	const int gau_num = (int)weight.size();
	for(int count1=0; count1<gau_num; count1++)
	{
		//--------------------------------------------------	
		//左端高斯点坐标
		Point_3D gaupoi_left(0, 0, 0);
		if(flag_left==0)
		{
			gaupoi_left.x = gauss_po[count1][0] + elec_left[0];
			gaupoi_left.y = gauss_po[count1][1] + elec_left[1];
		}
		else
		{
			for(int i=0; i<4; i++) 
			{
				gaupoi_left.x += gauss_ns[count1][i]*elenodes_left[i][0];
				gaupoi_left.y += gauss_ns[count1][i]*elenodes_left[i][1];
			}
		}

		//--------------------------------------------------
		//The displacements at Gaussian points
		double gaudis_left[3] = {0.0};
		for(int i=0; i<4; i++) 
		{
			gaudis_left[0] += gauss_ns[count1][i]*U_ele_l[3*i];
			gaudis_left[1] += gauss_ns[count1][i]*U_ele_l[3*i+1];
		}

		//--------------------------------------------------
		//Evaluate Jacobi matix
		double Jleft_val;
		if(flag_left==0) Jleft_val = Jacobi;
		else
		{
			//Jacobi matix
			double Jmatrix_left[2][2];
			//以上两个矩阵的积
			for(int i=0; i<2; i++)
				for(int j=0; j<2; j++)
				{
					Jmatrix_left[i][j]=0;
					for(int k=0; k<4; k++)
						Jmatrix_left[i][j] += gauss_dif[count1][i][k]*elenodes_left[k][j];
				}

			//求出J矩阵的行列式
			Jleft_val = Jmatrix_left[0][0]*Jmatrix_left[1][1]-Jmatrix_left[0][1]*Jmatrix_left[1][0];	
		}

		//------------------------------------------------------------------------------------------------------------------------
		//循环外单元高斯点计算积分
		for(int count2=0; count2<gau_num; count2++)
		{
			cogau++;
			//--------------------------------------------
			//计算权重函数值(此段程序放在计算bond的长度之前为节省计算量)
			double alpha = 0.5*(Gp_val[pri_ele][count1] + Gp_val[ass_ele][count2]);
			if (fabs(alpha - 1.0) <= Zero) alpha = 1.0;  //为下面精确赋值
			else if (alpha <= Zero) alpha = 0.0;

			if (alpha <= Zero) continue;  //等于零值，没有长程效果
			//--------------------------------------------
			//右端高斯点坐标
			Point_3D gaupoi_right(0, 0, 0);
			if(flag_right==0)
			{
				gaupoi_right.x = gauss_po[count2][0] + elec_right[0];
				gaupoi_right.y = gauss_po[count2][1] + elec_right[1];
			}
			else
			{
				for(int i=0; i<4; i++) 
				{
					gaupoi_right.x +=  gauss_ns[count2][i]*elenodes_right[i][0];
					gaupoi_right.y +=  gauss_ns[count2][i]*elenodes_right[i][1];
				}
			}

			//--------------------------------------------------
			//The displacements at Gaussian points
			double gaudis_right[3] = {0.0};
			for(int i=0; i<4; i++) 
			{
				gaudis_right[0] += gauss_ns[count2][i]*U_ele_r[3*i];
				gaudis_right[1] += gauss_ns[count2][i]*U_ele_r[3*i+1];
			}

			//--------------------------------------------
			//计算长程作用衰减函数值
			const double x = gaupoi_right.x-gaupoi_left.x;
			const double y = gaupoi_right.y-gaupoi_left.y;

			const double dis_squr = x*x + y*y;
			const double poi_dis = sqrt(dis_squr);		//两点之间的距离
			
			if(poi_dis>peri_para.horizon_R||poi_dis<Zero) continue;		//圆形积分域

			//--------------------------------------------
			//(u_{\xi}(p)-u_{\xi}(x))^2
			const double diff_U_xi = (gaudis_right[0]-gaudis_left[0])*x+(gaudis_right[1]-gaudis_left[1])*y;

			//计算Ｊ矩阵；
			double Jright_val;
			if(flag_right==0) Jright_val = Jacobi;
			else
			{
				//--------------------------------------------------
				//J矩阵
				double Jmatrix_right[2][2];
				//以上两个矩阵的积
				for(int i=0; i<2; i++)
					for(int j=0; j<2; j++)
					{
						Jmatrix_right[i][j]=0;
						for(int k=0; k<4; k++)
							Jmatrix_right[i][j]=Jmatrix_right[i][j] + gauss_dif[count2][i][k]*elenodes_right[k][j];
					}

				//求出J矩阵的行列式
				Jright_val = Jmatrix_right[0][0]*Jmatrix_right[1][1]-Jmatrix_right[0][1]*Jmatrix_right[1][0];
			}
	
			//--------------------------------------------
			//计算长程力矩阵
			double cos2pha;
			if(fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
			else cos2pha = x*x/(x*x + y*y);					//cos(pha)^2

			double sum = peri_para.acoe[0] + peri_para.acoe[1]*(2*cos2pha-1) + peri_para.acoe[2]*(8*cos2pha*cos2pha-8*cos2pha+1);

			Y_nonlocal_energy[count1] += exp(-poi_dis/peri_para.intrinsic_L)*sum*alpha*weight[count2]*0.25*Jright_val*diff_U_xi*diff_U_xi/dis_squr;
//			Y_nonlocal_energy[count1] += sum*alpha*weight[count2]*0.25*Jright_val*diff_U_xi*diff_U_xi/dis_squr;  //constant micro-modulus

			//--------------------------------------------------------------------------------------------
			//计算新断裂bond释放的能量
//			int btco = cogau - 1;
//			if(break_table[0][btco]) continue;   //上一步已经断裂
//			else if(!break_table[1][btco]) continue;		//这一步还没有断裂

//			Bond_dissip_en[count1] += exp(-poi_dis / peri_para.intrinsic_L)*sum*alpha*weight[count2]*(peri_para.broken_factor - 1.0)*(peri_para.broken_factor - 1.0)*0.25*weight[count1]*dis_squr*Jright_val*Jleft_val;
//			Bond_dissip_en[count1] += sum*alpha*weight[count2]*(peri_para.broken_factor-1.0)*(peri_para.broken_factor-1.0)*0.25*weight[count1]*dis_squr*Jright_val*Jleft_val;  //constant micro-modulus
		}
	}
}
//===============================================================
