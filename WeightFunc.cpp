//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	WeightFunc.cpp
//OBJECTIVE:	Estimate the values of points by weighting function
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================

#include "WeightFunc.h"

//---------------------------------------------------------------------------
//Value weight function
double WeightFunc::Value_weight_function(const Point_3D &poi, const struct Weight_func &weight_func)const
{
	double value = 0.0;
	double max_val=0.0;
	for(int i=0; i<weight_func.num; i++)
	{	
		if(weight_func.shape[i]=="Null") continue;
		
		//计算参数
		const double x1 = weight_func.r0[i];
		const double x2 = weight_func.r1[i];
		const double r = weight_func.ratio[i];

		double dx=0, dy=0, dz=0, xx=0;

		//同心球面
		if(weight_func.shape[i]=="Sphere")
		{
			dx = poi.x-weight_func.center[i].x;
			dy = poi.y-weight_func.center[i].y;
			dz = poi.z-weight_func.center[i].z;
			xx = sqrt(dx*dx+dy*dy+dz*dz);
		}
		else if(weight_func.shape[i]=="Cylinder_x")
		{
			dy = poi.y-weight_func.center[i].y;
			dz = poi.z-weight_func.center[i].z;
			xx = sqrt(dy*dy+dz*dz/(r*r));
		}
		else if(weight_func.shape[i]=="Cylinder_y")
		{
			dx = poi.x-weight_func.center[i].x;
			dz = poi.z-weight_func.center[i].z;
			xx = sqrt(dx*dx+dz*dz/(r*r));
		}
		else if(weight_func.shape[i]=="Cylinder_z")
		{
			dx = poi.x-weight_func.center[i].x;
			dy = poi.y-weight_func.center[i].y;
			xx = sqrt(dx*dx+dy*dy/(r*r));
		}
		else 
		{
			cout << "Error: i=" << i << " the shape of weighting area " << weight_func.shape[i] << " is not defined!" << endl;
			hout << "Error: i=" << i << " the shape of weighting area " << weight_func.shape[i] << " is not defined!" << endl;
		}

		if(xx<=x1) value = 1.0;
		else if(xx>=x2) value = 0.0;
		else 
		{
			if(weight_func.func_order[i]=="Constant") value = weight_func.func_constant[i];		//常数权重函数值	
			else if(weight_func.func_order[i]=="Linear") value = 1.0-(xx-x1)/(x2-x1);					//线性权重函数值
			else if(weight_func.func_order[i]=="Cubic") value = 1.0+(xx-x1)*(xx-x1)*(2*xx-3*x2+x1)/((x2-x1)*(x2-x1)*(x2-x1)); //三次幂权重函数值
			else 
			{
				cout << "Error: i=" << i << " the order of weighting function " << weight_func.func_order[i] << " is not defined!" << endl;
				hout << "Error: i=" << i << " the order of weighting function " << weight_func.func_order[i] << " is not defined!" << endl;
			}
		}
		if(value>max_val) max_val=value;
	}
	return max_val;
}
//-----------------------------------------------------------------------------------------------
//Output weighting function value
void WeightFunc::Output_weight_function_value(const string &output_file_name, const struct Weight_func &weight_func)const
{
	ofstream otec(output_file_name.c_str());

	otec << weight_func.keywords << endl;
	otec << weight_func.mark << endl;
	otec << weight_func.num << endl;
	//Attention: i start from number 1, skip over the initial weight_function value
	for(int i=1; i<weight_func.num; i++)
	{
		otec << weight_func.shape[i] << "  ";
		otec << weight_func.center[i].x << "  ";
		otec << weight_func.center[i].y << "  ";
		otec << weight_func.center[i].z << "  ";
		otec << weight_func.r0[i] << "  ";
		otec << weight_func.r1[i] << "  ";
		otec << weight_func.ratio[i] << "  ";
		otec << weight_func.func_order[i] << "  ";
		otec << weight_func.func_constant[i] << endl;
	}
	otec.close(); 
}
//-----------------------------------------------------------------------------------------------
//Intput weighting function value
void WeightFunc::Input_weight_function_value(const string &input_file_name, struct Weight_func &weight_func)
{
	ifstream infile(input_file_name.c_str());
	
	Input IPt;
	istringstream iteck(IPt.Get_Line(infile));
	string keywords;
	iteck >> keywords;

	istringstream itecm(IPt.Get_Line(infile));
	bool mark;
	itecm >> mark;

	istringstream itecw(IPt.Get_Line(infile));
	int wnum;
	itecw >> wnum;
	weight_func.num = wnum;

	string shape, func_order;
	Point_3D center;
	double r0, r1, ratio, func_constant;
	//Attention: i start from number 1, skip over the initial weight_function value
	for(int i=1; i<wnum; i++)
	{
		istringstream itec(IPt.Get_Line(infile));
		itec >> shape >> center.x >> center.y >> center.z >> r0 >> r1 >> ratio >> func_order >> func_constant;
		weight_func.shape.push_back(shape);
		weight_func.center.push_back(center);
		weight_func.r0.push_back(r0);
		weight_func.r1.push_back(r1);
		weight_func.ratio.push_back(ratio);
		weight_func.func_order.push_back(func_order);
		weight_func.func_constant.push_back(func_constant);
	}

	infile.close(); 
}
//-----------------------------------------------------------------------------------------------
//Output the values of weighting function at every node while updating
void WeightFunc::Output_weight_function_contour(const string &output_file_name, vector<Node> &nodes, const vector<Element> &elements, const vector<Node> &gauss,
																						  const vector<double> &weight, const struct Weight_func &weight_func)const
{
	//Define the relative elements of every node
	for(int i=0; i<(int)nodes.size(); i++)
	{
		nodes[i].relative_eles.clear();
	}
	for(int i=0; i<(int)elements.size(); i++)
	{
		int node_size = int(elements[i].nodes_id.size());
		for(int j=0; j<node_size; j++)
		{
			nodes[elements[i].nodes_id[j]].relative_eles.push_back(i);
		}
	}

	//-----------------------------------------------------------------------
	Mesher Mesh;
	//Estimate the weighting value of every elements
	vector<double> weight_val(elements.size(), 0.0);
	for(int i=0; i<(int)elements.size(); i++)
	{
		//---------------------------------------------------
		//单元权重值
		Node elenod[8];
		for(int j=0; j<8; j++) elenod[j] = nodes[elements[i].nodes_id[j]];
		double ele_vol = Mesh.Calculate_brick_volume(elenod);

		double elenodes[8][3];
		for(int j=0; j<8; j++) 
		{
			elenodes[j][0] = elenod[j].x;
			elenodes[j][1] = elenod[j].y;
			elenodes[j][2] = elenod[j].z;
		}

		//---------------------------------------------------
		for(int count=0; count<(int)gauss.size(); count++)
		{
			double Nshape[8] = {0};
			Nshape[0]=0.125*(1.0-gauss[count].x)*(1.0-gauss[count].y)*(1.0-gauss[count].z);
			Nshape[1]=0.125*(1.0+gauss[count].x)*(1.0-gauss[count].y)*(1.0-gauss[count].z);
			Nshape[2]=0.125*(1.0+gauss[count].x)*(1.0+gauss[count].y)*(1.0-gauss[count].z);
			Nshape[3]=0.125*(1.0-gauss[count].x)*(1.0+gauss[count].y)*(1.0-gauss[count].z);
			Nshape[4]=0.125*(1.0-gauss[count].x)*(1.0-gauss[count].y)*(1.0+gauss[count].z);
			Nshape[5]=0.125*(1.0+gauss[count].x)*(1.0-gauss[count].y)*(1.0+gauss[count].z);
			Nshape[6]=0.125*(1.0+gauss[count].x)*(1.0+gauss[count].y)*(1.0+gauss[count].z);
			Nshape[7]=0.125*(1.0-gauss[count].x)*(1.0+gauss[count].y)*(1.0+gauss[count].z);

			//--------------------------------------------------	
			//Coordinates of gaussian points
			Point_3D gaupoi(0, 0, 0);
			for(int j=0; j<8; j++) 
			{
				gaupoi.x += Nshape[j]*elenodes[j][0];
				gaupoi.y += Nshape[j]*elenodes[j][1];
				gaupoi.z += Nshape[j]*elenodes[j][2];
			}

			//--------------------------------------------
			//The weighting function values of gaussian points
			double gauss_wei = Value_weight_function(gaupoi, weight_func);

			//--------------------------------------------
			//EstimateＪmatrix
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
			//J matrix
			double Jmatrix[3][3];
			for(int j=0; j<3; j++)
				for(int k=0; k<3; k++)
				{
					Jmatrix[j][k]=0;
					for(int m=0; m<8; m++)
					Jmatrix[j][k] += diff[j][m]*elenodes[m][k];
				}
			//--------------------------------------------------
			//The determinant of J matrix
			double Jac_val = Jmatrix[0][0]*(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])
										-Jmatrix[0][1]*(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])
										+Jmatrix[0][2]*(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0]);			

			weight_val[i] += gauss_wei*Jac_val*weight[count];
		}
		weight_val[i] = weight_val[i]/ele_vol;
	}

	//Average value of weighting function at every node
	vector<double> nod_val(nodes.size(), 0);
	for(int i=0; i<(int)nodes.size(); i++)
	{
		const int nres = (int)nodes[i].relative_eles.size();
		for(int j=0; j<nres; j++)
			nod_val[i] += weight_val[nodes[i].relative_eles[j]];
		nod_val[i] = nod_val[i]/nres;
	}

	//Clear relative elements info of every node
	for(int i=0; i<(int)nodes.size(); i++)	nodes[i].relative_eles.clear(); 

	//---------------------------------------------------------------------------
	//Outpur resutls
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Weighting_Functing_Contour" << endl;
	otec << "VARIABLES = X, Y, Z, Weight" << endl;
	
	otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size() << ", F=FEPOINT, ET=BRICK" << endl;
	for(int i=0; i<(int)nodes.size(); i++)	
		otec << nodes[i].x << "  " << nodes[i].y << "  " << nodes[i].z << "  " << nod_val[i] << endl;
	otec << endl;
	for (int i=0; i < (int)elements.size(); i++)
	{
		otec	 << elements[i].nodes_id[0]+1 << "  " << elements[i].nodes_id[1]+1 << "  " 
				 << elements[i].nodes_id[2]+1 << "  " << elements[i].nodes_id[3]+1 << "  " 
				 << elements[i].nodes_id[4]+1 << "  " << elements[i].nodes_id[5]+1 << "  " 
				 << elements[i].nodes_id[6]+1 << "  " << elements[i].nodes_id[7]+1 << endl;
	}
	otec.close(); 
}
//-----------------------------------------------------------------------------------------------
//Export the values of weighting function at every node for testing
void WeightFunc::Export_weight_function_contour(const string &output_file_name, vector<Node> &nodes, const vector<Element> &elements, const double (*gauss_ns)[8], const double (*gauss_dif)[3][8], 
																						  const double &Jacobi, const double (*gauss_po)[3], const vector<double> &weight, const struct Weight_func &weight_func, const double (*ele_cent)[4])const
{
	//决定所以节点的相关单元
	for(int i=0; i<(int)nodes.size(); i++)
	{
		nodes[i].relative_eles.clear();
	}
	for(int i=0; i<(int)elements.size(); i++)
	{
		int node_size = int(elements[i].nodes_id.size());
		for(int j=0; j<node_size; j++)
		{
			nodes[elements[i].nodes_id[j]].relative_eles.push_back(i);
		}
	}

	Mesher Mesh;
	//计算每个单元的权重值
	vector<double> weight_val(elements.size(), 0);
	for(int i=0; i<(int)elements.size(); i++)
	{
		//-----------------------------------------------------------------------
		//单元权重值
		Node elenod[8];
		for(int j=0; j<8; j++) elenod[j] = nodes[elements[i].nodes_id[j]];
		double ele_vol = Mesh.Calculate_brick_volume(elenod);

		double elenodes[8][3];
		for(int j=0; j<8; j++) 
		{
			elenodes[j][0] = elenod[j].x;
			elenodes[j][1] = elenod[j].y;
			elenodes[j][2] = elenod[j].z;
		}

		for(int count=0; count<(int)weight.size(); count++)
		{
			//--------------------------------------------------	
			//左端高斯点坐标
			Point_3D gaupoi(0, 0, 0);
			if(elements[i].flag==0)
			{
				gaupoi.x = gauss_po[count][0] + ele_cent[i][0];
				gaupoi.y = gauss_po[count][1] + ele_cent[i][1];
				gaupoi.z = gauss_po[count][2] + ele_cent[i][2];
			}
			else
			{
				for(int i=0; i<8; i++) 
				{
					gaupoi.x += gauss_ns[count][i]*elenodes[i][0];
					gaupoi.y += gauss_ns[count][i]*elenodes[i][1];
					gaupoi.z += gauss_ns[count][i]*elenodes[i][2];
				}
			}

			//--------------------------------------------
			//计算左端高斯点的权重
			double gauss_wei;
			if(ele_cent[i][3]<=Zero) gauss_wei = 0.0;
			else if(fabs(ele_cent[i][3]-1.0)<=Zero) gauss_wei = 1.0;
			else gauss_wei = Value_weight_function(gaupoi, weight_func);

			//计算Ｊ矩阵；
			//--------------------------------------------------
			double Jac_val;
			if(elements[i].flag==0) Jac_val = Jacobi;
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
							Jmatrix_left[i][j] += gauss_dif[count][i][k]*elenodes[k][j];
					}

				//求出J矩阵的行列式
				Jac_val = Jmatrix_left[0][0]*(Jmatrix_left[1][1]*Jmatrix_left[2][2]-Jmatrix_left[1][2]*Jmatrix_left[2][1])
									-Jmatrix_left[0][1]*(Jmatrix_left[1][0]*Jmatrix_left[2][2]-Jmatrix_left[1][2]*Jmatrix_left[2][0])
									+Jmatrix_left[0][2]*(Jmatrix_left[1][0]*Jmatrix_left[2][1]-Jmatrix_left[1][1]*Jmatrix_left[2][0]);	
			}		
			weight_val[i] += gauss_wei*Jac_val*weight[count];
		}
		weight_val[i] = weight_val[i]/ele_vol;
	}
	//计算每个节点的平均权重值
	vector<double> nod_val(nodes.size(), 0);
	for(int i=0; i<(int)nodes.size(); i++)
	{
		const int nres = (int)nodes[i].relative_eles.size();
		for(int j=0; j<nres; j++)
			nod_val[i] += weight_val[nodes[i].relative_eles[j]];
		nod_val[i] = nod_val[i]/nres;
	}

	//清除网格节点的相关单元信息
	for(int i=0; i<(int)nodes.size(); i++)	nodes[i].relative_eles.clear(); 

	//---------------------------------------------------------------------------
	//输出结果
	ofstream otec(output_file_name.c_str());
	otec << "TITLE =" << output_file_name << endl;
	otec << "VARIABLES = X, Y, Z, Weight" << endl;
	
	otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size() << ", F=FEPOINT, ET=BRICK" << endl;
	for(int i=0; i<(int)nodes.size(); i++)	
		otec << nodes[i].x << "  " << nodes[i].y << "  " << nodes[i].z << "  " << nod_val[i] << endl;
	otec << endl;
	for (int i=0; i < (int)elements.size(); i++)
	{
		otec	 << elements[i].nodes_id[0]+1 << "  " << elements[i].nodes_id[1]+1 << "  " 
				 << elements[i].nodes_id[2]+1 << "  " << elements[i].nodes_id[3]+1 << "  " 
				 << elements[i].nodes_id[4]+1 << "  " << elements[i].nodes_id[5]+1 << "  " 
				 << elements[i].nodes_id[6]+1 << "  " << elements[i].nodes_id[7]+1 << endl;
	}
	otec.close(); 
}
//-----------------------------------------------------------------------------------------------
//Calculate the values of weighting function of gaussian points in all elements
void WeightFunc::Value_weightfunc_gausspois_elements(const struct Weight_func &weight_func, const vector<Node> &nodes, const vector<Element> &elements, const int &GS,
																									const double (*gauss_ns)[8], const double (*gauss_po)[3], double (*ele_cent)[4], vector<vector<double> > &Gp_val)const
{
	for(int i=0; i<(int)elements.size(); i++)
	{
		//-----------------------------------------------------------------------
		Node elenod[8];
		for(int j=0; j<8; j++) elenod[j] = nodes[elements[i].nodes_id[j]];
		double elenodes[8][3];
		for(int j=0; j<8; j++) 
		{
			elenodes[j][0] = elenod[j].x;
			elenodes[j][1] = elenod[j].y;
			elenodes[j][2] = elenod[j].z;
		}
		vector<double> gauss_temp;
		for(int count=0; count<GS; count++)
		{
			//--------------------------------------------------	
			//Coordinates of gaussian points
			Point_3D gaupoi(0, 0, 0);
			if(elements[i].flag==0)
			{
				gaupoi.x = gauss_po[count][0] + ele_cent[i][0];
				gaupoi.y = gauss_po[count][1] + ele_cent[i][1];
				gaupoi.z = gauss_po[count][2] + ele_cent[i][2];
			}
			else
			{
				for(int i=0; i<8; i++) 
				{
					gaupoi.x += gauss_ns[count][i]*elenodes[i][0];
					gaupoi.y += gauss_ns[count][i]*elenodes[i][1];
					gaupoi.z += gauss_ns[count][i]*elenodes[i][2];
				}
			}

			//--------------------------------------------
			//The values of weighting function of every gaussian points
			gauss_temp.push_back(Value_weight_function(gaupoi, weight_func));
		}
		Gp_val.push_back(gauss_temp);
	}
}
//-------------------------------------------------------------------------
//Update weighting function: alpha
int WeightFunc::Update_weighting_function_alpha(const vector<Node> &nodes, const vector<Element> &elements, const struct Peri_para &peri_para, const vector<bool> &damaged_iter,
																						vector<bool> &dam_ele_id, struct Weight_func &weight_func, bool &alpha_change)const
{
	//Judging the size of variables
	if(weight_func.num==0||
		(int)weight_func.shape.size()!=weight_func.num||
		(int)weight_func.center.size()!=weight_func.num||
		(int)weight_func.r0.size()!=weight_func.num||
		(int)weight_func.r1.size()!=weight_func.num||
		(int)weight_func.ratio.size()!=weight_func.num||
		(int)weight_func.func_order.size()!=weight_func.num||
		(int)weight_func.func_constant.size()!=weight_func.num)
	{
		cout << "Error: the size of variables in weighting function are different!" << endl; 
		hout << "Error: the size of variables in weighting function are different!" << endl; 
		return 0; 
	}

	for(int i=0; i<(int)elements.size(); i++)
	{
		//Update damage element table
		if(damaged_iter[i]&&!dam_ele_id[i])
		{
			dam_ele_id[i] = true;
			alpha_change = true;

			//Identify element center in the original configuration
			Point_3D scenter(0,0,0);
			for(int j=0; j<8; j++) 
			{
				scenter.x += nodes[elements[i].nodes_id[j]].x;
				scenter.y += nodes[elements[i].nodes_id[j]].y;
				scenter.z += nodes[elements[i].nodes_id[j]].z;
			}
			scenter = scenter/8;

			//Define the parameters of the new weighting function alpha
			weight_func.num++;
			weight_func.shape.push_back("Sphere");
			weight_func.center.push_back(scenter);
			weight_func.r0.push_back(2*peri_para.horizon_R);
			weight_func.r1.push_back(4*peri_para.horizon_R);
			weight_func.ratio.push_back(1.0);
			weight_func.func_order.push_back(weight_func.func_order[0]);
			weight_func.func_constant.push_back(weight_func.func_constant[0]);
		}
	}
	return 1;
}
//-------------------------------------------------------------------------
//Update weighting function and full_damaged_elements
int WeightFunc::Update_weighting_function_alpha(const vector<Node> &nodes, const vector<Element> &elements, const struct Peri_para &peri_para, const struct Damage damages, const vector<bool> &dam_iter,
																			const vector<vector<double> > &damage_table,  vector<bool> &full_dam_eles, struct Weight_func &weight_func, bool &alpha_change)const
{
	//Judging the size of variables
	if(weight_func.num==0||
		(int)weight_func.shape.size()!=weight_func.num||
		(int)weight_func.center.size()!=weight_func.num||
		(int)weight_func.r0.size()!=weight_func.num||
		(int)weight_func.r1.size()!=weight_func.num||
		(int)weight_func.ratio.size()!=weight_func.num||
		(int)weight_func.func_order.size()!=weight_func.num||
		(int)weight_func.func_constant.size()!=weight_func.num)
	{
		cout << "Error: the size of variables in weighting function are different!" << endl; 
		hout << "Error: the size of variables in weighting function are different!" << endl; 
		return 0; 
	}

	for(int i=0; i<(int)elements.size(); i++)
	{
		//Update damaged element table
		if(dam_iter[i]&&!full_dam_eles[i])
		{
			bool temp_dam = true;
			for(int j=0; j<(int)damage_table[i].size(); j++)
			{
				if(damages.d_crit - damage_table[i][j]>Zero)  
				{
					temp_dam = false;
					break;
				}
			}
			if(temp_dam) 
			{
				full_dam_eles[i] = true;
				alpha_change = true;

				//Identify element center in the original configuration
				Point_3D scenter(0,0,0);
				for(int j=0; j<(int)elements[i].nodes_id.size(); j++) 
				{
					scenter.x += nodes[elements[i].nodes_id[j]].x;
					scenter.y += nodes[elements[i].nodes_id[j]].y;
					scenter.z += nodes[elements[i].nodes_id[j]].z;
				}
				scenter = scenter/((int)elements[i].nodes_id.size());

				//Calculate the maximum distance between the center and the vertexes of the element
				double mdver = 0.0;
				for(int j=0; j<(int)elements[i].nodes_id.size(); j++) 
				{
					double disver = scenter.distance_to(nodes[elements[i].nodes_id[j]].x, nodes[elements[i].nodes_id[j]].y, nodes[elements[i].nodes_id[j]].z);
					if(mdver<disver) mdver = disver;	
				}
				
				double dist = peri_para.horizon_R + mdver;
				//Define the parameters of the new weighting function alpha
				weight_func.num++;
				weight_func.shape.push_back("Sphere");
				weight_func.center.push_back(scenter);
				weight_func.r0.push_back(dist);
				weight_func.r1.push_back(dist+peri_para.horizon_R);
				weight_func.ratio.push_back(1.0);
				weight_func.func_order.push_back(weight_func.func_order[0]);
				weight_func.func_constant.push_back(weight_func.func_constant[0]);
			}
		}
	}

	return 1;
}
//-------------------------------------------------------------------------
//Update weighting function: alpha with damaged and broken bonds elements
int WeightFunc::Update_weighting_function_alpha(const vector<Node> &nodes, const vector<Element> &elements, const struct Peri_para &peri_para, const vector<bool> &break_iter,
																						const vector<bool> &dam_iter, vector<bool> &dam_ele_id, struct Weight_func &weight_func, bool &alpha_change)const
{
	//Judging the size of variables
	if(weight_func.num==0||
		(int)weight_func.shape.size()!=weight_func.num||
		(int)weight_func.center.size()!=weight_func.num||
		(int)weight_func.r0.size()!=weight_func.num||
		(int)weight_func.r1.size()!=weight_func.num||
		(int)weight_func.ratio.size()!=weight_func.num||
		(int)weight_func.func_order.size()!=weight_func.num||
		(int)weight_func.func_constant.size()!=weight_func.num)
	{
		cout << "Error: the size of variables in weighting function are different!" << endl; 
		hout << "Error: the size of variables in weighting function are different!" << endl; 
		return 0; 
	}

	for(int i=0; i<(int)elements.size(); i++)
	{
		//Update damaged element table with damged element and broken bonds
		if((dam_iter[i]||break_iter[i])&&!dam_ele_id[i])
		{
			dam_ele_id[i] = true;
			alpha_change = true;

			//Identify element center in the original configuration
			Point_3D scenter(0,0,0);
			for(int j=0; j<8; j++) 
			{
				scenter.x += nodes[elements[i].nodes_id[j]].x;
				scenter.y += nodes[elements[i].nodes_id[j]].y;
				scenter.z += nodes[elements[i].nodes_id[j]].z;
			}
			scenter = scenter/8;

			//Define the parameters of the new weighting function alpha
			weight_func.num++;
			weight_func.shape.push_back("Sphere");
			weight_func.center.push_back(scenter);
			weight_func.r0.push_back(2*peri_para.horizon_R);
			weight_func.r1.push_back(4*peri_para.horizon_R);
			weight_func.ratio.push_back(1.0);
			weight_func.func_order.push_back(weight_func.func_order[0]);
			weight_func.func_constant.push_back(weight_func.func_constant[0]);
		}
	}
	return 1;
}
//-------------------------------------------------------------------------
//Update full_damaged_elements and put weighting function inside completely damaged zone
int WeightFunc::Update_weighting_func_in_damaged_zone(const vector<Node> &nodes, const vector<Element> &elements, const struct Peri_para &peri_para, const struct Damage &damages, const vector<bool> &dam_iter,
																										const vector<vector<double> > &damage_table,  vector<bool> &full_dam_eles, struct Weight_func &weight_func, vector<bool> &weighted_eles, bool &alpha_change)const
{
	//Judging the size of variables
	if(weight_func.num==0||
		(int)weight_func.shape.size()!=weight_func.num||
		(int)weight_func.center.size()!=weight_func.num||
		(int)weight_func.r0.size()!=weight_func.num||
		(int)weight_func.r1.size()!=weight_func.num||
		(int)weight_func.ratio.size()!=weight_func.num||
		(int)weight_func.func_order.size()!=weight_func.num||
		(int)weight_func.func_constant.size()!=weight_func.num)
	{
		cout << "Error: the size of variables in weighting function are different!" << endl; 
		hout << "Error: the size of variables in weighting function are different!" << endl; 
		return 0; 
	}

	vector<Point_3D> centrele(elements.size());
	//Update full damaged element table
	for(int i=0; i<(int)elements.size(); i++)
	{
		if(dam_iter[i]&&!full_dam_eles[i])
		{
			bool temp_dam = true;
			for(int j=0; j<(int)damage_table[i].size(); j++)
			{
				if(damages.d_crit - damage_table[i][j]>Zero)  
				{
					temp_dam = false;
					break;
				}
			}
			if(temp_dam) full_dam_eles[i] = true;
		}

		Point_3D poi_tem(0,0,0);
		int ensize = elements[i].nodes_id.size();
		for(int j=0; j<ensize; j++)
		{
			poi_tem.x += nodes[elements[i].nodes_id[j]].x;
			poi_tem.y += nodes[elements[i].nodes_id[j]].y;
			poi_tem.z += nodes[elements[i].nodes_id[j]].z;
		}
		centrele[i] = poi_tem/ensize;
	}

	//Update weighting function in the full damaged zone
	double dist = 5*peri_para.horizon_R; //4 used before
	for(int i=0; i<(int)elements.size(); i++)
	{
		if(full_dam_eles[i]&&(!weighted_eles[i]))
		{
			bool temp_key = true;
			int count = 0;
			for(int j=0; j<(int)elements.size(); j++)
			{
//				if(i!=j&&centrele[i].distance_to(centrele[j])<dist+Zero)
				if(i!=j&&
					fabs(centrele[j].x-centrele[i].x)<dist+Zero&&
					fabs(centrele[j].y-centrele[i].y)<dist+Zero&&
					fabs(centrele[j].z-centrele[i].z)<dist+Zero)
				{
					count++;
					if(!full_dam_eles[j])
					{
						temp_key = false;
						break;
					}
				}
			}
			if(count>0&&temp_key)
			{
				alpha_change = true;
				weighted_eles[i] = true;
				//Define the parameters of the new weighting function alpha
				weight_func.num++;
				weight_func.shape.push_back("Sphere");
				weight_func.center.push_back(centrele[i]);
				weight_func.r0.push_back(2.5*peri_para.horizon_R); //2 used before
				weight_func.r1.push_back(4*peri_para.horizon_R); //3 used before
				weight_func.ratio.push_back(1.0);
				weight_func.func_order.push_back(weight_func.func_order[0]);
				weight_func.func_constant.push_back(weight_func.func_constant[0]);
			}
		}
	}
	
	return 1;
}
//-------------------------------------------------------------------------
//This is a temporary test for pure damge model, at the beginning, dcirt=0.3 for example, when the full_damge_zone propagate to 5delt (for example), change dcrit=1.0 and continue pure damge model
int WeightFunc::Update_dcrit_pure_damage_test(const vector<Node> &nodes, const vector<Element> &elements, const struct Peri_para &peri_para, struct Damage &damages, const vector<bool> &dam_iter,
	const vector<vector<double> > &damage_table, vector<bool> &full_dam_eles, struct Weight_func &weight_func, vector<bool> &weighted_eles, bool &alpha_change)
{
	if (damages.d_crit==1.0) return 1;

	//Judging the size of variables
	if (weight_func.num == 0 ||
		(int)weight_func.shape.size() != weight_func.num ||
		(int)weight_func.center.size() != weight_func.num ||
		(int)weight_func.r0.size() != weight_func.num ||
		(int)weight_func.r1.size() != weight_func.num ||
		(int)weight_func.ratio.size() != weight_func.num ||
		(int)weight_func.func_order.size() != weight_func.num ||
		(int)weight_func.func_constant.size() != weight_func.num)
	{
		cout << "Error: the size of variables in weighting function are different!" << endl;
		hout << "Error: the size of variables in weighting function are different!" << endl;
		return 0;
	}

	vector<Point_3D> centrele(elements.size());
	//Update full damaged element table
	for (int i = 0; i<(int)elements.size(); i++)
	{
		if (dam_iter[i] && !full_dam_eles[i])
		{
			bool temp_dam = true;
			for (int j = 0; j<(int)damage_table[i].size(); j++)
			{
				if (damages.d_crit - damage_table[i][j]>Zero)
				{
					temp_dam = false;
					break;
				}
			}
			if (temp_dam) full_dam_eles[i] = true;
		}

		Point_3D poi_tem(0, 0, 0);
		int ensize = elements[i].nodes_id.size();
		for (int j = 0; j<ensize; j++)
		{
			poi_tem.x += nodes[elements[i].nodes_id[j]].x;
			poi_tem.y += nodes[elements[i].nodes_id[j]].y;
			poi_tem.z += nodes[elements[i].nodes_id[j]].z;
		}
		centrele[i] = poi_tem / ensize;
	}

	//Update weighting function in the full damaged zone
	double dist = 5 * peri_para.horizon_R;
	for (int i = 0; i<(int)elements.size(); i++)
	{
		if (full_dam_eles[i] && (!weighted_eles[i]))
		{
			bool temp_key = true;
			int count = 0;
			for (int j = 0; j<(int)elements.size(); j++)
			{
				//				if(i!=j&&centrele[i].distance_to(centrele[j])<dist+Zero)
				if (i != j&&
					fabs(centrele[j].x - centrele[i].x)<dist + Zero&&
					fabs(centrele[j].y - centrele[i].y)<dist + Zero&&
					fabs(centrele[j].z - centrele[i].z)<dist + Zero)
				{
					count++;
					if (!full_dam_eles[j])
					{
						temp_key = false;
						break;
					}
				}
			}
			if (count>0 && temp_key)
			{
				damages.d_crit = 1.0;
				hout << "    Attation: damages.d_crit has been changed to 1.0!" << endl;
				break;
			}
		}
	}

	return 1;
}
//-----------------------------------------------------------------------------------------------
//Calculate the values of weighting function of gaussian points in all elements for 2D problem
void WeightFunc::Value_weightfunc_gausspois_elements_2D(const struct Weight_func &weight_func, const vector<Node> &nodes, const vector<Element> &elements, const int &GS,
																									const double (*gauss_ns)[4], const double (*gauss_po)[3], const double (*ele_cent)[4], vector<vector<double> > &Gp_val)const
{
	for(int i=0; i<(int)elements.size(); i++)
	{
		//-----------------------------------------------------------------------
		Node elenod[4];
		for(int j=0; j<4; j++) elenod[j] = nodes[elements[i].nodes_id[j]];
		double elenodes[4][2];
		for(int j=0; j<4; j++) 
		{
			elenodes[j][0] = elenod[j].x;
			elenodes[j][1] = elenod[j].y;
		}
		vector<double> gauss_temp;
		for(int count=0; count<GS; count++)
		{
			//--------------------------------------------------	
			//Coordinates of gaussian points
			Point_3D gaupoi(0, 0, 0);
			if(elements[i].flag==0)
			{
				gaupoi.x = gauss_po[count][0] + ele_cent[i][0];
				gaupoi.y = gauss_po[count][1] + ele_cent[i][1];
			}
			else
			{
				for(int i=0; i<4; i++) 
				{
					gaupoi.x += gauss_ns[count][i]*elenodes[i][0];
					gaupoi.y += gauss_ns[count][i]*elenodes[i][1];
				}
			}

			//--------------------------------------------
			//The values of weighting function of every gaussian points
			gauss_temp.push_back(Value_weight_function(gaupoi, weight_func));
		}
		Gp_val.push_back(gauss_temp);
	}
}
//-----------------------------------------------------------------------------------------------
//Output the values of weighting function at every node while updating for 2D problem
void WeightFunc::Output_weight_function_contour_2D(const string &output_file_name, vector<Node> &nodes, const vector<Element> &elements, vector<bool> &full_dam_eles,
																								 const vector<Node> &gauss, const vector<double> &weight, const struct Weight_func &weight_func)const
{
	//Define the relative elements of every node
	for(int i=0; i<(int)nodes.size(); i++)
	{
		nodes[i].relative_eles.clear();
	}
	for(int i=0; i<(int)elements.size(); i++)
	{
		int node_size = int(elements[i].nodes_id.size());
		for(int j=0; j<node_size; j++)
		{
			nodes[elements[i].nodes_id[j]].relative_eles.push_back(i);
		}
	}

	//-----------------------------------------------------------------------
	Mesher Mesh;
	//Estimate the weighting value of every elements
	vector<double> weight_val(elements.size(), 0.0);
	for(int i=0; i<(int)elements.size(); i++)
	{
		//---------------------------------------------------
		//单元权重值
		Node elenod[4];
		for(int j=0; j<4; j++) elenod[j] = nodes[elements[i].nodes_id[j]];
		double ele_area = Mesh.Calculate_quadri_area(elenod);

		double elenodes[4][2];
		for(int j=0; j<4; j++) 
		{
			elenodes[j][0] = elenod[j].x;
			elenodes[j][1] = elenod[j].y;
		}

		//---------------------------------------------------
		for(int count=0; count<(int)gauss.size(); count++)
		{
			double Nshape[4] = {0};
			Nshape[0]=0.25*(1.0-gauss[count].x)*(1.0-gauss[count].y);
			Nshape[1]=0.25*(1.0+gauss[count].x)*(1.0-gauss[count].y);
			Nshape[2]=0.25*(1.0+gauss[count].x)*(1.0+gauss[count].y);
			Nshape[3]=0.25*(1.0-gauss[count].x)*(1.0+gauss[count].y);

			//--------------------------------------------------	
			//Coordinates of gaussian points
			Point_3D gaupoi(0, 0, 0);
			for(int j=0; j<4; j++) 
			{
				gaupoi.x += Nshape[j]*elenodes[j][0];
				gaupoi.y += Nshape[j]*elenodes[j][1];
			}

			//--------------------------------------------
			//The weighting function values of gaussian points
			double gauss_wei = Value_weight_function(gaupoi, weight_func);

			//--------------------------------------------
			//EstimateＪmatrix
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
			//J matrix
			double Jmatrix[2][2];
			for(int j=0; j<2; j++)
				for(int k=0; k<2; k++)
				{
					Jmatrix[j][k]=0;
					for(int m=0; m<4; m++)
					Jmatrix[j][k] += diff[j][m]*elenodes[m][k];
				}
			//--------------------------------------------------
			//The determinant of J matrix
			double Jac_val = Jmatrix[0][0]*Jmatrix[1][1]-Jmatrix[0][1]*Jmatrix[1][0];			

			weight_val[i] += gauss_wei*Jac_val*weight[count];
		}
		weight_val[i] = weight_val[i]/ele_area;
	}

	//Average value of weighting function at every node
	vector<double> nod_val(nodes.size(), 0);
	for(int i=0; i<(int)nodes.size(); i++)
	{
		const int nres = (int)nodes[i].relative_eles.size();
		for(int j=0; j<nres; j++)
			nod_val[i] += weight_val[nodes[i].relative_eles[j]];
		nod_val[i] = nod_val[i]/nres;
	}

	//Clear relative elements info of every node
	for(int i=0; i<(int)nodes.size(); i++)	nodes[i].relative_eles.clear(); 

	//---------------------------------------------------------------------------
	//Outpur resutls
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Weighting_Functing_Contour" << endl;
	otec << "VARIABLES = X, Y, Weight" << endl;
	
	otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size() << ", F=FEPOINT, ET=QUADRILATERAL" << endl;
	for(int i=0; i<(int)nodes.size(); i++)	
		otec << nodes[i].x << "  " << nodes[i].y << "  " << nod_val[i] << endl;
	otec << endl;
	for (int i=0; i<(int)elements.size(); i++)
	{
		otec	 << elements[i].nodes_id[0]+1 << "  " << elements[i].nodes_id[1]+1 << "  " 
				 << elements[i].nodes_id[2]+1 << "  " << elements[i].nodes_id[3]+1 << endl;
	}

	//---------------------------------------------------------------------------
	//Outpur full damaged elements
	int count =0;
	for (int i=0; i<(int)elements.size(); i++) 
		if(full_dam_eles[i]) count++;
	if(count>0)
	{
		otec << "ZONE N=" << (int)nodes.size() << ", E=" << count << ", F=FEPOINT, ET=QUADRILATERAL" << endl;
		for(int i=0; i<(int)nodes.size(); i++)	
			otec << nodes[i].x << "  " << nodes[i].y << "  " << nod_val[i] << endl;
		otec << endl;
		for (int i=0; i<(int)elements.size(); i++)
		{
			if(full_dam_eles[i])
			{
				otec	 << elements[i].nodes_id[0]+1 << "  " << elements[i].nodes_id[1]+1 << "  " 
						 << elements[i].nodes_id[2]+1 << "  " << elements[i].nodes_id[3]+1 << endl;
			}
		}
	}

	otec.close(); 
}
//-----------------------------------------------------------------------------------------------
//Output the weighting function vectors for 2D problem
int WeightFunc::Output_weight_function_vectors(const string &output_file_name, const struct Weight_func &weight_func)const
{
	//Judging the size of variables
	if(weight_func.num==0||
		(int)weight_func.shape.size()!=weight_func.num||
		(int)weight_func.center.size()!=weight_func.num||
		(int)weight_func.r0.size()!=weight_func.num||
		(int)weight_func.r1.size()!=weight_func.num||
		(int)weight_func.ratio.size()!=weight_func.num||
		(int)weight_func.func_order.size()!=weight_func.num||
		(int)weight_func.func_constant.size()!=weight_func.num)
	{
		cout << "Error: the size of variables in weighting function vectors are different!" << endl; 
		hout << "Error: the size of variables in weighting function vectors are different!" << endl; 
		return 0; 
	}

	//output
	ofstream owfv(output_file_name.c_str());
	owfv << "Weighting function vectors:" << weight_func.num << endl;
	for(int i=0; i<(int)weight_func.num; i++)
	{
		owfv << i << "  " << weight_func.shape[i];
		owfv << "  ("  << weight_func.center[i].x << ", " << weight_func.center[i].y << ", " << weight_func.center[i].z << ")  ";
		owfv << weight_func.r0[i] << "  " << weight_func.r1[i] << "  " << weight_func.ratio[i] << "  ";
		owfv << weight_func.func_order[i] << "  " << weight_func.func_constant[i] << endl;
	}
	owfv.close();
	
	return 1;
}
//===========================================================================
