//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Postprocessor.cpp
//OBJECTIVE:	Implement the postprocessor (the visual output)
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================

#include "Postprocessor.h"

//--------------------------------------------------------------------------
Postprocessor::Postprocessor(const vector<double> &post_u, const vector<Node> &post_nodes, const vector<Element> &post_elements)
{
	u = post_u;
	nodes = post_nodes;
	elements = post_elements;
}
//--------------------------------------------------------------------------
Postprocessor::Postprocessor(const vector<double> &post_u, const vector<double> &post_dissi, const vector<Node> &post_nodes, const vector<Element> &post_elements)
{
	u = post_u;
	dissi = post_dissi;
	nodes = post_nodes;
	elements = post_elements;
}
//--------------------------------------------------------------------------
//Implement postprocessor for 2D problem
int Postprocessor::Treatment_2D(Input *Init, const vector<vector<Node> > &nodes_tot, const vector<vector<Element> > &elements_tot, const vector<MatPro> &mats,
														   const vector<vector<double> > &U_tot, const string &wr_mod)
{
	int Unum = (int)U_tot.size();
	stringstream  wrfw;
	wrfw << "./Results/Data_Backup_Files_Number.dat";
	//Attention: "Init->weight_func" is only used for "write", on the other hand, "postp_weight_func" is used for "read" inside
	if(wr_files_num_weight_funcs(wrfw.str(), wr_mod, Init->weight_func, Unum)==0) return 0;

	if(wr_mod=="Write") 
	{		
		//-------------------------------------------------------------------------------
		//Pre-read for calvulating the total energy
		postp_weight_func = Init->weight_func;
		elements = elements_tot[Unum-1];
		nodes = nodes_tot[Unum-1];
	}
	else if(wr_mod=="Read")
	{
		//-------------------------------------------------------------------------------
		//Pre-read for calvulating the total energy
		stringstream redus;
		if(Unum-1<10)	redus << "./Results/Data_Backup_000" << Unum-1 << ".dat";
		else if (Unum-1<100)	redus << "./Results/Data_Backup_00" << Unum-1 << ".dat";
		else if (Unum-1<1000)	redus << "./Results/Data_Backup_0" << Unum-1 << ".dat";
		else	redus << "./Results/Data_Backup_" << Unum-1 << ".dat";
		if(write_or_read_data(redus.str(), wr_mod)==0) return 0;
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Calculate the total energy available for break in every element
	total_break_energy.assign(elements.size(), 0.0);
	if (Init->mod_disc.mod == "Hybrid")
	{
		//update of non-local element neighbours
		if(wr_mod=="Read")
		{
			Mesher MSR;		//Define the class of Mesher
			if(MSR.Deter_relative_elements(nodes, elements, Init->peri_para.horizon_R, Init->mod_disc.mod, Init->cracks, postp_weight_func)==0) return 0;
		}
		Global_Stiff_Matrix GSM;    //Define the class of Global_Stiff_Matrix
		if(GSM.Available_total_break_energy(Init->gauss.num, Init->peri_para, Init->mod_disc.mod, postp_weight_func, nodes, elements, total_break_energy) == 0) return 0;
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	vector<double> ave_strain_xx(Unum+1, 0.0);
	vector<double> ave_stress_xx(Unum+1, 0.0);
	vector<double> TPB_force_yy(Unum+1, 0.0);
	vector<double> damage_ener(Unum, 0.0);
	vector<double> break_ener(Unum, 0.0);
	for(int i=0; i<Unum; i++)
	{
		cout << endl;
		cout << "Dealing with data from loop " << i << endl;
		cout << "-------------------------------------" << endl;
		hout << endl;
		hout << "Dealing with data from loop " << i << endl;
		hout << "-------------------------------------" << endl;
		if(wr_mod=="Write")	
		{
			u = U_tot[i];
			nodes = nodes_tot[i];
			elements = elements_tot[i];
		}

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//输出或读取位移解
		stringstream worus;
		if(i<10)	worus << "./Results/Data_Backup_000" << i << ".dat";
		else if (i<100)	worus << "./Results/Data_Backup_00" << i << ".dat";
		else if (i<1000)	worus << "./Results/Data_Backup_0" << i << ".dat";
		else	worus << "./Results/Data_Backup_" << i << ".dat";
		if(write_or_read_data(worus.str(), wr_mod)==0) return 0;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//输出位移场云图
		//cout << endl;
		//cout << "    Exporting displacement field..." << endl;
		//hout << endl;
		//hout << "    Exporting displacement field..." << endl;
		//stringstream EDC;
		//if(i<10)	EDC << "./Results/Displacement_Field_Contour_000" << i << ".dat";
		//else if (i<100)	EDC << "./Results/Displacement_Field_Contour_00" << i << ".dat";
		//else if (i<1000)	EDC << "./Results/Displacement_Field_Contour_0" << i << ".dat";
		//else EDC << "./Results/Displacement_Field_Contour_" << i << ".dat";
		//Export_Displacement_Contour_2D(EDC.str(), nodes, elements);

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//计算单元形心处应变向量
		cout << "    Calculating strain..." << endl;
		hout << "    Calculating strain..." << endl;
		Calculate_Elements_Strain_2D(nodes, elements);

		Calculate_Elements_Stress_Energy_2D(mats, elements); //计算单元形心处应力向量和应变能密度

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//Estimate the average stress and strain at the boundary layer elements
		if(Estimate_boundary_stress_strain(nodes, elements, mats, i, Init->geom_rve, ave_strain_xx, ave_stress_xx)==0) return 0;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//Estimate the average force and average displacements at the top center of beam under three point bending
		if(Estimate_boundary_force_2D(nodes, elements, mats, i, Init->force_disp, TPB_force_yy)==0) return 0;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//计算节点处应变向量(Include Damage_energy, Dissipative_energy)
//		Calculate_Nodes_Strain(nodes, elements);
		Calculate_Nodes_Strain_Stress_Energy(Init->mod_disc.mod, Init->damages.d_crit, nodes, elements);  //计算节点处应变向量、应力向量和应变能密度(Include Damage_energy, Dissipative_energy)

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//输出应变场云图(Include Damage_energy, Dissipative_energy)
		cout << "    Exporting strain field..." << endl;
		hout << "    Exporting strain field..." << endl;
		stringstream ESM;
		if(i<10) ESM << "./Results/Strain_Field_Contour_000" << i << ".dat";
		else if (i<100) ESM << "./Results/Strain_Field_Contour_00" << i << ".dat";
		else if (i<1000) ESM << "./Results/Strain_Field_Contour_0" << i << ".dat";
		else ESM << "./Results/Strain_Field_Contour_" << i << ".dat";
//		Export_Strain_Mesh_2D(ESM.str(), nodes, elements);
		Export_Strain_Stress_Energy_Mesh_2D(ESM.str(), nodes, elements);   //输出应变场、应力场和应变能密度云图

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//An accumulation of dissipated energy
		for(int j=0; j<(int)elements.size(); j++)
		{
			damage_ener[i] += elements[j].Damage_energy;
			break_ener[i] += elements[j].Dissipative_energy;
		}
	}
	
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Calculate the sum of average_stress
	double sum_ave_stress = 0.0;
	for(int i=0; i<(int)ave_stress_xx.size(); i++) sum_ave_stress += ave_stress_xx[i];
	hout << endl << "***** The sum of boundary stress: " << sum_ave_stress << " *****" << endl <<endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Export stress strain curve
//	cout<<"    Exporting stess-strain curve..."<<endl;
//	hout<<"    Exporting stess-strain curve..."<<endl;
//	Export_Stress_Strain_Curve("./Results/Stress_strain_curve.dat", ave_strain_xx, ave_stress_xx);

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Export boundary force displacement curve for 2D three point bending
	cout<<"    Exporting force-displacement curve..."<<endl;
	hout<<"    Exporting force-displacement curve..."<<endl;
//	Export_Force_Displacement_Curve("./Results/Force_displacement_curve.dat", TPB_force_yy, Init->force_disp.delta_disp);

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Export dissipative energy
	cout<<"    Exporting dissipative energy..."<<endl;
	hout<<"    Exporting dissipative energy..."<<endl;
	Export_Dissipative_Energy_2D("./Results/Dissipative_energy_plot.dat", damage_ener, break_ener);
	
	return 1;
}
//--------------------------------------------------------------------------
//执行后处理过程
int Postprocessor::Treatment(const vector<vector<Node> > &nodes_tot, const vector<vector<Element> > &elements_tot, const vector<MatPro> &mats,
												   const vector<vector<double> > &U_tot, const string &wr_mod, const vector<vector<double> > &Dissipation)
{
	int Unum = (int)U_tot.size();
	if(wr_mod=="Write") 
	{
		ofstream onum("./Results/Data_Backup_Files_Number.dat");
		onum << Unum << endl;
		onum.close();
	}
	else if(wr_mod=="Read")
	{
		ifstream inum("./Results/Data_Backup_Files_Number.dat");
		inum >> Unum;
		inum.close();
	}

	vector<double> ave_strain_xx(Unum+1, 0.0);
	vector<double> ave_stress_xx(Unum+1, 0.0);
	//-----------------------------------------------------------------------------------------------------------------------------------------
	for(int i=0; i<Unum; i++)
	{
		cout << endl;
		cout << "Dealing with data from loop " << i << endl;
		cout << "-------------------------------------" << endl;
		hout << endl;
		hout << "Dealing with data from loop " << i << endl;
		hout << "-------------------------------------" << endl;
		if(wr_mod=="Write")	
		{
			u = U_tot[i];
			nodes = nodes_tot[i];
			elements = elements_tot[i];
			dissi = Dissipation[i];
		}

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//输出或读取位移解
		stringstream worus;
		if(i<10)	worus << "./Results/Data_Backup_000" << i << ".dat";
		else if (i<100)	worus << "./Results/Data_Backup_00" << i << ".dat";
		else if (i<1000)	worus << "./Results/Data_Backup_0" << i << ".dat";
		else	worus << "./Results/Data_Backup_" << i << ".dat";
		if(write_or_read_data(worus.str(), wr_mod)==0) return 0;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//输出位移场云图
		cout << endl;
		cout << "    Exporting displacement field..." << endl;
		hout << endl;
		hout << "    Exporting displacement field..." << endl;
		stringstream EDC;
		if(i<10)	EDC << "./Results/Displacement_Field_Contour_000" << i << ".dat";
		else if (i<100)	EDC << "./Results/Displacement_Field_Contour_00" << i << ".dat";
		else if (i<1000)	EDC << "./Results/Displacement_Field_Contour_0" << i << ".dat";
		else EDC << "./Results/Displacement_Field_Contour_" << i << ".dat";
		Export_Displacement_Contour(EDC.str(), nodes, elements);

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//输出边界面上位移场云图（内点处位移为0）
		cout << "    Exporting boundarydisplacement field..." << endl;
		hout << "    Exporting boundarydisplacement field..." << endl;
		stringstream EBDC;
		if(i<10) EBDC << "./Results/Boundary_Displacement_Field_Contour_000" << i << ".dat";
		else if (i<100) EBDC << "./Results/Boundary_Displacement_Field_Contour_00" << i << ".dat";
		else if (i<1000) EBDC << "./Results/Boundary_Displacement_Field_Contour_0" << i << ".dat";
		else EBDC <<"./Results/Boundary_Displacement_Field_Contour_"<< i <<".dat";
//		Export_Boundary_Displacement_Contour(EBDC.str(), nodes, elements, 0);

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//输出网格变形图
		cout << "    Exporting deformed mesh..." << endl;
		hout << "    Exporting deformed mesh..." << endl;
		stringstream EDM;
		if(i<10) EDM <<"./Results/Deformed_Mesh_000"<< i <<".dat";
		else if (i<100) EDM <<"./Results/Deformed_Mesh_00"<< i <<".dat";
		else if (i<1000) EDM <<"./Results/Deformed_Mesh_0"<< i <<".dat";
		else EDM <<"./Results/Deformed_Mesh_"<< i <<".dat";
		Export_Deformed_Mesh(EDM.str(), nodes, elements);

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//计算单元形心处应变向量
		cout << "    Calculating strain..." << endl;
		hout << "    Calculating strain..." << endl;
		Calculate_Elements_Strain(nodes, elements);

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//计算节点处应变向量
		Calculate_Nodes_Strain(nodes, elements);

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//输出应变场云图
		cout << "    Exporting strain field..." << endl;
		hout << "    Exporting strain field..." << endl;
		stringstream ESM;
		if(i<10) ESM << "./Results/Strain_Field_Contour_000" << i << ".dat";
		else if (i<100) ESM << "./Results/Strain_Field_Contour_00" << i << ".dat";
		else if (i<1000) ESM << "./Results/Strain_Field_Contour_0" << i << ".dat";
		else ESM << "./Results/Strain_Field_Contour_" << i << ".dat";
		Export_Strain_Mesh(ESM.str(), nodes, elements);

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//坐标变换并输出应变场云图
		//Export_Strain_After_Coordinates_Transformation(nodes, elements);
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Export total dissipative energy
	cout<<"    Exporting dissipative energy..."<<endl;
	hout<<"    Exporting dissipative energy..."<<endl;
	Export_Total_Dissipative_Energy("./Results/Dissipative_energy_plot.dat", Dissipation,  wr_mod);
	
	return 1;
}
//--------------------------------------------------------------------------
//执行后处理过程
int Postprocessor::Treatment(Input *Init, const vector<vector<Node> > &nodes_tot, const vector<vector<Element> > &elements_tot, const vector<MatPro> &mats,
												   const vector<vector<double> > &U_tot, const string &wr_mod, const vector<vector<double> > &Dissipation)
{
	int Unum = (int)U_tot.size();
	if(wr_mod=="Write") 
	{
		ofstream onum("./Results/Data_Backup_Files_Number.dat");
		onum << Unum << endl;
		onum.close();
	}
	else if(wr_mod=="Read")
	{
		ifstream inum("./Results/Data_Backup_Files_Number.dat");
		inum >> Unum;
		inum.close();
	}

	vector<double> ave_strain_xx(Unum+1, 0.0);
	vector<double> ave_stress_xx(Unum+1, 0.0);
	//-----------------------------------------------------------------------------------------------------------------------------------------
	for(int i=0; i<Unum; i++)
	{
		cout << endl;
		cout << "Dealing with data from loop " << i << endl;
		cout << "-------------------------------------" << endl;
		hout << endl;
		hout << "Dealing with data from loop " << i << endl;
		hout << "-------------------------------------" << endl;
		if(wr_mod=="Write")	
		{
			u = U_tot[i];
			nodes = nodes_tot[i];
			elements = elements_tot[i];
			dissi = Dissipation[i];
		}

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//输出或读取位移解
		stringstream worus;
		if(i<10)	worus << "./Results/Data_Backup_000" << i << ".dat";
		else if (i<100)	worus << "./Results/Data_Backup_00" << i << ".dat";
		else if (i<1000)	worus << "./Results/Data_Backup_0" << i << ".dat";
		else	worus << "./Results/Data_Backup_" << i << ".dat";
		if(write_or_read_data(worus.str(), wr_mod)==0) return 0;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//输出位移场云图
		cout << endl;
		cout << "    Exporting displacement field..." << endl;
		hout << endl;
		hout << "    Exporting displacement field..." << endl;
		stringstream EDC;
		if(i<10)	EDC << "./Results/Displacement_Field_Contour_000" << i << ".dat";
		else if (i<100)	EDC << "./Results/Displacement_Field_Contour_00" << i << ".dat";
		else if (i<1000)	EDC << "./Results/Displacement_Field_Contour_0" << i << ".dat";
		else EDC << "./Results/Displacement_Field_Contour_" << i << ".dat";
		Export_Displacement_Contour(EDC.str(), nodes, elements);

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//输出边界面上位移场云图（内点处位移为0）
		cout << "    Exporting boundarydisplacement field..." << endl;
		hout << "    Exporting boundarydisplacement field..." << endl;
		stringstream EBDC;
		if(i<10) EBDC << "./Results/Boundary_Displacement_Field_Contour_000" << i << ".dat";
		else if (i<100) EBDC << "./Results/Boundary_Displacement_Field_Contour_00" << i << ".dat";
		else if (i<1000) EBDC << "./Results/Boundary_Displacement_Field_Contour_0" << i << ".dat";
		else EBDC <<"./Results/Boundary_Displacement_Field_Contour_"<< i <<".dat";
//		Export_Boundary_Displacement_Contour(EBDC.str(), nodes, elements, 0);

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//输出网格变形图
		cout << "    Exporting deformed mesh..." << endl;
		hout << "    Exporting deformed mesh..." << endl;
		stringstream EDM;
		if(i<10) EDM <<"./Results/Deformed_Mesh_000"<< i <<".dat";
		else if (i<100) EDM <<"./Results/Deformed_Mesh_00"<< i <<".dat";
		else if (i<1000) EDM <<"./Results/Deformed_Mesh_0"<< i <<".dat";
		else EDM <<"./Results/Deformed_Mesh_"<< i <<".dat";
//		Export_Deformed_Mesh(EDM.str(), nodes, elements);

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//计算单元形心处应变向量
		cout << "    Calculating strain..." << endl;
		hout << "    Calculating strain..." << endl;
		Calculate_Elements_Strain(nodes, elements);

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//Estimate the average stress and strain at the boundary layer elements
		Estimate_boundary_stress_strain(nodes, elements, mats, i, Init->geom_rve, Init->iter, Init->displace, ave_strain_xx, ave_stress_xx);

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//计算节点处应变向量
		Calculate_Nodes_Strain(nodes, elements);

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//输出应变场云图
		cout << "    Exporting strain field..." << endl;
		hout << "    Exporting strain field..." << endl;
		stringstream ESM;
		if(i<10) ESM << "./Results/Strain_Field_Contour_000" << i << ".dat";
		else if (i<100) ESM << "./Results/Strain_Field_Contour_00" << i << ".dat";
		else if (i<1000) ESM << "./Results/Strain_Field_Contour_0" << i << ".dat";
		else ESM << "./Results/Strain_Field_Contour_" << i << ".dat";
		Export_Strain_Mesh(ESM.str(), nodes, elements);

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//坐标变换并输出应变场云图
		//Export_Strain_After_Coordinates_Transformation(nodes, elements);
	}
	
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Export stress strain curve
	cout<<"    Exporting stess-strain curve..."<<endl;
	hout<<"    Exporting stess-strain curve..."<<endl;
	Export_Stress_Strain_Curve("./Results/Stress_strain_curve.dat", ave_strain_xx, ave_stress_xx);

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Export total dissipative energy
	cout<<"    Exporting dissipative energy..."<<endl;
	hout<<"    Exporting dissipative energy..."<<endl;
	Export_Total_Dissipative_Energy("./Results/Dissipative_energy_plot.dat", Dissipation,  wr_mod);
	
	return 1;
}
//---------------------------------------------------------------------------
//Export total dissipative energy
int Postprocessor::Export_Total_Dissipative_Energy(const string &output_file_name, const vector<vector<double> > &Dissipation, const string &wr_mod)
{
	if(wr_mod=="Write")
	{
		total_ener.assign((int)Dissipation.size(), 0.0);
		for(int i=0; i<(int)Dissipation.size();i++)
			for(int j=0; j<(int)Dissipation[i].size(); j++) 
				total_ener[i] += Dissipation[i][j];
	}

	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Dissipative_energy_plot" << endl;
	otec << "VARIABLES = Steps, Energy_dissip" << endl;
	otec << "ZONE I=" << (int)total_ener.size() << ", F=POINT" << endl;

	for(int i=0; i<(int)total_ener.size(); i++) otec << setprecision(15) << i << " " << total_ener[i] << endl;

	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//Export 2D dissipatve energy in tecplot
int Postprocessor::Export_Dissipative_Energy_2D(const string &output_file_name, const vector<double> &damage_ener, const vector<double> &break_ener)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Dissipative_energy_plot" << endl;
	otec << "VARIABLES = Steps, dissipative_energy" << endl;
	
	//Total dissipative energy
	otec << "ZONE I=" << (int)damage_ener.size() << ", F=POINT" << endl;
	for(int i=0; i<(int)damage_ener.size(); i++) 
	otec << setprecision(15) << i << " " <<  damage_ener[i]+break_ener[i] << endl;

	//Damaged energy
	otec << "ZONE I=" << (int)damage_ener.size() << ", F=POINT" << endl;
	for(int i=0; i<(int)damage_ener.size(); i++) 
	otec << setprecision(15) << i << " " << damage_ener[i] << endl;

	//Bonds broken energy
	otec << "ZONE I=" << (int)break_ener.size() << ", F=POINT" << endl;
	for(int i=0; i<(int)break_ener.size(); i++) 
	otec << setprecision(15) << i << " " << break_ener[i] << endl;

	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//Export stress strain curve
int Postprocessor::Export_Stress_Strain_Curve(const string &output_file_name, const vector<double> &ave_strain_xx, const vector<double> &ave_stress_xx)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Stress_strain_curve" << endl;
	otec << "VARIABLES = Strain_Steps, Stress" << endl;
	otec << "ZONE I=" << (int)ave_strain_xx.size() << ", F=POINT" << endl;

//	for(int i=0; i<(int)ave_strain_xx.size(); i++) otec << setprecision(15) << ave_strain_xx[i] << " " << ave_stress_xx[i] << endl;
	for(int i=0; i<(int)ave_strain_xx.size(); i++) otec << i << setprecision(15) << " " << ave_stress_xx[i] << endl;

	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//Export boundary force displacement curve for 2D three point bending
int Postprocessor::Export_Force_Displacement_Curve(const string &output_file_name, const vector<double> &TPB_force_yy, double &delta_disp)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Force_displacement_curve" << endl;
	otec << "VARIABLES =Displaement, Force" << endl;
	otec << "ZONE I=" << (int)TPB_force_yy.size() << ", F=POINT" << endl;

//	for(int i=0; i<(int)TPB_force_yy.size(); i++) otec << i*delta_disp << setprecision(15) << " " << TPB_force_yy[i] << endl;
	for(int i=0; i<(int)TPB_force_yy.size(); i++) otec << i*delta_disp << setprecision(15) << " " << -TPB_force_yy[i]/1000 << endl;

	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//输出或读取位移解
int Postprocessor::write_or_read_data(const string &output_file_name, const string &wr_mod)
{
	//输出数据
	if(wr_mod=="Write")
	{
		ofstream odata(output_file_name.c_str());
		odata << "u_data:" << endl;
		odata << (int)u.size() << endl;
		for(int i=0; i<(int)u.size(); i++)
		{
			odata << setprecision(18) << setw(24) << setiosflags(ios::left) << u[i] << endl;
		}
		//---------------------------------------------------------------------------
		odata << "element_data:" << endl;
		odata << (int)elements.size() << endl;
		for(int i=0; i<(int)elements.size(); i++)
		{
			odata << elements[i].type << endl; 
			odata << elements[i].mat << endl;
			odata << elements[i].flag << endl;
			odata << elements[i].DG_elem << endl;
			odata << elements[i].dam << endl;
			odata << elements[i].Damage_energy << endl;
			odata << elements[i].Dissipative_energy << endl;
			odata << (int)elements[i].nodes_id.size() << endl;
			for(int j=0; j<(int)elements[i].nodes_id.size(); j++) odata << elements[i].nodes_id[j] << endl;
		}
		//---------------------------------------------------------------------------
		odata << "node_data:" << endl;
		odata << (int)nodes.size() << endl;
		for(int i=0; i<(int)nodes.size(); i++)
		{
			odata << nodes[i].type << endl;
			odata << nodes[i].DG_node << endl;
			odata << nodes[i].x << endl;
			odata << nodes[i].y << endl;
			odata << nodes[i].z << endl;
		}
		odata.close();
	}
	else if(wr_mod=="Read") //读入数据
	{
		//---------------------------------------------------------------------------
		u.clear();
		elements.clear();
		nodes.clear();

		//---------------------------------------------------------------------------
		ifstream idata(output_file_name.c_str());
		string u_comments;
		idata >> u_comments;
		int u_size;
		idata >> u_size;
		for(int i=0; i<u_size; i++)
		{
			double num;
			idata >> num;
			u.push_back(num);
		}
		//---------------------------------------------------------------------------
		string element_comments;
		idata >> element_comments;
		int eles_size;
		idata >> eles_size;
		for(int i=0; i<eles_size; i++)
		{
			Element eles;
			idata >> eles.type;
			idata >> eles.mat;
			idata >> eles.flag;
			idata >> eles.DG_elem;
			idata >> eles.dam;
			idata >> eles.Damage_energy;
			idata >> eles.Dissipative_energy;
			int ele_nod_size;
			idata >> ele_nod_size;
			int nodid = 0;
			for(int j=0; j<ele_nod_size; j++) 
			{
				idata >> nodid;
				eles.nodes_id.push_back(nodid);
			}
			elements.push_back(eles);
		}
		//---------------------------------------------------------------------------
		string node_comments;
		idata >> node_comments;
		int nods_size;
		idata >> nods_size;
		for(int i=0; i<nods_size; i++)
		{
			Node nods;
			idata >> nods.type;
			idata >> nods.DG_node;
			idata >> nods.x;
			idata >> nods.y;
			idata >> nods.z;
			nodes.push_back(nods);
		}
		idata.close();
	}
	else
	{
		hout << "注意！读写指令即不是读也不是写指令，请检查(write_or_read_data)！" << endl;
		return 0;
	}
	return 1;
}
//---------------------------------------------------------------------------
//计算单元形心处应变向量
int Postprocessor::Calculate_Elements_Strain(const vector<Node> &nodes, const vector<Element> &elements)
{
	//---------------------------------------------------------------------------
	//计算单元形心点处的应变向量
	vector<double> temp_strain(6,0.0);
	ele_strain.assign((int)elements.size(), temp_strain);
	for(int i=0; i<(int)elements.size(); i++)
	{
		//计算B矩阵；
		//--------------------------------------------
		//形函数N对正六面体中心点(0,0,0)的偏导矩阵	
		double diff[3][8];
		diff[0][0]=-0.125;
		diff[0][1]=-diff[0][0];                         
		diff[0][2]=0.125;
		diff[0][3]=-diff[0][2];
		diff[0][4]=-0.125;
		diff[0][5]=-diff[0][4];
		diff[0][6]=0.125;
		diff[0][7]=-diff[0][6];

		diff[1][0]=-0.125;
		diff[1][1]=-0.125;
		diff[1][2]=-diff[1][1];
		diff[1][3]=-diff[1][0];
		diff[1][4]=-0.125;
		diff[1][5]=-0.125;
		diff[1][6]=-diff[1][5];
		diff[1][7]=-diff[1][4];

		diff[2][0]=-0.125;
		diff[2][1]=-0.125;
		diff[2][2]=-0.125;
		diff[2][3]=-0.125;
		diff[2][4]=-diff[2][0];
		diff[2][5]=-diff[2][1];
		diff[2][6]=-diff[2][2];
		diff[2][7]=-diff[2][3];		
		//--------------------------------------------------
		//单元节点坐标矩阵
		double elenode[8][3];
		for(int j=0; j<8; j++)
		{
			elenode[j][0]=nodes[elements[i].nodes_id[j]].x;
			elenode[j][1]=nodes[elements[i].nodes_id[j]].y;
			elenode[j][2]=nodes[elements[i].nodes_id[j]].z;
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
		for(int j=0; j<3; j++)
			for(int k=0; k<8; k++)
			{
				diffxy[j][k]=0;
				for(int m=0; m<3; m++)
					diffxy[j][k] += Jinverse[j][m]*diff[m][k];
			}
		//--------------------------------------------------------
		//求出B矩阵
		double B[6][24];
		for(int j=0; j<6; j++)
			for(int k=0; k<24; k++)
				B[j][k]=0;
		for(int j=0; j<8; j++)
		{
			B[0][j*3+0]=diffxy[0][j];
			B[1][j*3+1]=diffxy[1][j];
			B[2][j*3+2]=diffxy[2][j];
			B[3][j*3+0]=diffxy[1][j];
			B[3][j*3+1]=diffxy[0][j];
			B[4][j*3+1]=diffxy[2][j];
			B[4][j*3+2]=diffxy[1][j];
			B[5][j*3+0]=diffxy[2][j];
			B[5][j*3+2]=diffxy[0][j];
		}		
		//--------------------------------------------
		//计算单元形心点应变值
		vector<double> ele_u(24,0.0);
		for(int j=0; j<8; j++)
		{
			ele_u[3*j] = u[3*elements[i].nodes_id[j]];
			ele_u[3*j+1] = u[3*elements[i].nodes_id[j]+1];
			ele_u[3*j+2] = u[3*elements[i].nodes_id[j]+2];
		}
		for(int j=0; j<6; j++)
		{
			for(int k=0; k<24; k++)
			{
				ele_strain[i][j] = ele_strain[i][j] + B[j][k]*ele_u[k];
			}
		}
	}
		
	return 1;
}
//---------------------------------------------------------------------------
//计算节点处应变向量
int Postprocessor::Calculate_Nodes_Strain(vector<Node> &nodes, const vector<Element> &elements)
{
	//计算节点处的应变向量
	vector<double> temp_strain(6,0.0);
	nod_strain.assign((int)nodes.size(), temp_strain);
	nod_damage.assign((int)nodes.size(),0.0);
	nod_dissip.assign((int)nodes.size(),0.0);

	//---------------------------------------------------------------------------
	//统计节点的相关单元信息
	for(int i=0; i<(int)nodes.size(); i++)  //相关单元清零
	{
		nodes[i].relative_eles.clear();
	}
	for(int i=0; i<(int)elements.size(); i++)
	{
		int node_size = (int)elements[i].nodes_id.size();
		for(int j=0; j<node_size; j++)
		{
			nodes[elements[i].nodes_id[j]].relative_eles.push_back(i);
		}
	}

	//---------------------------------------------------------------------------
	//Calculate the area of every element
	vector<double> ele_area(elements.size(), 0.0);
	Mesher Mesh;
	Node elenod[4];
	for(int i=0; i<(int)elements.size(); i++)
	{
		for(int j=0; j<4; j++) elenod[j] = nodes[elements[i].nodes_id[j]];
		ele_area[i] = Mesh.Calculate_quadri_area(elenod);
	}

	//---------------------------------------------------------------------------
	//计算节点处的应变向量
	for(int i=0; i<(int)nodes.size(); i++)
	{
		vector<int> same_nods;
		Find_same_DG_nodes(i, nodes, same_nods);
		int ncount = 0;
		for(int j=0; j<(int)same_nods.size(); j++)
		{
			int nrev_size = (int)nodes[same_nods[j]].relative_eles.size();
			ncount += nrev_size;
			for(int k=0; k<nrev_size; k++)
			{
				for(int m=0; m<6; m++)
				{
					nod_strain[i][m] += ele_strain[nodes[same_nods[j]].relative_eles[k]][m];
				}
				int Nsnjrek = nodes[same_nods[j]].relative_eles[k];
				nod_damage[i] += elements[Nsnjrek].Damage_energy/ele_area[Nsnjrek];
				nod_dissip[i] += elements[Nsnjrek].Dissipative_energy/ele_area[Nsnjrek];
			}
		}
		for(int j=0; j<6; j++)	nod_strain[i][j] = nod_strain[i][j]/ncount;
		nod_damage[i] = nod_damage[i]/ncount;
		nod_dissip[i] = nod_dissip[i]/ncount;
	}

	for(int i=0; i<(int)nodes.size(); i++) //相关单元清零
	{
		nodes[i].relative_eles.clear();
	}

	return 1;
}
//---------------------------------------------------------------------------
//计算节点处应变向量、应力向量和应变能密度
int Postprocessor::Calculate_Nodes_Strain_Stress_Energy(string mod, double d_crit, vector<Node> &nodes, const vector<Element> &elements)
{
	//计算节点处的应变向量
	vector<double> temp(6, 0.0);
	nod_strain.assign((int)nodes.size(), temp);
	nod_stress.assign((int)nodes.size(), temp);
	nod_local_dam.assign((int)nodes.size(), 0.0);
	nod_nonlocal_dam.assign((int)nodes.size(), 0.0);
	nod_damage.assign((int)nodes.size(), 0.0);
	nod_dissip.assign((int)nodes.size(), 0.0);
	nod_dissip_ratio.assign((int)nodes.size(), 0.0);
	nod_energy.assign((int)nodes.size(), 0.0);

	//---------------------------------------------------------------------------
	//统计节点的相关单元信息
	for (int i = 0; i<(int)nodes.size(); i++)  //相关单元清零
	{
		nodes[i].relative_eles.clear();
	}
	for (int i = 0; i<(int)elements.size(); i++)
	{
		int node_size = (int)elements[i].nodes_id.size();
		for (int j = 0; j<node_size; j++)
		{
			nodes[elements[i].nodes_id[j]].relative_eles.push_back(i);
		}
	}

	//---------------------------------------------------------------------------
	//Calculate the area of every element
	vector<double> ele_area(elements.size(), 0.0);
	Mesher Mesh;
	Node elenod[4];
	for (int i = 0; i<(int)elements.size(); i++)
	{
		for (int j = 0; j<4; j++) elenod[j] = nodes[elements[i].nodes_id[j]];
		ele_area[i] = Mesh.Calculate_quadri_area(elenod);
	}

	//---------------------------------------------------------------------------
	//计算节点处的应变向量
	for (int i = 0; i<(int)nodes.size(); i++)
	{
		vector<int> same_nods;
		Find_same_DG_nodes(i, nodes, same_nods);
		int ncount = 0;
		for (int j = 0; j<(int)same_nods.size(); j++)
		{
			int nrev_size = (int)nodes[same_nods[j]].relative_eles.size();
			ncount += nrev_size;
			for (int k = 0; k<nrev_size; k++)
			{
				for (int m = 0; m<6; m++)
				{
					nod_strain[i][m] += ele_strain[nodes[same_nods[j]].relative_eles[k]][m];
					nod_stress[i][m] += ele_stress[nodes[same_nods[j]].relative_eles[k]][m];
				}
				int Nsnjrek = nodes[same_nods[j]].relative_eles[k];
				nod_damage[i] += elements[Nsnjrek].Damage_energy/ele_area[Nsnjrek];
				nod_dissip[i] += elements[Nsnjrek].Dissipative_energy/ele_area[Nsnjrek];
				if (mod == "Hybrid")
				{
					nod_local_dam[i] += elements[Nsnjrek].dam;
					nod_nonlocal_dam[i] += (1.0-d_crit)*elements[Nsnjrek].Dissipative_energy/total_break_energy[Nsnjrek];
					nod_dissip_ratio[i] += elements[Nsnjrek].Dissipative_energy/total_break_energy[Nsnjrek];
				}
				else
				{
					nod_local_dam[i] += elements[Nsnjrek].dam;
					nod_dissip_ratio[i] += 0.0;
				}
				nod_energy[i] += ele_energy[Nsnjrek];
			}
		}
		for (int j = 0; j < 6; j++)
		{
			nod_strain[i][j] = nod_strain[i][j] / ncount;
			nod_stress[i][j] = nod_stress[i][j] / ncount;
		}
		nod_local_dam[i] = nod_local_dam[i] / ncount;
		nod_nonlocal_dam[i] = nod_nonlocal_dam[i] / ncount;
		nod_damage[i] = nod_damage[i] / ncount;
		nod_dissip[i] = nod_dissip[i] / ncount;
		nod_dissip_ratio[i] = nod_dissip_ratio[i] / ncount;
		nod_energy[i] = nod_energy[i] / ncount;
	}

	for (int i = 0; i<(int)nodes.size(); i++) //相关单元清零
	{
		nodes[i].relative_eles.clear();
	}

	return 1;
}
//---------------------------------------------------------------------------
//输出Tecplot可视化位移场云图
int Postprocessor::Export_Displacement_Contour(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Displacement_Field_Contour" << endl;
	otec << "VARIABLES = X, Y, Z";
	otec << ", Ux, Uy, Uz" << endl;
	
	otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size() << ", F=FEPOINT, ET=BRICK" << endl;
	for(int i=0; i<(int)nodes.size(); i++)
	{
		otec << nodes[i].x << "  " << nodes[i].y << "  " << nodes[i].z;
		otec << "  " << u[3*i] << "  " << u[3*i+1] << "  " << u[3*i+2];
		otec << endl;
	}
	otec << endl;
	for (int i=0; i < (int)elements.size(); i++)
	{
		otec	 << elements[i].nodes_id[0]+1 << "  " << elements[i].nodes_id[1]+1 << "  " 
				 << elements[i].nodes_id[2]+1 << "  " << elements[i].nodes_id[3]+1 << "  " 
				 << elements[i].nodes_id[4]+1 << "  " << elements[i].nodes_id[5]+1 << "  " 
				 << elements[i].nodes_id[6]+1 << "  " << elements[i].nodes_id[7]+1 << endl;
	}
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//输出Tecplot可视化边界面位移场云图
int Postprocessor::Export_Boundary_Displacement_Contour(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements, const int key)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Boundary_Displacement_Field_Contour" << endl;
	otec << "VARIABLES = X, Y, Z";
	otec << ", Ux, Uy, Uz" << endl;
	
	otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size() << ", F=FEPOINT, ET=BRICK" << endl;
	for(int i=0; i<(int)nodes.size(); i++)
	{
		otec << nodes[i].x << "  " << nodes[i].y << "  " << nodes[i].z;
		if(nodes[i].type==0)	otec << "  " << 0 << "  " << 0 << "  " << 0; //内点
		else if(key==1) //纯剪切
		{
			if(nodes[i].z==1||nodes[i].x==1||nodes[i].y==0) otec << "  " << u[3*i] << "  " << u[3*i+1] << "  " << u[3*i+2]; //输出三个面
			else otec << "  " << 0 << "  " << 0 << "  " << 0;
		}
		else if(key==0||key==2) //拉伸或裂纹例子 
		{	
			if(nodes[i].z==1) otec << "  " << u[3*i] << "  " << u[3*i+1] << "  " << u[3*i+2];  //输出一个面
			else otec << "  " << 0 << "  " << 0 << "  " << 0; 
		}
		else { hout<<"输入的key没有找到对应操作！（Export_Boundary_Displacement_Contour）"; }

		otec << endl;
	}
	otec << endl;
	for (int i=0; i < (int)elements.size(); i++)
	{
		otec	 << elements[i].nodes_id[0]+1 << "  " << elements[i].nodes_id[1]+1 << "  " 
				 << elements[i].nodes_id[2]+1 << "  " << elements[i].nodes_id[3]+1 << "  " 
				 << elements[i].nodes_id[4]+1 << "  " << elements[i].nodes_id[5]+1 << "  " 
				 << elements[i].nodes_id[6]+1 << "  " << elements[i].nodes_id[7]+1 << endl;
	}
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//输出网格变形图
int Postprocessor::Export_Deformed_Mesh(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Deformed_Mesh" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size()  << ", F=FEPOINT, ET=BRICK" << endl;
	for(int i=0; i<(int)nodes.size(); i++)
	{
		otec << setprecision(15) << nodes[i].x+u[3*i] << "  " << nodes[i].y+u[3*i+1] << "  " << nodes[i].z+u[3*i+2] << endl;
	}
	for (int i=0; i<(int)elements.size(); i++)
	{
		otec	 << elements[i].nodes_id[0]+1 << "  " << elements[i].nodes_id[1]+1 << "  " 
					<< elements[i].nodes_id[2]+1 << "  " << elements[i].nodes_id[3]+1 << "  " 
					<< elements[i].nodes_id[4]+1 << "  " << elements[i].nodes_id[5]+1 << "  " 
					<< elements[i].nodes_id[6]+1 << "  " << elements[i].nodes_id[7]+1 << endl;
	}
	otec << endl;
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//Export deformed mesh for testing
int Postprocessor::Testing_export_deformed_mesh(const string &output_file_name)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE =" << output_file_name << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size()  << ", F=FEPOINT, ET=BRICK" << endl;
	for(int i=0; i<(int)nodes.size(); i++)
	{
		otec << setprecision(15) << nodes[i].x+u[3*i] << "  " << nodes[i].y+u[3*i+1] << "  " << nodes[i].z+u[3*i+2] << endl;
	}
	for (int i=0; i<(int)elements.size(); i++)
	{
		otec	 << elements[i].nodes_id[0]+1 << "  " << elements[i].nodes_id[1]+1 << "  " 
					<< elements[i].nodes_id[2]+1 << "  " << elements[i].nodes_id[3]+1 << "  " 
					<< elements[i].nodes_id[4]+1 << "  " << elements[i].nodes_id[5]+1 << "  " 
					<< elements[i].nodes_id[6]+1 << "  " << elements[i].nodes_id[7]+1 << endl;
	}
	otec << endl;
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//输出应变场云图
int Postprocessor::Export_Strain_Mesh(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Strain_Field_Contour" << endl;
	otec << "VARIABLES = X, Y, Z";
	otec << ", dissip_en";
	otec << ", <greek>e</greek><sub>xx</sub>, <greek>e</greek><sub>yy</sub>, <greek>e</greek><sub>zz</sub>";
	otec << ", <greek>e</greek><sub>xy</sub>, <greek>e</greek><sub>yz</sub>, <greek>e</greek><sub>xz</sub>";
	otec << ", <greek>e</greek><sub>xx</sub>error(%), <greek>e</greek><sub>yy</sub>error(%), <greek>e</greek><sub>zz</sub>error(%)";
	otec << ", <greek>e</greek><sub>xy</sub>error(%), <greek>e</greek><sub>yz</sub>error(%), <greek>e</greek><sub>xz</sub>error(%)" << endl;

	//------------------------------------------------------------------------------------------------------------------------------------
	//平均应变
	double strain_average[6] = {0};
	for(int i=0; i<6; i++)
	{
		for(int j=0; j<(int)nodes.size(); j++)
			strain_average[i] +=  nod_strain[j][i];
		strain_average[i] = strain_average[i]/(int)nodes.size();
	}
	//------------------------------------------------------------------------------------------------------------------------------------
	//输出
	otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size()  << ", F=FEPOINT, ET=BRICK" << endl;
	for(int i=0; i<(int)nodes.size(); i++)
	{
		otec << setprecision(15) << nodes[i].x+u[3*i] << "  " << nodes[i].y+u[3*i+1] << "  " << nodes[i].z+u[3*i+2];
		otec << setprecision(15) << "  " << nod_dissip[i];
		otec << setprecision(15) << "  " << nod_strain[i][0] << "  " << nod_strain[i][1] << "  " << nod_strain[i][2];
		otec << setprecision(15) << "  " << nod_strain[i][3] << "  " << nod_strain[i][4] << "  " << nod_strain[i][5];
		otec << setprecision(15) << "  " << 100*(nod_strain[i][0]-strain_average[0])/strain_average[0] << "  " << 100*(nod_strain[i][1]-strain_average[1])/strain_average[1];
		otec << setprecision(15) << "  " << 100*(nod_strain[i][2]-strain_average[2])/strain_average[2] << "  " << 100*(nod_strain[i][3]-strain_average[3])/strain_average[3];
		otec << setprecision(15) << "  " << 100*(nod_strain[i][4]-strain_average[4])/strain_average[4] << "  " << 100*(nod_strain[i][5]-strain_average[5])/strain_average[5] << endl;
	}
	for (int i=0; i<(int)elements.size(); i++)
	{
		otec	 << elements[i].nodes_id[0]+1 << "  " << elements[i].nodes_id[1]+1 << "  " 
					<< elements[i].nodes_id[2]+1 << "  " << elements[i].nodes_id[3]+1 << "  " 
					<< elements[i].nodes_id[4]+1 << "  " << elements[i].nodes_id[5]+1 << "  " 
					<< elements[i].nodes_id[6]+1 << "  " << elements[i].nodes_id[7]+1 << endl;
	}
	otec << endl;

	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//坐标变换并输出应变场云图
int Postprocessor::Export_Strain_After_Coordinates_Transformation(const vector<Node> &nodes, const vector<Element> &elements)const
{
	//------------------------------------------------------------------------------------------------------------------------------------
	//平均应变
	double strain_average[6] = {0};
	for(int i=0; i<6; i++)
	{
		for(int j=0; j<(int)nodes.size(); j++)
			strain_average[i] +=  nod_strain[j][i];
		strain_average[i] = strain_average[i]/(int)nodes.size();
	}

	//------------------------------------------------------------------------------------------------------------------------------------
	//求变换后的坐标系
	const double original_coord[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };  //原点在(0, 0, 0)
	const double sval = 0.25;  //shear value
	const double smat[3][3] = { {1, sval, sval}, {sval, 1, sval}, {sval, sval, 1} };  //shear_deformed_matrix

	double transf_coord[3][3] = {{0}, {0}, {0}};  //原点在(0, 0, 0)
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			for(int k=0; k<3; k++)
				transf_coord[i][j] += smat[i][k]*original_coord[j][k];

	Point_3D coord[3];
	for(int i=0; i<3; i++) 
	{
		coord[i].x = transf_coord[i][0];
		coord[i].y = transf_coord[i][1];
		coord[i].z = transf_coord[i][2];
	}

	Point_3D normal[9];
	normal[0] = Vector_Products(coord[0], coord[1]);
	normal[1] = Vector_Products(coord[0], coord[2]);
	normal[2] = Vector_Products(coord[1], coord[2]);
	normal[3] = Vector_Products(coord[0]-coord[1], coord[2]);
	normal[4] = Vector_Products(coord[0]-coord[2], coord[1]);
	normal[5] = Vector_Products(coord[1]-coord[2], coord[0]);
	normal[6] = Vector_Products(coord[0]+coord[1], coord[2]);
	normal[7] = Vector_Products(coord[0]+coord[2], coord[1]);
	normal[8] = Vector_Products(coord[1]+coord[2], coord[0]);

	//---------------------------------------------------------------------------
	//输出图形
	for(int i=0; i<9; i++)
	{
		MathMatrix transf(3,3);
		Coordinates_Transformation_Matrix(normal[i], transf);  //法向向量旋转至Z轴变换矩阵
	
		string str[9] = {"XY", "XZ", "YZ", "XmYZ", "XmZY", "YmZX", "XpYZ", "XpZY", "YpZX"};
		string output_file_name = "Strain_Transf_Coord＿" + str[i] +".dat";

		ofstream otec(output_file_name.c_str());
		otec << "TITLE =" << output_file_name << endl;
		otec << "VARIABLES = X, Y, Z";
		otec << ", <greek>e</greek><sub>xx</sub>, <greek>e</greek><sub>yy</sub>, <greek>e</greek><sub>zz</sub>";
		otec << ", <greek>e</greek><sub>xy</sub>, <greek>e</greek><sub>yz</sub>, <greek>e</greek><sub>xz</sub>";
		otec << ", <greek>e</greek><sub>xx</sub>error(%), <greek>e</greek><sub>yy</sub>error(%), <greek>e</greek><sub>zz</sub>error(%)";
		otec << ", <greek>e</greek><sub>xy</sub>error(%), <greek>e</greek><sub>yz</sub>error(%), <greek>e</greek><sub>xz</sub>error(%)" << endl;

		//记录8个角点的位置以及应变值
		double eit_pois[8][9];
		for(int j=0; j<3; j++) eit_pois[0][j] = 0;
		for(int j=0; j<3; j++) eit_pois[1][j] = transf_coord[0][j];
		for(int j=0; j<3; j++) eit_pois[2][j] = transf_coord[0][j] + transf_coord[1][j];
		for(int j=0; j<3; j++) eit_pois[3][j] = transf_coord[1][j];
		for(int j=0; j<3; j++) eit_pois[4][j] = transf_coord[2][j];
		for(int j=0; j<3; j++) eit_pois[5][j] = transf_coord[0][j] + transf_coord[2][j];
		for(int j=0; j<3; j++) eit_pois[6][j] = transf_coord[0][j] + transf_coord[1][j] + transf_coord[2][j];
		for(int j=0; j<3; j++) eit_pois[7][j] = transf_coord[1][j] + transf_coord[2][j];
		int mark_pois[8] = {0};

		//------------------------------------------------------------------------------------------------------------------------------------
		//输出
		otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size()  << ", F=FEPOINT, ET=BRICK" << endl;
		for(int j=0; j<(int)nodes.size(); j++)
		{
			const double pois[3] = { nodes[j].x+u[3*j], nodes[j].y+u[3*j+1], nodes[j].z+u[3*j+2] }; 
			
			for(int k=0; k<8; k++)  //记录8个角点的应变值
			{
				if(mark_pois[k]==0&&fabs(eit_pois[k][0]-pois[0])<Zero&&fabs(eit_pois[k][1]-pois[1])<Zero&&fabs(eit_pois[k][2]-pois[2])<Zero)
				{
					for(int m=0; m<6; m++) eit_pois[k][m+3] = nod_strain[j][m];
					mark_pois[k] = 1;
					break;
				}
			}

			double new_pois[3] = { 0 };
			for(int k=0; k<3; k++)
				for(int m=0; m<3; m++)
					new_pois[k] += transf.element[k][m]*pois[m];

			otec << setprecision(15) << new_pois[0] << "  " << new_pois[1] << "  " << new_pois[2];
			otec << setprecision(15) << "  " << nod_strain[j][0] << "  " << nod_strain[j][1] << "  " << nod_strain[j][2];
			otec << setprecision(15) << "  " << nod_strain[j][3] << "  " << nod_strain[j][4] << "  " << nod_strain[j][5];
			otec << setprecision(15) << "  " << 100*(nod_strain[j][0]-strain_average[0])/strain_average[0] << "  " << 100*(nod_strain[j][1]-strain_average[1])/strain_average[1];
			otec << setprecision(15) << "  " << 100*(nod_strain[j][2]-strain_average[2])/strain_average[2] << "  " << 100*(nod_strain[j][3]-strain_average[3])/strain_average[3];
			otec << setprecision(15) << "  " << 100*(nod_strain[j][4]-strain_average[4])/strain_average[4] << "  " << 100*(nod_strain[j][5]-strain_average[5])/strain_average[5] << endl;
		}
		for(int j=0; j<(int)elements.size(); j++)
		{
			otec	 << elements[j].nodes_id[0]+1 << "  " << elements[j].nodes_id[1]+1 << "  " 
						<< elements[j].nodes_id[2]+1 << "  " << elements[j].nodes_id[3]+1 << "  " 
						<< elements[j].nodes_id[4]+1 << "  " << elements[j].nodes_id[5]+1 << "  " 
						<< elements[j].nodes_id[6]+1 << "  " << elements[j].nodes_id[7]+1 << endl;
		}
		otec << endl;

		//------------------------------------------------------------------------------------------------------------------------------------
		for(int j=0; j<8; j++)
		{
			if(mark_pois[j]!=1)
			{
				hout << "寻找8个角点的应变值时出错！" << endl;
				return 0;
			}
		}

		//------------------------------------------------------------------------------------------------------------------------------------
		//输出
		otec << "ZONE N=" << 8 << ", E=" << 1  << ", F=FEPOINT, ET=BRICK" << endl;
		for(int j=0; j<8; j++)
		{
			double new_pois[3] = { 0 };
			for(int k=0; k<3; k++)
				for(int m=0; m<3; m++)
					new_pois[k] += transf.element[k][m]*eit_pois[j][m];
			for(int k=0; k<3; k++) otec << setprecision(15) << new_pois[k] << "  ";
			for(int k=0; k<6; k++) otec << setprecision(15) << eit_pois[j][k+3] << "  ";
			for(int k=0; k<6; k++) otec << setprecision(15) << 100*(eit_pois[j][k+3]-strain_average[k])/strain_average[k] << "  ";
			otec << endl;
		}
		otec	 << "1 2 3 4 5 6 7 8" << endl;		

		otec.close();
	}

	return 1;
}
//---------------------------------------------------------------------------
//向量叉积
Point_3D Postprocessor::Vector_Products(const Point_3D &p1, const Point_3D &p2)const
{
	Point_3D temp;
	temp.x = p1.y*p2.z - p1.z*p2.y;
	temp.y = p1.z*p2.x - p1.x*p2.z;
	temp.z = p1.x*p2.y - p1.y*p2.x;

	return temp;
}
//---------------------------------------------------------------------------
//法向向量旋转至Z轴变换矩阵
void Postprocessor::Coordinates_Transformation_Matrix(const Point_3D &poi, MathMatrix &transf)const
{
	transf.element[0][0] = poi.x*poi.z/(sqrt(poi.x*poi.x+poi.y*poi.y+poi.z*poi.z)*sqrt(poi.x*poi.x+poi.y*poi.y));
	transf.element[0][1] = poi.y*poi.z/(sqrt(poi.x*poi.x+poi.y*poi.y+poi.z*poi.z)*sqrt(poi.x*poi.x+poi.y*poi.y));
	transf.element[0][2] = sqrt(poi.x*poi.x+poi.y*poi.y)/sqrt(poi.x*poi.x+poi.y*poi.y+poi.z*poi.z);
	transf.element[1][0] = -poi.y/sqrt(poi.x*poi.x+poi.y*poi.y);
	transf.element[1][1] = poi.x/sqrt(poi.x*poi.x+poi.y*poi.y);
	transf.element[1][2] = 0;
	transf.element[2][0] = poi.x/sqrt(poi.x*poi.x+poi.y*poi.y+poi.z*poi.z);
	transf.element[2][1] = poi.y/sqrt(poi.x*poi.x+poi.y*poi.y+poi.z*poi.z);
	transf.element[2][2] = poi.z/sqrt(poi.x*poi.x+poi.y*poi.y+poi.z*poi.z);
}
//---------------------------------------------------------------------------
//Immediately write result data
int Postprocessor::Immedi_write_data(const int &Unum, const string &wr_mod)
{
	cout << endl;
	cout << "Immediately output result data: loop  " << Unum << endl;
	hout << endl;
	hout << "Immediately output result data: loop  " << Unum << endl;
	hout << "-------------------------------------" << endl;

	if(wr_mod=="Write") 
	{
		ofstream onum("./Results/Data_Backup_Files_Number.dat");
		onum << Unum+1 << endl;
		onum.close();
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Output result data
	stringstream worus;
	if(Unum<10)	worus << "./Results/Data_Backup_000" << Unum << ".dat";
	else if (Unum<100)	worus << "./Results/Data_Backup_00" << Unum << ".dat";
	else if (Unum<1000)	worus << "./Results/Data_Backup_0" << Unum << ".dat";
	else	worus << "./Results/Data_Backup_" << Unum << ".dat";
	if(write_or_read_data(worus.str(), wr_mod)==0) return 0;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	cout << "The loop  " << Unum << " result data output successfully!" << endl;
	cout << "-------------------------------------" << endl << endl;
	hout << "The loop  " << Unum << " result data output successfully!" << endl;
	hout << "-------------------------------------" << endl << endl;

	return 1;
}
//---------------------------------------------------------------------------
//Immediately write result data with weight functions
int Postprocessor::Immedi_write_data(const int &Unum, const string &wr_mod, const struct Weight_func &weight_func)
{
	cout << endl;
	cout << "Immediately output result data: loop  " << Unum << endl;
	hout << endl;
	hout << "Immediately output result data: loop  " << Unum << endl;
	hout << "-------------------------------------" << endl;

	if(wr_mod=="Write") 
	{
		ofstream onum("./Results/Data_Backup_Files_Number.dat");
		onum << Unum+1 << endl;
		onum << weight_func.keywords << endl;
		onum << weight_func.mark << endl;
		onum << weight_func.num << endl;
		for(int i=0; i<weight_func.num; i++) 
		{
			onum << weight_func.shape[i] << endl;
			onum << weight_func.center[i].x << endl;
			onum << weight_func.center[i].y << endl;
			onum << weight_func.center[i].z << endl;
			onum << weight_func.r0[i] << endl;
			onum << weight_func.r1[i] << endl;
			onum << weight_func.ratio[i] << endl;
			onum << weight_func.func_order[i] << endl;
			onum << weight_func.func_constant[i] << endl;
		}
		onum.close();
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Output result data
	stringstream worus;
	if(Unum<10)	worus << "./Results/Data_Backup_000" << Unum << ".dat";
	else if (Unum<100)	worus << "./Results/Data_Backup_00" << Unum << ".dat";
	else if (Unum<1000)	worus << "./Results/Data_Backup_0" << Unum << ".dat";
	else	worus << "./Results/Data_Backup_" << Unum << ".dat";
	if(write_or_read_data(worus.str(), wr_mod)==0) return 0;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	cout << "The loop  " << Unum << " result data output successfully!" << endl;
	cout << "-------------------------------------" << endl << endl;
	hout << "The loop  " << Unum << " result data output successfully!" << endl;
	hout << "-------------------------------------" << endl << endl;

	return 1;
}
//---------------------------------------------------------------------------
//Estimate the average stress and strain at the boundary layer elements
void Postprocessor::Estimate_boundary_stress_strain(const vector<Node> &nodes, const vector<Element> &elements, const vector<MatPro> &mats, const int &inum, const struct Geom_RVE &geom_rve,
																							  const struct Iterative &iter, const struct Displace &displace, vector<double> &ave_strain_xx, vector<double> &ave_stress_xx)const
{
	//-----------------------------------------------------------------
	//Updata strain
	double total_strain = 2.0*displace.value[1][0]/geom_rve.len_x;
	ave_strain_xx[inum+1] = total_strain*(inum+1)/iter.ramp_para;

	//-----------------------------------------------------------------
	//Updata stress
	int ele_num = 0;
	double stress_xx = 0.0;
	//-----------------------------------------------------------------
	//Stiffness matrix
	double ele_elas[6][6];		
	for(int i=0; i<6; i++)
		for(int j=0; j<6; j++)
			ele_elas[i][j] = mats[1].elas_matrix[i][j];

	//-----------------------------------------------------------------
	for(int i=0; i<(int)elements.size(); i++)
		for(int j=0; j<(int)elements[i].nodes_id.size(); j++)
		{
			//-------------------------------------------------------------------
			if(nodes[elements[i].nodes_id[j]].x==geom_rve.origin.x+geom_rve.len_x)
			{
				ele_num ++;
				for(int k=0; k<6; k++) stress_xx += ele_elas[0][k]*ele_strain[i][k];
				break;
			}
		}
	
	//-----------------------------------------------------------------
	ave_stress_xx[inum+1] = 2.0*stress_xx/ele_num;
}
//---------------------------------------------------------------------------
//Estimate the average stress and strain at the boundary layer elements
int Postprocessor::Estimate_boundary_stress_strain(const vector<Node> &nodes, const vector<Element> &elements, const vector<MatPro> &mats, const int &inum, const struct Geom_RVE &geom_rve,
																				vector<double> &ave_strain_xx, vector<double> &ave_stress_xx)const
{
	double boudary_length = 0.0;
	double deform_xx=0.0, force_xx=0.0;
	//-----------------------------------------------------------------
	//查找边界节点
	for(int i=0; i<(int)elements.size(); i++)
	{
		vector<vector<double> > nod_strain;
		vector<int> nodj;
		for(int j=0; j<(int)elements[i].nodes_id.size(); j++)
		{
			if(nodes[elements[i].nodes_id[j]].x==geom_rve.origin.x+geom_rve.len_x)
			{
				//--------------------------------------------
				//仅适用于四边形
				if((int)elements[i].nodes_id.size()!=4) 
				{
					hout << "Attention: the number of nodes in the element which is on the boundary isn't equal to 4." << endl;
					return 0;
				}

				//计算这一节点在此单元中的应变
				//--------------------------------------------

				//计算B矩阵；
				//--------------------------------------------
				//形函数N对某节点的偏导矩阵
				double nx=0.0, ny=0.0;
				switch(j)
				{
					case 0:{ nx=-1.0; ny=-1.0; break; }
					case 1:{ nx=1.0; ny=-1.0; break; }
					case 2:{ nx=1.0; ny=1.0; break; }
					case 3:{ nx=-1.0; ny=1.0; break; }
					default: hout << "Attention: the number of nodes isn't from 0 to 3" << endl;
				}
				double diff[2][4];
				diff[0][0]=-0.25*(1.0-ny);
				diff[0][1]=-diff[0][0];
				diff[0][2]=0.25*(1.0+ny);
				diff[0][3]=-diff[0][2];

				diff[1][0]=-0.25*(1.0-nx);
				diff[1][1]=-0.25*(1.0+nx);
				diff[1][2]=-diff[1][1];
				diff[1][3]=-diff[1][0];

				//--------------------------------------------------
				//单元节点坐标矩阵
				double elenode[4][2];
				for(int k=0; k<4; k++)
				{
					elenode[k][0]=nodes[elements[i].nodes_id[k]].x;
					elenode[k][1]=nodes[elements[i].nodes_id[k]].y;
				}
				//--------------------------------------------------
				//J矩阵
				double Jmatrix[2][2];
				//以上两个矩阵的积
				for(int k=0; k<2; k++)
					for(int m=0; m<2; m++)
					{
						Jmatrix[k][m]=0;
						for(int p=0; p<4; p++)
						Jmatrix[k][m] += diff[k][p]*elenode[p][m];
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
				//求出N对x,y,的偏导
				double diffxy[2][4];
				for(int k=0; k<2; k++)
					for(int m=0; m<4; m++)
					{
						diffxy[k][m]=0;
						for(int p=0; p<2; p++)
							diffxy[k][m] += Jinverse[k][p]*diff[p][m];
					}

				//--------------------------------------------------------
				//求出B矩阵
				double B[6][12];
				for(int k=0; k<6; k++)
					for(int m=0; m<12; m++)
						B[k][m]=0;
				for(int k=0; k<4; k++)
				{
					B[0][k*3+0]=diffxy[0][k];
					B[1][k*3+1]=diffxy[1][k];
					B[3][k*3+0]=diffxy[1][k];
					B[3][k*3+1]=diffxy[0][k];
					B[4][k*3+2]=diffxy[1][k];
					B[5][k*3+2]=diffxy[0][k];
				}		
				//--------------------------------------------
				//计算单元形心点应变值
				vector<double> ele_u(12,0.0);
				for(int k=0; k<4; k++)
				{
					ele_u[3*k] = u[3*elements[i].nodes_id[k]];
					ele_u[3*k+1] = u[3*elements[i].nodes_id[k]+1];
				}
				vector<double> temp_strain(6,0.0);
				for(int k=0; k<6; k++)
					for(int m=0; m<12; m++)
						temp_strain[k] += B[k][m]*ele_u[m];

				nod_strain.push_back(temp_strain);
				nodj.push_back(j);
			}
		}
		if((int)nodj.size()>2)
		{
			hout << "Attention: more than 2 nodes of the element are in this boundary!" << endl;
			return 0;
		}
		else if((int)nodj.size()==2)
		{
			//计算边长
			Point_3D p0(nodes[elements[i].nodes_id[nodj[0]]].x, nodes[elements[i].nodes_id[nodj[0]]].y, nodes[elements[i].nodes_id[nodj[0]]].z);
			Point_3D p1(nodes[elements[i].nodes_id[nodj[1]]].x, nodes[elements[i].nodes_id[nodj[1]]].y, nodes[elements[i].nodes_id[nodj[1]]].z);
			double ele_bs =p0.distance_to(p1);

			//计算弹性矩阵
			double ele_elas[6][6] = { {0}, {0}, {0}, {0}, {0}, {0} };
			double damage_differ = 1.0 - elements[i].dam;
			if(damage_differ>Zero)
			{
				for(int j=0; j<6; j++)
					for(int k=0; k<6; k++)
						ele_elas[j][k] = mats[1].elas_matrix[j][k]*damage_differ;
			}

			//计算单元节点平均应变和应力
			double strain_xx=0.0, stress_xx=0.0;
			for(int j=0; j<(int)nod_strain.size(); j++)
			{
				strain_xx += nod_strain[j][0];
				for(int k=0; k<6; k++) stress_xx += ele_elas[0][k]*nod_strain[j][k];
			}
			deform_xx += 0.5*strain_xx*ele_bs;
			force_xx += 0.5*stress_xx*ele_bs;
			boudary_length += ele_bs; //总边长
		}
	}

	//-----------------------------------------------------------------
	ave_strain_xx[inum+1] = deform_xx/boudary_length;
	ave_stress_xx[inum+1] = force_xx/boudary_length;

	return 1;
}
//---------------------------------------------------------------------------
//Estimate the average force at the top center of beam under three point bending
int Postprocessor::Estimate_boundary_force_2D(const vector<Node> &nodes, const vector<Element> &elements, const vector<MatPro> &mats, const int &inum, const struct Force_Disp_TPB &force_disp, vector<double> &TPB_force_yy)const
{
	int boudary_eles_num = 0;
	double force_yy=0.0;
	//-----------------------------------------------------------------
	//查找边界节点
	for(int i=0; i<(int)elements.size(); i++)
	{
		vector<vector<double> > nod_strain;
		vector<int> nodj;
		for(int j=0; j<(int)elements[i].nodes_id.size(); j++)
		{
			if(nodes[elements[i].nodes_id[j]].x>=force_disp.cx0&&
			   nodes[elements[i].nodes_id[j]].x<=force_disp.cx1&&
			   nodes[elements[i].nodes_id[j]].y==force_disp.cy0)
			{
				//--------------------------------------------
				//Only for quadrangle
				if((int)elements[i].nodes_id.size()!=4)
				{
					hout << "Attention: the number of nodes in the element which is on the boundary isn't equal to 4." << endl;
					return 0;
				}

				//计算这一节点在此单元中的应变
				//--------------------------------------------

				//计算B矩阵；
				//--------------------------------------------
				//形函数N对某节点的偏导矩阵
				double nx=0.0, ny=0.0;
				switch(j)
				{
					case 0:{ nx=-1.0; ny=-1.0; break; }
					case 1:{ nx=1.0; ny=-1.0; break; }
					case 2:{ nx=1.0; ny=1.0; break; }
					case 3:{ nx=-1.0; ny=1.0; break; }
					default: hout << "Attention: the number of nodes isn't from 0 to 3" << endl;
				}
				double diff[2][4];
				diff[0][0]=-0.25*(1.0-ny);
				diff[0][1]=-diff[0][0];
				diff[0][2]=0.25*(1.0+ny);
				diff[0][3]=-diff[0][2];

				diff[1][0]=-0.25*(1.0-nx);
				diff[1][1]=-0.25*(1.0+nx);
				diff[1][2]=-diff[1][1];
				diff[1][3]=-diff[1][0];

				//--------------------------------------------------
				//单元节点坐标矩阵
				double elenode[4][2];
				for(int k=0; k<4; k++)
				{
					elenode[k][0]=nodes[elements[i].nodes_id[k]].x;
					elenode[k][1]=nodes[elements[i].nodes_id[k]].y;
				}
				//--------------------------------------------------
				//J矩阵
				double Jmatrix[2][2];
				//以上两个矩阵的积
				for(int k=0; k<2; k++)
					for(int m=0; m<2; m++)
					{
						Jmatrix[k][m]=0;
						for(int p=0; p<4; p++)
						Jmatrix[k][m] += diff[k][p]*elenode[p][m];
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
				//求出N对x,y,的偏导
				double diffxy[2][4];
				for(int k=0; k<2; k++)
					for(int m=0; m<4; m++)
					{
						diffxy[k][m]=0;
						for(int p=0; p<2; p++)
							diffxy[k][m] += Jinverse[k][p]*diff[p][m];
					}

				//--------------------------------------------------------
				//求出B矩阵
				double B[6][12];
				for(int k=0; k<6; k++)
					for(int m=0; m<12; m++)
						B[k][m]=0;
				for(int k=0; k<4; k++)
				{
					B[0][k*3+0]=diffxy[0][k];
					B[1][k*3+1]=diffxy[1][k];
					B[3][k*3+0]=diffxy[1][k];
					B[3][k*3+1]=diffxy[0][k];
					B[4][k*3+2]=diffxy[1][k];
					B[5][k*3+2]=diffxy[0][k];
				}		
				//--------------------------------------------
				//计算单元形心点应变值
				vector<double> ele_u(12,0.0);
				for(int k=0; k<4; k++)
				{
					ele_u[3*k] = u[3*elements[i].nodes_id[k]];
					ele_u[3*k+1] = u[3*elements[i].nodes_id[k]+1];
				}
				vector<double> temp_strain(6,0.0);
				for(int k=0; k<6; k++)
					for(int m=0; m<12; m++)
						temp_strain[k] += B[k][m]*ele_u[m];

				nod_strain.push_back(temp_strain);
				nodj.push_back(j);
			}
		}
		if((int)nodj.size()>2)
		{
			hout << "Attention: more than 2 nodes of the element are in this boundary!" << endl;
			return 0;
		}
		else if((int)nodj.size()==2)
		{
			//计算边长
			Point_3D p0(nodes[elements[i].nodes_id[nodj[0]]].x, nodes[elements[i].nodes_id[nodj[0]]].y, nodes[elements[i].nodes_id[nodj[0]]].z);
			Point_3D p1(nodes[elements[i].nodes_id[nodj[1]]].x, nodes[elements[i].nodes_id[nodj[1]]].y, nodes[elements[i].nodes_id[nodj[1]]].z);
			double ele_bs =p0.distance_to(p1);

			//计算弹性矩阵
			double ele_elas[6][6] = { {0}, {0}, {0}, {0}, {0}, {0} };
			double damage_differ = 1.0 - elements[i].dam;
			if(damage_differ>Zero)
			{
				for(int j=0; j<6; j++)
					for(int k=0; k<6; k++)
						ele_elas[j][k] = mats[1].elas_matrix[j][k]*damage_differ;
			}

			//计算单元节点平均应力
			double stress_yy=0.0;
			for(int j=0; j<(int)nod_strain.size(); j++)
				for(int k=0; k<6; k++) stress_yy += ele_elas[1][k]*nod_strain[j][k];
			force_yy += 0.5*stress_yy*ele_bs;
			boudary_eles_num++;
		}
	}

	TPB_force_yy[inum+1] = force_yy;

	return 1;
}
//---------------------------------------------------------------------------
//Find the same nodes seperated by DG elements 
int Postprocessor::Find_same_DG_nodes(const int &nnum, const vector<Node> &nodes, vector<int> &same_nods)const
{
	same_nods.push_back(nnum);

	if(nodes[nnum].DG_node)
	{
		Node Cnod = nodes[nnum];
		for(int i=0; i<(int)nodes.size(); i++)
		{
			if(nodes[i].DG_node&&i!=nnum)
			{
				if(Cnod.x==nodes[i].x&&Cnod.y==nodes[i].y&&Cnod.z==nodes[i].z)
				{
					same_nods.push_back(i);
				}
			}
		}
	}

	return 1;
}
//---------------------------------------------------------------------------
//输出Tecplot可视化位移场云图 for 2D problem
int Postprocessor::Export_Displacement_Contour_2D(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Displacement_Field_Contour" << endl;
	otec << "VARIABLES = X, Y";
	otec << ", Ux, Uy" << endl;
	
	otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size() << ", F=FEPOINT, ET=QUADRILATERAL" << endl;
	for(int i=0; i<(int)nodes.size(); i++)
	{
		otec << nodes[i].x << "  " << nodes[i].y;
		otec << "  " << u[3*i] << "  " << u[3*i+1];
		otec << endl;
	}
	otec << endl;
	for (int i=0; i < (int)elements.size(); i++)
	{
		otec	 << elements[i].nodes_id[0]+1 << "  " << elements[i].nodes_id[1]+1 << "  " 
				 << elements[i].nodes_id[2]+1 << "  " << elements[i].nodes_id[3]+1 << endl;
	}
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//计算单元形心处应变向量for 2D problem (Plane Strain Assumption)
int Postprocessor::Calculate_Elements_Strain_2D(const vector<Node> &nodes, const vector<Element> &elements)
{
	//---------------------------------------------------------------------------
	//计算单元形心点处的应变向量
	vector<double> temp_strain(6,0.0);
	ele_strain.assign((int)elements.size(), temp_strain);
	for(int i=0; i<(int)elements.size(); i++)
	{
		//计算B矩阵；
		//--------------------------------------------
		//形函数N对正方形中心点(0,0)的偏导矩阵	
		double diff[2][4];
		diff[0][0]=-0.25;
		diff[0][1]=-diff[0][0];                         
		diff[0][2]=0.25;
		diff[0][3]=-diff[0][2];

		diff[1][0]=-0.25;
		diff[1][1]=-0.25;
		diff[1][2]=-diff[1][1];
		diff[1][3]=-diff[1][0];

		//--------------------------------------------------
		//单元节点坐标矩阵
		double elenode[4][2];
		for(int j=0; j<4; j++)
		{
			elenode[j][0]=nodes[elements[i].nodes_id[j]].x;
			elenode[j][1]=nodes[elements[i].nodes_id[j]].y;
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
		double J_val = Jmatrix[0][0]*Jmatrix[1][1]-Jmatrix[0][1]*Jmatrix[1][0];
		//----------------------------------------------------
		//求出J矩阵的逆矩阵
		double Jinverse[2][2];
		Jinverse[0][0]=Jmatrix[1][1]/J_val;
		Jinverse[1][1]=Jmatrix[0][0]/J_val;
		Jinverse[0][1]=-Jmatrix[0][1]/J_val;
		Jinverse[1][0]=-Jmatrix[1][0]/J_val;

		//-------------------------------------------------------
		//求出N对x,y,的偏导
		double diffxy[2][4];
		for(int j=0; j<2; j++)
			for(int k=0; k<4; k++)
			{
				diffxy[j][k]=0;
				for(int m=0; m<2; m++)
					diffxy[j][k] += Jinverse[j][m]*diff[m][k];
			}

		//--------------------------------------------------------
		//求出B矩阵
		double B[6][12];
		for(int j=0; j<6; j++)
			for(int k=0; k<12; k++)
				B[j][k]=0;
		for(int j=0; j<4; j++)
		{
			B[0][j*3+0]=diffxy[0][j];
			B[1][j*3+1]=diffxy[1][j];
			B[3][j*3+0]=diffxy[1][j];
			B[3][j*3+1]=diffxy[0][j];
			B[4][j*3+2]=diffxy[1][j];
			B[5][j*3+2]=diffxy[0][j];
		}		
		//--------------------------------------------
		//计算单元形心点应变值
		vector<double> ele_u(12,0.0);
		for(int j=0; j<4; j++)
		{
			ele_u[3*j] = u[3*elements[i].nodes_id[j]];
			ele_u[3*j+1] = u[3*elements[i].nodes_id[j]+1];
		}
		for(int j=0; j<6; j++)
		{
			for(int k=0; k<12; k++)
			{
				ele_strain[i][j] = ele_strain[i][j] + B[j][k]*ele_u[k];
			}
		}
	}
		
	return 1;
}
//---------------------------------------------------------------------------
//计算单元形心处应力向量和应变能密度for 2D problem (Plane Strain Assumption)
int Postprocessor::Calculate_Elements_Stress_Energy_2D(const vector<MatPro> &mats, const vector<Element> &elements)
{
	//---------------------------------------------------------------------------
	//计算单元形心点处的应力向量和应变能密度
	vector<double> temp_stress(6, 0.0);
	ele_stress.assign((int)elements.size(), temp_stress);
	ele_energy.assign((int)elements.size(), 0.0);
	for (int i = 0; i < (int)elements.size(); i++)
	{
		for (int j = 0; j < 6; j++)
		{
			for (int k = 0; k < 6; k++)
			{
				ele_stress[i][j] += (1.0-elements[i].dam) * mats[elements[i].mat+1].elas_matrix[j][k] * ele_strain[i][k];
			}
		}
		for (int j = 0; j < 3; j++)	ele_energy[i] += 0.5*ele_stress[i][j] * ele_strain[i][j];
	}

	return 1;
}
//---------------------------------------------------------------------------
//输出应变场云图 for 2D problem
int Postprocessor::Export_Strain_Mesh_2D(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Strain_Field_Contour" << endl;
	otec << "VARIABLES = X, Y";
	otec << ", dam_en, break_en, total_dissip";
	otec << ", <greek>e</greek><sub>xx</sub>, <greek>e</greek><sub>yy</sub>, <greek>e</greek><sub>xy</sub>" << endl;

	//------------------------------------------------------------------------------------------------------------------------------------
	//输出
	otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size()  << ", F=FEPOINT, ET=QUADRILATERAL" << endl;
	for(int i=0; i<(int)nodes.size(); i++)
	{
		otec << setprecision(15) << nodes[i].x+u[3*i] << "  " << nodes[i].y+u[3*i+1];
		otec << setprecision(15) << "  " << nod_damage[i] << "  " << nod_dissip[i] << "  " << nod_damage[i]+nod_dissip[i];
		otec << setprecision(15) << "  " << nod_strain[i][0] << "  " << nod_strain[i][1] << "  " << nod_strain[i][3] << endl;
	}
	for (int i=0; i<(int)elements.size(); i++)
	{
		otec	 << elements[i].nodes_id[0]+1 << "  " << elements[i].nodes_id[1]+1 << "  " 
					<< elements[i].nodes_id[2]+1 << "  " << elements[i].nodes_id[3]+1 << endl;
	}
	otec << endl;

	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//输出应变场、应力场和应变能密度云图 for 2D problem
int Postprocessor::Export_Strain_Stress_Energy_Mesh_2D(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Strain_Field_Contour" << endl;
	otec << "VARIABLES = X, Y";
	otec << ", dam_en, break_en, total_dissip";
	otec << ", local_dam, nonlocal_dam, total_dam, nonloc_dissip_ratio";
	otec << ", <greek>e</greek><sub>xx</sub>, <greek>e</greek><sub>yy</sub>, <greek>e</greek><sub>xy</sub>" << endl;
	otec << ", <greek>s</greek><sub>xx</sub>, <greek>s</greek><sub>yy</sub>, <greek>s</greek><sub>xy</sub>" << endl;
	otec << ", e_d" << endl;

	//------------------------------------------------------------------------------------------------------------------------------------
	//输出
	otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size() << ", F=FEPOINT, ET=QUADRILATERAL" << endl;
	for (int i = 0; i<(int)nodes.size(); i++)
	{
		otec << setprecision(15) << nodes[i].x + u[3 * i] << "  " << nodes[i].y + u[3 * i + 1];
		otec << setprecision(15) << "  " << nod_damage[i] << "  " << nod_dissip[i] << "  " << nod_damage[i] + nod_dissip[i];
		otec << setprecision(15) << "  " << nod_local_dam[i] << "  " << nod_nonlocal_dam[i] << "  " << nod_local_dam[i]+nod_nonlocal_dam[i] << "  " << nod_dissip_ratio[i];
		otec << setprecision(15) << "  " << nod_strain[i][0] << "  " << nod_strain[i][1] << "  " << nod_strain[i][3];
		otec << setprecision(15) << "  " << nod_stress[i][0] << "  " << nod_stress[i][1] << "  " << nod_stress[i][3];
		otec << setprecision(15) << "  " << nod_energy[i] << endl;
	}
	for (int i = 0; i<(int)elements.size(); i++)
	{
		otec << elements[i].nodes_id[0] + 1 << "  " << elements[i].nodes_id[1] + 1 << "  "
			<< elements[i].nodes_id[2] + 1 << "  " << elements[i].nodes_id[3] + 1 << endl;
	}
	otec << endl;

	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//输出或读取文件数和加权函数
int Postprocessor::wr_files_num_weight_funcs(const string &output_file_name, const string &wr_mod, const struct Weight_func &weight_func, int &Unum)
{
	if(wr_mod=="Write") //输出数据
	{
		ofstream onum(output_file_name.c_str());
		onum << Unum << endl;
		onum << weight_func.keywords << endl;
		onum << weight_func.mark << endl;
		onum << weight_func.num << endl;
		for(int i=0; i<weight_func.num; i++) 
		{
			onum << weight_func.shape[i] << endl;
			onum << weight_func.center[i].x << endl;
			onum << weight_func.center[i].y << endl;
			onum << weight_func.center[i].z << endl;
			onum << weight_func.r0[i] << endl;
			onum << weight_func.r1[i] << endl;
			onum << weight_func.ratio[i] << endl;
			onum << weight_func.func_order[i] << endl;
			onum << weight_func.func_constant[i] << endl;
		}
		onum.close();
	}
	else if(wr_mod=="Read") //读入数据
	{
		ifstream inum(output_file_name.c_str());
		inum >> Unum;
		
		string keywords;
		inum >> keywords;
		postp_weight_func.keywords = keywords;
		
		bool mark;
		inum >> mark;
		postp_weight_func.mark = mark;

		int num;
		inum >> num;
		postp_weight_func.num = num;

		for(int i=0; i<postp_weight_func.num; i++) 
		{
			string strshape;
			inum >> strshape;
			postp_weight_func.shape.push_back(strshape);
			Point_3D poicenter;
			inum >> poicenter.x;
			inum >> poicenter.y;
			inum >> poicenter.z;
			postp_weight_func.center.push_back(poicenter);
			double r0, r1, ratio;
			inum >> r0;
			inum >> r1;
			inum >> ratio;
			postp_weight_func.r0.push_back(r0);
			postp_weight_func.r1.push_back(r1);
			postp_weight_func.ratio.push_back(ratio);
			string strforder;
			inum >> strforder;
			postp_weight_func.func_order.push_back(strforder);
			double funconst;
			inum >> funconst;
			postp_weight_func.func_constant.push_back(funconst);
		}
		inum.close();
	}
	else
	{
		hout << "注意！读写指令即不是读也不是写指令，请检查(wr_files_num_weight_funcs)！" << endl;
		return 0;
	}
	return 1;
}
//===========================================================================
