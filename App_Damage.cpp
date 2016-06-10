//============================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	App_Damage.cpp
//OBJECTIVE:	An adaptive algorithm between local and non-local models for the simulation of static damage and fracture propagation.
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//============================================================================================
#include "App_Damage.h"

int App_Damage::Application_damage_2D(Input *Init)const
{
	clock_t ct0,ct1; //time markers for cpu time of local operations

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Generate material database
	ct0 = clock();
	cout << "-_- Creating material database..." << endl;
	hout << "-_- Creating material database..." << endl;

	MatBase Material;
	if(Material.Generate_matbase_2D(Init->stif_nonloc, Init->nonloc_gsize, Init->nonloc_gau.num, Init->damages, Init->peri_para)==0) return 0;

	ct1 = clock();
	hout << "    The material database built in " << (double)(ct1-ct0)/CLOCKS_PER_SEC << "secs." << endl;
	hout << "^_^ The material database built successfully!" << endl << endl;
	cout << "^_^ The material database built successfully!" << endl << endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Define displacement solution of every computational step
	vector<vector<double> > U_tot;
	//Define nodes and elements solutions of every computational step
	vector<vector<Node> > nodes_tot;
	vector<vector<Element> > elements_tot;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//if "Read" mod, jump over Solving
	if(Init->rw_mod.type!="Read"&&Init->rw_mod.type=="Write")
	{
		//-----------------------------------------------------------------------------------------------------------------------------------------
		//Define the class of Global_Stiff_Matrix for the element matrices and total matrix
		Global_Stiff_Matrix Glosmat;
		//Define a table list of the broken bonds
		vector<bool> break_table[2];
		//Damaged energy of every element for continuum damage mechanics
		vector<double> damage_en[2];
		//Energy dissipation of every element for bond broken
		vector<double> dissip_en[2];
		//Define a table list of the damaged gaussian points
		vector<vector<double> > damage_table[2];
		//Define a table list of the damage force at gaussian points
		vector<vector<double> > Y_force[2];
		//Full damaged element
		vector<bool> full_dam_eles;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		time_t it0 = time(NULL); //Standard time
	
		cout<<endl;
		cout<<"------------------------------------"<<endl;
		cout<<"|            Iterations            |"<<endl;
		cout<<"------------------------------------"<<endl;
		cout<<endl;
		cout<<"-_- Starting iterative algorithm..."<<endl;
		cout<<endl;
		hout<<endl;
		hout<<"------------------------------------"<<endl;
		hout<<"|            Iterations            |"<<endl;
		hout<<"------------------------------------"<<endl;
		hout<<endl;
		hout<<"-_- Starting interative algorithm..."<<endl;
		hout<<endl;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//Loops
		int totalcount = 0;
		int count_ramp = 0;
		int ramp_num = 0;  //记录同一载荷步的循环次数(载荷步变化, 即count_ramp++, ramp_num清零)
		bool alpha_change = false;
		bool key_break_tab = false;
		while(count_ramp<Init->iter.ramp_para)
		{
			hout<<"--------------------------------------------------------------"<<endl;
			hout<<"---Ramp on boundary condition, loop "<<count_ramp << " --- " << ramp_num <<" on "<<Init->iter.ramp_para<<"---"<<endl;
			hout<<"--------------------------------------------------------------"<<endl;

			cout<<"--------------------------------------------------------------"<<endl;
			cout<<"---Ramp on boundary condition, loop "<<count_ramp << " --- " << ramp_num <<" on "<<Init->iter.ramp_para<<"---"<<endl;
			cout<<"--------------------------------------------------------------"<<endl;

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//The loading ratio of boundary conditions
			const double rmp = (double)(count_ramp+1)/Init->iter.ramp_para;

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Generate grids of RVE
			ct0 = clock();
			cout << "-_- Generating grids of RVE..." << endl;
			hout << "-_- Generating grids of RVE..." << endl;

			Mesher Mesh;
			if(Mesh.Import_mesh_reconfiguration_2D("mesh.dat", Init->mod_disc.disc, full_dam_eles, Init->weight_func) == 0) return 0; //Import mesh data from a file and reconfiguration mesh data(DGFEM)

			ct1 = clock();
			hout << "    Grids generated in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs."  << endl;
			hout << "^_^ Grids generated successfully!" << endl << endl;
			cout << "^_^ Grids generated successfully!" << endl << endl;

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Seek non-local neightbours of elements
			ct0 = clock();
			cout << "-_- Seeking non-local neighbours of elements..." << endl;
			hout << "-_- Seeking non-local neighbours of elements..." << endl;

			if(Mesh.Deter_relative_nodes_elements(Mesh.nodes, Mesh.elements, Init->peri_para.horizon_R, Init->mod_disc.mod, Init->cracks, Init->weight_func)==0) return 0; //update of non-local element neighbours

			SolveEqu Solv;
			vector<int> Iz, Ig;			//Iz and Ig for sparce storage and calculation
			if(Solv.izig(Mesh.nodes, Iz, Ig)==0) return 0;

			for(int i=0; i<(int)Mesh.nodes.size(); i++) Mesh.nodes[i].relative_nods.clear();		//Cleaning the relative nodes after function izig

			ct1 = clock();
			hout << "    Seeking neighbours done in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs."  << endl;
			hout << "^_^ Non-local neighbours calculated successfully!" << endl << endl;
			cout << "^_^ Non-local neighbours calculated successfully!" << endl << endl;

			if(count_ramp==0) //First time
			{
				//-----------------------------------------------------------------------------------------------------------------------------------------
				//Initialization of  damaged energy of every element
				damage_en[0].assign(Mesh.elements.size(), 0.0);
				//Initialization of  dissipative energy of every element
				dissip_en[0].assign(Mesh.elements.size(), 0.0);
				//Initialization of  full damaged element
				full_dam_eles.assign(Mesh.elements.size(), false);

				//-----------------------------------------------------------------------------------------------------------------------------------------
				//Initialization of break_table
				ct0 = clock();
				cout << "-_- Initialization of break_table..." << endl;
				hout << "-_- Initialization of break_table..." << endl;
	
				long int total_bonds=0;
				for(int i=0; i<(int)Mesh.elements.size(); i++)     total_bonds += (long int)Mesh.elements[i].relative_eles.size();
				int num_gauss = Init->gauss.num*Init->gauss.num; //the number of gaussian points in one element
				total_bonds = total_bonds*num_gauss*num_gauss;

				cout << "This simulation involves " << total_bonds << " bonds."<< endl;
				hout << "This simulation involves " << total_bonds << " bonds."<< endl;
				cout << "Memory required for the break_table is " << (double)total_bonds*32/1024/1024/1024/8 << "GB." << endl;
				hout << "Memory required for the break_table is " << (double)total_bonds*32/1024/1024/1024/8 << "GB." << endl;
	
				//---------------------------------------------------------------------------
				//assign break_table
				break_table[0].assign(total_bonds, false);

				ct1 = clock();
				hout << "    The break_table initialized in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
				hout << "^_^ The break_table initialized successfully!" << endl << endl;
				cout << "^_^ The break_table initialized successfully!" << endl << endl;

				//-----------------------------------------------------------------------------------------------------------------------------------------
				//Initialization of damage_table and damage_force
				ct0 = clock();
				cout << "-_- Initialization of damage_table and damage_force..." << endl;
				hout << "-_- Initialization of damage_table and damage_force..." << endl;
	
				vector<double> dam_temp;
				int num_gpoi = Init->gauss.num*Init->gauss.num; //the number of gaussian points in one element
				for(int i=0; i<(int)Mesh.elements.size(); i++)		
				{
					if(Mesh.elements[i].type==241)
					{
						dam_temp.assign(num_gpoi, 0.0);
						damage_table[0].push_back(dam_temp);
						Y_force[0].push_back(dam_temp);
					}
					else 
					{
						hout << "Error: the element type " << Mesh.elements[i].type << " is not defined!" << endl;
						cout << "Error: the element type " << Mesh.elements[i].type << " is not defined!" << endl;
						return 0;
					}
				}

				ct1 = clock();
				hout << "    The damage_table initialized in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
				hout << "^_^ The damage_table initialized successfully!" << endl << endl;
				cout << "^_^ The damage_table initialized successfully!" << endl << endl;
			}
			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Reassign break_table
			if(Init->weight_func.num!=1&&Init->weight_func.shape[0]=="Null"&&!key_break_tab)
			{
				key_break_tab = true; //Initialized only at the first time
				//-----------------------------------------------------------------------------------------------------------------------------------------
				//Reinitialization of break_table
				ct0 = clock();
				cout << "-_- Reinitialization of break_table..." << endl;
				hout << "-_- Reinitialization of break_table..." << endl;
	
				long int total_bonds=0;
				for(int i=0; i<(int)Mesh.elements.size(); i++)     total_bonds += (long int)Mesh.elements[i].relative_eles.size();
				int num_gauss = Init->gauss.num*Init->gauss.num; //the number of gaussian points in one element
				total_bonds = total_bonds*num_gauss*num_gauss;

				cout << "This simulation involves " << total_bonds << " bonds."<< endl;
				hout << "This simulation involves " << total_bonds << " bonds."<< endl;
				cout << "Memory required for the break_table is " << (double)total_bonds*32/1024/1024/1024/8 << "GB." << endl;
				hout << "Memory required for the break_table is " << (double)total_bonds*32/1024/1024/1024/8 << "GB." << endl;
	
				//---------------------------------------------------------------------------
				//assign break_table
				break_table[0].assign(total_bonds, false);

				ct1 = clock();
				hout << "    The break_table reinitialized in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
				hout << "^_^ The break_table reinitialized successfully!" << endl << endl;
				cout << "^_^ The break_table reinitialized successfully!" << endl << endl;
			}

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Define the matrices of elements
			MathMatrix tem_mat(12,12);
			vector<MathMatrix> ele_self_matrix[2];
			ele_self_matrix[0].assign(Mesh.elements.size(), tem_mat);
			vector<vector<MathMatrix> > ele_relative_matrix[2];
			for(int i=0; i<(int)Mesh.elements.size(); i++)
			{
				vector<MathMatrix> tem_vec_mat(Mesh.elements[i].relative_eles.size(), tem_mat);
				ele_relative_matrix[0].push_back(tem_vec_mat);
			}

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Estimate element matrices
			ct0 = clock();
			cout <<"-_- Calculating new element matrices..." << endl;
			hout <<"-_- Calculating new element matrices..." << endl;

			vector<vector<double> > Gp_val; //The values of weighting function of gaussian points in all elements
			vector<double> temp_ak_vec(Init->gauss.num*Init->gauss.num, -1.0);
			vector<vector<double> > Ak_val(Mesh.elements.size(), temp_ak_vec); //The values of alpha_key of gaussian points in all elements
			vector<double> temp_stiff(36*Init->gauss.num*Init->gauss.num, 0.0);     //The elastic matrix 6*6 at every gaussian point in a element
			vector<vector<double> > Ele_local_stiff(Mesh.elements.size(), temp_stiff);
			if(Glosmat.Gen_element_matrices_damage_2D(Init->gauss.num, Init->peri_para, Init->mod_disc.mod, Init->weight_func, Init->damages, Material.mats_vec, Mesh.nodes, Mesh.elements, 
																					full_dam_eles, damage_table[0], break_table[0], Gp_val, Ak_val, Ele_local_stiff, ele_self_matrix[0], ele_relative_matrix[0])==0) return 0;
			ct1 = clock();
			hout << "    Element matrices calculated in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
			hout << "^_^ Element matrices calculated successfully!" << endl << endl;
			cout << "^_^ Element matrices calculated successfully!" << endl << endl;

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Build load vector
			ct0 = clock();
			cout <<"-_- Building load vector..." << endl;
			hout <<"-_- Building load vector..." << endl;

			vector<double> equright(3*Mesh.nodes.size(), 0.0);
			Global_Load_Vector Gloload;
			if(Gloload.Gen_global_load_vector_2D(Init->load, rmp, Mesh.nodes, Mesh.elements, equright)==0) return 0;

			ct1 = clock();
			hout << "    Load vector built in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
			hout << "^_^ Load vector built successfully!" << endl << endl;
			cout << "^_^ Load vector built successfully!" << endl << endl;

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Build displacement constraints
			ct0 = clock();
			cout <<"-_- Building displacement constraints..." << endl;
			hout <<"-_- Building displacement constraints..." << endl;

			int bnod_num = 0;		//the number of nodes on the boundary
			vector<int> ip;			//the sign of boudary nodes
			vector<double> vp;	//the value of boudary nodes
			if(Solv.Fixed_displacement_constraints_2D(Init->displace, rmp, Mesh.nodes, bnod_num, ip, vp)==0) return 0;
			ct1 = clock();

			hout << "    Displacement constraints built in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs."  << endl;
			hout << "^_^ Displacement constraints built successfully!" << endl << endl;
			cout << "^_^ Displacement constraints built successfully!" << endl << endl;

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Copy variables
			ele_self_matrix[1] = ele_self_matrix[0];
			ele_relative_matrix[1] = ele_relative_matrix[0];
			damage_table[1] = damage_table[0];
			Y_force[1] = Y_force[0];
			break_table[1] = break_table[0];
			damage_en[1] = damage_en[0];
			dissip_en[1] = dissip_en[0];

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Solving equations and Updating element matrics
			bool flag_break = true;
			bool flag_dam = true;
			alpha_change=false;
			int count_iter = 0;
			int dam_break_order = 0;  //For judging to update damaged elements or broken bonds
			vector<double> U_Solution;
			vector<bool> dam_iter(Mesh.elements.size(), false);
			vector<bool> break_iter(Mesh.elements.size(), false);

			//-----------------------------------------------------------------------------------------------------------------------------------------
			while((flag_dam||flag_break)&&count_iter<Init->iter.max_iter)
			{
				hout<<endl;
				hout <<"-_- iteration "<<count_iter << " updated item " << dam_break_order%2 <<endl;
				hout <<"----------------"<<endl;
				hout<<endl;
				cout<<endl;
				cout <<"iteration "<<count_iter << " updated item " << dam_break_order%2<<endl;
				cout <<"------------"<<endl;
				cout<<endl;
				//-------------------------------------------------------------------------------------------------------

				totalcount++; //to record the total count of iterations
				//-------------------------------------------------------------------------------------------------------
				ct0 = clock();
				cout << "-_- Assembling global stiffness matrix..." << endl;
				hout << "-_- Assembling global stiffness matrix..." << endl;

				vector<double> total_matrix(6*(int)Mesh.nodes.size()+9*Iz.back(), 0);
				if(Glosmat.Update_global_matrix(ele_self_matrix[1], ele_relative_matrix[1], Mesh.elements, Iz, Ig, total_matrix)==0) return 0;

				ct1 = clock();
				hout << "    Global matrix assembled in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
				hout << "^_^ Global matrix assembled successfully!" << endl << endl;
				cout << "^_^ Global matrix assembled successfully!" << endl << endl;

				//-------------------------------------------------------------------------------------------------------
				ct0 = clock();
				cout << "-_- Solving equations..." << endl;
				hout << "-_- Solving equations..." << endl;

				Solv.Deal_with_displacement_zero_value(bnod_num, (int)Mesh.nodes.size(), Iz, Ig, ip, vp, equright, total_matrix);

				U_Solution.assign(3*Mesh.nodes.size(), 0.0);
				Solv.Solve_linear_equations(bnod_num, (int)Mesh.nodes.size(), Iz, Ig, ip, vp, total_matrix, equright, U_Solution);
//				Solv.Solve_linear_equations_omp(bnod_num, (int)Mesh.nodes.size(), Iz, Ig, ip, vp, total_matrix, equright, U_Solution);

				ct1 = clock();
				hout << "    Linear equations solved in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
				hout << "^_^ Linear equations solved successfully!" << endl << endl;
				cout << "^_^ Linear equations solved successfully!" << endl << endl;

				//-------------------------------------------------------------------------------------------------------
				ct0 = clock();
				cout << "-_- Updating element matrices..." << endl;
				hout << "-_- Updating element matrices..." << endl;

				//-------------------------------------------------------------------------------------------------------
				if (dam_break_order%2==0)
				{
					int damaged_sum = 0;
					vector<MathMatrix> temp_self_matrix(ele_self_matrix[0]);
					vector<double> temp_damage_en(damage_en[0]);

					if (Glosmat.Update_damaged_matrices_2D(Init->gauss.num, Init->damages, Init->peri_para, Init->mod_disc.mod, Init->weight_func, Material.mats_vec, Mesh.nodes, Mesh.elements, U_Solution, 
																		full_dam_eles, Ak_val, Ele_local_stiff, break_table, Gp_val, temp_self_matrix, damage_table, Y_force, temp_damage_en, damaged_sum, dam_iter) == 0) return 0;
				
					if(Init->mod_disc.mod == "Local")  flag_break = false;
					else if(Init->mod_disc.mod == "Hybrid")
					{
						if(Init->weight_func.num==1&&Init->weight_func.shape[0]=="Null") flag_break = false;
						else dam_break_order++;
					}
					else{ hout << "Attention: a wrong word of mod is used." << endl;  return 0; }

					if(damaged_sum==0)
					{
						flag_dam = false;
						hout << "    !!damcov " << damaged_sum << " gaussian points were damaged." << endl;
					}
					else
					{
						flag_dam = true;
						hout << "    !!G " << damaged_sum << " gaussian points were damaged." << endl;

						ele_self_matrix[1] = temp_self_matrix;
						ele_relative_matrix[1] = ele_relative_matrix[0];
						damage_en[1] = temp_damage_en;

						if(dam_break_order%2==1) continue;
					}
				}

				if(dam_break_order%2==1)
				{
					int broken_sum=0;
					vector<MathMatrix> temp_self_matrix(ele_self_matrix[0]);
					vector<vector<MathMatrix> > temp_relative_matrix(ele_relative_matrix[0]);
					vector<double> temp_dissip_en(dissip_en[0]);

					if(Glosmat.Update_nonlocal_matrices_2D(Init->gauss.num, Init->damages, Init->peri_para, Init->mod_disc.mod, Init->weight_func, Mesh.nodes, Gp_val, Ak_val, U_Solution, 
																	Mesh.elements, damage_table[0], temp_self_matrix, temp_relative_matrix, break_table, temp_dissip_en, broken_sum, break_iter)==0) return 0;
					
					dam_break_order++;
					if(broken_sum==0)
					{
						flag_break = false;
						hout << "    !!breakcov " << broken_sum << " bonds were broken." << endl;
					}
					else 
					{
						flag_break = true;
						hout << "    !!B " << broken_sum << " bonds were broken." << endl;
	
						ele_self_matrix[1] = temp_self_matrix;
						ele_relative_matrix[1] = temp_relative_matrix;
						dissip_en[1] = temp_dissip_en;
					}
				}

				//-------------------------------------------------------------------------------------------------------
				ct1 = clock();
				hout << "    " << "Elment matrices update in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
				hout << "^_^ Element matrices updated successfully!" << endl << endl;
				cout << "^_^ Element matrices updated successfully!" << endl << endl;
				//-------------------------------------------------------------------------------------------------------
				count_iter++;
			}

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Copy variables
			damage_table[0] = damage_table[1];
			Y_force[0] = Y_force[1];
			break_table[0] = break_table[1];
			damage_en[0] = damage_en[1];
			dissip_en[0] = dissip_en[1];

			//-------------------------------------------------------------------------------------------------------
			//Record dissipative energy of every element
			for(int i=0; i<(int)Mesh.elements.size(); i++)
			{
				Mesh.elements[i].Damage_energy = damage_en[1][i];
				Mesh.elements[i].Dissipative_energy = dissip_en[1][i];
			}

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Updating weighting function: alpha, and updating full damaged elements 
			if(Init->mod_disc.mod=="Hybrid"&&(count_iter!=1||(Init->iter.max_iter==1&&(flag_break||flag_dam))))
			{
				cout << "-_- Updating alpha weighting function and full_damaged_eles..." << endl;
				hout << "-_- Updating alpha weighting function and full_damaged_eles..." << endl;

				WeightFunc Wei_Fun;
				if(Wei_Fun.Update_weighting_function_alpha(Mesh.nodes, Mesh.elements, Init->peri_para, Init->damages, dam_iter, damage_table[0], full_dam_eles, Init->weight_func, alpha_change)==0) return 0;

				if(alpha_change)
				{
					ramp_num++;
					hout << "    Weighting functions alpha and full_damaged_eles are updated." << endl;
					hout << "^_^ Weighting functions alpha and full_damaged_eles are updated successfully!" << endl << endl;
					cout << "^_^ Weighting functions alpha and full_damaged_eles are updated successfully!" << endl << endl;
				}
				else
				{
					hout << "    Weighting functions alpha and full_damaged_eles did not change." << endl;
					hout << "^_^ Weighting functions alpha and full_damaged_eles tested successfully!" << endl << endl;
					cout << "^_^ Weighting functions alpha and full_damaged_eles tested successfully!" << endl << endl;
				}
			}

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Record values of every step
			if(!alpha_change)
			{
				//---------------------------------------------------------------------------
				//Output data
				Gauss gau;		//gaussian potins
				if(gau.Generate_gauss_2D(Init->gauss.num)==0) return 0;

				int utsize = (int)U_tot.size();
				//---------------------------------------------------------------------------
				//Output weighting function contour
				stringstream owfc;
				if(utsize<10)	owfc << "./Results/Weight_Func_Contour_000" << utsize << ".dat";
				else if (utsize<100)	owfc << "./Results/Weight_Func_Contour_00" << utsize << ".dat";
				else if (utsize<1000)	owfc << "./Results/Weight_Func_Contour_0" << utsize << ".dat";
				else	owfc << "./Results/Weight_Func_Contour_" << utsize << ".dat";

//				WeightFunc Wei_Fun;
//				Wei_Fun.Output_weight_function_contour_2D(owfc.str(), Mesh.nodes, Mesh.elements, full_dam_eles, gau.gauss, gau.weight, Init->weight_func);

				//Output weighting function values
				//stringstream owfv;				
				//if(utsize<10)	owfv << "./Results/Weight_Func_Vector_000" << utsize << ".dat";
				//else if (utsize<100)	owfv << "./Results/Weight_Func_Vector_00" << utsize << ".dat";
				//else if (utsize<1000)	owfv << "./Results/Weight_Func_Vector_0" << utsize << ".dat";
				//else	owfv << "./Results/Weight_Func_Vector_" << utsize << ".dat";
				//Wei_Fun.Output_weight_function_vectors(owfv.str(), Init->weight_func);

				//Output damage contour(attention: we add an average damage value of element into Mesh.elements)
				stringstream odac;
				if(utsize<10)	odac << "./Results/Damage_Contour_000" << utsize << ".dat";
				else if (utsize<100)	odac << "./Results/Damage_Contour_00" << utsize << ".dat";
				else if (utsize<1000)	odac << "./Results/Damage_Contour_0" << utsize << ".dat";
				else	odac << "./Results/Damage_Contour_" << utsize << ".dat";
				Output_damage_contour_2D(odac.str(), Mesh.nodes, Mesh.elements, gau.gauss, gau.weight, damage_table[0]);     //注意此单元中计算了单元的平均损伤值

				//Immediately output result data
				Postprocessor Output(U_Solution, Mesh.nodes, Mesh.elements);
				if (Output.Immedi_write_data(utsize, Init->rw_mod.type, Init->weight_func) == 0) return 0;

				//---------------------------------------------------------------------------
				//storage data
				U_tot.push_back(U_Solution);
				nodes_tot.push_back(Mesh.nodes);
				elements_tot.push_back(Mesh.elements);
				count_ramp++;
				ramp_num = 0; //记录同一载荷步的循环次数(载荷步变化, 即count_ramp++, ramp_num清零)
			}
		}

		//-----------------------------------------------------------------------------------------------------------------------------------------
		time_t it1 = time(NULL); //Standard time
		cout << endl;
		cout << "----------------------------------------------------------------"<<endl;
		cout << "Iterative algorithm finished in " << totalcount <<" iterations."<<endl;
		cout << "It took about " << (int)(it1 - it0) << "secs."  << endl;
		cout << "Iterative algorithm done successfully!" <<endl;
		cout << endl;
		hout << endl;
		hout << "----------------------------------------------------------------"<<endl;
		hout << "Iterative algorithm finished in " << totalcount <<" iterations."<<endl;
		hout << "It took about " << (int)(it1 - it0) << "secs."  << endl;
		hout << "Iterative algorithm done successfully!" <<endl;
		hout << endl;
	}
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Postprocess
	ct0 = clock();
	cout << endl;
	cout << "-------------------------------------------------" << endl;
	cout << "|               Postprocessing                  |" << endl;
	cout << "-------------------------------------------------" << endl;
	cout << endl;
	cout << "-_- Starting postprocessor..." << endl;
	cout << endl;
	hout << endl;
	hout << "-------------------------------------------------" << endl;
	hout << "|               Postprocessing                  |" << endl;
	hout << "-------------------------------------------------" << endl;
	hout << endl;
	hout << "-_- Starting postprocessor..." << endl;
	hout << endl;
	
	//-------------------------------------------------------------------------------------------------------
	Postprocessor Post;
	if(Post.Treatment_2D(Init, nodes_tot, elements_tot, Material.mats_vec, U_tot, Init->rw_mod.type)==0) return 0; //执行后处理过程	
	
	//-------------------------------------------------------------------------------------------------------
	ct1 = clock();
	cout << endl;
	cout << "----------------------------------------------------------------"<<endl;
	cout << "    Postprocessor took about " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs."  << endl;
	cout << "^_^ Postprocessor done successfully!" <<endl;
	cout << endl;
	hout << endl;
	hout << "----------------------------------------------------------------"<<endl;
	hout << "    Postprocessor took about " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs."  << endl;
	hout << "^_^ Postprocessor done successfully!" <<endl;
	hout << endl;

	return 1;
}
int App_Damage::Application_damage(Input *Init)const
{
	clock_t ct0,ct1; //time markers for cpu time of local operations

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Generate material database
	ct0 = clock();
	cout << "-_- Creating material database..." << endl;
	hout << "-_- Creating material database..." << endl;

	MatBase Material;
	if(Material.Generate_matbase(Init->stif_nonloc, Init->nonloc_gsize, Init->nonloc_gau.num, Init->damages, Init->peri_para)==0) return 0;

	ct1 = clock();
	hout << "    The material database built in " << (double)(ct1-ct0)/CLOCKS_PER_SEC << "secs." << endl;
	hout << "^_^ The material database built successfully!" << endl << endl;
	cout << "^_^ The material database built successfully!" << endl << endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Define displacement solution of every computational step
	vector<vector<double> > U_tot;
	//Define dissipative energy of every computational step
	vector<vector<double> > Dissipation;
	//Define nodes and elements solutions of every computational step
	vector<vector<Node> > nodes_tot;
	vector<vector<Element> > elements_tot;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//if "Read" mod, jump over Solving
	if(Init->rw_mod.type!="Read"&&Init->rw_mod.type=="Write")
	{
		//-----------------------------------------------------------------------------------------------------------------------------------------
		//Define the class of Global_Stiff_Matrix for the element matrices and total matrix
		Global_Stiff_Matrix Glosmat;
		//Define a table list of the broken bonds
		vector<bool> break_table;
		//Energy dissipation of every element
		vector <double> dissip_en;
		//Define a table list of the damaged gaussian points
		vector<vector<double> > damage_table;
		//Full damaged element
		vector<bool> full_dam_eles;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		time_t it0 = time(NULL); //Standard time
	
		cout<<endl;
		cout<<"------------------------------------"<<endl;
		cout<<"|            Iterations            |"<<endl;
		cout<<"------------------------------------"<<endl;
		cout<<endl;
		cout<<"-_- Starting iterative algorithm..."<<endl;
		cout<<endl;
		hout<<endl;
		hout<<"------------------------------------"<<endl;
		hout<<"|            Iterations            |"<<endl;
		hout<<"------------------------------------"<<endl;
		hout<<endl;
		hout<<"-_- Starting interative algorithm..."<<endl;
		hout<<endl;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//Loops
		int totalcount = 0;
		int count_ramp = 0;
		while(count_ramp<Init->iter.ramp_para)
		{
			hout<<"--------------------------------------------------------------"<<endl;
			hout<<"---Ramp on boundary condition, loop "<<count_ramp <<" on "<<Init->iter.ramp_para<<"---"<<endl;
			hout<<"--------------------------------------------------------------"<<endl;

			cout<<"--------------------------------------------------------------"<<endl;
			cout<<"---Ramp on boundary condition, loop "<<count_ramp <<" on "<<Init->iter.ramp_para<<"---"<<endl;
			cout<<"--------------------------------------------------------------"<<endl;

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//The loading ratio of boundary conditions
			const double rmp = (double)(count_ramp+1)/Init->iter.ramp_para;

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Generate grids of RVE
			ct0 = clock();
			cout << "-_- Generating grids of RVE..." << endl;
			hout << "-_- Generating grids of RVE..." << endl;

			Mesher Mesh;
			if(Mesh.Import_mesh_reconfiguration("mesh.dat", Init->mod_disc.disc, full_dam_eles)==0) return 0; //Import mesh data from a file and reconfiguration mesh data(DGFEM)

//			if(Mesh.Import_mesh_reconfiguration("mesh.dat", Init->mod_disc.disc, Init->weight_func)==0) return 0; //Import mesh data from a file and reconfiguration mesh data(DGFEM) for pure nonlocal fracture propagation
//			if(Mesh.Generate_mesh(Init->grid_size, Init->ele_prop, Init->geom_rve, Init->mod_disc.disc, full_dam_eles)== 0) return 0;  //Generate background regular mesh
			
//			if(Mesh.Generate_mesh(Init->grid_size, Init->ele_prop, Init->cracks, Init->geom_rve, Init->mod_disc.disc, full_dam_eles)== 0) return 0;  //Generate background regular mesh
//			if(Mesh.Generate_mesh(Init->grid_size, Init->ele_prop, Init->cracks, Init->geom_rve, Init->mod_disc.disc, Init->weight_func)== 0) return 0;  //Generate background regular mesh

			ct1 = clock();
			hout << "    Grids generated in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs."  << endl;
			hout << "^_^ Grids generated successfully!" << endl << endl;
			cout << "^_^ Grids generated successfully!" << endl << endl;
			
			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Seek non-local neightbours of elements
			ct0 = clock();
			cout << "-_- Seeking non-local neighbours of elements..." << endl;
			hout << "-_- Seeking non-local neighbours of elements..." << endl;

			if(Mesh.Deter_relative_nodes_elements(Mesh.nodes, Mesh.elements, Init->peri_para.horizon_R, Init->mod_disc.mod, Init->cracks)==0) return 0; //update of non-local element neighbours

			SolveEqu Solv;
			vector<int> Iz, Ig;			//Iz and Ig for sparce storage and calculation
			if(Solv.izig(Mesh.nodes, Iz, Ig)==0) return 0;
	
			for(int i=0; i<(int)Mesh.nodes.size(); i++) Mesh.nodes[i].relative_nods.clear();		//Cleaning the relative nodes after function izig

			ct1 = clock();
			hout << "    Seeking neighbours done in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs."  << endl;
			hout << "^_^ Non-local neighbours calculated successfully!" << endl << endl;
			cout << "^_^ Non-local neighbours calculated successfully!" << endl << endl;
	
			if(count_ramp==0) //First time
			{
				//-----------------------------------------------------------------------------------------------------------------------------------------
				//Initialization of  dissipative energy of every element
				dissip_en.assign(Mesh.elements.size(), 0.0);
				//Initialization of  full damaged element
				full_dam_eles.assign(Mesh.elements.size(), false);

				//-----------------------------------------------------------------------------------------------------------------------------------------
				//Initialization of break_table
				ct0 = clock();
				cout << "-_- Initialization of break_table..." << endl;
				hout << "-_- Initialization of break_table..." << endl;
	
				long int total_bonds=0;
				for(int i=0; i<(int)Mesh.elements.size(); i++)		total_bonds += (long int)Mesh.elements[i].relative_eles.size();
				int num_gauss = Init->gauss.num*Init->gauss.num*Init->gauss.num; //the number of gaussian points in one element
				total_bonds = total_bonds*num_gauss*num_gauss;
	
				cout << "This simulation involves " << total_bonds << " bonds."<< endl;
				hout << "This simulation involves " << total_bonds << " bonds."<< endl;
				cout << "Memory required for the break_table is " << (double)total_bonds/1024/1024/1024/8 << "GB." << endl;
				hout << "Memory required for the break_table is " << (double)total_bonds/1024/1024/1024/8 << "GB." << endl;
	
				//---------------------------------------------------------------------------
				//assign break_table
				break_table.assign(total_bonds, false);

				ct1 = clock();
				hout << "    The break_table initialized in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
				hout << "^_^ The break_table initialized successfully!" << endl << endl;
				cout << "^_^ The break_table initialized successfully!" << endl << endl;

				//-----------------------------------------------------------------------------------------------------------------------------------------
				//Initialization of damage_table
				ct0 = clock();
				cout << "-_- Initialization of damage_table..." << endl;
				hout << "-_- Initialization of damage_table..." << endl;
	
				vector<double> dam_temp;
				int num_gpoi_brick = Init->gauss.num*Init->gauss.num*Init->gauss.num; //the number of gaussian points in brick element
				for(int i=0; i<(int)Mesh.elements.size(); i++)		
				{
					if(Mesh.elements[i].type==381)
					{
						dam_temp.assign(num_gpoi_brick, 0.0);
						damage_table.push_back(dam_temp);
					}
					else 
					{
						hout << "Error: the element type " << Mesh.elements[i].type << " is not defined!" << endl;
						cout << "Error: the element type " << Mesh.elements[i].type << " is not defined!" << endl;
						return 0;
					}
				}

				ct1 = clock();
				hout << "    The damage_table initialized in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
				hout << "^_^ The damage_table initialized successfully!" << endl << endl;
				cout << "^_^ The damage_table initialized successfully!" << endl << endl;
			}

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Assign values of dissipative energy of every element
			for(int i=0; i<(int)Mesh.elements.size(); i++) Mesh.elements[i].Dissipative_energy = dissip_en[i];

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Define the matrices of elements
			MathMatrix tem_mat(24,24);
			vector<MathMatrix> ele_self_matrix(Mesh.elements.size(), tem_mat);
			vector<vector<MathMatrix> > ele_relative_matrix;
			for(int i=0; i<(int)Mesh.elements.size(); i++)
			{
				vector<MathMatrix> tem_vec_mat(Mesh.elements[i].relative_eles.size(), tem_mat);
				ele_relative_matrix.push_back(tem_vec_mat);
			}

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Estimate element matrices
			ct0 = clock();
			cout <<"-_- Calculating new element matrices..." << endl;
			hout <<"-_- Calculating new element matrices..." << endl;

			vector<vector<double> > Gp_val; //The values of weighting function of gaussian points in all elements
			if(Glosmat.Gen_element_matrices_damage(Init->gauss.num, Init->peri_para, Init->mod_disc.mod, Init->weight_func, Init->damages, Material.mats_vec, 
																						Mesh.nodes, Mesh.elements, full_dam_eles, damage_table, break_table, Gp_val, ele_self_matrix, ele_relative_matrix)==0) return 0;

			ct1 = clock();
			hout << "    Element matrices calculated in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
			hout << "^_^ Element matrices calculated successfully!" << endl << endl;
			cout << "^_^ Element matrices calculated successfully!" << endl << endl;

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Build load vector
			ct0 = clock();
			cout <<"-_- Building load vector..." << endl;
			hout <<"-_- Building load vector..." << endl;

			vector<double> equright(3*Mesh.nodes.size(), 0.0);
			Global_Load_Vector Gloload;
			if(Gloload.Gen_global_load_vector(Init->load, rmp, Mesh.nodes, Mesh.elements, equright)==0) return 0;

			ct1 = clock();
			hout << "    Load vector built in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
			hout << "^_^ Load vector built successfully!" << endl << endl;
			cout << "^_^ Load vector built successfully!" << endl << endl;

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Build displacement constraints
			ct0 = clock();
			cout <<"-_- Building displacement constraints..." << endl;
			hout <<"-_- Building displacement constraints..." << endl;

			int bnod_num = 0;		//the number of nodes on the boundary
			vector<int> ip;			//the sign of boudary nodes
			vector<double> vp;	//the value of boudary nodes
			if(Solv.Fixed_displacement_constraints(Init->displace, rmp, Mesh.nodes, bnod_num, ip, vp)==0) return 0;
			ct1 = clock();

			hout << "    Displacement constraints built in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs."  << endl;
			hout << "^_^ Displacement constraints built successfully!" << endl << endl;
			cout << "^_^ Displacement constraints built successfully!" << endl << endl;

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Solving equations and Updating element matrics
			bool flag_break = true;
			bool flag_dam = true;
			bool alpha_change=false;
			bool full_damaged_change = false;
			int count_iter = 0;
			int dam_break_order = 0;  //For judging to update damaged elements or broken bonds
			vector<double> U_Solution;
			vector<bool> dam_iter(Mesh.elements.size(), false);
			vector<bool> break_iter(Mesh.elements.size(), false);
			while((flag_dam||flag_break)&&count_iter<Init->iter.max_iter)
			{
				hout <<"-_- iteration "<<count_iter << " updated item " << dam_break_order%2 <<endl;
				hout <<"----------------"<<endl;
				hout<<endl;
				cout <<"iteration "<<count_iter << " updated item " << dam_break_order%2<<endl;
				cout <<"------------"<<endl;
				cout<<endl;
				//-------------------------------------------------------------------------------------------------------

				totalcount++; //to record the total count of iterations
				//-------------------------------------------------------------------------------------------------------
				ct0 = clock();
				cout << "-_- Assembling global stiffness matrix..." << endl;
				hout << "-_- Assembling global stiffness matrix..." << endl;

				vector<double> total_matrix(6*(int)Mesh.nodes.size()+9*Iz.back(), 0);
				if(Glosmat.Update_global_matrix(ele_self_matrix, ele_relative_matrix, Mesh.elements, Iz, Ig, total_matrix)==0) return 0;
			
				ct1 = clock();
				hout << "    Global matrix assembled in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
				hout << "^_^ Global matrix assembled successfully!" << endl << endl;
				cout << "^_^ Global matrix assembled successfully!" << endl << endl;

				//-------------------------------------------------------------------------------------------------------
				ct0 = clock();
				cout << "-_- Solving equations..." << endl;
				hout << "-_- Solving equations..." << endl;

				Solv.Deal_with_displacement_zero_value(bnod_num, (int)Mesh.nodes.size(), Iz, Ig, ip, vp, equright, total_matrix);

				U_Solution.assign(3*Mesh.nodes.size(), 0.0);
				Solv.Solve_linear_equations(bnod_num, (int)Mesh.nodes.size(), Iz, Ig, ip, vp, total_matrix, equright, U_Solution);

				ct1 = clock();
				hout << "    Linear equations solved in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
				hout << "^_^ Linear equations solved successfully!" << endl << endl;
				cout << "^_^ Linear equations solved successfully!" << endl << endl;

				//-------------------------------------------------------------------------------------------------------
				ct0 = clock();
				cout << "-_- Updating element matrices..." << endl;
				hout << "-_- Updating element matrices..." << endl;

				//-------------------------------------------------------------------------------------------------------
				if(dam_break_order%2==0)
				{
					int damaged_sum=0;
					if(Glosmat.Update_damaged_element_matrices(Init->gauss.num, Init->damages, Material.mats_vec, Mesh.nodes, Mesh.elements,
																							U_Solution, full_dam_eles, ele_self_matrix, damage_table, damaged_sum, dam_iter)==0) return 0;

					dam_break_order++;
					if(damaged_sum==0)
					{
						hout << "    No gaussian point was damaged." << endl;
						flag_dam = false;
					}
					else	
					{
						//-------------------------------------------------------------------------------------------------------
						//Clear ele_self_matrix and ele_relative_matrix
						//MathMatrix mat_zero(24,24);
						//for(int i=0; i<(int)Mesh.elements.size(); i++)
						//{
						//	ele_self_matrix[i] = mat_zero;
						//	for(int j=0; j<(int)Mesh.elements[i].relative_eles.size(); j++) ele_relative_matrix[i][j] = mat_zero;
						//}

						//Gp_val.clear();
						//if(Glosmat.Gen_element_matrices_damage(Init->gauss.num, Init->peri_para, Init->mod_disc.mod, Init->weight_func, Init->damages, Material.mats_vec, 
						//																			Mesh.nodes, Mesh.elements, full_dam_eles, damage_table, break_table, Gp_val, ele_self_matrix, ele_relative_matrix)==0) return 0;

						hout << "    !!G " << damaged_sum << " gaussian points were damaged." << endl;
						continue;
					}
				}

				if(dam_break_order%2==1)
				{
					int broken_sum=0;
					if(Glosmat.Update_nonlocal_element_matrices(Init->gauss.num, Init->damages, Init->peri_para, Init->mod_disc.mod, Init->weight_func, Material.mats_vec, Mesh.nodes,
																									   Gp_val, U_Solution, Mesh.elements, damage_table, ele_self_matrix, ele_relative_matrix, break_table, broken_sum, break_iter)==0) return 0;
					
					dam_break_order++;
					if(broken_sum==0)
					{
						hout << "    No bond was broken." << endl;
						flag_break = false;
					}
					else 
					{
						//-------------------------------------------------------------------------------------------------------
						//Clear ele_self_matrix and ele_relative_matrix
						//MathMatrix mat_zero(24,24);
						//for(int i=0; i<(int)Mesh.elements.size(); i++)
						//{
						//	ele_self_matrix[i] = mat_zero;
						//	for(int j=0; j<(int)Mesh.elements[i].relative_eles.size(); j++) ele_relative_matrix[i][j] = mat_zero;
						//}

						//Gp_val.clear();
						//if(Glosmat.Gen_element_matrices_damage(Init->gauss.num, Init->peri_para, Init->mod_disc.mod, Init->weight_func, Init->damages, Material.mats_vec, 
						//																			Mesh.nodes, Mesh.elements, full_dam_eles, damage_table, break_table, Gp_val, ele_self_matrix, ele_relative_matrix)==0) return 0;

						hout << "    !!B " << broken_sum << " bonds were broken." << endl;
					}
				}

				ct1 = clock();
				hout << "    " << "Ement matrices update in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
				hout << "^_^ Element matrices updated successfully!" << endl << endl;
				cout << "^_^ Element matrices updated successfully!" << endl << endl;
				//-------------------------------------------------------------------------------------------------------
				count_iter++;
			}

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Updating weighting function: alpha, and updating full damaged elements 
			if(count_iter!=1||(Init->iter.max_iter==1&&(flag_break||flag_dam)))
			{
				cout << "-_- Updating alpha weighting function and full_damaged_eles..." << endl;
				hout << "-_- Updating alpha weighting function and full_damaged_eles..." << endl;

				WeightFunc Wei_Fun;
				if(Wei_Fun.Update_weighting_function_alpha(Mesh.nodes, Mesh.elements, Init->peri_para, Init->damages, dam_iter, damage_table, full_dam_eles, Init->weight_func, alpha_change)==0) return 0;

				if(alpha_change)
				{
					hout << "    Weighting functions alpha and full_damaged_eles are updated." << endl;
					hout << "^_^ Weighting functions alpha and full_damaged_eles are updated successfully!" << endl << endl;
					cout << "^_^ Weighting functions alpha and full_damaged_eles are updated successfully!" << endl << endl;
				}
				else
				{
					hout << "    Weighting functions alpha and full_damaged_eles did not change." << endl;
					hout << "^_^ Weighting functions alpha and full_damaged_eles tested successfully!" << endl << endl;
					cout << "^_^ Weighting functions alpha and full_damaged_eles tested successfully!" << endl << endl;
				}
			}

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Record dissipative energy of every element
			for(int i=0; i<(int)Mesh.elements.size(); i++) dissip_en[i] = Mesh.elements[i].Dissipative_energy;

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Record values of every step
			if(!alpha_change||!full_damaged_change)
			{
				count_ramp++;
				U_tot.push_back(U_Solution);
				Dissipation.push_back(dissip_en);
				nodes_tot.push_back(Mesh.nodes);
				elements_tot.push_back(Mesh.elements);

				//---------------------------------------------------------------------------
				//Output data
				Gauss gau;		//gaussian potins
				if(gau.Generate_gauss(Init->gauss.num)==0) return 0;

				int utsize = (int)U_tot.size()-1;
				
				//Output weighting function contour
				stringstream owfc;
				if(utsize<10)	owfc << "./Results/Weight_Func_Contour_000" << utsize << ".dat";
				else if (utsize<100)	owfc << "./Results/Weight_Func_Contour_00" << utsize << ".dat";
				else if (utsize<1000)	owfc << "./Results/Weight_Func_Contour_0" << utsize << ".dat";
				else	owfc << "./Results/Weight_Func_Contour_" << utsize << ".dat";

				WeightFunc Wei_Fun;
				Wei_Fun.Output_weight_function_contour(owfc.str(), Mesh.nodes, Mesh.elements, gau.gauss, gau.weight, Init->weight_func);

				//Output damage contour
				stringstream odac;
				if(utsize<10)	odac << "./Results/Damage_Contour_000" << utsize << ".dat";
				else if (utsize<100)	odac << "./Results/Damage_Contour_00" << utsize << ".dat";
				else if (utsize<1000)	odac << "./Results/Damage_Contour_0" << utsize << ".dat";
				else	odac << "./Results/Damage_Contour_" << utsize << ".dat";
				Output_damage_contour(odac.str(), Mesh.nodes, Mesh.elements, gau.gauss, gau.weight, damage_table);


				//Immediately output result data
				Postprocessor Output(U_Solution, dissip_en, Mesh.nodes, Mesh.elements);
				if(Output.Immedi_write_data((int)U_tot.size()-1, Init->rw_mod.type)==0) return 0;
			}
		}

		//-----------------------------------------------------------------------------------------------------------------------------------------
		time_t it1 = time(NULL); //Standard time
		cout << endl;
		cout << "----------------------------------------------------------------"<<endl;
		cout << "Iterative algorithm finished in " << totalcount <<" iterations."<<endl;
		cout << "It took about " << (int)(it1 - it0) << "secs."  << endl;
		cout << "Iterative algorithm done successfully!" <<endl;
		cout << endl;
		hout << endl;
		hout << "----------------------------------------------------------------"<<endl;
		hout << "Iterative algorithm finished in " << totalcount <<" iterations."<<endl;
		hout << "It took about " << (int)(it1 - it0) << "secs."  << endl;
		hout << "Iterative algorithm done successfully!" <<endl;
		hout << endl;

	}
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Postprocess
	ct0 = clock();
	cout << endl;
	cout << "-------------------------------------------------" << endl;
	cout << "|               Postprocessing                  |" << endl;
	cout << "-------------------------------------------------" << endl;
	cout << endl;
	cout << "-_- Starting postprocessor..." << endl;
	cout << endl;
	hout << endl;
	hout << "-------------------------------------------------" << endl;
	hout << "|               Postprocessing                  |" << endl;
	hout << "-------------------------------------------------" << endl;
	hout << endl;
	hout << "-_- Starting postprocessor..." << endl;
	hout << endl;
	
	//-------------------------------------------------------------------------------------------------------
	Postprocessor Post;
	if(Post.Treatment(Init, nodes_tot, elements_tot, Material.mats_vec, U_tot, Init->rw_mod.type, Dissipation)==0) return 0; //执行后处理过程	
	
	//-------------------------------------------------------------------------------------------------------
	ct1 = clock();
	cout << endl;
	cout << "----------------------------------------------------------------"<<endl;
	cout << "    Postprocessor took about " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs."  << endl;
	cout << "^_^ Postprocessor done successfully!" <<endl;
	cout << endl;
	hout << endl;
	hout << "----------------------------------------------------------------"<<endl;
	hout << "    Postprocessor took about " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs."  << endl;
	hout << "^_^ Postprocessor done successfully!" <<endl;
	hout << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------
//Output the damage values at every node while updating
void App_Damage::Output_damage_contour(const string &output_file_name, vector<Node> &nodes, const vector<Element> &elements, const vector<Node> &gauss,
																	const vector<double> &weight, const vector<vector<double> > &damage_table)const
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
	vector<double> damage_val(elements.size(), 0);
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

			damage_val[i] += damage_table[i][count]*Jac_val*weight[count];
		}
		damage_val[i] = damage_val[i]/ele_vol;
	}
	//Average value of weighting function at every node
	vector<double> nod_val(nodes.size(), 0);
	for(int i=0; i<(int)nodes.size(); i++)
	{
		const int nres = (int)nodes[i].relative_eles.size();
		for(int j=0; j<nres; j++)
			nod_val[i] += damage_val[nodes[i].relative_eles[j]];
		nod_val[i] = nod_val[i]/nres;
	}

	//Clear relative elements info of every node
	for(int i=0; i<(int)nodes.size(); i++)	nodes[i].relative_eles.clear(); 

	//---------------------------------------------------------------------------
	//Outpur resutls
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Damage_Contour" << endl;
	otec << "VARIABLES = X, Y, Z, Damage" << endl;
	
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
//Output the damage values at every node while updating for 2D problem
void App_Damage::Output_damage_contour_2D(const string &output_file_name, vector<Node> &nodes, vector<Element> &elements, const vector<Node> &gauss,
																		   const vector<double> &weight, const vector<vector<double> > &damage_table)const
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
	double max_dam_value = 0.0;
	Mesher Mesh;
	//Estimate the weighting value of every elements
	vector<double> damage_val(elements.size(), 0);
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
			if(max_dam_value<damage_table[i][count]) max_dam_value = damage_table[i][count];
			
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

			damage_val[i] += damage_table[i][count]*Jac_val*weight[count];
		}
		damage_val[i] = damage_val[i]/ele_area;

		elements[i].dam = damage_val[i];
	}
	//Average value of weighting function at every node
	vector<double> nod_val(nodes.size(), 0);
	for(int i=0; i<(int)nodes.size(); i++)
	{
		const int nres = (int)nodes[i].relative_eles.size();
		for(int j=0; j<nres; j++)
			nod_val[i] += damage_val[nodes[i].relative_eles[j]];
		nod_val[i] = nod_val[i]/nres;
	}

	//Clear relative elements info of every node
	for(int i=0; i<(int)nodes.size(); i++)	nodes[i].relative_eles.clear(); 

	//---------------------------------------------------------------------------
	//Outpur resutls
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Damage_Contour" << endl;
	otec << "VARIABLES = X, Y, Damage" << endl;
	
	otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size() << ", F=FEPOINT, ET=QUADRILATERAL" << endl;
	for(int i=0; i<(int)nodes.size(); i++)	
		otec << nodes[i].x << "  " << nodes[i].y << "  " << nod_val[i] << endl;
	otec << endl;
	for (int i=0; i < (int)elements.size(); i++)
	{
		otec	 << elements[i].nodes_id[0]+1 << "  " << elements[i].nodes_id[1]+1 << "  " 
				 << elements[i].nodes_id[2]+1 << "  " << elements[i].nodes_id[3]+1 << endl;
	}
	otec.close(); 

	//---------------------------------------------------------------------------
	//Output the maximum damage value
cout << "***Maximu damage value = " << max_dam_value << endl;
hout << "***Maximu damage value = " << max_dam_value << endl;
}
//===========================================================================
