//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	App_Fracture.cpp
//OBJECTIVE:	An adaptive algorithm between local and non-local models for the simulation of static fracture problems.
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================
#include "App_Fracture.h"

int App_Fracture::Application_fracture(Input *Init)const
{
	clock_t ct0,ct1; //time markers for cpu time of local operations

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Generate material database
	ct0 = clock();
	cout << "-_- Creating material database..." << endl;
	hout << "-_- Creating material database..." << endl;

	MatBase Material;
	if(Material.Generate_matbase(Init->stif_nonloc, Init->nonloc_gsize, Init->nonloc_gau.num, Init->peri_para)==0) return 0;
	
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

	//-------------------------------------------------------------------------------------------------------
	//Input weighting function value
	//cout << "-_- Totally updating alpha weighting function values..." << endl;
	//hout << "-_- Totally updating alpha weighting function values..." << endl;

	//WeightFunc Wei_Fun_In;
	//Wei_Fun_In.Input_weight_function_value("Total_weighting_function_value.dat", Init->weight_func);

	//cout << "-_- Alpha weighting function values updated successfully!" << endl;
	//hout << "-_- Alpha weighting function values updated successfully!" << endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//if "Read" mod, jump over Solving
	if(Init->rw_mod.type!="Read"&&Init->rw_mod.type=="Write")
	{
		//-----------------------------------------------------------------------------------------------------------------------------------------
		//Define the class of Global_Stiff_Matrix for the element matrices and total matrix
		Global_Stiff_Matrix Glosmat;
		//Define a table list of the broken bonds
		vector<bool> break_table;
		//Damaged element Id records
		vector<bool> dam_ele_id;
		//Energy dissipation of every element
		vector <double> dissip_en;

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
			if(Mesh.Generate_mesh(Init->grid_size, Init->ele_prop, Init->cracks, Init->geom_rve, Init->mod_disc.disc, Init->weight_func)== 0) return 0;  //Generate background regular mesh

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
				//Initialization of  damged element id
				dam_ele_id.assign(Mesh.elements.size(), false);
				//Initialization of  dissipative energy of every element
				dissip_en.assign(Mesh.elements.size(), 0.0);

				//-----------------------------------------------------------------------------------------------------------------------------------------
				//Initialization of  break_table
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
	
				break_table.assign(total_bonds, false);

				ct1 = clock();
				hout << "    The break_table initialized in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
				hout << "^_^ The break_table initialized successfully!" << endl << endl;
				cout << "^_^ The break_table initialized successfully!" << endl << endl;
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
			if(Glosmat.Gen_element_matrices(Init->gauss.num, Init->peri_para, Init->mod_disc.mod, Init->weight_func, Material.mats_vec, 
																			Mesh.nodes, Mesh.elements, break_table, Gp_val, ele_self_matrix, ele_relative_matrix)==0) return 0;

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
			bool alpha_change=false;
			int count_iter = 0;
			vector<double> U_Solution;
			vector<bool> dam_iter(Mesh.elements.size(), false);
			while (flag_break&&count_iter<Init->iter.max_iter)
			{
				hout <<"-_- iteration "<<count_iter<<endl;
				hout <<"----------------"<<endl;
				hout<<endl;
				cout <<"iteration "<<count_iter<<endl;
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

				int broken_sum=0;
				if(Glosmat.Update_nonlocal_element_matrices(Init->gauss.num, Init->peri_para, Init->mod_disc.mod, Init->weight_func, Material.mats_vec, Mesh.nodes,
																									Gp_val, U_Solution, Mesh.elements, ele_self_matrix, ele_relative_matrix, break_table, broken_sum, dam_iter)==0) return 0;

				ct1 = clock();
				hout << "    " << broken_sum << " bonds were broken, element matrices update in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
				hout << "^_^ Element matrices updated successfully!" << endl << endl;
				cout << "^_^ Element matrices updated successfully!" << endl << endl;
				//-------------------------------------------------------------------------------------------------------
			
				if(broken_sum==0) flag_break = false;
				count_iter++;
			}

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Updating weighting function: alpha
			if(count_iter!=1||(flag_break&&Init->iter.max_iter==1))
			{
				ct0 = clock();
				cout << "-_- Updating alpha weighting function..." << endl;
				hout << "-_- Updating alpha weighting function..." << endl;

				WeightFunc Wei_Fun;
				if(Wei_Fun.Update_weighting_function_alpha(Mesh.nodes, Mesh.elements, Init->peri_para, dam_iter, dam_ele_id, Init->weight_func, alpha_change)==0) return 0;

				ct1 = clock();
				if(alpha_change)
				{
					hout << "    Weighting functions alpha update in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
					hout << "^_^ Weighting functions alpha updated successfully!" << endl << endl;
					cout << "^_^ Weighting functions alpha updated successfully!" << endl << endl;
				}
				else
				{
					hout << "    Weighting functions alpha did not change, test finished in " << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "secs." << endl;
					hout << "^_^ Weighting functions alpha tested successfully!" << endl << endl;
					cout << "^_^ Weighting functions alpha tested successfully!" << endl << endl;
				}
			}

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Record dissipative energy of every element
			for(int i=0; i<(int)Mesh.elements.size(); i++) dissip_en[i] = Mesh.elements[i].Dissipative_energy;

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//Record values of every step
			if(!alpha_change)
			{
				count_ramp++;
				U_tot.push_back(U_Solution);
				Dissipation.push_back(dissip_en);
				nodes_tot.push_back(Mesh.nodes);
				elements_tot.push_back(Mesh.elements);
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
	//Output weighting function value
	cout << "-_- Outputing all weighting function values..." << endl;
	hout << "-_- Outputing all weighting function values..." << endl;

	WeightFunc Wei_Fun_Out;
	Wei_Fun_Out.Output_weight_function_value("Total_weighting_function_value.dat", Init->weight_func);

	cout << "-_- All weighting function values output successfully!" << endl;
	hout << "-_- All weighting function values output successfully!" << endl;

	//-------------------------------------------------------------------------------------------------------
	Postprocessor Post;
	if(Post.Treatment(nodes_tot, elements_tot, Material.mats_vec, U_tot, Init->rw_mod.type, Dissipation)==0) return 0; //执行后处理过程	
	
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
//===========================================================================
