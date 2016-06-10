//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Input_Reader.cpp
//OBJECTIVE:	Reading the input data
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================
extern int myid;
extern int group_size;
#include"mpi.h"
#include "Input_Reader.h"

//---------------------------------------------------------------------------
//Read data
int Input::Read_Infile(ifstream &infile)
{
  if (myid==0)
    {
	cout << "Reading input file..." << endl;
	hout << "Reading input file..." << endl;
    }
	while(!infile.eof())
	{
		istringstream istr(Get_Line(infile));
		if(infile.eof()) break;
		string str_temp;
		istr >> str_temp;

		if(str_temp.empty()) continue;  //skip the space or Enter key after every input item
		else if(str_temp=="Application_Name") { if(Read_application_name(app_name, infile)==0) return 0; }
		else if(str_temp=="Simulation_Name")	{ if(Read_simulation_name(simu_name, infile)==0) return 0; }
		else if(str_temp=="Local_Stiffness")	{ if(Read_local_stiffness(stif_loc, infile)==0) return 0; }
		else if(str_temp=="Local_Stiffness_2D")	{ if(Read_local_stiffness_2D(stif_loc, infile)==0) return 0; }
		else if(str_temp=="Nonlocal_Stiffness"){ if(Read_nonlocal_stiffness(stif_nonloc, infile)==0) return 0; }
		else if(str_temp=="Nonlocal_Stiffness_2D"){ if(Read_nonlocal_stiffness_2D(stif_nonloc, infile)==0) return 0; }
		else if(str_temp=="Peridynamic_Parameters")	{ if(Read_peridynamic_parameters(peri_para, infile)==0) return 0; }
		else if(str_temp=="Nonlocal_Grid_Size")	{ if(Read_grid_size(nonloc_gsize, infile)==0) return 0; }
		else if(str_temp=="Nonlocal_Gauss")	{ if(Read_gauss(nonloc_gau, infile)==0) return 0; }
		else if(str_temp=="RVE_Geometry")	{ if(Read_rve_geometry(geom_rve, infile)==0) return 0; }
		else if(str_temp=="Grid_Size")	{ if(Read_grid_size(grid_size, infile)==0) return 0; }
		else if(str_temp=="Weight_Function")	{ if(Read_weight_function(weight_func, infile)==0) return 0; }
		else if(str_temp=="Element_Material_Properties")	{ if(Read_element_properties(ele_prop, infile)==0) return 0; }
		else if(str_temp=="Crack")	{ if(Read_crack(cracks, infile)==0) return 0; }
		else if(str_temp=="Damage") { if(Read_damage(damages, infile)==0) return 0; }
		else if(str_temp=="Model_Discretization")	{ if(Read_model_discret(mod_disc, infile)==0) return 0; }
		else if(str_temp=="Iterative")	{ if(Read_iterative(iter, infile)==0) return 0; }
		else if(str_temp=="Read_Write")		{ if(Read_rw_mod(rw_mod, infile)==0) return 0; }
		else if(str_temp=="Gauss") { if(Read_gauss(gauss, infile)==0) return 0; }
		else if(str_temp=="Load")	{ if(Read_load(load, infile)==0) return 0; }
		else if(str_temp=="Load_2D")	{ if(Read_load_2D(load, infile)==0) return 0; }
		else if(str_temp=="Displacement")	{ if(Read_displacement(displace, infile)==0) return 0; }
		else if(str_temp=="Displacement_2D")	{ if(Read_displacement_2D(displace, infile)==0) return 0; }
		else if(str_temp=="Force_Disp_TPB_2D")  { if(Read_force_disp_TPB_2D(force_disp, geom_rve, iter.ramp_para, infile)==0) return 0; }
		else 
		  { if(myid==0)
		      {
                        cout << "Error: the keywords \"" << str_temp << "\" is not defined!" << endl; 
			hout << "Error: the keywords \"" << str_temp << "\" is not defined!" << endl; 
		      }
			return 0; 
		}

	}
	if(myid==0)
	  {
	cout << "Reading the keywords is finished!" << endl;
	hout << "Reading the keywords is finished!" << endl;
	  }
	if(!app_name.mark) { if(myid==0){cout << "Attention: \"Application_Name\" will use default parameters!" << endl; hout << "Attention: \"Application_Name\" will use default parameters!" << endl;} }
	if(!simu_name.mark) {if(myid==0){ cout << "Attention: \"Simulation_Name\" will use default parameters!" << endl; hout << "Attention: \"Simulation_Name\" will use default parameters!" << endl; }}
	if(!stif_loc.mark)	{ if(myid==0){cout << "Attention: \"Local_Stiffness\" will use default parameters!" << endl; hout << "Attention: \"Local_Stiffness\" will use default parameters!" << endl; }}
	if(!stif_nonloc.mark) {	if(myid==0){cout << "Attention: \"Nonlocal_Stiffness\" will use default parameters!" << endl; hout << "Attention: \"Nonlocal_Stiffness\" will use default parameters!" << endl; }}
	if(!peri_para.mark) {if(myid==0){  cout << "Attention: \"Peridynamic_Parameters\" will use default parameters!" << endl; hout << "Attention: \"Peridynamic_Parameters\" will use default parameters!" << endl; }}
	  if(!nonloc_gsize.mark) {  if(myid==0){cout << "Attention: \"Nonlocal_Grid_Size\" will use default parameters!" << endl; hout << "Attention: \"Nonlocal_Grid_Size\" will use default parameters!" << endl; }}
	  if(!nonloc_gau.mark) {	if(myid==0){cout << "Attention: \"Nonloca_Gauss\" will use default parameters!" << endl; hout << "Attention: \"Nonloca_Gauss\" will use default parameters!" << endl; }}
	  if(!geom_rve.mark) {if(myid==0){ cout << "Attention: \"RVE_Geometry\" will use default parameters!" << endl; hout << "Attention: \"RVE_Geometry\" will use default parameters!" << endl; }}
	  if(!grid_size.mark) {if(myid==0){ cout << "Attention: \"Grid_Size\" will use default parameters!" << endl; hout << "Attention: \"Grid_Size\" will use default parameters!" << endl; }}
	  if(!weight_func.mark) {if(myid==0){ cout << "Attention: \"Weight_Function\" will use default parameters!" << endl; hout << "Attention: \"Weight_Function\" will use default parameters!" << endl; }}
	  if(!ele_prop.mark) {if(myid==0){ cout << "Attention: \"Element_Material_Properties\" will use default parameters!" << endl; hout << "Attention: \"Element_Material_Properties\" will use default parameters!" << endl; }}
	  if(!cracks.mark) {if(myid==0){	cout << "Attention: \"Crack\" will use default parameters!" << endl; hout << "Attention: \"Crack\" will use default parameters!" << endl; }}
	  if(!damages.mark) {if(myid==0){	cout << "Attention: \"Damage\" will use default parameters!" << endl; hout << "Attention: \"Damage\" will use default parameters!" << endl; }}
	  if(!mod_disc.mark) {if(myid==0){ cout << "Attention: \"Model_Discretization\" will use default parameters!" << endl; hout << "Attention: \"Model_Discretization\" will use default parameters!" << endl; }}
	  if(!iter.mark) {if(myid==0){	cout << "Attention: \"Iterative\" will use default parameters!" << endl; hout << "Attention: \"Iterative\" will use default parameters!" << endl; }}
	  if(!rw_mod.mark) {if(myid==0){ cout << "Attention: \"Read_Write\" will use default parameters!" << endl; hout << "Attention: \"Read_Write\" will use default parameters!" << endl; }}
	  if(!gauss.mark) {if(myid==0){ cout << "Attention: \"Gauss\" will use default parameters!" << endl; hout << "Attention: \"Gauss\" will use default parameters!" << endl; }}
	  if(!load.mark) {if(myid==0){	cout << "Attention: \"Load\" will use default parameters!" << endl; hout << "Attention: \"Load\" will use default parameters!" << endl; }}
	  if(!displace.mark) {if(myid==0){ cout << "Attention: \"Displacement\" will use default parameters!" << endl; hout << "Attention: \"Displacement\" will use default parameters!" << endl; }}
	  if(!force_disp.mark) { if(myid==0){cout << "Attention: \"Force_Disp_TPB (Three Point Bending)\" will use default parameters!" << endl; hout << "Attention: \"Force_Disp_TPB (Three Point Bending)\" will use default parameters!" << endl; }}

	return 1;
}
//---------------------------------------------------------------------------
//Initialize data
int Input::Data_Initialization()
{
	//Initialize name of simulation
	app_name.keywords = "Application_Name";
	app_name.mark = false;
	app_name.str = "App_Fracture";

	//Initialize name of simulation
	simu_name.keywords = "Simulation_Name";
	simu_name.mark = false;
	simu_name.str = "Test";

	//Initialize material stiffness of local continuum model
	stif_loc.keywords = "Local_Stiffness";
	stif_loc.mark = false;
	stif_loc.type.push_back("Isotropic");
	stif_loc.E11.push_back(100);
	stif_loc.E22.push_back(100);
	stif_loc.E33.push_back(100);
	stif_loc.Nu12.push_back(0.25);
	stif_loc.Nu23.push_back(0.25);
	stif_loc.Nu13.push_back(0.25);
	stif_loc.G12.push_back(40);
	stif_loc.G23.push_back(40);
	stif_loc.G13.push_back(40);

	//Initialize effective stiffness of non-local continuum model
	stif_nonloc.keywords = "Nonlocal_Stiffness";
	stif_nonloc.mark = false;
	stif_nonloc.type = "Isotropic";
	stif_nonloc.E11 = stif_nonloc.E22 = stif_nonloc.E33 = 100;
	stif_nonloc.Nu12 = stif_nonloc.Nu23 = stif_nonloc.Nu13 = 0.25;
	stif_nonloc.G12 = stif_nonloc.G23 = stif_nonloc.G13 = 40;	

	//Initialize peridynamic(non-local continuum) parameters
	peri_para.keywords = "Peridynamic_Parameters";
	peri_para.mark = false;
	peri_para.horizon_R = 1.0;
	peri_para.intrinsic_L = 0.01;
	peri_para.broken_factor = 1.1;

	//Initialize grid's size for computing coefficients of non-local effective stiffness
	nonloc_gsize.keywords = "Nonlocal_Grid_Size";
	nonloc_gsize.mark = false;
	nonloc_gsize.delta_x = 0.1;
	nonloc_gsize.delta_y = 0.1;
	nonloc_gsize.delta_z = 0.1;

	//Initialize the number gauss points for computing coefficients of non-local effective stiffness
	nonloc_gau.keywords = "Nonloca_Gauss";
	nonloc_gau.mark = false;
	nonloc_gau.num = 4;

	//Initialize the geometry of the RVE
	geom_rve.keywords = "RVE_Geometry";
	geom_rve.mark = false;
	geom_rve.origin.x = 0.0;
	geom_rve.origin.y = 0.0;
	geom_rve.origin.z = 0.0;
	geom_rve.origin.flag = 0;
	geom_rve.len_x = 1.0;
	geom_rve.wid_y = 1.0;
	geom_rve.hei_z = 1.0;
	geom_rve.volume = geom_rve.len_x*geom_rve.wid_y*geom_rve.hei_z;

	//Initialize grid's size
	grid_size.keywords = "Grid_Size";
	grid_size.mark = false;
	grid_size.delta_x = 0.1;
	grid_size.delta_y = 0.1;
	grid_size.delta_z = 0.1;

	//Initialize weighting function
	weight_func.keywords = "Weight_Function";
	weight_func.mark = false;
	weight_func.num = 0;

	//Initialize material properties of elements
	ele_prop.keywords = "Element_Material_Properties";
	ele_prop.mark = false;
	ele_prop.type = "Pure";
	ele_prop.radius = 0.0;
	ele_prop.zone_xmin = 0.0;
	ele_prop.zone_xmax = 0.0;
	ele_prop.zone_ymin = 0.0;
	ele_prop.zone_ymax = 0.0;
	ele_prop.zone_zmin = 0.0;
	ele_prop.zone_zmax = 0.0;

	//Initialize crack information
	cracks.keywords = "Crack";
	cracks.mark = false;
	cracks.num = 0;

	//Initialize damage parameters
	damages.keywords = "Damge";
	damages.mark = false;
	damages.d_crit = 0;
	damages.k1 = 1.0;
	damages.k0 = 0.0;

	//Initialize the types of model and  discretization
	mod_disc.keywords = "Model_Discretization";
	mod_disc.mark = false;
	mod_disc.mod = "Local";
	mod_disc.disc = "FEM";

	//Initialize numerical iterative number
	iter.keywords = "Iterative_Number";
	iter.mark = false;
	iter.type = false;
	iter.max_iter = 1;
	iter.ramp_para = 1;

	//Initialize read or write mode for data backup
	rw_mod.keywords = "Read_Write";
	rw_mod.mark = false;
	rw_mod.type = "Write";

	//Initialize the number of gauss points
	gauss.keywords = "Gauss";
	gauss.mark = false;
	gauss.num = 4;

	//Initialize the load condition
	load.keywords = "Load";
	load.mark = false;
	load.num = 0;

	//Initialize the displacement condition
	displace.keywords = "Displacement";
	displace.mark = false;
	displace.num = 0;

	//Initialize the Force_Disp_TPB condition for three point bending
	force_disp.keywords = "Force_Disp_TPB";
	force_disp.mark = false;
	force_disp.cx0 = 0.0;
	force_disp.cx1 = 0.0;
	force_disp.cy0 = 0.0;
	force_disp.cy1 = 0.0;
	force_disp.czz = 0.0;
	force_disp.delta_disp = 0.0;
if (myid==0)
    {
	cout << "^_^ Data initialization achieves" <<endl<<endl;
	hout << "^_^ Data initialization achieves" <<endl<<endl;
    }
	return 1;
}
//---------------------------------------------------------------------------
//Reading the name of application case
int Input::Read_application_name(struct App_name &app_name, ifstream &infile)
{
	if(app_name.mark)
	{
if (myid==0)
    {
		cout << "Attention: \"" << app_name.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << app_name.keywords << "\" has been input!" << endl;
    }
		return 0;
	}
	else app_name.mark = true;

	istringstream istr(Get_Line(infile));
	istr >> app_name.str;

	return 1;
}
//---------------------------------------------------------------------------
//Reading the name of simulation
int Input::Read_simulation_name(struct Simu_name &simu_name, ifstream &infile)
{
	if(simu_name.mark)
	  {
         if(myid==0)
	      {
		cout << "Attention: \"" << simu_name.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << simu_name.keywords << "\" has been input!" << endl;
	      }
		return 0;
	}
	else simu_name.mark = true;

	istringstream istr(Get_Line(infile));
	istr >> simu_name.str;

	return 1;
}
//---------------------------------------------------------------------------
//Reading material stiffness of local continuum model
int Input::Read_local_stiffness(struct Stif_loc &stif_loc, ifstream &infile)
{
	if(stif_loc.mark)
	  {
           if(myid==0)
	      {
		cout << "Attention: \"" << stif_loc.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << stif_loc.keywords << "\" has been input!" << endl;
	      }
		return 0;
	}
	else stif_loc.mark = true;

	istringstream istr0(Get_Line(infile));
	istr0 >> stif_loc.num;

	for(int i=0; i<stif_loc.num; i++)
	{
		istringstream istr(Get_Line(infile));
		string type_temp;
		istr >> type_temp;
		stif_loc.type.push_back(type_temp);

		if(type_temp=="Isotropic")
		{
			double E, Nu;
			istr >> E >> Nu;

			stif_loc.E11.push_back(E);
			stif_loc.E22.push_back(E);
			stif_loc.E33.push_back(E);
			stif_loc.Nu12.push_back(Nu);
			stif_loc.Nu23.push_back(Nu);
			stif_loc.Nu13.push_back(Nu);
			stif_loc.G12.push_back(0.5*E/(1.0+Nu));
			stif_loc.G23.push_back(0.5*E/(1.0+Nu));
			stif_loc.G13.push_back(0.5*E/(1.0+Nu));
		}
		else if(type_temp=="Transverse")
		{
			double E1, E2, Nu1, Nu2, G2;
			istr >> E1 >> E2 >> Nu1 >> Nu2 >> G2;

			stif_loc.E11.push_back(E1);
			stif_loc.E22.push_back(E1);
			stif_loc.E33.push_back(E2);
			stif_loc.Nu12.push_back(Nu1); 
			stif_loc.Nu23.push_back(Nu2);
			stif_loc.Nu13.push_back(Nu2);
			stif_loc.G12.push_back(0.5*E1/(1.0+Nu1)); 
			stif_loc.G23.push_back(G2);
			stif_loc.G13.push_back(G2);
		}
		else if(type_temp=="Orthotropic")
		{
			double E11, E22, E33, Nu12, Nu23, Nu13, G12, G23, G13;

			istr >> E11 >> E22 >> E33;
			istr >> Nu12 >> Nu23 >> Nu13;
			istr >> G12 >> G23 >> G13;

			stif_loc.E11.push_back(E11);
			stif_loc.E22.push_back(E22);
			stif_loc.E33.push_back(E33);
			stif_loc.Nu12.push_back(Nu12); 
			stif_loc.Nu23.push_back(Nu23);
			stif_loc.Nu13.push_back(Nu13);
			stif_loc.G12.push_back(G12); 
			stif_loc.G23.push_back(G23);
			stif_loc.G13.push_back(G13);
		}
		else
		  {if(myid==0)
		      {
			hout << "Error: the material type of local model " << type_temp << " is not defined!" << endl;
			cout << "Error: the material type of local model " << type_temp << " is not defined!" << endl;
		      }
			return 0;
		}
	}

	return 1;
}
//---------------------------------------------------------------------------
//Reading material stiffness of local continuum model for 2D problem
int Input::Read_local_stiffness_2D(struct Stif_loc &stif_loc, ifstream &infile)
{
	if(stif_loc.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << stif_loc.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << stif_loc.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else stif_loc.mark = true;

	istringstream istr0(Get_Line(infile));
	istr0 >> stif_loc.num;

	for(int i=0; i<stif_loc.num; i++)
	{
		istringstream istr(Get_Line(infile));
		string type_temp;
		istr >> type_temp;
		stif_loc.type.push_back(type_temp);

		if(type_temp=="Isotropic")
		{
			double E, Nu;
			istr >> E >> Nu;

			stif_loc.E11.push_back(E);
			stif_loc.E22.push_back(E);
			stif_loc.E33.push_back(E);
			stif_loc.Nu12.push_back(Nu);
			stif_loc.Nu23.push_back(Nu);
			stif_loc.Nu13.push_back(Nu);
			stif_loc.G12.push_back(0.5*E/(1.0+Nu));
			stif_loc.G23.push_back(0.5*E/(1.0+Nu));
			stif_loc.G13.push_back(0.5*E/(1.0+Nu));
		}
		else if(type_temp=="Orthotropic")
		{
			double E11, E22, Nu12, G12;

			istr >> E11 >> E22 >> Nu12 >> G12;

			stif_loc.E11.push_back(E11);
			stif_loc.E22.push_back(E22);
			stif_loc.E33.push_back(1.0E+30);
			stif_loc.Nu12.push_back(Nu12); 
			stif_loc.Nu23.push_back(0.0);
			stif_loc.Nu13.push_back(0.0);
			stif_loc.G12.push_back(G12); 
			stif_loc.G23.push_back(1.0E+30);
			stif_loc.G13.push_back(1.0E+30);
		}
		else
		{if(myid==0)
		      {
			hout << "Error: the material type of local model " << type_temp << " is not defined!" << endl;
			cout << "Error: the material type of local model " << type_temp << " is not defined!" << endl;
		      }
			return 0;
		}
	}

	return 1;
}
//---------------------------------------------------------------------------
//Reading effective stiffness of non-local continuum model
int Input::Read_nonlocal_stiffness(struct Stif_nonloc &stif_nonloc, ifstream &infile)
{
	if(stif_nonloc.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << stif_nonloc.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << stif_nonloc.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else stif_nonloc.mark = true;

	istringstream istr(Get_Line(infile));
	istr >> stif_nonloc.type;

	double E11, E22, E33, Nu12, Nu23, Nu13;
	if(stif_nonloc.type=="Isotropic")
	{
		double E;
		istr >> E;

		E11 = E22 = E33 = E;
		Nu12 = Nu23 = Nu13 = 0.25;
	}
	else if(stif_nonloc.type=="Transverse")
	{
				
		double E1, E2, Nu2;
		istr >> E1 >> E2 >> Nu2;

		E11 = E22 = E1; E33 = E2;
		Nu23 = Nu2; Nu13 = Nu2;   //注意: Nu31=Nu13*E33/E11
		Nu12 = (E1-4*Nu2*Nu2*E2)/(3*E1); 
	}
	else if(stif_nonloc.type=="Orthotropic")
	{
		istr >> E11 >> E22 >> E33 >> Nu12 >> Nu23 >> Nu13;
	}
	else
	{if(myid==0)
		      {
		cout << "Error: the material type of nonlocal model " << stif_nonloc.type << " is not defined!" << endl;
		hout << "Error: the material type of nonlocal model " << stif_nonloc.type << " is not defined!" << endl;
		      }
		return 0;
	}

	stif_nonloc.E11 = E11;
	stif_nonloc.E22 = E22;
	stif_nonloc.E33 = E33;
	stif_nonloc.Nu12 = Nu12;
	stif_nonloc.Nu23 = Nu23;
	stif_nonloc.Nu13 = Nu13;
	stif_nonloc.G12 = E11*E22*(Nu12*E22+Nu23*Nu13*E33)/(E11*E22-Nu12*Nu12*E22*E22-Nu23*Nu23*E11*E33-Nu13*Nu13*E22*E33-2*Nu12*Nu23*Nu13*E22*E33);
	stif_nonloc.G23 = E22*E33*(Nu23*E11+Nu12*Nu13*E22)/(E11*E22-Nu12*Nu12*E22*E22-Nu23*Nu23*E11*E33-Nu13*Nu13*E22*E33-2*Nu12*Nu23*Nu13*E22*E33);
	stif_nonloc.G13 = E11*E22*E33*(Nu13+Nu12*Nu23)/(E11*E22-Nu12*Nu12*E22*E22-Nu23*Nu23*E11*E33-Nu13*Nu13*E22*E33-2*Nu12*Nu23*Nu13*E22*E33);

	return 1;
}
//---------------------------------------------------------------------------
//Reading effective stiffness of non-local continuum model for 2D problem
int Input::Read_nonlocal_stiffness_2D(struct Stif_nonloc &stif_nonloc, ifstream &infile)
{
	if(stif_nonloc.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << stif_nonloc.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << stif_nonloc.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else stif_nonloc.mark = true;

	istringstream istr(Get_Line(infile));
	istr >> stif_nonloc.type;

	double E11, E22, E33, Nu12, Nu23, Nu13;
	if(stif_nonloc.type=="Isotropic")
	{
		double E;
		istr >> E;

		E11 = E22 = E33 = E;
		Nu12 = Nu23 = Nu13 = 0.33333333;
	}
	else if(stif_nonloc.type=="Orthotropic")
	{
		istr >> E11 >> E22 >> Nu12;
	}
	else
	{if(myid==0)
		      {
		cout << "Error: the material type of nonlocal model " << stif_nonloc.type << " is not defined!" << endl;
		hout << "Error: the material type of nonlocal model " << stif_nonloc.type << " is not defined!" << endl;
		      }
		return 0;
	}

	stif_nonloc.E11 = E11;
	stif_nonloc.E22 = E22;
	stif_nonloc.E33 = 1.0E+30;
	stif_nonloc.Nu12 = Nu12;
	stif_nonloc.Nu23 = 0.0;
	stif_nonloc.Nu13 = 0.0;
	stif_nonloc.G12 = Nu12*E11*E22/(E11-Nu12*Nu12*E22);
	stif_nonloc.G23 = 1.0E+30;
	stif_nonloc.G13 = 1.0E+30;

	return 1;
}
//---------------------------------------------------------------------------
//Reading peridynamic(non-local continuum) parameters
int Input::Read_peridynamic_parameters(struct Peri_para &peri_para, ifstream &infile)
{
	if(peri_para.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << peri_para.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << peri_para.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else peri_para.mark = true;

	istringstream istr(Get_Line(infile));
	istr >> peri_para.horizon_R >> peri_para.intrinsic_L  >> peri_para.broken_factor;
	if(peri_para.horizon_R<0||peri_para.intrinsic_L<=0||peri_para.broken_factor<=1.0) 
	{ if(myid==0)
		      {
		cout <<"Error: selected parameter in \"Peridynamic_Parameters\" is wrong!" << endl; 
		hout <<"Error: selected parameter in \"Peridynamic_Parameters\" is wrong!" << endl; 
		      }
		return 0; 
	}	

	return 1;
}
//---------------------------------------------------------------------------
//Reading geometric information of the RVE
int Input::Read_rve_geometry(struct Geom_RVE &geom_rve, ifstream &infile)
{
	if(geom_rve.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << geom_rve.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << geom_rve.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else geom_rve.mark = true;

	istringstream istr(Get_Line(infile));
	istr >> geom_rve.origin.x >> geom_rve.origin.y >> geom_rve.origin.z;
	istr >> geom_rve.len_x >> geom_rve.wid_y >> geom_rve.hei_z;
	if(geom_rve.len_x<0||geom_rve.wid_y<0||geom_rve.hei_z<0)
	{if(myid==0)
		      {
		cout << "Error: the sizes of RVE should be positive!" << endl;
		hout << "Error: the sizes of RVE should be positive!" << endl;
		      }
		return 0;
	} 	

	return 1;
}
//---------------------------------------------------------------------------
//Reading grid's size
int Input::Read_grid_size(struct Grid_size &grid_size, ifstream &infile)
{
	if(grid_size.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << grid_size.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << grid_size.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else grid_size.mark = true;

	istringstream istr(Get_Line(infile));
	istr >> grid_size.delta_x >> grid_size.delta_y >> grid_size.delta_z;
	if(grid_size.delta_x<0||grid_size.delta_y<0||grid_size.delta_z<0)
	{if(myid==0)
		      {
		cout << "Error: size of grid in \"" << grid_size.keywords << "\" is less than 0!" << endl;
		hout << "Error: size of grid in \"" << grid_size.keywords << "\" is less than 0!" << endl;
		      }
		return 0;
	}

	return 1;
}
//---------------------------------------------------------------------------
//Reading weighting function information
int Input::Read_weight_function(struct Weight_func &weight_func, ifstream &infile)
{
	if(weight_func.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << weight_func.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << weight_func.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else weight_func.mark = true;

	istringstream iwf(Get_Line(infile));
	iwf >> weight_func.num;
	
	for(int i=0; i<weight_func.num; i++)
	{
		istringstream istr(Get_Line(infile));
		string str_temp;
		istr >> str_temp;
		weight_func.shape.push_back(str_temp);
		if(str_temp=="Sphere")
		{
			Point_3D poi_temp;
			istr >> poi_temp.x >> poi_temp.y >> poi_temp.z;
			weight_func.center.push_back(poi_temp);
			
			double r0, r1;
			istr >> r0 >> r1;
			weight_func.r0.push_back(r0);
			weight_func.r1.push_back(r1);
			weight_func.ratio.push_back(1.0);		//fullfilled data
			if(r0>r1) 
			{if(myid==0)
		      {
				cout << "Error: the inner radius is larger than the outer radius of homocentric sphere!" << endl;
				hout << "Error: the inner radius is larger than the outer radius of homocentric sphere!" << endl;
		      }
				return 0; 
			}
		}
		else if(str_temp=="Cylinder_x"||str_temp=="Cylinder_y"||str_temp=="Cylinder_z")
		{
			Point_3D poi_temp;
			if(str_temp=="Cylinder_x") istr >> poi_temp.y >> poi_temp.z;
			else if(str_temp=="Cylinder_y") istr >> poi_temp.x >> poi_temp.z;
			else if(str_temp=="Cylinder_z") istr >> poi_temp.x >> poi_temp.y;
			weight_func.center.push_back(poi_temp);

			double r0, r1, ratio;
			istr >> r0 >> r1 >> ratio;
			weight_func.r0.push_back(r0);
			weight_func.r1.push_back(r1);
			weight_func.ratio.push_back(ratio);
			if(r0>r1) 
			{if(myid==0)
		      {
				cout << "Error: the inner radius is larger than the outer radius of homocentric ellipse in the cross section!" << endl;
				hout << "Error: the inner radius is larger than the outer radius of homocentric ellipse in the cross section!" << endl;
		      }
				return 0;
			}
			if(ratio<=0) 
			{if(myid==0)
		      {
				cout << "Error: the ratio of axes in ellipse of the cylinder section is under 0!" << endl;
				hout << "Error: the ratio of axes in ellipse of the cylinder section is under 0!" << endl;
		      }
				return 0; 
			}
		}
		else if(str_temp=="Null")
		{
			Point_3D poi_temp(0,0,0);
			weight_func.center.push_back(poi_temp);
			weight_func.r0.push_back(0.0);
			weight_func.r1.push_back(0.0);
			weight_func.ratio.push_back(1.0);
		}
		else
		{if(myid==0)
		      {
			cout << "Error: the shape of weighting area " << str_temp << " is not defined!" << endl;
			hout << "Error: the shape of weighting area " << str_temp << " is not defined!" << endl;
		      }
			return 0;
		}
		
		string str_func;
		istr >> str_func;
		weight_func.func_order.push_back(str_func);

		if(str_func=="Constant") 
		{
			double val_temp;
			istr >> val_temp;
			weight_func.func_constant.push_back(val_temp);
			if(val_temp>1.0||val_temp<0.0) 
			{if(myid==0)
		      {
				cout << "Error: the specified value of constant weighting function is not between 0 and 1!" << endl;
				hout << "Error: the specified value of constant weighting function is not between 0 and 1!" << endl;
		      }
				return 0;
			}
		}
		if(str_func=="Linear"||str_func=="Cubic")
		{
			weight_func.func_constant.push_back(0.0);    //keep a correspondence between  "func_order" and "func_constant".
		}
		else
		{ if(myid==0)
		      {
			cout << "Error, the weighting function order " << str_func << " is not yet defined!" << endl; 
			hout << "Error, the weighting function order " << str_func << " is not yet defined!" << endl; 
		      }
			return 0; 
		}
	}

	return 1;
}
//---------------------------------------------------------------------------
//Reading material properties of elements
int Input::Read_element_properties(struct Ele_prop &ele_prop, ifstream &infile)
{
	if(ele_prop.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << ele_prop.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << ele_prop.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else ele_prop.mark = true;

	istringstream istr(Get_Line(infile));
	istr >> ele_prop.type;

	if(ele_prop.type=="Fiber")
	{
		istr >> ele_prop.radius;
		if(ele_prop.radius<0.0) 
		{ if(myid==0)
		      {
			cout <<"Error: the radius of fiber is less than 0!" << endl;	
			hout <<"Error: the radius of fiber is less than 0!" << endl;
		      }	
			return 0; 
		}
	}
	else if(ele_prop.type=="Zone") 
	{
		istr >> ele_prop.zone_xmin >> ele_prop.zone_xmax >> ele_prop.zone_ymin >> ele_prop.zone_ymax >> ele_prop.zone_zmin >> ele_prop.zone_zmax;
		if(ele_prop.zone_xmin>ele_prop.zone_xmax||
			ele_prop.zone_ymin>ele_prop.zone_ymax||
			ele_prop.zone_zmin>ele_prop.zone_zmax)
		{ if(myid==0)
		      {
			cout <<"Error: the parameters of zone in \" Element Material Properties\" are defined incorrectly!" << endl;	
			hout <<"Error: the parameters of zone in \" Element Material Properties\" are defined incorrectly!" << endl;
		      }	
			return 0; 
		}
	}
	else if(ele_prop.type!="Pure")
	{if(myid==0)
		      {
		cout << "The type of element material properties is defined incorrectly!" << endl; 
		hout << "The type of element material properties is defined incorrectly!" << endl; 
		      }
		return 0;
	}

	return 1;
}
//---------------------------------------------------------------------------
//Reading crack information
int Input::Read_crack(struct Crack &cracks, ifstream &infile)
{
	if(cracks.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << cracks.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << cracks.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else cracks.mark = true;

	istringstream istr(Get_Line(infile));
	istr >> cracks.num;
	
	for(int i=0; i<cracks.num; i++)
	{
		int n0, n1, n2, n3;
		char ch0, ch1, ch2;
		double d10, d20, d21, d30, d31;
		istringstream istr0(Get_Line(infile));
		istr0 >> n0 >> ch0;
		if(n0!=0||(ch0!='x'&&ch0!='y'&&ch0!='z'))
		{if(myid==0)
		      {
			cout << "Error: The parameter in the line 0 of  the crack " << i << " is defined incorrectly!" << endl;
			hout << "Error: The parameter in the line 0 of  the crack " << i << " is defined incorrectly!" << endl;
		      }
			return 0;
		}

		istringstream istr1(Get_Line(infile));
		istr1 >> n1 >> ch1 >> d10;
		if(n1!=1||(ch1!='x'&&ch1!='y'&&ch1!='z'))
		{if(myid==0)
		      {
			cout << "Error: The parameter in the line 1 of  the crack " << i << " is defined incorrectly!" << endl;
			hout << "Error: The parameter in the line 1 of  the crack " << i << " is defined incorrectly!" << endl;
		      }
			return 0;
		}

		istringstream istr2(Get_Line(infile));
		istr2 >> n2 >> ch2 >> d20 >> d21;
		if(n2!=2||(ch2!='x'&&ch2!='y'&&ch2!='z')||d20>=d21)
		{if(myid==0)
		      {
			cout << "Error: The parameter in the line 2 of  the crack " << i << " is defined incorrectly!" << endl;
			hout << "Error: The parameter in the line 2 of  the crack " << i << " is defined incorrectly!" << endl;
		      }
			return 0;
		}
		if((ch0!='x'&&ch1!='x'&&ch2!='x')||
			(ch0!='y'&&ch1!='y'&&ch2!='y')||
			(ch0!='z'&&ch1!='z'&&ch2!='z')) 
		{if(myid==0)
		      {
			cout << "Error: 'x', 'y' or 'z' is repeated in crack representation!" << endl;
			hout << "Error: 'x', 'y' or 'z' is repeated in crack representation!" << endl;
		      }
			return 0;
		}

		istringstream istr3(Get_Line(infile));
		istr3 >> n3 >> d30 >> d31;
		if(n3!=3||d30<0||d30>1||d31<0||d31>2)		
		{if(myid==0)
		      {
			cout << "Error: The parameter in the line 3 of  the crack " << i << " is defined incorrectly!" << endl;
			hout << "Error: The parameter in the line 3 of  the crack " << i << " is defined incorrectly!" << endl;
		      }
			return 0;
		}

		cracks.n0.push_back(n0);			cracks.n1.push_back(n1);			cracks.n2.push_back(n2);			cracks.n3.push_back(n3);
		cracks.ch0.push_back(ch0);		cracks.ch1.push_back(ch1);		cracks.ch2.push_back(ch2);
		cracks.d10.push_back(d10);		cracks.d20.push_back(d20);		cracks.d21.push_back(d21);		cracks.d30.push_back(d30);		cracks.d31.push_back(d31);
	}

	return 1;
}
//---------------------------------------------------------------------------
//Reading damage parameters
int Input::Read_damage(struct Damage &damages, ifstream &infile)
{
	if(damages.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << damages.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << damages.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else damages.mark = true;

	istringstream istr(Get_Line(infile));
	istr >> damages.d_crit >> damages.k1 >> damages.k0;

	if(damages.d_crit<0||damages.d_crit>1.0)
	{if(myid==0)
		      {
		cout << "The critical parameter in the item \"Damage\" is defined incorrectly!" << endl; 
		hout << "The critical parameter in the item \"Damage\" is defined incorrectly!" << endl; 
		      }
		return 0;
	}
	if(fabs(damages.k1)<=Zero)
	{if(myid==0)
		      {
		cout << "The coefficient of Pseudo potential in the item \"Damage\" equals to zero!" << endl; 
		hout << "The coefficient of Pseudo potential in the item \"Damage\" equals to zero!" << endl; 
		      }
		return 0;
	}

	return 1;
}
//---------------------------------------------------------------------------
//Reading the types of model and discretization
int Input::Read_model_discret(struct Model_Discret &mod_disc, ifstream &infile)
{
	if(mod_disc.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << mod_disc.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << mod_disc.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else mod_disc.mark = true;

	istringstream istr(Get_Line(infile));
	istr >> mod_disc.mod >> mod_disc.disc;
	if(mod_disc.mod!="Local"&&mod_disc.mod!="Nonlocal"&&mod_disc.mod!="Hybrid")
	{if(myid==0)
		      {
		cout << "Error: the type of model \"" << mod_disc.mod << "\" is not yet defined!" << endl;
		hout << "Error: the type of model \"" << mod_disc.mod << "\" is not yet defined!" << endl;
		      }
		return 0;
	}
	if(mod_disc.disc!="FEM"&&mod_disc.disc!="DGFEM")
	{if(myid==0)
		      {
		cout << "Error: the type of discretization \"" << mod_disc.disc << "\" is not yet defined!" << endl;
		hout << "Error: the type of discretization \"" << mod_disc.disc << "\" is not yet defined!" << endl;
		      }
		return 0;
	}

	return 1;
}
//---------------------------------------------------------------------------
//Reading the numerical iterative number
int Input::Read_iterative(struct Iterative &iter, ifstream &infile)
{
	if(iter.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << iter.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << iter.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else iter.mark = true;

	istringstream istr(Get_Line(infile));
	istr >> iter.type;

	if(iter.type!=0&&iter.type!=1)
	{if(myid==0)
		      {
		cout << "Error: " << iter.type << "=" << iter.type << " is not 0 or 1!" << endl;
		hout << "Error: " << iter.type << "=" << iter.type << " is not 0 or 1!" << endl;
		      }
		return 0;
	}

	if(iter.type)	 
	{
		istr >> iter.max_iter >> iter.ramp_para;

		if(iter.max_iter<=0||iter.ramp_para<=0)
		{if(myid==0)
		      {
			cout << "The parameter in the item \"Iterative\" is defined incorrectly!" << endl; 
			hout << "The parameter in the item \"Iterative\" is defined incorrectly!" << endl; 
		      }
			return 0; 
		}
	}

	return 1;
}
//---------------------------------------------------------------------------
//Reading read or write mode for data backup
int Input::Read_rw_mod(struct RW_mod &rw_mod, ifstream &infile)
{
	if(rw_mod.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << rw_mod.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << rw_mod.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else rw_mod.mark = true;

	istringstream istr(Get_Line(infile));
	istr >> rw_mod.type;

	return 1;
}
//---------------------------------------------------------------------------
//Reading number of gauss points
int Input::Read_gauss(struct Gauss_Point &gauss, ifstream &infile)
{
	if(gauss.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << gauss.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << gauss.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else gauss.mark = true;

	istringstream istr(Get_Line(infile));
	istr >> gauss.num;
	if(gauss.num<=0||gauss.num>=100)
	{if(myid==0)
		      {
		cout <<"Error: the number of gauss point in 1D is " << gauss.num << "!" << endl; 
		hout <<"Error: the number of gauss point in 1D is " << gauss.num << "!" << endl; 
		      }
		return 0; 
	}

	return 1;
}
//---------------------------------------------------------------------------
//Reading the load condition
int Input::Read_load(struct Load &load, ifstream &infile)
{
	if(load.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << load.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << load.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else load.mark = true;

	istringstream istr0(Get_Line(infile));
	istr0 >> load.num;
	
	for(int i=0; i<load.num; i++)
	{
		string domain_type;
		istringstream istr(Get_Line(infile));
		istr >> domain_type;
		load.domain_type.push_back(domain_type);

		vector<double> coef;
		if(domain_type=="Point")	{ coef.resize(3); istr >> coef[0] >> coef[1] >> coef[2]; }
		else if(domain_type=="Surface")	{ coef.resize(4); istr >> coef[0] >> coef[1] >> coef[2] >> coef[3]; }
		else if(domain_type=="Zone")	{ coef.resize(6); istr >> coef[0] >> coef[1] >> coef[2] >> coef[3] >> coef[4] >> coef[5]; }
		else
		{if(myid==0)
		      {
			cout << "Error: selected type for loaded area is not defined!" << endl; 
			hout << "Error: selected type for loaded area is not defined!" << endl;
		      }
			return 0; 
		}

		//coefficients of loaded area
		load.coef.push_back(coef);
		
		vector<string> load_type;
		vector<double> value;
		int sign = 0;
		double val_temp;
		while(!istr.eof())
		{
			string str_temp;
			istr >> str_temp;
			if(str_temp.empty()) continue;
			if(str_temp=="Force_x")
			{
				if(sign>=1) 
				{if(myid==0)
		      {
					cout << "Error: the force type is repeated or in a wrong order!" << endl; 
					hout << "Error: the force type is repeated or in a wrong order!" << endl; 
		      }
					return 0;
				}
				sign += 1;
			}
			else if(str_temp=="Force_y")
			{
				if(sign>=2) 
				{if(myid==0)
		      {
					cout << "Error: the force type is repeated or in a wrong order!" << endl; 
					hout << "Error: the force type is repeated or in a wrong order!" << endl; 
		      }
					return 0;
				}
				sign += 2;
			}
			else if(str_temp=="Force_z")
			{
				if(sign>=4) 
				{if(myid==0)
		      {
					cout << "Error: the force type is repeated or in a wrong order!" << endl;
					hout << "Error: the force type is repeated or in a wrong order!" << endl;
		      }
					return 0; 
				}
				sign += 4;
			}
			else 
			{if(myid==0)
		      {
				cout << "Error: the force type is not defined!" << endl;
				hout << "Error: the force type is not defined!" << endl;
		      }
				return 0; 
			}

			load_type.push_back(str_temp);
			istr >> val_temp;
			value.push_back(val_temp);
		}
		if((int)value.size()==0||(int)value.size()>3||(int)value.size()!=(int)load_type.size())
		{if(myid==0)
		      {
			cout << "Error: force conditions are not input or excessive!" << endl; 
			hout << "Error: force conditions are not input or excessive!" << endl; 
		      }
			return 0;
		}

		load.load_type.push_back(load_type);
		load.value.push_back(value);
	}

	return 1;
}
//---------------------------------------------------------------------------
//Reading the load condition for 2D problem
int Input::Read_load_2D(struct Load &load, ifstream &infile)
{
	if(load.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << load.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << load.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else load.mark = true;

	istringstream istr0(Get_Line(infile));
	istr0 >> load.num;
	
	for(int i=0; i<load.num; i++)
	{
		string domain_type;
		istringstream istr(Get_Line(infile));
		istr >> domain_type;
		load.domain_type.push_back(domain_type);

		vector<double> coef;
		if(domain_type=="Point")	{ coef.resize(2); istr >> coef[0] >> coef[1]; }
		else if(domain_type=="Line")	{ coef.resize(4); istr >> coef[0] >> coef[1] >> coef[2] >> coef[3]; }
		else if(domain_type=="Zone")	{ coef.resize(4); istr >> coef[0] >> coef[1] >> coef[2] >> coef[3]; }
		else
		{if(myid==0)
		      {
			cout << "Error: selected type for loaded area is not defined!" << endl; 
			hout << "Error: selected type for loaded area is not defined!" << endl;
		      }
			return 0; 
		}

		//coefficients of loaded area
		load.coef.push_back(coef);
		
		vector<string> load_type;
		vector<double> value;
		int sign = 0;
		double val_temp;
		while(!istr.eof())
		{
			string str_temp;
			istr >> str_temp;
			if(str_temp.empty()) continue;
			if(str_temp=="Force_x")
			{
				if(sign>=1) 
				{if(myid==0)
		      {
					cout << "Error: the force type is repeated or in a wrong order!" << endl; 
					hout << "Error: the force type is repeated or in a wrong order!" << endl; 
		      }
					return 0;
				}
				sign += 1;
			}
			else if(str_temp=="Force_y")
			{
				if(sign>=2) 
				{if(myid==0)
		      {
					cout << "Error: the force type is repeated or in a wrong order!" << endl; 
					hout << "Error: the force type is repeated or in a wrong order!" << endl; 
		      }
					return 0;
				}
				sign += 2;
			}
			else 
			{if(myid==0)
		      {
				cout << "Error: the force type is not defined!" << endl;
				hout << "Error: the force type is not defined!" << endl;
		      }
				return 0; 
			}

			load_type.push_back(str_temp);
			istr >> val_temp;
			value.push_back(val_temp);
		}
		if((int)value.size()==0||(int)value.size()>2||(int)value.size()!=(int)load_type.size())
		{if(myid==0)
		      {
			cout << "Error: force conditions are not input or excessive!" << endl; 
			hout << "Error: force conditions are not input or excessive!" << endl; 
		      }
			return 0;
		}

		load.load_type.push_back(load_type);
		load.value.push_back(value);
	}

	return 1;
}
//---------------------------------------------------------------------------
//Reading the displacement condition
int Input::Read_displacement(struct Displace &displace, ifstream &infile)
{
	if(displace.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << displace.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << displace.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else displace.mark = true;

	istringstream istr0(Get_Line(infile));
	istr0 >> displace.num;
	
	for(int i=0; i<displace.num; i++)
	{
		string domain_type;
		istringstream istr(Get_Line(infile));
		istr >> domain_type;
		displace.domain_type.push_back(domain_type);

		vector<double> coef;
		if(domain_type=="Point")	{ coef.resize(3); istr >> coef[0] >> coef[1] >> coef[2]; }
		else if(domain_type=="Surface")	{ coef.resize(4); istr >> coef[0] >> coef[1] >> coef[2] >> coef[3]; }
		else if(domain_type=="Zone")	{ coef.resize(6); istr >> coef[0] >> coef[1] >> coef[2] >> coef[3] >> coef[4] >> coef[5]; }
		else if(domain_type!="All_surfaces")
		{if(myid==0)
		      {
			cout << "Error: selected type for displacement area is not defined!" << endl; 
			hout << "Error: selected type for displacement area is not defined!" << endl;
		      }
			return 0; 
		}

		//coefficients of displacement area
		displace.coef.push_back(coef);
		
		vector<string> disp_type;
		vector<double> value;
		int sign = 0;
		double val_temp;
		while(!istr.eof())
		{
			string str_temp;
			istr >> str_temp;
			if(str_temp.empty()) continue;
			if(str_temp=="Fixed_displacement_x")
			{
				if(sign>=1) 
				{ if(myid==0)
		      {
					cout << "Error: the displacement type is repeated or in a wrong order!" << endl;
					hout << "Error: the displacement type is repeated or in a wrong order!" << endl;
		      }
					return 0; 
				}
				sign += 1;
			}
			else if(str_temp=="Fixed_displacement_y")
			{
				if(sign>=2) 
				{ if(myid==0)
		      {
					cout << "Error: the displacement type is repeated or in a wrong order!" << endl; 
					hout << "Error: the displacement type is repeated or in a wrong order!" << endl; 
		      }
					return 0; 
				}
				sign += 2;
			}
			else if(str_temp=="Fixed_displacement_z")
			{
				if(sign>=4) 
				{ if(myid==0)
		      {
					cout << "Error: the displacement type is repeated or in a wrong order!" << endl; 
					hout << "Error: the displacement type is repeated or in a wrong order!" << endl;
		      } 
					return 0; 
				}
				sign += 4;
			}
			else if(str_temp=="Pure_shear") 
			{
				if(sign>0) 
				{ if(myid==0)
		      {
					cout << "Error: \"Pure_shear\" is not the only displacement condition!" << endl; 
					hout << "Error: \"Pure_shear\" is not the only displacement condition!" << endl;
		      } 
					return 0; 
				}
			}
			else 
			{if(myid==0)
		      {
				cout << "Error: selected type for displacement condition is not defined!" << endl; 
				hout << "Error: selected type for displacement condition is not defined!" << endl; 
		      }
				return 0; 
			}

			disp_type.push_back(str_temp);
			istr >> val_temp;
			value.push_back(val_temp);
		}
		if((int)value.size()==0||(int)value.size()>3||(int)value.size()!=(int)disp_type.size())
		{if(myid==0)
		      {
			cout << "Error: force conditions are not input or excessive!" << endl; 
			hout << "Error: force conditions are not input or excessive!" << endl;
		      } 
			return 0;
		}

		displace.disp_type.push_back(disp_type);
		displace.value.push_back(value);
	}

	return 1;
}
//---------------------------------------------------------------------------
//Reading the displacement condition for 2D problem
int Input::Read_displacement_2D(struct Displace &displace, ifstream &infile)
{
	if(displace.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << displace.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << displace.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else displace.mark = true;

	istringstream istr0(Get_Line(infile));
	istr0 >> displace.num;
	
	for(int i=0; i<displace.num; i++)
	{
		string domain_type;
		istringstream istr(Get_Line(infile));
		istr >> domain_type;
		displace.domain_type.push_back(domain_type);

		vector<double> coef;
		if(domain_type=="Point")	{ coef.resize(2); istr >> coef[0] >> coef[1]; }
		else if(domain_type=="Line")	{ coef.resize(4); istr >> coef[0] >> coef[1] >> coef[2] >> coef[3]; }
		else if(domain_type=="Zone")	{ coef.resize(4); istr >> coef[0] >> coef[1] >> coef[2] >> coef[3]; }
		else if(domain_type!="All_boundaries")
		{if(myid==0)
		      {
			cout << "Error: selected type for displacement area is not defined!" << endl; 
			hout << "Error: selected type for displacement area is not defined!" << endl;
		      }
			return 0; 
		}

		//coefficients of displacement area
		displace.coef.push_back(coef);
		
		vector<string> disp_type;
		vector<double> value;
		int sign = 0;
		double val_temp;
		while(!istr.eof())
		{
			string str_temp;
			istr >> str_temp;
			if(str_temp.empty()) continue;
			if(str_temp=="Fixed_displacement_x")
			{
				if(sign>=1) 
				{ if(myid==0)
		      {
					cout << "Error: the displacement type is repeated or in a wrong order!" << endl;
					hout << "Error: the displacement type is repeated or in a wrong order!" << endl;
		      }
					return 0; 
				}
				sign += 1;
			}
			else if(str_temp=="Fixed_displacement_y")
			{
				if(sign>=2) 
				{ if(myid==0)
		      {
					cout << "Error: the displacement type is repeated or in a wrong order!" << endl; 
					hout << "Error: the displacement type is repeated or in a wrong order!" << endl; 
		      }
					return 0; 
				}
				sign += 2;
			}
			else if(str_temp=="Pure_shear") 
			{
				if(sign>0) 
				{ if(myid==0)
		      {
					cout << "Error: \"Pure_shear\" is not the only displacement condition!" << endl; 
					hout << "Error: \"Pure_shear\" is not the only displacement condition!" << endl; 
		      }
					return 0; 
				}
			}
			else 
			{if(myid==0)
		      {
				cout << "Error: selected type for displacement condition is not defined!" << endl; 
				hout << "Error: selected type for displacement condition is not defined!" << endl; 
		      }
				return 0; 
			}

			disp_type.push_back(str_temp);
			istr >> val_temp;
			value.push_back(val_temp);
		}
		if((int)value.size()==0||(int)value.size()>2||(int)value.size()!=(int)disp_type.size())
		{if(myid==0)
		      {
			cout << "Error: force conditions are not input or excessive!" << endl; 
			hout << "Error: force conditions are not input or excessive!" << endl; 
		      }
			return 0;
		}

		displace.disp_type.push_back(disp_type);
		displace.value.push_back(value);
	}

	return 1;
}
//---------------------------------------------------------------------------
//Reading the parameters for the force-displacement plot of three point bending
int Input::Read_force_disp_TPB_2D(struct Force_Disp_TPB &force_disp, const struct Geom_RVE &geom_rve, const int &iter_num, ifstream &infile)
{
	if(force_disp.mark)
	{if(myid==0)
		      {
		cout << "Attention: \"" << force_disp.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << force_disp.keywords << "\" has been input!" << endl;
		      }
		return 0;
	}
	else force_disp.mark = true;

	double total_displacement;
	istringstream istr0(Get_Line(infile));
	istr0 >> force_disp.cx0 >> force_disp.cx1 >> force_disp.cy0 >> total_displacement;
	if(force_disp.cx0>force_disp.cx1||
		force_disp.cx0<geom_rve.origin.x||
		force_disp.cx1>geom_rve.origin.x+geom_rve.len_x||
		force_disp.cy0<geom_rve.origin.y||
		force_disp.cy0>geom_rve.origin.y+geom_rve.wid_y)
	{if(myid==0)
		      {
		cout << "Error: the position parameters for force displacement plot of three point bending are not right!" << endl; 
		hout << "Error: the position parameters for force displacement plot of three point bending are not right!" << endl; 
		      }
		return 0;
	}
	if(total_displacement<0) 
	{if(myid==0)
		      {
		cout << "Error: the total displacement for force displacement plot of three point bending is less than 0!(Here an absolute value is needed!)" << endl; 
		hout << "Error: the total displacement for force displacement plot of three point bending are not right!(Here an absolute value is needed!)" << endl; 
		      }
		return 0;
	}

	force_disp.delta_disp = total_displacement/iter_num;

	return 1;
}
//---------------------------------------------------------------------------
//读入一行信息，并跳过注释行（以"%"开头）；
string Input::Get_Line(ifstream &infile)const
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
