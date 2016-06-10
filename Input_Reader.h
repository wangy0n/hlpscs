//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Input_Reader.h
//OBJECTIVE:	Reading the input data
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================

#ifndef INPUTREADER_H
#define INPUTREADER_H

#include<iomanip>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include "time.h"
#include "Hns.h"
using namespace hns;

#include "Geometry_3D.h"

//---------------------------------------------------------------------------
//Name of application case
struct App_name{
			string keywords;
			bool mark;
			string str;
		};
//Name of simulation
struct Simu_name{
			string keywords;
			bool mark;
			string str;
		};
//Material stiffness of local continuum model
struct Stif_loc{
			string keywords;
			bool mark;
			int num;
			vector<string> type;
			vector<double> E11, E22, E33;
			vector<double> Nu12, Nu23, Nu13;
			vector<double> G12, G23, G13;
		};
//Effective stiffness of non-local continuum model
struct Stif_nonloc{
			string keywords;
			bool mark;
			string type;
			double E11, E22, E33;
			double Nu12, Nu23, Nu13;
			double G12, G23, G13;
		};
//Peridynamic(non-local continuum) parameters
struct Peri_para{
			string keywords;
			bool mark;
			double horizon_R;
			double intrinsic_L;
			double acoe[6];				//coefficients of peridynamic parameters to an effective stiffness 
			double broken_factor;  //stretch
		};
//The geometry of the RVE
struct Geom_RVE{
			string keywords;
			bool mark;
			Point_3D origin;
			double len_x, wid_y, hei_z;
			double volume;
		};
//Grid's size
struct Grid_size{
			string keywords;
			bool mark;
			double delta_x, delta_y, delta_z;
		};
//Weighting function
struct Weight_func{
			string keywords;
			bool mark;
			int num;		//the number of the weighting areas
			vector<string> shape;
			vector<Point_3D> center;
			vector<double> r0, r1, ratio;
			vector<string> func_order;
			vector<double> func_constant;
		};
//Material properties of elements
struct Ele_prop{
			string keywords;
			bool mark;
			string type;
			double radius;
			double zone_xmin;
			double zone_xmax;
			double zone_ymin;
			double zone_ymax;
			double zone_zmin;
			double zone_zmax;
		};
//Crack
struct Crack{
			string keywords;
			bool mark;
			int num;
			vector<int> n0, n1, n2, n3;
			vector<char> ch0, ch1, ch2;
			vector<double> d10, d20, d21, d30, d31;
		};
//Damage parameters
struct Damage{
			string keywords;
			bool mark;
			double d_crit;
			double k1;				//The parameter in the condition of Pseudo potential
			double k0;				//The parameter in the condition of Pseudo potential
		};
//The types of model and discretization
struct Model_Discret{
			string keywords;
			bool mark;
			string mod;
			string disc;
		};
//Numerical iterative number
struct Iterative{
			string keywords;
			bool mark;
			bool type;
			int max_iter;
			int ramp_para;
		};
//Read or write mode for data backup
struct RW_mod{
			string keywords;
			bool mark;
			string type;
		};
//The number of gauss points
struct Gauss_Point{
			string keywords;
			bool mark;
			int num;
		};
//The load condition
struct Load{
			string keywords;
			bool mark;
			int num;
			vector<string> domain_type;
			vector<vector<double> > coef;
			vector<vector<string> > load_type;
			vector<vector<double> > value;
		};
//The displacement condition
struct Displace{
			string keywords;
			bool mark;
			int num;
			vector<string> domain_type;
			vector<vector<double> > coef;
			vector<vector<string> > disp_type;
			vector<vector<double> > value;
		};
//The parameters for the force-displacement plot of three point bending
struct Force_Disp_TPB{
			string keywords;
			bool mark;
			double cx0, cx1;
			double cy0, cy1;
			double czz;
			double delta_disp;
		};
//---------------------------------------------------------------------------
class Input 
{
	public:
		//Data members
		struct App_name app_name;
		struct Simu_name simu_name;
		struct Stif_loc stif_loc;
		struct Stif_nonloc stif_nonloc;
		struct Peri_para peri_para;
		struct Geom_RVE geom_rve;
		struct Grid_size grid_size, nonloc_gsize;
		struct Weight_func weight_func;
		struct Ele_prop ele_prop;
		struct Crack cracks;
		struct Damage damages;
		struct Model_Discret mod_disc;
		struct Iterative iter;
		struct RW_mod rw_mod;
		struct Gauss_Point gauss, nonloc_gau;
		struct Load load;
		struct Displace displace;
		struct Force_Disp_TPB force_disp;

		//Constructor
		Input(){};  

		//Member functions
		int Data_Initialization();		//Initialize data
		int Read_Infile(ifstream &infile);		//Read data
		string Get_Line(ifstream &infile)const;									//读入信息一行，跳过注释行（以%开头）

private:
		//Member functions
		int Read_application_name(struct App_name &app_name, ifstream &infile);
		int Read_simulation_name(struct Simu_name &simu_name, ifstream &infile);
		int Read_local_stiffness(struct Stif_loc &stif_loc, ifstream &infile);
		int Read_local_stiffness_2D(struct Stif_loc &stif_loc, ifstream &infile);
		int Read_nonlocal_stiffness(struct Stif_nonloc &stif_nonloc, ifstream &infile);
		int Read_nonlocal_stiffness_2D(struct Stif_nonloc &stif_nonloc, ifstream &infile);
		int Read_peridynamic_parameters(struct Peri_para &peri_para, ifstream &infile);
		int Read_rve_geometry(struct Geom_RVE &geom_rve, ifstream &infile);
		int Read_grid_size(struct Grid_size &grid_size, ifstream &infile);
		int Read_weight_function(struct Weight_func &weight_func, ifstream &infile);
		int Read_element_properties(struct Ele_prop &ele_prop, ifstream &infile);
		int Read_crack(struct Crack &cracks, ifstream &infile);
		int Read_damage(struct Damage &damages, ifstream &infile);
		int Read_model_discret(struct Model_Discret &mod_disc, ifstream &infile);
		int Read_iterative(struct Iterative &iter, ifstream &infile);
		int Read_rw_mod(struct RW_mod &rw_mod, ifstream &infile);
		int Read_gauss(struct Gauss_Point &gauss, ifstream &infile);
		int Read_load(struct Load &load, ifstream &infile);
		int Read_load_2D(struct Load &load, ifstream &infile);
		int Read_displacement(struct Displace &displace, ifstream &infile);
		int Read_displacement_2D(struct Displace &displace, ifstream &infile);
		int Read_force_disp_TPB_2D(struct Force_Disp_TPB &force_disp, const struct Geom_RVE &geom_rve, const int &iter_num, ifstream &infile);
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
