//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	WeightFunc.h
//OBJECTIVE:	Estimate the values of points by weighting function
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================
#ifndef WEIGHTFUN_H
#define WEIGHTFUN_H

#include "Input_Reader.h"
#include "Fem_3D.h"
#include "Geometry_3D.h"
#include "Mesher.h"
#include "Hns.h"
using namespace hns;

//---------------------------------------------------------------------------
class WeightFunc
{
	public:
		//Data Member

		//Constructor
		WeightFunc(){};

		//Value weight function
		double Value_weight_function(const Point_3D &poi, const struct Weight_func &weight_func)const;

		//Export the values of weighting function at every node for testing
		void Export_weight_function_contour(const string &output_file_name, vector<Node> &nodes, const vector<Element> &elements, const double (*gauss_ns)[8], const double (*gauss_dif)[3][8], 
																		   const double &Jacobi, const double (*gauss_po)[3], const vector<double> &weight, const struct Weight_func &weight_func, const double (*ele_cent)[4])const;

		//Output the values of weighting function at every node while updating
		void Output_weight_function_contour(const string &output_file_name, vector<Node> &nodes, const vector<Element> &elements, const vector<Node> &gauss,
																		   const vector<double> &wight, const struct Weight_func &weight_func)const;
		//Output the values of weighting function at every node while updating for 2D problem
		void Output_weight_function_contour_2D(const string &output_file_name, vector<Node> &nodes, const vector<Element> &elements, vector<bool> &full_dam_eles,
																				  const vector<Node> &gauss, const vector<double> &weight, const struct Weight_func &weight_func)const;
		//Output the weighting function vectors for 2D problem
		int Output_weight_function_vectors(const string &output_file_name, const struct Weight_func &weight_func)const;

		//Output weighting function value
		void Output_weight_function_value(const string &output_file_name, const struct Weight_func &weight_func)const;
		//Intput weighting function value
		void Input_weight_function_value(const string &input_file_name, struct Weight_func &weight_func);

		//Calculate the values of weighting function of gaussian points in all elements
		void Value_weightfunc_gausspois_elements(const struct Weight_func &weight_func, const vector<Node> &nodes, const vector<Element> &elements, const int &GS,
																					 const double (*gauss_ns)[8], const double (*gauss_po)[3], double (*ele_cent)[4], vector<vector<double> > &Gp_val)const;

		//Calculate the values of weighting function of gaussian points in all elements for 2D problem
		void Value_weightfunc_gausspois_elements_2D(const struct Weight_func &weight_func, const vector<Node> &nodes, const vector<Element> &elements, const int &GS,
																					 const double (*gauss_ns)[4], const double (*gauss_po)[3], const double (*ele_cent)[4], vector<vector<double> > &Gp_val)const;

		//Update weighting function: alpha
		int Update_weighting_function_alpha(const vector<Node> &nodes, const vector<Element> &elements, const struct Peri_para &peri_para, const vector<bool> &break_iter,
																vector<bool> &dam_ele_id, struct Weight_func &weight_func, bool &alpha_change)const;
		//Update weighting function: alpha with damaged and broken bonds elements
		int Update_weighting_function_alpha(const vector<Node> &nodes, const vector<Element> &elements, const struct Peri_para &peri_para, const vector<bool> &break_iter,
																const vector<bool> &damage_iter, vector<bool> &dam_ele_id, struct Weight_func &weight_func, bool &alpha_change)const;
		//Update weighting function (and full_damaged_elements)
		int Update_weighting_function_alpha(const vector<Node> &nodes, const vector<Element> &elements, const struct Peri_para &peri_para, const struct Damage damages, const vector<bool> &dam_iter,
																 const vector<vector<double> > &damage_table,  vector<bool> &full_dam_eles, struct Weight_func &weight_func, bool &alpha_change)const;
		//Update full_damaged_elements and put weighting function inside completely damaged zone
		int Update_weighting_func_in_damaged_zone(const vector<Node> &nodes, const vector<Element> &elements, const struct Peri_para &peri_para, const struct Damage &damages, const vector<bool> &dam_iter,
																			  const vector<vector<double> > &damage_table,  vector<bool> &full_dam_eles, struct Weight_func &weight_func, vector<bool> &weighted_eles, bool &alpha_change)const;
		//This is a temporary test for pure damge model, at the beginning, dcirt=0.3 for example, when the full_damge_zone propagate to 5delt (for example), change dcrit=1.0 and continue pure damge model
		int Update_dcrit_pure_damage_test(const vector<Node> &nodes, const vector<Element> &elements, const struct Peri_para &peri_para, struct Damage &damages, const vector<bool> &dam_iter,
														const vector<vector<double> > &damage_table, vector<bool> &full_dam_eles, struct Weight_func &weight_func, vector<bool> &weighted_eles, bool &alpha_change);
	
	private:

};
//-------------------------------------------------------
#endif
//=====================================================