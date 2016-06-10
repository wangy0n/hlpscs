//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Global_Stiff_Matrix.h
//OBJECTIVE:	Compute and assemble the global stiffness matrix
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================

#ifndef Global_Stiff_Matrix_H
#define Global_Stiff_Matrix_H

#include"Input_Reader.h"
#include"Fem_3D.h"
#include"MatPro.h"
#include"Gauss.h"
#include"Geometry_3D.h"
#include "WeightFunc.h"
#include"time.h"
#include "Hns.h"
using namespace hns;
#include <omp.h>  //openmp

#define CHUNKSIZE 1		//defines the chunk size as 1 contiguous iteration
//---------------------------------------------------------------------------
class Global_Stiff_Matrix
{
	public:
		//Constructor
		Global_Stiff_Matrix(){};
		
		//Member functions

		//Estimate element matrices
		int Gen_element_matrices(const int &gaussnum, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func,
													   const vector<MatPro> &mats, vector<Node> &nodes, const vector<Element> &elements, const vector<bool> &break_table, 
													   vector<vector<double> > &Gp_val, vector<MathMatrix> &ele_self_matrix, vector<vector<MathMatrix> > &ele_relative_matrix)const;
		//Estimate element matrices with local damage
		int Gen_element_matrices_damage(const int &gaussnum, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func, const struct Damage &damages,
																	  const vector<MatPro> &mats, vector<Node> &nodes, const vector<Element> &elements, const vector<bool> &full_dam_eles, const vector<vector<double> > &damage_table,
																	  const vector<bool> &break_table, vector<vector<double> > &Gp_val, vector<MathMatrix> &ele_self_matrix, vector<vector<MathMatrix> > &ele_relative_matrix)const;
		//Estimate element matrices with local damage for 2D problem
		int Gen_element_matrices_damage_2D(const int &gaussnum, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func, const struct Damage &damages,
							const vector<MatPro> &mats, vector<Node> &nodes, const vector<Element> &elements, const vector<bool> &full_dam_eles, const vector<vector<double> > &damage_table,
							const vector<bool> &break_table, vector<vector<double> > &Gp_val, vector<vector<double> > &Ak_val, vector<vector<double> > &Ele_local_stiff, vector<MathMatrix> &ele_self_matrix, 
							vector<vector<MathMatrix> > &ele_relative_matrix)const;
		//Assemble to the global stiffness matrix
		int Update_global_matrix(const vector<MathMatrix> &ele_self_matrix, const vector<vector<MathMatrix> > &ele_relative_matrix, const vector<Element> &elements,
													 const vector<int> &Iz, const vector<int> &Ig, vector<double> &total_matrix)const;
		//Update element matrices with broken bonds
		int Update_nonlocal_element_matrices(const int &gaussnum, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func, const vector<MatPro> &mats,
																			 const vector<Node> &nodes, const vector<vector<double> > &Gp_val, const vector<double> &U_s, vector<Element> &elements,
																			 vector<MathMatrix> &ele_self_matrix, vector<vector<MathMatrix> > &ele_relative_matrix, vector<bool> &break_table, int &broken_sum, vector<bool> &break_iter)const;		
		//Update element matrices with broken bonds (after up to damage criterion)
		int Update_nonlocal_element_matrices(const int &gaussnum, const struct Damage damages, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func, const vector<MatPro> &mats,
																			 const vector<Node> &nodes, const vector<vector<double> > &Gp_val, const vector<double> &U_s, vector<Element> &elements, const vector<vector<double> > &damage_table,
																			 vector<MathMatrix> &ele_self_matrix, vector<vector<MathMatrix> > &ele_relative_matrix, vector<bool> &break_table, int &broken_sum, vector<bool> &break_iter)const;
		//Update element matrices with broken bonds (after up to damage criterion) for 2D problem
		int Update_nonlocal_element_matrices_2D(const int &gaussnum, const struct Damage damages, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func, const vector<MatPro> &mats,
																					const vector<Node> &nodes, const vector<vector<double> > &Gp_val, const vector<double> &U_s, vector<Element> &elements, const vector<vector<double> > &damage_table,
																					vector<bool> &full_dam_eles,  vector<MathMatrix> &ele_self_matrix, vector<vector<MathMatrix> > &ele_relative_matrix, vector<bool> &break_table, int &broken_sum, vector<bool> &break_iter)const;
		//Update element matrices with broken bonds (after up to damage criterion) for 2D problem (New Algorithm)
		int Update_nonlocal_matrices_2D(const int &gaussnum, const struct Damage damages, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func,
																	const vector<Node> &nodes, const vector<vector<double> > &Gp_val, const vector<vector<double> > &Ak_val, const vector<double> &U_s, vector<Element> &elements,
																	const vector<vector<double> > &damage_table, vector<MathMatrix> &ele_self_matrix, vector<vector<MathMatrix> > &ele_relative_matrix, vector<bool> *break_table, vector<double> &dissip_en,
																	int &broken_sum, vector<bool> &break_iter)const;
		//Update damaged element matrices
		int Update_damaged_element_matrices(const int &gaussnum, const struct Damage damages, const vector<MatPro> &mats, const vector<Node> &nodes, const vector<Element> &elements, const vector<double> &U_s,
																			 const vector<bool> &full_dam_eles, vector<MathMatrix> &ele_self_matrix, vector<vector<double> > &damage_table, int &damaged_sum, vector<bool> &dam_iter)const;
		//Update damaged element matrices for 2D
		int Update_damaged_element_matrices_2D(const int &gaussnum, const struct Damage damages, const vector<MatPro> &mats, const vector<Node> &nodes, const vector<Element> &elements, const vector<double> &U_s,
																					 const vector<bool> &full_dam_eles, vector<MathMatrix> &ele_self_matrix, vector<vector<double> > &damage_table, int &damaged_sum, vector<bool> &dam_iter)const;
		//Update damaged element matrices for 2D (New Algorithm)
		int Update_damaged_matrices_2D(const int &gaussnum, const struct Damage damages, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func, 
															const vector<MatPro> &mats, const vector<Node> &nodes, const vector<Element> &elements, const vector<double> &U_s, const vector<bool> &full_dam_eles,
															const vector<vector<double> > &Ak_val, const vector<vector<double> > &Ele_local_stiff, const vector<bool> *break_table, const vector<vector<double> > &Gp_val, 
															vector<MathMatrix> &ele_self_matrix, vector<vector<double> > *damage_table, vector<vector<double> > *Y_force, vector<double> &damage_en, int &damaged_sum, vector<bool> &dam_iter)const;
		//Calculate the total energy available for break in every element in Postprocessor
		int Available_total_break_energy(const int &gaussnum, const struct Peri_para &peri_para, const string &com_mod, const struct Weight_func &weight_func,
			                                              const vector<Node> &nodes, const vector<Element> &elements, vector<double> &dissip_en)const;
private:

		//Member functions

		//Calculate data of every gaussian point of element
		void Generate_element_gauss_data(const vector<Element> &elements, const vector<Node> &nodes, const string &com_mod, const struct Weight_func &weight_func, const vector<Node> &gauss, 
																	  const vector<double> &weight, double (*gauss_ns)[8], double (*gauss_nw)[8], double (*gauss_dif)[3][8], double &Jacobi,  double (*gauss_po)[3], double (*ele_cent)[4])const;
		//Calculate data of every gaussian point of element for 2D problem
		void Generate_element_gauss_data_2D(const vector<Element> &elements, const vector<Node> &nodes, const string &com_mod, const struct Weight_func &weight_func, const vector<Node> &gauss, 
																	  const vector<double> &weight, double (*gauss_ns)[4], double (*gauss_nw)[4], double (*gauss_dif)[2][4], double &Jacobi,  double (*gauss_po)[3], double (*ele_cent)[4])const;
		//Generate the element matrix based on local continuum model (contact forces)
		void Generate_Contacforce_Elestiff(double (*element_stiff_matrix)[24], const double *alpha_key, const double (*Long_range_stiffness)[6][6], const double (*elenodes)[3], 
																		const vector<MatPro> &mats, const int &flag, const double &Jacobi, const double (*gauss_dif)[3][8], const vector<double> &weight)const;
		//Generate the element matrix based on local continuum model (coupled part + damaged part)
		void Generate_Contacforce_Elestiff_Coupled_Damaged(double (*element_stiff_matrix)[24], const double *alpha_key, const double (*Long_range_stiffness)[6][6], const Damage &damages, const bool full_dam_i,
				const vector<double> &ele_dam_table, const double (*elenodes)[3], const vector<MatPro> &mats, const int &flag, const int &elemat, const double &Jacobi, const double (*gauss_dif)[3][8], const vector<double> &weight)const;
		//Generate the element matrix based on local continuum model (coupled part + damaged part) for 2D problem
		void Generate_Contacforce_Elestiff_Coupled_Damaged_2D(double(*element_stiff_matrix)[12], const double *alpha_key, const double(*Long_range_stiffness)[6][6], double(*temp_elas)[6][6], const Damage &damages, const bool full_dam_i,
				const vector<double> &ele_dam_table, const double (*elenodes)[2], const vector<MatPro> &mats, const int &flag, const int &elemat, const double &Jacobi, const double (*gauss_dif)[2][4], const vector<double> &weight)const;
		//Generate the element matrix based on non-local continuum model (long-range forces) without broken bonds
		void Generate_Longforce_Elestiff_Brokeless(double (*element_stiff_matrix1)[24], double (*element_stiff_matrix2)[24], double *alpha_key, double (*Long_range_stiffness)[6][6], const struct Weight_func &weight_func,
																						const int &flag_left, const int &flag_right, const double (*elenodes_left)[3], const double (*elenodes_right)[3], const double (*gauss_ns)[8], const double (*gauss_nw)[8],
																						const double (*gauss_dif)[3][8], const vector<double> &weight, const struct Peri_para &peri_para, const double &Jacobi, const double (*gauss_po)[3], const double elec_left[],
																						const double elec_right[], const vector<bool> &break_table, const long int &pos_bt, const vector<vector<double> > &Gp_val, const int &pri_ele, const int &ass_ele)const;
		//Generate the element matrix based on non-local continuum model (long-range forces) without broken bonds for 2D problem
		void Generate_Longforce_Elestiff_Brokeless_2D(double(*element_stiff_matrix1)[12], double(*element_stiff_matrix2)[12], double *alpha_key, double(*Long_range_stiffness)[6][6], const bool full_dam_i,
																						const int &flag_left, const int &flag_right, const double (*elenodes_left)[2], const double (*elenodes_right)[2], const double (*gauss_ns)[4], const double (*gauss_nw)[4],
																						const double (*gauss_dif)[2][4], const vector<double> &weight, const struct Peri_para &peri_para, const double &Jacobi, const double (*gauss_po)[3], const double elec_left[],
																						const double elec_right[], const vector<bool> &break_table, const long int &pos_bt, const vector<vector<double> > &Gp_val, const int &pri_ele, const int &ass_ele)const;
		//(Constant stress assumption in a horizon) generate the element matrix based on non-local continuum model (long-range forces) without broken bonds for 2D problem
		void Generate_Longforce_Elestiff_Brokeless_2D_Stress_Assump(double(*element_stiff_matrix1)[12], double(*element_stiff_matrix2)[12], double *alpha_key, double(*Long_range_stiffness)[6][6], const struct Weight_func &weight_func,
																											const int &flag_left, const int &flag_right, const double(*elenodes_left)[2], const double(*elenodes_right)[2], const double(*gauss_ns)[4], const double(*gauss_nw)[4],
																											const double(*gauss_dif)[2][4], const vector<double> &weight, const struct Peri_para &peri_para, const double &Jacobi, const double(*gauss_po)[3], const double elec_left[],
																											const double elec_right[], const vector<vector<double> > &damage_table, const vector<bool> &break_table, const long int &pos_bt, const vector<vector<double> > &Gp_val, 
																											const int &pri_ele, const int &ass_ele)const;
		//Generate the element matrix based on non-local continuum model (long-range forces) with broken bonds
		void Generate_Longforce_Elestiff_Breaks(double (*element_stiff_matrix1)[24], double (*element_stiff_matrix2)[24], double *alpha_key, vector<bool> &break_table, int &broken_num, double &Dissip_ener,
																				  const struct Weight_func &weight_func, const int &flag_left, const int &flag_right, const double (*elenodes_left)[3], const double(*n_ele_l)[3], const double (*elenodes_right)[3],
																				   const double(*n_ele_r)[3], const double (*gauss_ns)[8], const double (*gauss_nw)[8], const double (*gauss_dif)[3][8], const vector<double> &weight, const struct Peri_para &peri_para,
																				   const double &Jacobi, const double (*gauss_po)[3], const double elec_left[], const double elec_right[], const long int &pos_bt, const vector<vector<double> > &Gp_val, const int &pri_ele,
																				   const int &ass_ele)const;
		//Generate the element matrix based on non-local continuum model (long-range forces) with broken bonds (after up to damage criterion)
		void Generate_Longforce_Elestiff_Breaks(double (*element_stiff_matrix1)[24], double (*element_stiff_matrix2)[24], double *alpha_key, vector<bool> &break_table, int &broken_num, const struct Damage damages,  const vector<vector<double> > &damage_table,
																				   double &Dissip_ener, const struct Weight_func &weight_func, const int &flag_left, const int &flag_right, const double (*elenodes_left)[3], const double(*n_ele_l)[3], const double (*elenodes_right)[3],
																				   const double(*n_ele_r)[3], const double (*gauss_ns)[8], const double (*gauss_nw)[8], const double (*gauss_dif)[3][8], const vector<double> &weight, const struct Peri_para &peri_para,
																				   const double &Jacobi, const double (*gauss_po)[3], const double elec_left[], const double elec_right[], const long int &pos_bt, const vector<vector<double> > &Gp_val, const int &pri_ele,
																				   const int &ass_ele)const;
		//Generate the element matrix based on non-local continuum model (long-range forces) with broken bonds (after up to damage criterion) for 2D problem
		void Generate_Longforce_Elestiff_Breaks_2D(double (*element_stiff_matrix1)[12], double (*element_stiff_matrix2)[12], double *alpha_key, vector<bool> &break_table, int &broken_num, const struct Damage damages,  const vector<vector<double> > &damage_table,
																				   double &Dissip_ener, const struct Weight_func &weight_func, const int &flag_left, const int &flag_right, const double (*elenodes_left)[2], const double(*n_ele_l)[2], const double (*elenodes_right)[2],
																				   const double(*n_ele_r)[2], const double (*gauss_ns)[4], const double (*gauss_nw)[4], const double (*gauss_dif)[2][4], const vector<double> &weight, const struct Peri_para &peri_para,
																				   const double &Jacobi, const double (*gauss_po)[3], const double elec_left[], const double elec_right[], const long int &pos_bt, const vector<vector<double> > &Gp_val, const int &pri_ele,
																				   const int &ass_ele)const;
		//Generate the element matrix based on non-local continuum model (long-range forces) with broken bonds (after up to damage criterion) for 2D problem (New Algorithm)
		void Longforce_Elestiff_Breaks_2D(double (*element_stiff_matrix1)[12], double (*element_stiff_matrix2)[12], const vector<bool> &break_table0, vector<bool> &break_table1, int &increased_num, int &broken_num,
																const struct Damage damages, const vector<vector<double> > &damage_table, double &Dissip_ener, const struct Weight_func &weight_func, const int &flag_left, const int &flag_right,
																const double (*elenodes_left)[2], const double(*n_ele_l)[2], const double (*elenodes_right)[2], const double(*n_ele_r)[2], const double (*gauss_ns)[4], const double (*gauss_nw)[4], 
																const double (*gauss_dif)[2][4], const vector<double> &weight, const struct Peri_para &peri_para,const double &Jacobi, const double (*gauss_po)[3], const double elec_left[], 
																const double elec_right[], const long int &pos_bt, const vector<vector<double> > &Gp_val, const vector<vector<double> > &Ak_val, const int &pri_ele,const int &ass_ele)const;

		//Add element matrix to globla matrix(in one element)
		void Add_to_gsmatrix(const double (*element_stiff_matrix)[24], const vector<int> &Iz, const vector<int> &Ig, vector<double> &total_matrix, const Element &element)const;
		void Add_to_gsmatrix(const MathMatrix &element_stiff_matrix, const vector<int> &Iz, const vector<int> &Ig, vector<double> &total_matrix, const Element &element)const;
		//Add element matrix to globla matrix(between different elements)
		void Add_to_gsmatrix(const double (*element_stiff_matrix)[24], const vector<int> &Iz, const vector<int> &Ig, vector<double> &total_matrix, const Element &ele_row, const Element &ele_col)const;
		void Add_to_gsmatrix(const MathMatrix &element_stiff_matrix, const vector<int> &Iz, const vector<int> &Ig, vector<double> &total_matrix, const Element &ele_row, const Element &ele_col)const;
		//Update the damaged brick element
		void Brick_line_update_damaged(double (*element_stiff_matrix)[24], vector<double> &ele_dam_table, int &damaged_num, const struct Damage damages, const double (*elenodes)[3],
																  const double *U_ele_nod, const vector<MatPro> &mats, const int &elemat, const vector<Node> &gauss, const vector<double> &wight)const;
		//Update the damaged quadrilateral element
		void Quadri_line_update_damaged(double (*element_stiff_matrix)[12], vector<double> &ele_dam_table, int &damaged_num, const struct Damage damages, const double (*elenodes)[2],
																  const double *U_ele_nod, const vector<MatPro> &mats, const int &elemat, const vector<Node> &gauss, const vector<double> &wight)const;
		//Update the damaged brick element (New Algorithm)
		void Quadri_damage_update(double (*element_stiff_matrix)[12], const vector<double> &Ele_stiff, const vector<double> Y_nonlocal_energy, vector<double> &Bond_dissip_en, 
						const vector<double> &ele_damtab0, const vector<double> &Yfors0, vector<double> &ele_damtab1, vector<double> &Yfors1, double &Ymax, double &Damage_ener, 
						int &increased_num, int &damaged_num, const struct Damage damages, const double (*elenodes)[2], const double *U_ele_nod, const vector<MatPro> &mats, const int &elemat, 
						const vector<Node> &gauss, const vector<double> &wight)const;
		//Update the damaged brick element (Last Algorithm)
		void Quadri_damage_update(double (*element_stiff_matrix)[12], const vector<double> &Ele_stiff, const vector<double> Y_nonlocal_energy, const vector<double> &ele_damtab0, 
						const vector<double> &Yfors0, vector<double> &ele_damtab1, vector<double> &Yfors1, vector<double> &ele_damtab2, vector<double> &Yfors2, double &Ymax, 
						double &Damage_ener, int &increased_num, int &damaged_num, const struct Damage damages, const double (*elenodes)[2], const double *U_ele_nod, const vector<MatPro> &mats, 
						const int &elemat, const vector<Node> &gauss, const vector<double> &wight)const;
		//Available break energy in every element
		void Available_Breaks_Energy_2D(double &Ava_dissip_ener, const int &flag_left, const int &flag_right, const double(*elenodes_left)[2], const double(*elenodes_right)[2], 
															const double(*gauss_ns)[4], const double(*gauss_nw)[4], const double(*gauss_dif)[2][4], const vector<double> &weight,
															const struct Peri_para &peri_para, const double &Jacobi, const double(*gauss_po)[3], const double elec_left[], const double elec_right[])const;
		//Generate the energy density at every gaussian point based on non-local continuum model judging broken bonds for 2D problem
		void Longforce_Energy_Density_Brokeless_2D(vector<double> &Y_nonlocal_energy, vector<double> &Bond_dissip_en, const double *U_ele_l, const double *U_ele_r, const int &flag_left, const int &flag_right, 
																				const double (*elenodes_left)[2], const double (*elenodes_right)[2], const double (*gauss_ns)[4], const double (*gauss_dif)[2][4], const vector<double> &weight, 
																				const struct Peri_para &peri_para, const double &Jacobi, const double (*gauss_po)[3], const double elec_left[], const double elec_right[], const vector<bool> *break_table,
																				const long int &pos_bt, const vector<vector<double> > &Gp_val, const int &pri_ele, const int &ass_ele)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
