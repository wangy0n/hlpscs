//============================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	App_Damage.h
//OBJECTIVE:	An adaptive algorithm between local and non-local models for the simulation of static damage and fracture propagation.
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//============================================================================================

#ifndef APP_DAMAGE_H
#define APP_DAMAGE_H

#include "time.h"
#include "Hns.h"
using namespace hns;

#include "Input_Reader.h"
#include "MatBase.h"
#include "SolveEqu.h"
#include "Global_Stiff_Matrix.h"
#include "Global_Load_Vector.h"
#include "Postprocessor.h"
#include "WeightFunc.h"

//---------------------------------------------------------------------------
class App_Damage
{
	public:
		//Data Member

		//Constructor
		App_Damage(){};

		//Member Functions
		int Application_damage(Input *Init)const;
		int Application_damage_2D(Input *Init)const;

	private:
		//Output the damage values at every node while updating
		void Output_damage_contour(const string &output_file_name, vector<Node> &nodes, const vector<Element> &elements, const vector<Node> &gauss,
													  const vector<double> &weight, const vector<vector<double> > &damage_table)const;
		//Output the damage values at every node while updating for 2D problem
		void Output_damage_contour_2D(const string &output_file_name, vector<Node> &nodes, vector<Element> &elements, const vector<Node> &gauss,
															 const vector<double> &weight, const vector<vector<double> > &damage_table)const;
};
//------------------------------------------------
#endif
//===========================================================================
