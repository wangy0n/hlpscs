//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Global_Load_Vector.h
//OBJECTIVE:	Estimate the global load vector
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#ifndef GLOBAL_LOAD_VECTOR_H
#define GLOBAL_LOAD_VECTOR_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<cmath>
#include<string>
#include"Input_Reader.h"
#include "Geometry_2D.h"
#include "Geometry_3D.h"

#include "Fem_3D.h"
#include "Hns.h"
using namespace hns;
//----------------------------------------------------------

class Global_Load_Vector
{
	public:
	
		//Constructor
		Global_Load_Vector(){};
		
		//Member functions

		//Generate global load vector
		int Gen_global_load_vector(const struct Load &load, const double &rmp, const vector<Node> &nodes, const vector<Element> &elements, vector<double> &equright)const;
		//Generate global load vector for 2D problem
		int Gen_global_load_vector_2D(const struct Load &load, const double &rmp, const vector<Node> &nodes, const vector<Element> &elements, vector<double> &equright)const;

	private:
};

//---------------------------------------------------------------------------
#endif
//========================================================================
