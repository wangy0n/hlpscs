//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	App_Fracture.h
//OBJECTIVE:	An adaptive algorithm between local and non-local models for the simulation of static fracture problems.
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================

#ifndef APP_FRACTURE_H
#define APP_FRACTURE_H

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
class App_Fracture
{
	public:
		//Data Member

		//Constructor
		App_Fracture(){};

		//Member Functions
		int Application_fracture(Input *Init)const;

	private:
};
//------------------------------------------------
#endif
//===========================================================================
