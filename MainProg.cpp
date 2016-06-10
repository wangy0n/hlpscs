//===========================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	MainProg.cpp
//OBJECTIVE:	Main program beginning
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//===========================================================================
#include "mpi.h"
#include <string>
#include "time.h"
#include "Hns.h"
using namespace hns;
int group_size;//The number of processes.
int myid;// The id of one process.
#include "Input_Reader.h"
#include "App_Fracture.h"
#include "App_Damage.h"
int main(int argc, char** argv)
{

  //-------------------------MPI initialize------------------------------------
  MPI_Init (&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&group_size);//get total number of the processes
  MPI_Comm_rank(MPI_COMM_WORLD,&myid); // get the id of one process

	
//-----------------------------------------------------------------------------------------------------------------------------------------
	//Read input file name into in_file
	string in_file;
	if(argc > 1)		in_file = argv[1];
	else
	{
	       in_file = "input.dat";
   if(myid==0)//The 1st if(myid==0).
	 {
		cout << "The input file name is:  ";		    
		cout << in_file << endl;
	  }  //The end of the 1st if(myid==0).
	};
	//Open the input file
	ifstream infile;
	infile.open(in_file.c_str());
	if(!infile) {
  if(myid==0)//The 2nd if(myid==0).
     {
 hout << "Failed to open input file: "  << in_file << endl; 
      }  //The end of the 2nd if(myid==0). 
       return 0;   
                     }

	//Read output file name into out_file
	string out_file;
	if(argc > 2)		out_file = argv[2];
	else
	{
                out_file = "output.dat";
  if(myid==0)//The 3rd if(myid==0).
     {
		cout << "The output file name is:  ";		     	       
		cout << out_file << endl;
      }  //The end of the 3rd if(myid==0)
	      
	};
	//Open the output stream
	if(out_file.size()>0) open_deffo_stream( (char*)out_file.c_str() );

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Identification Tag

     if(myid==0)//The 4th if(myid==0).
     {
	cout<<endl;
	cout<<"*************************************************"<<endl;
	cout<<"*                  NECA   v3.0                  *"<<endl;
	cout<<"*                                               *"<<endl;
	cout<<"*          Fei Han            Yan Azdoud        *"<<endl;
	cout<<"*                                               *"<<endl;
	cout<<"*     Propriety of KAUST/COHMAS Laboratory      *"<<endl;
	cout<<"* fei.han@kaust.edu.sa  yan.azdoud@kaust.edu.sa *"<<endl;
	cout<<"*************************************************"<<endl;
	cout<<endl;
	cout<<endl;
     }
     if(myid==0)
       {hout<<endl;
	hout<<endl;
	hout<<"*************************************************"<<endl;
	hout<<"*                  NECA   v3.0                  *"<<endl;
	hout<<"*                                               *"<<endl;
	hout<<"*          Fei Han            Yan Azdoud        *"<<endl;
	hout<<"*                                               *"<<endl;
	hout<<"*     Propriety of KAUST/COHMAS Laboratory      *"<<endl;
	hout<<"* fei.han@kaust.edu.sa  yan.azdoud@kaust.edu.sa *"<<endl;
	hout<<"*************************************************"<<endl;
	hout<<endl;
	hout<<endl;
     }
 //The end of the 4th if(myid==0). 
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Call for application cases

	//Time markers for total simulation
	time_t it_begin, it_end;
	it_begin = time(NULL);

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Input reader
   if(myid==0)//The 5th if(myid==0).
     {
	hout<<"======================================================" << endl;
	hout<<"-_- Input file reader......"<<endl;
	cout<<endl;
	cout<<"-------------------------------------------------"<<endl;
	cout<<"|                Data input                     |"<<endl;
	cout<<"-------------------------------------------------"<<endl;
	cout<<endl;
  }  
//The end of the 5th if(myid==0).  
	Input *Init = new Input;
	if(Init->Data_Initialization())
	{ 
		if(Init->Read_Infile(infile)==0) return 0;
	}
	else return 0;
	it_end= time(NULL);
   if(myid==0)//The 6th if(myid==0).
     {
	hout<<"    Operation done in "<<(int)(it_end-it_begin)<<"secs."<<endl;
	hout<<"^_^ Input achieves"<<endl<<endl;
	cout<<"^_^ Input achieves"<<endl<<endl;
	cout<<"Operation done in"<<(int)(it_end-it_begin)<<"secs."<<endl;//王增加的句子
  }  //The end of the 6th if(myid==0).	  
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Implementation
	if(Init->app_name.str=="App_Fracture")
	{
		//Definition
		App_Fracture *Compute =  new  App_Fracture;
		if(Compute->Application_fracture(Init)==0) return 0;
		delete Compute;
	}
	else if(Init->app_name.str=="App_Damage")
	{
		//Definition
		App_Damage *Compute =  new  App_Damage;
		if(Compute->Application_damage(Init)==0) return 0;
		delete Compute;
	}
	else if(Init->app_name.str=="App_Damage_2D")
	{
		//Definition
		App_Damage *Compute =  new  App_Damage;
		if(Compute->Application_damage_2D(Init)==0) return 0;
		delete Compute;
	}
    //-----------------------------------------------------------------------------------------------------------------------------------------
	//Delete the pointer of classes
	delete Init;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Time markers for total simulation
	it_end = time(NULL);
   if(myid==0)//The 7th if(myid==0).
     {
	cout<<endl;
	cout<<"*******************************************************************"<<endl;
	cout<<"    The simulation took "<<(int)(it_end-it_begin)<<"secs."<<endl;
	cout<<"^_^ End of simulation "<<endl;
	cout<<"*******************************************************************"<<endl;
	hout<<endl;
	hout<<"*******************************************************************"<<endl;
	hout<<"    The simulation took "<<(int)(it_end-it_begin)<<"secs."<<endl;
	hout<<"^_^ End of simulation "<<endl;
	hout<<"*******************************************************************"<<endl;
	  
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//End comments
	cout<< endl;
	cout<<"*************************************************"<<endl;
	cout<<"*             Hope it helped ;)                 *"<<endl;
	cout<<"*                                               *"<<endl;
	cout<<"*          Fei Han            Yan Azdoud        *"<<endl;
	cout<<"*                                               *"<<endl;
	cout<<"*     Propriety of KAUST/COHMAS Laboratory      *"<<endl;
	cout<<"* fei.han@kaust.edu.sa  yan.azdoud@kaust.edu.sa *"<<endl;
	cout<<"*************************************************"<<endl;
	cout<<endl;
	cout<<endl;

	hout<<endl;
	hout<<endl;
	hout<<"*************************************************"<<endl;
	hout<<"*             Hope it helped ;)                 *"<<endl;
	hout<<"*                                               *"<<endl;
	hout<<"*          Fei Han            Yan Azdoud        *"<<endl;
	hout<<"*                                               *"<<endl;
	hout<<"*     Propriety of KAUST/COHMAS Laboratory      *"<<endl;
	hout<<"* fei.han@kaust.edu.sa  yan.azdoud@kaust.edu.sa *"<<endl;
	hout<<"*************************************************"<<endl;
	hout<<endl;
	hout<<endl;
     }//The end of the 7th if(myid==0).
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Close the output stream
	close_deffo_stream();
	MPI_Finalize();
	return 1;

}
//===========================================================================
