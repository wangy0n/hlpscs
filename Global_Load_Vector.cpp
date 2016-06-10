//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Global_Load_Vector.cpp
//OBJECTIVE:	Estimate the global load vector
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#include "Global_Load_Vector.h"

//-------------------------------------------------------------------------------------------------------------------------------------------------------
//Generate global load vector
int Global_Load_Vector::Gen_global_load_vector(const struct Load &load, const double &rmp, const vector<Node> &nodes, const vector<Element> &elements, vector<double> &equright)const
{
	for(int i=0; i<load.num; i++)
	{
		//-----------------------------------------------------------
		//Define the value of displacement constraint by the ratio
		vector<double> value;
		for(int j=0; j<(int)load.value[i].size(); j++)  value.push_back(load.value[i][j]*rmp);

		//-----------------------------------------------------------
		if(load.domain_type[i]=="Point")
		{
			Point_3D fpoi(load.coef[i][0], load.coef[i][1], load.coef[i][2]);
			for(int j=0; j<(int)nodes.size(); j++)
			{
				if(fpoi.distance_to(nodes[j].x, nodes[j].y, nodes[j].z)<=Zero)
				{
					for(int k=0; k<(int)load.load_type[i].size(); k++)
					{
						if(load.load_type[i][k]=="Force_x") equright[3*j] += value[k];
						else if(load.load_type[i][k]=="Force_y") equright[3*j+1] += value[k];
						else if(load.load_type[i][k]=="Force_z") equright[3*j+2] += value[k];
					}
				}
			}
		}
		else if(load.domain_type[i]=="Surface")
		{
			Plane_3D fpla(load.coef[i]);
			for(int j=0; j<(int)elements.size(); j++)
			{
				vector<int> nodn;
				for(int k=0; k<(int)elements[j].nodes_id.size(); k++)
				{
					int enid =  elements[j].nodes_id[k];
					if(fpla.contain(nodes[enid].x, nodes[enid].y, nodes[enid].z)) nodn.push_back(enid);
				}
				if((int)nodn.size()==4)  //A surface of element in this plane //并没有排除是对角面的可能性
				{
					//calculate area (Wiki chinese)
					double sa=(nodes[nodn[0]].x-nodes[nodn[1]].x)*(nodes[nodn[0]].x-nodes[nodn[1]].x)
									+(nodes[nodn[0]].y-nodes[nodn[1]].y)*(nodes[nodn[0]].y-nodes[nodn[1]].y)
									+(nodes[nodn[0]].z-nodes[nodn[1]].z)*(nodes[nodn[0]].z-nodes[nodn[1]].z);

					double sb=(nodes[nodn[1]].x-nodes[nodn[2]].x)*(nodes[nodn[1]].x-nodes[nodn[2]].x)
									+(nodes[nodn[1]].y-nodes[nodn[2]].y)*(nodes[nodn[1]].y-nodes[nodn[2]].y)
									+(nodes[nodn[1]].z-nodes[nodn[2]].z)*(nodes[nodn[1]].z-nodes[nodn[2]].z);

					double sc=(nodes[nodn[2]].x-nodes[nodn[3]].x)*(nodes[nodn[2]].x-nodes[nodn[3]].x)
									+(nodes[nodn[2]].y-nodes[nodn[3]].y)*(nodes[nodn[2]].y-nodes[nodn[3]].y)
									+(nodes[nodn[2]].z-nodes[nodn[3]].z)*(nodes[nodn[2]].z-nodes[nodn[3]].z);
	
					double sd=(nodes[nodn[3]].x-nodes[nodn[0]].x)*(nodes[nodn[3]].x-nodes[nodn[0]].x)
									+(nodes[nodn[3]].y-nodes[nodn[0]].y)*(nodes[nodn[3]].y-nodes[nodn[0]].y)
									+(nodes[nodn[3]].z-nodes[nodn[0]].z)*(nodes[nodn[3]].z-nodes[nodn[0]].z);

					double se=(nodes[nodn[0]].x-nodes[nodn[2]].x)*(nodes[nodn[0]].x-nodes[nodn[2]].x)
									+(nodes[nodn[0]].y-nodes[nodn[2]].y)*(nodes[nodn[0]].y-nodes[nodn[2]].y)
									+(nodes[nodn[0]].z-nodes[nodn[2]].z)*(nodes[nodn[0]].z-nodes[nodn[2]].z);

					double sf=(nodes[nodn[1]].x-nodes[nodn[3]].x)*(nodes[nodn[1]].x-nodes[nodn[3]].x)
									+(nodes[nodn[1]].y-nodes[nodn[3]].y)*(nodes[nodn[1]].y-nodes[nodn[3]].y)
									+(nodes[nodn[1]].z-nodes[nodn[3]].z)*(nodes[nodn[1]].z-nodes[nodn[3]].z);
					
					double s_area = 0.25*sqrt(4*se*sf-(sb+sd-sa-sc)*(sb+sd-sa-sc));

					for(int m=0; m<(int)nodn.size(); m++)
					{
						for(int n=0; n<(int)load.load_type[i].size(); n++)
						{
							if(load.load_type[i][n]=="Force_x") equright[3*nodn[m]] += 0.25*value[n]*s_area;
							else if(load.load_type[i][n]=="Force_y") equright[3*nodn[m]+1] += 0.25*value[n]*s_area;
							else if(load.load_type[i][n]=="Force_z") equright[3*nodn[m]+2] += 0.25*value[n]*s_area;
						}
					}
				}
			}
		}
		else if(load.domain_type[i]=="Zone")
		{
			double	xmin	=	load.coef[i][0];
			double	xmax	=	load.coef[i][1];
			double	ymin	=	load.coef[i][2];
			double	ymax	=	load.coef[i][3];
			double	zmin	=	load.coef[i][4];
			double	zmax	=	load.coef[i][5];
			for(int j=0; j<(int)elements.size(); j++)
			{
				int ncount = 0;
				for(int k=0; k<(int)elements[j].nodes_id.size(); k++)
				{
					int enid =  elements[j].nodes_id[k];
					if((nodes[enid].x-xmin>=Zero)&&(xmax-nodes[enid].x>=Zero)&&
						(nodes[enid].y-ymin>=Zero)&&(ymax-nodes[enid].y>=Zero)&&
						(nodes[enid].z-zmin>=Zero)&&(zmax-nodes[enid].z>=Zero))
					{
						ncount++;
					}
				}
				if(ncount==8)  //An element in this zone
				{
					//calculate volume
					Hexahedron hex(elements[j].nodes_id[0], elements[j].nodes_id[1], elements[j].nodes_id[2], elements[j].nodes_id[3],
												 elements[j].nodes_id[4], elements[j].nodes_id[5], elements[j].nodes_id[6], elements[j].nodes_id[7]);
					double s_vol = hex.cal_hex_volume(nodes);
					for(int m=0; m<(int)elements[j].nodes_id.size(); m++)
					{
						for(int n=0; n<(int)load.load_type[i].size(); n++)
						{
							if(load.load_type[i][n]=="Force_x") equright[3*elements[j].nodes_id[m]] += 0.125*value[n]*s_vol;
							else if(load.load_type[i][n]=="Force_y") equright[3*elements[j].nodes_id[m]+1] += 0.125*value[n]*s_vol;
							else if(load.load_type[i][n]=="Force_z") equright[3*elements[j].nodes_id[m]+2] += 0.125*value[n]*s_vol;
						}
					}
				}
			}
		}
		else
		{
			cout << "Error: Failed to find the type of loaded domain!" << endl; 
			hout << "Error: Failed to find the type of loaded domain!" << endl;
			return 0; 
		}
	}

	return 1;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------
//Generate global load vector for 2D prolem
int Global_Load_Vector::Gen_global_load_vector_2D(const struct Load &load, const double &rmp, const vector<Node> &nodes, const vector<Element> &elements, vector<double> &equright)const
{
	for(int i=0; i<load.num; i++)
	{
		//-----------------------------------------------------------
		//Define the value of displacement constraint by the ratio
		vector<double> value;
//		for(int j=0; j<(int)load.value[i].size(); j++)  value.push_back(load.value[i][j]*rmp); 
		for (int j = 0; j<(int)load.value[i].size(); j++)  value.push_back(load.value[i][j]); //用于固定加载

		//-----------------------------------------------------------
		if(load.domain_type[i]=="Point")
		{
			Point_2D fpoi(load.coef[i][0], load.coef[i][1]);
			for(int j=0; j<(int)nodes.size(); j++)
			{
				if(fpoi.distance_to(nodes[j].x, nodes[j].y)<=Zero)
				{
					for(int k=0; k<(int)load.load_type[i].size(); k++)
					{
						if(load.load_type[i][k]=="Force_x") equright[3*j] += value[k];
						else if(load.load_type[i][k]=="Force_y") equright[3*j+1] += value[k];
					}
				}
			}
		}
		else if(load.domain_type[i]=="Line")  //注意这里的Line其实是线段的意思
		{
			vector<Side> sides;
			
			Point_2D poi[2];
			poi[0].x = load.coef[i][0];
			poi[0].y = load.coef[i][1];
			poi[1].x = load.coef[i][2];
			poi[1].y = load.coef[i][3];

			Line_2D line(poi[0], poi[1]);

			for(int j=0; j<(int)elements.size(); j++)
			{
				for(int k=0; k<(int)elements[j].nodes_id.size(); k++)
				{
					//线段在力边界线上
					vector<int> noid;
					if(k==(int)elements[j].nodes_id.size()-1)
					{
						noid.push_back(elements[j].nodes_id[k]);
						noid.push_back(elements[j].nodes_id[0]);
					}
					else
					{
						noid.push_back(elements[j].nodes_id[k]); 
						noid.push_back(elements[j].nodes_id[k+1]); 
					}
					//点到线段的两个端点的距离等于线段的长度，那么这个点在此线段内
					Point_2D poi_nod0(nodes[noid[0]].x, nodes[noid[0]].y);
					Point_2D poi_nod1(nodes[noid[1]].x, nodes[noid[1]].y);
					if(line.contain(poi_nod0)&&line.contain(poi_nod1))
					{
						Side side_temp;
						side_temp.nodes_id = noid;
						side_temp.type = 1;
						sides.push_back(side_temp);
						//注意，这个地方没有查是否有重复的边
					}
				}
			}
			//Go through sides
			if((int)sides.size()==0) { hout << "注意！载荷边界线位置输入有误（可能不经过任何网格节点），请重新输入载荷参数！" << endl; return 0; }
			for(int j=0; j<(int)sides.size(); j++)
			{
				double s_length = sqrt((nodes[sides[j].nodes_id[1]].x-nodes[sides[j].nodes_id[0]].x)*(nodes[sides[j].nodes_id[1]].x-nodes[sides[j].nodes_id[0]].x)
													 +(nodes[sides[j].nodes_id[1]].y-nodes[sides[j].nodes_id[0]].y)*(nodes[sides[j].nodes_id[1]].y-nodes[sides[j].nodes_id[0]].y));				
				
				for(int m=0; m<2; m++)
				{
					for(int n=0; n<(int)load.load_type[i].size(); n++)
					{
						if(load.load_type[i][n]=="Force_x") equright[3*sides[j].nodes_id[m]] += 0.5*value[n]*s_length;
						else if(load.load_type[i][n]=="Force_y") equright[3*sides[j].nodes_id[m]+1] += 0.5*value[n]*s_length;
					}
				}
			}
		}
		else if(load.domain_type[i]=="Zone")
		{
			double	xmin	=	load.coef[i][0];
			double	xmax	=	load.coef[i][1];
			double	ymin	=	load.coef[i][2];
			double	ymax	=	load.coef[i][3];
			for(int j=0; j<(int)elements.size(); j++)
			{
				vector<int> nodn;
				for(int k=0; k<(int)elements[j].nodes_id.size(); k++)
				{
					int enid =  elements[j].nodes_id[k];
					if((nodes[enid].x-xmin>=Zero)&&(xmax-nodes[enid].x>=Zero)&&
						(nodes[enid].y-ymin>=Zero)&&(ymax-nodes[enid].y>=Zero))
					{
						nodn.push_back(enid);
					}
				}
				if((int)nodn.size()==4)  //An element in this zone
				{
					//calculate area (Wiki chinese)
					double sa=(nodes[nodn[0]].x-nodes[nodn[1]].x)*(nodes[nodn[0]].x-nodes[nodn[1]].x)
									+(nodes[nodn[0]].y-nodes[nodn[1]].y)*(nodes[nodn[0]].y-nodes[nodn[1]].y)
									+(nodes[nodn[0]].z-nodes[nodn[1]].z)*(nodes[nodn[0]].z-nodes[nodn[1]].z);

					double sb=(nodes[nodn[1]].x-nodes[nodn[2]].x)*(nodes[nodn[1]].x-nodes[nodn[2]].x)
									+(nodes[nodn[1]].y-nodes[nodn[2]].y)*(nodes[nodn[1]].y-nodes[nodn[2]].y)
									+(nodes[nodn[1]].z-nodes[nodn[2]].z)*(nodes[nodn[1]].z-nodes[nodn[2]].z);

					double sc=(nodes[nodn[2]].x-nodes[nodn[3]].x)*(nodes[nodn[2]].x-nodes[nodn[3]].x)
									+(nodes[nodn[2]].y-nodes[nodn[3]].y)*(nodes[nodn[2]].y-nodes[nodn[3]].y)
									+(nodes[nodn[2]].z-nodes[nodn[3]].z)*(nodes[nodn[2]].z-nodes[nodn[3]].z);
	
					double sd=(nodes[nodn[3]].x-nodes[nodn[0]].x)*(nodes[nodn[3]].x-nodes[nodn[0]].x)
									+(nodes[nodn[3]].y-nodes[nodn[0]].y)*(nodes[nodn[3]].y-nodes[nodn[0]].y)
									+(nodes[nodn[3]].z-nodes[nodn[0]].z)*(nodes[nodn[3]].z-nodes[nodn[0]].z);

					double se=(nodes[nodn[0]].x-nodes[nodn[2]].x)*(nodes[nodn[0]].x-nodes[nodn[2]].x)
									+(nodes[nodn[0]].y-nodes[nodn[2]].y)*(nodes[nodn[0]].y-nodes[nodn[2]].y)
									+(nodes[nodn[0]].z-nodes[nodn[2]].z)*(nodes[nodn[0]].z-nodes[nodn[2]].z);

					double sf=(nodes[nodn[1]].x-nodes[nodn[3]].x)*(nodes[nodn[1]].x-nodes[nodn[3]].x)
									+(nodes[nodn[1]].y-nodes[nodn[3]].y)*(nodes[nodn[1]].y-nodes[nodn[3]].y)
									+(nodes[nodn[1]].z-nodes[nodn[3]].z)*(nodes[nodn[1]].z-nodes[nodn[3]].z);
					
					double s_area = 0.25*sqrt(4*se*sf-(sb+sd-sa-sc)*(sb+sd-sa-sc));

					for(int m=0; m<(int)nodn.size(); m++)
					{
						for(int n=0; n<(int)load.load_type[i].size(); n++)
						{
							if(load.load_type[i][n]=="Force_x") equright[3*nodn[m]] += 0.25*value[n]*s_area;
							else if(load.load_type[i][n]=="Force_y") equright[3*nodn[m]+1] += 0.25*value[n]*s_area;
						}
					}
				}
			}
		}
		else
		{
			cout << "Error: Failed to find the type of loaded domain!" << endl; 
			hout << "Error: Failed to find the type of loaded domain!" << endl;
			return 0; 
		}
	}

	return 1;
}
//===============================================================
