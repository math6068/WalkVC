/************************************************
** This is a local search solver for Minimum Vertex Cover.                                                       
************************************************/

/************************************************
** Date:	2015.2.2  
** FastVC
** Author: Shaowei Cai, caisw@ios.ac.cn   
**		   Key Laboratory of Computer Science,
**		   Institute of Software, Chinese Academy of Sciences, 
**		   Beijing, China     
** Revison: Revise the choose_remove_v() funtion, i.e., with 
** probability p follow BMS, with probability 1-p remove a random 
** vertex.
** By Zongjie Ma, Yi Fan
** Institute for Integrated and Intelligent Systems, Griffith ** University, Brisbane, Australia.
** Date:  2016.9.10                                                                
************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <sys/times.h>
#include <cmath>

using namespace std;


#define pop(stack) stack[--stack ## _fill_pointer]
#define push(item, stack) stack[stack ## _fill_pointer++] = item

/*max vertex count and max edge count*/
#define	MAXV	15000000
#define MAXE	80000000

#define individual_analysis_on_init_sls_mode //for step speed statistics

#ifdef individual_analysis_on_init_sls_mode
double	init_time;
double	sls_time;
double	sls_step_speed_per_ms;
#endif

tms start, finish;
int start_time;

struct Edge{
	int v1;
	int v2;
};

/*parameters of algorithm*/
long long	max_steps;			//step limit
int			cutoff_time;			//time limit
long long	step;
int			optimal_size;			//terminate the algorithm before step limit if it finds a vertex cover of optimal_size

/*parameters of the instance*/
int		v_num;//|V|: 1...v
int		e_num;//|E|: 0...e-1

/*structures about edge*/
Edge	edge[MAXE];  

/*structures about vertex*/
int		dscore[MAXV];			//dscore of v
long long	time_stamp[MAXV];


//from vertex to it's edges and neighbors
int*	v_edges[MAXV];	//edges related to v, v_edges[i][k] means vertex v_i's k_th edge
int*	v_adj[MAXV];		//v_adj[v_i][k] = v_j(actually, that is v_i's k_th neighbor)
int		v_degree[MAXV];	//amount of edges (neighbors) related to v


/* structures about solution */
//current candidate solution
int		c_size;						//cardinality of C
bool	v_in_c[MAXV];				//a flag indicates whether a vertex is in C
int		remove_cand[MAXV];			//remove candidates, an array consists of only vertices in C, not including tabu_remove
int		index_in_remove_cand[MAXV];
int		remove_cand_size;

//best solution found
int		best_c_size;
bool	best_v_in_c[MAXV];			//a flag indicates whether a vertex is in best solution
double  best_comp_time;
long    best_step;


//uncovered edge stack
int		uncov_stack[MAXE];		//store the uncov edge number
int		uncov_stack_fill_pointer;
int		index_in_uncov_stack[MAXE];//which position is an edge in the uncov_stack


//CC and taboo
//int 	conf_change[MAXV];
//int		tabu_remove=0;


/* functions declaration */
int build_instance(char *filename);
void init_sol();
void cover_LS();
void add(int v);
void remove(int v);
void update_edge_weight();
void cover_rest_edges();
int check_solution();



void update_best_sol()
{
	int i;

	for (i=1;i<=v_num;i++)
		best_v_in_c[i] = v_in_c[i];
	
	best_c_size = c_size;
	times(&finish);
	best_comp_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
	best_comp_time = round(best_comp_time * 100)/100.0;
	best_step = step;
}



int v_degree_tmp[MAXV];

int build_instance(char *filename)
{
	char line[1024];
	char tempstr1[10];
	char tempstr2[10];
	int  v,e;
	
	char	tmp;
	int		v1,v2;
	
	ifstream infile(filename);
    if(infile==NULL) return 0;

	/*** build problem data structures of the instance ***/
	infile.getline(line,1024);
	while (line[0] != 'p') infile.getline(line,1024);
	sscanf(line, "%s %s %d %d", tempstr1, tempstr2, &v_num, &e_num);
	
	if (v_num > MAXV) {
		cout<<"the number of vertices ("<<v_num<<") exceeds MAXV ("<<MAXV<<")."<<endl;
		exit(0);
	}
	if (e_num > MAXE) {
		cout<<"the number of vertices ("<<e_num<<") exceeds MAXV ("<<MAXE<<")."<<endl;
		exit(0);
	}

	/* read edges and compute v_degree */
	for (v=1; v<=v_num; v++) v_degree[v] = 0;
	
	for (e=0; e<e_num; e++)
	{
		infile>>tmp>>v1>>v2;
		v_degree[v1]++;
		v_degree[v2]++;
		
		edge[e].v1 = v1;
		edge[e].v2 = v2;
	}
	infile.close();
	
	/* build v_adj and v_edges arrays */
	for (v=1; v<=v_num; v++)
	{
		v_adj[v] = new int[v_degree[v]];
		v_edges[v] = new int[v_degree[v]];
	}
	
	//for(v=1; v<=v_num; v++) v_degree_tmp[v]=0;
	
	for (e=0; e<e_num; e++)
	{
		v1=edge[e].v1;
		v2=edge[e].v2;

		v_edges[v1][v_degree_tmp[v1]] = e;
		v_edges[v2][v_degree_tmp[v2]] = e;

		v_adj[v1][v_degree_tmp[v1]] = v2;
		v_adj[v2][v_degree_tmp[v2]] = v1;

		v_degree_tmp[v1]++;
		v_degree_tmp[v2]++;
	}
	
	return 1;
}


void free_memory()
{
	for (int v=1; v<=v_num; v++)
	{
		delete[] v_adj[v];
		delete[] v_edges[v];
	}
}

void reset_remove_cand()
{
	int v,j;
	j=0;
	for (v=1; v<=v_num; v++)
	{
		if(v_in_c[v]==1)
		{
			remove_cand[j] = v;
			index_in_remove_cand[v]=j;
			j++;
		}
		else index_in_remove_cand[v]=0;
	}
	
	remove_cand_size = j;
}



void update_target_size()
{
	c_size--;

	/*int cand_counts=50;
	int i,v;

	int best_remove_v = remove_cand[rand()%remove_cand_size];
	
	for (i=1; i<cand_counts; ++i)
	{
		v = remove_cand[rand()%remove_cand_size];
	
		if( dscore[v] < dscore[best_remove_v])
			continue;
		else if( dscore[v]> dscore[best_remove_v] )
			best_remove_v = v;
		else if (time_stamp[v]<time_stamp[best_remove_v])
			best_remove_v = v;
	}  */
	

	int v,i;
	int best_dscore;
	int best_remove_v;//vertex with the highest improvement in C

	best_remove_v = remove_cand[0];
	best_dscore = dscore[best_remove_v];
	
	if(dscore[best_remove_v]!=0)
	{
		for (i=1; i<remove_cand_size; ++i)
		{
			v = remove_cand[i];
			
			if(dscore[v]==0) break;

			if (dscore[v] > dscore[best_remove_v])
				best_remove_v = v;
		}
	} 
	
	remove(best_remove_v);
	
	//remove best_remove_v from remove_cand, and move the last vertex in remove_cand to the position
	int last_remove_cand_v = remove_cand[--remove_cand_size];
	int index = index_in_remove_cand[best_remove_v];
	remove_cand[index] = last_remove_cand_v;
	index_in_remove_cand[last_remove_cand_v] = index;

	//reset_remove_cand();
}




//update the best vertex in C 
int cand_count=50;
int choose_remove_v()
{   
   
// Random Walk with BMS mode  
   int i,v;
    int best_v;
	if (rand() % 1000 < 600 )
    {    best_v = remove_cand[rand()%remove_cand_size];
         for (i=1; i<cand_count; ++i)
	     {
		v = remove_cand[rand()%remove_cand_size];
	
		if( dscore[v] < dscore[best_v])
			continue;
		else if( dscore[v]> dscore[best_v] )
			best_v = v;
		else if (time_stamp[v]<time_stamp[best_v])
			best_v = v;
	      }
    }
    else 
         best_v = remove_cand[rand()%remove_cand_size];

  /*  int i,v;

	int best_v = remove_cand[rand()%remove_cand_size]; //the first v selected randomly.
	
	for (i=1; i<cand_count; ++i)
	{
		v = remove_cand[rand()%remove_cand_size];
	
		if( dscore[v] < dscore[best_v])
			continue;
		else if( dscore[v]> dscore[best_v] )
			best_v = v;
		else if (time_stamp[v]<time_stamp[best_v])
			best_v = v;
	}  */
	
	return best_v;   //notice that here doesnot need to update best_v's position information, because that when the vertex selected to add will take up its position, and renew the potition information ater that.
} 





inline
void uncover(int e) 
{
	index_in_uncov_stack[e] = uncov_stack_fill_pointer;
	push(e,uncov_stack);
}


inline
void cover(int e)
{
	int index,last_uncov_edge;

	//since the edge is satisfied, its position can be reused to store the last_uncov_edge
	last_uncov_edge = pop(uncov_stack);
	index = index_in_uncov_stack[e];
	uncov_stack[index] = last_uncov_edge;
	index_in_uncov_stack[last_uncov_edge] = index;
}


void init_sol()    //Al.2.
{
	int i,v,e;
	int v1, v2;

	/*** build solution data structures of the instance ***/
	for (v=1; v<=v_num; v++)
	{
		v_in_c[v] = 0;
		dscore[v] = 0;
		//conf_change[v] = 1;     //cc
		time_stamp[v]= 0; 
	}
	
	//construct a vertex cover
	c_size = 0;
	for (e=0; e<e_num; e++)
	{	
		v1=edge[e].v1;
		v2=edge[e].v2;
		
		if (v_in_c[v1]==0 && v_in_c[v2]==0)//if uncovered, choose the endpoint with higher degree
		{
			if(v_degree[v1] > v_degree[v2]) 
			{
				v_in_c[v1]=1;
			}
			else{
				v_in_c[v2]=1;
			}
			c_size++;
		}
	}

	//init uncovered edge stack
	uncov_stack_fill_pointer = 0;
	
	//calculate dscores
	for (e=0; e<e_num; e++)
	{
		v1=edge[e].v1;
		v2=edge[e].v2;
		
		if (v_in_c[v1]==1 && v_in_c[v2]==0) dscore[v1]--;
		else if (v_in_c[v2]==1 && v_in_c[v1]==0) dscore[v2]--;
	}
	

	//remove redundent vertices
	for (v=1; v<=v_num; v++)
	{
		if (v_in_c[v]==1 && dscore[v]==0) 
		{
			remove(v);
			c_size--;
		}
	}

	update_best_sol();//initialize the best found solution

	//cout<<c_size<<' '<<best_comp_time<<' ';
	
	reset_remove_cand();

}


void add(int v)
{
	v_in_c[v] = 1;
	dscore[v] = -dscore[v];
	
	int i,e,n;

	int edge_count = v_degree[v];
	
	for (i=0; i<edge_count; ++i)
	{
		e = v_edges[v][i];// v's i'th edge
		n = v_adj[v][i];//v's i'th neighbor

		if (v_in_c[n]==0)//this adj isn't in cover set
		{
			dscore[n]--;
			//conf_change[n] = 1;

			cover(e);
		}
		else
		{
			dscore[n]++; 
		}
	}
	
}

void remove(int v)
{
	v_in_c[v] = 0;
	dscore[v] = -dscore[v];
	//conf_change[v] = 0;

	int i,e,n;

	int edge_count = v_degree[v];
	for (i=0; i<edge_count; ++i)
	{
		e = v_edges[v][i];
		n = v_adj[v][i];

		if (v_in_c[n]==0)//this adj isn't in cover set
		{
			dscore[n]++;  
			//conf_change[n] = 1;

			uncover(e);
		}
		else
		{
			dscore[n]--; 
		}
	}
}


void print_solution()
{
	for (int i=1; i<=v_num; i++)
	{
		if (best_v_in_c[i]==1)//output vertex cover
			cout<<i<<'\t';
	}
	cout<<endl;
}


int check_solution()
{
	for(int e=0; e<e_num; ++e)
	{
		if(best_v_in_c[edge[e].v1]!=1 && best_v_in_c[edge[e].v2]!=1)
		{
			cout<<"c error: uncovered edge "<<e<<endl;
			return 0;
		}
	}
	
	int verified_vc_size=0;
	for (int i=1; i<=v_num; i++)
	{
		if (best_v_in_c[i]==1)
			verified_vc_size++;
	}
	
	if(best_c_size==verified_vc_size) return 1;
	
	else{
		cout<<"c error: claimed best found vc size!=verified vc size"<<endl;
		cout<<"c claimed best found vc size="<<best_c_size<<endl;
		cout<<"c verified vc size="<<verified_vc_size<<endl;
		return 0;
	}
}

