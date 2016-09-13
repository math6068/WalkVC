#include "WalkVC.h"

int try_step=100;


int edge_cand;

void cover_LS()
{
	int		remove_v, add_v;
	int		remove_dscore, add_dscore;
	int		e,v1,v2;
	int		i;

	step = 1;

	while(1)
	{
		if (uncov_stack_fill_pointer == 0)//update best solution if needed
		{
			update_best_sol();
			
			//if (c_size==optimal_size) return;
				
			update_target_size();
			
			continue;
		}
		
		if(step%try_step==0) //check cutoff
		{
			times(&finish);
			double elap_time = (finish.tms_utime + finish.tms_stime - start_time)/sysconf(_SC_CLK_TCK);
			if(elap_time >= cutoff_time)return;
		}
		
		remove_v = choose_remove_v();
		
		remove(remove_v);

		
		e = uncov_stack[rand()%uncov_stack_fill_pointer];
		v1 = edge[e].v1;
		v2 = edge[e].v2;

//FastVC mode 
		if(dscore[v1]>dscore[v2] || (dscore[v1]==dscore[v2] && time_stamp[v1]<time_stamp[v2]) )
			add_v=v1;
		else add_v=v2;
     
        //if(rand()%2==1)
			 //add_v=v1;
		//else add_v=v2;

// Random Walk with randuncoveredge mode; can be useful
       /* if (rand() % 1000 < 980 )
		{
        if(dscore[v1]>dscore[v2] || (dscore[v1]==dscore[v2] && time_stamp[v1]<time_stamp[v2]) )
			 add_v=v1;
		else add_v=v2;
		}
        else 
        {
        if(rand()%2==1)
			 add_v=v1;
		else add_v=v2;
        } */

		add(add_v);
		
		int index = index_in_remove_cand[remove_v];
		index_in_remove_cand[remove_v] = 0;
		
		remove_cand[index] = add_v;
		index_in_remove_cand[add_v] = index;
		
		time_stamp[add_v]=time_stamp[remove_v]=step;

		step++;
	}

}


int main(int argc, char* argv[])
{
	int seed,i;

	//cout<<"c This is a local search solver for the Minimum Vertex Cover problem."<<endl;
	
	if(build_instance(argv[1])!=1){
		cout<<"can't open instance file"<<endl;
		return -1;
	}
	
	//optimal_size=0;
	i=2;
// run WalkVC input seed and cutoff in this order
	sscanf(argv[i++],"%d",&seed);
	sscanf(argv[i++],"%d",&cutoff_time);

//generate random number sequence...
	srand(seed);
	//cout<<seed<<' ';
	//cout<<"c This is WalkVC, solving instnce "<<argv[1]<<endl;
		
	times(&start);
	start_time = start.tms_utime + start.tms_stime;

   	init_sol();   //constructVC in Al.1.
#ifdef individual_analysis_on_init_sls_mode
		times(&finish);
		init_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
		init_time = round(init_time * 100)/100.0;
#endif    
	//if(c_size + uncov_stack_fill_pointer > optimal_size ) 
	//{
		//cout<<"c Start local search..."<<endl;
		cover_LS();
	//}
#ifdef individual_analysis_on_init_sls_mode
		sls_time = best_comp_time - init_time;
		sls_step_speed_per_ms = double(best_step) * 0.001 / sls_time;
#endif			
	//check solution
	if(check_solution()==1)
	{
		//print_solution();
		cout<<"o "<<best_c_size<<endl;
		cout<<"c searchSteps "<<best_step<<endl;
		cout<<"c solveTime "<<best_comp_time<<endl;
#ifdef 	individual_analysis_on_init_sls_mode
		//cout<<"c initTime " << init_time << endl;
		//cout<<"c slsTime " << sls_time << endl;
		cout<<"c stepSpeed(/ms) "<< sls_step_speed_per_ms << endl;
#endif			
		//cout<<best_c_size<<' '<<best_comp_time<<endl;
	}
	
	free_memory();

	return 0;
}
