#include"func_state.h"
#include"global_variables.h"
#include <iostream>
#include<math.h>
#include<string.h>
#include <time.h>
using namespace std;


void memetic()
{
	Pe_factor = 1;
	Lambda = 5;
	Incre = 2;
	Best_sol_one_run_obj = -MAXNUM;
	initial_population();
	int generations = 0;
	//while (1.0*(clock() - Start_time) / CLOCKS_PER_SEC < Time_limit)
	while(generations < Max_generations)
	{
		knapsack_cross_over(Cur_sol);
		build_delta_matrix(Cur_sol);
		build_global_data(Cur_sol);
		
		double f_part1 = compute_obj_part1(Cur_sol);
		double f_part2 = compute_obj_part2(Cur_sol);
		double f_sol = f_part1 + f_part2;
		update_population(Cur_sol, f_sol);
		generations++;
		cout << "generations=" << generations << ", Best_sol_one_run_obj=" << Best_sol_one_run_obj << ", time=" << 1.0*(clock() - Start_time) / CLOCKS_PER_SEC << endl;
	}
}



int main(int argc, char *argv[])
{
	if (argc < 5)
	{
		cout << "HESA usage: input_file output_stat_file output_sol_file seed" << endl;
		cout << "input_file is the instance name, output_stat_file is a file used to store the running information, output_sol_file stores the solution information, and seed is the random seed, such as 1, 2, ...";
		exit(-1);
	}
	//srand(unsigned(time(NULL)));
	Instance_name = argv[1];
	char *out_stat_name = argv[2];
	char *out_sol_name = argv[3];
	int seed = atoi(argv[4]);
	srand(seed);
	
	//the following are used parameters
	Pop_size = 10;
	Num_greedy_ini = 1000;
	Pert_str = 0.1;
	Tabu_tenure = 50;
	Tabu_depth = 2000;

	Time_limit = 3600.0; 		//it is unused
	Max_generations = 100;


	read_instance();
	allocate_memory();

	int runs = 1;
	double ff_best = -MAXNUM;
	double ff_worst = MAXNUM;
	int hit = 0;
	double ff_avg = 0;
	double avg_time = 0;
	build_neighbors();


	for (int i = 0; i < runs; i++)
	{
		Start_time = clock();
		memetic();	
		//memetic_mutationRL();
		//memetic_crossover_RL();
		//local_search_tps();
		//adaptive_fits();
		if (Best_sol_one_run_obj > ff_best + PRECISION)
		{
			hit = 1;
			ff_best = Best_sol_one_run_obj;
			memcpy(Global_best_sol, Best_sol_one_run, sizeof(int)*Instance.items);
		}
		else if (fabs(Best_sol_one_run_obj - ff_best) <= PRECISION)
			hit++;
		if (Best_sol_one_run_obj < ff_worst - PRECISION)
			ff_worst = Best_sol_one_run_obj;
		ff_avg += Best_sol_one_run_obj;
		avg_time += Run_time;
		//out_results_one_run(out_sol_name, Instance_name, Best_k, Run_time);
	}
	ff_avg = ff_avg / runs;
	avg_time /= runs;
	proof(Global_best_sol, ff_best);
	out_results_best_sol(out_sol_name, Instance_name, ff_best, Global_best_sol);
	out_results_stat(ff_best, ff_avg, ff_worst, avg_time, hit, out_stat_name, Instance_name);
	free_memory();
}
