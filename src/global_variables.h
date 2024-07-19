#ifndef _GLOBAL_VARIABLES_H
#define _GLOBAL_VARIABLES_H
#include<bitset>

#define MAXCAN 200
#define MAXNUM 9999999999
#define PRECISION 1.0e-8
#define MUTATION_INI_POP
#define DEBUG_INITIAL_SOL111
#define DEBUG_READ_INSTANCE111

typedef struct instance_data
{
	int items;				//the number of items
	int knaps;				//the number of knapsacks
	int classes;			//the number of class
	int u;					//a positive large number
	double *we;				//weight of each item
	double capac;			//capacity for each knapsack
	double *po;
	double **q_ij;			//profit between item i and item j
	int **t_rj;
	int *t_jr;
	double *sr;
	int *nr;
	double **psi;
	double **p_jk;

}instance_data;

typedef struct neighborhood{
	int  type;
	int  v;
	int  g;
	int  x;
	int  y;
}neighborhood;

extern char *Instance_name;
extern instance_data Instance;
extern int**Size_kr;			//Size_kr[k][r] indiactes the number of items in knapsack k with the class r
extern int*Size_r;
extern double*Capacity1_k;
extern double*Capacity2_k;
extern double*Exceed_cap_k;
extern int *Exceed_class_r;
extern double**Delta_matrix_jk;

extern int Num_greedy_ini;

//for infeasible tabu search
extern int Pe_factor;			//penalty factor
extern int Lambda;				//every Lambda consecutive iterations 
extern int Incre;

//for tabu search
extern int Tabu_depth;
extern int Tabu_tenure;
extern int**Tabu_list;

//for perturbation
extern neighborhood*Neighbors;
extern double Pert_str;				//perturbation strength


//for reinforcement learning
extern double RL_omega;
extern double RL_alpha;
extern double RL_beta;
extern double RL_gamma;
extern double RL_smoothing_coeff;
extern double RL_proba_threshold;
extern double **Proba_matrix_jk;	//probability learning matrix

//for memetic search
extern int Pop_size;				//population size
extern int **Pop_sol;
extern double*Pop_obj;

extern int *Cur_sol;
extern int *Local_best_sol;
extern int *Best_sol_one_run;
extern int *Global_best_sol;
extern double Best_sol_one_run_obj;

extern double Time_limit, Start_time, Run_time;
extern int Max_generations;
#endif
