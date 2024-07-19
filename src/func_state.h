#ifndef _STATEMENT_H
#define _STATEMENT_H


void read_instance();		//read instance
void initial_global_data();
void build_global_data(int *sol);
void update_global_data(int item, int old_knap, int new_knap);
void update_global_data_exceed(int old_knap, int new_knap, int r1, int r2);			//for infeasible tabu search

void random_initial_sol(int *sol);	//unused
void greedy_helpfunc_value_density(int *sol, int *can_item, int *can_knaps, int *address_knaps, int &len_knap);
void greedy_initial_sol(int *sol);
void initial_sol(int *sol);
void initial_sol_helpfunc(int *sol);

bool satisfy_all_constraints_relocate(int item_j, int old_knap, int new_knap);
bool satisfy_all_constraints_swap(int item_j1, int knap_j1, int item_j2, int knap_j2);
bool satisfy_all_constraints_2_1_exchange(int item_x, int item_y, int knap_x_y, int item_z, int knap_z);
bool satisfy_nr_constraints_relocate(int item_j, int old_knap, int new_knap);
bool satisfy_nr_constraints_swap(int item_j1, int knap_j1, int item_j2, int knap_j2);
bool satisfy_cap_constraints_swap(int item_j1, int knap_j1, int item_j2, int knap_j2);
bool satisfy_cap_constraints_relocate(int item_j, int old_knap, int new_knap);

void build_delta_matrix(int *sol);
void update_delta_matrix(int item, int old_knap, int new_knap);
double compute_obj_part1(int *sol);
double compute_obj_part2(int *sol);
double compute_delta_relocate(int *sol, int item_j, int knap);
double compute_delta_swap(int *sol, int item_j, int item_j2);
double compute_delta_2_1_exchange(int *sol, int item_x, int item_y, int item_z);
double compute_delta_relocate_g1_g2(int *sol, int item_j, int new_knap);
double compute_delta_swap_g1_g2(int *sol, int item_j, int item_j2);
int compute_delta_relocate_g2(int *sol, int item_j, int new_knap);
int compute_delta_swap_g2(int *sol, int item_j, int item_j2);
double compute_g1();
double compute_delta_relocate_g1(int *sol, int item_j, int new_knap);
double compute_delta_swap_g1(int *sol, int item_j, int item_j2);

void relocate_move(int *sol, int item, int knap);
void swap_move(int *sol, int item1, int item2);
void two_one_exchange_move(int *sol, int item_x, int item_y, int item_z);

void build_neighbors();
void random_perturbation(int *sol, double &f_sol, int pert_len);
void mutation_RL(int *sol, double &f_sol, int pert_len);
void greedy_helpfunc_proba_matrix(int *sol, int *can_item, int *can_knaps, int *address_knaps, int &len_knap);

//for local search
double descent_local_search(int *sol, double &f_sol);
double feasible_tabu_search(int *sol, double &f_sol);
double infeasible_tabu_search(int *sol, double &f_sol);

//for reinforcement learning
void probability_updateing(int *pre_sol, int *cur_sol);
void probability_smoothing();

//for memetic search
void initial_population();										//initial the population
void knapsack_cross_over(int *offspring);						//knapsack-based crossover operator
void backbone_cross_over(int *offspring);						//backbone based crossover operator, according to yuning chen ITEC paper
void uniform_cross_over(int *offspring);						//uniform crossover operator
void update_population(int *offspring, double f_off);			//update the population
bool is_same_solution(int *sol_1, int *sol_2);


void update_best_sol_one_run(double f_sol, int *sol);
void allocate_memory();
void free_memory();
void out_results_best_sol(char *out_filename, char *instance_name, double ff, int *sol);
void out_results_stat(double best, double ave, double worst, double avg_time, int hit, char *stat_filename, char instance[]);
void proof(int *sol, double ff);
void check_move(int *sol, double fs, double gs);
void verify_sol_from_file(int *sol, double ff, int type);


#endif