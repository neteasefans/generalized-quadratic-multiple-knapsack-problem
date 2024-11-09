#include"func_state.h"
#include"global_variables.h"
#include<time.h>
#include<string.h>
#include<math.h>
#include<iostream>
using namespace std;

//initial the population
void initial_population()
{			
	int *sol = new int[Instance.items];
	int **sol_tmp = new int*[Pop_size];
	for (int p = 0; p < Pop_size; p++)
		sol_tmp[p] = new int[Instance.items];
	double *obj_tmp = new double[Pop_size];
	for (int i = 0; i < Pop_size; i++)
		obj_tmp[i] = -MAXNUM;
	double time1 = clock();
	for (int i = 0; i < Num_greedy_ini; i++)
	{
		greedy_initial_sol(sol);
		double f_sol = 0.0;
		for (int j = 0; j < Instance.items; j++)
			if (sol[j] < Instance.knaps)
				f_sol += Instance.p_jk[j][sol[j]];
		for (int j = 0; j < Instance.items; j++)
			for (int j2 = j + 1; j2 < Instance.items; j2++)
				if (sol[j] == sol[j2] && sol[j] != Instance.knaps)
					f_sol += Instance.q_ij[j][j2];	
		int sel_pos = Pop_size;
		for (int i = 0; i < Pop_size; i++)
		{
			if (f_sol > obj_tmp[i] + PRECISION)
			{
				sel_pos = i;
				break;
			}
		}
		for (int i = Pop_size - 1; i > sel_pos; i--)
		{
			obj_tmp[i] = obj_tmp[i - 1];
			memcpy(sol_tmp[i], sol_tmp[i - 1], sizeof(int)*Instance.items);
		}
		if (sel_pos != Pop_size)
		{
			obj_tmp[sel_pos] = f_sol;
			memcpy(sol_tmp[sel_pos], sol, sizeof(int)*Instance.items);
		}
		//cout << "i=" << i << ", f_sol=" << f_sol << ", sel_pos=" << sel_pos << endl;
		//for (int i = 0; i < Pop_size; i++)
		//	cout << "i=" << i << ", obj_tmp=" << obj_tmp[i] << endl;
	}
	double initial_time = 1.0*(clock() - time1) / CLOCKS_PER_SEC;
	cout << "in initial_sol func, initial_time=" << initial_time << endl;

	for (int pop_len = 0; pop_len < Pop_size; pop_len++)
	//int pop_len = 0;
	//while (pop_len < Pop_size)
	{
		memcpy(Pop_sol[pop_len], sol_tmp[pop_len], sizeof(int)*Instance.items);
		//greedy_initial_sol(Pop_sol[pop_len]);
		build_delta_matrix(Pop_sol[pop_len]);
		build_global_data(Pop_sol[pop_len]);
		double f_part1 = compute_obj_part1(Pop_sol[pop_len]);
		double f_part2 = compute_obj_part2(Pop_sol[pop_len]);
		Pop_obj[pop_len] = f_part1 + f_part2;
		descent_local_search(Pop_sol[pop_len], Pop_obj[pop_len]);
		infeasible_tabu_search(Pop_sol[pop_len], Pop_obj[pop_len]);		
		
#ifdef MUTATION_INI_POP
		while (1)
		{
			bool is_exist = false;
			for (int i = 0; i < pop_len; i++)
			{
				is_exist = is_same_solution(Pop_sol[i], Pop_sol[pop_len]);
				if (is_exist)
					break;
			}			
			if (!is_exist)
				break;
			else
				random_perturbation(Pop_sol[pop_len], Pop_obj[pop_len], int(Pert_str*Instance.items));
			cout << "pop_len=" << pop_len << ", is_exist=" << is_exist << endl;
		}	
#endif	
	}
	for (int i = 0; i < Pop_size; i++)	
		proof(Pop_sol[i], Pop_obj[i]);	
	cout << "after initial_population func, " << ", Best_sol_one_run_obj=" << Best_sol_one_run_obj << endl;

	
	for (int i = 0; i < Pop_size; i++)
	{
		delete[]sol_tmp[i]; sol_tmp[i] = NULL;
	}
	delete[]sol_tmp; sol_tmp = NULL;
	delete[]obj_tmp; obj_tmp = NULL;
	delete[]sol; sol = NULL;
}

void calculate_kanp_obj(double *knap_obj, int *parent, int *offspring)
{
	for (int k = 0; k < Instance.knaps + 1; k++)
		knap_obj[k] = 0.0;
	for (int j = 0; j < Instance.items; j++)
		if (offspring[j] == Instance.knaps)			//item j is not allocated to any knap of the offspring 
			knap_obj[parent[j]] += Instance.p_jk[j][parent[j]];
	for (int j = 0; j < Instance.items; j++)
	{
		if (offspring[j] == Instance.knaps)
		{
			for (int j2 = j + 1; j2 < Instance.items; j2++)
				if (parent[j] == parent[j2] && offspring[j2] == Instance.knaps)
					knap_obj[parent[j]] += Instance.q_ij[j][j2];
		}
	}
}

//knapsack-based crossover operator
void knapsack_cross_over(int *offspring)
{
	double *knap_obj_par1 = new double[Instance.knaps + 1];
	double *knap_obj_par2 = new double[Instance.knaps + 1];
	int *can_item = new int[Instance.items];
	int *can_knaps = new int[Instance.knaps];
	int *address_knaps = new int[Instance.knaps];
	int*flag_knap = new int[Instance.knaps];
	int best_knap_par1[MAXCAN], best_knap_par2[MAXCAN];
	int best_knap_num_par1 = 0, best_knap_num_par2 = 0;
	int parent1 = rand() % Pop_size;
	int parent2 = rand() % Pop_size;
	while (parent1 == parent2)
		parent2 = rand() % Pop_size;
	for (int j = 0; j < Instance.items; j++)
		offspring[j] = Instance.knaps;			//item j is not allocated to any knap of the offspring 
	initial_global_data();
	int len_knap = Instance.knaps;
	for (int k = 0; k < Instance.knaps; k++)	//knapsack k can be selected
		flag_knap[k] = 1;
	//first step: knapsack based preserving items
	int step_order = 0;
	while (len_knap > 0)
	{	
		if (step_order % 2 == 0)
		{
			calculate_kanp_obj(knap_obj_par1, Pop_sol[parent1], offspring);
			double knap_obj_maximum_par1 = -MAXNUM;
			for (int k = 0; k < Instance.knaps; k++)
			{
				if (knap_obj_par1[k] > knap_obj_maximum_par1 + PRECISION && flag_knap[k] == 1)
				{
					best_knap_num_par1 = 0;
					knap_obj_maximum_par1 = knap_obj_par1[k];
					best_knap_par1[best_knap_num_par1++] = k;
				}
				else if (fabs(knap_obj_par1[k] - knap_obj_maximum_par1) <= PRECISION && best_knap_num_par1 < MAXCAN)
					best_knap_par1[best_knap_num_par1++] = k;
			}
			int rx = rand() % best_knap_num_par1;
			int sel_knap = best_knap_par1[rx];
			for (int j = 0; j < Instance.items; j++)
			{
				if (Pop_sol[parent1][j] == sel_knap && offspring[j] == Instance.knaps && satisfy_all_constraints_relocate(j, offspring[j], sel_knap))
				{
					update_global_data(j, offspring[j], sel_knap);
					offspring[j] = sel_knap;
					flag_knap[sel_knap] = 0;
				}
			}

		}
		else
		{
			calculate_kanp_obj(knap_obj_par2, Pop_sol[parent2], offspring);
			double knap_obj_maximum_par2 = -MAXNUM;
			for (int k = 0; k < Instance.knaps; k++)
			{
				if (knap_obj_par2[k] > knap_obj_maximum_par2 + PRECISION && flag_knap[k] == 1)
				{
					best_knap_num_par2 = 0;
					knap_obj_maximum_par2 = knap_obj_par2[k];
					best_knap_par2[best_knap_num_par2++] = k;
				}
				else if (fabs(knap_obj_par2[k] - knap_obj_maximum_par2) <= PRECISION && best_knap_num_par2 < MAXCAN)
					best_knap_par2[best_knap_num_par2++] = k;
			}
			int rx = rand() % best_knap_num_par2;
			int sel_knap = best_knap_par2[rx];
			for (int j = 0; j < Instance.items; j++)
			{
				if (Pop_sol[parent2][j] == sel_knap && offspring[j] == Instance.knaps && satisfy_all_constraints_relocate(j, offspring[j], sel_knap))
				{
					update_global_data(j, offspring[j], sel_knap);
					offspring[j] = sel_knap;
					flag_knap[sel_knap] = 0;
				}
			}
		}		
		len_knap--;
		step_order++;
		//proof(offspring, 0.0);
	}	
	
	//second step: complete the partial solution in a greedy manner
	len_knap = 0;
	for (int k = 0; k < Instance.knaps; k++)
	{
		address_knaps[k] = len_knap;
		can_knaps[len_knap] = k;
		len_knap++;
	}
	while (len_knap > 0)
	{
		greedy_helpfunc_value_density(offspring, can_item, can_knaps, address_knaps, len_knap);
	}

	delete[]can_item; can_item = NULL;
	delete[]can_knaps; can_knaps = NULL;
	delete[]address_knaps; address_knaps = NULL;

	delete[]knap_obj_par1; knap_obj_par1 = NULL;
	delete[]knap_obj_par2; knap_obj_par2 = NULL;
	delete[]flag_knap; flag_knap = NULL;
}

void uniform_cross_over(int *offspring)
{
	int *can_knaps = new int[Instance.knaps];
	int *address_knaps = new int[Instance.knaps];
	int *can_item = new int[Instance.items];

	int parent1 = rand() % Pop_size;
	int parent2 = rand() % Pop_size;
	while (parent1 == parent2)
		parent2 = rand() % Pop_size;
	for (int j = 0; j < Instance.items; j++)
		offspring[j] = Instance.knaps;			//item j is not allocated to any knap of the offspring 
	initial_global_data();
	for (int j = 0; j < Instance.items; j++)
	{
		int rx = rand() % 2;
		int new_knap = -1;
		if (rx == 0)
			new_knap = Pop_sol[parent1][j];
		else
			new_knap = Pop_sol[parent2][j];
		if (new_knap != Instance.knaps  && satisfy_all_constraints_relocate(j, offspring[j], new_knap))
		{
			update_global_data(j, offspring[j], new_knap);
			offspring[j] = new_knap;
		}
	}
	int len_knap = 0;
	for (int k = 0; k < Instance.knaps; k++)
	{
		address_knaps[k] = len_knap;
		can_knaps[len_knap] = k;
		len_knap++;
	}
	while (len_knap > 0)
	{
		greedy_helpfunc_value_density(offspring, can_item, can_knaps, address_knaps, len_knap);
	}
	delete[]can_knaps; can_knaps = NULL;
	delete[]address_knaps; address_knaps = NULL;
	delete[]can_item; can_item = NULL;
}


int common_vertices_between_two_knaps(int *sol1, int k1, int *sol2, int k2, int *off)
{
	int common_nodes = 0;
	for (int j1 = 0; j1 < Instance.items; j1++)
	{
		if (sol1[j1] == k1 && off[j1] == Instance.knaps)
		{
			for (int j2 = 0; j2 < Instance.items; j2++)
			{
				if (sol2[j2] == k2 && off[j2] == Instance.knaps && j1 == j2)
				{
					common_nodes++;
					break;
				}
			}
		}
	}
	return common_nodes;
}

double value_density_item(int item, int knap, int *off)
{
	double vd_jk = Instance.p_jk[item][knap];
	for (int j = 0; j < Instance.items; j++)
		if (off[j] == knap)
			vd_jk += Instance.q_ij[item][j];
	if (fabs(Instance.we[item] - 0) <= PRECISION)
		vd_jk = MAXNUM;
	else
	{
		int class_j_r = Instance.t_jr[item];
		if (Size_kr[knap][class_j_r] == 0)
			vd_jk /= (Instance.we[item] + Instance.sr[class_j_r]);
		else
			vd_jk /= Instance.we[item];
	}
	return vd_jk;
}

//backbone-based crossover operator
void backbone_cross_over(int *offspring)
{
	int *can_item = new int[Instance.items];
	int *can_knaps = new int[Instance.knaps];
	int *address_knaps = new int[Instance.knaps];
	int *flag_knap1 = new int[Instance.knaps];
	int *flag_knap2 = new int[Instance.knaps];
	int best_knap_par1[MAXCAN], best_knap_par2[MAXCAN];
	int best_knap_num = 0;
	int parent1 = rand() % Pop_size;
	int parent2 = rand() % Pop_size;
	while (parent1 == parent2)
		parent2 = rand() % Pop_size;
	for (int j = 0; j < Instance.items; j++)
		offspring[j] = Instance.knaps;			//item j is not allocated to any knap of the offspring 
	initial_global_data();
	int len_knap = 0;
	for (int k = 0; k < Instance.knaps; k++)	//knapsack k can be selected
	{
		flag_knap1[k] = 1;
		flag_knap2[k] = 1;
	}
	//first step: backbone based preserving items

	while (len_knap < Instance.knaps)
	{
		int cn_maximum = -1;
		for (int k1 = 0; k1 < Instance.knaps; k1++)
		{
			if (flag_knap1[k1])
			{
				for (int k2 = 0; k2 < Instance.knaps; k2++)
				{
					if (flag_knap2[k2])
					{
						int common_nodes = common_vertices_between_two_knaps(Pop_sol[parent1], k1, Pop_sol[parent2], k2, offspring);
						if (common_nodes > cn_maximum)
						{
							cn_maximum = common_nodes;
							best_knap_num = 0;
							best_knap_par1[best_knap_num] = k1;
							best_knap_par2[best_knap_num] = k2;
							best_knap_num++;
						}
						else if (common_nodes == cn_maximum && best_knap_num < MAXCAN)
						{
							best_knap_par1[best_knap_num] = k1;
							best_knap_par2[best_knap_num] = k2;
							best_knap_num++;
						}
					}
				}
			}
		}
		
		if (cn_maximum > 0)
		{
			int rx = rand() % best_knap_num;
			int knap_par1 = best_knap_par1[rx];
			int knap_par2 = best_knap_par2[rx];
			for (int j = 0; j < Instance.items; j++)			//preserve the backbone 
			{
				if (Pop_sol[parent1][j] == knap_par1 && Pop_sol[parent2][j] == knap_par2 && offspring[j] == Instance.knaps && satisfy_all_constraints_relocate(j, offspring[j], len_knap))
				{
					update_global_data(j, offspring[j], len_knap);
					offspring[j] = len_knap;
					flag_knap1[knap_par1] = 0;
					flag_knap2[knap_par2] = 0;
				}
			}
			//cout << "cn_maximum=" << cn_maximum << endl;
			//assign some vertices in an alternating way
			int flag_par1 = 1;
			int flag_par2 = 1;
			int step = 1;
			while (flag_par1 || flag_par2)
			{
				//cout << "step=" << step << endl;
				int tar_par = -1;
				int tar_knap = -1;
				int tar_item = -1;
				if (step % 2 == 0)		//select vertices from parent1
				{
					tar_par = parent1;
					tar_knap = knap_par1;
					flag_par1 = 0;
				}
				else
				{
					tar_par = parent2;
					tar_knap = knap_par2;
					flag_par2 = 0;

				}
				double vd_maximum = -MAXNUM;
				for (int j = 0; j < Instance.items; j++)
				{
					if (Pop_sol[tar_par][j] == tar_knap && offspring[j] == Instance.knaps && satisfy_all_constraints_relocate(j, offspring[j], len_knap))
					{
						double vd = value_density_item(j, len_knap, offspring);
						if (vd > vd_maximum + PRECISION)
						{
							vd_maximum = vd;
							tar_item = j;
							if (step % 2 == 0)
								flag_par1 = 1;
							else
								flag_par2 = 1;
						}
					}
				}
				if ((flag_par1 && (step % 2 == 0)) || (flag_par2 && (step % 2 != 0)))
				{
					update_global_data(tar_item, offspring[tar_item], len_knap);
					offspring[tar_item] = len_knap;
				}
				step++;
			}
		}
		len_knap++;
	}

	//second step: complete the partial solution in a greedy manner
	len_knap = 0;
	for (int k = 0; k < Instance.knaps; k++)
	{
		address_knaps[k] = len_knap;
		can_knaps[len_knap] = k;
		len_knap++;
	}
	while (len_knap > 0)
	{
		greedy_helpfunc_value_density(offspring, can_item, can_knaps, address_knaps, len_knap);
	}

	delete[]can_item; can_item = NULL;
	delete[]can_knaps; can_knaps = NULL;
	delete[]address_knaps; address_knaps = NULL;

	delete[]flag_knap1; flag_knap1 = NULL;
	delete[]flag_knap2; flag_knap2 = NULL;
}

//update the population, update the worst solution in terms of the objective value
void update_population(int *offspring, double f_off)
{
	double obj_worst = MAXNUM;
	int index_worst = -1;
	for (int i = 0; i < Pop_size; i++)
	{
		if (Pop_obj[i] < obj_worst - PRECISION)
		{
			obj_worst = Pop_obj[i];
			index_worst = i;
		}
	}
	if (f_off > obj_worst + PRECISION)
	{
		bool is_exist = is_same_solution(Pop_sol[index_worst], offspring);
		if (!is_exist)
		{
			memcpy(Pop_sol[index_worst], offspring, sizeof(int)*Instance.items);
			Pop_obj[index_worst] = f_off;
		}
	}
}

