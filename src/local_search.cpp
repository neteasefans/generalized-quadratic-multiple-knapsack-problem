#include"func_state.h"
#include"global_variables.h"
#include <iostream>
#include<string.h>
#include <time.h>
#include<math.h>
using namespace std;


double descent_local_search(int *sol, double &f_sol)
{
	int *rand_item = new int[Instance.items];
	int *rand_knap = new int[Instance.knaps + 1];
	int *flag_item = new int[Instance.items];
	int *flag_knap = new int[Instance.knaps + 1];
	int flag_improve = 1;
	//while (flag_improve && 1.0*(clock() - Start_time) / CLOCKS_PER_SEC < Time_limit)
	while (flag_improve)
	{		
		flag_improve = 0;
		for (int j = 0; j < Instance.items; j++)
			flag_item[j] = 0;
		for (int k = 0; k < Instance.knaps + 1; k++)
			flag_knap[k] = 0;
		int len1 = 0, len2 = 0;
		while (len1 < Instance.items)
		{
			int r_item = rand() % Instance.items;
			if (flag_item[r_item] == 0)
			{
				rand_item[len1++] = r_item;
				flag_item[r_item] = 1;
			}
		}
		while (len2 < Instance.knaps + 1)
		{
			int r_k = rand() % (Instance.knaps + 1);
			if (flag_knap[r_k] == 0)
			{
				rand_knap[len2++] = r_k;
				flag_knap[r_k] = 1;
			}
		}
		
		//for relocate operator
		for (int j = 0; j < Instance.items; j++)
		{
			int item = rand_item[j];
			for (int k = 0; k < Instance.knaps + 1; k++)
			{
				int knap = rand_knap[k];
				if (sol[item] != knap && satisfy_all_constraints_relocate(item, sol[item], knap))
				{
					double delta = compute_delta_relocate(sol, item, knap);
					if (delta > 0 + PRECISION)
					{
						relocate_move(sol, item, knap);
						flag_improve = 1;
						f_sol += delta;
						if (f_sol > Best_sol_one_run_obj + PRECISION)
							update_best_sol_one_run(f_sol, sol);
					}
				}
			}			
		}
		
		
		//for swap operator
		for (int j1 = 0; j1 < Instance.items; j1++)
		{
			int item_j1 = rand_item[j1];
			for (int j2 = j1 + 1; j2 < Instance.items; j2++)
			{
				int item_j2 = rand_item[j2];
				if (sol[item_j1] != sol[item_j2] && satisfy_all_constraints_swap(item_j1, sol[item_j1], item_j2, sol[item_j2]))
				{
					double delta = compute_delta_swap(sol, item_j1, item_j2);
					if (delta > 0 + PRECISION)
					{
						swap_move(sol, item_j1, item_j2);
						flag_improve = 1;
						f_sol += delta;
						if (f_sol > Best_sol_one_run_obj + PRECISION)
							update_best_sol_one_run(f_sol, sol);
					}
				}
			}			
		}	

		/*
		//for 2-1 exchange operator
		for (int j1 = 0; j1 < Instance.items; j1++)
		{
			int item_j1 = rand_item[j1];
			for (int j2 = j1 + 1; j2 < Instance.items; j2++)
			{
				int item_j2 = rand_item[j2];
				for (int j3 = 0; j3 < Instance.items; j3++)
				{
					int item_j3 = rand_item[j3];
					if (sol[item_j1] == sol[item_j2] && sol[item_j1] != sol[item_j3] && satisfy_all_constraints_2_1_exchange(item_j1, item_j2, sol[item_j1], item_j3, sol[item_j3]))
					{
						double delta = compute_delta_2_1_exchange(sol, item_j1, item_j2, item_j3);
						if (delta > 0 + PRECISION)
						{
							two_one_exchange_move(sol, item_j1, item_j2, item_j3);
							flag_improve = 1;
							f_sol += delta;
							if (f_sol > Best_sol_one_run_obj + PRECISION)
								update_best_sol_one_run(f_sol, sol);
						}
					}
				}
			}
		}*/

	}
	delete[]rand_item; rand_item = NULL;
	delete[]rand_knap; rand_knap = NULL;
	delete[]flag_item; flag_item = NULL;
	delete[]flag_knap; flag_knap = NULL;
	return f_sol;
}


double feasible_tabu_search(int *sol, double &f_sol)
{	
	int best_re_item[MAXCAN], best_re_knap[MAXCAN];
	int tabu_best_re_item[MAXCAN], tabu_best_re_knap[MAXCAN];
	int best_swap_item1[MAXCAN], best_swap_item2[MAXCAN];
	int tabu_best_swap_item1[MAXCAN], tabu_best_swap_item2[MAXCAN];
	//int best_2_1_exchange_item1[MAXCAN], best_2_1_exchange_item2[MAXCAN], best_2_1_exchange_item3[MAXCAN];
	//int tabu_best_2_1_exchange_item1[MAXCAN], tabu_best_2_1_exchange_item2[MAXCAN], tabu_best_2_1_exchange_item3[MAXCAN];

	int best_re_num, tabu_best_re_num, best_swap_num, tabu_best_swap_num, best_2_1_exchange_num, tabu_best_2_1_exchange_num;
	double delta, delta_re_best, tabu_delta_re_best, delta_swap_best, tabu_delta_swap_best, delta_2_1_exchange_best, tabu_delta_2_1_exchange_best;
	
	for (int j = 0; j < Instance.items; j++)
		for (int k = 0; k < Instance.knaps + 1; k++)
			Tabu_list[j][k] = 0;
	double f_local_best = f_sol;
	for (int j = 0; j < Instance.items; j++)
		Local_best_sol[j] = sol[j];

	int non_improve = 0, iter = 0;
	//while (non_improve < Tabu_depth && 1.0*(clock() - Start_time) / CLOCKS_PER_SEC < Time_limit)
	while (non_improve < Tabu_depth)
	{
		delta_re_best = -MAXNUM;
		tabu_delta_re_best = -MAXNUM;
		delta_swap_best = -MAXNUM;
		tabu_delta_swap_best = -MAXNUM;
		delta_2_1_exchange_best = -MAXNUM;
		tabu_delta_2_1_exchange_best = -MAXNUM;
		best_re_num = 0;
		tabu_best_re_num = 0;
		best_swap_num = 0;
		tabu_best_swap_num = 0;
		best_2_1_exchange_num = 0;
		tabu_best_2_1_exchange_num = 0;

		//for relocate operator
		for (int j = 0; j < Instance.items; j++)
		{
			for (int k = 0; k < Instance.knaps + 1; k++)
			{
				if (sol[j] != k && satisfy_all_constraints_relocate(j, sol[j], k))
				{					
					delta = compute_delta_relocate(sol, j, k);
					if (Tabu_list[j][k] <= iter)
					{
						if (delta > delta_re_best + PRECISION)
						{
							best_re_num = 0;
							best_re_item[best_re_num] = j;
							best_re_knap[best_re_num] = k;
							delta_re_best = delta;
							best_re_num++;
						}
						else if (fabs(delta - delta_re_best) <= PRECISION && best_re_num < MAXCAN)
						{
							best_re_item[best_re_num] = j;
							best_re_knap[best_re_num] = k;
							best_re_num++;
						}
					}
					else
					{
						if (delta > tabu_delta_re_best + PRECISION)
						{
							tabu_best_re_num = 0;
							tabu_best_re_item[tabu_best_re_num] = j;
							tabu_best_re_knap[tabu_best_re_num] = k;
							tabu_delta_re_best = delta;
							tabu_best_re_num++;
						}
						else if (fabs(delta - tabu_delta_re_best) <= PRECISION && tabu_best_re_num < MAXCAN)
						{
							tabu_best_re_item[tabu_best_re_num] = j;
							tabu_best_re_knap[tabu_best_re_num] = k;
							tabu_best_re_num++;
						}
					}
				}
			}
		}
				
		//for swap operator
		for (int j = 0; j < Instance.items; j++)
		{
			for (int j2 = j + 1; j2 < Instance.items; j2++)
			{
				if (sol[j] != sol[j2] && satisfy_all_constraints_swap(j, sol[j], j2, sol[j2]))
				{					
					delta = compute_delta_swap(sol, j, j2);
					if (Tabu_list[j][sol[j2]] <= iter && Tabu_list[j2][sol[j]] <= iter)
					{
						if (delta > delta_swap_best + PRECISION)
						{
							best_swap_num = 0;
							best_swap_item1[best_swap_num] = j;
							best_swap_item2[best_swap_num] = j2;
							delta_swap_best = delta;
							best_swap_num++;
						}
						else if (fabs(delta - delta_swap_best) <= PRECISION && best_swap_num < MAXCAN)
						{
							best_swap_item1[best_swap_num] = j;
							best_swap_item2[best_swap_num] = j2;
							best_swap_num++;
						}
					}
					else
					{
						if (delta > tabu_delta_swap_best + PRECISION)
						{
							tabu_best_swap_num = 0;
							tabu_best_swap_item1[tabu_best_swap_num] = j;
							tabu_best_swap_item2[tabu_best_swap_num] = j2;
							tabu_delta_swap_best = delta;
							tabu_best_swap_num++;
						}
						else if (fabs(delta - tabu_delta_swap_best) <= PRECISION && tabu_best_swap_num < MAXCAN)
						{
							tabu_best_swap_item1[tabu_best_swap_num] = j;
							tabu_best_swap_item2[tabu_best_swap_num] = j2;
							tabu_best_swap_num++;
						}
					}
				}
			}
		}

		/*
		//for 2-1 exchange operator
		for (int j = 0; j < Instance.items; j++)
		{
			for (int j2 = j + 1; j2 < Instance.items; j2++)
			{
				for (int j3 = 0; j3 < Instance.items; j3++)
				{
					if (sol[j] == sol[j2] && sol[j] != sol[j3] && satisfy_all_constraints_2_1_exchange(j, j2, sol[j], j3, sol[j3]))
					{
						delta = compute_delta_2_1_exchange(sol, j, j2, j3);
						if (Tabu_list[j][sol[j3]] <= iter && Tabu_list[j2][sol[j3]] <= iter && Tabu_list[j3][sol[j]] <= iter)
						{
							if (delta > delta_2_1_exchange_best + PRECISION)
							{
								best_2_1_exchange_num = 0;
								best_2_1_exchange_item1[best_2_1_exchange_num] = j;
								best_2_1_exchange_item2[best_2_1_exchange_num] = j2;
								best_2_1_exchange_item3[best_2_1_exchange_num] = j3;
								delta_2_1_exchange_best = delta;
								best_2_1_exchange_num++;
							}
							else if (fabs(delta - delta_2_1_exchange_best) <= PRECISION && best_2_1_exchange_num < MAXCAN)
							{
								best_2_1_exchange_item1[best_2_1_exchange_num] = j;
								best_2_1_exchange_item2[best_2_1_exchange_num] = j2;
								best_2_1_exchange_item3[best_2_1_exchange_num] = j3;
								best_2_1_exchange_num++;
							}
						}
						else
						{
							if (delta > tabu_delta_2_1_exchange_best + PRECISION)
							{
								tabu_best_2_1_exchange_num = 0;
								tabu_best_2_1_exchange_item1[tabu_best_2_1_exchange_num] = j;
								tabu_best_2_1_exchange_item2[tabu_best_2_1_exchange_num] = j2;
								tabu_best_2_1_exchange_item3[tabu_best_2_1_exchange_num] = j3;							
								tabu_delta_2_1_exchange_best = delta;
								tabu_best_2_1_exchange_num++;
							}
							else if (fabs(delta - tabu_delta_2_1_exchange_best) <= PRECISION && tabu_best_2_1_exchange_num < MAXCAN)
							{
								tabu_best_2_1_exchange_item1[tabu_best_2_1_exchange_num] = j;
								tabu_best_2_1_exchange_item2[tabu_best_2_1_exchange_num] = j2;
								tabu_best_2_1_exchange_item3[tabu_best_2_1_exchange_num] = j3;
								tabu_best_2_1_exchange_num++;
							}
						}
					}
				}
			}
		}*/

		int move_type = -1;
		double delta_maximum = -MAXNUM;
		if (best_re_num > 0 && delta_re_best > delta_maximum + PRECISION)
		{
			move_type = 1;
			delta_maximum = delta_re_best;
		}
		if ((tabu_best_re_num > 0 && tabu_delta_re_best > delta_maximum + PRECISION && f_sol + tabu_delta_re_best > f_local_best + PRECISION)
			|| (best_re_num == 0 && tabu_best_re_num > 0))
		{
			move_type = 2;
			delta_maximum = tabu_delta_re_best;
		}
		if (best_swap_num > 0 && delta_swap_best > delta_maximum + PRECISION)
		{
			move_type = 3;
			delta_maximum = delta_swap_best;
		}
		if ((tabu_best_swap_num > 0 && tabu_delta_swap_best > delta_maximum + PRECISION && f_sol + tabu_delta_swap_best > f_local_best + PRECISION) 
			|| (best_swap_num == 0 && tabu_best_swap_num > 0))
		{
			move_type = 4;
			delta_maximum = tabu_delta_swap_best;
		}
		if (best_2_1_exchange_num > 0 && delta_2_1_exchange_best > delta_maximum + PRECISION)
		{
			move_type = 5;
			delta_maximum = delta_2_1_exchange_best;
		}
		if ((tabu_best_2_1_exchange_num > 0 && tabu_delta_2_1_exchange_best > delta_maximum + PRECISION && f_sol + tabu_delta_2_1_exchange_best > f_local_best + PRECISION)
			|| (best_2_1_exchange_num == 0 && tabu_best_2_1_exchange_num > 0))
		{
			move_type = 6;
			delta_maximum = tabu_delta_2_1_exchange_best;
		}
		int item = -1, knap = -1;
		int item1 = -1, item2 = -1, item3 = -1;
		int rx = -1;
		switch(move_type)
		{
			case 1:
				rx = rand() % best_re_num;
				item = best_re_item[rx];
				knap = best_re_knap[rx];
				break;
			case 2:
				rx = rand() % tabu_best_re_num;
				item = tabu_best_re_item[rx];
				knap = tabu_best_re_knap[rx];
				break;
			case 3:
				rx = rand() % best_swap_num;
				item1 = best_swap_item1[rx];
				item2 = best_swap_item2[rx];
				break;
			case 4:
				rx = rand() % tabu_best_swap_num;
				item1 = tabu_best_swap_item1[rx];
				item2 = tabu_best_swap_item2[rx];
				break;
			/*case 5:
				rx = rand() % best_2_1_exchange_num;
				item1 = best_2_1_exchange_item1[rx];
				item2 = best_2_1_exchange_item2[rx];
				item3 = best_2_1_exchange_item3[rx];
				break;
			case 6:
				rx = rand() % tabu_best_2_1_exchange_num;
				item1 = tabu_best_2_1_exchange_item1[rx];
				item2 = tabu_best_2_1_exchange_item2[rx];
				item3 = tabu_best_2_1_exchange_item3[rx];
				break;*/
		}
	
		if (move_type == 1 || move_type == 2)
		{
			relocate_move(sol, item, knap);			
			Tabu_list[item][knap] = Tabu_tenure + iter;			
		}
		else if (move_type == 3 || move_type == 4)
		{
			int knap1 = sol[item1];
			int knap2 = sol[item2];
			swap_move(sol, item1, item2);
			Tabu_list[item1][knap1] = Tabu_tenure + iter;
			Tabu_list[item2][knap2] = Tabu_tenure + iter;			
		}
		else if (move_type == 5 || move_type == 6)
		{
			int knap1 = sol[item1];
			int knap2 = sol[item3];
			two_one_exchange_move(sol, item1, item2, item3);
			Tabu_list[item1][knap1] = Tabu_tenure + iter;
			Tabu_list[item2][knap1] = Tabu_tenure + iter;
			Tabu_list[item3][knap2] = Tabu_tenure + iter;
		}
		if (move_type != -1)
		{
			f_sol += delta_maximum;
			if (f_sol > f_local_best + PRECISION)
			{
				non_improve = 0;
				for (int i = 0; i < Instance.items; i++)
					Local_best_sol[i] = sol[i];
				f_local_best = f_sol;
			}
			else
				non_improve++;
			if (f_sol > Best_sol_one_run_obj + PRECISION)
				update_best_sol_one_run(f_sol, sol);
		}	
		iter++;		
		//if (iter % 1000 == 0)
		//cout << "iter=" << iter << ", move_type=" << move_type << ", f_sol=" << f_sol << ", f_loca_best=" << f_local_best << ", Best_sol_one_run_obj=" << Best_sol_one_run_obj << endl;
		//proof(sol, f_sol);
		//check_move(sol, f_sol, 0.0);
	}
	for (int j = 0; j < Instance.items; j++)				//return the best local solution
		sol[j] = Local_best_sol[j];
	f_sol = f_local_best;
	build_global_data(sol);	
	build_delta_matrix(sol);
	return f_sol;
}

//guided by the extend evaluation function F(s) = f(s) - penelay_factor * g1(s) and penelay_factor is adjusted adaptively
double infeasible_tabu_search(int *sol, double &f_sol)
{
	int best_re_item[MAXCAN], best_re_knap[MAXCAN];	
	int tabu_best_re_item[MAXCAN], tabu_best_re_knap[MAXCAN];
	int best_swap_item1[MAXCAN], best_swap_item2[MAXCAN];
	int tabu_best_swap_item1[MAXCAN], tabu_best_swap_item2[MAXCAN];

	int inf_count = 0;

	double best_re_f_s[MAXCAN], best_re_g_s[MAXCAN];
	double tabu_best_re_f_s[MAXCAN], tabu_best_re_g_s[MAXCAN];
	double best_swap_f_s[MAXCAN], best_swap_g_s[MAXCAN];
	double tabu_best_swap_f_s[MAXCAN], tabu_best_swap_g_s[MAXCAN];

	int best_re_num, tabu_best_re_num, best_swap_num, tabu_best_swap_num;
	double delta_f_s, delta_g_s, delta_fg_s;
	double delta_re_best, tabu_delta_re_best, delta_swap_best, tabu_delta_swap_best;

	for (int j = 0; j < Instance.items; j++)
		for (int k = 0; k < Instance.knaps + 1; k++)
			Tabu_list[j][k] = 0;	
	double g_sol = 0.0;
	double fg_sol = f_sol - Pe_factor *g_sol;
	double fg_local_best = fg_sol;
	for (int j = 0; j < Instance.items; j++)
		Local_best_sol[j] = sol[j];
	for (int j = 0; j < Instance.items; j++)
		Global_best_sol[j] = sol[j];				//perserve the best feasible solution in this procedure
	double f_global_best = f_sol;
	int non_improve = 0, iter = 0;
	//while (non_improve < Tabu_depth && 1.0*(clock() - Start_time) / CLOCKS_PER_SEC < Time_limit)
	while (non_improve < Tabu_depth)
	{
		delta_re_best = -MAXNUM;
		tabu_delta_re_best = -MAXNUM;
		delta_swap_best = -MAXNUM;
		tabu_delta_swap_best = -MAXNUM;
		best_re_num = 0;
		tabu_best_re_num = 0;
		best_swap_num = 0;
		tabu_best_swap_num = 0;
		
		
		//for relocate operator
		for (int j = 0; j < Instance.items; j++)
		{
			for (int k = 0; k < Instance.knaps + 1; k++)
			{
				if (sol[j] != k && satisfy_nr_constraints_relocate(j, sol[j], k))			//relax constraint of capacity
				//if (sol[j] != k && satisfy_cap_constraints_relocate(j, sol[j], k))		//relax constraint of class
				//if (sol[j] != k)															//relax constraints of capacity and class
				{
					delta_f_s = compute_delta_relocate(sol, j, k);
					delta_g_s = compute_delta_relocate_g1(sol, j, k);
					//delta_g_s = compute_delta_relocate_g2(sol, j, k);
					//delta_g_s = compute_delta_relocate_g1(sol, j, k) + compute_delta_relocate_g2(sol, j, k);
					delta_fg_s = delta_f_s - Pe_factor*delta_g_s;
					if (Tabu_list[j][k] <= iter)
					{
						if (delta_fg_s > delta_re_best + PRECISION)
						{
							best_re_num = 0;
							best_re_item[best_re_num] = j;
							best_re_knap[best_re_num] = k;
							best_re_f_s[best_re_num] = delta_f_s;
							best_re_g_s[best_re_num] = delta_g_s;
							delta_re_best = delta_fg_s;							
							best_re_num++;
						}
						else if (fabs(delta_fg_s - delta_re_best) <= PRECISION && best_re_num < MAXCAN)
						{
							best_re_item[best_re_num] = j;
							best_re_knap[best_re_num] = k;
							best_re_f_s[best_re_num] = delta_f_s;
							best_re_g_s[best_re_num] = delta_g_s;
							best_re_num++;
						}
					}
					else
					{
						if (delta_fg_s > tabu_delta_re_best + PRECISION)
						{
							tabu_best_re_num = 0;
							tabu_best_re_item[tabu_best_re_num] = j;
							tabu_best_re_knap[tabu_best_re_num] = k;
							tabu_best_re_f_s[tabu_best_re_num] = delta_f_s;
							tabu_best_re_g_s[tabu_best_re_num] = delta_g_s;
							tabu_delta_re_best = delta_fg_s;						
							tabu_best_re_num++;
						}
						else if (fabs(delta_fg_s - tabu_delta_re_best) <= PRECISION && tabu_best_re_num < MAXCAN)
						{
							tabu_best_re_item[tabu_best_re_num] = j;
							tabu_best_re_knap[tabu_best_re_num] = k;
							tabu_best_re_f_s[tabu_best_re_num] = delta_f_s;
							tabu_best_re_g_s[tabu_best_re_num] = delta_g_s;
							tabu_best_re_num++;
						}
					}
				}
			}
		}				
		
		
		//for swap operator
		for (int j = 0; j < Instance.items; j++)
		{
			for (int j2 = j + 1; j2 < Instance.items; j2++)
			{
				if (sol[j] != sol[j2] && satisfy_nr_constraints_swap(j, sol[j], j2, sol[j2]))	
				//if (sol[j] != sol[j2] && satisfy_cap_constraints_swap(j, sol[j], j2, sol[j2]))	
				//if (sol[j] != sol[j2])
				{
					delta_f_s = compute_delta_swap(sol, j, j2);
					delta_g_s = compute_delta_swap_g1(sol, j, j2);
					//delta_g_s = compute_delta_swap_g2(sol, j, j2);
					//delta_g_s = compute_delta_swap_g1(sol, j, j2) + compute_delta_swap_g2(sol, j, j2);
					delta_fg_s = delta_f_s - Pe_factor*delta_g_s;
					if (Tabu_list[j][sol[j2]] <= iter && Tabu_list[j2][sol[j]] <= iter)
					{
						if (delta_fg_s > delta_swap_best + PRECISION)
						{							
							best_swap_num = 0;
							best_swap_item1[best_swap_num] = j;
							best_swap_item2[best_swap_num] = j2;
							best_swap_f_s[best_swap_num] = delta_f_s;
							best_swap_g_s[best_swap_num] = delta_g_s;
							delta_swap_best = delta_fg_s;
							best_swap_num++;
						}
						else if (fabs(delta_fg_s - delta_swap_best) <= PRECISION && best_swap_num < MAXCAN)
						{
							best_swap_item1[best_swap_num] = j;
							best_swap_item2[best_swap_num] = j2;
							best_swap_f_s[best_swap_num] = delta_f_s;
							best_swap_g_s[best_swap_num] = delta_g_s;
							best_swap_num++;
						}
					}
					else
					{
						if (delta_fg_s > tabu_delta_swap_best + PRECISION)
						{
							tabu_best_swap_num = 0;
							tabu_best_swap_item1[tabu_best_swap_num] = j;
							tabu_best_swap_item2[tabu_best_swap_num] = j2;
							tabu_best_swap_f_s[tabu_best_swap_num] = delta_f_s;
							tabu_best_swap_g_s[tabu_best_swap_num] = delta_g_s;
							tabu_delta_swap_best = delta_fg_s;
							tabu_best_swap_num++;
						}
						else if (fabs(delta_fg_s - tabu_delta_swap_best) <= PRECISION && tabu_best_swap_num < MAXCAN)
						{
							tabu_best_swap_item1[tabu_best_swap_num] = j;
							tabu_best_swap_item2[tabu_best_swap_num] = j2;
							tabu_best_swap_f_s[tabu_best_swap_num] = delta_f_s;
							tabu_best_swap_g_s[tabu_best_swap_num] = delta_g_s;
							tabu_best_swap_num++;
						}
					}
				}
			}
		}
		
		int move_type = -1;
		double delta_maximum = -MAXNUM;		
				
		if (best_re_num > 0 && delta_re_best > delta_maximum + PRECISION)
		{
			move_type = 1;
			delta_maximum = delta_re_best;		
		}
		if ((tabu_best_re_num > 0 && tabu_delta_re_best > delta_maximum + PRECISION && fg_sol + tabu_delta_re_best > fg_local_best + PRECISION)
			|| (best_re_num == 0 && tabu_best_re_num > 0))
		{
			move_type = 2;
			delta_maximum = tabu_delta_re_best;
		}
		
		if (best_swap_num > 0 && delta_swap_best > delta_maximum + PRECISION)
		{
			move_type = 3;
			delta_maximum = delta_swap_best;			
		}
		if ((tabu_best_swap_num > 0 && tabu_delta_swap_best > delta_maximum + PRECISION && fg_sol + tabu_delta_swap_best > fg_local_best + PRECISION)
			|| (best_swap_num == 0 && tabu_best_swap_num > 0))
		{
			move_type = 4;
			delta_maximum = tabu_delta_swap_best;			
		}

		int item = -1, knap = -1;
		int item1 = -1, item2 = -1;
		int rx = -1;
		double delta_fs_sel = -MAXNUM;
		double delta_gs_sel = -MAXNUM;
		switch (move_type)
		{
		case 1:
			rx = rand() % best_re_num;
			item = best_re_item[rx];
			knap = best_re_knap[rx];		
			delta_fs_sel = best_re_f_s[rx];
			delta_gs_sel = best_re_g_s[rx];
			break;
		case 2:
			rx = rand() % tabu_best_re_num;
			item = tabu_best_re_item[rx];
			knap = tabu_best_re_knap[rx];
			delta_fs_sel = tabu_best_re_f_s[rx];
			delta_gs_sel = tabu_best_re_g_s[rx];
			break;
		case 3:
			rx = rand() % best_swap_num;
			item1 = best_swap_item1[rx];
			item2 = best_swap_item2[rx];
			delta_fs_sel = best_swap_f_s[rx];
			delta_gs_sel = best_swap_g_s[rx];
			break;
		case 4:
			rx = rand() % tabu_best_swap_num;
			item1 = tabu_best_swap_item1[rx];
			item2 = tabu_best_swap_item2[rx];
			delta_fs_sel = tabu_best_swap_f_s[rx];
			delta_gs_sel = tabu_best_swap_g_s[rx];
			break;
		}

		if (move_type == 1 || move_type == 2)
		{
			relocate_move(sol, item, knap);			
			Tabu_list[item][knap] = Tabu_tenure + iter;
		}
		else if (move_type == 3 || move_type == 4)
		{
			int knap1 = sol[item1];
			int knap2 = sol[item2];
			swap_move(sol, item1, item2);
			Tabu_list[item1][knap1] = Tabu_tenure + iter;
			Tabu_list[item2][knap2] = Tabu_tenure + iter;
		}
		if (move_type != -1)
		{
			f_sol += delta_fs_sel;
			g_sol += delta_gs_sel;
			fg_sol = f_sol - Pe_factor*g_sol;
			if (fg_sol > fg_local_best + PRECISION)
			{
				non_improve = 0;
				for (int i = 0; i < Instance.items; i++)
					Local_best_sol[i] = sol[i];
				fg_local_best = fg_sol;
			}
			else
				non_improve++;
			if (f_sol > Best_sol_one_run_obj + PRECISION && fabs(g_sol - 0) <= PRECISION)			
				update_best_sol_one_run(f_sol, sol);			
			if (f_sol > f_global_best + PRECISION && fabs(g_sol - 0) <= PRECISION)
			{
				memcpy(Global_best_sol, sol, sizeof(int)*Instance.items);
				f_global_best = f_sol;
			}
		}						
		iter++;
		//if (iter % 1000 == 0)
		//cout << "iter=" << iter << ", move_type=" << move_type << ", fg_sol=" << fg_sol << ", fg_loca_best=" << fg_local_best << ", f_sol=" << f_sol <<
		//	", g_sol=" << g_sol << ", Best_sol_one_run_obj=" << Best_sol_one_run_obj << ", Pe_factor=" << Pe_factor << endl;
		//proof(sol, f_sol);
		//check_move(sol, f_sol, g_sol);
		if (fabs(g_sol - 0) < PRECISION && fabs(fg_sol - f_sol) > 0 + PRECISION)
		{
			cout << "in infeasible_tabu_search func, error, g_sol=" << g_sol << ", f_sol=" << f_sol << ", but fg_sol=" << fg_sol << endl;
			getchar();
		}
		if (g_sol > 0 + PRECISION)
			inf_count++;
		if (iter % Lambda == 0)
		{
			if (inf_count == Lambda)
				Pe_factor *= Incre;
			if (inf_count == 0)
				Pe_factor /= Incre;
			if (Pe_factor < 1)
				Pe_factor = 1;		
			inf_count = 0;
		}
	
	}
	for (int j = 0; j < Instance.items; j++)				//return the best local solution
		sol[j] = Global_best_sol[j];
	f_sol = f_global_best;
	build_global_data(sol);
	build_delta_matrix(sol);
	return f_sol;
}


