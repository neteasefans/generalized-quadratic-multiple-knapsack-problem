#include"func_state.h"
#include"global_variables.h"
#include<time.h>
#include<string.h>
#include<math.h>
#include <iostream>

#include <fstream>
using namespace std;


char *Instance_name;
instance_data Instance;
int**Size_kr;		//Size_kr[k][r] indiactes the number of items in knapsack k with the class r, k \in {0, 1,..., Instance.knaps}
int*Size_r;			//Size_r[r] indicates the number of knapsacks that involves the item with the class r, and excludes the knapsack labeling Instance.knaps
double*Capacity1_k;
double*Capacity2_k;
double*Exceed_cap_k;
int *Exceed_class_r;
double**Delta_matrix_jk;
int Num_greedy_ini;

//for infeasible tabu search
int Pe_factor;			//penalty factor
int Lambda;				//every Lambda consecutive iterations 
int Incre;

//for tabu search
int Tabu_depth;
int Tabu_tenure;
int**Tabu_list;

//for perturbation
neighborhood*Neighbors;
double Pert_str;			//perturbation strength


//for reinforcement learning
double RL_omega;
double RL_alpha;
double RL_beta;
double RL_gamma;
double RL_smoothing_coeff;
double RL_proba_threshold;
double **Proba_matrix_jk;	//probability learning matrix

//for memetic search
int Pop_size;				//population size
int **Pop_sol;
double*Pop_obj;

int *Cur_sol;
int *Local_best_sol;
int *Best_sol_one_run;
int *Global_best_sol;
double Best_sol_one_run_obj;



double Time_limit, Start_time, Run_time;
int Max_generations;
//read instance
void read_instance()
{
	ifstream FIC;
	FIC.open(Instance_name);
	if (FIC.fail())
	{
		cout << "### Error open, Instance_name " << Instance_name << endl;
		exit(0);
	}
	char str_reading[MAXCAN];
	char str_last[MAXCAN];
	char str_tmp1[MAXCAN], str_tmp2[MAXCAN];
	FIC >> str_reading;	
	while (1)
	{	
		strcpy(str_last, str_reading);
		FIC >> str_reading;
		if ((strcmp(str_last, "siparis") == 0 && strcmp(str_reading, "turu") == 0) || 
			(strcmp(str_last, "knapsack") == 0 && strcmp(str_reading, "indisi") == 0) ||
			(strcmp(str_last, "kalip") == 0 && strcmp(str_reading, "indisi") == 0) ||
			(strcmp(str_last, "sayi") == 0 && strcmp(str_reading, "scalar") == 0))				//for items, knapsacks, classes  and U
		{
			strcpy(str_tmp1, str_last);
			strcpy(str_tmp2, str_reading);
			FIC >> str_reading;
			int data = 0;
			int num_len = 0;
			for (int i = strlen(str_reading) - 1; i >= 0; i--)
				if (str_reading[i] >= '0' && str_reading[i] <= '9')
					data += (str_reading[i] - '0') * int(pow(10, num_len++));
				else if (str_reading[i] == '*')
					break;				
		
			if (strcmp(str_tmp1, "siparis") == 0 && strcmp(str_tmp2, "turu") == 0)
			{
				Instance.items = data;
				Instance.we = new double[Instance.items];
				Instance.po = new double[Instance.items];
				for (int i = 0; i < Instance.items; i++)
				{
					Instance.we[i] = 0.0;
					Instance.po[i] = 0.0;
				}
				Instance.q_ij = new double*[Instance.items];
				for (int i = 0; i < Instance.items; i++)
					Instance.q_ij[i] = new double[Instance.items];
				for (int i = 0; i < Instance.items; i++)
					for (int j = 0; j < Instance.items; j++)
						Instance.q_ij[i][j] = 0.0;
			}
			else if (strcmp(str_tmp1, "knapsack") == 0 && strcmp(str_tmp2, "indisi") == 0)
			{
				Instance.knaps = data;
				Instance.p_jk = new double*[Instance.items];
				for (int i = 0; i < Instance.items; i++)
					Instance.p_jk[i] = new double[Instance.knaps];
			}
			else if (strcmp(str_tmp1, "kalip") == 0 && strcmp(str_tmp2, "indisi") == 0)
			{
				Instance.classes = data;
				Instance.t_rj = new int*[Instance.classes];
				for (int i = 0; i < Instance.classes; i++)
					Instance.t_rj[i] = new int[Instance.items];
				for (int i = 0; i < Instance.classes; i++)
					for (int j = 0; j < Instance.items; j++)
						Instance.t_rj[i][j] = 0;
				Instance.t_jr = new int[Instance.items];
				for (int i = 0; i < Instance.items; i++)
					Instance.t_jr[i] = -1;
				Instance.psi = new double*[Instance.classes];
				for (int i = 0; i < Instance.classes; i++)
					Instance.psi[i] = new double[Instance.knaps];
				for (int i = 0; i < Instance.classes; i++)
					for (int j = 0; j < Instance.knaps; j++)
						Instance.psi[i][j] = 0.0;
				Instance.sr = new double[Instance.classes];
				Instance.nr = new int[Instance.classes];
				for (int i = 0; i < Instance.classes; i++)
				{
					Instance.sr[i] = 0.0;
					Instance.nr[i] = 0;
				}
			}
			else if (strcmp(str_tmp1, "sayi") == 0 && strcmp(str_tmp2, "scalar") == 0)
				Instance.u = data;			
			cout<<"Instance.items="<<Instance.items<<endl;
			cout<<"Instance.knapsacks="<<Instance.knapsacks<<endl;
			cout<<"Instance.classes="<<Instance.classes<<endl;
			cout<<"Instance.U="<<Instance.u<<endl;
		}		
		else if ((strcmp(str_last, "parameter") == 0 && strcmp(str_reading, "w(j)/") == 0) ||
			(strcmp(str_last, "parameter") == 0 && strcmp(str_reading, "po(j)/") == 0) ||
			(strcmp(str_last, "parameter") == 0 && strcmp(str_reading, "s(r)/") == 0) || 
			(strcmp(str_last, "parameter") == 0 && strcmp(str_reading, "nr(r)/") == 0))					//for we, po, sr and nr
		{
			strcpy(str_tmp1, str_last);
			strcpy(str_tmp2, str_reading);
			while (1)
			{
				FIC >> str_last >> str_reading;
				if (strcmp(str_tmp1, "parameter") == 0 && strcmp(str_tmp2, "po(j)/") == 0)
				{
					if (strcmp(str_last, "/") == 0)
						break;
				}
				else
					if (strcmp(str_last, "/;") == 0)
						break;
				int item = atoi(str_last);				
				if (strcmp(str_tmp1, "parameter") == 0 && strcmp(str_tmp2, "w(j)/") == 0)
					Instance.we[item - 1] = atof(str_reading);
				else if (strcmp(str_tmp1, "parameter") == 0 && strcmp(str_tmp2, "po(j)/") == 0)
					Instance.po[item - 1] = atof(str_reading);		
				else if (strcmp(str_tmp1, "parameter") == 0 && strcmp(str_tmp2, "s(r)/") == 0)
					Instance.sr[item - 1] = atof(str_reading);
				else if (strcmp(str_tmp1, "parameter") == 0 && strcmp(str_tmp2, "nr(r)/") == 0)
					Instance.nr[item - 1] = atoi(str_reading);
			}			
			for(int i=0;i<Instance.items;i++)
			cout<<"i="<<i+1<<" "<< Instance.we[i]<<endl;						
		}		
		else if (strcmp(str_last, "cap(k);") == 0 && strcmp(str_reading, "cap(k)=") == 0)				//for capacity
		{
			FIC >> str_reading;
			Instance.capac = atof(str_reading);		
		}		
		else if ((strcmp(str_last, "parameter") == 0 && strcmp(str_reading, "pp(i,j)/") == 0) ||		//for q_ij, t_rj and psi(r, k)
			(strcmp(str_last, "parameter") == 0 && strcmp(str_reading, "t(r,j)/") == 0) || 
			(strcmp(str_last, "parameter") == 0 && strcmp(str_reading, "psi(r,k)/") == 0))
		{
			strcpy(str_tmp1, str_last);
			strcpy(str_tmp2, str_reading);
			while (1)
			{
				FIC >> str_last >> str_reading;
				if (strcmp(str_last, "/;") == 0)
					break;
				int item1 = 0, item2 = 0;
				int num_len = 0;
				int pos_point = -1;
				for (int i = strlen(str_last) - 1; i >= 0; i--)
				{
					if (str_last[i] >= '0' && str_last[i] <= '9')
						item2 += (str_last[i] - '0') * int(pow(10, num_len++));
					else if (str_last[i] == '.')
					{
						pos_point = i;
						break;
					}
				}
				num_len = 0;
				for (int i = pos_point - 1; i >= 0; i--)
					if (str_last[i] >= '0' && str_last[i] <= '9')
						item1 += (str_last[i] - '0') * int(pow(10, num_len++));
				if (strcmp(str_tmp1, "parameter") == 0 && strcmp(str_tmp2, "pp(i,j)/") == 0)
					Instance.q_ij[item1 - 1][item2 - 1] = Instance.q_ij[item2 - 1][item1 - 1] = atof(str_reading);
				else if (strcmp(str_tmp1, "parameter") == 0 && strcmp(str_tmp2, "t(r,j)/") == 0)
				{
					Instance.t_rj[item1 - 1][item2 - 1] = atoi(str_reading);
					Instance.t_jr[item2 - 1] = item1 - 1;
				}
				else if (strcmp(str_tmp1, "parameter") == 0 && strcmp(str_tmp2, "psi(r,k)/") == 0)
					Instance.psi[item1 - 1][item2 - 1] = atof(str_reading);							
			}
		}	
		if (strcmp(str_reading, "epsilon(j,k)=") == 0)
			break;
	}
	//calculate p(j,k) according to p(j,k)= po(j)*sum(r,t(r,j)*psi(r,k));	
	for (int j = 0; j < Instance.items; j++)
	{
		for (int k = 0; k < Instance.knaps; k++)
		{
			double sum = 0;
			for (int r = 0; r < Instance.classes; r++)			
				sum += Instance.t_rj[r][j] * Instance.psi[r][k];
			Instance.p_jk[j][k] = Instance.po[j] * sum;
		}		
	}
	FIC.close();
#ifdef DEBUG_READ_INSTANCE
	cout << "items=" << Instance.items << ", knapsacks=" << Instance.knaps << ", class=" << Instance.classes << ",u=" << Instance.u << endl;
	cout << "capacity=" << Instance.capac << endl;
	cout << "we===================" << endl;
	for (int j = 0; j < Instance.items; j++)
		cout << "j=" << j + 1 << ", we=" << Instance.we[j] << endl;
	cout << "po===================" << endl;
	for (int j = 0; j < Instance.items; j++)
		cout << "j=" << j + 1 << ", po=" << Instance.po[j] << endl;
	cout << "qij===================" << endl;
	for (int i = 0; i < Instance.items; i++)
		for (int j = 0; j < Instance.items; j++)
			cout << i + 1 << "." << j + 1 << "= " << Instance.q_ij[i][j] << endl;
	cout << "trj===================" << endl;
	for (int i = 0; i < Instance.classes; i++)
		for (int j = 0; j < Instance.items;j++)
			cout << i + 1 << "." << j + 1 << "= " << Instance.t_rj[i][j] << endl;
	cout << "sr===================" << endl;
	for (int i = 0; i < Instance.classes; i++)
		cout << i + 1 << " " << Instance.sr[i] << endl;
	cout << "nr===================" << endl;
	for (int i = 0; i < Instance.classes; i++)
		cout << i + 1 << " " << Instance.nr[i] << endl;
	cout << "psi_rk===================" << endl;
	for (int i = 0; i < Instance.classes; i++)
		for (int j = 0; j < Instance.knaps; j++)
			cout << i + 1 << "." << j + 1 << "= " << Instance.psi[i][j] << endl;
	cout << "p_jk===================" << endl;
	for (int j = 0; j < Instance.items; j++)
		for (int k = 0; k < Instance.knaps; k++)
			cout << j + 1 << "." << k + 1 << "= " << Instance.p_jk[j][k] << endl;
#endif
}

//Complexity: 0(1) time
bool satisfy_all_constraints_relocate(int item_j, int old_knap, int new_knap)
{
	if (satisfy_cap_constraints_relocate(item_j, old_knap, new_knap) && satisfy_nr_constraints_relocate(item_j, old_knap, new_knap))
		return true;
	else
		return false;
}

//Complexity: 0(1) time
bool satisfy_nr_constraints_relocate(int item_j, int old_knap, int new_knap)
{
	int class_j_r = Instance.t_jr[item_j];
	int size_r = Size_r[class_j_r];
	if (Size_kr[old_knap][class_j_r] == 1 && old_knap < Instance.knaps)
		size_r--;
	if (Size_kr[new_knap][class_j_r] == 0 && new_knap < Instance.knaps)
		size_r++;
	if (size_r > Instance.nr[class_j_r])
		return false;	
	return true;
}

//Complexity: 0(1) time
bool satisfy_cap_constraints_relocate(int item_j, int old_knap, int new_knap)
{
	int class_j_r = Instance.t_jr[item_j];
	double sum1 = Capacity1_k[new_knap] + Instance.we[item_j];
	double sum2 = Capacity2_k[new_knap];
	if (Size_kr[new_knap][class_j_r] == 0)
		sum2 += Instance.sr[class_j_r];
	if (sum1 + sum2 > Instance.capac + PRECISION && new_knap < Instance.knaps)
		return false;
	return true;
}

//Complexity: 0(1) time
bool satisfy_all_constraints_swap(int item_j1, int knap_j1, int item_j2, int knap_j2)
{
	if (satisfy_nr_constraints_swap(item_j1, knap_j1, item_j2, knap_j2) && satisfy_cap_constraints_swap(item_j1, knap_j1, item_j2, knap_j2))
		return true;
	else
		return false;
}

//Complexity: 0(1) time
bool satisfy_nr_constraints_swap(int item_j1, int knap_j1, int item_j2, int knap_j2)
{
	int class_j1_r = Instance.t_jr[item_j1];
	int class_j2_r = Instance.t_jr[item_j2];
	int size_r1 = Size_r[class_j1_r];
	int size_r2 = Size_r[class_j2_r];
	if (Size_kr[knap_j1][class_j1_r] == 1 && class_j2_r != class_j1_r && knap_j1 < Instance.knaps)
		size_r1--;
	if (Size_kr[knap_j2][class_j1_r] == 0 && knap_j2 < Instance.knaps)
		size_r1++;

	if (Size_kr[knap_j2][class_j2_r] == 1 && class_j1_r != class_j2_r  && knap_j2 < Instance.knaps)
		size_r2--;
	if (Size_kr[knap_j1][class_j2_r] == 0 && knap_j1 < Instance.knaps)
		size_r2++;
	if (size_r1 > Instance.nr[class_j1_r] || size_r2 > Instance.nr[class_j2_r])
		return false;
	return true;
}

//Complexity: 0(1) time
bool satisfy_cap_constraints_swap(int item_j1, int knap_j1, int item_j2, int knap_j2)
{
	int class_j1_r = Instance.t_jr[item_j1];
	int class_j2_r = Instance.t_jr[item_j2];
	double sum1_k1 = Capacity1_k[knap_j1] + Instance.we[item_j2] - Instance.we[item_j1];
	double sum1_k2 = Capacity1_k[knap_j2] + Instance.we[item_j1] - Instance.we[item_j2];
	double sum2_k1 = Capacity2_k[knap_j1];
	double sum2_k2 = Capacity2_k[knap_j2];

	if (Size_kr[knap_j1][class_j1_r] == 1 && class_j2_r != class_j1_r)
		sum2_k1 -= Instance.sr[class_j1_r];
	if (Size_kr[knap_j1][class_j2_r] == 0)
		sum2_k1 += Instance.sr[class_j2_r];

	if (Size_kr[knap_j2][class_j2_r] == 1 && class_j1_r != class_j2_r)
		sum2_k2 -= Instance.sr[class_j2_r];
	if (Size_kr[knap_j2][class_j1_r] == 0)
		sum2_k2 += Instance.sr[class_j1_r];

	if ((sum1_k1 + sum2_k1 > Instance.capac + PRECISION && knap_j1 < Instance.knaps) || (sum1_k2 + sum2_k2 > Instance.capac + PRECISION && knap_j2 < Instance.knaps))
		return false;
	return true;
}

//Complexity: 0(1) time, item x and y exchange with item z, and x, y belong the same knapsack
bool satisfy_all_constraints_2_1_exchange(int item_x, int item_y, int knap_x_y, int item_z, int knap_z)
{
	int class_x_r = Instance.t_jr[item_x];
	int class_y_r = Instance.t_jr[item_y];
	int class_z_r = Instance.t_jr[item_z];

	int size_x_r = Size_r[class_x_r];
	int size_y_r = Size_r[class_y_r];
	int size_z_r = Size_r[class_z_r];

	if (class_x_r != class_y_r)
	{
		if (Size_kr[knap_x_y][class_x_r] == 1 && class_z_r != class_x_r && knap_x_y < Instance.knaps)
			size_x_r--;
		if (Size_kr[knap_z][class_x_r] == 0 && knap_z < Instance.knaps)
			size_x_r++;

		if (Size_kr[knap_x_y][class_y_r] == 1 && class_z_r != class_y_r && knap_x_y < Instance.knaps)
			size_y_r--;
		if (Size_kr[knap_z][class_y_r] == 0 && knap_z < Instance.knaps)
			size_y_r++;
	}
	else
	{
		if (Size_kr[knap_x_y][class_x_r] == 2 && class_z_r != class_x_r && knap_x_y < Instance.knaps)
			size_x_r--;
		if (Size_kr[knap_z][class_x_r] == 0 && knap_z < Instance.knaps)
			size_x_r++;
		size_y_r = size_x_r;
	}
	if (Size_kr[knap_z][class_z_r] == 1 && class_x_r != class_z_r && class_y_r != class_z_r  && knap_z < Instance.knaps)
		size_z_r--;
	if (Size_kr[knap_x_y][class_z_r] == 0 && knap_x_y < Instance.knaps)
		size_z_r++;
	if (size_x_r > Instance.nr[class_x_r] || size_y_r > Instance.nr[class_y_r] || size_z_r > Instance.nr[class_z_r])
		return false;

	double sum1_k1 = Capacity1_k[knap_x_y] + Instance.we[item_z] - Instance.we[item_x] - Instance.we[item_y];
	double sum1_k2 = Capacity1_k[knap_z] + Instance.we[item_x] + Instance.we[item_y] - Instance.we[item_z];
	double sum2_k1 = Capacity2_k[knap_x_y];
	double sum2_k2 = Capacity2_k[knap_z];


	if (class_x_r != class_y_r)
	{		
		if (Size_kr[knap_x_y][class_x_r] == 1 && class_z_r != class_x_r)
			sum2_k1 -= Instance.sr[class_x_r];
		if (Size_kr[knap_x_y][class_y_r] == 1 && class_z_r != class_y_r)
			sum2_k1 -= Instance.sr[class_y_r];
		if (Size_kr[knap_z][class_x_r] == 0)
			sum2_k2 += Instance.sr[class_x_r];
		if (Size_kr[knap_z][class_y_r] == 0)
			sum2_k2 += Instance.sr[class_y_r];
	}
	else
	{
		if (Size_kr[knap_x_y][class_x_r] == 2 && class_z_r != class_x_r)
			sum2_k1 -= Instance.sr[class_x_r];	
		if (Size_kr[knap_z][class_x_r] == 0)
			sum2_k2 += Instance.sr[class_x_r];	
	}
	if (Size_kr[knap_z][class_z_r] == 1 && class_x_r != class_z_r && class_y_r != class_z_r)
		sum2_k2 -= Instance.sr[class_z_r];
	if (Size_kr[knap_x_y][class_z_r] == 0)
		sum2_k1 += Instance.sr[class_z_r];

	if ((sum1_k1 + sum2_k1 > Instance.capac + PRECISION && knap_x_y < Instance.knaps) || (sum1_k2 + sum2_k2 > Instance.capac + PRECISION && knap_z < Instance.knaps))
		return false;
	return true;
}

void initial_global_data()
{
	for (int k = 0; k < Instance.knaps + 1; k++)
		for (int r = 0; r < Instance.classes; r++)
			Size_kr[k][r] = 0;
	for (int k = 0; k < Instance.knaps + 1; k++)
	{
		Capacity1_k[k] = 0.0;
		Capacity2_k[k] = 0.0;
		Exceed_cap_k[k] = 0.0;
	}
	for (int r = 0; r < Instance.classes; r++)
		Exceed_class_r[r] = 0;
	for (int r = 0; r < Instance.classes; r++)	
		Size_r[r] = 0;	
}

void build_global_data(int *sol)
{
	initial_global_data();
	for (int j = 0; j < Instance.items; j++)
		Size_kr[sol[j]][Instance.t_jr[j]]++;	
	for (int j = 0; j < Instance.items; j++)
		Capacity1_k[sol[j]] += Instance.we[j];
	for (int k = 0; k < Instance.knaps + 1; k++)
		for (int r = 0; r < Instance.classes; r++)
			if (Size_kr[k][r]>0)
				Capacity2_k[k] += Instance.sr[r];
	for (int r = 0; r < Instance.classes; r++)
		Size_r[r] = 0;
	for (int r = 0; r < Instance.classes; r++)	
		for (int k = 0; k < Instance.knaps; k++)	
			if (Size_kr[k][r] > 0)
				Size_r[r]++;		
}

void update_global_data(int item, int old_knap, int new_knap)
{
	int item_r = Instance.t_jr[item];
	Size_kr[old_knap][item_r]--;
	Size_kr[new_knap][item_r]++;
	Capacity1_k[old_knap] -= Instance.we[item];
	Capacity1_k[new_knap] += Instance.we[item];	
	if (Size_kr[old_knap][item_r] == 0)
	{
		if (old_knap < Instance.knaps)
			Size_r[item_r]--;
		Capacity2_k[old_knap] -= Instance.sr[item_r];
	}
	if (Size_kr[new_knap][item_r] == 1)
	{		
		if (new_knap < Instance.knaps)
			Size_r[item_r]++;
		Capacity2_k[new_knap] += Instance.sr[item_r];
	}	
}

void update_global_data_exceed(int old_knap, int new_knap, int r1, int r2)
{	
	if (Capacity1_k[old_knap] + Capacity2_k[old_knap] > Instance.capac + PRECISION)
		Exceed_cap_k[old_knap] = Capacity1_k[old_knap] + Capacity2_k[old_knap] - Instance.capac;
	else
		Exceed_cap_k[old_knap] = 0.0;
	if (Capacity1_k[new_knap] + Capacity2_k[new_knap] > Instance.capac + PRECISION)
		Exceed_cap_k[new_knap] = Capacity1_k[new_knap] + Capacity2_k[new_knap] - Instance.capac;
	else
		Exceed_cap_k[new_knap] = 0.0;

	if (Size_r[r1] > Instance.nr[r1])
		Exceed_class_r[r1] = Size_r[r1] - Instance.nr[r1];
	else
		Exceed_class_r[r1] = 0;
	if (Size_r[r2] > Instance.nr[r2])
		Exceed_class_r[r2] = Size_r[r2] - Instance.nr[r2];
	else
		Exceed_class_r[r2] = 0;
}


void initial_sol_helpfunc(int *sol)
{
	initial_global_data();
	for (int j = 0; j < Instance.items; j++)
		sol[j] = Instance.knaps;					//indicate item j is unallocated
	for (int j = 0; j < Instance.items; j++)
	{
		int j_r = Instance.t_jr[j];
		Size_kr[Instance.knaps][j_r]++;
	}
	for (int j = 0; j < Instance.items; j++)
		if (sol[j] == Instance.knaps)
			Capacity1_k[Instance.knaps] += Instance.we[j];
	for (int r = 0; r < Instance.classes; r++)
		if (Size_kr[Instance.knaps][r] > 0)
			Capacity2_k[Instance.knaps] += Instance.sr[r];
}


//unused
void random_initial_sol(int *sol)
{
	int *can_item = new int[Instance.items*Instance.knaps];
	int *can_knaps = new int[Instance.items*Instance.knaps];
	initial_sol_helpfunc(sol);
	
	while (1)
	{
		int can_len = 0;
		for (int j = 0; j < Instance.items; j++)
		{
			if (sol[j] == Instance.knaps)
			{
				int j_r = Instance.t_jr[j];				
				for (int k = 0; k < Instance.knaps; k++)
				{					
					if (satisfy_all_constraints_relocate(j, sol[j], k))
					{
						can_item[can_len] = j;
						can_knaps[can_len] = k;
						can_len++;
					}
				}
			}
			//cout << "can_len=" << can_len << endl;
			//getchar();
		}
		cout << "can_len=" << can_len << endl;
		if (can_len == 0)
			break;
		else
		{
			int rx = rand() % can_len;
			int item = can_item[rx];
			int knap = can_knaps[rx];
			int old_knap = sol[item];
			update_global_data(item, old_knap, knap);
			sol[item] = knap;
		}
	}
	delete[]can_item; can_item = NULL;
	delete[]can_knaps; can_knaps = NULL;
#ifdef DEBUG_INITIAL_SOL
	cout << "items=" << Instance.items << ", k=" << Instance.knaps << ", r=" << Instance.classes << ", Capacity=" << Instance.capac << endl;
	cout << "nr============" << endl;
	for (int r = 0; r < Instance.classes; r++)
		cout << Instance.nr[r] << " ";
	cout << endl;
	cout << "Size_kr=====================" << endl;
	for (int k = 0; k < Instance.knaps + 1; k++)
	{
		cout << "k=" << k << "===" << endl;
		for (int r = 0; r < Instance.classes; r++)		
			cout << Size_kr[k][r] << " ";
		cout << endl;
	}
	cout << "Capacity1_k============" << endl;
	for (int k = 0; k < Instance.knaps + 1; k++)
		cout << "k=" << k << " " << Capacity1_k[k] << endl;
	cout << "Capacity2_k============" << endl;
	for (int k = 0; k < Instance.knaps + 1; k++)
		cout << "k=" << k << " " << Capacity2_k[k] << endl;
	cout << "Size_r================" << endl;
	for (int r = 0; r < Instance.classes; r++)
		cout << "r=" << r << " " << Size_r[r] << endl;
#endif
}

//based on value density
void greedy_helpfunc_value_density(int *sol, int *can_item, int *can_knaps, int *address_knaps, int &len_knap)
{
	int best_vd_item[MAXCAN];
	int rx = rand() % len_knap;
	int sel_knap = can_knaps[rx];			//randomly select a knapsack from the list can_knaps (KL)
	int pos = address_knaps[sel_knap];
	while (1)
	{
		int len_item = 0;
		for (int j = 0; j < Instance.items; j++)
		{
			if (sol[j] == Instance.knaps && satisfy_all_constraints_relocate(j, sol[j], sel_knap))
				can_item[len_item++] = j;
		}
		if (len_item == 0)					//the list can_item (FL) is empty
		{
			for (int k = pos; k < len_knap - 1; k++)
			{
				can_knaps[k] = can_knaps[k + 1];
				address_knaps[can_knaps[k]] = k;
			}
			len_knap--;
			break;
		}
		else
		{
			int vd_len = 0;
			double vd_maximum = -MAXNUM;
			for (int jnd = 0; jnd < len_item; jnd++)
			{
				int item_j = can_item[jnd];
				double vd_jk = Instance.p_jk[item_j][sel_knap];
				for (int j2 = 0; j2 < Instance.items; j2++)
					if (sol[j2] == sel_knap)
						vd_jk += Instance.q_ij[item_j][j2];
				//if (fabs(Instance.we[item_j] - 0) > PRECISION)
				//	vd_jk /= Instance.we[item_j];
				//else
				//	vd_jk = MAXNUM;

				if (fabs(Instance.we[item_j] - 0) <= PRECISION)
					vd_jk = MAXNUM;
				else
				{
					int class_j_r = Instance.t_jr[item_j];
					if (Size_kr[sel_knap][class_j_r] == 0)
						vd_jk /= (Instance.we[item_j] + Instance.sr[class_j_r]);
					else
						vd_jk /= Instance.we[item_j];
				}

				if (vd_jk > vd_maximum + PRECISION)
				{
					vd_len = 0;
					best_vd_item[vd_len] = item_j;
					vd_maximum = vd_jk;
					vd_len++;
				}
				else if (fabs(vd_jk - vd_maximum) <= PRECISION && vd_len < MAXCAN)
				{
					best_vd_item[vd_len] = item_j;
					vd_len++;
				}
			}
			int rx_item = rand() % vd_len;
			int sel_item = best_vd_item[rx_item];
			update_global_data(sel_item, sol[sel_item], sel_knap);
			sol[sel_item] = sel_knap;
			//cout << "sel_item=" << sel_item << ", sel_knap=" << sel_knap << ", len_knap=" << len_knap << endl;
		}
	}
}

//according to the paper "A multi-start iterated local search algorithm for the generalized quadratic multiple knapsack problem"
void greedy_initial_sol(int *sol)
{	
	int *can_item = new int[Instance.items];
	int *can_knaps = new int[Instance.knaps];
	int *address_knaps = new int[Instance.knaps];	
	initial_sol_helpfunc(sol);
	int len_knap = 0;
	for (int k = 0; k < Instance.knaps; k++)
	{
		address_knaps[k] = len_knap;
		can_knaps[len_knap] = k;
		len_knap++;
	}
	while (len_knap > 0)
	{		
		greedy_helpfunc_value_density(sol, can_item, can_knaps, address_knaps, len_knap);
	}
	delete[]can_item; can_item = NULL;
	delete[]can_knaps; can_knaps = NULL;
	delete[]address_knaps;address_knaps=NULL;
#ifdef DEBUG_INITIAL_SOL
	cout << "items=" << Instance.items << ", k=" << Instance.knaps << ", r=" << Instance.classes << ", Capacity=" << Instance.capac << endl;
	cout << "nr============" << endl;
	for (int r = 0; r < Instance.classes; r++)
		cout << Instance.nr[r] << " ";
	cout << endl;
	cout << "Size_kr=====================" << endl;
	for (int k = 0; k < Instance.knaps + 1; k++)
	{
		cout << "k=" << k << "===" << endl;
		for (int r = 0; r < Instance.classes; r++)
			cout << Size_kr[k][r] << " ";
		cout << endl;
	}
	cout << "Capacity1_k============" << endl;
	for (int k = 0; k < Instance.knaps + 1; k++)
		cout << "k=" << k << " " << Capacity1_k[k] << endl;
	cout << "Capacity2_k============" << endl;
	for (int k = 0; k < Instance.knaps + 1; k++)
		cout << "k=" << k << " " << Capacity2_k[k] << endl;
	cout << "Size_r================" << endl;
	for (int r = 0; r < Instance.classes; r++)
		cout << "r=" << r << " " << Size_r[r] << endl;
#endif
}

//call the greedy_intial_solution func for  Num_greedy_ini times, and choose the solution with the best objective value
void initial_sol(int *sol)
{
	int *sol_tmp = new int[Instance.items];
	double obj_maximum = -MAXNUM;
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
		if (f_sol > obj_maximum + PRECISION)
		{
			memcpy(sol_tmp, sol, sizeof(int)*Instance.items);
			obj_maximum = f_sol;
		}
		//cout << "i=" << i << ", f_sol=" << f_sol << endl;
	}
	memcpy(sol, sol_tmp, sizeof(int)*Instance.items);	
	double initial_time = 1.0*(clock() - time1) / CLOCKS_PER_SEC;
	cout << "in initial_sol func, initial_time=" << initial_time << endl;
	delete[]sol_tmp; sol_tmp = NULL;
}

void build_delta_matrix(int *sol)
{
	for (int j = 0; j < Instance.items; j++)
		for (int k = 0; k < Instance.knaps + 1; k++)
			Delta_matrix_jk[j][k] = 0.0;
	for (int j = 0; j < Instance.items; j++)	
		for (int j2 = 0; j2 < Instance.items; j2++)		
			Delta_matrix_jk[j][sol[j2]] += Instance.q_ij[j][j2];	
}

void update_delta_matrix(int item, int old_knap, int new_knap)
{
	for (int j = 0; j < Instance.items; j++)
	{
		if( j != item)
		{
			Delta_matrix_jk[j][old_knap] -= Instance.q_ij[item][j];
			Delta_matrix_jk[j][new_knap] += Instance.q_ij[item][j];
		}
	}
}

double compute_obj_part1(int *sol)
{
	double f_part1 = 0;
	for (int j = 0; j < Instance.items; j++)
		if (sol[j] < Instance.knaps)
			f_part1 += Instance.p_jk[j][sol[j]];
	return f_part1;
}

double compute_obj_part2(int *sol)
{
	double f_part2 = 0;
	for (int j = 0; j < Instance.items; j++)
		if (sol[j] < Instance.knaps)
			f_part2 += Delta_matrix_jk[j][sol[j]];
	f_part2 /= 2;
	return f_part2;
}

double compute_delta_relocate(int *sol, int item_j, int knap)
{
	double delta_part1, delta_part2;
	if (sol[item_j] == Instance.knaps)
	{
		delta_part1 = Instance.p_jk[item_j][knap];
		delta_part2 = Delta_matrix_jk[item_j][knap];
	}
	else if (knap == Instance.knaps)
	{
		delta_part1 = -Instance.p_jk[item_j][sol[item_j]];
		delta_part2 = -Delta_matrix_jk[item_j][sol[item_j]];
	}
	else
	{
		delta_part1 = Instance.p_jk[item_j][knap] - Instance.p_jk[item_j][sol[item_j]];
		delta_part2 = Delta_matrix_jk[item_j][knap] - Delta_matrix_jk[item_j][sol[item_j]];
	}
	return (delta_part1 + delta_part2);
}

double compute_delta_swap(int *sol, int item_j, int item_j2)
{
	double delta_part1, delta_part2;
	if (sol[item_j] == Instance.knaps)
	{
		delta_part1 = Instance.p_jk[item_j][sol[item_j2]] - Instance.p_jk[item_j2][sol[item_j2]];
		delta_part2 = Delta_matrix_jk[item_j][sol[item_j2]] - Delta_matrix_jk[item_j2][sol[item_j2]] - Instance.q_ij[item_j][item_j2];
	}
	else if (sol[item_j2] == Instance.knaps)
	{
		delta_part1 = -Instance.p_jk[item_j][sol[item_j]] + Instance.p_jk[item_j2][sol[item_j]];
		delta_part2 = -Delta_matrix_jk[item_j][sol[item_j]] + Delta_matrix_jk[item_j2][sol[item_j]] - Instance.q_ij[item_j][item_j2];
	}
	else
	{
		delta_part1 = Instance.p_jk[item_j][sol[item_j2]] - Instance.p_jk[item_j][sol[item_j]] + Instance.p_jk[item_j2][sol[item_j]] - Instance.p_jk[item_j2][sol[item_j2]];
		delta_part2 = Delta_matrix_jk[item_j][sol[item_j2]] - Delta_matrix_jk[item_j][sol[item_j]] +
			Delta_matrix_jk[item_j2][sol[item_j]] - Delta_matrix_jk[item_j2][sol[item_j2]] - 2 * Instance.q_ij[item_j][item_j2];
	}
	return (delta_part1 + delta_part2);
}

double compute_delta_2_1_exchange(int *sol, int item_x, int item_y, int item_z)
{
	double delta_part1, delta_part2;
	int knap_x_y = sol[item_x];
	int knap_z = sol[item_z];
	if (knap_x_y == Instance.knaps)
	{
		delta_part1 = Instance.p_jk[item_x][knap_z] + Instance.p_jk[item_y][knap_z] - Instance.p_jk[item_z][knap_z];
		delta_part2 = Delta_matrix_jk[item_x][knap_z] + Delta_matrix_jk[item_y][knap_z] -
			Delta_matrix_jk[item_z][knap_z] - Instance.q_ij[item_x][item_z] - Instance.q_ij[item_y][item_z] + Instance.q_ij[item_x][item_y];
	}
	else if (knap_z == Instance.knaps)
	{
		delta_part1 = -Instance.p_jk[item_x][knap_x_y] - Instance.p_jk[item_y][knap_x_y] + Instance.p_jk[item_z][knap_x_y];
		delta_part2 = -Delta_matrix_jk[item_x][knap_x_y] - Delta_matrix_jk[item_y][knap_x_y] +
			Delta_matrix_jk[item_z][knap_x_y] - Instance.q_ij[item_x][item_z] - Instance.q_ij[item_y][item_z] +  Instance.q_ij[item_x][item_y];
	}
	else
	{
		delta_part1 = Instance.p_jk[item_x][knap_z] - Instance.p_jk[item_x][knap_x_y] + Instance.p_jk[item_y][knap_z] -	Instance.p_jk[item_y][knap_x_y] + 
			Instance.p_jk[item_z][knap_x_y] - Instance.p_jk[item_z][knap_z];
		delta_part2 = Delta_matrix_jk[item_x][knap_z] - Delta_matrix_jk[item_x][knap_x_y] + Delta_matrix_jk[item_y][knap_z] - Delta_matrix_jk[item_y][knap_x_y] +
			Delta_matrix_jk[item_z][knap_x_y] - Delta_matrix_jk[item_z][knap_z] - 2 * Instance.q_ij[item_x][item_z] - 2 * Instance.q_ij[item_y][item_z] + 2 * Instance.q_ij[item_x][item_y];
	}
	return (delta_part1 + delta_part2);
}

//g1 (s) sums the exceeded weights of the knapsacks 
//g2 (s) sums the exceeded numbers of the knapsacks for each class r, unused
double compute_g1()
{
	double gs = 0.0;	
	for (int k = 0; k < Instance.knaps; k++)
		gs += Exceed_cap_k[k];
	return gs;
}

double compute_delta_relocate_g1(int *sol, int item_j, int new_knap)
{
	int class_j_r = Instance.t_jr[item_j];
	int old_knap = sol[item_j];
	double sum1 = 0.0, sum2 = 0.0;
	double delta_g1 = 0.0, delta_g1_old_k = 0.0, delta_g1_new_k = 0.0;
	if (old_knap != Instance.knaps)
	{
		sum1 = Capacity1_k[old_knap] - Instance.we[item_j];
		sum2 = Capacity2_k[old_knap];
		if (Size_kr[old_knap][class_j_r] == 1)		
			sum2 -= Instance.sr[class_j_r];		
		double diff = sum1 + sum2 - Instance.capac;
		if (diff > 0 + PRECISION)
			delta_g1_old_k = diff - Exceed_cap_k[old_knap];
		else
			delta_g1_old_k = -Exceed_cap_k[old_knap];
	}
	if (new_knap != Instance.knaps)
	{
		sum1 = Capacity1_k[new_knap] + Instance.we[item_j];
		sum2 = Capacity2_k[new_knap];
		if (Size_kr[new_knap][class_j_r] == 0)		
			sum2 += Instance.sr[class_j_r];	
		
		double diff = sum1 + sum2 - Instance.capac;
		if (diff > 0 + PRECISION)
			delta_g1_new_k = diff - Exceed_cap_k[new_knap];
		else
			delta_g1_new_k = -Exceed_cap_k[new_knap];
	}
	delta_g1 = delta_g1_old_k + delta_g1_new_k;	
	return delta_g1;
}

int compute_delta_relocate_g2(int *sol, int item_j, int new_knap)
{
	int class_j_r = Instance.t_jr[item_j];
	int old_knap = sol[item_j];
	int size_r = Size_r[class_j_r];
	int delta_g2_r = 0;
	if (old_knap != Instance.knaps && Size_kr[old_knap][class_j_r] == 1)
		size_r--;
	if (Size_kr[new_knap][class_j_r] == 0 && new_knap != Instance.knaps)
		size_r++;
	int diff = size_r - Instance.nr[class_j_r];
	if (diff > 0)
		delta_g2_r = diff - Exceed_class_r[class_j_r];
	else
		delta_g2_r = -Exceed_class_r[class_j_r];
	return delta_g2_r;
}

double compute_delta_swap_g1(int *sol, int item_j, int item_j2)
{
	int class_j_r = Instance.t_jr[item_j];
	int class_j2_r = Instance.t_jr[item_j2];
	int j_knap = sol[item_j];
	int j2_knap = sol[item_j2];
	double sum1 = 0.0, sum2 = 0.0;
	double delta_g1 = 0.0, delta_g1_j_knap = 0.0, delta_g1_j2_knap = 0.0;
	if (j_knap != Instance.knaps)
	{
		sum1 = Capacity1_k[j_knap] + Instance.we[item_j2] - Instance.we[item_j];
		sum2 = Capacity2_k[j_knap];
		if (Size_kr[j_knap][class_j_r] == 1 && class_j2_r != class_j_r)
			sum2 -= Instance.sr[class_j_r];
		if (Size_kr[j_knap][class_j2_r] == 0)
			sum2 += Instance.sr[class_j2_r];
		double diff = sum1 + sum2 - Instance.capac;
		if (diff > 0 + PRECISION)
			delta_g1_j_knap = diff - Exceed_cap_k[j_knap];
		else
			delta_g1_j_knap = -Exceed_cap_k[j_knap];
	}
	if (j2_knap != Instance.knaps)
	{
		sum1 = Capacity1_k[j2_knap] - Instance.we[item_j2] + Instance.we[item_j];
		sum2 = Capacity2_k[j2_knap];
		if (Size_kr[j2_knap][class_j2_r] == 1 && class_j_r != class_j2_r)
			sum2 -= Instance.sr[class_j2_r];
		if (Size_kr[j2_knap][class_j_r] == 0)
			sum2 += Instance.sr[class_j_r];

		double diff = sum1 + sum2 - Instance.capac;
		if (diff > 0 + PRECISION)
			delta_g1_j2_knap = diff - Exceed_cap_k[j2_knap];
		else
			delta_g1_j2_knap = -Exceed_cap_k[j2_knap];
	}
	delta_g1 = delta_g1_j_knap + delta_g1_j2_knap;
	return delta_g1;
}

int compute_delta_swap_g2(int *sol, int item_j, int item_j2)
{
	int class_j_r = Instance.t_jr[item_j];
	int class_j2_r = Instance.t_jr[item_j2];
	int j_knap = sol[item_j];
	int j2_knap = sol[item_j2];
	int size_r = Size_r[class_j_r];
	int size_r2 = Size_r[class_j2_r];
	int delta_g2_r = 0;
	int delta_g2_r2 = 0;
	if (j_knap != Instance.knaps && Size_kr[j_knap][class_j_r] == 1 && class_j2_r != class_j_r)
		size_r--;
	if (Size_kr[j2_knap][class_j_r] == 0 && j2_knap != Instance.knaps)
		size_r++;
	if (j2_knap != Instance.knaps && Size_kr[j2_knap][class_j2_r] == 1 && class_j2_r != class_j_r)
		size_r2--;
	if (Size_kr[j_knap][class_j2_r] == 0 && j_knap != Instance.knaps)
		size_r2++;
	int diff = size_r - Instance.nr[class_j_r];
	int diff2 = size_r2 - Instance.nr[class_j2_r];
	if (diff > 0)
		delta_g2_r = diff - Exceed_class_r[class_j_r];
	else
		delta_g2_r = -Exceed_class_r[class_j_r];
	if (diff2 > 0)
		delta_g2_r2 = diff2 - Exceed_class_r[class_j2_r];
	else
		delta_g2_r2 = -Exceed_class_r[class_j2_r];
	return delta_g2_r + delta_g2_r2;
}

//relocate move
void relocate_move(int *sol, int item, int knap)
{
	int old_knap = sol[item];
	int new_knap = knap;
	int r = Instance.t_jr[item];
	update_global_data(item, old_knap, new_knap);
	update_global_data_exceed(old_knap, new_knap, r, r);
	update_delta_matrix(item, old_knap, new_knap);
	sol[item] = knap;
}

//swap move
void swap_move(int *sol, int item1, int item2)
{
	int knap1 = sol[item1];
	int knap2 = sol[item2];
	int r1 = Instance.t_jr[item1];
	int r2 = Instance.t_jr[item2];
	update_global_data(item1, knap1, knap2);
	update_global_data(item2, knap2, knap1);
	update_global_data_exceed(knap1, knap2, r1, r2);
	update_delta_matrix(item1, knap1, knap2);
	update_delta_matrix(item2, knap2, knap1);
	sol[item1] = knap2;
	sol[item2] = knap1;
}

//2-1 exchange move
void two_one_exchange_move(int *sol, int item_x, int item_y, int item_z)
{
	int knap1 = sol[item_x];
	int knap2 = sol[item_z];
	update_global_data(item_x, knap1, knap2);
	update_global_data(item_y, knap1, knap2);
	update_global_data(item_z, knap2, knap1);
	update_global_data_exceed(knap1, knap2, -1, -1);
	update_delta_matrix(item_x, knap1, knap2);
	update_delta_matrix(item_y, knap1, knap2);
	update_delta_matrix(item_z, knap2, knap1);
	sol[item_x] = knap2;
	sol[item_y] = knap2;
	sol[item_z] = knap1;
}

//for the perturbation
void build_neighbors()
{
	int count = 0;
	//int numbers = N*(N - 1) / 2 + N*K;
	//relocate neighborhood
	for (int i = 0; i < Instance.items; i++)
	{
		for (int g = 0; g < Instance.knaps + 1; g++)
		{
			Neighbors[count].type = 1;
			Neighbors[count].v = i;
			Neighbors[count].g = g;
			count++;
		}
	}
	//swap neighborhood
	for (int i = 0; i < Instance.items; i++)
	{
		for (int j = i + 1; j < Instance.items; j++)
		{
			Neighbors[count].type = 2;
			Neighbors[count].x = i;
			Neighbors[count].y = j;
			count++;
		}
	}
}

void update_best_sol_one_run(double f_sol, int *sol)
{
	Best_sol_one_run_obj = f_sol;
	for (int j = 0; j < Instance.items; j++)
		Best_sol_one_run[j] = sol[j];
	Run_time = 1.0*(clock() - Start_time) / CLOCKS_PER_SEC;
}

bool is_same_solution(int *sol_1, int *sol_2)
{
	bool is_exist = true;
	for (int j = 0; j < Instance.items; j++)
	{
		if (sol_1[j] != sol_2[j])
		{
			is_exist = false;
			break;
		}
	}
	return is_exist;
}

void random_perturbation(int *sol, double &f_sol, int pert_len)
{
	int count = 0;
	int number_neighbors = Instance.items*(Instance.items - 1) / 2 + Instance.items*(Instance.knaps + 1);
	do
	{
		int cur_index = rand() % number_neighbors;
		if (Neighbors[cur_index].type == 1)
		{
			int v = Neighbors[cur_index].v;
			int g = Neighbors[cur_index].g;
			if ((sol[v] != g) && satisfy_all_constraints_relocate(v, sol[v], g))
			{
				int old_g = sol[v];
				relocate_move(sol, v, g);
				count++;
			}
		}
		else if (Neighbors[cur_index].type == 2)
		{
			int x = Neighbors[cur_index].x;
			int y = Neighbors[cur_index].y;
			if ((sol[x] != sol[y]) && satisfy_all_constraints_swap(x, sol[x], y, sol[y]))
			{
				int old_g = sol[x];
				swap_move(sol, x, y);
				count++;
			}
		}
	} while (count < pert_len);
	double f_part1 = compute_obj_part1(sol);
	double f_part2 = compute_obj_part2(sol);
	f_sol = f_part1 + f_part2;
}

void allocate_memory()
{
	Size_kr = new int*[Instance.knaps + 1];
	for (int k = 0; k < Instance.knaps + 1; k++)
		Size_kr[k] = new int[Instance.classes];
	Size_r = new int[Instance.classes];
	Capacity1_k = new double[Instance.knaps + 1];
	Capacity2_k = new double[Instance.knaps + 1];
	Exceed_cap_k = new double[Instance.knaps + 1];
	Exceed_class_r = new int[Instance.classes];
	Delta_matrix_jk = new double*[Instance.items];
	Tabu_list = new int*[Instance.items];
	Proba_matrix_jk = new double*[Instance.items];
	for (int j = 0; j < Instance.items; j++)
	{
		Delta_matrix_jk[j] = new double[Instance.knaps + 1];
		Tabu_list[j] = new int[Instance.knaps + 1];
		Proba_matrix_jk[j] = new double[Instance.knaps + 1];
	}	
	Pop_sol = new int*[Pop_size];
	for (int p = 0; p < Pop_size; p++)
		Pop_sol[p] = new int[Instance.items];
	Pop_obj = new double[Pop_size];
	int number_neighbors = Instance.items*(Instance.items - 1) / 2 + Instance.items*(Instance.knaps + 1);
	Neighbors = new neighborhood[number_neighbors];
	Cur_sol = new int[Instance.items];
	Local_best_sol = new int[Instance.items];
	Global_best_sol = new int[Instance.items];
	Best_sol_one_run = new int[Instance.items];

}

void free_memory()
{
	for (int i = 0; i < Instance.knaps + 1; i++)
	{
		delete[]Size_kr[i]; Size_kr[i] = NULL;

	}
	delete[]Size_kr; Size_kr = NULL;
	delete[]Size_r; Size_r = NULL;
	delete[]Capacity1_k; Capacity1_k = NULL;
	delete[]Capacity2_k; Capacity2_k = NULL;
	delete[]Exceed_cap_k; Exceed_cap_k = NULL;
	delete[]Exceed_class_r; Exceed_class_r = NULL;
	for (int j = 0; j < Instance.items; j++)
	{
		delete[]Delta_matrix_jk[j]; Delta_matrix_jk[j] = NULL;
		delete[]Tabu_list[j]; Tabu_list[j] = NULL;
		delete[]Proba_matrix_jk[j]; Proba_matrix_jk[j] = NULL;
	}
	delete[]Delta_matrix_jk; Delta_matrix_jk = NULL;
	delete[]Tabu_list; Tabu_list = NULL;
	delete[]Proba_matrix_jk; Proba_matrix_jk = NULL;
	for (int p = 0; p < Pop_size; p++)
	{
		delete[]Pop_sol[p]; Pop_sol[p] = NULL;
	}
	delete[]Pop_sol; Pop_sol = NULL;
	delete[]Pop_obj; Pop_obj = NULL;
	delete[]Neighbors; Neighbors = NULL;
	delete[]Cur_sol; Cur_sol = NULL;
	delete[]Local_best_sol; Local_best_sol = NULL;
	delete[]Global_best_sol; Global_best_sol = NULL;
	delete[]Best_sol_one_run; Best_sol_one_run = NULL;
}

void out_results_best_sol(char *out_filename, char *instance_name, double ff, int *sol)
{
	FILE *fp;
	char buff[MAXCAN];
	sprintf(buff, "%s", out_filename);
	fp = fopen(buff, "a+");
	fprintf(fp, "%s %.6f\n", instance_name, ff);
	for (int i = 0; i < Instance.items; i++)
		fprintf(fp, "%d ", sol[i]);
	fprintf(fp, "\n");
	fclose(fp);
}

void out_results_stat(double best, double ave, double worst, double avg_time, int hit, char *stat_filename, char instance[])
{
	FILE *fp;
	char buff[MAXCAN];
	sprintf(buff, "%s", stat_filename);
	fp = fopen(buff, "a+");
	fprintf(fp, "%s %f %f %f %d %f\n", instance, best, ave, worst, hit, avg_time);
	fclose(fp);
}

/*  check the sol reported in the literature and the sol obtained by our proposed algorithm */
void verify_sol_from_file(int *sol, double ff, int type)
{
	char sol_file[] = "E:\\G-QMKP\\sol.txt";
	ifstream FIC;
	FIC.open(sol_file);
	if (FIC.fail())
	{
		cout << "### Error open, sol_file " << sol_file << endl;
		exit(0);
	}
	int item = 0, knap;
	while (!FIC.eof() && item < Instance.items)
	{		
		FIC >> knap;
		//cout << "item=" << item << ", knap=" << knap << endl;
		sol[item++] = knap;		
	}
	if (type == 1)
	{
		for (int j = 0; j < Instance.items; j++)
		{
			sol[j]--;
			if (sol[j] < 0)
				sol[j] = Instance.knaps;
		}		
	}

	for (int k = 0; k < Instance.knaps + 1; k++)
		for (int r = 0; r < Instance.classes; r++)
			Size_kr[k][r] = 0;
	for (int j = 0; j < Instance.items; j++)
		Size_kr[sol[j]][Instance.t_jr[j]]++;

	for (int k = 0; k < Instance.knaps + 1; k++)
		Capacity1_k[k] = 0.0;
	for (int j = 0; j < Instance.items; j++)
		Capacity1_k[sol[j]] += Instance.we[j];

	for (int k = 0; k < Instance.knaps + 1; k++)
		Capacity2_k[k] = 0.0;
	for (int k = 0; k < Instance.knaps + 1; k++)
		for (int r = 0; r < Instance.classes; r++)
			if (Size_kr[k][r]>0)
				Capacity2_k[k] += Instance.sr[r];

	for (int r = 0; r < Instance.classes; r++)
		Size_r[r] = 0;
	for (int r = 0; r < Instance.classes; r++)
	{
		for (int k = 0; k < Instance.knaps; k++)
		{
			if (Size_kr[k][r] > 0)
				Size_r[r]++;
		}
	}
	for (int k = 0; k < Instance.knaps; k++)
	{
		double capacity = Capacity1_k[k] + Capacity2_k[k];
		if (capacity > Instance.capac + PRECISION)
		{
			cout << "in verify func, an error is detected, capacity is exceeded, k=" << k << ",cap_real=" << capacity << ",Capacity=" << Instance.capac << endl;
			getchar();
		}
		cout << "k=" << k << ", capacity=" << capacity << endl;
	}
	for (int r = 0; r < Instance.classes; r++)
	{
		if (Size_r[r] > Instance.nr[r])
		{
			cout << "in verify func, an error is detected, nr is exceeded, r=" << r << ", Size_r=" << Size_r[r] << ", nr=" << Instance.nr[r] << endl;
			getchar();
		}
		cout << "r=" << r + 1 << ", size=" << Size_r[r] << endl;
	}
	double f_check = 0.0;
	for (int j = 0; j < Instance.items; j++)
	{
		if (sol[j] < Instance.knaps)
		{
			f_check += Instance.p_jk[j][sol[j]];
		}
	}
	for (int j = 0; j < Instance.items; j++)
	{
		for (int j2 = j + 1; j2 < Instance.items; j2++)
		{
			if (sol[j] == sol[j2] && sol[j] != Instance.knaps)
			{
				f_check += Instance.q_ij[j][j2];
			}
		}
	}
	if (fabs(f_check - ff) > PRECISION)
	{
		cout << "in verify func, an error is detected, f_check != ff, & f_check=" << f_check << ", ff=" << ff << endl;
		getchar();
	}	
	cout << "in verify_sol_from_file func, and sol is feasible, & f_sol = " << f_check << endl;
	getchar();
}

//verify whether the obtained solution is feasible or not
void proof(int *sol, double ff)
{
	for (int k = 0; k < Instance.knaps + 1; k++)
		for (int r = 0; r < Instance.classes; r++)
			Size_kr[k][r] = 0;
	for (int j = 0; j < Instance.items; j++)
		Size_kr[sol[j]][Instance.t_jr[j]]++;

	for (int k = 0; k < Instance.knaps + 1; k++)
		Capacity1_k[k] = 0.0;
	for (int j = 0; j < Instance.items; j++)
		Capacity1_k[sol[j]] += Instance.we[j];
	
	for (int k = 0; k < Instance.knaps + 1; k++)
		Capacity2_k[k] = 0.0;
	for (int k = 0; k < Instance.knaps + 1; k++)
		for (int r = 0; r < Instance.classes; r++)
			if (Size_kr[k][r]>0)
				Capacity2_k[k] += Instance.sr[r];

	for (int r = 0; r < Instance.classes; r++)
		Size_r[r] = 0;
	for (int r = 0; r < Instance.classes; r++)
	{
		for (int k = 0; k < Instance.knaps; k++)
		{
			if (Size_kr[k][r] > 0)
				Size_r[r]++;
		}
	}
	
	for (int k = 0; k < Instance.knaps; k++)
	{
		double capacity = Capacity1_k[k] + Capacity2_k[k] - Instance.capac;
		if (capacity > 0 + PRECISION)
		{
			cout << "in proof, an error is detected, capacity is exceeded, k=" << k << ",cap_real=" << Capacity1_k[k] + Capacity2_k[k] << ",Capacity=" << Instance.capac << endl;
			//getchar();
		}
		else
			capacity = 0.0;
		//cout << "k=" << k << ", capacity=" << capacity << endl;
		if (fabs(capacity - Exceed_cap_k[k]) > 0 + PRECISION)
		{
			cout << "in proof, k=" << k << ", error exceed_cap" << Exceed_cap_k[k] << endl;
			//getchar();
		}
	}	

	for (int r = 0; r < Instance.classes; r++)
	{
		int diff = Size_r[r] - Instance.nr[r];
		if (diff > 0)
		{
			cout << "in proof, an error is detected, nr is exceeded, r=" << r << ", Size_r=" << Size_r[r] << ", nr=" << Instance.nr[r] << endl;
			//getchar();
		}
		else
			diff = 0;		
	}
	double f_check = 0.0;
	for (int j = 0; j < Instance.items; j++)
	{
		if (sol[j] < Instance.knaps)
		{
			f_check += Instance.p_jk[j][sol[j]];
		}
	}
	for (int j = 0; j < Instance.items; j++)
	{
		for (int j2 = j + 1; j2 < Instance.items; j2++)
		{
			if (sol[j] == sol[j2] && sol[j] != Instance.knaps)
			{
				f_check += Instance.q_ij[j][j2];
			}
		}
	}
	if (fabs(f_check - ff) > PRECISION)
	{
		cout << "in proof, an error is detected, f_check != ff, & f_check=" << f_check << ", ff=" << ff << endl;
		getchar();
	}
}

void check_move(int *sol, double fs, double gs)
{
	double fs_check = 0.0;
	double gs_check = 0.0;
		
	for (int j = 0; j < Instance.items; j++)
	{
		if (sol[j] < Instance.knaps)
		{
			fs_check += Instance.p_jk[j][sol[j]];
		}
	}
	for (int j = 0; j < Instance.items; j++)
	{
		for (int j2 = j + 1; j2 < Instance.items; j2++)
		{
			if (sol[j] == sol[j2] && sol[j] != Instance.knaps)
			{
				fs_check += Instance.q_ij[j][j2];
			}
		}
	}
	for (int r = 0; r < Instance.classes; r++)
		if (Size_r[r] > Instance.nr[r])		
		gs_check += Size_r[r] - Instance.nr[r];		
	for (int k = 0; k < Instance.knaps; k++)
	{
		if (Capacity1_k[k] + Capacity2_k[k] > Instance.capac + PRECISION)
		{
			gs_check += Capacity1_k[k] + Capacity2_k[k] - Instance.capac;
		}
	}
	if (fabs(fs_check - fs) > PRECISION)
	{
		cout << "in check_move, an error is detected, f_check != ff, & f_check=" << fs_check << ", fs=" << fs << endl;
		getchar();
	}
	if (fabs(gs_check - gs) > PRECISION)
	{
		cout << "in check_move, an error is detected, gs_check != gs, & gs_check=" << gs_check << ", gs=" << gs << endl;
		getchar();
	}	
}
