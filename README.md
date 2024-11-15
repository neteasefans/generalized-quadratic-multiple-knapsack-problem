# A hybrid evolutionary search for the generalized quadratic multiple knapsack problem (GQMKP)
This repository includes the source code and the best solution certificates of the proposed HESA algorithm published in an EJOR paper titled with "A hybrid evolutionary search for the generalized quadratic multiple knapsack problem".

The two sets of 96 benchmark instances were first introduced in Saraçand Sipahioglu (2014) and later used in Adouani et al. (2019) , Avci and Topaloglu (2017) , Chen and Hao (2016). To facilitate the further research, we upload the instances here.

We made comparisons between HESA and some state-of-the-art methods from the following related GQMKP papers:

[1] Saraç, T. , & Sipahioglu, A. (2014). Generalized quadratic multiple knapsack problem and two solution approaches. Computers & Operations Research, 43, 78–89.

[2].Chen, Y. N. , & Hao, J. K. (2016). Memetic search for the generalized quadratic mul- tiple knapsack problem. IEEE Transactions on Evolutionary Computation, 20(6), 908–923.

[3].Avci, M. , & Topaloglu, S. (2017). A multi-start iterated local search algorithm for the generalized quadratic multiple knapsack problem. Computers & Operations Research, 83 , 54–65.

[4].Adouani, Y., Jarboui, B., & Masmoudi, M. (2019). A matheuristic for the 0–1 generalized quadratic multiple knapsack problem. In Optimization letters . https: //doi.org/10.1007/s11590-019-01503-z. In press

Please cite our work as:
Zhou, Q., Hao, J. K., & Wu, Q. (2022). A hybrid evolutionary search for the generalized quadratic multiple knapsack problem. European Journal of Operational Research, 296(3), 788-803.

** Instructions to use the souce code of HESA

*** To compile:

q.zhou$ make

q.zhou$

*** To run:

q.zhou$ ./HESA_GQMKP ./instance_file ./output_stat_file ./output_sol_file seed

(where "instance_file is the instance name, output_stat_file is a file used to store the running information, output_sol_file stores the solution information, and seed is the random seed, such as 1, 2, ...")

*** To clean

q.zhou$ make clean

q.zhou$
