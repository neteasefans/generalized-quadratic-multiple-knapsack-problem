sets	
	j siparis turu   /1*30/
	k knapsack indisi /1*3/
	r kalip indisi /1*15/;
alias(i,j);	
*buyuk pozitif sayi	
scalar U/30/; 	
*kapasite harcama miktarlari	
parameter w(j)/	
1	29
2	12
3	46
4	17
5	15
6	38
7	22
8	34
9	46
10	6
11	2
12	15
13	8
14	23
15      0
16      0
17      0 
18      0
19      0
20      0 
21	0
22	10
23	11
24	11
25	32
26	5
27	35
28	21
29	14
30	29
/;	
parameter cap(k);	
cap(k)= 79;	
	
*orjinal p degerleri	
parameter po(j)/	
1	77
2	294
3	151
4	318
5	131
6	5
7	118
8	72
9	367
10	128
11	172
12	119
13	369
14	401
15	476
16	436
17	444
18	106
19	46
20	165
21	379
22	309
23	1
24	365
25	185
26	119
27	307
28	256
29	351
30	433
/ ;	
parameter pp(i,j)/	
1.2= 2	
1.3= 19	
1.4= 4	
1.5= 7	
1.6= 5	
1.7= 13	
1.8= 4	
1.9= 4	
1.10= 2	
1.11= 20	
1.12= 12	
1.13= 18	
1.14= 10	
1.15= 16	
1.16= 2	
1.17= 11	
1.18= 6	
1.19= 12	
1.20= 1	
1.21= 7	
1.22= 5	
1.23= 6	
1.24= 12	
1.25= 1	
1.26= 2	
1.27= 14	
1.28= 8	
1.29= 10	
1.30= 10	
2.3= 17	
2.4= 1	
2.5= 3	
2.6= 4	
2.7= 12	
2.8= 11	
2.9= 9	
2.10= 13	
2.11= 3	
2.12= 1	
2.13= 9	
2.14= 11	
2.15= 13	
2.16= 17	
2.17= 8	
2.18= 4	
2.19= 9	
2.20= 9	
2.21= 19	
2.22= 6	
2.23= 14	
2.24= 15	
2.25= 18	
2.26= 1	
2.27= 19	
2.28= 10	
2.29= 14	
2.30= 18	
3.4= 13	
3.5= 20	
3.6= 15	
3.7= 2	
3.8= 8	
3.9= 9	
3.10= 4	
3.11= 1	
3.12= 14	
3.13= 12	
3.14= 1	
3.15= 17	
3.16= 10	
3.17= 1	
3.18= 1	
3.19= 15	
3.20= 8	
3.21= 16	
3.22= 19	
3.23= 15	
3.24= 9	
3.25= 1	
3.26= 3	
3.27= 2	
3.28= 9	
3.29= 3	
3.30= 2	
4.5= 13	
4.6= 11	
4.7= 12	
4.8= 15	
4.9= 12	
4.10= 15	
4.11= 7	
4.12= 5	
4.13= 13	
4.14= 8	
4.15= 4	
4.16= 4	
4.17= 18	
4.18= 12	
4.19= 9	
4.20= 8	
4.21= 1	
4.22= 18	
4.23= 7	
4.24= 3	
4.25= 20	
4.26= 10	
4.27= 14	
4.28= 18	
4.29= 7	
4.30= 1	
5.6= 1	
5.7= 11	
5.8= 1	
5.9= 1	
5.10= 11	
5.11= 4	
5.12= 9	
5.13= 6	
5.14= 6	
5.15= 10	
5.16= 9	
5.17= 7	
5.18= 16	
5.19= 3	
5.20= 8	
5.21= 3	
5.22= 11	
5.23= 9	
5.24= 19	
5.25= 11	
5.26= 9	
5.27= 16	
5.28= 10	
5.29= 11	
5.30= 8	
6.7= 2	
6.8= 10	
6.9= 3	
6.10= 15	
6.11= 8	
6.12= 6	
6.13= 11	
6.14= 15	
6.15= 6	
6.16= 5	
6.17= 17	
6.18= 3	
6.19= 9	
6.20= 13	
6.21= 9	
6.22= 14	
6.23= 4	
6.24= 14	
6.25= 12	
6.26= 10	
6.27= 16	
6.28= 11	
6.29= 18	
6.30= 14	
7.8= 5	
7.9= 10	
7.10= 14	
7.11= 16	
7.12= 10	
7.13= 9	
7.14= 15	
7.15= 6	
7.16= 2	
7.17= 11	
7.18= 12	
7.19= 4	
7.20= 15	
7.21= 15	
7.22= 2	
7.23= 8	
7.24= 13	
7.25= 3	
7.26= 20	
7.27= 6	
7.28= 2	
7.29= 8	
7.30= 6	
8.9= 3	
8.10= 1	
8.11= 13	
8.12= 17	
8.13= 16	
8.14= 9	
8.15= 9	
8.16= 18	
8.17= 19	
8.18= 15	
8.19= 15	
8.20= 8	
8.21= 15	
8.22= 14	
8.23= 3	
8.24= 16	
8.25= 7	
8.26= 8	
8.27= 8	
8.28= 3	
8.29= 7	
8.30= 4	
9.10= 2	
9.11= 1	
9.12= 6	
9.13= 18	
9.14= 17	
9.15= 16	
9.16= 3	
9.17= 9	
9.18= 2	
9.19= 1	
9.20= 18	
9.21= 7	
9.22= 9	
9.23= 6	
9.24= 18	
9.25= 5	
9.26= 5	
9.27= 7	
9.28= 10	
9.29= 6	
9.30= 19	
10.11= 17	
10.12= 11	
10.13= 13	
10.14= 17	
10.15= 17	
10.16= 16	
10.17= 5	
10.18= 6	
10.19= 1	
10.20= 9	
10.21= 19	
10.22= 9	
10.23= 12	
10.24= 13	
10.25= 5	
10.26= 13	
10.27= 9	
10.28= 1	
10.29= 13	
10.30= 18	
11.12= 19	
11.13= 3	
11.14= 7	
11.15= 9	
11.16= 3	
11.17= 1	
11.18= 11	
11.19= 5	
11.20= 11	
11.21= 14	
11.22= 1	
11.23= 8	
11.24= 7	
11.25= 10	
11.26= 14	
11.27= 16	
11.28= 4	
11.29= 18	
11.30= 19	
12.13= 9	
12.14= 9	
12.15= 14	
12.16= 7	
12.17= 3	
12.18= 20	
12.19= 14	
12.20= 19	
12.21= 16	
12.22= 6	
12.23= 2	
12.24= 1	
12.25= 2	
12.26= 6	
12.27= 3	
12.28= 6	
12.29= 17	
12.30= 15	
13.14= 10	
13.15= 11	
13.16= 14	
13.17= 14	
13.18= 14	
13.19= 6	
13.20= 7	
13.21= 1	
13.22= 16	
13.23= 19	
13.24= 6	
13.25= 1	
13.26= 18	
13.27= 13	
13.28= 4	
13.29= 15	
13.30= 8	
14.15= 11	
14.16= 8	
14.17= 10	
14.18= 14	
14.19= 6	
14.20= 11	
14.21= 10	
14.22= 14	
14.23= 16	
14.24= 18	
14.25= 6	
14.26= 4	
14.27= 8	
14.28= 19	
14.29= 5	
14.30= 9	
15.16= 11	
15.17= 17	
15.18= 2	
15.19= 15	
15.20= 6	
15.21= 20	
15.22= 6	
15.23= 15	
15.24= 16	
15.25= 10	
15.26= 15	
15.27= 1	
15.28= 14	
15.29= 13	
15.30= 2	
16.17= 15	
16.18= 17	
16.19= 10	
16.20= 7	
16.21= 18	
16.22= 11	
16.23= 3	
16.24= 7	
16.25= 1	
16.26= 5	
16.27= 12	
16.28= 4	
16.29= 10	
16.30= 19	
17.18= 18	
17.19= 20	
17.20= 5	
17.21= 15	
17.22= 18	
17.23= 20	
17.24= 3	
17.25= 6	
17.26= 12	
17.27= 6	
17.28= 11	
17.29= 8	
17.30= 15	
18.19= 10	
18.20= 12	
18.21= 6	
18.22= 4	
18.23= 15	
18.24= 1	
18.25= 10	
18.26= 13	
18.27= 16	
18.28= 12	
18.29= 18	
18.30= 17	
19.20= 1	
19.21= 6	
19.22= 11	
19.23= 9	
19.24= 6	
19.25= 6	
19.26= 3	
19.27= 1	
19.28= 4	
19.29= 2	
19.30= 8	
20.21= 4	
20.22= 12	
20.23= 3	
20.24= 11	
20.25= 9	
20.26= 7	
20.27= 8	
20.28= 1	
20.29= 4	
20.30= 10	
21.22= 15	
21.23= 15	
21.24= 1	
21.25= 1	
21.26= 10	
21.27= 14	
21.28= 2	
21.29= 3	
21.30= 19	
22.23= 12	
22.24= 18	
22.25= 19	
22.26= 17	
22.27= 1	
22.28= 13	
22.29= 17	
22.30= 18	
23.24= 6	
23.25= 12	
23.26= 1	
23.27= 5	
23.28= 15	
23.29= 18	
23.30= 4	
24.25= 20	
24.26= 6	
24.27= 15	
24.28= 12	
24.29= 14	
24.30= 18	
25.26= 4	
25.27= 1	
25.28= 17	
25.29= 4	
25.30= 7	
26.27= 2	
26.28= 8	
26.29= 16	
26.30= 2	
27.28= 11	
27.29= 6	
27.30= 19	
28.29= 11	
28.30= 4	
29.30= 18	
/;	
parameter t(r,j)/	
1.1= 1	
2.2= 1	
3.3= 1	
4.4= 1	
5.5= 1	
6.6= 1	
7.7= 1	
8.8= 1	
9.9= 1	
10.10= 1	
11.11= 1	
12.12= 1	
13.13= 1	
14.14= 1	
15.15= 1	
13.16= 1	
14.17= 1	
3.18= 1	
8.19= 1	
4.20= 1	
1.21= 1	
11.22= 1	
8.23= 1	
13.24= 1	
10.25= 1	
2.26= 1	
12.27= 1	
5.28= 1	
13.29= 1	
5.30= 1	
/;	
parameter s(r)/	
1	24
2	13
3	13
4	13
5	16
6	22
7	12
8	18
9	23
10	22
11	22
12	13
13	22
14	19
15	17
/;	
parameter nr(r)/	
1	1
2	1
3	1
4	1
5	1
6	1
7	1
8	1
9	1
10	1
11	1
12	1
13	1
14	1
15	1
/;	
parameter psi(r,k)/	
1.1= 0.94	
1.2= 0.76	
1.3= 0	
2.1= 0.94	
2.2= 0.76	
2.3= 0	
3.1= 0.94	
3.2= 0.76	
3.3= 0	
4.1= 0.94	
4.2= 0.76	
4.3= 0	
5.1= 0.94	
5.2= 0.76	
5.3= 0.92	
6.1= 0.94	
6.2= 0.76	
6.3= 0	
7.1= 0.94	
7.2= 0	
7.3= 0	
8.1= 0.94	
8.2= 0	
8.3= 0.92	
9.1= 0.94	
9.2= 0	
9.3= 0	
10.1= 0.94	
10.2= 0.76	
10.3= 0.92	
11.1= 0.94	
11.2= 0.76	
11.3= 0.92	
12.1= 0.94	
12.2= 0	
12.3= 0.92	
13.1= 0.94	
13.2= 0.76	
13.3= 0.92	
14.1= 0.94	
14.2= 0.76	
14.3= 0.92	
15.1= 0.94	
15.2= 0.76	
15.3= 0.92	
/;	
parameter sigma(r,k)/	
1.1= 1	
1.2= 1	
1.3= 0	
2.1= 1	
2.2= 1	
2.3= 0	
3.1= 1	
3.2= 1	
3.3= 0	
4.1= 1	
4.2= 1	
4.3= 0	
5.1= 1	
5.2= 1	
5.3= 1	
6.1= 1	
6.2= 1	
6.3= 0	
7.1= 1	
7.2= 0	
7.3= 0	
8.1= 1	
8.2= 0	
8.3= 1	
9.1= 1	
9.2= 0	
9.3= 0	
10.1= 1	
10.2= 1	
10.3= 1	
11.1= 1	
11.2= 1	
11.3= 1	
12.1= 1	
12.2= 0	
12.3= 1	
13.1= 1	
13.2= 1	
13.3= 1	
14.1= 1	
14.2= 1	
14.3= 1	
15.1= 1	
15.2= 1	
15.3= 1	
/;	
parameter p(j,k);	
p(j,k)= po(j)*sum(r,t(r,j)*psi(r,k));	
parameter epsilon(j,k);	
epsilon(j,k)= sum(r,t(r,j)*sigma(r,k)); 	