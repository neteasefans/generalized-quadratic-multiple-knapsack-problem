sets	
	j siparis turu   /1*30/
	k knapsack indisi /1*1/
	r kalip indisi /1*3/;
alias(i,j);	
*buyuk pozitif sayi	
scalar U/30/; 	
*kapasite harcama miktarlari	
parameter w(j)/	
1	38
2	49
3	18
4	7
5	20
6	10
7	34
8	27
9	13
10	30
11	33
12	7
13	49
14	9
15	1
16	10
17	39
18	3
19	43
20	26
21	17
22	1
23	18
24	5
25	40
26	18
27	13
28	44
29	12
30	40
/;	
parameter cap(k);	
cap(k)= 539;	
	
*orjinal p degerleri	
parameter po(j)/	
1	61
2	407
3	298
4	53
5	157
6	198
7	200
8	217
9	116
10	133
11	418
12	373
13	262
14	151
15	392
16	477
17	211
18	290
19	477
20	214
21	119
22	442
23	189
24	244
25	428
26	99
27	462
28	453
29	122
30	104
/ ;	
parameter pp(i,j)/	
1.2= 16	
1.3= 8	
1.4= 19	
1.5= 16	
1.6= 19	
1.7= 4	
1.8= 12	
1.9= 12	
1.10= 17	
1.11= 7	
1.12= 7	
1.13= 7	
1.14= 11	
1.15= 10	
1.16= 10	
1.17= 12	
1.18= 13	
1.19= 4	
1.20= 19	
1.21= 1	
1.22= 9	
1.23= 8	
1.24= 12	
1.25= 1	
1.26= 11	
1.27= 2	
1.28= 19	
1.29= 5	
1.30= 11	
2.3= 12	
2.4= 13	
2.5= 1	
2.6= 3	
2.7= 7	
2.8= 15	
2.9= 16	
2.10= 2	
2.11= 3	
2.12= 12	
2.13= 16	
2.14= 6	
2.15= 14	
2.16= 1	
2.17= 10	
2.18= 12	
2.19= 4	
2.20= 17	
2.21= 19	
2.22= 2	
2.23= 18	
2.24= 12	
2.25= 12	
2.26= 16	
2.27= 8	
2.28= 18	
2.29= 7	
2.30= 18	
3.4= 4	
3.5= 15	
3.6= 10	
3.7= 12	
3.8= 7	
3.9= 10	
3.10= 1	
3.11= 6	
3.12= 18	
3.13= 17	
3.14= 19	
3.15= 5	
3.16= 4	
3.17= 14	
3.18= 5	
3.19= 2	
3.20= 13	
3.21= 10	
3.22= 18	
3.23= 5	
3.24= 2	
3.25= 16	
3.26= 13	
3.27= 9	
3.28= 7	
3.29= 13	
3.30= 10	
4.5= 6	
4.6= 12	
4.7= 6	
4.8= 1	
4.9= 13	
4.10= 1	
4.11= 10	
4.12= 3	
4.13= 7	
4.14= 9	
4.15= 1	
4.16= 9	
4.17= 8	
4.18= 5	
4.19= 8	
4.20= 13	
4.21= 13	
4.22= 12	
4.23= 1	
4.24= 4	
4.25= 8	
4.26= 10	
4.27= 1	
4.28= 17	
4.29= 14	
4.30= 8	
5.6= 11	
5.7= 18	
5.8= 18	
5.9= 14	
5.10= 6	
5.11= 20	
5.12= 3	
5.13= 6	
5.14= 3	
5.15= 7	
5.16= 2	
5.17= 19	
5.18= 8	
5.19= 9	
5.20= 1	
5.21= 1	
5.22= 11	
5.23= 11	
5.24= 14	
5.25= 14	
5.26= 14	
5.27= 1	
5.28= 12	
5.29= 13	
5.30= 18	
6.7= 19	
6.8= 11	
6.9= 10	
6.10= 18	
6.11= 12	
6.12= 4	
6.13= 8	
6.14= 8	
6.15= 20	
6.16= 13	
6.17= 6	
6.18= 17	
6.19= 6	
6.20= 14	
6.21= 15	
6.22= 12	
6.23= 10	
6.24= 9	
6.25= 3	
6.26= 15	
6.27= 14	
6.28= 19	
6.29= 9	
6.30= 6	
7.8= 13	
7.9= 10	
7.10= 10	
7.11= 1	
7.12= 8	
7.13= 5	
7.14= 11	
7.15= 11	
7.16= 13	
7.17= 14	
7.18= 18	
7.19= 6	
7.20= 7	
7.21= 2	
7.22= 6	
7.23= 12	
7.24= 2	
7.25= 1	
7.26= 6	
7.27= 10	
7.28= 15	
7.29= 2	
7.30= 3	
8.9= 20	
8.10= 17	
8.11= 20	
8.12= 17	
8.13= 16	
8.14= 2	
8.15= 7	
8.16= 18	
8.17= 11	
8.18= 18	
8.19= 1	
8.20= 17	
8.21= 11	
8.22= 20	
8.23= 19	
8.24= 11	
8.25= 5	
8.26= 13	
8.27= 3	
8.28= 19	
8.29= 10	
8.30= 8	
9.10= 10	
9.11= 7	
9.12= 8	
9.13= 11	
9.14= 6	
9.15= 1	
9.16= 10	
9.17= 11	
9.18= 13	
9.19= 16	
9.20= 11	
9.21= 11	
9.22= 19	
9.23= 6	
9.24= 15	
9.25= 17	
9.26= 8	
9.27= 14	
9.28= 9	
9.29= 12	
9.30= 12	
10.11= 1	
10.12= 12	
10.13= 1	
10.14= 19	
10.15= 19	
10.16= 4	
10.17= 2	
10.18= 3	
10.19= 18	
10.20= 16	
10.21= 1	
10.22= 15	
10.23= 18	
10.24= 1	
10.25= 12	
10.26= 13	
10.27= 16	
10.28= 3	
10.29= 19	
10.30= 12	
11.12= 20	
11.13= 10	
11.14= 17	
11.15= 11	
11.16= 16	
11.17= 18	
11.18= 7	
11.19= 8	
11.20= 20	
11.21= 7	
11.22= 1	
11.23= 7	
11.24= 6	
11.25= 8	
11.26= 16	
11.27= 8	
11.28= 19	
11.29= 14	
11.30= 7	
12.13= 19	
12.14= 20	
12.15= 7	
12.16= 8	
12.17= 2	
12.18= 1	
12.19= 15	
12.20= 9	
12.21= 8	
12.22= 14	
12.23= 13	
12.24= 4	
12.25= 3	
12.26= 17	
12.27= 1	
12.28= 2	
12.29= 2	
12.30= 11	
13.14= 10	
13.15= 14	
13.16= 12	
13.17= 2	
13.18= 9	
13.19= 7	
13.20= 12	
13.21= 19	
13.22= 4	
13.23= 12	
13.24= 13	
13.25= 16	
13.26= 19	
13.27= 12	
13.28= 15	
13.29= 10	
13.30= 10	
14.15= 16	
14.16= 15	
14.17= 4	
14.18= 16	
14.19= 7	
14.20= 14	
14.21= 1	
14.22= 6	
14.23= 11	
14.24= 13	
14.25= 15	
14.26= 11	
14.27= 8	
14.28= 16	
14.29= 16	
14.30= 15	
15.16= 4	
15.17= 5	
15.18= 6	
15.19= 1	
15.20= 18	
15.21= 11	
15.22= 3	
15.23= 5	
15.24= 11	
15.25= 18	
15.26= 15	
15.27= 5	
15.28= 2	
15.29= 18	
15.30= 12	
16.17= 9	
16.18= 10	
16.19= 19	
16.20= 9	
16.21= 9	
16.22= 11	
16.23= 18	
16.24= 3	
16.25= 6	
16.26= 4	
16.27= 5	
16.28= 18	
16.29= 4	
16.30= 3	
17.18= 18	
17.19= 17	
17.20= 19	
17.21= 3	
17.22= 13	
17.23= 15	
17.24= 7	
17.25= 17	
17.26= 10	
17.27= 3	
17.28= 16	
17.29= 1	
17.30= 8	
18.19= 2	
18.20= 16	
18.21= 2	
18.22= 11	
18.23= 13	
18.24= 10	
18.25= 15	
18.26= 1	
18.27= 4	
18.28= 19	
18.29= 7	
18.30= 3	
19.20= 4	
19.21= 19	
19.22= 10	
19.23= 4	
19.24= 10	
19.25= 6	
19.26= 10	
19.27= 14	
19.28= 19	
19.29= 7	
19.30= 15	
20.21= 19	
20.22= 16	
20.23= 12	
20.24= 11	
20.25= 14	
20.26= 5	
20.27= 1	
20.28= 4	
20.29= 9	
20.30= 1	
21.22= 18	
21.23= 17	
21.24= 1	
21.25= 16	
21.26= 16	
21.27= 15	
21.28= 3	
21.29= 12	
21.30= 11	
22.23= 19	
22.24= 16	
22.25= 9	
22.26= 9	
22.27= 13	
22.28= 17	
22.29= 17	
22.30= 20	
23.24= 14	
23.25= 8	
23.26= 5	
23.27= 3	
23.28= 5	
23.29= 9	
23.30= 12	
24.25= 1	
24.26= 3	
24.27= 7	
24.28= 14	
24.29= 7	
24.30= 7	
25.26= 4	
25.27= 18	
25.28= 16	
25.29= 12	
25.30= 11	
26.27= 7	
26.28= 14	
26.29= 16	
26.30= 12	
27.28= 18	
27.29= 9	
27.30= 13	
28.29= 1	
28.30= 3	
29.30= 1	
/;	
parameter t(r,j)/	
1.1= 1	
2.2= 1	
3.3= 1	
3.4= 1	
3.5= 1	
3.6= 1	
3.7= 1	
3.8= 1	
3.9= 1	
3.10= 1	
3.11= 1	
3.12= 1	
3.13= 1	
3.14= 1	
3.15= 1	
3.16= 1	
3.17= 1	
3.18= 1	
3.19= 1	
3.20= 1	
3.21= 1	
3.22= 1	
3.23= 1	
3.24= 1	
3.25= 1	
3.26= 1	
3.27= 1	
3.28= 1	
3.29= 1	
3.30= 1	
/;	
parameter s(r)/	
1	11
2	8
3	5
/;	
parameter nr(r)/	
1	5
2	3
3	9
/;	
parameter psi(r,k)/	
1.1= 0.61	
2.1= 0.61	
3.1= 0.61	
/;	
parameter sigma(r,k)/	
1.1= 1	
2.1= 1	
3.1= 1	
/;	
parameter p(j,k);	
p(j,k)= po(j)*sum(r,t(r,j)*psi(r,k));	
parameter epsilon(j,k);	
epsilon(j,k)= sum(r,t(r,j)*sigma(r,k)); 	