sets	
	j siparis turu   /1*30/
	k knapsack indisi /1*3/
	r kalip indisi /1*15/;
alias(i,j);	
*buyuk pozitif sayi	
scalar U/30/; 	
*kapasite harcama miktarlari	
parameter w(j)/	
1	6
2	4
3	8
4	2
5	8
6	4
7	6
8	4
9	6
10	2
11	9
12	3
13	8
14	5
15	2
16	7
17	1
18	3
19	4
20	7
21	2
22	6
23	2
24	6
25	3
26	4
27	8
28	3
29	2
30	8
/;	
parameter cap(k);	
cap(k)= 38;	

*orjinal p degerleri	
parameter po(j)/	
1	0
2	0
3	0
4	0
5	0
6	0
7	70
8	50
9	0
10	0
11	0
12	0
13	0
14	0
15	0
16	80
17	0
18	0
19	0
20	0
21	0
22	0
23	0
24	0
25	40
26	50
27	90
28	0
29	0
30	0
/ ;	
parameter pp(i,j)/	
1.2= 0	
1.3= 0	
1.4= 0	
1.5= 0
1.6= 0
1.7= 6
1.8= 0
1.9= 0
1.10= 0
1.11= 0
1.12= 0
1.13= 0
1.14= 0
1.15= 6
1.16= 8
1.17= 0
1.18= 0
1.19= 0
1.20= 0
1.21= 0
1.22= 0
1.23= 0
1.24= 6
1.25= 0
1.26= 0
1.27= 8
1.28= 0
1.29= 0
1.30= 0
2.3= 0
2.4= 0
2.5= 6
2.6= 0
2.7= 6
2.8= 0
2.9= 0
2.10= 0
2.11= 0
2.12= 0
2.13= 6
2.14= 0
2.15= 0
2.16= 0
2.17= 0
2.18= 0
2.19= 0
2.20= 0
2.21= 0
2.22= 6
2.23= 4
2.24= 0
2.25= 4
2.26= 6
2.27= 0
2.28= 0
2.29= 4
2.30= 0
3.4= 0
3.5= 0
3.6= 0
3.7= 0
3.8= 6
3.9= 0
3.10= 0
3.11= 0
3.12= 0
3.13= 8
3.14= 8
3.15= 6
3.16= 8
3.17= 0
3.18= 0
3.19= 6
3.20= 8
3.21= 6
3.22= 8
3.23= 0
3.24= 0
3.25= 0
3.26= 0
3.27= 0
3.28= 6
3.29= 0
3.30= 0
4.5= 0
4.6= 0
4.7= 6
4.8= 0
4.9= 0
4.10= 0
4.11= 0
4.12= 0
4.13= 0
4.14= 0
4.15= 0
4.16= 6
4.17= 0
4.18= 4
4.19= 0
4.20= 6
4.21= 4
4.22= 0
4.23= 0
4.24= 0
4.25= 4
4.26= 0
4.27= 0
4.28= 4
4.29= 0
4.30= 0
5.6= 6
5.7= 0
5.8= 6
5.9= 0
5.10= 0
5.11= 0
5.12= 0
5.13= 0
5.14= 0
5.15= 6
5.16= 8
5.17= 6
5.18= 0
5.19= 0
5.20= 0
5.21= 0
5.22= 0
5.23= 0
5.24= 0
5.25= 0
5.26= 0
5.27= 0
5.28= 0
5.29= 0
5.30= 0
6.7= 6
6.8= 0
6.9= 0
6.10= 4
6.11= 8
6.12= 0
6.13= 6
6.14= 0
6.15= 4
6.16= 0
6.17= 0
6.18= 4
6.19= 6
6.20= 6
6.21= 4
6.22= 0
6.23= 4
6.24= 6
6.25= 0
6.26= 6
6.27= 0
6.28= 0
6.29= 0
6.30= 0
7.8= 6
7.9= 6
7.10= 0
7.11= 8
7.12= 0
7.13= 0
7.14= 0
7.15= 0
7.16= 8
7.17= 0
7.18= 0
7.19= 0
7.20= 0
7.21= 0
7.22= 0
7.23= 6
7.24= 6
7.25= 6
7.26= 6
7.27= 0
7.28= 0
7.29= 0
7.30= 0
8.9= 0
8.10= 0
8.11= 0
8.12= 0
8.13= 0
8.14= 0
8.15= 0
8.16= 6
8.17= 0
8.18= 0
8.19= 0
8.20= 0
8.21= 0
8.22= 0
8.23= 4
8.24= 0
8.25= 0
8.26= 0
8.27= 0
8.28= 0
8.29= 0
8.30= 0
9.10= 0
9.11= 0
9.12= 0
9.13= 8
9.14= 0
9.15= 0
9.16= 0
9.17= 0
9.18= 0
9.19= 6
9.20= 8
9.21= 0
9.22= 6
9.23= 0
9.24= 6
9.25= 6
9.26= 0
9.27= 0
9.28= 0
9.29= 0
9.30= 0
10.11= 6
10.12= 0
10.13= 0
10.14= 0
10.15= 0
10.16= 6
10.17= 0
10.18= 0
10.19= 4
10.20= 6
10.21= 0
10.22= 0
10.23= 0
10.24= 0
10.25= 0
10.26= 0
10.27= 0
10.28= 0
10.29= 0
10.30= 0
11.12= 6
11.13= 0
11.14= 0
11.15= 6
11.16= 0
11.17= 0
11.18= 0
11.19= 0
11.20= 8
11.21= 0
11.22= 0
11.23= 0
11.24= 0
11.25= 6
11.26= 8
11.27= 8
11.28= 0
11.29= 0
11.30= 0
12.13= 0
12.14= 0
12.15= 0
12.16= 0
12.17= 0
12.18= 0
12.19= 0
12.20= 6
12.21= 0
12.22= 0
12.23= 0
12.24= 6
12.25= 0
12.26= 0
12.27= 0
12.28= 0
12.29= 0
12.30= 0
13.14= 8
13.15= 0
13.16= 8
13.17= 6
13.18= 0
13.19= 0
13.20= 8
13.21= 0
13.22= 0
13.23= 0
13.24= 0
13.25= 6
13.26= 0
13.27= 8
13.28= 0
13.29= 0
13.30= 8
14.15= 0
14.16= 0
14.17= 0
14.18= 0
14.19= 6
14.20= 0
14.21= 4
14.22= 0
14.23= 0
14.24= 0
14.25= 0
14.26= 0
14.27= 0
14.28= 0
14.29= 4
14.30= 0
15.16= 0
15.17= 4
15.18= 0
15.19= 0
15.20= 0
15.21= 4
15.22= 6
15.23= 0
15.24= 0
15.25= 0
15.26= 0
15.27= 0
15.28= 0
15.29= 0
15.30= 0
16.17= 0
16.18= 0
16.19= 0
16.20= 0
16.21= 0
16.22= 8
16.23= 6
16.24= 0
16.25= 6
16.26= 6
16.27= 8
16.28= 6
16.29= 0
16.30= 0
17.18= 4
17.19= 0
17.20= 0
17.21= 0
17.22= 0
17.23= 4
17.24= 0
17.25= 0
17.26= 0
17.27= 6
17.28= 0
17.29= 0
17.30= 0
18.19= 0
18.20= 0
18.21= 0
18.22= 0
18.23= 0
18.24= 6
18.25= 0
18.26= 0
18.27= 6
18.28= 0
18.29= 4
18.30= 0
19.20= 0
19.21= 0
19.22= 0
19.23= 0
19.24= 6
19.25= 4
19.26= 0
19.27= 0
19.28= 0
19.29= 0
19.30= 0
20.21= 6
20.22= 0
20.23= 0
20.24= 8
20.25= 6
20.26= 0
20.27= 8
20.28= 0
20.29= 6
20.30= 0
21.22= 6
21.23= 4
21.24= 0
21.25= 4
21.26= 4
21.27= 0
21.28= 4
21.29= 4
21.30= 0
22.23= 0
22.24= 0
22.25= 0
22.26= 6
22.27= 0
22.28= 0
22.29= 0
22.30= 0
23.24= 0
23.25= 0
23.26= 0
23.27= 0
23.28= 0
23.29= 0
23.30= 0
24.25= 0
24.26= 0
24.27= 0
24.28= 0
24.29= 6
24.30= 0
25.26= 0
25.27= 6
25.28= 0
25.29= 0
25.30= 0
26.27= 0
26.28= 0
26.29= 0
26.30= 0
27.28= 0
27.29= 0
27.30= 8
28.29= 0
28.30= 6
29.30= 6
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
11.16= 1
6.17= 1
11.18= 1
4.19= 1
7.20= 1
7.21= 1
5.22= 1
15.23= 1
3.24= 1
8.25= 1
12.26= 1
11.27= 1
11.28= 1
10.29= 1
11.30= 1
/;	
parameter s(r)/	
1	1
2	1
3	2
4	2
5	2
6	2
7	2
8	2
9	2
10	1
11	1
12	2
13	2
14	2
15	2
/;	
parameter nr(r)/	
1	2
2	2
3	2
4	2
5	2
6	2
7	2
8	2
9	2
10	2
11	2
12	2
13	2
14	2
15	2
/;	
parameter psi(r,k)/	
1.1= 0.93	
1.2= 0.07	
1.3= 0	
2.1= 0.01	
2.2= 0	
2.3= 0	
3.1= 0.44	
3.2= 0.18	
3.3= 0.25	
4.1= 0.08	
4.2= 0	
4.3= 0	
5.1= 1
5.2= 0.6
5.3= 0.03
6.1= 0.49
6.2= 0.56
6.3= 0.92
7.1= 0.29
7.2= 0.28
7.3= 0.75
8.1= 0.94
8.2= 0.75
8.3= 0
9.1= 0.24
9.2= 0.9
9.3= 0.79
10.1= 0.02
10.2= 0.34
10.3= 0.26
11.1= 0.01
11.2= 0
11.3= 0.84
12.1= 0.9
12.2= 0.92
12.3= 0.73
13.1= 0.32
13.2= 0
13.3= 0.7
14.1= 0.57
14.2= 0.53
14.3= 0.88
15.1= 0.3
15.2= 0.11
15.3= 0.14
/;
parameter sigma(r,k)/
1.1= 1
1.2= 1
1.3= 0
2.1= 1
2.2= 0
2.3= 0
3.1= 1
3.2= 1
3.3= 1
4.1= 1
4.2= 0
4.3= 0
5.1= 1
5.2= 1
5.3= 1
6.1= 1
6.2= 1
6.3= 1
7.1= 1
7.2= 1
7.3= 1
8.1= 1
8.2= 1
8.3= 0
9.1= 1
9.2= 1
9.3= 1
10.1= 1
10.2= 1
10.3= 1
11.1= 1
11.2= 0
11.3= 1
12.1= 1
12.2= 1
12.3= 1
13.1= 1
13.2= 0
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