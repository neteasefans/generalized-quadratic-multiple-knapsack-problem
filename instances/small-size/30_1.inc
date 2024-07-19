sets	
	j siparis turu   /1*30/
	k knapsack indisi /1*3/
	r kalip indisi /1*3/;
alias(i,j);	
*buyuk pozitif sayi	
scalar U/30/; 	
*kapasite harcama miktarlari	
parameter w(j)/	
1	6
2	3
3	8
4	2
5	6
6	2
7	9
8	1
9	4
10	5
11	3
12	9
13	10
14	5
15	7
16	4
17	7
18	9
19	1
20	6
21	2
22	8
23	5
24	2
25	7
26	9
27	5
28	8
29	5
30	5
/;	
parameter cap(k);	
cap(k)= 21;	

*orjinal p degerleri	
parameter po(j)/	
1	58
2	80
3	93
4	76
5	88
6	7
7	93
8	69
9	28
10	14
11	38
12	20
13	72
14	73
15	80
16	87
17	40
18	39
19	13
20	36
21	17
22	1
23	30
24	89
25	13
26	27
27	67
28	77
29	5
30	23
/ ;	
parameter pp(i,j)/	
1.2= 2	
1.3= 4	
1.4= 4	
1.5= 3
1.6= 3
1.7= 1
1.8= 2
1.9= 1
1.10= 4
1.11= 1
1.12= 2
1.13= 1
1.14= 2
1.15= 3
1.16= 2
1.17= 1
1.18= 1
1.19= 3
1.20= 3
1.21= 3
1.22= 3
1.23= 3
1.24= 1
1.25= 4
1.26= 4
1.27= 1
1.28= 3
1.29= 1
1.30= 2
2.3= 2
2.4= 4
2.5= 2
2.6= 2
2.7= 3
2.8= 1
2.9= 4
2.10= 2
2.11= 4
2.12= 4
2.13= 2
2.14= 1
2.15= 4
2.16= 2
2.17= 3
2.18= 1
2.19= 1
2.20= 2
2.21= 3
2.22= 4
2.23= 1
2.24= 4
2.25= 1
2.26= 2
2.27= 1
2.28= 1
2.29= 2
2.30= 1
3.4= 3
3.5= 2
3.6= 1
3.7= 2
3.8= 1
3.9= 3
3.10= 3
3.11= 1
3.12= 2
3.13= 2
3.14= 2
3.15= 1
3.16= 1
3.17= 2
3.18= 2
3.19= 1
3.20= 3
3.21= 4
3.22= 2
3.23= 1
3.24= 2
3.25= 2
3.26= 2
3.27= 1
3.28= 2
3.29= 2
3.30= 2
4.5= 1
4.6= 2
4.7= 2
4.8= 3
4.9= 4
4.10= 2
4.11= 2
4.12= 1
4.13= 2
4.14= 3
4.15= 1
4.16= 3
4.17= 2
4.18= 4
4.19= 2
4.20= 3
4.21= 1
4.22= 4
4.23= 3
4.24= 3
4.25= 1
4.26= 1
4.27= 3
4.28= 2
4.29= 2
4.30= 1
5.6= 1
5.7= 4
5.8= 2
5.9= 3
5.10= 3
5.11= 2
5.12= 2
5.13= 2
5.14= 2
5.15= 1
5.16= 3
5.17= 4
5.18= 1
5.19= 3
5.20= 1
5.21= 3
5.22= 1
5.23= 1
5.24= 3
5.25= 2
5.26= 4
5.27= 1
5.28= 2
5.29= 3
5.30= 2
6.7= 1
6.8= 4
6.9= 1
6.10= 1
6.11= 3
6.12= 2
6.13= 1
6.14= 2
6.15= 1
6.16= 4
6.17= 2
6.18= 2
6.19= 2
6.20= 1
6.21= 3
6.22= 1
6.23= 1
6.24= 3
6.25= 1
6.26= 3
6.27= 3
6.28= 2
6.29= 3
6.30= 1
7.8= 2
7.9= 1
7.10= 3
7.11= 1
7.12= 2
7.13= 4
7.14= 3
7.15= 1
7.16= 1
7.17= 1
7.18= 1
7.19= 4
7.20= 1
7.21= 1
7.22= 1
7.23= 1
7.24= 1
7.25= 3
7.26= 3
7.27= 2
7.28= 4
7.29= 1
7.30= 1
8.9= 1
8.10= 2
8.11= 1
8.12= 3
8.13= 1
8.14= 3
8.15= 4
8.16= 1
8.17= 1
8.18= 1
8.19= 2
8.20= 2
8.21= 3
8.22= 2
8.23= 1
8.24= 1
8.25= 1
8.26= 2
8.27= 4
8.28= 3
8.29= 3
8.30= 1
9.10= 1
9.11= 1
9.12= 3
9.13= 3
9.14= 1
9.15= 3
9.16= 3
9.17= 4
9.18= 1
9.19= 1
9.20= 1
9.21= 3
9.22= 1
9.23= 2
9.24= 1
9.25= 3
9.26= 2
9.27= 1
9.28= 2
9.29= 1
9.30= 1
10.11= 1
10.12= 3
10.13= 1
10.14= 4
10.15= 2
10.16= 2
10.17= 1
10.18= 2
10.19= 2
10.20= 4
10.21= 3
10.22= 2
10.23= 2
10.24= 3
10.25= 3
10.26= 3
10.27= 3
10.28= 2
10.29= 3
10.30= 3
11.12= 1
11.13= 3
11.14= 1
11.15= 2
11.16= 4
11.17= 3
11.18= 1
11.19= 2
11.20= 2
11.21= 2
11.22= 3
11.23= 1
11.24= 1
11.25= 2
11.26= 2
11.27= 1
11.28= 3
11.29= 2
11.30= 4
12.13= 1
12.14= 1
12.15= 3
12.16= 3
12.17= 2
12.18= 1
12.19= 3
12.20= 3
12.21= 1
12.22= 1
12.23= 1
12.24= 2
12.25= 3
12.26= 1
12.27= 1
12.28= 1
12.29= 3
12.30= 1
13.14= 1
13.15= 2
13.16= 2
13.17= 1
13.18= 1
13.19= 1
13.20= 2
13.21= 3
13.22= 4
13.23= 1
13.24= 1
13.25= 1
13.26= 2
13.27= 1
13.28= 1
13.29= 3
13.30= 3
14.15= 3
14.16= 1
14.17= 3
14.18= 3
14.19= 2
14.20= 1
14.21= 1
14.22= 3
14.23= 3
14.24= 3
14.25= 2
14.26= 4
14.27= 2
14.28= 3
14.29= 3
14.30= 3
15.16= 2
15.17= 1
15.18= 1
15.19= 1
15.20= 1
15.21= 1
15.22= 3
15.23= 3
15.24= 2
15.25= 1
15.26= 4
15.27= 3
15.28= 1
15.29= 4
15.30= 2
16.17= 4
16.18= 2
16.19= 1
16.20= 2
16.21= 3
16.22= 3
16.23= 2
16.24= 3
16.25= 2
16.26= 2
16.27= 4
16.28= 2
16.29= 3
16.30= 2
17.18= 1
17.19= 2
17.20= 3
17.21= 2
17.22= 1
17.23= 1
17.24= 3
17.25= 1
17.26= 2
17.27= 3
17.28= 3
17.29= 2
17.30= 3
18.19= 2
18.20= 2
18.21= 1
18.22= 1
18.23= 4
18.24= 2
18.25= 2
18.26= 2
18.27= 1
18.28= 2
18.29= 2
18.30= 1
19.20= 3
19.21= 4
19.22= 1
19.23= 3
19.24= 1
19.25= 1
19.26= 4
19.27= 2
19.28= 4
19.29= 1
19.30= 1
20.21= 3
20.22= 1
20.23= 3
20.24= 1
20.25= 1
20.26= 4
20.27= 2
20.28= 1
20.29= 3
20.30= 1
21.22= 3
21.23= 2
21.24= 4
21.25= 4
21.26= 4
21.27= 1
21.28= 1
21.29= 3
21.30= 4
22.23= 4
22.24= 3
22.25= 3
22.26= 1
22.27= 3
22.28= 4
22.29= 1
22.30= 2
23.24= 3
23.25= 1
23.26= 1
23.27= 2
23.28= 1
23.29= 1
23.30= 2
24.25= 2
24.26= 3
24.27= 3
24.28= 2
24.29= 1
24.30= 2
25.26= 2
25.27= 4
25.28= 1
25.29= 1
25.30= 3
26.27= 1
26.28= 1
26.29= 2
26.30= 1
27.28= 3
27.29= 1
27.30= 2
28.29= 3
28.30= 2
29.30= 1
/;
parameter t(r,j)/
1.1= 1
2.2= 1
3.3= 1
2.4= 1
3.5= 1
3.6= 1
1.7= 1
2.8= 1
2.9= 1
2.10= 1
3.11= 1
1.12= 1
3.13= 1
1.14= 1
2.15= 1
1.16= 1
2.17= 1
3.18= 1
3.19= 1
1.20= 1
2.21= 1
2.22= 1
2.23= 1
3.24= 1
2.25= 1
2.26= 1
2.27= 1
2.28= 1
3.29= 1
3.30= 1
/;	
parameter s(r)/	
1	2
2	2
3	4
/;	
parameter nr(r)/	
1	5
2	3
3	9
/;	
parameter psi(r,k)/	
1.1= 0.68	
1.2= 0.2	
1.3= 0.42	
2.1= 0.41	
2.2= 0.65
2.3= 0.77
3.1= 0.54
3.2= 0.77
3.3= 0.62
/;
parameter sigma(r,k)/
1.1= 1
1.2= 1
1.3= 1
2.1= 1
2.2= 1
2.3= 1
3.1= 1
3.2= 1
3.3= 1
/;
parameter p(j,k);
p(j,k)= po(j)*sum(r,t(r,j)*psi(r,k));
parameter epsilon(j,k);
epsilon(j,k)= sum(r,t(r,j)*sigma(r,k)); 