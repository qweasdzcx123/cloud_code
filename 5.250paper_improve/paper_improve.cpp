/*************************************************************************************\
 *                                                                                   *
 *  This procedure demonstrates the use of the DQPSO* algorithm for solving          * 
 *  the multidimensional knapsack problem.                                           *
 *  Reference: https://doi.org/10.1016/j.eswa.2020.113310                            *
 *                                                                                   *
\*************************************************************************************/ 

/*************************************************************************************/
/****  0. Header files, data structures, and global varialbes  ***********************/
/*************************************************************************************/ 
//改进算法 100-500-3 测试
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <ctime>
#include <vector>
#include<math.h>
#include<ctype.h>
using namespace std;

typedef struct{
	int *X; 
	int *S;
	int *NS;
	int IN, ON;
	int *IW; 
	long int f; 
} Solution; // Representation of solutions

typedef struct Neighbor{
	int f; 
    int type;
    int IO; 
    int k; 
    int x;
    int x1;
    int y1; 
} Neighbor; // Representation of neighbors

char * File_Name; // File of benchmark instance 
char * outfilename; // Output file 
double time_limit, time_to_target, AvgTime; 
double time_one_run, starting_time;

int f, f_best;
double BestResult, AvgResult, WorResult;  
double sigmma;   
int Nhit;  
int *P;
int *P1;  
int *B;        
int **R;  // weights matrix
int **R1; 
int N, M; // number of terms, number of dimensions
Solution SC;
Solution S_BEST; 
Solution G_BEST; 
Solution SS_BEST; 

const int np = 1000;  // Size of population, which is set according to the size of instances
double Prob_LS;       // Probability of applying the VND procedure after the repair operator
Solution Pbest[np];   // Personal historical best positions of discrete particle swarm 

double *PR; 
double *PR1; 
int *order; 
double *dataM;
double *dataNM;
double **QPOP; //Quantum Particle Swarm
double *PY;
double *GY; 
double *PFit;
int *orderPSO; 

/***********************************************************************************/
/***********************          1. Initializing           ************************/
/***********************************************************************************/ 
//a. Inputing 
void Initializing()
{
     int i,j, x1, x2; 
     int count = 0;  
            
	 ifstream FIC; 
	 ofstream FIC2;  
     FIC.open(File_Name); 
     if ( FIC.fail() )
     {
           cout << "### Erreur open, File_Name " << File_Name << endl;
           exit(0);
     } 
	 FIC >> N >> M;
 
	 P = new int [N];
	 P1= new int [N]; 
	 B = new int [N];
	 
	 R = new int *[M];
	 for(i=0;i<M;i++) 
		 R[i] = new int [N];
	 R1 = new int *[M];
	 for(i=0;i<M;i++) 
		 R1[i] = new int [N];	
		  
     while ( ! FIC.eof() )
     {
         for(i=0;i<N;i++) FIC >> P[i];
		   
		 for(i=0;i<M;i++) 
		    for(j=0;j<N;j++)FIC >> R[i][j];
	     for(j=0;j<M;j++) FIC >> B[j];  
     }	  
     FIC.close();
}
void Initializing1()
{
     int i,j, x1, x2; 
     int count=0; 

	 ifstream FIC;
	 ofstream FIC2; 
     FIC.open(File_Name);
     if ( FIC.fail() )
     {
           cout << "### Erreur open, File_Name " << File_Name << endl;
           exit(0);
     } 
	 FIC >> N >> M;
 
	 P = new int [N];
	 P1= new int [N]; 
	 B = new int [N];
	 
	 R = new int *[M];
	 for(i=0;i<M;i++) 
		 R[i] = new int [N];
	 R1 = new int *[M];
	 for(i=0;i<M;i++) 
		 R1[i] = new int [N];	
		  
     while ( ! FIC.eof() )
     {
         for(i=0;i<N;i++) FIC >> P[i];
         for(j=0;j<M;j++) FIC >> B[j];   
		 for(i=0;i<M;i++) 
		    for(j=0;j<N;j++)FIC >> R[i][j];    
     }
     
     FIC.close();
}
// b. Assignment of memery
void AssignMemery()
{
     int i, j, swap;
	 SC.X      = new int [N];
	 SC.S      = new int [N];
	 SC.NS     = new int [N]; 
	 SC.IW     = new int [M];
	 S_BEST.X  = new int [N];
     S_BEST.S  = new int [N];
	 S_BEST.NS = new int [N];
	 S_BEST.IW = new int [M]; 
	 SS_BEST.X  = new int [N];
     SS_BEST.S  = new int [N];
	 SS_BEST.NS = new int [N];
	 SS_BEST.IW = new int [M]; 
	 G_BEST.X  = new int [N];
     G_BEST.S  = new int [N];
     G_BEST.NS = new int [N];
     G_BEST.IW = new int [M];
     
     QPOP = new double *[np];
     for(i=0;i<np;i++) QPOP[i] = new double [N]; 
     PY = new double [N];
     GY = new double [N]; 
     PFit = new double [np]; 
     orderPSO = new int [np]; 
      
	 for(i=0;i<np;i++)
     {
     	Pbest[i].X  = new int [N];
     	Pbest[i].S  = new int [N];
     	Pbest[i].NS = new int [N];
     	Pbest[i].IW = new int [M]; 
	 } 
	 
     PR = new double [N]; 
     PR1= new double [N]; 
     order = new int [N];      
	 
	 dataM  = new double [N];
     dataNM = new double [N]; 
}

/*****************************************************************************/
/***** 2. Hamming distance  between two discrete solutions   *****************/
/*****************************************************************************/ 
int distance(int S1[], int S2[])
{
	int j, count;
	count =0; 
	for(j=0;j<N;j++)
	if(S1[j] != S2[j]) count ++; 
    return count; 
}

/*****************************************************************************/
/***************************      3. QuickSort     ***************************/
/*****************************************************************************/
void q_sort(double numbers[], int left, int right, int index[]) 
{
  int pivot, pivot_index, l_hold, r_hold;
  l_hold = left;
  r_hold = right;
  pivot = numbers[left];
  pivot_index = index[left] ;
  while (left < right)
  {
    while ((numbers[right] >= pivot) && (left < right))
      right--;
    if (left != right)
    {
      numbers[left] = numbers[right];
      index[left] = index[right];
      left++;
    }
    while ((numbers[left] <= pivot) && (left < right))
      left++;
    if (left != right)
    {
      numbers[right] = numbers[left];
      index[right] = index[left];
      right--;
    }
  }
  numbers[left] = pivot;
  index[left] = pivot_index;
  pivot = left;
  left = l_hold;
  right = r_hold;
  if (left < pivot)
    q_sort(numbers, left, pivot-1, index);
  if (right > pivot)
    q_sort(numbers, pivot+1, right, index);
}
void Quick_Sort(double numbers[], int array_size, int index[])
{
  q_sort(numbers, 0, array_size - 1, index);
}

/*****************************************************************************/
/*******************        4. Proof of the solution            **************/
/*****************************************************************************/ 
void proof(Solution &S)
{
   int i,j; 
   int count,p; 
   printf("\n %d \n",S.f) ;
// for(i=0;i<M;i++) printf("%d ",S.IW[i]); printf("\n"); 
// for(i=0;i<M;i++) printf("%d ",B[i]); printf("\n"); 
// printf("%d   %d   %d\n",S.IN,S.ON,S.f);
   for(i=0;i<N;i++) printf("%d ",S.X[i]); printf("\n");  
   for(i=0;i<M;i++) S.IW[i] = 0;
   for(i=0;i<M;i++)
   for(j=0;j<N;j++)
   S.IW[i] += S.X[j]*R[i][j];
  // printf("\n"); 
   for(i=0;i<M;i++)
   {  
      if(S.IW[i]>B[i]) printf("errer !") ;
	  //printf("%d  %d\n", S.IW[i],B[i]); 
   } 
   printf("\n"); 
   count = 0;
   p=0;
   for(i=0;i<N;i++) if(S.X[i]==1) count++;
   for(i=0;i<N;i++) if(S.X[i]==1) p+=P[i];
   printf("count = %d   p=%d \n",count,p);
}

/*****************************************************************************/
/****************************   5. Repair Operator  **************************/
/*****************************************************************************/
void RepairOperator(Solution &S)
{
	int i,j;
	int assign; 
	int c1,c2; 
	
	for(i=0;i<M;i++) S.IW[i] = 0;
    for(i=0;i<M;i++)
       for(j=0;j<N;j++)
         S.IW[i] += S.X[j]*R[i][j];
    
    for(j=0;j<N;j++) 
    {
       if(S.X[j] == 1)
       {
		assign = 0; 
    	for(i=0;i<M;i++)
    	{
		   if(S.IW[i] > B[i]) { assign = 1; break; }
	    }
		if(assign == 1)
		{
			S.X[j] = 0;
			for(i=0;i<M;i++) S.IW[i] -= R[i][j]; 
		}
	  }
	}
	
	for(j=N-1;j>=0;j--)
	{
		if(S.X[j] == 0)
		{
		  assign = 0;
		  for(i=0;i<M;i++)
		  {
		    if(S.IW[i] + R[i][j] > B[i] ) { assign = 1;  break;}
	      }
		  
		  if(assign == 0) 
		  {
		 	S.X[j] = 1; 
		 	for(i=0;i<M;i++) S.IW[i] += R[i][j]; 
		  } 
		   
	    }
	}
	
	S.f = 0;
	for(j=0;j<N;j++) if( S.X[j]==1 ) S.f += P[j];  
	
	c1 = 0;
	c2 = 0;
	S.IN = 0; 
	S.ON = 0; 
	for(j=0;j<N;j++)
	{
		if(S.X[j] == 1) 
		{
			S.S[c1] = j;
			c1++;
			S.IN ++; 
		}
		else
		{
			S.NS[c2] = j;
			c2++;
			S.ON ++; 
		} 
	}
 
} 

int CheckOperator(Solution &S)
{
	int i,j;
	int assign; 
	int c1,c2; 
	
	for(i=0;i<M;i++) S.IW[i] = 0;
    for(i=0;i<M;i++)
       for(j=0;j<N;j++)
         S.IW[i] += S.X[j]*R[i][j];
    
    for(i=0;i<M;i++) 
	if(S.IW[i] > B[i]) 
	{
		S.f = -99999; printf("error !\n"); return S.f; 
	} 
	
	S.f = 0;
	for(j=0;j<N;j++) if( S.X[j]==1 ) S.f += P[j];  

	c1 = 0;
	c2 = 0;
	S.IN = 0; 
	S.ON = 0; 
	for(j=0;j<N;j++)
	{
		if(S.X[j] == 1) 
		{
			S.S[c1] = j;
			c1++;
			S.IN ++; 
		}
		else
		{
			S.NS[c2] = j;
			c2++;
			S.ON ++; 
		} 
	}
 return S.f; 
} 

/*****************************************************************************/
/*********************     6. Local Search Method          *******************/
/*****************************************************************************/ 
//-----------------------------------------------------------------------------
// Variable neighoborhood search (VND) method
//-----------------------------------------------------------------------------
void VND(Solution &S)
{
	
	 int i, j, m, ci, assign, swap;
	 int Flag_VND, Flag_Swap;
	 int x,y, k, countor; 
	 int f_c;
	
		 
	Flag_VND = 1; 
	while(Flag_VND) 
	 {
	 	  Flag_VND = 0; 
	 	  for(i=0; i<S.IN; i++)
	      {
		    dataM[i] = 1.0*PR[S.S[i]]; 
	      }
	      Quick_Sort(dataM, S.IN, S.S);
		  
	      for(i=0;i<S.ON;i++)
	      {
	    	dataNM[i] = -1.0*PR[S.NS[i]]; 
	      }
	      Quick_Sort(dataNM, S.ON, S.NS); 
	      
	 	  //a.  neighborhood search for the ADD neighborhod N1
	 	  label:
	 	  for(i=0;i<S.ON;i++)
	 	  	{
			   ci = S.NS[i]; 
	 	  	   assign = 0; 
	 	  	   for(m=0;m<M;m++)
	 	  	   {
	 	  		 if(S.IW[m] + R[m][ci] > B[m] ) { assign = 1;  break;}
			   }
			   
	 	       if(assign == 0) 
		       {
		 	      S.X[ci] = 1; 
		 	      S.f += P[ci];
		 	      for(m=0;m<M;m++) S.IW[m] += R[m][ci]; 
		 	      swap = S.NS[i]; 
		 	      S.NS[i] = S.NS[S.ON - 1];
		 	      S.NS[S.ON - 1] = swap ; 
		 	      S.S[S.IN] = S.NS[S.ON - 1] ;
		 	      S.IN ++;
		 	      S.ON --; 
		 	      Flag_VND = 1;
		 	     
		 	   }
		    }
		    
	 	 // neighborhood search for the SWAP neighborhood  N2
	 	 Flag_Swap = 0; 
         for(x = 0; x < S.IN; x++)
	     { 
           for(y = 0; y < S.ON; y ++)	
			  {
				
				f_c = S.f + (P[S.NS[y]]-P[S.S[x]]) ; 
			    if(f_c <= S.f) continue; 
			    
			    assign = 0; 
			    for(m=0; m<M; m++) if(S.IW[m]+(R[m][S.NS[y]]-R[m][S.S[x]]) > B[m]) { assign  =  1;  break;  }
			    if(assign == 1) continue; 
			    
			    for(m=0;m<M;m++) 
                 {
                    S.IW[m] += (R[m][S.NS[y]] - R[m][S.S[x]]);   
				 } 
				 
                S.X[S.S[x]]  = 1 - S.X[S.S[x]]; 
                S.X[S.NS[y]] = 1 - S.X[S.NS[y]];
                       
                swap    = S.S[x] ; 
                S.S[x]  = S.NS[y];
                S.NS[y] = swap ;  
				 
                S.f = f_c;  
                Flag_Swap = 1;
                Flag_VND  = 1; 
               
                goto label; 
	          }
       	} 	
	      	    
    }
}

/*****************************************************************************/
/************************     7. The DQPSO  Algorithm    *********************/
/*****************************************************************************/
void Initialization()
{
	int i,j,k;
	int p_max; 
	double r;
	for(i=0;i<np;i++)
	  for(j=0;j<N;j++)
	  {
		 QPOP[i][j] = 1.0*(rand()%RAND_MAX)/RAND_MAX; 
	  }
	
	 for(i=0;i<np;i++)
	 {
	 	for(j=0;j<N;j++) 
	 	{
	 	   	r = 1.0*(rand()%RAND_MAX)/RAND_MAX; 
	 	   	if(r > QPOP[i][j])  Pbest[i].X[j] = 1;  
	 	   	else Pbest[i].X[j] = 0; 
		}
	 }
	 
	 for(i=0;i<np;i++) RepairOperator(Pbest[i]);
	 
	 p_max = -99999;
	 for(i=0;i<np;i++)
	 {
	 	if(Pbest[i].f > p_max)
	 	{
	 		for(k=0;k<M;k++) G_BEST.IW[k] = Pbest[i].IW[k]; 
	 		for(j=0;j<N;j++) G_BEST.S[j]  = Pbest[i].S[j]; 
			for(j=0;j<N;j++) G_BEST.NS[j] = Pbest[i].NS[j]; 
			for(j=0;j<N;j++) G_BEST.X[j]  = Pbest[i].X[j]; 
			G_BEST.f  = Pbest[i].f;
			G_BEST.IN = Pbest[i].IN; 
			G_BEST.ON = Pbest[i].ON; 	 
			p_max = Pbest[i].f; 
		}
	 }
}// Initialization of population

int pop_updating(Solution & off_spring, int theta,double T) 
{
	int i,j;
	int f_max = -9999999; 
	int f_min = 9999999;
	int k_max, k_min;
	int k_closest;  
	int count; 
	int similarity = 999999; 
	float random;
	// cout<<T<<endl;
	//cout << T<<endl;
    for(i=0;i<np;i++)
     {
        if(Pbest[i].f > f_max) { k_max = i; f_max = Pbest[i].f; } 
        if(Pbest[i].f < f_min) { k_min = i; f_min = Pbest[i].f; }             
     } 
     
	if((off_spring).f < f_min) return 0;  // don't need to update  
	
	for(i=0;i<np;i++) 
	{
		count = 0;
		for(j=0;j<N;j++) if(Pbest[i].X[j] != (off_spring).X[j]) count ++;  
		if( count < similarity ) { similarity = count; k_closest = i; } 	   
	} // find out the most similar solution 
	random = 1.0*(rand()%RAND_MAX)/RAND_MAX;
	if(off_spring.f > f_max ) 
	{
		//cout << exp((off_spring.f - f_max)/T)<<endl;
	    for(j=0;j<N;j++) Pbest[k_closest].X[j]= off_spring.X[j]; 
		Pbest[k_closest].f =  off_spring.f; 
        return 1;              
    }// deplace the most simlilar solution by the offspring if the offspring is better than the best solution in the population
    random = 1.0*(rand()%RAND_MAX)/RAND_MAX;
   /* if (off_spring.f < Pbest[k_closest].f)
    {   if (random < exp((off_spring.f - Pbest[k_closest].f)/T*1.0))
    {   cout<<"random: "<<random <<endl;
        cout<< exp((off_spring.f - Pbest[k_closest].f)/T*1.0)<<endl;
    }
    
       
    }*/
    if((off_spring.f > Pbest[k_closest].f) && (similarity <= theta) || random < exp((off_spring.f - Pbest[k_closest].f)/T*1.0))  
    {
        for(j=0;j<N;j++) Pbest[k_closest].X[j]= off_spring.X[j]; 
		Pbest[k_closest].f =  off_spring.f;     
        return 1;                                               
    }
	random = 1.0*(rand()%RAND_MAX)/RAND_MAX;
    if( (off_spring.f > f_min) && (similarity > theta) ||random < exp((off_spring.f - f_min)/T*1.0) )
    {
        for(j=0;j<N;j++) Pbest[k_min].X[j]= off_spring.X[j];  
		Pbest[k_min].f =  off_spring.f;  
        return 1;                                                
    }
   return 0;  
} // updating the personal historical best positions of discrete particle swarm (D^{lb})

// The DQPSO Algorithm
void DQPSO( int num)
{
	int i, j, k, m, iter;  
	int lbest, count;  
	double r, e1, e2, e3, a, b,e1_start,e2_start,e1_end,e2_end,sita_start,sita_end,sita_k,sita;  
    int K;   
	int F1, F2; 
	int iterMax;
	double starting_time;  
	double p; 
	starting_time = clock();   
	Initialization(); 
	int answer[30] = {117779,119206,119215,118813,116509,119470,119827,118320,117781,119212,217365,};
	double T_start,T_end;
	double T ;
	// parameters used in the evolution formulas of DQPSO
	a  = 0.000;   
	b  = 1.000;  
	e1 = 0.2;  
	e2 = 0.4;   
	e3 = 0.4;   
	e1_start = 0.3;
	e1_end =0.25;
	e2_start = 0.4;
	e2_end  = 0.25;
	sita_start = 5;
	sita_end = 0;
	sita_k = 1;
	iterMax = 5000;
	T_start = 200;
	T_end = 0;
	
	iter = 0;  // iter denotes the current number of iterations
	while(iter < 5000) // 5000 is the maximum number of iterations which is set according to the size of instances
	{
        e1 = e1_start - (iter*(e1_start - e1_end))/iterMax*1.0;
        e2 = e2_start - (iter*(e2_start - e2_end))/iterMax*1.0;
        //cout << "e1"<<e1<<endl;
	    for(m=0; m<np; m++) 
	    {
           T = (T_start - T_end)* tan(0.875 * pow((1 - 1.0*m/np),2))*( iterMax  - iter )/iterMax + T_end;
           sita = (sita_start - sita_end)* tan(0.875 * pow((1 - 1.0*m/np),2))*( iterMax -iter )/iterMax + sita_end;
		   lbest = -99999;
		   count = 0; 
		   while(1)
		   {
		   	  k = rand()%np;
		   	  if(Pbest[k].f > lbest && k != m) 
				 {
				 	K = k; 
				 	lbest = Pbest[k].f; 
				 }
			  count++; 
			  if(count > 10) break; // 10 represents the number of neighbors of particles
		   }
		   
		   for(j=0;j<N;j++)  PY[j]  =  a*Pbest[m].X[j]  +  b*(1 - Pbest[m].X[j]) ; 
		   for(j=0;j<N;j++)  GY[j]  =  a*Pbest[K].X[j]  +  b*(1 - Pbest[K].X[j]) ;
		   for(j=0;j<N;j++)  QPOP[m][j] = (e1*QPOP[m][j] + e2*PY[j] + (1-e1-e2)*GY[j] ); 
		   
	       for(j=0;j<N;j++)
	    	{
	    		r = 1.0*(rand()%RAND_MAX)/RAND_MAX;  
				if(r > QPOP[m][j]) SC.X[j] = 1;
				else SC.X[j] = 0; 
			} 
			
			RepairOperator(SC); 
			
			p = 1.0*(rand()%RAND_MAX)/RAND_MAX; // p denotes a real-valued random number in [0,1]
			if(p < Prob_LS) VND(SC);
            //cout<<T<<endl;
			pop_updating(SC,sita,T);
			if(SC.f > G_BEST.f)
			{
				for(j=0;j<N;j++) G_BEST.X[j] = SC.X[j];
				G_BEST.f = SC.f; 
				time_one_run = (double) (1.0*(clock()-starting_time)/CLOCKS_PER_SEC); 
				printf("improved objetive value f = %d!\n", G_BEST.f);
                printf("iter= %d!\n", iter );
                cout<<time_one_run<<endl;
			}	
	    }
	    if(iter%200 == 0)
        {
            putchar('\r');
			cout<<"Num:    "<<num<<"  ";
            cout<<"iter:    "<<iter<<endl;
			putchar('\r');
        }
		iter++; 
	}
	
    for(i=0;i<np;i++)
    {
    	RepairOperator(Pbest[i]);
    	VND(Pbest[i]); 
    	if(Pbest[i].f > G_BEST.f)
			{
				for(j=0;j<N;j++) G_BEST.X[j] = Pbest[i].X[j];
				G_BEST.f = Pbest[i].f; 
				time_one_run = (double) (1.0*(clock()-starting_time)/CLOCKS_PER_SEC); 
			}	
	}
	RepairOperator(G_BEST);
}

/*****************************************************************************/
/**** 8. Outputing  results and preprocessing of instances      **************/
/*****************************************************************************/ 
// a. Computing the standard deviation
double Deviation(int arr[], int n)
{
    int i;
    double sum = 0, tmp = 0, x_avg;
    for(i = 0; i < n; ++i) sum += arr[i];
    x_avg = 1.0*sum / n;
    for(i = 0; i < n; ++i)
        tmp += (arr[i] - x_avg)*(arr[i] - x_avg);
    return sqrt(tmp/n); 
} 
//b. Output the best solution found
void OutSol(Solution &S, char *filename)
{
    int i;
    int r;
	FILE *fp; 
	char buff[80];
    sprintf(buff,"%s.sol",filename); 
    fp=fopen(buff,"a+");
    fprintf(fp,"N = %d  M = %d  f= %d\n", N, M, S.f); 
    for(i=0;i<N;i++)
    fprintf(fp,"%d\n", S.X[i]); 
	fclose(fp);
}

//c. Output the statistics information over multipe runs of algorithm
void Outresulting(char *filename, char *outfile)
{
    int i,j;
    FILE *fp; 
   	char buff[80];
    sprintf(buff,"%s",outfile);   
    fp=fopen(buff,"a+");
    fprintf(fp,"%s  %d  %lf  %lf  %lf  %lf  %lf   %d\n",filename,N,BestResult,AvgResult,WorResult,sigmma,AvgTime,Nhit);  
	fclose(fp);         
}

void Outresulting1(char *filename, char *outfile, double F[], double T[], int Number)
{
    int i,j;
    FILE *fp; 
   	char buff[80];
    sprintf(buff,"%s",outfile);   
    fp=fopen(buff,"a+");
    fprintf(fp,"%s  %d  %lf  %lf  %lf  %lf  %lf   %d\n",filename,N,BestResult,AvgResult,WorResult,sigmma,AvgTime,Nhit);  
    for(i=0;i<101;i++)
    {
    	 fprintf(fp,"%d  %lf  %lf \n",i,F[i]/Number, T[i]/Number);  
	}
	fclose(fp);     
}


// d. preprocessing of instances
void ComputeValue()
{
	int i,j;
	for(j=0;j<N;j++) order[j] = j; 
	for(j=0;j<N;j++)
	{ 
	    PR[j] = 0.0; 
		for(i=0;i<M;i++)
		{
			PR[j] += ( (1.0*P[j]*B[i])/(1.0*R[i][j]) ); 
		}
	}
    Quick_Sort(PR, N, order); 
} 

void ComputeValue1()
{
	int i,j;
	for(j=0;j<N;j++) order[j] = j; 
	for(j=0;j<N;j++)
	{ 
	    PR[j] = 0.0; 
		for(i=0;i<M;i++)
		{
			PR[j] +=  (1.0*R[i][j]) / (1.0*B[i]) ; 
		}
		
		PR[j] = 1.0*P[j]/PR[j]; 
	}
	
    Quick_Sort(PR, N, order); 
} 

void Preprocessing()
{
	int i,j; 
	ComputeValue1();
	for(j=0;j<N;j++) P1[j] = P[order[j]];
	for(i=0;i<M;i++)
	{
		for(j=0;j<N;j++) R1[i][j] = R[i][order[j]];
	}
	for(j=0;j<N;j++) P[j] = P1[j]; 
	for(i=0;i<M;i++)
	{
		for(j=0;j<N;j++) R[i][j] = R1[i][j];  
	}   
} 

/*****************************************************************************/
/*****************          9. Main Framework         ************************/
/*****************************************************************************/ 
int main(int argc, char **argv)
{ 
	char *file_name[30] = {"5.250.0.txt","5.250.1.txt", "5.250.2.txt","5.250.3.txt","5.250.4.txt",
                           "5.250.5.txt","5.250.6.txt","5.250.7.txt","5.250.8.txt","5.250.9.txt",
                           "5.250.10.txt","5.250.11.txt","5.250.12.txt","5.250.13.txt","5.250.14.txt",
                           "5.250.15.txt","5.250.16.txt","5.250.17.txt","5.250.18.txt","5.250.19.txt",
                           "5.250.20.txt","5.250.21.txt","5.250.22.txt","5.250.23.txt","5.250.24.txt",
                           "5.250.25.txt","5.250.26.txt","5.250.27.txt","5.250.28.txt","5.250.29.txt",
                        };
	for(int t = 0;t < 30;t++){
    //改进代码  sita 在种群循环中
     int i,j ; 
     int Nruns; // number of runs
     int Result[101];
     srand( time(NULL) ) ; // random seed 
     Nruns = 100;  // number of runs
     Prob_LS = 0.01; 
     File_Name  = file_name[t];
	 cout<<File_Name<<endl;
     //File_Name  = argv[1];
    // outfilename =  argv[2]; 
     outfilename =  "Results.txt"; 
     BestResult = -99999999; 
     WorResult  = 99999999; 
     Nhit = 0; 
     AvgResult = 0.0; 
     AvgTime = 0.0; 
	 

						
     Initializing(); 
     AssignMemery(); 
     Preprocessing();
    
     for(i=0;i<Nruns;i++) // perform the DQPSO algorithm Nruns times
     {
	   DQPSO(i) ; 
       Result[i] =  G_BEST.f; 
       AvgResult += G_BEST.f; 
       AvgTime   += time_one_run;  
       if(G_BEST.f > BestResult)
       {    
           BestResult = G_BEST.f; 
           Nhit = 1;  
           for(j=0;j<N;j++)   SS_BEST.S[j]  = G_BEST.S[j];
	       for(j=0;j<N;j++)   SS_BEST.NS[j] = G_BEST.NS[j]; 
	       for(j=0;j<N;j++)   SS_BEST.X[j]  = G_BEST.X[j] ;
	       for(j=0;j<M;j++)   SS_BEST.IW[j] = G_BEST.IW[j]; 
           SS_BEST.f  = G_BEST.f;
           SS_BEST.IN = G_BEST.IN;
           SS_BEST.ON = G_BEST.ON;  
       }
       else if( abs(G_BEST.f - BestResult) < 1.0e-7)  
       { 
           Nhit++;  
       }
       if(G_BEST.f < WorResult)
       {
           WorResult = G_BEST.f;          
       }
       
     }
     AvgResult /= Nruns; 
     AvgTime   /= Nruns; 
     sigmma = Deviation(Result, Nruns); 
     Outresulting(File_Name,outfilename); 
    // Outresulting1(File_Name,outfilename, FF, TT, Nruns);
     OutSol(SS_BEST, File_Name); 
    // proof(SS_BEST); 
     cout<<File_Name<<endl;
	}
	return 1;
}
