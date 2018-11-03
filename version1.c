//
//  version1.c
//  
//
//  Created by Enrico Fiorenza
//

#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#define N 12
#define NSPIN 4
#define NCOL 3
double **matrix();
double *LU_decomposition_solver(double **S, double *b);
double **matrix_inverse(double **S  );
int main(void)
{
    FILE *prop;
    double **S,**D;
    double v_5[N][N],v_4[N][N],t_3,t_2,t_1;
    int i,j,k;
    //inizializzo a zero tutte le variabili e matrici
    t_3=0;
    t_2=0;
    t_1=0;
    S=(double**) malloc(sizeof(double*)*N);
    for (i=0;i<N;i++)
    {
        S[i]=(double*) calloc(N,sizeof(double));
    }
    D=(double**) malloc(sizeof(double*)*N);
    for (i=0;i<N;i++)
    {
        D[i]=(double*) calloc(N,sizeof(double));
    }
    for (k=0;k<N;k++)
    {
        for (j=0;j<N;j++)
        {
            v_5[k][j]=0;
            v_4[k][j]=0;
        }
    }
    //definisco e stampo a schermo un propagatore fittizio 2*id
    for (j=0;j<N;j++){
        S[j][j]=2;
    }
    for (k=0;k<N;k++)
    {
        for (j=0;j<N;j++)
        {
            S[k][j]=matrix()[k][j];
        }
        //printf("\n");
    }
    printf("Il propagatore S:\n");
    for (k=0;k<N;k++)
    {
        for (j=0;j<N;j++)
        {
            printf("\t%lf", S[k][j]);
        }
        printf("\n");
    }
    //definisco e stampo gli elementi della matrice v_5 definita come (gamma5 prod.tens id(3x3))
        printf("\n");
    for (i=0;i<=5;i++){
        v_5[i][i]=1;
    }
    for (i=6;i<N;i++){
        v_5[i][i]=-1;
    }
    printf("la matrice v_5:\n");
    for (k=0;k<N;k++)
    {
        for (j=0;j<N;j++)
        {
            printf("\t%lf", v_5[k][j]);
        }
    printf("\n");
    }
    printf("\n");
    //definisco e stampo gli elementi della matrice v_4 definita come (gamma4 prod.tens id(3x3))
    for (i=0;i<6;i++){
        v_4[i][i+6]=1;
    }
    for (i=6;i<N;i++){
        v_4[i][i-6]=1;
    }
    printf("la matrice v_4:\n");
    for (k=0;k<N;k++)
    {
        for (j=0;j<N;j++)
        {
            printf("\t%lf", v_4[k][j]);
        }
        printf("\n");
    }
    /* Reading part, that gives the segmentation fault
    prop=fopen("matrix.dat","r");
    for (k=0;k<N;k++)
    {
        for (j=0;j<N;j++)
        {
            fscanf(prop,"%lf",&S[i][j]);
        }
        
    }
    */
    //inverto la matrice S e stampo a schermo il risultato per controllo
    printf("\n");
    for (k=0;k<N;k++)
    {
        for (j=0;j<N;j++)
        {
            D[k][j]=matrix_inverse(S)[k][j];
        }
        //printf("\n");
    }
    printf("l'operatore di Dirac D:\n");
    for (k=0;k<N;k++)
    {
        for (j=0;j<N;j++)
        {
            printf("\t%lf", D[k][j]);
        }
        printf("\n");
    }
    //calcolo e stampo a schermo la traccia relativa al fattore di forma 3
    for (i=0;i<N;i++){
        for(j=0;j<N;j++){
            t_3=t_3+(v_5[i][j]*D[j][i]);
        }
        //printf("traccia_3 è %lf \n",t_3);
    }
    printf("traccia_3 è %lf \n",t_3);
    //calcolo e stampo a schermo la traccia relativa al fattore di forma 2
    for (i=0;i<N;i++){
        t_2=t_2+D[i][i];
    }
    printf("traccia_2 è %lf \n",t_2);
    //calcolo e stampo a schermo la traccia relativa al fattore di forma 1
    for (i=0;i<N;i++){
        for(j=0;j<N;j++){
            t_1=t_1+(v_4[i][j]*D[j][i]);
        }
        //printf("traccia_1 è %lf \n",t_1);
    }
    printf("traccia_1 è %lf \n",t_1);
    
    
    
    
    
    
    
}



//LU matrix decomposition (Numerical Recipes, Marco)
double *LU_decomposition_solver(double **S, double *b){
    double **U,**L,*y,*x;
    int i,j,k;
    
    x=(double*) malloc(sizeof(double)*N);
    y=(double*) malloc(sizeof(double)*N);
    L=(double**) malloc(sizeof(double*)*N);
    U=(double**) malloc(sizeof(double*)*N);
    for (i=0;i<N;i++){
        U[i]=(double*) calloc(N,sizeof(double));
        L[i]=(double*) calloc(N,sizeof(double));
    }
    
    for (i=0;i<N;i++)
        L[i][i]=1;
    
    for (j=0;j<N;j++){
        for (i=0;i<=j;i++){
            U[i][j]=S[i][j];
            for (k=0;k<i;k++)
                U[i][j]-=L[i][k]*U[k][j];
        }
        for (i=j+1;i<N;i++){
            L[i][j]=S[i][j];
            for (k=0;k<j;k++)
                L[i][j]-=L[i][k]*U[k][j];
            L[i][j]/=U[j][j];
        }
    }
    
    y[0]=b[0]/L[0][0];
    for (i=0;i<N;i++){
        y[i]=b[i];
        for (k=0;k<i;k++)
            y[i]-=L[i][k]*y[k];
        y[i]/=L[i][i];
    }
    
    x[N-1]=y[N-1]/U[N-1][N-1];
    for (i=N-2;i>=0;i--){
        x[i]=y[i];
        for (k=i+1;k<N;k++)
            x[i]-=U[i][k]*x[k];
        x[i]/=U[i][i];
    }
    
    free(y);
    for (i=0;i<N;i++){
        free(L[i]);
        free(U[i]);
    }
    free(U);
    free(L);
    
    return x;
    
}

//return the inverse matrix of M
double **matrix_inverse(double **S ){
    double *b,**r,*a;
    int i,j;
    
    b=(double*) calloc(N,sizeof(double));
    r=(double**) malloc(sizeof(double*)*N);
    for (i=0;i<N;i++){
        r[i]=(double*) malloc(N*sizeof(double));
    }
    
    for (i=0;i<N;i++){
        b[i]=1.;
        a=LU_decomposition_solver(S, b);
        for (j=0;j<N;j++)
            r[j][i]=a[j];
        free(a);
        b[i]=0;
    }
    
    free(b);
    
    return r;
}

double **matrix()
{
    FILE *prop;
    int i,j,k,is,ic,is_i,ic_i;
    double **S;
    S=(double**) malloc(sizeof(double*)*N);
    for(i=0;i<N;i++)
    {
        S[i]=(double*) calloc(N,sizeof(double));
    }
    printf("il propagatore prima della lettura\n");
    for(j=0;j<N;j++)
    {
        for(k=0;k<N;k++){
            printf("%f\t",S[j][k]);
        }
        printf("\n");
    }
    
    prop=fopen("fft_S_M0_R0_0","r");
    if(prop==NULL)
    {
        printf("Error opening fft_S_M0_R0_0\n");
        exit(1);
    }
    i=0;
    /*for(is=0;is<NSPIN;is++)
     for(ic=0;ic<NCOL;ic++)
     for(is_i=0;is_i<NSPIN;is_i++)
     for(ic_i=0;ic_i<NCOL;ic_i++)
     {
     printf("hello %d\n",k);
     k=k+1;
     i=ic_i*is_i;
     j=ic*is;
     S[i][j]=fgetc(prop);
     //fscanf(prop,"%lf",&S[i][j]);
     }
     */
    for(j=0;j<N;j++)
    {
        for(k=0;k<N;k++){
            printf("hello %d\n",i);
            i=i+1;
            S[j][k]=fgetc(prop);
        }
        printf("\n");
    }
    fclose(prop);
    /*printf("il propagatore dopo la lettura\n");
    for(j=0;j<N;j++)
    {
        for(k=0;k<N;k++){
            printf("%lf\t",S[j][k]);
        }
        printf("\n");
    }
    */
    return S;
}


