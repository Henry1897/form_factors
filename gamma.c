//
//  gamma.c
//  
//
//  Created by Enrico Fiorenza on 05/11/18.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#define N 12
#define NSPIN 4
#define NCOL 3
int main()
{
    double complex v_0[N][N],v_1[N][N],v_2[N][N],v_3[N][N],v_5[N][N];
    int i,j,k;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            v_0[i][j]=0+0*I;
            v_1[i][j]=0+0*I;
            v_2[i][j]=0+0*I;
            v_3[i][j]=0+0*I;
            v_5[i][j]=0+0*I;
        }
    }
    for (i=0;i<6;i++){
        v_0[i][i+6]=1+0*I;
    }
    for (i=6;i<N;i++){
        v_0[i][i-6]=1+0*I;
    }
    for(i=9;i<12;i++){
        v_1[i][i-9]=0+1*I;
        v_1[i-9][i]=0-1*I;
    }
    for(i=6;i<9;i++){
        v_1[i][i-3]=0+1*I;
        v_1[i-3][i]=0-1*I;
    }
    for(i=9;i<12;i++){
        v_2[i][i-9]=-1+0*I;
        v_2[i-9][i]=-1+0*I;
    }
    for(i=6;i<9;i++){
        v_2[i][i-3]=1+0*I;
        v_2[i-3][i]=1+0*I;
    }
    for (i=6;i<9;i++){
        v_3[i][i-6]=0+1*I;
        v_3[i-6][i]=0-1*I;
    }
    for (i=3;i<6;i++){
        v_3[i][i+6]=0+1*I;
        v_3[i+6][i]=0-1*I;
    }
    for (i=0;i<=5;i++){
        v_5[i][i]=1+0*I;
    }
    for (i=6;i<N;i++){
        v_5[i][i]=-1+0*I;
    }
    printf("la matrice v_0:\n");
    for (k=0;k<N;k++)
    {
        for (j=0;j<N;j++)
        {
            printf("\t%lf", creal(v_0[k][j]));
        }
        printf("\n");
    }
    printf("la matrice v_1:\n");
    for (k=0;k<N;k++)
    {
        for (j=0;j<N;j++)
        {
            printf("\t%lfi", cimag(v_1[k][j]));
        }
        printf("\n");
    }
    printf("la matrice v_2:\n");
    for (k=0;k<N;k++)
    {
        for (j=0;j<N;j++)
        {
            printf("\t%lf", creal(v_2[k][j]));
        }
        printf("\n");
    }
    printf("la matrice v_3:\n");
    for (k=0;k<N;k++)
    {
        for (j=0;j<N;j++)
        {
            printf("\t%lfi", cimag(v_3[k][j]));
        }
        printf("\n");
    }
    printf("la matrice v_5:\n");
    for (k=0;k<N;k++)
    {
        for (j=0;j<N;j++)
        {
            printf("\t%lf", creal(v_5[k][j]));
        }
        printf("\n");
    }
}

