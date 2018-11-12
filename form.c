//
//  form.c
//  Programma che legge il propagatore medio calcolato da seek3.c e calcola le tracce e i fattori di forma. Output: 4 file in cui sono stampati (ap)^2 ed il corrispondente sigma in out 4 stampo sigma_1 calcolato come suggerito in 1004.1115 formula (33)
//
//  Created by Enrico Fiorenza on 07/11/18.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#define NMOM 4508
#define N 12
#define NSPIN 4
#define NCOL 3
#define V 663552
double complex *LU_decomposition_solver(double complex **S, double complex *b);
double complex **matrix_inverse(double complex **S);
int isc(int is,int ic);
int main()
{
    FILE *impulse,*prop,*out1,*out2,*out3;
    double complex **S,**D;
    double r=0,im=0,s_1,s_2,s_3,p_0,p_1,p_2,p_3,p_q,p,p_4,cut,pt_0,pt_1,pt_2,pt_3,pt_q,pt,pt_4,cutt,sm_1;
    double complex t_1,t_2,t_3,g_0,g_1,g_2,g_3;
    double complex v_0[N][N],v_1[N][N],v_2[N][N],v_3[N][N],v_5[N][N],ps[N][N];
    int i,j,k,j_0,j_1,j_2,j_3,imom,counter;
    //inizializzo a zero tutte le variabili
    t_3=0+0*I;
    t_2=0+0*I;
    t_1=0+0*I;
    g_0=0+0*I;
    g_1=0+0*I;
    g_2=0+0*I;
    g_3=0+0*I;
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
    S=(complex double**) malloc(sizeof(complex double*)*N);
    for (i=0;i<N;i++)
    {
        S[i]=(complex double*) calloc(N,sizeof(complex double));
    }
    D=(complex double** ) malloc(sizeof(complex double* )*N);
    for (i=0;i<N;i++)
    {
        D[i]=(complex double* ) calloc(N,sizeof(complex double));
    }
    ///apro i file necessari
    impulse=fopen("/Users/enricofiorenza/Desktop/Tesi/fac_form//mom_list.txt","r");
    if(impulse==NULL)
    {
        printf("Error opening mom_list.txt\n");
        exit(1);
    }
    out1=fopen("form_1cut.txt","w");
    if(out1==NULL)
    {
        printf("Error opening form_1\n");
        exit(1);
    }
    out2=fopen("form_2cut.txt","w");
    if(out2==NULL)
    {
        printf("Error opening form_2\n");
        exit(1);
    }
    out3=fopen("form_3cut.txt","w");
    if(out3==NULL)
    {
        printf("Error opening form_3\n");
        exit(1);
    }
    out4=fopen("form_1modcut.txt","w");
    if(out4==NULL)
    {
        printf("Error opening form_1mod\n");
        exit(1);
    }
    for(jack=100;jack<115;jack++){
    for(imom=0;imom<NMOM;imom++)
    {
        ///////lettura propagatore
        for(i=0;i<12;i++)
        {
            for(j=0;j<12;j++)
            {
                fscanf(prop,"%lf %lf",&r,&im);
                S[i][j]=r+im*I;
            }
        }
        /*
        printf("propagatore\n");
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                printf("\t%e+%ei",creal(S[i][j]),cimag(S[i][j]));
            }
            printf("\n");
        }
        printf("\n");
        */
        
        //////lettura e calcolo impulso
        fscanf(impulse,"%d %d %d %d",&j_0,&j_1,&j_2,&j_3);
        p_0=j_0*2.0*M_PI/48.0+(M_PI/48.0);
        p_1=j_1*2.0*M_PI/24.0;
        p_2=j_2*2.0*M_PI/24.0;
        p_3=j_3*2.0*M_PI/24.0;
        pt_0=sin(p_0);
        pt_1=sin(p_1);
        pt_2=sin(p_2);
        pt_3=sin(p_3);
        p_q=pow(p_0,2)+pow(p_1,2)+pow(p_2,2)+pow(p_3,2);
        p=sqrt(p_q);
        p_4=pow(p_0,4)+pow(p_1,4)+pow(p_2,4)+pow(p_3,4);
        cut=p_4/(pow(p_q,2));
        pt_q=pow(pt_0,2)+pow(pt_1,2)+pow(pt_2,2)+pow(pt_3,2);
        pt=sqrt(pt_q);
        pt_4=pow(pt_0,4)+pow(pt_1,4)+pow(pt_2,4)+pow(pt_3,4);
        cutt=pt_4/(pow(pt_q,2));
        printf("impulso numero %d p=%f\n cut=%f",imom+1,p,cut);
        if(pt_q<=2.2)
        {
        //////inverto il propagatore
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                D[i][j]=matrix_inverse(S)[i][j];
            }
        }
        /*
        printf("operatore di Dirac");
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                printf("\t%e+%ei", creal(D[i][j]),cimag(D[i][j]));
            }
            printf("\n");
        }
        */
        printf("\n");
        //////calcolo pslash
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                ps[i][j]=(p_0*v_0[i][j])+(p_1*v_1[i][j])+(p_2*v_2[i][j])+(p_3*v_3[i][j]);
            }
        }
        /*
        printf("pslash\n");
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                printf("\t%e", creal(ps[i][j]));
            }
            printf("\n");
        }
        */
        printf("\n");
        ////calcolo le tracce
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                t_1=t_1+(ps[i][j]*D[j][i]);
            }
        }
        for(i=0;i<N;i++)
        {
            t_2=t_2+D[i][i];
        }
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                t_3=t_3+(v_5[i][j]*D[j][i]);
            }
        }
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                g_0=g_0+(v_0[i][j]*D[j][i]);
                g_1=g_1+(v_1[i][j]*D[j][i]);
                g_2=g_2+(v_2[i][j]*D[j][i]);
                g_3=g_3+(v_3[i][j]*D[j][i]);
            }
        }
        ///////calcolo fattori di forma
        if(p_q!=0)
        {
            
            s_1=(cimag(t_1))/(12.0*p_q*V);
            s_2=(creal(t_2))/(12.0*V);
            s_3=(cimag(t_3))/((-12.0)*V);
            fprintf(out1,"%f %f\n",pt_q,s_1);
            fprintf(out2,"%f %f\n",pt_q,s_2);
            fprintf(out3,"%f %f\n",pt_q,s_3);
            ///metodo alternativo s_1 pag 17 articolo 1004.1115
            t_1=0+0*I;
            counter=0;
            if(p_0!=0)
            {
                t_1=t_1+(g_0/pt_0);
                counter=counter+1;
            }
            if(p_1!=0)
            {
                t_1=t_1+(g_1/pt_1);
                counter=counter+1;
            }
            if(p_2!=0)
            {
                t_1=t_1+(g_2/pt_2);
                counter=counter+1;
            }
            if(p_3!=0)
            {
                t_1=t_1+(g_3/pt_3);
                counter=counter+1;
            }
            sm_1=cimag(t_1)/(12*V*counter);
            fprintf(out4,"%f %f\n",pt_q,sm_1);
        }
        g_0=0+0*I;
        g_1=0+0*I;
        g_2=0+0*I;
        g_3=0+0*I;
        t_3=0+0*I;
        t_2=0+0*I;
        t_1=0+0*I;
        }
////////////fine loop impulsi
    }
///////////fine loop jackknife
    }
    fclose(impulse);
    fclose(prop);
    fclose(out1);
    fclose(out2);
    fclose(out3);
    fclose(out4);
    return(0);
    
    
}
//////funzioni usate
int isc(int is,int ic)
{
    return(ic+NCOL*is);
}
//LU matrix decomposition (Numerical Recipes, Marco)
double complex *LU_decomposition_solver(double complex **S, double complex *b){
    double complex **U,**L,*y,*x;
    int i,j,k;
    
    x=(complex double*) malloc(sizeof(complex double)*N);
    y=(complex double*) malloc(sizeof(complex double)*N);
    L=(complex double**) malloc(sizeof(complex double*)*N);
    U=(complex double**) malloc(sizeof(complex double*)*N);
    for (i=0;i<N;i++){
        U[i]=(complex double*) calloc(N,sizeof(complex double));
        L[i]=(complex double*) calloc(N,sizeof(complex double));
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
double complex **matrix_inverse(double complex **S ){
    double complex *b,**r,*a;
    int i,j;
    
    b=(complex double*) calloc(N,sizeof(complex double));
    r=(complex double**) malloc(sizeof(complex double*)*N);
    for (i=0;i<N;i++){
        r[i]=(complex double*) malloc(N*sizeof(complex double));
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





