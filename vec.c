//
//  vec.c
//  
//
//  Created by Enrico Fiorenza on 25/11/18.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#define NMOM 4508
#define N 12
#define NSPIN 4
#define NCOL 3
#define Vol 663552
#define NJACK 15
#define tiny 0.000000000000001
double complex *LU_decomposition_solver(double complex **S, double complex *b);
double complex **matrix_inverse(double complex **S);
void swap(double *xp,double *yp);
int main()
{
    bool swapped;
    FILE *impulse,*prop,*out,*vecT,*vecX,*vecY,*vecZ;
    char name1[80],name2[80];
    double complex **S,**D,VT[N][N],VX[N][N],VY[N][N],VZ[N][N];
    double r,im,pt_0,pt_1,pt_2,pt_3,pt_q,pt_4,cutt,p_0,p_1,p_2,p_3,p_q;
    double eT,eX,eY,eZ,gammaT,gammaX,gammaY,gammaZ,e,gamma;
    double imp[NMOM][6],impm[150];
    double complex t,x,y,z;
    double complex v_0[N][N],v_1[N][N],v_2[N][N],v_3[N][N],v_5[N][N];
    double res[NMOM][NJACK*4],s_T[150][NJACK],s_X[150][NJACK],s_Y[150][NJACK],s_Z[150][NJACK];
    int i,j,k,j_0,j_1,j_2,j_3,imom,jack,counter,l;
    //inizializzo a zero tutte le variabili
    t=0+0*I;
    gammaT=0;
    eT=0;
    gammaX=0;
    eX=0;
    gammaY=0;
    eY=0;
    gammaZ=0;
    eZ=0;
    for(i=0;i<150;i++)
    {
        for(j=0;j<NJACK;j++)
        {
            s_T[i][j]=0;
            s_X[i][j]=0;
            s_Y[i][j]=0;
            s_Z[i][j]=0;
        }
    }
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
        v_0[i][i+6]=-1+0*I;
    }
    for (i=6;i<N;i++){
        v_0[i][i-6]=-1+0*I;
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
    impulse=fopen("/Users/enricofiorenza/Desktop/Tesi/fac_form//mom_list.txt","r");
    if(impulse==NULL){
        printf("error opening mom_list.txt\n");
        exit(1);
    }
    for(imom=0;imom<NMOM;imom++)
    {
        fscanf(impulse,"%d %d %d %d",&j_0,&j_1,&j_2,&j_3);
        p_0=2.0*M_PI*(j_0+0.5)/48.0;
        p_1=2.0*M_PI*(j_1)/24.0;
        p_2=2.0*M_PI*(j_2)/24.0;
        p_3=2.0*M_PI*(j_3)/24.0;
        p_q=pow(p_0,2.0)+pow(p_1,2.0)+pow(p_2,2.0)+pow(p_3,2.0);
        pt_0=sin(2.0*M_PI*(j_0+0.5)/48.0);
        pt_1=sin(2.0*M_PI*(j_1)/24.0);
        pt_2=sin(2.0*M_PI*(j_2)/24.0);
        pt_3=sin(2.0*M_PI*(j_3)/24.0);
        pt_q=pow(pt_0,2.0)+pow(pt_1,2.0)+pow(pt_2,2.0)+pow(pt_3,2.0);
        pt_4=pow(pt_0,4.0)+pow(pt_1,4.0)+pow(pt_2,4.0)+pow(pt_3,4.0);
        cutt=pt_4/(pow(pt_q,2));
        imp[imom][0]=pt_0;
        imp[imom][1]=pt_1;
        imp[imom][2]=pt_2;
        imp[imom][3]=pt_3;
        imp[imom][4]=pt_q;
        imp[imom][5]=cutt;
    }
    for(jack=100;jack<100+NJACK;jack++)
    {
        sprintf(name1,"/Users/enricofiorenza/Desktop/vector/15A/file_jack//prop%d",jack);
        prop=fopen(name1,"r");
        if(prop==NULL)
        {
            printf("error opening prop j=%d\n",jack);
            exit(1);
        }
        sprintf(name2,"/Users/enricofiorenza/Desktop/vector/15A/file_jack//vecT%d",jack);
        vecT=fopen(name2,"r");
        if(vecT==NULL)
        {
            printf("error opening vecT j=%d\n",jack);
            exit(1);
        }
        sprintf(name2,"/Users/enricofiorenza/Desktop/vector/15A/file_jack//vecX%d",jack);
        vecX=fopen(name2,"r");
        if(vecX==NULL)
        {
            printf("error opening vecX j=%d\n",jack);
            exit(1);
        }
        sprintf(name2,"/Users/enricofiorenza/Desktop/vector/15A/file_jack//vecY%d",jack);
        vecY=fopen(name2,"r");
        if(vecY==NULL)
        {
            printf("error opening vecY j=%d\n",jack);
            exit(1);
        }
        sprintf(name2,"/Users/enricofiorenza/Desktop/vector/15A/file_jack//vecZ%d",jack);
        vecZ=fopen(name2,"r");
        if(vecZ==NULL)
        {
            printf("error opening vecZ j=%d\n",jack);
            exit(1);
        }
        for(imom=0;imom<NMOM;imom++)
        {
            for(i=0;i<N;i++)
            {
                for(j=0;j<N;j++)
                {
                    fscanf(prop,"%lf %lf",&r,&im);
                    S[i][j]=r+im*I;
                    fscanf(vecT,"%lf %lf",&r,&im);
                    VT[i][j]=r+im*I;
                    fscanf(vecX,"%lf %lf",&r,&im);
                    VX[i][j]=r+im*I;
                    fscanf(vecY,"%lf %lf",&r,&im);
                    VY[i][j]=r+im*I;
                    fscanf(vecZ,"%lf %lf",&r,&im);
                    VZ[i][j]=r+im*I;
                }
            }
            for(i=0;i<12;i++)
            {
                for(j=0;j<12;j++)
                {
                    D[i][j]=matrix_inverse(S)[i][j];
                    //D[i][j]=S[i][j];
                }
            }
            t=0+0*I;
            x=0+0*I;
            y=0+0*I;
            z=0+0*I;
            for(i=0;i<N;i++)
                for(j=0;j<N;j++)
                    for(k=0;k<N;k++)
                        for(l=0;l<N;l++)
                        {
                            t+=v_0[i][j]*D[j][k]*VT[k][l]*D[l][i];
                            x+=v_1[i][j]*D[j][k]*VX[k][l]*D[l][i];
                            y+=v_2[i][j]*D[j][k]*VY[k][l]*D[l][i];
                            z+=v_3[i][j]*D[j][k]*VZ[k][l]*D[l][i];
                        }
            res[imom][jack-100]=creal(t/(Vol*12.0));
            res[imom][jack-100+NJACK]=creal(x/(Vol*12.0));
            res[imom][jack-100+(2*NJACK)]=creal(y/(Vol*12.0));
            res[imom][jack-100+(3*NJACK)]=creal(z/(Vol*12.0));
        }
        fclose(prop);
        fclose(vecT);
        fclose(vecX);
        fclose(vecY);
        fclose(vecZ);
        printf("fine configurazione %d\n",jack-100);
    }
    //sorting
    for(imom=0;imom<NMOM-1;imom++)
    {
        swapped=false;
        for(j=0;j<NMOM-imom-1;j++)
        {
            if(imp[j][4]>imp[j+1][4])
            {
                for(k=0;k<(4*NJACK);k++)
                {
                    swap(&res[j][k],&res[j+1][k]);
                }
                for(k=0;k<6;k++)
                {
                    swap(&imp[j][k],&imp[j+1][k]);
                }
                swapped=true;
            }
        }
        if(swapped==false)
            break;
    }
    printf("sorting done \n");
    //media H(4)
    for(imom=0;imom<NMOM;imom++)
    {
        if(imom==0)
        {
            for(j=0;j<NJACK;j++){
                s_T[0][j]+=res[imom][j];
                s_X[0][j]+=res[imom][j+NJACK];
                s_Y[0][j]+=res[imom][j+(2*NJACK)];
                s_Z[0][j]+=res[imom][j+(3*NJACK)];
            }
            k=1;
            counter=0;
        }
        else
        {
            if(imp[imom][4]-imp[imom-1][4]>tiny)
            {
                for(j=0;j<NJACK;j++){
                    s_T[counter][j]=s_T[counter][j]/k;
                    s_X[counter][j]=s_X[counter][j]/k;
                    s_Y[counter][j]=s_Y[counter][j]/k;
                    s_Z[counter][j]=s_Z[counter][j]/k;
                }
                k=1;
                impm[counter]=imp[imom-1][4];
                counter+=1;
                for(j=0;j<NJACK;j++){
                    s_T[counter][j]+=res[imom][j];
                    s_X[counter][j]+=res[imom][j+NJACK];
                    s_Y[counter][j]+=res[imom][j+(2*NJACK)];
                    s_Z[counter][j]+=res[imom][j+(3*NJACK)];
                }
            }
            if(imom==4507 && imp[imom][4]-imp[imom-1][4]>tiny)
            {
                for(j=0;j<NJACK;j++){
                    s_T[counter][j]=res[imom][j];
                    s_X[counter][j]=res[imom][j+NJACK];
                    s_Y[counter][j]=res[imom][j+(2*NJACK)];
                    s_Z[counter][j]=res[imom][j+(3*NJACK)];
                }
                impm[counter]=imp[imom][4];
            }
            if(imp[imom][4]-imp[imom-1][4]<=tiny)
            {
                for(j=0;j<NJACK;j++){
                    s_T[counter][j]+=res[imom][j];
                    s_X[counter][j]+=res[imom][j+NJACK];
                    s_Y[counter][j]+=res[imom][j+(2*NJACK)];
                    s_Z[counter][j]+=res[imom][j+(3*NJACK)];
                }
                k+=1;
                if(imom==4507)
                {
                    for(j=0;j<NJACK;j++){
                        s_T[counter][j]=s_T[counter][j]/k;
                        s_X[counter][j]=s_X[counter][j]/k;
                        s_Y[counter][j]=s_Y[counter][j]/k;
                        s_Z[counter][j]=s_Z[counter][j]/k;
                    }
                    impm[counter]=imp[imom][4];
                }
            }
        }
    }
    out=fopen("gamma.txt","w");
    printf("fine medie H4\n");
    //stampa risultati
    for(i=0;i<150;i++)
    {
        if(impm[i]!=0)
        {
            for(j=0;j<NJACK;j++)
            {
                gammaT+=(s_T[i][j]/(NJACK));
                gammaX+=(s_X[i][j]/(NJACK));
                gammaY+=(s_Y[i][j]/(NJACK));
                gammaZ+=(s_Z[i][j]/(NJACK));
            }
            for(j=0;j<NJACK;j++){
                eT+=pow((s_T[i][j]-gammaT),2.0);
                eX+=pow((s_X[i][j]-gammaX),2.0);
                eY+=pow((s_Y[i][j]-gammaY),2.0);
                eZ+=pow((s_Z[i][j]-gammaZ),2.0);
            }
            eT=sqrt(((NJACK-1.0)/(NJACK))*eT);
            eX=sqrt(((NJACK-1.0)/(NJACK))*eX);
            eY=sqrt(((NJACK-1.0)/(NJACK))*eY);
            eZ=sqrt(((NJACK-1.0)/(NJACK))*eZ);
            gamma=(gammaT+gammaX+gammaY+gammaZ)/(4.0);
            e=(eT+eX+eY+eZ)/(4.0);
            fprintf(out,"%.16g   %.16g   %.16g\n",impm[i],gamma,e);
        }
        eT=0;
        eX=0;
        eY=0;
        eZ=0;
        gammaT=0;
        gammaX=0;
        gammaY=0;
        gammaZ=0;
    }
    fclose(out);
    fclose(impulse);
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

void swap(double *xp,double *yp)
{
    double temp=*xp;
    *xp=*yp;
    *yp=temp;
}
