//
//  H4form.c
//  
//
//  Created by Enrico Fiorenza on 15/11/18.
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
#define V 663552
#define NJACK 15
#define tiny 0.000000000000001
double complex *LU_decomposition_solver(double complex **S, double complex *b);
double complex **matrix_inverse(double complex **S);
void swap(double *xp,double *yp);
int main()
{
    bool swapped;
    FILE *impulse,*prop,*out1,*out2,*out3;
    char outname[80];
    double complex **S,**D;
    double r=0,im=0,pt_0,pt_1,pt_2,pt_3,pt_q,pt_4,cutt;
    double sigma_1,sigma_2,sigma_3,sigmam_1,e_1,e_2,e_3;
    double s[NMOM][6+(NJACK*3)],s_1[150][1+NJACK],s_2[150][1+NJACK],s_3[150][1+NJACK];
    double complex t_1,t_2,t_3,g_0,g_1,g_2,g_3;
    double complex v_0[N][N],v_1[N][N],v_2[N][N],v_3[N][N],v_5[N][N],ps[N][N];
    int i,j,k,j_0,j_1,j_2,j_3,imom,counter,jack;
    //inizializzo a zero tutte le variabili
    t_3=0+0*I;
    t_2=0+0*I;
    t_1=0+0*I;
    g_0=0+0*I;
    g_1=0+0*I;
    g_2=0+0*I;
    g_3=0+0*I;
    sigma_1=0;
    e_1=0;
    sigma_2=0;
    e_2=0;
    sigma_3=0;
    e_3=0;
    for(i=0;i<150;i++)
    {
        for(j=0;j<NJACK+1;j++)
        {
            s_1[i][j]=0;
            s_2[i][j]=0;
            s_3[i][j]=0;
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
    impulse=fopen("/Users/enricofiorenza/Desktop/Tesi/fac_form//mom_list.txt","r");
    if(impulse==NULL)
    {
        printf("Error opening mom_list.txt\n");
        exit(1);
    }
    out1=fopen("sigma_1.txt","w");
    if(out1==NULL){
        printf("Error opening out1\n");
        exit(1);
    }
    out2=fopen("sigma_2.txt","w");
    if(out2==NULL){
        printf("Error opening out2\n");
        exit(1);
    }
    out3=fopen("sigma_3.txt","w");
    if(out3==NULL){
        printf("Error opening out3\n");
        exit(1);
    }
    for(imom=0;imom<NMOM;imom++)
    {
        fscanf(impulse,"%d %d %d %d",&j_0,&j_1,&j_2,&j_3);
        pt_0=sin(2.0*M_PI*(j_0+0.5)/48.0);
        pt_1=sin(2.0*M_PI*(j_1)/24.0);
        pt_2=sin(2.0*M_PI*(j_2)/24.0);
        pt_3=sin(2.0*M_PI*(j_3)/24.0);
        pt_q=pow(pt_0,2.0)+pow(pt_1,2.0)+pow(pt_2,2.0)+pow(pt_3,2.0);
        pt_4=pow(pt_0,4.0)+pow(pt_1,4.0)+pow(pt_2,4.0)+pow(pt_3,4.0);
        cutt=pt_4/(pow(pt_q,2));
        s[imom][0]=pt_0;
        s[imom][1]=pt_1;
        s[imom][2]=pt_2;
        s[imom][3]=pt_3;
        s[imom][4]=pt_q;
        s[imom][5]=cutt;
    }
    for(jack=100;jack<115;jack++)
    {
        sprintf(outname,"/Users/enricofiorenza/Desktop/jackknife/H4average//prop%d",jack);
        prop=fopen(outname,"r");
        if(prop==NULL)
        {
            printf("error opening prop j=%d",jack);
            exit(1);
        }
        for(imom=0;imom<NMOM;imom++)
        {
            for(i=0;i<12;i++)
            {
                for(j=0;j<12;j++)
                {
                    fscanf(prop,"%lf %lf",&r,&im);
                    S[i][j]=r+im*I;
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
            //calcolo prima traccia
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
            counter=0;
            if(s[imom][0]!=0)
            {
                t_1=t_1+(g_0/s[imom][0]);
                counter=counter+1;
            }
            if(s[imom][1]!=0)
            {
                t_1=t_1+(g_1/s[imom][1]);
                counter=counter+1;
            }
            if(s[imom][2]!=0)
            {
                t_1=t_1+(g_2/s[imom][2]);
                counter=counter+1;
            }
            if(s[imom][3]!=0)
            {
                t_1=t_1+(g_3/s[imom][3]);
                counter=counter+1;
            }
            //calcolo secondo traccia
            for(i=0;i<N;i++)
            {
                t_2=t_2+D[i][i];
            }
            //calcolo terza traccia
            for(i=0;i<N;i++)
            {
                for(j=0;j<N;j++)
                {
                    t_3=t_3+(v_5[i][j]*D[j][i]);
                }
            }
            s[imom][6+3*(jack-100)]=cimag(t_1)/(12.0*V*counter);
            s[imom][7+3*(jack-100)]=(creal(t_2))/(12.0*V);
            s[imom][8+3*(jack-100)]=(cimag(t_3))/((-12.0)*V);
            g_0=0+0*I;
            g_1=0+0*I;
            g_2=0+0*I;
            g_3=0+0*I;
            t_1=0+0*I;
            t_2=0+0*I;
            t_3=0+0*I;
            
        }
        fclose(prop);
        rewind(impulse);
        printf("fine configurazione %d\n",jack-100);
    }
    //sorting
    for(imom=0;imom<NMOM-1;imom++)
    {
        swapped=false;
        for(j=0;j<NMOM-imom-1;j++)
        {
            if(s[j][4]>s[j+1][4])
            {
                for(k=0;k<(6+3*NJACK);k++)
                {
                    swap(&s[j][k],&s[j+1][k]);
                }
                swapped=true;
            }
        }
        if(swapped==false)
            break;
    }
    printf("sorting done\n");
    //medie H(4)
    for(imom=0;imom<NMOM;imom++)
    {
        if(imom==0)
        {
            for(j=0;j<NJACK;j++){
                s_1[0][j]+=s[imom][6+(3*j)];
                s_2[0][j]+=s[imom][7+(3*j)];
                s_3[0][j]+=s[imom][8+(3*j)];
            }
            k=1;
            counter=0;
        }
        else
        {
            if(s[imom][4]-s[imom-1][4]>tiny)
            {
                for(j=0;j<NJACK;j++){
                     s_1[counter][j]=s_1[counter][j]/k;
                     s_2[counter][j]=s_2[counter][j]/k;
                     s_3[counter][j]=s_3[counter][j]/k;
                }
                k=1;
                s_1[counter][NJACK]=s[imom-1][4];
                s_2[counter][NJACK]=s[imom-1][4];
                s_3[counter][NJACK]=s[imom-1][4];
                counter+=1;
                for(j=0;j<NJACK;j++){
                    s_1[counter][j]+=s[imom][6+(3*j)];
                    s_2[counter][j]+=s[imom][7+(3*j)];
                    s_3[counter][j]+=s[imom][8+(3*j)];
                }
            }
            if(imom==4507 && s[imom][4]-s[imom-1][4]>tiny)
            {
                for(j=0;j<NJACK;j++){
                    s_1[counter][j]=s[imom][6+(3*j)];
                    s_2[counter][j]=s[imom][7+(3*j)];
                    s_3[counter][j]=s[imom][8+(3*j)];
                }
                s_1[counter][NJACK]=s[imom][4];
                s_2[counter][NJACK]=s[imom][4];
                s_3[counter][NJACK]=s[imom][4];
            }
            if(s[imom][4]-s[imom-1][4]<=tiny)
            {
                for(j=0;j<NJACK;j++){
                    s_1[counter][j]+=s[imom][6+(3*j)];
                    s_2[counter][j]+=s[imom][7+(3*j)];
                    s_3[counter][j]+=s[imom][8+(3*j)];
                }
                k+=1;
                if(imom==4507)
                {
                    for(j=0;j<NJACK;j++){
                        s_1[counter][j]=s_1[counter][j]/k;
                        s_2[counter][j]=s_2[counter][j]/k;
                        s_3[counter][j]=s_3[counter][j]/k;
                    }
                    s_1[counter][NJACK]=s[imom][4];
                    s_2[counter][NJACK]=s[imom][4];
                    s_3[counter][NJACK]=s[imom][4];
                }
            }
        }
    }
    printf("fine medie H4\n");
    //stampa risultati
    for(i=0;i<150;i++)
    {
        if(s_1[i][NJACK]!=0)
        {
            for(j=0;j<NJACK;j++){
                sigma_1+=(s_1[i][j]/15.0);
                sigma_2+=(s_2[i][j]/15.0);
                sigma_3+=(s_3[i][j]/15.0);
            }
            for(j=0;j<NJACK;j++){
                e_1+=pow((s_1[i][j]-sigma_1),2.0);
                e_2+=pow((s_2[i][j]-sigma_2),2.0);
                e_3+=pow((s_3[i][j]-sigma_3),2.0);
            }
            e_1=sqrt(((14.0)/15.0)*e_1);
            e_2=sqrt(((14.0)/15.0)*e_2);
            e_3=sqrt(((14.0)/15.0)*e_3);
            fprintf(out1,"%.16g   %.16g   %.16g\n",s_1[i][NJACK],sigma_1,e_1);
            fprintf(out2,"%.16g   %.16g   %.16g\n",s_2[i][NJACK],sigma_2,e_2);
            fprintf(out3,"%.16g   %.16g   %.16g\n",s_1[i][NJACK],sigma_3,e_3);
        }
        e_1=0;
        e_2=0;
        e_3=0;
        sigma_1=0;
        sigma_2=0;
        sigma_3=0;
    }
    fclose(out1);
    fclose(out2);
    fclose(out3);
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

