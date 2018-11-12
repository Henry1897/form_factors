//
//  jackform.c
//  
//
//  Created by Enrico Fiorenza on 09/11/18.
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
#define NJACK 15
double complex *LU_decomposition_solver(double complex **S, double complex *b);
double complex **matrix_inverse(double complex **S);
int isc(int is,int ic);
int main ()
{
    FILE *impulse,*prop,*out1,*out2,*out3,*out4;
    char outname[80];
    double complex **S,**D;
    double r=0,im=0,p_0,p_1,p_2,p_3,p_q,p,p_4,cut,pt_0,pt_1,pt_2,pt_3,pt_q,pt,pt_4,cutt;
    double sigma_1,sigma_2,sigma_3,sigmam_1,e_1,e_2,e_3,em_1;
    double s_1[NMOM][NJACK],s_2[NMOM][NJACK],s_3[NMOM][NJACK],s_4[NMOM][NJACK],sm_1[NMOM][NJACK];
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
    for(jack=100;jack<115;jack++)
    {
        sprintf(outname,"/Users/enricofiorenza/Desktop/jackknife//prop%d",jack);
        prop=fopen(outname,"r");
        if(prop==NULL)
        {
            printf("error opening prop");
            exit(1);
        }
        for(imom=0;imom<NMOM;imom++)
        {
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
            printf("impulso numero=%d jack=%d\n",imom+1,jack-99);
            for(i=0;i<12;i++)
            {
                for(j=0;j<12;j++)
                {
                    fscanf(prop,"%lf %lf",&r,&im);
                    S[i][j]=r+im*I;
                }
            }
            if(pt_q<=2.2 && cutt<0.28)
            {
                for(i=0;i<N;i++)
                {
                    for(j=0;j<N;j++)
                    {
                        D[i][j]=matrix_inverse(S)[i][j];
                    }
                }
                for(i=0;i<N;i++)
                {
                    for(j=0;j<N;j++)
                    {
                        ps[i][j]=(p_0*v_0[i][j])+(p_1*v_1[i][j])+(p_2*v_2[i][j])+(p_3*v_3[i][j]);
                    }
                }
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
                s_1[imom][jack-100]=(cimag(t_1))/(12.0*p_q*V);
                s_2[imom][jack-100]=(creal(t_2))/(12.0*V);
                s_3[imom][jack-100]=(cimag(t_3))/((-12.0)*V);
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
                sm_1[imom][jack-100]=cimag(t_1)/(12*V*counter);
                g_0=0+0*I;
                g_1=0+0*I;
                g_2=0+0*I;
                g_3=0+0*I;
                t_1=0+0*I;
                t_2=0+0*I;
                t_3=0+0*I;
                ///fine if democratic cut
            }
            /////fine loop impulsi
        }
        fclose(prop);
        rewind(impulse);
        /////fine loop jack
    }
    out1=fopen("jackform1.txt","w");
    if(out1==NULL)
    {
        printf("error opening out1");
        exit(1);
    }
    out2=fopen("jackform2.txt","w");
    if(out2==NULL)
    {
        printf("error opening out2");
        exit(1);
    }
    out3=fopen("jackform3.txt","w");
    if(out3==NULL)
    {
        printf("error opening out3");
        exit(1);
    }
    out4=fopen("jackform1m.txt","w");
    if(out4==NULL)
    {
        printf("error opening out4");
        exit(1);
    }
    sigma_1=0;
    sigma_2=0;
    sigma_3=0;
    sigmam_1=0;
    e_1=0;
    e_2=0;
    e_3=0;
    em_1=0;
    rewind(impulse);
    for(imom=0;imom<NMOM;imom++)
    {
        fscanf(impulse,"%d %d %d %d",&j_0,&j_1,&j_2,&j_3);
        p_0=j_0*2.0*M_PI/48.0+(M_PI/48.0);
        p_1=j_1*2.0*M_PI/24.0;
        p_2=j_2*2.0*M_PI/24.0;
        p_3=j_3*2.0*M_PI/24.0;
        pt_0=sin(p_0);
        pt_1=sin(p_1);
        pt_2=sin(p_2);
        pt_3=sin(p_3);
        pt_q=pow(pt_0,2)+pow(pt_1,2)+pow(pt_2,2)+pow(pt_3,2);
        pt=sqrt(pt_q);
        pt_4=pow(pt_0,4)+pow(pt_1,4)+pow(pt_2,4)+pow(pt_3,4);
        cutt=pt_4/(pow(pt_q,2));
        if(pt_q<=2.2 && cutt<0.28)
        {
            for(jack=0;jack<NJACK;jack++)
            {
                sigma_1=sigma_1+s_1[imom][jack];
                sigma_2=sigma_2+s_2[imom][jack];
                sigma_3=sigma_3+s_3[imom][jack];
                sigmam_1=sigmam_1+sm_1[imom][jack];
            }
            sigma_1=sigma_1/NJACK;
            sigma_2=sigma_2/NJACK;
            sigma_3=sigma_3/NJACK;
            sigmam_1=sigmam_1/NJACK;
            for(jack=0;jack<NJACK;jack++)
            {
                e_1=e_1+pow((s_1[imom][jack]-sigma_1),2);
                e_2=e_2+pow((s_2[imom][jack]-sigma_2),2);
                e_3=e_3+pow((s_3[imom][jack]-sigma_3),2);
                em_1=em_1+pow((sm_1[imom][jack]-sigmam_1),2);
            }
            e_1=sqrt((14.0/15.0)*e_1);
            e_2=sqrt((14.0/15.0)*e_2);
            e_3=sqrt((14.0/15.0)*e_3);
            em_1=sqrt((14.0/15.0)*em_1);
            fprintf(out1,"%.15g    %.15g    %.15g\n",pt_q,sigma_1,e_1);
            fprintf(out2,"%.15g    %.15g    %.15g\n",pt_q,sigma_2,e_2);
            fprintf(out3,"%.15g    %.15g    %.15g\n",pt_q,sigma_3,e_3);
            fprintf(out4,"%.15g    %.15g    %.15g\n",pt_q,sigmam_1,em_1);
        }
        sigma_1=0;
        sigma_2=0;
        sigma_3=0;
        sigmam_1=0;
        e_1=0;
        e_2=0;
        e_3=0;
        em_1=0;
    }
    fclose(impulse);
    fclose(out1);
    fclose(out2);
    fclose(out3);
    fclose(out4);
    free(S);
    free(D);
    
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

