//
//  nobin_test.c
//  
//
//  Created by Enrico Fiorenza on 18/02/2019.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#define NMOM 4508
#define N 12
#define NSPIN 4
#define NCOL 3
#define NCONF 30
#define vol 663552
int isc(int is,int ic);
double complex *LU_decomposition_solver(double complex **S, double complex *b);
double complex **matrix_inverse(double complex **S);
double Z_q_function(double complex **D,double pt0, double pt1, double pt2, double pt3);
int main()
{
    double complex ***Su,***Sd,***GA,***GV,***GP,***GS,***GT;
    FILE *impulse,*prop_u,*prop_d,*out_zq_u,*out_zq_d,*out_V,*out_A,*out_P,*out_S,*out_T;
    char prop_u_name[100],prop_d_name[100],out_zq_uname[100],out_zq_dname[100],out_Vname[100],out_Aname[100];
    double complex **Sum,**Sdm,**Dmu,**Dmd;
    double complex arrayu[144],arrayd[144],GmA[N][N],GmV[N][N],GmP[N][N],GmS[N][N],GmT[N][N],tV,tA,tP,tS,tT;
    double complex v_0[N][N],v_1[N][N],v_2[N][N],v_3[N][N],v_5[N][N],Vec_0[N][N],Id[N][N],s03[N][N];
    int jack,i,j,k,imom,is,ic,is_i,ic_i,m,l;
    double pt_0,pt_1,pt_2,pt_3,pt_q,pt_4,cut,p_0,p_1,p_2,p_3,p_q,imp[NMOM][6];
    double Zqu,Zqd,ZV,ZA,ZP,ZS,ZT;
    int j_0,j_1,j_2,j_3,iconf;
    Su=(complex double***) malloc(sizeof(complex double**)*N);
    Sd=(complex double***) malloc(sizeof(complex double**)*N);
    GA=(complex double***) malloc(sizeof(complex double**)*N);
    GV=(complex double***) malloc(sizeof(complex double**)*N);
    GP=(complex double***) malloc(sizeof(complex double**)*N);
    GS=(complex double***) malloc(sizeof(complex double**)*N);
    GT=(complex double***) malloc(sizeof(complex double**)*N);
    for(i=0;i<N;i++)
    {
        Su[i]=(complex double**) malloc(sizeof(complex double*)*N);
        Sd[i]=(complex double**) malloc(sizeof(complex double*)*N);
        GV[i]=(complex double**) malloc(sizeof(complex double*)*N);
        GA[i]=(complex double**) malloc(sizeof(complex double*)*N);
        GP[i]=(complex double**) malloc(sizeof(complex double*)*N);
        GS[i]=(complex double**) malloc(sizeof(complex double*)*N);
        GT[i]=(complex double**) malloc(sizeof(complex double*)*N);
        for(j=0;j<N;j++)
        {
            Su[i][j]=(complex double*) calloc(NMOM,sizeof(complex double));
            Sd[i][j]=(complex double*) calloc(NMOM,sizeof(complex double));
            GV[i][j]=(complex double*) calloc(NMOM,sizeof(complex double));
            GA[i][j]=(complex double*) calloc(NMOM,sizeof(complex double));
            GP[i][j]=(complex double*) calloc(NMOM,sizeof(complex double));
            GS[i][j]=(complex double*) calloc(NMOM,sizeof(complex double));
            GT[i][j]=(complex double*) calloc(NMOM,sizeof(complex double));
        }
    }
    Sum=(complex double**) malloc(sizeof(complex double*)*N);
    Sdm=(complex double**) malloc(sizeof(complex double*)*N);
    for (i=0;i<N;i++)
    {
        Sum[i]=(complex double*) calloc(N,sizeof(complex double));
        Sdm[i]=(complex double*) calloc(N,sizeof(complex double));
    }
    //matrici gamma
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            Id[i][j]=0+0*I;
            v_0[i][j]=0+0*I;
            v_1[i][j]=0+0*I;
            v_2[i][j]=0+0*I;
            v_3[i][j]=0+0*I;
            v_5[i][j]=0+0*I;
            Vec_0[i][j]=0+0*I;
            s03[i][j]=0+0*I;
        }
    }
    for(i=0;i<3;i++){
        s03[i][i]=1+0*I;
    }
    for(i=3;i<9;i++){
        s03[i][i]=-1+0*I;
    }
    for(i=9;i<N;i++){
        s03[i][i]=1+0*I;
    }
    for(i=0;i<12;i++)
    {
        Id[i][i]=1;
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
    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
            for(l=0;l<N;l++)
            {
                Vec_0[i][j]+=v_0[i][l]*v_5[l][j];
            }
    //impulsi
    impulse=fopen("/Users/enricofiorenza/Desktop/out/mom_list.txt","r");
    if(impulse==NULL)
    {
        printf("error opening impulse file\n");
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
        cut=pt_4/(pow(pt_q,2));
        imp[imom][0]=pt_0;
        imp[imom][1]=pt_1;
        imp[imom][2]=pt_2;
        imp[imom][3]=pt_3;
        imp[imom][4]=pt_q;
        imp[imom][5]=cut;
    }
    fclose(impulse);
    for(iconf=0;iconf<NCONF;iconf++)
    {
        sprintf(prop_u_name,"/Users/enricofiorenza/Desktop/out/0%d/fft_S_M0_R0_0",100+iconf);
        sprintf(prop_d_name,"/Users/enricofiorenza/Desktop/out/0%d/fft_S_M0_R1_0",100+iconf);
        prop_u=fopen(prop_u_name,"r");
        prop_d=fopen(prop_d_name,"r");
        if(prop_u==NULL)
        {
            printf("error opening file fft_S_M0_R0_0 conf=%d\n",100+iconf);
            exit(1);
        }
        if(prop_d==NULL)
        {
            printf("error opening file fft_S_M0_R1_0 conf=%d\n",100+iconf);
            exit(1);
        }
        for(imom=0;imom<NMOM;imom++)
        {
            fread(arrayu,sizeof(complex double),144,prop_u);
            fread(arrayd,sizeof(complex double),144,prop_d);
            for(is=0;is<NSPIN;is++)
                for(ic=0;ic<NCOL;ic++)
                    for(is_i=0;is_i<NSPIN;is_i++)
                        for(ic_i=0;ic_i<NCOL;ic_i++)
                        {
                            k=((isc(is_i,ic_i))+(N*isc(is,ic)));
                            //fread(&arrayu,sizeof(complex double),1,prop_u);
                            //fread(&arrayd,sizeof(complex double),1,prop_d);
                            Sum[isc(is_i,ic_i)][isc(is,ic)]=arrayu[k];
                            Sdm[isc(is_i,ic_i)][isc(is,ic)]=arrayd[k];
                        }
            for(i=0;i<N;i++)
                for(j=0;j<N;j++)
                {
                    GmV[i][j]=0;
                    GmA[i][j]=0;
                    GmP[i][j]=0;
                    GmS[i][j]=0;
                    GmT[i][j]=0;
                }
            for(i=0;i<N;i++)
                for(j=0;j<N;j++)
                    for(k=0;k<N;k++)
                        for(l=0;l<N;l++)
                        {
                            GmV[i][j]+=Sum[i][k]*Vec_0[k][l]*Sum[l][j];
                            GmA[i][j]+=Sum[i][k]*v_0[k][l]*Sum[l][j];
                            GmP[i][j]+=Sum[i][k]*v_5[k][l]*Sum[l][j];
                            GmS[i][j]+=Sum[i][k]*Id[k][l]*Sum[l][j];
                            GmT[i][j]+=Sum[i][k]*s03[k][l]*Sum[l][j];
                            
                        }
            for(i=0;i<N;i++)
                for(j=0;j<N;j++)
                {
                    Su[i][j][imom]+=Sum[i][j]/(NCONF);
                    Sd[i][j][imom]+=Sdm[i][j]/(NCONF);
                    GV[i][j][imom]+=GmV[i][j]/(NCONF);
                    GA[i][j][imom]+=GmA[i][j]/(NCONF);
                    GP[i][j][imom]+=GmP[i][j]/(NCONF);
                    GS[i][j][imom]+=GmS[i][j]/(NCONF);
                    GT[i][j][imom]+=GmT[i][j]/(NCONF);
                }
        }
        fclose(prop_u);
        fclose(prop_d);
    }
    out_A=fopen("/Users/enricofiorenza/Desktop/Vertici_Nf_2/nobinVuu.txt","w");
    out_V=fopen("/Users/enricofiorenza/Desktop/Vertici_Nf_2/nobinAuu.txt","w");
    out_P=fopen("/Users/enricofiorenza/Desktop/Vertici_Nf_2/nobinPuu.txt","w");
    out_S=fopen("/Users/enricofiorenza/Desktop/Vertici_Nf_2/nobinSuu.txt","w");
    out_T=fopen("/Users/enricofiorenza/Desktop/Vertici_Nf_2/nobinTuu.txt","w");
    out_zq_u=fopen("/Users/enricofiorenza/Desktop/Vertici_Nf_2/nobinZqu.txt","w");
    out_zq_d=fopen("/Users/enricofiorenza/Desktop/Vertici_Nf_2/nobinZqd.txt","w");
    for(imom=0;imom<NMOM;imom++)
    {
        for(i=0;i<N;i++)
            for(j=0;j<N;j++)
            {
                Sum[i][j]=Su[i][j][imom];
                Sdm[i][j]=Sd[i][j][imom];
            }
        Dmu=matrix_inverse(Sum);
        Dmd=matrix_inverse(Sdm);
        Zqu=Z_q_function(Dmu,imp[imom][0],imp[imom][1],imp[imom][2],imp[imom][3]);
        Zqd=Z_q_function(Dmd,imp[imom][0],imp[imom][1],imp[imom][2],imp[imom][3]);
        tA=0;
        tV=0;
        tP=0;
        tS=0;
        tT=0;
        for(i=0;i<N;i++)
            for(j=0;j<N;j++)
                for(k=0;k<N;k++)
                    for(l=0;l<N;l++)
                    {
                        tA+=Dmu[i][j]*GA[j][k][imom]*Dmu[k][l]*v_0[l][i];
                        tV+=Dmu[i][j]*GV[j][k][imom]*Dmu[k][l]*Vec_0[l][i];
                        tP+=Dmu[i][j]*GP[j][k][imom]*Dmu[k][l]*v_5[l][i];
                        tS+=Dmu[i][j]*GS[j][k][imom]*Dmu[k][l]*Id[l][i];
                        tT+=Dmu[i][j]*GT[j][k][imom]*Dmu[k][l]*s03[l][i];
                        
                    }
        ZA=creal(tA)/(12.0);
        ZV=creal(tV)/(-12.0);
        ZP=creal(tP)/(12.0);
        ZS=creal(tS)/(12.0);
        ZT=creal(tT)/(12.0);
        /*
        ZA=ZA/sqrt(Zqu*Zqd);
        ZV=ZV/sqrt(Zqu*Zqd);
        ZP=ZP/sqrt(Zqu*Zqd);
        ZA=1/ZA;
        ZP=1/ZP;
        ZV=(-1)/ZV;
         */
        fprintf(out_zq_u,"%.16g %.16g\n",imp[imom][4],Zqu);
        fprintf(out_zq_d,"%.16g %.16g\n",imp[imom][4],Zqd);
        fprintf(out_A,"%.16g %.16g\n",imp[imom][4],ZA);
        fprintf(out_V,"%.16g %.16g\n",imp[imom][4],ZV);
        fprintf(out_P,"%.16g %.16g\n",imp[imom][4],ZP);
        fprintf(out_S,"%.16g %.16g\n",imp[imom][4],ZS);
        fprintf(out_T,"%.16g %.16g\n",imp[imom][4],ZT);
    }
    fclose(out_A);
    fclose(out_V);
    fclose(out_P);
    fclose(out_S);
    fclose(out_T);
    fclose(out_zq_u);
    fclose(out_zq_d);
    free(Su);
    free(Sd);
    free(GA);
    free(GV);
    free(GP);
    free(GS);
    free(GT);
    free(Sum);
    free(Sdm);
}
//funzione per combinare gli indici dei tensori
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
//funzione che calcola Zq  RI'MOM (\Sigma_1)
double Z_q_function(double complex **D,double pt0, double pt1, double pt2, double pt3){
    int i,j,k,counter;
    double Z_q;
    double complex v_0[N][N],v_1[N][N],v_2[N][N],v_3[N][N],v_5[N][N];
    double complex t_1,g_0,g_1,g_2,g_3;
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
    if(pt0!=0)
    {
        t_1=t_1+(g_0/pt0);
        counter=counter+1;
    }
    if(pt1!=0)
    {
        t_1=t_1+(g_1/pt1);
        counter=counter+1;
    }
    if(pt2!=0)
    {
        t_1=t_1+(g_2/pt2);
        counter=counter+1;
    }
    if(pt3!=0)
    {
        t_1=t_1+(g_3/pt3);
        counter=counter+1;
    }
    Z_q=cimag(t_1)/(12.0*(counter)*663552.0);
    return Z_q;
}
