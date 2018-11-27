//
//  jackseek.c
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
int isc(int is,int ic);
int main()
{
    FILE *propbin,*out,*vecbinT,*vecbinX,*vecbinY,*vecbinZ,*outvecT,*outvecX,*outvecY,*outvecZ;
    char filename[80],outname[80];
    double complex S[12][12],Sm[12][12],T[12][12],Tm[12][12],X[12][12],Xm[12][12],Y[12][12],Ym[12][12],Z[12][12],Zm[12][12];
    double complex *array;
    int iconf,i,j,k,p_0,p_1,p_2,p_3,imom,is,ic,is_i,ic_i,jack;
    array=(complex double*) calloc(144,sizeof(complex double));
    for(i=0;i<12;i++)
    {
        for(j=0;j<12;j++)
        {
            S[i][j]=0+0*I;
            Sm[i][j]=0+0*I;
            T[i][j]=0+0*I;
            Tm[i][j]=0+0*I;
            X[i][j]=0+0*I;
            Xm[i][j]=0+0*I;
            Y[i][j]=0+0*I;
            Ym[i][j]=0+0*I;
            Z[i][j]=0+0*I;
            Zm[i][j]=0+0*I;
        }
    }
    for(jack=100;jack<115;jack++)
    {
        sprintf(outname,"/Users/enricofiorenza/Desktop/vector/15A/file_jack//prop%d",jack);
        out=fopen(outname,"w");
        if(out==NULL)
        {
            printf("Error opening output\n");
            exit(1);
        }
        sprintf(outname,"/Users/enricofiorenza/Desktop/vector/15A/file_jack//vecT%d",jack);
        outvecT=fopen(outname,"w");
        if(outvecT==NULL)
        {
            printf("Error opening outputvecT\n");
            exit(1);
        }
        sprintf(outname,"/Users/enricofiorenza/Desktop/vector/15A/file_jack//vecX%d",jack);
        outvecX=fopen(outname,"w");
        if(outvecX==NULL)
        {
            printf("Error opening outputvecX\n");
            exit(1);
        }
        sprintf(outname,"/Users/enricofiorenza/Desktop/vector/15A/file_jack//vecY%d",jack);
        outvecY=fopen(outname,"w");
        if(outvecY==NULL)
        {
            printf("Error opening outputvecY\n");
            exit(1);
        }
        sprintf(outname,"/Users/enricofiorenza/Desktop/vector/15A/file_jack//vecZ%d",jack);
        outvecZ=fopen(outname,"w");
        if(outvecZ==NULL)
        {
            printf("Error opening outputvecZ\n");
            exit(1);
        }
        for(imom=0;imom<NMOM;imom++)
        {
            for(iconf=100;iconf<115;iconf++)
            {
                sprintf(filename,"/Users/enricofiorenza/Desktop/out/0%d//fft_S_M0_R0_0",iconf);
                propbin=fopen(filename,"r");
                if(propbin==NULL)
                {
                    printf("Error opening input\n");
                    exit(1);
                }
                sprintf(filename,"/Users/enricofiorenza/Desktop/out/0%d//fft_S_M0_R0_RI_VT",iconf);
                vecbinT=fopen(filename,"r");
                if(vecbinT==NULL)
                {
                    printf("Error opening inputvecT\n");
                    exit(1);
                }
                sprintf(filename,"/Users/enricofiorenza/Desktop/out/0%d//fft_S_M0_R0_RI_VX",iconf);
                vecbinX=fopen(filename,"r");
                if(vecbinX==NULL)
                {
                    printf("Error opening inputvecX\n");
                    exit(1);
                }
                sprintf(filename,"/Users/enricofiorenza/Desktop/out/0%d//fft_S_M0_R0_RI_VY",iconf);
                vecbinY=fopen(filename,"r");
                if(vecbinY==NULL)
                {
                    printf("Error opening inputvecY\n");
                    exit(1);
                }
                sprintf(filename,"/Users/enricofiorenza/Desktop/out/0%d//fft_S_M0_R0_RI_VZ",iconf);
                vecbinZ=fopen(filename,"r");
                if(vecbinZ==NULL)
                {
                    printf("Error opening inputvecZ\n");
                    exit(1);
                }
                fseek(propbin,sizeof(double)*288*imom,SEEK_SET);
                fread(array,sizeof(complex double),144,propbin);
                for(is=0;is<NSPIN;is++)
                    for(ic=0;ic<NCOL;ic++)
                        for(is_i=0;is_i<NSPIN;is_i++)
                            for(ic_i=0;ic_i<NCOL;ic_i++)
                            {
                                k=((isc(is_i,ic_i))+(N*isc(is,ic)));
                                S[isc(is_i,ic_i)][isc(is,ic)]=array[k];
                            }
                fseek(vecbinT,sizeof(double)*288*imom,SEEK_SET);
                fread(array,sizeof(complex double),144,vecbinT);
                for(is=0;is<NSPIN;is++)
                    for(ic=0;ic<NCOL;ic++)
                        for(is_i=0;is_i<NSPIN;is_i++)
                            for(ic_i=0;ic_i<NCOL;ic_i++)
                            {
                                k=((isc(is_i,ic_i))+(N*isc(is,ic)));
                                T[isc(is_i,ic_i)][isc(is,ic)]=array[k];
                            }
                fseek(vecbinX,sizeof(double)*288*imom,SEEK_SET);
                fread(array,sizeof(complex double),144,vecbinX);
                for(is=0;is<NSPIN;is++)
                    for(ic=0;ic<NCOL;ic++)
                        for(is_i=0;is_i<NSPIN;is_i++)
                            for(ic_i=0;ic_i<NCOL;ic_i++)
                            {
                                k=((isc(is_i,ic_i))+(N*isc(is,ic)));
                                X[isc(is_i,ic_i)][isc(is,ic)]=array[k];
                            }
                fseek(vecbinY,sizeof(double)*288*imom,SEEK_SET);
                fread(array,sizeof(complex double),144,vecbinY);
                for(is=0;is<NSPIN;is++)
                    for(ic=0;ic<NCOL;ic++)
                        for(is_i=0;is_i<NSPIN;is_i++)
                            for(ic_i=0;ic_i<NCOL;ic_i++)
                            {
                                k=((isc(is_i,ic_i))+(N*isc(is,ic)));
                                Y[isc(is_i,ic_i)][isc(is,ic)]=array[k];
                            }
                fseek(vecbinZ,sizeof(double)*288*imom,SEEK_SET);
                fread(array,sizeof(complex double),144,vecbinZ);
                for(is=0;is<NSPIN;is++)
                    for(ic=0;ic<NCOL;ic++)
                        for(is_i=0;is_i<NSPIN;is_i++)
                            for(ic_i=0;ic_i<NCOL;ic_i++)
                            {
                                k=((isc(is_i,ic_i))+(N*isc(is,ic)));
                                Z[isc(is_i,ic_i)][isc(is,ic)]=array[k];
                            }
                if(iconf!=jack)
                {
                    for(i=0;i<12;i++)
                    {
                        for(j=0;j<12;j++)
                        {
                            Sm[i][j]=Sm[i][j]+S[i][j];
                            Tm[i][j]=Tm[i][j]+T[i][j];
                            Xm[i][j]=Xm[i][j]+X[i][j];
                            Ym[i][j]=Ym[i][j]+Y[i][j];
                            Zm[i][j]=Zm[i][j]+Z[i][j];
                        }
                    }
                }
                fclose(propbin);
                fclose(vecbinT);
                fclose(vecbinX);
                fclose(vecbinY);
                fclose(vecbinZ);
            }
            for(i=0;i<12;i++)
            {
                for(j=0;j<12;j++)
                {
                    Sm[i][j]=Sm[i][j]/14.0;
                    Tm[i][j]=Tm[i][j]/14.0;
                    Xm[i][j]=Xm[i][j]/14.0;
                    Ym[i][j]=Ym[i][j]/14.0;
                    Zm[i][j]=Zm[i][j]/14.0;
                    fprintf(out,"%.16g\n%.16g\n",creal(Sm[i][j]),cimag(Sm[i][j]));
                    fprintf(outvecT,"%.16g\n%.16g\n",creal(Tm[i][j]),cimag(Tm[i][j]));
                    fprintf(outvecX,"%.16g\n%.16g\n",creal(Xm[i][j]),cimag(Xm[i][j]));
                    fprintf(outvecY,"%.16g\n%.16g\n",creal(Ym[i][j]),cimag(Ym[i][j]));
                    fprintf(outvecZ,"%.16g\n%.16g\n",creal(Zm[i][j]),cimag(Zm[i][j]));
                }
            }
            //printf("fine impulso %d\n",imom+1);
            for(i=0;i<12;i++)
            {
                for(j=0;j<12;j++)
                {
                    Sm[i][j]=0;
                    S[i][j]=0;
                    Tm[i][j]=0;
                    T[i][j]=0;
                    Xm[i][j]=0;
                    X[i][j]=0;
                    Ym[i][j]=0;
                    Y[i][j]=0;
                    Zm[i][j]=0;
                    Z[i][j]=0;
                }
            }
        }
        printf("fine file %d\n",jack-99);
        fclose(out);
        fclose(outvecT);
        fclose(outvecX);
        fclose(outvecY);
        fclose(outvecZ);
    }
}

int isc(int is,int ic)
{
    return(ic+NCOL*is);
}

