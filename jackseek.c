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
    FILE *propbin,*out;
    char filename[80],outname[80];
    double complex S[12][12],Sm[12][12];
    double complex *array;
    int iconf,i,j,k,p_0,p_1,p_2,p_3,imom,is,ic,is_i,ic_i,jack;
    array=(complex double*) calloc(144,sizeof(complex double));
    for(i=0;i<12;i++)
    {
        for(j=0;j<12;j++)
        {
            S[i][j]=0+0*I;
            Sm[i][j]=0+0*I;
        }
    }
    for(jack=100;jack<115;jack++)
    {
        sprintf(outname,"/Users/enricofiorenza/Desktop/jackknife//prop%d",jack);
        out=fopen(outname,"w");
        if(out==NULL)
        {
            printf("Error opening output\n");
            exit(1);
        }
        for(imom=0;imom<NMOM;imom++)
        {
            for(iconf=100;iconf<115;iconf++)
            {
                sprintf(filename,"/Users/enricofiorenza/Desktop/Tesi/fac_form/out/0%d//fft_S_M0_R0_0",iconf);
                propbin=fopen(filename,"r");
                if(out==NULL)
                {
                    printf("Error opening input\n");
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
                if(iconf!=jack)
                {
                    for(i=0;i<12;i++)
                    {
                        for(j=0;j<12;j++)
                        {
                            Sm[i][j]=Sm[i][j]+S[i][j];
                        }
                    }
                }
                fclose(propbin);
            }
            for(i=0;i<12;i++)
            {
                for(j=0;j<12;j++)
                {
                    Sm[i][j]=Sm[i][j]/14;
                    fprintf(out,"%.16g\n%.16g\n",creal(Sm[i][j]),cimag(Sm[i][j]));
                }
            }
            printf("fine impulso %d\n",imom+1);
            for(i=0;i<12;i++)
            {
                for(j=0;j<12;j++)
                {
                    Sm[i][j]=0;
                    S[i][j]=0;
                }
            }
        }
        fclose(out);
    }
}

int isc(int is,int ic)
{
    return(ic+NCOL*is);
}

