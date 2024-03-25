/*********************************************************************
 * model.cpp
 * Keep in mind:
 * <> Use 0-based indexing as always in C or C++
 * <> Indexing is column-based as in Matlab (not row-based as in C)
 * <> Use linear indexing.  [x*dimy+y] instead of [x][y]
 * Adapted from the code by Sen Pei (http://www.columbia.edu/~sp3449/)
 * and Shawn Lankton (http://www.shawnlankton.com/2008/03/getting-started-with-mex-a-short-tutorial/)
 ********************************************************************/
#include <matrix.h>
#include <mex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <random>
#include <vector>
#include <algorithm> 
using namespace std; 

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

void mexFunction(int nlmxhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    mxArray *mxnl, *mxpart, *mxC, *mxCave, *mxSm, *mxSh, *mxIm, *mxIhr, *mxIhu, *mxpara, *mxtaos, *mxbetahtom, *mxbetamtoh, *mxbirth_rate, *mxbetamap, *mxalphamap;//input
    mxArray *mxnewSm, *mxnewSh, *mxnewIm, *mxnewIhr, *mxnewIhu, *mxdailyIm, *mxdailyIhr, *mxdailyIhu;//output
    const mwSize *dims;
    double *nl, *part, *C, *Cave, *Sm, *Sh, *Im, *Ihr, *Ihu, *para, *taos, *betahtom, *betamtoh, *birth_rate, *betamap, *alphamap;//input
    double *newSm, *newSh, *newIm, *newIhr, *newIhu, *dailyIm, *dailyIhr, *dailyIhu;//output
    int num_mp, num_loc, num_para;

//associate inputs
    mxnl = mxDuplicateArray(prhs[0]);
    mxpart = mxDuplicateArray(prhs[1]);
    mxC = mxDuplicateArray(prhs[2]);
    mxCave = mxDuplicateArray(prhs[3]);
    mxSm = mxDuplicateArray(prhs[4]);
    mxSh = mxDuplicateArray(prhs[5]);
    mxIm = mxDuplicateArray(prhs[6]);
    mxIhr = mxDuplicateArray(prhs[7]);
    mxIhu = mxDuplicateArray(prhs[8]);
    mxpara = mxDuplicateArray(prhs[9]);
    mxtaos = mxDuplicateArray(prhs[10]);
    mxbetahtom = mxDuplicateArray(prhs[11]);
    mxbetamtoh = mxDuplicateArray(prhs[12]);
    mxbirth_rate = mxDuplicateArray(prhs[13]);
    mxbetamap = mxDuplicateArray(prhs[14]);
    mxalphamap = mxDuplicateArray(prhs[15]);
    
//figure out dimensions
    dims = mxGetDimensions(prhs[0]);//number of subpopulation
    num_mp = (int)dims[0];
    dims = mxGetDimensions(prhs[1]);//number of locations
    num_loc = (int)dims[0]-1;
    dims = mxGetDimensions(prhs[9]);//number of parameters
    num_para = (int)dims[0];

//associate outputs
    mxnewSm = plhs[0] = mxCreateDoubleMatrix(num_loc,1,mxREAL);
    mxnewSh = plhs[1] = mxCreateDoubleMatrix(num_mp,1,mxREAL);
    mxnewIm = plhs[2] = mxCreateDoubleMatrix(num_loc,1,mxREAL);
    mxnewIhr = plhs[3] = mxCreateDoubleMatrix(num_mp,1,mxREAL);
    mxnewIhu = plhs[4] = mxCreateDoubleMatrix(num_mp,1,mxREAL);
    mxdailyIm = plhs[5] = mxCreateDoubleMatrix(num_loc,1,mxREAL);
    mxdailyIhr = plhs[6] = mxCreateDoubleMatrix(num_mp,1,mxREAL);
    mxdailyIhu = plhs[7] = mxCreateDoubleMatrix(num_mp,1,mxREAL);

//associate pointers
    nl = mxGetPr(mxnl);
    part = mxGetPr(mxpart);
    C = mxGetPr(mxC);
    Cave = mxGetPr(mxCave);
    Sm = mxGetPr(mxSm);
    Sh = mxGetPr(mxSh);
    Im = mxGetPr(mxIm);
    Ihr = mxGetPr(mxIhr);
    Ihu = mxGetPr(mxIhu);
    para = mxGetPr(mxpara);
    taos = mxGetPr(mxtaos);
    betahtom = mxGetPr(mxbetahtom);
    betamtoh = mxGetPr(mxbetamtoh);
    birth_rate = mxGetPr(mxbirth_rate);
    betamap = mxGetPr(mxbetamap);
    alphamap = mxGetPr(mxalphamap);
    newSm = mxGetPr(mxnewSm);
    newSh = mxGetPr(mxnewSh);
    newIm = mxGetPr(mxnewIm);
    newIhr = mxGetPr(mxnewIhr);
    newIhu = mxGetPr(mxnewIhu);
    dailyIm = mxGetPr(mxdailyIm);    
    dailyIhr = mxGetPr(mxdailyIhr);
    dailyIhu = mxGetPr(mxdailyIhu);

    ////////////////////////////////////////
    //do something
    default_random_engine generator((unsigned)time(NULL));
    //initialize auxillary variables
    int i, j;
    //para=U,death,D,mu,theta,alpha1,alpha2,alpha3,beta0,beta1,...
    double U=para[0],death=para[1],D=para[2],mu=para[3],theta=para[4];
    double dt1=(double)1/3, dt2=1-dt1;
    //daytime population,Im,Ihr,Ihu,Sm,Sh,R
    vector<double> ND(num_loc),ImD(num_loc),IhrD(num_loc),IhuD(num_loc),SmD(num_loc),ShD(num_loc),RD(num_loc);
    //daytime enter Sh, Ihu
    vector<double> ShentD(num_loc),IhuentD(num_loc);
    //nighttime population,Im,Ir,Iu,Sm,Sh,R
    vector<double> NN(num_loc),ImN(num_loc),IhrN(num_loc),IhuN(num_loc),SmN(num_loc),ShN(num_loc),RN(num_loc);
    //nighttime enter Sh, Ihu
    vector<double> ShentN(num_loc),IhuentN(num_loc);
    //total outgoing population
    vector<double> popleft(num_loc);
    //intermediate Sm,Sh, Im, Ihu, Ihr, R
    vector<double> tempSm(num_loc),tempSh(num_mp),tempIm(num_loc),tempIhu(num_mp),tempIhr(num_mp),tempR(num_mp);
    //R and newR
    vector<double> R(num_mp), newR(num_mp);
    
    /////////////////////////////
    //change index in nl and part (0-based index)
    for (i=0; i<num_mp; i++)
        nl[i]=nl[i]-1;
    for (i=0; i<num_loc+1; i++)
        part[i]=part[i]-1;
    for (i=0; i<num_loc; i++){
        betamap[i]=betamap[i]-1;
        alphamap[i]=alphamap[i]-1;
    }
    //compute popleft
    for (i=0; i<num_loc; i++){
        for (j=part[i]+1; j<part[i+1]; j++){
            popleft[i]=popleft[i]+Cave[j];
        }
    }
    
    
    //assgn intermediate Sm, Sh, Im, Ihu, Ihr, R
    for (i=0; i<num_mp; i++){
        R[i] = max(C[i]-Sh[i]-Ihr[i]-Ihu[i],0.0);
        tempSh[i] = Sh[i];
        tempIhr[i] = Ihr[i];
        tempIhu[i] = Ihu[i];
        tempR[i] = R[i];
    }
    
    for (i=0; i<num_loc; i++){
        tempSm[i] = Sm[i];
        tempIm[i] = Im[i];
    }

    ///////////////////////////
    //daytime transmission
    //compute ND
    for (i=0; i<num_loc; i++){
        ND[i]=C[(int)(part[i])];//i<-i
        for (j=part[i]+1; j<part[i+1]; j++){
            ND[i]=ND[i]+Ihr[j];//reported infections (no mobility)
        }
    }
    for (i=0; i<num_loc; i++){
        for (j=part[i]+1; j<part[i+1]; j++){
            ND[(int)nl[j]]=ND[(int)nl[j]]+C[j]-Ihr[j];//commuting with reported infections removed
        }
    }

    
    //comput IhrD,IhuD,ShD,RD,SmD,ImD
    for (i=0; i<num_loc; i++){
        for (j=part[i]; j<part[i+1]; j++){
            IhrD[i]=IhrD[i]+Ihr[j];
            IhuD[(int)nl[j]]=IhuD[(int)nl[j]]+Ihu[j];
            ShD[(int)nl[j]]=ShD[(int)nl[j]]+Sh[j];
            RD[(int)nl[j]]=RD[(int)nl[j]]+R[j];
        }
        SmD[i]=Sm[i];
        ImD[i]=Im[i];
    }
    
    //compute ShentD and IuentD
    for (i=0; i<num_loc; i++){
        for (j=part[i]+1; j<part[i+1]; j++){
            ShentD[(int)nl[j]]=ShentD[(int)nl[j]]+Cave[j]*ShD[i]/(ND[i]-IhrD[i]);
            IhuentD[(int)nl[j]]=IhuentD[(int)nl[j]]+Cave[j]*IhuD[i]/(ND[i]-IhrD[i]);
        }
    }
    //compute for each subpopulation
    for (i=0; i<num_loc; i++){
        for (j=part[i]; j<part[i+1]; j++){
            //////////////////////////////
            double beta=para[(int)betamap[(int)nl[j]]];
            double alpha=para[(int)alphamap[(int)nl[j]]];
            double EImtoIh=beta*taos[(int)nl[j]]*betamtoh[(int)nl[j]]*Sh[j]*ImD[(int)nl[j]]/ND[(int)nl[j]]*dt1;//
            double EImtoIhr=alpha*beta*taos[(int)nl[j]]*betamtoh[(int)nl[j]]*Sh[j]*ImD[(int)nl[j]]/ND[(int)nl[j]]*dt1;
            double EImtoIhu=(1-alpha)*beta*taos[(int)nl[j]]*betamtoh[(int)nl[j]]*Sh[j]*ImD[(int)nl[j]]/ND[(int)nl[j]]*dt1;
            double EIhrtor=Ihr[j]/D*dt1;
            double EIhutor=Ihu[j]/D*dt1;
            ///////////////////////////
            double EShenter=theta*dt1*(C[j]-Ihr[j])/ND[(int)nl[j]]*ShentD[(int)nl[j]];//incoming S
            double EShleft=theta*dt1*Sh[j]/(ND[(int)nl[j]]-IhrD[(int)nl[j]])*popleft[(int)nl[j]];//outgoing S
            double EIhuenter=theta*dt1*(C[j]-Ihr[j])/ND[(int)nl[j]]*IhuentD[(int)nl[j]];//incoming Iu
            double EIhuleft=theta*dt1*Ihu[j]/(ND[(int)nl[j]]-IhrD[(int)nl[j]])*popleft[(int)nl[j]];//outgoing Iu
            
            ////////////////////
            //stochastic: poisson
            poisson_distribution<int> distribution1(EImtoIh);
            int eImtoIh=min(distribution1(generator),(int)(Sh[j]*dt1));
            poisson_distribution<int> distribution2(EImtoIhr);
            int eImtoIhr=min(distribution2(generator),(int)(Sh[j]*dt1));
            poisson_distribution<int> distribution3(EImtoIhu);
            int eImtoIhu=min(distribution3(generator),(int)(Sh[j]*dt1));
            poisson_distribution<int> distribution4(EIhrtor);
            int eIhrtor=min(distribution4(generator),(int)(Ihr[j]*dt1));
            poisson_distribution<int> distribution5(EIhutor);
            int eIhutor=min(distribution5(generator),(int)(Ihu[j]*dt1));
            //////////////////////////
            poisson_distribution<int> distribution6(EShenter);
            int eShenter=distribution6(generator);
            poisson_distribution<int> distribution7(EShleft);
            int eShleft=min(distribution7(generator),(int)(Sh[j]*dt1));
            poisson_distribution<int> distribution8(EIhuenter);
            int eIhuenter=distribution8(generator);
            poisson_distribution<int> distribution9(EIhuleft);
            int eIhuleft=min(distribution9(generator),(int)(Ihr[j]*dt1));

            
            /////////////////////
            tempSh[j]=max((int)(tempSh[j]-eImtoIh+eShenter-eShleft),0);
            tempIhr[j]=max((int)(tempIhr[j]+eImtoIhr-eIhrtor),0);
            tempIhu[j]=max((int)(tempIhu[j]+eImtoIhu-eIhutor+eIhuenter-eIhuleft),0);
            dailyIhr[j]=max((int)(dailyIhr[j]+eImtoIhr),0);
            dailyIhu[j]=max((int)(dailyIhu[j]+eImtoIhu),0);
            tempR[j]=max((int)(C[j]-tempIhr[j]-tempIhu[j]-tempSh[j]),0);
        }
    }

    for (i=0; i<num_loc; i++){
        ///////////////////////////////
        double beta=para[(int)betamap[i]];
        double EIhrtoSm=beta*taos[i]*betahtom[i]*SmD[i]*IhrD[i]/ND[i]*dt1;
        double EIhutoSm=mu*beta*taos[i]*betahtom[i]*SmD[i]*IhuD[i]/ND[i]*dt1;
        double Ebirthmu=birth_rate[i]*SmD[i]*dt1;
        double EbirthmU=birth_rate[i]*(1-U)*ImD[i]*dt1;
        double EdeathSm=death*SmD[i]*dt1;
        double EbirthmUU=birth_rate[i]*U*ImD[i]*dt1;
        double EdeathIm=death*ImD[i]*dt1;

        ///////////////////////////////////
        poisson_distribution<int> distribution1(EIhrtoSm);
        int eIhrtoSm=min(distribution1(generator),(int)(Sm[i]*dt1));
        poisson_distribution<int> distribution2(EIhutoSm);
        int eIhutoSm=min(distribution2(generator),(int)(Sm[i]*dt1));
        poisson_distribution<int> distribution3(Ebirthmu);
        int ebirthmu=min(distribution3(generator),(int)(Sm[i]*dt1));
        poisson_distribution<int> distribution4(EbirthmU);
        int ebirthmU=min(distribution4(generator),(int)(Im[i]*dt1));
        poisson_distribution<int> distribution5(EdeathSm);
        int edeathSm=min(distribution5(generator),(int)(Sm[i]*dt1));
        poisson_distribution<int> distribution6(EbirthmUU);
        int ebirthmUU=min(distribution6(generator),(int)(Im[i]*dt1));
        poisson_distribution<int> distribution7(EdeathIm);
        int edeathIm=min(distribution7(generator),(int)(Im[i]*dt1));

        /////////////////////
        tempSm[i]=max((int)(tempSm[i]-eIhrtoSm-eIhutoSm+ebirthmu+ebirthmU-edeathSm),0);
        tempIm[i]=max((int)(tempIm[i]+eIhrtoSm+eIhutoSm+ebirthmUU-edeathIm),0);
        dailyIm[i]=max((int)(dailyIm[i]+eIhrtoSm+eIhutoSm+ebirthmUU),0);
    }
    ////////////////////////////////
    //nighttime transmission
    //assgn new Sh, Ihu, Ihr, R, Sm, Im
    for (i=0; i<num_mp; i++){
        newSh[i] = tempSh[i];
        newIhr[i] = tempIhr[i];
        newIhu[i] = tempIhu[i];
        newR[i] = tempR[i];
    }

    for (i=0; i<num_loc; i++){
        newSm[i] = tempSm[i];
        newIm[i] = tempIm[i];
    }
    //compute NN
    for (i=0; i<num_loc; i++){
        for (j=part[i]; j<part[i+1]; j++){
            NN[i]=NN[i]+C[j];
        }
    }
    //comput IhrN, IhuN, ShN, RN, Sm, Im
    for (i=0; i<num_loc; i++){
        for (j=part[i]; j<part[i+1]; j++){
            IhrN[i]=IhrN[i]+tempIhr[j];
            IhuN[i]=IhuN[i]+tempIhu[j];
            ShN[i]=ShN[i]+tempSh[j];
            RN[i]=RN[i]+tempR[j];
        }
        SmN[i]=tempSm[i];
        ImN[i]=tempIm[i];
    }
    //compute SentN and IhuentN
    for (i=0; i<num_loc; i++){
        for (j=part[i]+1; j<part[i+1]; j++){
            ShentN[(int)nl[j]]=ShentN[(int)nl[j]]+Cave[j]*ShN[i]/(NN[i]-IhrN[i]);
            IhuentN[(int)nl[j]]=IhuentN[(int)nl[j]]+Cave[j]*IhuN[i]/(NN[i]-IhrN[i]);
        }
    }
    //compute for each subpopulation
    for (i=0; i<num_loc; i++){
        for (j=part[i]; j<part[i+1]; j++){
            //////////////////////////////
            double beta=para[(int)betamap[i]];
            double alpha=para[(int)alphamap[i]];
            double EImtoIh=beta*taos[i]*betamtoh[i]*tempSh[j]*tempIm[i]/NN[i]*dt2;//
            double EImtoIhr=alpha*beta*taos[i]*betamtoh[i]*tempSh[j]*tempIm[i]/NN[i]*dt2;
            double EImtoIhu=(1-alpha)*beta*taos[i]*betamtoh[i]*tempSh[j]*tempIm[i]/NN[i]*dt2;
            double EIhrtor=Ihr[j]/D*dt2;
            double EIhutor=Ihu[j]/D*dt2;
            ///////////////////////////
            double EShenter=theta*dt2*C[j]/NN[i]*ShentN[i];//incoming S
            double EShleft=theta*dt2*tempSh[j]/(NN[i]-IhrN[i])*popleft[i];//outgoing S
            double EIhuenter=theta*dt2*C[j]/NN[i]*IhuentN[i];//incoming Iu
            double EIhuleft=theta*dt2*tempIhu[j]/(NN[i]-IhrN[i])*popleft[i];//outgoing Iu
            
            ////////////////////
            //stochastic: poisson
            poisson_distribution<int> distribution1(EImtoIh);
            int eImtoIh=min(distribution1(generator),(int)(tempSh[j]*dt2));
            poisson_distribution<int> distribution2(EImtoIhr);
            int eImtoIhr=min(distribution2(generator),(int)(tempSh[j]*dt2));
            poisson_distribution<int> distribution3(EImtoIhu);
            int eImtoIhu=min(distribution3(generator),(int)(tempSh[j]*dt2));
            poisson_distribution<int> distribution4(EIhrtor);
            int eIhrtor=min(distribution4(generator),(int)(tempIhr[j]*dt2));
            poisson_distribution<int> distribution5(EIhutor);
            int eIhutor=min(distribution5(generator),(int)(tempIhu[j]*dt2));
            //////////////////////////
            poisson_distribution<int> distribution6(EShenter);
            int eShenter=distribution6(generator);
            poisson_distribution<int> distribution7(EShleft);
            int eShleft=min(distribution7(generator),(int)(tempSh[j]*dt2));
            poisson_distribution<int> distribution8(EIhuenter);
            int eIhuenter=distribution8(generator);
            poisson_distribution<int> distribution9(EIhuleft);
            int eIhuleft=min(distribution9(generator),(int)(tempIhr[j]*dt2));

            
            /////////////////////
            newSh[j]=max((int)(newSh[j]-eImtoIh+eShenter-eShleft),0);
            newIhr[j]=max((int)(newIhr[j]+eImtoIhr-eIhrtor),0);
            newIhu[j]=max((int)(newIhu[j]+eImtoIhu-eIhutor+eIhuenter-eIhuleft),0);
            dailyIhr[j]=max((int)(dailyIhr[j]+eImtoIhr),0);
            dailyIhu[j]=max((int)(dailyIhu[j]+eImtoIhu),0);
            newR[j]=max((int)(C[j]-newIhr[j]-newIhu[j]-newSh[j]),0);
        }
    }

    for (i=0; i<num_loc; i++){
        ///////////////////////////////
        double beta=para[(int)betamap[i]];
        double EIhrtoSm=beta*taos[i]*betahtom[i]*SmN[i]*IhrN[i]/NN[i]*dt2;
        double EIhutoSm=mu*beta*taos[i]*betahtom[i]*SmN[i]*IhuN[i]/NN[i]*dt2;
        double Ebirthmu=birth_rate[i]*SmN[i]*dt2;
        double EbirthmU=birth_rate[i]*(1-U)*ImN[i]*dt2;
        double EdeathSm=death*SmN[i]*dt2;
        double EbirthmUU=birth_rate[i]*U*ImN[i]*dt2;
        double EdeathIm=death*ImN[i]*dt2;

        ///////////////////////////////////
        poisson_distribution<int> distribution1(EIhrtoSm);
        int eIhrtoSm=min(distribution1(generator),(int)(tempSm[i]*dt2));
        poisson_distribution<int> distribution2(EIhutoSm);
        int eIhutoSm=min(distribution2(generator),(int)(tempSm[i]*dt2));
        poisson_distribution<int> distribution3(Ebirthmu);
        int ebirthmu=min(distribution3(generator),(int)(tempSm[i]*dt2));
        poisson_distribution<int> distribution4(EbirthmU);
        int ebirthmU=min(distribution4(generator),(int)(tempIm[i]*dt2));
        poisson_distribution<int> distribution5(EdeathSm);
        int edeathSm=min(distribution5(generator),(int)(tempSm[i]*dt2));
        poisson_distribution<int> distribution6(EbirthmUU);
        int ebirthmUU=min(distribution6(generator),(int)(tempIm[i]*dt2));
        poisson_distribution<int> distribution7(EdeathIm);
        int edeathIm=min(distribution7(generator),(int)(tempIm[i]*dt2));

        /////////////////////
        newSm[i]=max((int)(newSm[i]-eIhrtoSm-eIhutoSm+ebirthmu+ebirthmU-edeathSm),0);
        newIm[i]=max((int)(newIm[i]+eIhrtoSm+eIhutoSm+ebirthmUU-edeathIm),0);
        dailyIm[i]=max((int)(dailyIm[i]+eIhrtoSm+eIhutoSm+ebirthmUU),0);
    }
    /////////////////////
    return;
}
