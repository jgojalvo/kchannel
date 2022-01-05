#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
  int i,j,seed,N,nst,nme,ij,i2,r2,npert,nwarm;
  int *n1,*n2;
  double tt,t,st,dx,dt;
  double gk,gl,Vk,Vl,Sth,Vth,a0,b0,m,sa,ge,dl,dk,bs,gs,gt,Sthm,D,Spert,F,fgk;
  double rnti,Tpert,epert,Twarm,dasym;
  double *V,*ng,*E,*S,*th,*Vi,*ngi,*Ei,*Si,*thi;
  double *f_V,*f_ng,*f_E,*f_S,*f_th,*fi_V,*fi_ng,*fi_E,*fi_S,*fi_th,*rn;
  double dth,dx2,a1,a2;
  char dum[60],name1[60],name2[60],name3[60];
  int ncord(int, int, int);
  int dran_ini(int, int);
  void dran_gv(double *, int);

  FILE *input,*outdat;

/**********  INPUT DATA  **********/

  strcpy(name1,argv[1]);

  strcpy(name3,name1);
  sprintf(name2,".par");
  strcat(name3,name2);
  input = fopen(name3,"r");

  fscanf(input,"%lg %s",&gk,dum);
  fscanf(input,"%lg %s",&gl,dum);
  fscanf(input,"%lg %s",&Vk,dum);
  fscanf(input,"%lg %s",&Vl,dum);
  fscanf(input,"%lg %s",&Sth,dum);
  fscanf(input,"%lg %s",&Vth,dum);
  fscanf(input,"%lg %s",&a0,dum);
  fscanf(input,"%lg %s",&b0,dum);
  fscanf(input,"%lg %s",&m,dum);
  fscanf(input,"%lg %s",&sa,dum);
  fscanf(input,"%lg %s",&ge,dum);
  fscanf(input,"%lg %s",&dl,dum);
  fscanf(input,"%lg %s",&dk,dum);
  fscanf(input,"%lg %s",&bs,dum);
  fscanf(input,"%lg %s",&gs,dum);
  fscanf(input,"%lg %s",&gt,dum);
  fscanf(input,"%lg %s",&F,dum);
  fscanf(input,"%lg %s",&D,dum);
  fscanf(input,"%lg %s",&dasym,dum);
  fscanf(input,"%lg %s",&Twarm,dum);
  fscanf(input,"%lg %s",&Tpert,dum);
  fscanf(input,"%lg %s",&Spert,dum);
  fscanf(input,"%lg %s",&rnti,dum);
  fscanf(input,"%i %s",&N,dum);
  fscanf(input,"%lg %s",&tt,dum);
  fscanf(input,"%lg %s",&st,dum);
  fscanf(input,"%lg %s",&dx,dum);
  fscanf(input,"%lg %s",&dt,dum);
  fscanf(input,"%i %s",&seed,dum);

  fclose(input);

  strcpy(name3,name1);
  sprintf(name2,".dat");
  strcat(name3,name2);
  outdat = fopen(name3,"w");
  fclose(outdat);
  
/**********  SIMULATION CONSTANTS  **********/
  nst=(int)(st/dt+0.5);
  nme=(int)(tt/dt+0.5);
  npert=(int)(Tpert/dt+0.5);
  nwarm=(int)(Twarm/dt+0.5);

  dx2 = 1.0/pow(dx,2.0);
  dth = dt*0.5;
  a1 = -D*dx2;
  a2 = D*dx2/2.0;
  Sthm = pow(Sth,m);
  fgk = F*gk;

  V=(double *)calloc(N,sizeof(double));
  ng=(double *)calloc(N,sizeof(double));
  E=(double *)calloc(N+1,sizeof(double));
  S=(double *)calloc(N,sizeof(double));
  th=(double *)calloc(N,sizeof(double));
  Vi=(double *)calloc(N,sizeof(double));
  ngi=(double *)calloc(N,sizeof(double));
  Ei=(double *)calloc(N+1,sizeof(double));
  Si=(double *)calloc(N,sizeof(double));
  thi=(double *)calloc(N,sizeof(double));
  f_V=(double *)calloc(N,sizeof(double));
  f_ng=(double *)calloc(N,sizeof(double));
  f_E=(double *)calloc(N,sizeof(double));
  f_S=(double *)calloc(N,sizeof(double));
  f_th=(double *)calloc(N,sizeof(double));
  fi_V=(double *)calloc(N,sizeof(double));
  fi_ng=(double *)calloc(N,sizeof(double));
  fi_E=(double *)calloc(N,sizeof(double));
  fi_S=(double *)calloc(N,sizeof(double));
  fi_th=(double *)calloc(N,sizeof(double));
  rn=(double *)calloc(N,sizeof(double));

  n1=(int *)calloc(N,sizeof(int));
  n2=(int *)calloc(N,sizeof(int));

  for (i=0;i<N;i++)
  {
    n1[i]=ncord(N,i,1);
    n2[i]=ncord(N,i,-1);
  }

/**********  INITIAL CONDITIONS  **********/

  if (dran_ini(seed,N)==1) return 1;

  for (i=0;i<N;i++)
  {
     V[i]=0;
     ng[i]=0;
     E[i]=0;
     S[i]=0;
     th[i]=0;
  }
  E[N] = 0;
  Ei[N] = 0;
  
  // ek[0]=2.8;
  // ik[0]=50;
  // th[0]=at/(gi+gt);

/**** I.C.: random ****/

  dran_gv(rn,N);
  for (i=0;i<N;i++)
      V[i]+=exp(rn[i])*rnti;
  dran_gv(rn,N);
  for (i=0;i<N;i++)
      ng[i]+=exp(rn[i])*rnti;
  dran_gv(rn,N);
  for (i=0;i<N;i++)
      E[i]+=exp(rn[i])*rnti;
  dran_gv(rn,N);
  for (i=0;i<N;i++)
      S[i]+=exp(rn[i])*rnti;
  dran_gv(rn,N);
  for (i=0;i<N;i++)
      th[i]+=exp(rn[i])*rnti;

/**********  SIMULATION  **********/

  for (ij=1;ij<=nme;ij++)
  {
     if ((ij-nwarm)%npert==0 && ij>nwarm)
        for (i=0;i<1;i++)
           S[i] += Spert;
     for (i=0;i<N;i++)
        f_V[i] = -gk*pow(ng[i],4)*(V[i]-Vk-dk*E[i])-gl*(V[i]-Vl-dl*E[i]);
     for (i=0;i<N;i++)
        f_ng[i] = a0*pow(S[i],m)/(Sthm+pow(S[i],m))*(1-ng[i])-b0*ng[i];
     for (i=0;i<N;i++)
        f_E[i] = fgk*pow(ng[i],4.0)*(V[i]-Vk-dk*E[i])-ge*E[i]
                 +a1*E[i]+a2*(E[n1[i]]*(1-dasym)+E[n2[i]]*(1+dasym));
     for (i=0;i<N;i++)
        f_S[i]=bs*(Vth-V[i])/(exp((Vth-V[i])/sa)-1)-gs*S[i];              
     for (i=0;i<N;i++)
        f_th[i]=bs*(Vl-V[i])-gt*th[i]; 
     for (i=0;i<N;i++)
     {
        Vi[i]=V[i]+dt*f_V[i];
        ngi[i]=ng[i]+dt*f_ng[i];
        Ei[i]=E[i]+dt*f_E[i];
        Si[i]=S[i]+dt*f_S[i];
        thi[i]=th[i]+dt*f_th[i];
     }
        
     for (i=0;i<N;i++)
        fi_V[i] = -gk*pow(ngi[i],4)*(Vi[i]-Vk-dk*Ei[i])-gl*(Vi[i]-Vl-dl*Ei[i]);
     for (i=0;i<N;i++)
        fi_ng[i] = a0*pow(Si[i],m)/(Sthm+pow(Si[i],m))*(1-ngi[i])-b0*ngi[i];
     for (i=0;i<N;i++)
        fi_E[i] = fgk*pow(ngi[i],4)*(Vi[i]-Vk-dk*Ei[i])-ge*Ei[i]
                  +a1*Ei[i]+a2*(Ei[n1[i]]*(1-dasym)+Ei[n2[i]]*(1+dasym));
     for (i=0;i<N;i++)
        fi_S[i]=bs*(Vth-Vi[i])/(exp((Vth-Vi[i])/sa)-1)-gs*Si[i];              
     for (i=0;i<N;i++)
        fi_th[i]=bs*(Vl-Vi[i])-gt*thi[i]; 
     for (i=0;i<N;i++)
     {
        V[i]=V[i]+dth*(f_V[i]+fi_V[i]);
        ng[i]=ng[i]+dth*(f_ng[i]+fi_ng[i]);
        E[i]=E[i]+dth*(f_E[i]+fi_E[i]);
        S[i]=S[i]+dth*(f_S[i]+fi_S[i]);
        th[i]=th[i]+dth*(f_th[i]+fi_th[i]);
     }

     if (ij%nst==0)
     {
        t=dt*(double)ij;

        strcpy(name3,name1);
        sprintf(name2,".dat");
        strcat(name3,name2);
        outdat = fopen(name3,"a");
        fprintf(outdat,"%f\t",t);
        for (i=0;i<N;i++)
          fprintf(outdat,"%f\t%f\t%f\t%f\t%f\t",V[i],ng[i],E[i],S[i],th[i]);
        fprintf(outdat,"\n");
        //printf("%f\t%f\t%f\t%f\t%f\t%f\n",t,V[0],ng[0],E[0],S[0],th[0]);
        fclose(outdat);
     }
  } 

  printf("\nQuitting.\n\n");

  return 0;
}


int ncord(int l, int i, int ix)
{
   int nc,i1;

   i1=i+ix;

/* Periodic boundary conditions */
//   if (i1>=l) i1-=l;
//   if (i1<0) i1+=l;

/* No-flux boundary conditions */
//   if (i1>=l) i1=l-1;
//   if (i1<0) i1=0;

/* Free boundary conditions */
   if (i1>=l) i1=l;
   if (i1<0) i1=l;

   nc=i1;

   return nc;
}
