////////////////   poiseuille ， nonconvex scheme

#include<math.h> 
#include<cstdlib> 
//#include<iomanip.h> 

#include<string.h> 
#include<stdio.h> 
 

const int Q=9; 

const int NY=11; 
const int NX=20; 

 
const int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};    // 9 direction
const int ne[Q]={0,3,4,1,2,7,8,5,6};                                                // back direction
const double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36}; 
double rho[NX+1][NY+1],u[NX+1][NY+1][2],f[NX+1][NY+1][Q],F[NX+1][NY+1][Q],ff[NX+1][NY+1][Q],xlabel[NX+1][NY+1],ylabel[NX+1][NY+1],u_exact[NY+1]; 
double u0[NX+1][NY+1][2];


int i,j,k,ip,jp,n; 
double c,Re,dx,dy,Lx,Ly,dt,rho0,P0,s_nu,s_q,SS,niu,error,Force,U_max,U_s,cs_2,rela_error,ERROR; 
double a1,a2,a3,a4,a5,q_up,q_down;
 
void init(); 
double feq(int k,double rho,double u[2]); 
void evolution(); 
void output(int m); 
void Error(); 
int flag[NX+1][NY+1];
void Outdata(int m);
 
int main() 
{ 
  for (i=0;i<=NX;i++)    //tag different boundary
	  for (j=0;j<=NY;j++)
	  {
		  if      (j==0)   flag[i][j]=1;  //Up and Down Boundary
		  else if (j==NY)  flag[i][j]=5;  //Up and Down Boundary

		  else if (i==0)   flag[i][j]=2;
          else if (i==NX)  flag[i][j]=3; 
		  else  		   flag[i][j]=4;

	  }

         printf("flag [0][NY]=%d\n",flag[0][NY]);
         printf("flag [NX][NY]=%d\n",flag[NX][NY]);
         printf("flag [0][0]=%d\n",flag[0][0]);
         printf("flag [NX][0]=%d\n",flag[NX][0]);
         printf("flag [0][3]=%d\n",flag[0][3]);
         printf("flag [NX][3]=%d\n",flag[NX][3]);


  init();  // initate
  printf("U_max/c=%f\n",U_max/c);

 for(n=0;;n++) 
  { 

    evolution(); 

    if(n%1000==0) 
    { 
      Error(); 
     
      printf("The %d th computationa result:\n",n);
      
      printf("The relative error of 1000 steps is: %.13f ,  error=%.10f  \n",ERROR, error);
      
    }
      
    if(n>0 && fabs(ERROR)<1.0e-10) 
     { 
	   Outdata(n);
       printf("NY=%d, q_up=%f , tau=%f , error=%.14f \n",NX,q_up,SS,error);

       break; 
     }   
  } 

  return 0;      
} 



void init() 
{ 
  
  //------------------------------------------------
  //------------------------------------These are for dt=dx
  //dx=1.0;
  //dy=dx;
  //c=1.0; 
  //dt=dx;
  //Lx=NX*dx;
  //Ly=(NY-2)*dy+(q_up+q_down)*dy;
  //s_nu=-1.0/0.55;
  //s_q=-8.0*(2+s_nu)/(8+s_nu); //!!!!!
  //SS=-1.0/s_nu; //!!!!!!!
  //niu=(-1.0/s_nu-0.5)/3.0; //!!!!!!
  //------------------------------------------------


  //-----------------------------------These are for dt=a*dx^2
  q_up=0.25;     //上边界的网格步长
  q_down=q_up;
  SS=3.0;            //tau



  Ly=1.0;
  dy=Ly/( NY-2.0 + q_up + q_down  ); 
  dx=dy;
  Lx=NX*dx;



  niu=0.003;
  s_nu=-1.0/SS;
  s_q=-8.0*(2.0+s_nu)/(8.0+s_nu);             
  dt=(SS -0.5)/3.0 *dx*dx /niu;
  c=dx/dt;





  rho0=1.0;  
  U_max=0.1;


 
  Re=U_max*Ly/niu;  
  
  cs_2=c*c/3.0;

  Force=8.0*niu*rho0*U_max/Ly/Ly;



  printf("s_nu=%f, s_q=%f \n",s_nu,s_q);

  printf("U/c=%f",U_max/c);

  //cout<<"tau_f= "<<tau_f<<endl; 

  

  for(i=0;i<=NX;i++)
	  for(j=0;j<=NY;j++)
	  {

		   xlabel[i][j]=i*dx;

		   ylabel[i][j]=(j-1)*dy+q_down*dy;

		   ylabel[i][0]=0.0;
		   ylabel[i][NY] = Ly;

	  }

   for(j=0;j<=NY;j++)
   {
	   u_exact[j]=4.0*U_max*(1.0-ylabel[0][j]/Ly)*ylabel[0][j]/Ly;  //analytical solution
   }
 


 

  for(i=0;i<=NX;i++)    //  initialization of DF
    for(j=0;j<=NY;j++) 
    { 
      u[i][j][0]=0; 
      u[i][j][1]=0; 
   
      rho[i][j]=1.0;
     
      for(k=0;k<Q;k++) 
       { 
         f[i][j][k]=feq(k,rho[i][j],u[i][j]); 
       }
 
    } 
} 

 
double feq(int k,double rho,double u[2])   //  equilibrium distribution
{ 
  double eu,uv,feq; 
  eu=(e[k][0]*u[0]+e[k][1]*u[1]); 
  uv=(u[0]*u[0]+u[1]*u[1]); 
  feq=w[k]*rho*(1.0+3.0*eu/c+4.5*eu*eu/(c*c)-1.5*uv/(c*c)); 
  return feq; 
} 

 

void evolution() 
{ 

   for(j=NY-1;j>0;j--) 
	{
		for(i=0;i<=NX;i++) 
	    {
			for(k=0;k<Q;k++) 
				F[i][j][k]=f[i][j][k]
				           +
				           s_nu*(  0.5*(f[i][j][k]+f[i][j][ne[k]]) - 0.5*(feq(k,rho[i][j], u[i][j])+feq(ne[k],rho[i][j], u[i][j]))   )
						   +
				           s_q*(  0.5*(f[i][j][k]-f[i][j][ne[k]]) - 0.5*(feq(k,rho[i][j], u[i][j])-feq(ne[k],rho[i][j], u[i][j]))   )				
				
				           +dt*w[k]*Force*c*e[k][0]/cs_2;
				           
	    }
	}


    for(j=NY-1;j>0;j--) 
		for(i=0;i<=NX;i++) 
		{
			
				for(k=0;k<Q;k++) 
				{
					ip=i-e[k][0]; 
					ip=(ip+NX+1)%(NX+1);

					jp=j-e[k][1];
					
					if(flag[ip][jp]==4||flag[ip][jp]==2||flag[ip][jp]==3)    //inner fluid
						ff[i][j][k]=F[ip][jp][k];



					
					
					else if(flag[ip][jp]==5)

					{
						a1=2.0*q_up;
						a2=1.0-a1;

						ff[i][j][k] = a1 * ( F[i][j][ne[k]]-dt*w[ne[k]]*Force*c*e[ne[k]][0]/cs_2)  + a2*f[i][j][ne[k]];

                    }


					else if(flag[ip][jp]==1)

					{
						a1=2.0*q_down;
						a2=1.0-a1;

						ff[i][j][k] = a1* ( F[i][j][ne[k]]-dt*w[ne[k]]*Force*c*e[ne[k]][0]/cs_2)  + a2*f[i][j][ne[k]];

                    }
                                               

                              
				}
		}
   



   for(i=0;i<=NX;i++) 
        for(j=1;j<NY;j++) 
        { 

          rho[i][j]=0; 
          u[i][j][0]=0; 
          u[i][j][1]=0; 
          

          for(k=0;k<Q;k++) 
          { 
			f[i][j][k]=ff[i][j][k];
            rho[i][j]+=f[i][j][k]; 
            u[i][j][0]+=c*e[k][0]*f[i][j][k]; /////////////////////////////////////////  c!=1
            u[i][j][1]+=c*e[k][1]*f[i][j][k]; 
          } 
          
      
    

          u[i][j][0]=u[i][j][0]/rho[i][j] ; //+0.5*dt*Force/rho[i][j]; //   
          u[i][j][1]/=rho[i][j]; 

       } 

}


void Outdata(int m)
{
	
	FILE  *fm;
	char filename1[50];
	sprintf(filename1,"%s%d%s%d.dat","NX=",NX,"times=",m);

	fm=fopen(filename1,"w");
        
        if(fm==NULL)
           {
             printf("failed to open! \n");
             exit(1);
            }
	fprintf(fm,"TITLE = \"stream\"\n");
	fprintf(fm,"VARIABLES = \"X\",\"Y\",\"U\",\"V\",\"p\"\n");
	fprintf(fm,"ZONE I=%d J=%d  F=POINT\n",NX+1,NY+1);


    for(i=0;i<=NX;i++) 
	{
		rho[i][0]=-q_down*rho[i][2]+(1.0+q_down)*rho[i][1];
		rho[i][NY]=-q_up*rho[i][NY-2]+(1.0+q_up)*rho[i][NY-1];
	}
	
	
	for(j=0;j<=NY;j++) 
              for(i=0;i<=NX;i++) 
			{
				

				fprintf(fm,"%f  %f   %lf  %lf  %lf  \n",xlabel[i][j], ylabel[i][j], u[i][j][0], u[i][j][1], rho[i][j]);
			}

	
       	fclose(fm);
}






 void Error()   //compute error
 { 

          double temp1,temp2,temp3,temp4; 

          temp1=temp2=0; 
		  temp3=temp4=0;
		  
          for(i=0; i<=NX; i++)
          for(j=1;j<NY;j++) 
            { 

             

              temp1+=(  (u[i][j][0]-u_exact[j] )*(u[i][j][0]-u_exact[j])  +  u[i][j][1]*u[i][j][1]  )  ; 
              temp2+=(u_exact[j]*u_exact[j])  ;


			  temp3+=(  (u[i][j][0]-u0[i][j][0] )*(u[i][j][0]-u0[i][j][0])  +  (u[i][j][1]-u0[i][j][1])*(u[i][j][1]-u0[i][j][1])  )  ; 
              temp4+=(u[i][j][0]*u[i][j][0] + u[i][j][1]*u[i][j][1])  ;

			  u0[i][j][0] = u[i][j][0];
			  u0[i][j][1] = u[i][j][1];
            } 


          temp1=sqrt(temp1); 
          temp2=sqrt(temp2); 
          error=temp1/(temp2+1e-30); 


		  temp3=sqrt(temp3); 
          temp4=sqrt(temp4); 
          ERROR=temp3/(temp4+1e-30);

 } 





 
