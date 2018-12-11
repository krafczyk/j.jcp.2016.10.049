////////////////   poiseuille , convex scheme stability test

#include<math.h> 
#include<cstdlib> 
//#include<iomanip.h> 

#include<string.h> 
#include<stdio.h> 
#include <string>
#include "ArgParseStandalone.h"
#include "utilities.h"
#include <time.h>
 

const int Q=9; 

int NY=0; 
int NX=0; 

 
const int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};    // 9 direction
const int ne[Q]={0,3,4,1,2,7,8,5,6};                                                // back direction
const double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36}; 

Array2D<double> rho;
Array3D<double> u;
Array3D<double> f;
Array3D<double> F;
Array3D<double> ff;
Array2D<double> xlabel;
Array2D<double> ylabel;
double* u_exact = NULL;
Array3D<double> u0;
Array2D<int> flag;

int i,j,k,ip,jp,n; 
double c,Re,dx,dy,Lx,Ly,dt,rho0,P0,s_nu,s_q,SS,niu,error,Force,U_max,U_s,cs_2,rela_error,ERROR; 
double a1,a2,a3,a4,a5,q_up,q_down;
 
void init(double gamma, double tau); 
double feq(int k,double rho,double u_0, double u_1); 
void evolution(); 
void output(int m); 
void Error(); 
void Outdata(int m);
bool verbose = false;

const double tau_low = 0.5001;
const double tau_high = 3.0;

int main(int argc, char** argv) 
{ 
  double gamma = 0.;
  double tau = 0.;
  int Ny = 0;
  bool dump_solution_passed = false;
  std::string solution_filepath;
  bool header = false;
  time_t time_limit = 10*60;

  ArgParse::ArgParser Parser("Poiseuille Non Convex Simulation");
  Parser.AddArgument("--gamma", "Set the value for gamma", &gamma, ArgParse::Argument::Required);
  Parser.AddArgument("--tau", "Set the value for tau", &tau, ArgParse::Argument::Required);
  Parser.AddArgument("--Ny", "Set the y resolution Ny", &Ny, ArgParse::Argument::Required);
  Parser.AddArgument("--dump-solution", "Filepath to dump solution at.", &solution_filepath, ArgParse::Argument::Optional, &dump_solution_passed);
  Parser.AddArgument("--header", "Whether or not to include column headers in the output", &header, ArgParse::Argument::Optional);
  Parser.AddArgument("--verbose", "Whether to print extra stuff", &verbose, ArgParse::Argument::Optional);

  if(Parser.ParseArgs(argc, argv) < 0) {
	  printf("Problem parsing arguments!");
	  return -1;
  }

  if(Parser.HelpPrinted()) {
	  return 0;
  }

  // Set the resolution
  NY=Ny; 
  NX=2*(NY-1); 

  rho.init(NX+1,NY+1);
  u.init(NX+1,NY+1,2);
  f.init(NX+1,NY+1,Q);
  F.init(NX+1,NY+1,Q);
  ff.init(NX+1,NY+1,Q);
  xlabel.init(NX+1,NY+1);
  ylabel.init(NX+1,NY+1);
  u_exact = new double[NY+1]; 
  u0.init(NX+1,NY+1,2);
  flag.init(NX+1,NY+1);

  for (i=0;i<=NX;i++)    //tag different boundary
	  for (j=0;j<=NY;j++)
	  {
		  if      (j==0)   flag.assign(i,j)=1;  //Up and Down Boundary
		  else if (j==NY)  flag.assign(i,j)=5;  //Up and Down Boundary

		  else if (i==0)   flag.assign(i,j)=2;
          else if (i==NX)  flag.assign(i,j)=3; 
		  else  		   flag.assign(i,j)=4;

	  }

         //printf("flag [0][NY]=%d\n",flag[0][NY]);
         //printf("flag [NX][NY]=%d\n",flag[NX][NY]);
         //printf("flag [0][0]=%d\n",flag[0][0]);
         //printf("flag [NX][0]=%d\n",flag[NX][0]);
         //printf("flag [0][3]=%d\n",flag[0][3]);
         //printf("flag [NX][3]=%d\n",flag[NX][3]);


  init(gamma, tau);  // initate
  if(verbose) {
  printf("U_max/c=%f\n",U_max/c);
  }

  time_t init_time = time(NULL);
 for(n=0;;n++) 
  { 

    evolution(); 

    if(n%1000==0) 
    { 
      Error(); 
      if((time(NULL)-init_time) > time_limit) {
      	printf("%.14f, %i, %.14f, %.14f, False\n",dx, NY,gamma,tau);
      	break;
      }
     
      if(verbose) {
      printf("The %d th computationa result:\n",n);
      
      printf("The relative error of 1000 steps is: %.13f ,  error=%.10f  \n",ERROR, error);
      }
      
    }
      
    if(n>0 && fabs(ERROR)<1.0e-10) 
     { 
       if(dump_solution_passed) {
	   Outdata(n);
       }
       if(header) {
         printf("\"Lattice Size\", \"NY\", \"Gamma\", \"Tau\", \"Stable\"\n");
       }
       printf("%.14f, %i, %.14f, %.14f, True\n",dx, NY,gamma,tau);

       break; 
     } else if (fabs(ERROR) > 10) {
       printf("%.14f, %i, %.14f, %.14f, False\n",dx, NY,gamma,tau);
       break;
     }
  } 

  delete[] u_exact; 

  return 0;      
} 



void init(double gamma, double tau) 
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
  q_up=gamma;     //ÉÏ±ß½çµÄÍø¸ñ²½³¤
  q_down=q_up;
  SS=tau;     //tau



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



  if(verbose) {
  printf("s_nu=%f, s_q=%f \n",s_nu,s_q);

  printf("U/c=%f",U_max/c);
  }

  //cout<<"tau_f= "<<tau_f<<endl; 

  

  for(i=0;i<=NX;i++)
	  for(j=0;j<=NY;j++)
	  {

		   xlabel.assign(i,j)=i*dx;

		   ylabel.assign(i,j)=(j-1)*dy+q_down*dy;

		   ylabel.assign(i,0)=0.0;
		   ylabel.assign(i,NY) = Ly;

	  }

   for(j=0;j<=NY;j++)
   {
	   u_exact[j]=4.0*U_max*(1.0-ylabel(0,j)/Ly)*ylabel(0,j)/Ly;  //anlytical solution
   }
 


 

  for(i=0;i<=NX;i++)    //  initialization of DF
    for(j=0;j<=NY;j++) 
    { 
      u.assign(i,j,0)=0; 
      u.assign(i,j,1)=0; 
   
      rho.assign(i,j)=1.0;
     
      for(k=0;k<Q;k++) 
       { 
         f.assign(i,j,k)=feq(k,rho(i,j),u(i,j,0), u(i,j,1)); 
       }
 
    } 
} 

 
double feq(int k,double rho,double u_0, double u_1)   //  equilibrium distribution
{ 
  double eu,uv,feq; 
  eu=(e[k][0]*u_0+e[k][1]*u_1); 
  uv=(u_0*u_0+u_1*u_1); 
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
				F.assign(i,j,k)=f(i,j,k)
				           +
				           s_nu*(  0.5*(f(i,j,k)+f(i,j,ne[k])) - 0.5*(feq(k,rho(i,j), u(i,j,0), u(i,j,1))+feq(ne[k],rho(i,j), u(i,j,0), u(i,j,1)))   )
						   +
				           s_q*(  0.5*(f(i,j,k)-f(i,j,ne[k])) - 0.5*(feq(k,rho(i,j), u(i,j,0), u(i,j,1))-feq(ne[k],rho(i,j), u(i,j,0), u(i,j,1)))   )
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
					
					if(flag(ip,jp)==4||flag(ip,jp)==2||flag(ip,jp)==3)    //inner fluid
						ff.assign(i,j,k)=F(ip,jp,k);



					
					
					else if(flag(ip,jp)==5)

					{
					    a1=2.0*q_up/(1.0+2.0*q_up);
						a2=1.0-a1;
						ff.assign(i,j,k) = a1 * ( F(i,j,k)-dt*w[k]*Force*c*e[k][0]/cs_2)  + a2*f(i,j,ne[k]);
                    }


					else if(flag(ip,jp)==1)

					{
						a1=2.0*q_down/(1.0+2.0*q_down);
						a2=1.0-a1;
						ff.assign(i,j,k) = a1 * ( F(i,j,k)-dt*w[k]*Force*c*e[k][0]/cs_2)  + a2*f(i,j,ne[k]);

                    }
                                               

                              
				}
		}
   



   for(i=0;i<=NX;i++) 
        for(j=1;j<NY;j++) 
        { 

          rho.assign(i,j)=0; 
          u.assign(i,j,0)=0; 
          u.assign(i,j,1)=0; 
          

          for(k=0;k<Q;k++) 
          { 
			f.assign(i,j,k)=ff(i,j,k);
            rho.assign(i,j)+=f(i,j,k); 
            u.assign(i,j,0)+=c*e[k][0]*f(i,j,k); /////////////////////////////////////////  c!=1
            u.assign(i,j,1)+=c*e[k][1]*f(i,j,k); 
          } 
          
      
    

          u.assign(i,j,0)=u(i,j,0)/rho(i,j) ; //+0.5*dt*Force/rho[i][j]; //   
          u.assign(i,j,1)/=rho(i,j); 

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
		rho.assign(i,0)=-q_down*rho(i,2)+(1.0+q_down)*rho(i,1);
		rho.assign(i,NY)=-q_up*rho(i,NY-2)+(1.0+q_up)*rho(i,NY-1);
	}
	
	
	for(j=0;j<=NY;j++) 
              for(i=0;i<=NX;i++) 
			{
				

				fprintf(fm,"%f  %f   %lf  %lf  %lf  \n",xlabel(i,j), ylabel(i,j), u(i,j,0), u(i,j,1), rho(i,j));
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

             

              temp1+=(  (u(i,j,0)-u_exact[j] )*(u(i,j,0)-u_exact[j])  +  u(i,j,1)*u(i,j,1)  )  ; 
              temp2+=(u_exact[j]*u_exact[j])  ;


			  temp3+=(  (u(i,j,0)-u0(i,j,0) )*(u(i,j,0)-u0(i,j,0))  +  (u(i,j,1)-u0(i,j,1))*(u(i,j,1)-u0(i,j,1))  )  ; 
              temp4+=(u(i,j,0)*u(i,j,0) + u(i,j,1)*u(i,j,1))  ;

			  u0.assign(i,j,0) = u(i,j,0);
			  u0.assign(i,j,1) = u(i,j,1);
            } 


          temp1=sqrt(temp1); 
          temp2=sqrt(temp2); 
          error=temp1/(temp2+1e-30); 


		  temp3=sqrt(temp3); 
          temp4=sqrt(temp4); 
          ERROR=temp3/(temp4+1e-30);

 } 




 
