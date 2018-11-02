///////////////// Taylor-Green vortex flow with curved boundaries, convex scheme


#include<iostream> 
#include<cmath> 
#include<cstdlib> 
#include<iomanip> 
#include<fstream> 
#include<sstream> 
#include<string> 
#include<stdio.h> 
#include "utilities.h"
 
using namespace std; 
const int Q=9; 
const int NX=39; 
const int NY=39; 
const double U=0.05; 
const double pi=3.1415926;


 
int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};    //9 directions
int ne[Q]={0,3,4,1,2,7,8,5,6};                                                //back directions
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36}; 
//double rho[NX+1][NY+1],u[NX+1][NY+1][2],u0[NX+1][NY+1][2],f[NX+1][NY+1][Q],ff[NX+1][NY+1][Q],F[NX+1][NY+1][Q],xlabel[NX+1][NY+1],ylabel[NX+1][NY+1]; 
Array2D<double> rho;
Array3D<double> u;
Array3D<double> u0;
Array3D<double> f;
Array3D<double> ff;
Array3D<double> F;
Array2D<double> xlabel;
Array2D<double> ylabel; 
int i,j,k,ip,jp,n,q_flag; 
double c,Re,dx,dy,Lx,Ly,D,dt,rho0,p0,tau_f,niu,error,y,yy1,yy2,kk,b,cc,x1,x2,x,q; 

double iq,jq,AA,BB,CC,DD,EE,rr,uu,vv,Center_x,Center_y;
double R;         //radius of the circle 
double ell=1.0;  // a parameter to control the shape of the boundary: ell*(x-x0)^2+(y-y0)^2=R^2, ell=1 for circle
double s_nu,s_q,SS,cs_2;

//double abs( double i);
void comput_q (int i, int j, int ip, int jp);




 
void init(); 
double feq(int k,double rho,double u[2]); 
void evolution(); 
void output(int m); 
void Error(); 
//int flag[NX+1][NY+1];
Array2D<int> flag;


 
int main() 
{ 
  using namespace std; 

  rho.init(NX+1,NY+1);
  u.init(NX+1,NY+1,2);
  u0.init(NX+1,NY+1,2);
  f.init(NX+1,NY+1,Q);
  ff.init(NX+1,NY+1,Q);
  F.init(NX+1,NY+1,Q);
  xlabel.init(NX+1,NY+1);
  ylabel.init(NX+1,NY+1); 
  flag.init(NX+1,NY+1);

  init();  //initiate


  for(n=0; ;n++) 
  { 
    evolution(); 
    if(n%100==0) 
    { 
      Error(); 
	}

	if(n%100==0)
	{
      cout<<"The"<<n<<"th computation result:"<<endl<<"The u,v of point (NX/2,NY/2)is : " 
     
      <<setprecision(6)<<u(NX/2,NY/2,0)<<","<<u(NX/2,NY/2,1)<<endl; 
      cout<<"The max relative error of uv is:" 
        <<setiosflags(ios::scientific)<<error<<endl; 
        
	} 

//	if(n%1000==0)  
//		output(n);

	if(n==int(1.0*Lx/U/dt)) 
	{
        Error();
		cout<<"The max relative error of uv is:" 
        <<setiosflags(ios::scientific)<<error<<endl; 

		output(n+1);
		break;
	}

  } 
 
  return 0; 
} 



void init() 
{ 
  
  Lx=1.0; 
  Ly=1.0; 
  dx=Lx/(NX+1); 
  dy=dx; 
  niu=0.002;
  SS=2.0;      //tau
  s_nu=-1.0/SS;
  //s_q = s_nu;
  s_q=-8.0*(2+s_nu)/(8+s_nu);             
  dt=(SS -0.5)/3.0 *dx*dx /niu;
  c=dx/dt;

  cout<<"U/c = "<<U/c<<"\n";




  R=Lx/4.0;  // radius 
  rho0=1.0;
  cs_2=c*c/3.0;
  
  for(i=0;i<=NX;i++)    //corrdinates of the lattice nodes
  for(j=0;j<=NY;j++) 
  {
	  xlabel.assign(i,j)=i*dx+0.5*dx;
	  ylabel.assign(i,j)=j*dy+0.5*dy;

  }
  Center_x=Lx/2.0;  //the centre of the circle
  Center_y=Ly/2.0;

  

  for (i=0;i<=NX;i++)    //tag different lattice nodes 
	  for (j=0;j<=NY;j++)
	  {
		  flag.assign(i,j)=0;           

		  if( (  ell*(xlabel(i,j)-Center_x)*(xlabel(i,j)-Center_x) + (ylabel(i,j)-Center_y)*(ylabel(i,j)-Center_y) ) < R*R ) 
		  {
			 
		        flag.assign(i,j)=1;  //internal fluid lattice nodes	  

		  }

	  }
  
 
  for(i=0;i<=NX;i++)    //initialization of DF
    for(j=0;j<=NY;j++) 
    { 
      u.assign(i,j,0)=-U*cos(2.0*pi*xlabel(i,j))*sin(2.0*pi*ylabel(i,j)); 
      u.assign(i,j,1)=U*cos(2.0*pi*ylabel(i,j))*sin(2.0*pi*xlabel(i,j)); 
      rho.assign(i,j)= rho0 -3.0*U*U/4.0/c/c*( cos(4.0*pi*xlabel(i,j))+cos(4.0*pi*ylabel(i,j)) ); 


      for(k=0;k<Q;k++) 
      { 
        f.assign(i,j,k)=feq(k,rho(i,j),u(i,j)); 
      } 


    }

} 



 
double feq(int k,double rho,double u_0, double u_1)   //  equilibrium distribution
{ 
  double eu,uv,feq; 
  eu=(e[k][0]*u_0+e[k][1]*u_1); 
  uv=(u_0*u_0+u_1*u_1); 
  feq=w[k]*(rho + rho0* (3.0*eu/c+4.5*eu*eu/c/c-1.5*uv/c/c) ); 
  return feq; 
} 


 
void evolution() 
{ 
  
	 for(i=(NX+1)/4-3;i<=(NX+1)/4*3+3;i++) 
		for(j=(NY+1)/4-3;j<=(NY+1)/4*3+3;j++) 
		{
			for(k=0;k<Q;k++) 
				F.assign(i,j,k)=f(i,j,k)
				           +
				           s_nu*(  0.5*(f(i,j,k)+f(i,j,ne[k])) - 0.5*(feq(k,rho(i,j), u(i,j,0), u(i,j,1))+feq(ne[k],rho(i,j), u(i,j,0), u(i,j,1)))   )
						   +
				           s_q*(  0.5*(f(i,j,k)-f(i,j,ne[k])) - 0.5*(feq(k,rho(i,j), u(i,j,0), u(i,j,1))-feq(ne[k],rho(i,j), u(i,j,0), u(i,j,1)))   );    //collision
		

		}
	
	
		
for(i=(NX+1)/4-3;  i<=(NX+1)/4*3+3;  i++) 
		for(j=(NY+1)/4-3;   j<=(NY+1)/4*3+3;  j++) 
			{

			if(flag(i,j)==1)
			{
				for(k=0;k<Q;k++) 
				{
					ip=i-e[k][0]; 
					jp=j-e[k][1];
				

					if( flag(i,j)==1 && flag(ip,jp)==0 )
					{
					   
					   comput_q(i,j,ip,jp); // compute the ratio $gamma$, the computed value is 'q', which is equal to gamma * dx

					   iq=xlabel(i,j)-q*double(e[k][0]);
					   jq=ylabel(i,j)-q*double(e[k][1]);

					    
					    rr= rho0 - 3.0*U*U/4.0/c/c*( cos(iq*4.0*pi)+cos(jq*4.0*pi) )*exp(-16.0*niu*pi*pi*(n)*dt) ;

						uu=-U*cos(iq*2.0*pi)*sin(jq*2.0*pi)*exp(-8.0*niu*pi*pi*(n)*dt) ; 

                        vv=U*sin(iq*2.0*pi)*cos(jq*2.0*pi)*exp(-8.0*niu*pi*pi*(n)*dt)  ; 
					

						q=q/dx;

					
                        AA= 2.0*q;
						BB = 1.0-AA;
						CC = 2.0;

						ff.assign(i,j,k) = AA*F(i,j,ne[k]) +BB*f(i,j,ne[k]) + CC*w[k]*rho0*3.0/c*(e[k][0]*uu+e[k][1]*vv);

					}
					
					
					else
					
					ff.assign(i,j,k)=F(ip,jp,k);			
				}

			}
			
			}

    

     for(i=(NX+1)/4-3;i<=(NX+1)/4*3+3;i++) 
		for(j=(NY+1)/4-3;j<=(NY+1)/4*3+3;j++) 
		if(flag(i,j)==1)
        { 
			
          rho.assign(i,j)=0; 
          u.assign(i,j,0)=0; 
          u.assign(i,j,1)=0; 
		
          for(k=0;k<Q;k++) 
          { 
            f.assign(i,j,k) = ff(i,j,k);
            rho.assign(i,j)+=f(i,j,k); 
            u.assign(i,j,0)+=c*e[k][0]*f(i,j,k); 
            u.assign(i,j,1)+=c*e[k][1]*f(i,j,k); 
          } 


		  u.assign(i,j,0)/=rho0; 
          u.assign(i,j,1)/=rho0; 

        
		} 

}





//double abs( double i)
//{
//	if(i>=0.0) return i;
//	else 
//		return -i;
//}



void comput_q (int i, int j, int ip, int jp)  // compute the ratio $gamma$, the computed value is 'q', which is equal to gamma * dx
{
     if (ip==i)
	 {   
		 yy1  = fabs( Center_y+sqrt( (R*R-(xlabel(i,j)-Center_x)*(xlabel(i,j)-Center_x))/ell )-ylabel(i,j) );
		 yy2  = fabs( Center_y-sqrt( (R*R-(xlabel(i,j)-Center_x)*(xlabel(i,j)-Center_x))/ell )-ylabel(i,j) );

		 if(yy1<=yy2) q=yy1;
		 else q=yy2;
	 }


	 else
	 {
         kk  =  (ylabel(ip,jp)-ylabel(i,j))/(xlabel(ip,jp)-xlabel(i,j));
		 
		 b   =  (  2.0*Center_x - 2.0*ell*kk*(ylabel(i,j)-kk*xlabel(i,j)-Center_y)  ) / (ell*kk*kk+1.0);

		 cc  =  ( ell*(ylabel(i,j)-kk*xlabel(i,j)-Center_y)*(ylabel(i,j)-kk*xlabel(i,j)-Center_y)+Center_x*Center_x-R*R  ) / (ell*kk*kk+1.0);  

         x1  =  fabs( (b+sqrt(b*b-4.0*cc))/2.0-xlabel(i,j) );
         x2  =  fabs( (b-sqrt(b*b-4.0*cc))/2.0-xlabel(i,j) );

		 if(x1<=x2) q=x1;
		 else q=x2;
		 
	 }
}




 void output(int m)    //output the data
 { 
          ostringstream name; 
          name<<"TaylorGreen"<<"Mesh"<<NX<<"_"<<NY<<"compute_"<<m<<"times"".dat"; 
          ofstream out(name.str().c_str()); 
         out<<"Title= \"TaylorGreen\"\n"<<"VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"U0\",\"V0\",\"p\" \n"<<"ZONE T=\"BOX\",I=" 
            <<NX+1<<",J="<<NY+1<<",F=POINT"<<endl;
            for(j=0;j<=NY;j++) 
              for(i=0;i<=NX;i++) 
              { 
                out<<setprecision(15)<<xlabel(i,j)<<" "<<ylabel(i,j)<<" "<<u(i,j,0)<<" "<<  u(i,j,1)<<" "<<u0(i,j,0)<<" "<<  u0(i,j,1)<<" "<<rho(i,j)<<endl; 
			//	cout<<rho[i][j]<<endl;
              } 
 } 
 
      
		
void Error()   //compute error
{ 
          double temp1,temp2; 
          temp1=0; 
          temp2=0;
		  
          for(i=1;i<NX;i++) 
            for(j=1;j<NY;j++)
			{
			   if(flag.assign(i,j)==1)
			   {
				 
				
                  u0.assign(i,j,0)=-U*cos(xlabel(i,j)*2.0*pi)*sin(ylabel(i,j)*2.0*pi)*exp(-8.0*niu*pi*pi*(n)*dt);
				  u0.assign(i,j,1)= U*cos(ylabel(i,j)*2.0*pi)*sin(xlabel(i,j)*2.0*pi)*exp(-8.0*niu*pi*pi*(n)*dt);


                  temp1+=( (u(i,j,0)-u0(i,j,0))*(u(i,j,0)-u0(i,j,0))+(u(i,j,1)-u0(i,j,1))*(u(i,j,1)-u0(i,j,1))); 
                  temp2+=(u0(i,j,0)*u0(i,j,0)+u0(i,j,1)*u0(i,j,1)); 
			      		
				
			   }
			}


            temp1=sqrt(temp1); 
            temp2=sqrt(temp2); 
            error=temp1/(temp2+1e-30); 
 } 
 
