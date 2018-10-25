// Taylor Green vortex flow, convex scheme


#include<iostream> 
#include<cmath> 
#include<cstdlib> 
#include<iomanip> 
#include<fstream> 
#include<sstream> 
#include<string> 
#include<stdio.h> 
#include "ArgParseStandalone.h"
#include "utilities.h"
 
using namespace std; 
const int Q=9; 
//const int NX = 49; 
//const int NY = 49; 
int NY=0; 
int NX=0; 
const double U=0.05; 
const double pi=3.1415926;


 
int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};       //9个方向
int ne[Q]={0,3,4,1,2,7,8,5,6};                                                //用于标记反方向
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36}; 
double** rho = NULL;
double*** u = NULL;
double*** u0 = NULL;
double*** f = NULL;
double*** ff = NULL;
double*** F = NULL;
double** xlabel = NULL;
double** ylabel = NULL;
int i,j,k,ip,jp,n,q_flag; 
double c,Re,dx,dy,Lx,Ly,D,dt,rho0,p0,tau_f,niu,error,y,yy1,yy2,kk,b,cc,x1,x2,x,q; 

double iq,jq,AA,BB,CC,DD,EE,rr,uu,vv,Center_x,Center_y;
double R;  //圆半径 
double ell=1.0;  //椭圆参数，即ell*(x-x0)^2+(y-y0)^2=R^2 , ell控制椭圆的扁平
double s_nu,s_q,SS,cs_2;

//double abs( double i);
void comput_q (int i, int j, int ip, int jp);




 
void init(double tau); 
double feq(int k,double rho,double u[2]); 
void evolution(); 
void output(int m); 
void Error(); 
int** flag = NULL;
bool verbose = false;


// void outdata();
 
int main(int argc, char** argv) 
{ 
  using namespace std; 

  double tau = 0.;
  int Ny = 0;
  bool dump_solution_passed = false;
  std::string solution_filepath;
  bool header = false;

  ArgParse::ArgParser Parser("Poiseuille Non Convex Simulation");
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

  NY = Ny;
  NX = NY;

  rho = New2DArray<double>(NX+1,NY+1); 
  u = New3DArray<double>(NX+1,NY+1,2);
  u0 = New3DArray<double>(NX+1,NY+1,2);
  f = New3DArray<double>(NX+1,NY+1,Q);
  ff = New3DArray<double>(NX+1,NY+1,Q);
  F = New3DArray<double>(NX+1,NY+1,Q);
  xlabel = New2DArray<double>(NX+1,NY+1);
  ylabel = New2DArray<double>(NX+1,NY+1); 
  flag = New2DArray<int>(NX+1,NY+1);
  

  init(tau);  //初始化


  for(n=0; ;n++) 
  { 
    evolution(); 
    if(n%100==0) 
    { 
      Error(); 
	}

	if(n%100==0)
	{
		if(verbose) {
      cout<<"The"<<n<<"th computation result:"<<endl<<"The u,v of point (NX/2,NY/2)is : " 
     
      <<setprecision(6)<<u[NX/2][NY/2][0]<<","<<u[NX/2][NY/2][1]<<endl; 
      cout<<"The max relative error of uv is:" 
        <<setiosflags(ios::scientific)<<error<<endl; 
		}
       
	//  output(n);
        
	} 

	if(n==int(1.0*Lx/U/dt)) 
	{
        Error();
        if (verbose) {
		cout<<"The max relative error of uv is:" 
        <<setiosflags(ios::scientific)<<error<<endl; 
	}
	if(dump_solution_passed) {
		output(n+1);
	}
        if(header) {
          printf("\"Lattice Size\", \"NY\", \"Tau\", \"Error\"\n");
       	}
       	printf("%.14f, %i, %.14f, %.14f\n",dx, NY,tau,error);
		break;
	}

  } 

  Delete2DArray<double>(rho,NX+1,NY+1); 
  Delete3DArray<double>(u,NX+1,NY+1,2);
  Delete3DArray<double>(u0,NX+1,NY+1,2);
  Delete3DArray<double>(f,NX+1,NY+1,Q);
  Delete3DArray<double>(ff,NX+1,NY+1,Q);
  Delete3DArray<double>(F,NX+1,NY+1,Q);
  Delete2DArray<double>(xlabel,NX+1,NY+1);
  Delete2DArray<double>(ylabel,NX+1,NY+1); 
  Delete2DArray<int>(flag,NX+1,NY+1);

  return 0; 
} 



void init(double tau) 
{ 
  
  Lx=1.0; 
  Ly=1.0; 
  dx=Lx/(NX+1); 
  dy=dx; 
  niu=0.002;
  SS=tau;                    //tau
  s_nu=-1.0/SS;
  //s_q = s_nu;
  //s_q=8.0*(2+s_nu)/(8+7.0*s_nu); 
  s_q=8.0*(2+s_nu)/(8+s_nu); 
  dt=(SS -0.5)/3.0 *dx*dx /niu;
  c=dx/dt;

  if(verbose) {
  cout<<"U/c = "<<U/c<<"\n";
  }




  R=Lx/4.0;  // 圆半径  
  rho0=1.0;
  cs_2=c*c/3.0;
  
  for(i=0;i<=NX;i++)    //格点坐标
  for(j=0;j<=NY;j++) 
  {
	  xlabel[i][j]=i*dx+0.5*dx;
	  ylabel[i][j]=j*dy+0.5*dy;

  }
  Center_x=Lx/2.0;
  Center_y=Ly/2.0;

  

  for (i=0;i<=NX;i++)    //对flag进行赋值，将格点分类
	  for (j=0;j<=NY;j++)
	  {
		  flag[i][j]=0;           //内部流体节点

		  if( (  ell*(xlabel[i][j]-Center_x)*(xlabel[i][j]-Center_x) + (ylabel[i][j]-Center_y)*(ylabel[i][j]-Center_y) ) < R*R )//圆内部
		  {
			 
		        flag[i][j]=1;	  

		  }

	  }
  
 
  for(i=0;i<=NX;i++)    //速度初始化
    for(j=0;j<=NY;j++) 
    { 
      u[i][j][0]=-U*cos(2.0*pi*xlabel[i][j])*sin(2.0*pi*ylabel[i][j]); 
      u[i][j][1]=U*cos(2.0*pi*ylabel[i][j])*sin(2.0*pi*xlabel[i][j]); 
      rho[i][j]= rho0 -3.0*U*U/4.0/c/c*( cos(4.0*pi*xlabel[i][j])+cos(4.0*pi*ylabel[i][j]) ); 


      for(k=0;k<Q;k++) 
      { 
        f[i][j][k]=feq(k,rho[i][j],u[i][j]); 
      } 


    }


} //结束初始化



 
double feq(int k,double rho,double u[2])   // 计算平衡态分布函数
{ 
  double eu,uv,feq; 
  eu=(e[k][0]*u[0]+e[k][1]*u[1]); 
  uv=(u[0]*u[0]+u[1]*u[1]); 
  //feq=w[k]*rho*(1.0+3.0*eu/c+4.5*eu*eu/c/c-1.5*uv/c/c); 
  feq=w[k]*(rho + rho0* (3.0*eu/c+4.5*eu*eu/c/c-1.5*uv/c/c) ); 
  if(isinf(feq)) {
	  std::cout << w[k] << " " << eu << " " << uv << " " << u[0] << " " << u[1] << std::endl;
  }
  return feq; 
} 


 
void evolution() 
{ 
	std::cout << "evolution" << std::endl;


  
	 for(i=(NX+1)/4-3;i<=(NX+1)/4*3+3;i++) 
		for(j=(NY+1)/4-3;j<=(NY+1)/4*3+3;j++) 
		{
			for(k=0;k<Q;k++)  {
				bool thenan_before = false;
				bool thenan_after = false;
				if(isnan(F[i][j][k])) {
					thenan_before = true;
				}
				F[i][j][k]=f[i][j][k]
				           +
				           s_nu*(  0.5*(f[i][j][k]+f[i][j][ne[k]]) - 0.5*(feq(k,rho[i][j], u[i][j])+feq(ne[k],rho[i][j], u[i][j]))   )
						   +
				           s_q*(  0.5*(f[i][j][k]-f[i][j][ne[k]]) - 0.5*(feq(k,rho[i][j], u[i][j])-feq(ne[k],rho[i][j], u[i][j]))   );    //碰撞
			        if(isnan(F[i][j][k])) {
				        thenan_after = true;
				}	
			        if((!thenan_before)&&(thenan_after)) {
				        if(isnan(f[i][j][k])) {
						std::cout << "condition met: first summation term" << std::endl;
					}
					if(isnan(s_nu*(  0.5*(f[i][j][k]+f[i][j][ne[k]]) - 0.5*(feq(k,rho[i][j], u[i][j])+feq(ne[k],rho[i][j], u[i][j]))))) {
						std::cout << "condition met: second summation term" << std::endl;
					}
					if(isnan(s_q*(  0.5*(f[i][j][k]-f[i][j][ne[k]]) - 0.5*(feq(k,rho[i][j], u[i][j])-feq(ne[k],rho[i][j], u[i][j]))   ))) {
						std::cout << "condition met: third summation term" << std::endl;
						std::cout << "printing factors" << std::endl;
						std::cout << s_q << " " << f[i][j][k] << " " << f[i][j][ne[k]] << " " << feq(k, rho[i][j], u[i][j]) << " " << feq(ne[k], rho[i][j], u[i][j]) << std::endl;
					}
				}
			}

		}
	

	
	
		
for(i=(NX+1)/4-3;  i<=(NX+1)/4*3+3;  i++) 
		for(j=(NY+1)/4-3;   j<=(NY+1)/4*3+3;  j++) 
			{
			if(flag[i][j]==1)
			{
				for(k=0;k<Q;k++) 
				{
					ip=i-e[k][0]; 
					jp=j-e[k][1];
				
					if( flag[i][j]==1 && flag[ip][jp]==0 )
					{
					   
						comput_q(i,j,ip,jp);

					   iq=xlabel[i][j]-q*double(e[k][0]);
					   jq=ylabel[i][j]-q*double(e[k][1]);

					    
					    rr= rho0 - 3.0*U*U/4.0/c/c*( cos(iq*4.0*pi)+cos(jq*4.0*pi) )*exp(-16.0*niu*pi*pi*(n)*dt) ;

						uu=-U*cos(iq*2.0*pi)*sin(jq*2.0*pi)*exp(-8.0*niu*pi*pi*(n)*dt) ; 

                        vv=U*sin(iq*2.0*pi)*cos(jq*2.0*pi)*exp(-8.0*niu*pi*pi*(n)*dt)  ; 
					

						q=q/dx;


                        AA= 2.0*q/(1.0+2.0*q);
						BB = 1.0-AA;
						CC = 1.0-AA+BB;

						ff[i][j][k] = AA*F[i][j][k] +BB*f[i][j][ne[k]] + CC*w[k]*rho0*3.0/c*(e[k][0]*uu+e[k][1]*vv);

					}
					
					
					else
					
					ff[i][j][k]=F[ip][jp][k];			
				}

			}
			
			}

    

     for(i=(NX+1)/4-3;i<=(NX+1)/4*3+3;i++) 
		for(j=(NY+1)/4-3;j<=(NY+1)/4*3+3;j++) 
		if(flag[i][j]==1)
        { 
			
          rho[i][j]=0; 
          u[i][j][0]=0; 
          u[i][j][1]=0; 
		
          for(k=0;k<Q;k++) 
          { 
            f[i][j][k] = ff[i][j][k];
            rho[i][j]+=f[i][j][k]; 
            if(isnan(e[k][0])) {
		    std::cout << "e[" << k << "][0] is nan!" << std::endl;
	    }
            if(isnan(e[k][1])) {
		    std::cout << "e[" << k << "][1] is nan!" << std::endl;
	    }
            if(isnan(f[i][j][k])) {
		    std::cout << "f[" << i << "][" << j << "][" << k << "] is nan!" << std::endl;
	    }
            u[i][j][0]+=c*e[k][0]*f[i][j][k]; 
            u[i][j][1]+=c*e[k][1]*f[i][j][k]; 
          } 

	if(isnan(u[i][j][0])) {
		std::cout << "u[" << i << "][" << j << "][0] is nan!" << std::endl;
	}
	if(isnan(u[i][j][1])) {
		std::cout << "u[" << i << "][" << j << "][1] is nan!" << std::endl;
	}
          u[i][j][0]/=rho[i][j]; 
          u[i][j][1]/=rho[i][j]; 

		  //u[i][j][0]/=rho0; 
          //u[i][j][1]/=rho0; 

        
		} 

}





//double abs( double i)
//{
//	if(i>=0.0) return i;
//	else 
//		return -i;
//}



void comput_q (int i, int j, int ip, int jp)  //1表示大于0.5，-1表示小于0.5,2为正常格点
{
     if (ip==i)
	 {   
		 yy1  = fabs( Center_y+sqrt( (R*R-(xlabel[i][j]-Center_x)*(xlabel[i][j]-Center_x))/ell )-ylabel[i][j] );
		 yy2  = fabs( Center_y-sqrt( (R*R-(xlabel[i][j]-Center_x)*(xlabel[i][j]-Center_x))/ell )-ylabel[i][j] );

		 if(yy1<=yy2) q=yy1;
		 else q=yy2;
	 }


	 else
	 {
         kk  =  (ylabel[ip][jp]-ylabel[i][j])/(xlabel[ip][jp]-xlabel[i][j]);
		 
		 b   =  (  2.0*Center_x - 2.0*ell*kk*(ylabel[i][j]-kk*xlabel[i][j]-Center_y)  ) / (ell*kk*kk+1.0);

		 cc  =  ( ell*(ylabel[i][j]-kk*xlabel[i][j]-Center_y)*(ylabel[i][j]-kk*xlabel[i][j]-Center_y)+Center_x*Center_x-R*R  ) / (ell*kk*kk+1.0);  

         x1  =  fabs( (b+sqrt(b*b-4.0*cc))/2.0-xlabel[i][j] );
         x2  =  fabs( (b-sqrt(b*b-4.0*cc))/2.0-xlabel[i][j] );

		 if(x1<=x2) q=x1;
		 else q=x2;
		 
	 }
}




        void output(int m)    //输出
        { 
          ostringstream name; 
          name<<"TaylorGreen"<<"Re_"<<Re<<"网格"<<NX<<"_"<<NY<<"迭代_"<<m<<"次"".dat"; 
          ofstream out(name.str().c_str()); 
         out<<"Title= \"TaylorGreen\"\n"<<"VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"U0\",\"V0\",\"p\" \n"<<"ZONE T=\"BOX\",I=" 
            <<NX+1<<",J="<<NY+1<<",F=POINT"<<endl;
            for(j=0;j<=NY;j++) 
              for(i=0;i<=NX;i++) 
              { 
                out<<setprecision(15)<<xlabel[i][j]<<" "<<ylabel[i][j]<<" "<<u[i][j][0]<<" "<<  u[i][j][1]<<" "<<u0[i][j][0]<<" "<<  u0[i][j][1]<<" "<<rho[i][j]<<endl; 
			//	cout<<rho[i][j]<<endl;
              } 

			  
        } 
 
      
		
		void Error()   //计算相对误差
        { 
          double temp1,temp2; 
          temp1=0; 
          temp2=0;
		  
          for(i=1;i<NX;i++) 
            for(j=1;j<NY;j++)
			{
			   if(flag[i][j]==1)
			   {
				 
				
                  u0[i][j][0]=-U*cos(xlabel[i][j]*2.0*pi)*sin(ylabel[i][j]*2.0*pi)*exp(-8.0*niu*pi*pi*(n)*dt);
				  u0[i][j][1]= U*cos(ylabel[i][j]*2.0*pi)*sin(xlabel[i][j]*2.0*pi)*exp(-8.0*niu*pi*pi*(n)*dt);


                  temp1+=( (u[i][j][0]-u0[i][j][0])*(u[i][j][0]-u0[i][j][0])+(u[i][j][1]-u0[i][j][1])*(u[i][j][1]-u0[i][j][1])); 
                  temp2+=(u0[i][j][0]*u0[i][j][0]+u0[i][j][1]*u0[i][j][1]); 
			      
		
				
				
			   }
			}


	  std::cout<< "NX: " << NX << std::endl;
	  std::cout << "NY: " << NY << std::endl;
	  std::cout << "temp1: " << temp1 << std::endl;
	  std::cout << "temp2: " << temp2 << std::endl;
            temp1=sqrt(temp1); 
            temp2=sqrt(temp2); 
            error=temp1/(temp2+1e-30); 
        } 
 
