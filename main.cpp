#include "eigen-git-mirror/Eigen/Dense"
#include <vector>
#include <random>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <ostream>
#include <ctime>
#include <cstdlib> 
#include <time.h> 
#include <math.h>
//#include "../matplotlib-cpp/matplotlibcpp.h"
#define PREC 1E-7;
//namespace plt =matplotlibcpp;


double Vnn=40.0; //meV
double Vp=0.0;//meV
double V0=-2.0*Vnn;//meV //negative or positive?
//Choose Starting microstate
Eigen::MatrixXd get_starting_microstate()  //In reality, this would need to be Xf
{
	//Initialize Matrix with all ones
	Eigen::MatrixXd starting_microstate = (-1.0)*Eigen::MatrixXd::Ones(80, 80); 
	return starting_microstate;
}


//member functions for site class
//get random x_index based off of Nx
int x_Index(int Nx) {
	double random_number=rand()%Nx;
	return floor(random_number);
}

//get random y_index based off of Ny
int y_Index(int Ny) {
	double random_number=rand()%Ny;
	return floor(random_number);
}


int make_index(int index, int N)
{
	int make_index=(N+(index%N))%N;
	return make_index;
}

//double random_num()
//{
//	std::random_device my_rand;
//        std::mt19937_64 gen(my_rand());
//	std::uniform_real_distribution<> dis(0,1);
//	return dis(gen);	
//}
//Flip site on previous microstate
Eigen::MatrixXd flip_site(int Nx, int Ny, Eigen::MatrixXd& previous_microstate, const int& x_index, const int& y_index) //change to MatrixXd
{
	previous_microstate(x_index, y_index)=-previous_microstate(x_index, y_index);
	return previous_microstate;
}


//get N
double get_N(Eigen::MatrixXd microstate, int Nx, int Ny)
{	
	auto config=microstate;
	auto act_config=config.array()+1.0;
	return act_config.sum()/2.0;
}
//Calculate total grand canonical energy
double total_energy(int Nx, int Ny, Eigen::MatrixXd& microstate, const double& mu, const double& T)
{
	double k=0.08617333262; //meV/K 
	double chem_pot;
	double energy= V0+microstate.sum()*Vp;
	for (int i=0; i<Nx; i++)
	{
		for (int j=0; j<Ny; j++)
		{
		//	double temp_eng=Vp*microstate(x_index, y_index)+Vnn*(microstate(x_index,y_index))*(((microstate((x_index+1)%Nx, y_index)+microstate(x_index, (y_index+1)%Ny))));
			//double temp_eng=Vp*microstate(i, j)+Vnn*(microstate(i,j))*(((microstate((i+1)%Nx, j)+microstate(i, (j+1)%Ny))));
			double temp_eng= Vnn*(microstate(i,j))*(((microstate((i+1)%Nx, j)+microstate(i, (j+1)%Ny))));
			energy+=temp_eng;
		}
			//double energy= Vp*microstate(x_index,y_index)+Vnn*(microstate(x_index, y_index))*(((microstate((x_index+1)%Nx), y_index))+microstate(x_index, (y_index+1)%Ny)+microstate((x_index-1)%Nx, y_index)+microstate(x_index, (y_index-1)%Ny));
	}
	chem_pot=get_N(microstate, Nx, Ny)*mu;
	double tot_energy=(energy-chem_pot);//   /k/T;
	return tot_energy;
}

//Calculate change in grand canonical energy
double change_in_energy(int Nx, int Ny, Eigen::MatrixXd& current_microstate, const int& x_index, const int& y_index)
{
 	int index_xf =make_index(x_index+1, Nx);
	int index_xb =make_index(x_index-1, Nx);
	int index_yf =make_index(y_index+1, Ny);
	int index_yb =make_index(y_index-1, Ny);

	Eigen::MatrixXd potential_site=current_microstate;
	//might have wrong equation
	//double energy_change=-2*potential_site(x_index, y_index)*(Vp+Vnn*(((potential_site((x_index+1)%Nx), y_index))+potential_site((x_index-1)%Nx, y_index)+potential_site(x_index, (y_index+1)%Ny)+potential_site(x_index, (y_index-1)%Ny)));
	
	double energy_change=-2.0*potential_site(x_index, y_index)*(Vp+Vnn*((potential_site(index_xf, y_index))+potential_site(index_xb, y_index)+potential_site(x_index, index_yf)+potential_site(x_index, index_yb)));
	//double energy_change=-2*potential_site(x_index, y_index)*(Vp+Vnn*(((potential_site((x_index+1)%Nx), y_index))+potential_site(x_index, (y_index+1)%Ny)));
	return energy_change;
}


//Conditionally accept new microstate
bool accept_microstate(int Nx, int Ny, Eigen::MatrixXd& microstate, const double& chemical_potential, const double& T, const int& x_index, const int& y_index, std::uniform_real_distribution<>& dis, std::mt19937_64& gen)
{
	
	double delta_energy=change_in_energy(Nx, Ny, microstate, x_index, y_index);
	double k=0.08617333262; //meV/K 
	//double test_number=rand()%1;
	double test_number=dis(gen);
	//int N;
	//auto N=get_N(testing_microstate, Nx, Ny);
	//if ((delta_energy-chemical_potential*(N))<0)
	if ((delta_energy-((-1.0)*microstate(x_index, y_index)*chemical_potential))<0)
	{
		return true;
	}
	else if (std::exp((delta_energy-(-1.0)*chemical_potential*microstate(x_index,y_index))*(-1.0/(k*T)))>=test_number)  //may have to multiply by avg comp rather than testmicrostate(x,y)
	{
		return true;
	}
	else{
		return false;
	}
}

//Get average thermoqualities using importance sampling
double A(std::vector<double>& total_energy, int i)
{
	double sum_of_energies;
	for (auto& n : total_energy)
    		sum_of_energies += n;
	auto average=(1.0/i)*sum_of_energies;
	return average;
}

//Convergence Test
bool converged(std::vector<double>& energy_container)
{
//	std::cout<<"not converged";
	int size=energy_container.size();
	if (size>=4)
		{
		if ((energy_container.back()-energy_container.at(size-2))< 0.01)
			return true;
		}
	else if (size<4)
		return false;
	else 
		return false;
}


//get comp average
double comp_average(std::vector<double>& comp_container) 
{
	//auto comp_array=microstate.array();
	//auto frac_comp=(comp_array.sum()/Nx/Ny);
	//return frac_comp;
	double sum=0.0;
	//for (auto& n:comp_container)
	for(double n:comp_container)
	{
		sum+=n;
	}
	double av_sum=sum/comp_container.size();
	return av_sum;
}



double LHS_Energy_mean(std::vector<double>& energy_container)
{
	double av_sum=0.0;
	auto N=energy_container.size();
	double average;
	for (auto& n : energy_container)
	{
		//sum+=n;
		av_sum+=pow(n,2.0)/N;
	}
	return av_sum;
}


double RHS_Energy_mean(std::vector<double>& energy_container)
{
	double sum=0.0;
	auto N=energy_container.size();
	double average;
	for (auto& n: energy_container)
	{
		sum+=n;
	}
	average=(pow(sum/N,2.0))/N;
	return average;
}


double Cp_average(std::vector<double>& energy_container, const double& T)
{	
	double k=0.08617333262; //meV/K 
	auto average=(LHS_Energy_mean(energy_container)-RHS_Energy_mean(energy_container))/(k*pow(T,2.0));
	return average;
}

//run MC and get thermo averages?
std::vector<double> do_MC(const double& chemical_potential, const double& T, Eigen::MatrixXd *microstate)
{
	int max_pass=2000;
	//auto current_microstate=get_starting_microstate();
	Eigen::MatrixXd& current_microstate=*microstate;
	auto Nx=current_microstate.cols();
	auto Ny=current_microstate.rows();
	std::random_device my_rand;
        std::mt19937_64 gen(my_rand());
	std::uniform_real_distribution<> dis(0,1);
	int x_index;
	int y_index;    	
	//auto x_index=int(random_num()*Nx);
	//auto y_index=int(random_num()*Ny);
	auto current_energy=total_energy(Nx, Ny, current_microstate, chemical_potential, T);	
	double comp=get_N(current_microstate, Nx, Ny);
	std::vector<double> comp_container;
	std::vector<double> new_energy_container;
	std::vector<double> average_thermodynamics;
	std::vector<double> Cp_container;
	int i=0;
	//while (!converged(new_energy_container))
	//while(i<max_pass)
	for ( int i=0; i<max_pass; i++)
	{
		for (int loop=0; loop<5000; loop++)
		{
			
			x_index=int(dis(gen)*Nx)%Nx;
			y_index=int(dis(gen)*Ny)%Ny;
			double site=(-1.0)*current_microstate(x_index,y_index);
			bool check=accept_microstate(Nx, Ny, current_microstate, chemical_potential, T, x_index, y_index, dis, gen);
			if (check==true)
			{
				current_energy+=change_in_energy(Nx, Ny, current_microstate, x_index, y_index)-(-1.0*current_microstate(x_index, y_index)*chemical_potential); 
				comp+=site;
				current_microstate(x_index,y_index)=site; 
			//	i++;
			}
		}
		
		if (i>=3)
		{
			new_energy_container.push_back(current_energy);
                        //current_microstate(x_index, y_index)=site;
			comp_container.push_back(comp);
		}
	}
	//		if (i<500)
//		{
//			if (check==true)
//			{	
			//current_energy=total_energy(Nx, Ny, current_microstate, x_index, y_index, chemical_potential, T);
//			current_energy+=change_in_energy(Nx, Ny, current_microstate, x_index, y_index)-(-1.0*current_microstate(x_index, y_index)*chemical_potential); 
//			comp+=site;
//			current_microstate(x_index,y_index)=site;
//			i++;
//			}
//			
//		}
//		if (i>=500)
//		{	
//			if (check==true)
//			{
			//std::cout<<"In if";
			//current_energy=total_energy(Nx, Ny, current_microstate, x_index, y_index, chemical_potential, T);
//			current_energy+=change_in_energy(Nx, Ny, current_microstate, x_index, y_index)-(-1.0*current_microstate(x_index, y_index)*chemical_potential); 
//		        new_energy_container.push_back(current_energy);
//			current_microstate(x_index, y_index)=site;
//			comp+=site;
//			comp_container.push_back(comp);
//	                i++;
//			}
		//}
	//		std::cout<<"\n\n"<<current_microstate;
	
        auto av_Cp=Cp_average(new_energy_container, T);
	auto av_comp=comp_average(comp_container);
	average_thermodynamics.push_back(av_comp);
	average_thermodynamics.push_back(av_Cp);
	
	return average_thermodynamics;
	
//	std::vector<double> whatever;
//	return whatever;
}

std::string to_string(int i){
   std::stringstream s;
   s << i;
   return s.str();
}

int main()
{
	//std::cout<<"check!";

	srand((unsigned) time(0));		

	int T_iterations=30;
	int mu_iterations=30;
	double max_T=3000;
	double min_T=0;
	double max_mu=500;
	double min_mu=-500;	
	std::ofstream outputfile;
	//outputfile.open("muvscomp.txt");
	Eigen::MatrixXd Cp(T_iterations, mu_iterations);
	Eigen::MatrixXd Energy(T_iterations, mu_iterations);
	Eigen::MatrixXd Comp(T_iterations, mu_iterations);
        Eigen::MatrixXd microstate=get_starting_microstate();
//	Eigen::ArrayXd Temporary_T[T_iterations-1];
//	Eigen::ArrayXd Temporary_mu[mu_iterations-1];


//        for (int n=0; n<T_iterations; n++)
//	{
//		Temporary_T[n]=min_T+((max_T-min_T)/T_iterations*n);	
//	}
 //       for (int m=0; m<mu_iterations; m++)
//	{
//		Temporary_mu[m]=min_mu+((max_mu-min_mu)/mu_iterations*m);	
//	}
//	Eigen::MatrixXd mu_comp(T_iterator-1, mu_iterator-1);
//	Eigen::MatrixXd mu(T_iterator-1, mu_iterator-1);
	Eigen::VectorXd my_comp=Comp.row(1);
	Eigen::ArrayXd Temporary_T=Eigen::ArrayXd::LinSpaced(T_iterations, max_T, min_T);
	Eigen::ArrayXd Temporary_mu=Eigen::ArrayXd::LinSpaced(mu_iterations, min_mu, max_mu);
	//Vector2d plot1;
	//Vector2d plot2;
	//Vector2d plot3;
	
	for (int i=0; i<mu_iterations; i++)
//	for (auto& T:Temporary_T)
	{
		double mu=Temporary_mu(i);
		for (int j=0; j<T_iterations; j++)
		{
	//		std::cout<<"In loop";
			double T=Temporary_T(j);

		//std::vector<double> MC_results=do_MC(mu, T, x_index, y_index);
			std::vector<double> MC_results=do_MC(mu, T, &microstate);
		
			Comp(j,i)=MC_results[0];
			Cp(j,i)=MC_results[1];

		}
        my_comp=Comp.col(i).transpose();
	for (int l=0; l<T_iterations; l++)
	{
		my_comp(l)=my_comp(l)/80/80;
		Comp(l,i)=my_comp(l);
	}

	auto comp_forT=Comp.transpose();		
//	Eigen::VectorXd TotalComp=Comp.row(mu_iterations);
//	Eigen::VectorXd TotalT=Temporary_T.row(mu_iterations);
	//for (int m=0; m<mu_iterations; m++)
	//{
//		std::cout<<Cp(m,i)<<" ";
	//	std::cout<<Comp(m,i)<< " ";
	//	std::cout<<Temporary_T(m)<<'\n';
	//}
//	for

	//outputfile.close();
	//std::cout<<my_comp<<'\n';
	//std::cout<<Temporary_mu<<'\n';
	//gather data for constant T, vary mu vs comp (hopefully...)
	//auto mu_comp=Comp.row(i);
	//auto mu=Temporary_mu;
	//Eigen::MatrixXd plot1=Eigen::MatrixXd(i,2);
	//plot1(i, 0)=Comp.row(i);
	//plot1(i, 1)=chem_pots;
	//outputfile<<plot1;
	//outputfile.close();



	}
	std::string front("muvscomp_");
	std::string base(".txt");
	for (int i=0; i<T_iterations; i++)	
	{
		
		outputfile.open(front+to_string(floor(Temporary_T(i)))+base);
		for (int k=0; k<mu_iterations;k++)
	
		{
			outputfile<<Comp(i, k)<<" ";
			outputfile<<Temporary_mu(k)<<'\n';
		}
		outputfile.close();
		//std::cout<<'\n';
	}

}	
