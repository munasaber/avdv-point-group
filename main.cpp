#include "eigen-git-mirror/Eigen/Dense"
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <ostream>
#include <ctime>
#include <cstdlib> 
#include <time.h> 
//#include "../matplotlib-cpp/matplotlibcpp.h"
#define PREC 1E-7;
//namespace plt =matplotlibcpp;


float Vnn=40.0; //meV
float Vp=0.0;//meV
float V0=2*Vnn;//meV
//Choose Starting microstate
Eigen::MatrixXf get_starting_microstate()  //In reality, this would need to be Xf
{
	//Initialize Matrix with all ones
	Eigen::MatrixXf starting_microstate = Eigen::MatrixXf::Ones(50,50); 
	return starting_microstate;
}


//member functions for site class
//get random x_index based off of Nx
int x_Index(int Nx) {
	float random_number=rand()%Nx;
	return floor(random_number);
}

//get random y_index based off of Ny
int y_Index(int Ny) {
	float random_number=rand()%Ny;
	return floor(random_number);
}


//Flip site on previous microstate
Eigen::MatrixXf flip_site(int Nx, int Ny, Eigen::MatrixXf previous_microstate, int x_index, int y_index) //change to MatrixXf
{
	previous_microstate(x_index, y_index)=-previous_microstate(x_index, y_index);
//	previous_microstate.coeffRef(5, 10)=-previous_microstate.coeffRef(5,10);
	//	for (int i=0; i<previous_microstate.rows(); i++)
//	{
//		for (j=0; j<previous_microstate.cols(); j++)
//		{
//			if (i==x_index) //might be remaking index with every 
//			{
//				if (j==y_index)
//				{
//					previous_microstate[i][j]=-previous_microstate[i][j]
	return previous_microstate;
}


//get N
float get_N(Eigen::MatrixXf microstate, int Nx, int Ny)
{	
	auto config=microstate;
	auto act_config=config.array()+1.000;
	return act_config.sum()/2;
}
//Calculate total grand canonical energy
float total_energy(int Nx, int Ny, Eigen::MatrixXf microstate, int x_index, int y_index, const float& mu, const float& T)
{
	float k=0.08617333262; //meV/K 
	float chem_pot;
	float energy;
	for (int i=0; i<Nx; i++)
	{
		for (int j=0; j<Ny; j++)
		{
		//	float temp_eng=Vp*microstate(x_index, y_index)+Vnn*(microstate(x_index,y_index))*(((microstate((x_index+1)%Nx, y_index)+microstate(x_index, (y_index+1)%Ny))));
			float temp_eng=Vp*microstate(i, j)+Vnn*(microstate(i,j))*(((microstate((i+1)%Nx, j)+microstate(i, (j+1)%Ny))));
			energy+=temp_eng;
		}
			//float energy= Vp*microstate(x_index,y_index)+Vnn*(microstate(x_index, y_index))*(((microstate((x_index+1)%Nx), y_index))+microstate(x_index, (y_index+1)%Ny)+microstate((x_index-1)%Nx, y_index)+microstate(x_index, (y_index-1)%Ny));
	}
	chem_pot=get_N(microstate, Nx, Ny)*mu;
	float tot_energy=(energy-chem_pot)/k/T;
	return tot_energy;
}

//Calculate change in grand canonical energy
float change_in_energy(int Nx, int Ny, Eigen::MatrixXf current_microstate, int x_index, int y_index)
{
 	Eigen::MatrixXf potential_site=current_microstate;
	float energy_change=-2*potential_site(x_index, y_index)*(Vp+Vnn*(potential_site(x_index, y_index))*(((potential_site((x_index+1)%Nx), y_index))+potential_site(x_index, (y_index+1)%Ny)));
        return energy_change;
}


//Conditionally accept new microstate
bool accept_microstate(int Nx, int Ny, Eigen::MatrixXf previous_microstate, Eigen::MatrixXf testing_microstate, const float& chemical_potential, const float& T, int x_index, int y_index)
{
//	std::cout<<"In accept microstate";
	float delta_energy=change_in_energy(Nx, Ny, testing_microstate, x_index, y_index);
	float k=0.08617333262; //meV/K 
	Eigen::MatrixXf ending_microstate;
	float test_number=rand()%1;
	//int N;
	auto N=get_N(testing_microstate, Nx, Ny);
	if ((delta_energy-chemical_potential*(N))<0)
	{
		return true;
	}
	else if (exp(-(1/(k*T))*(delta_energy-chemical_potential*(N)))>=test_number)
	{
		return true;
	}
	else{
		return false;
	}
}

//Get average thermoqualities using importance sampling
float A(std::vector<float>& total_energy, int i)
{
	float sum_of_energies;
	for (auto& n : total_energy)
    		sum_of_energies += n;
	auto average=(1/i)*sum_of_energies;
	return average;
}

//Convergence Test
bool converged(std::vector<float>& energy_container)
{
	int size=energy_container.size();
	if (size>2)
		{
		if ((energy_container.back()-energy_container.at(size-2))< 0.0001)
			return true;
		}
	else if (size<3)
		return false;
}


//get comp average
float comp_average(std::vector<float>& comp_container) 
{
	//auto comp_array=microstate.array();
	//auto frac_comp=(comp_array.sum()/Nx/Ny);
	//return frac_comp;
	float sum=0.0;
	for (auto& n:comp_container)
	{
		sum+=n;
	}
	float av_sum=sum/comp_container.size();
	return av_sum;
}



float LHS_Energy_mean(std::vector<float>& energy_container)
{
	float av_sum=0.0;
	auto N=energy_container.size();
	float average;
	for (auto& n : energy_container)
	{
		//sum+=n;
		av_sum+=pow(n,2.0)/N;
	}
	return av_sum;
}


float RHS_Energy_mean(std::vector<float>& energy_container)
{
	float sum=0.0;
	auto N=energy_container.size();
	float average;
	for (auto& n: energy_container)
	{
		sum+=n;
	}
	average=(pow(sum/N,2.0))/N;
	return average;
}


float Cp_average(std::vector<float>& energy_container, const float& T)
{	
	float k=0.08617333262; //meV/K 
	auto average=(LHS_Energy_mean(energy_container)-RHS_Energy_mean(energy_container))/(k*pow(T,2.0));
	return average;
}

//run MC and get thermo averages?
std::vector<float> do_MC(const float& chemical_potential, const float& T)
{
	//std::cout<<"do_MC";
	//put these in the arguments of the function when finished debugging
	//float chemical_potential=50.0;
	//float T=50.0;
	//Site site;	
	auto current_microstate=get_starting_microstate();
	auto Nx=current_microstate.cols();
	auto Ny=current_microstate.rows();
	auto x_index=x_Index(50);
	auto y_index=y_Index(50);
	auto starting_energy=total_energy(Nx, Ny, current_microstate, x_index, y_index, chemical_potential, T);
	
	//auto starting_comp=comp(Nx,Ny, current_microstate);
	std::vector<float> comp_container;
	std::vector<float> new_energy_container;
	new_energy_container.push_back(starting_energy);
	//comp_container.push_back(starting_comp);
	std::vector<float> average_thermodynamics;
	std::vector<float> Cp_container;
	while (!converged(new_energy_container))
	{
		//std::cout<<"In converged!";
   	
		auto x_index=x_Index(50);
		auto y_index=y_Index(50);
		float current_energy=total_energy(Nx,Ny,current_microstate, x_index, y_index, chemical_potential, T);
		auto prospective_microstate=flip_site(Nx,Ny, current_microstate, x_index, y_index); 
		//float energy_change=change_in_energy(Nx, Ny, current_microstate, x_index, y_index);
		bool check=accept_microstate(Nx, Ny, current_microstate, prospective_microstate, chemical_potential, T, x_index, y_index);
		if (check==true)
		{
			//std::cout<<"In if";
			current_microstate=prospective_microstate;
			float current_energy=total_energy(Nx, Ny, current_microstate, x_index, y_index, chemical_potential, T);
		        new_energy_container.push_back(current_energy);
			//auto N=get_N(current_microstate, Nx, Ny);
			//comp_container.push_back(N);

  		       	//auto flipped=flip_site(50, 50, current_microstate, x_index, y_index);
		}

			//std::cout<<"\n\n"<<current_microstate;
			auto N=get_N(current_microstate, Nx, Ny);
			comp_container.push_back(N);
	}
        auto av_Cp=Cp_average(new_energy_container, T);
	auto av_comp=comp_average(comp_container);
	average_thermodynamics.push_back(av_Cp);
	average_thermodynamics.push_back(av_comp);
	return average_thermodynamics;
	
//	std::vector<float> whatever;
//	return whatever;
}

int main()
{
	//std::cout<<"check!";

	srand((unsigned) time(0));		

	int T_iterations=20;
	int mu_iterations=20;
	float max_T=1000;
	float min_T=100;
	float max_mu=0.5;
	float min_mu=-0.5;	
	std::ofstream outputfile;
	outputfile.open("muvscomp.txt");
	Eigen::MatrixXf Cp(T_iterations, mu_iterations);
	Eigen::MatrixXf Energy(T_iterations, mu_iterations);
	Eigen::MatrixXf Comp(T_iterations, mu_iterations);

//	Eigen::ArrayXf Temporary_T[T_iterations-1];
//	Eigen::ArrayXf Temporary_mu[mu_iterations-1];


//        for (int n=0; n<T_iterations; n++)
//	{
//		Temporary_T[n]=min_T+((max_T-min_T)/T_iterations*n);	
//	}
 //       for (int m=0; m<mu_iterations; m++)
//	{
//		Temporary_mu[m]=min_mu+((max_mu-min_mu)/mu_iterations*m);	
//	}
//	Eigen::MatrixXf mu_comp(T_iterator-1, mu_iterator-1);
//	Eigen::MatrixXf mu(T_iterator-1, mu_iterator-1);
	Eigen::VectorXf my_comp=Comp.row(1);
	Eigen::ArrayXf Temporary_T=Eigen::ArrayXf::LinSpaced(T_iterations, min_T, max_T);
	Eigen::ArrayXf Temporary_mu=Eigen::ArrayXf::LinSpaced(mu_iterations, min_mu, max_mu);
	//Vector2f plot1;
	//Vector2f plot2;
	//Vector2f plot3;
	
	for (int i=0; i<T_iterations; i++)
//	for (auto& T:Temporary_T)
	{
	float T=Temporary_T(i);
		for (int j=0; j<mu_iterations; j++)
		{
	//		std::cout<<"In loop";
		float mu=Temporary_mu(j);

		//std::vector<float> MC_results=do_MC(mu, T, x_index, y_index);
		std::vector<float> MC_results=do_MC(mu, T);
		
		Comp(i,j)=MC_results[0];
		Cp(i,j)=MC_results[1];
		}
        my_comp=Comp.row(i);
	for (int l=0; l<mu_iterations; l++)
	{
		my_comp(l)=my_comp(l)/50/50;
		Comp(i,l)=my_comp(l);
	}
	auto comp_forT=Comp.transpose();		
//	Eigen::VectorXf TotalComp=Comp.row(mu_iterations);
//	Eigen::VectorXf TotalT=Temporary_T.row(mu_iterations);
	for (int m=0; m<mu_iterations; m++)
	{
//		std::cout<<Cp(m,i)<<" ";
		//TotalComp(m)<<Comp(m,i);
		//Total_T(m)<<Temporary_T(m);
		std::cout<<Comp(m,i)<< " ";
		std::cout<<Temporary_T(m)<<'\n';
	}
//	for

	//outputfile.close();
	//std::cout<<my_comp<<'\n';
	//std::cout<<Temporary_mu<<'\n';
	//gather data for constant T, vary mu vs comp (hopefully...)
	//auto mu_comp=Comp.row(i);
	//auto mu=Temporary_mu;
	//Eigen::MatrixXf plot1=Eigen::MatrixXf(i,2);
	//plot1(i, 0)=Comp.row(i);
	//plot1(i, 1)=chem_pots;
	//outputfile<<plot1;
	//outputfile.close();



	}
}	
