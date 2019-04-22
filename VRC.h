#pragma once
#include<cmath>
#include<vector>
#include<ctime>
#include<random>
#include<iostream>

using namespace std;

void print(vector<double>pos)
{
	for (int i = 0; i < pos.size(); i++)
		cout << pos[i] << " ";
	cout << endl;
}


bool better(double a, double b)
{
	if (a < b)
		return true;
	else
		return false;
}


double randDouble(double min, double max)
{
	static default_random_engine engine(time(nullptr));
	uniform_real_distribution<double> dis(min, max);
	return dis(engine);
}





class Peakfitness
{
public:

	Peakfitness()
	{
		peaks = 10;
		h_max = 70;
		h_min = 30;
		w_max = 12;
		w_min = 1;
		h_s = 7;
		w_s = 1;

		fi = 0;

		for (int i = 0; i < peaks; i++)
		{
			height.push_back(randDouble(h_min, h_max));
			width.push_back(randDouble(w_min, w_max));
		}
		print(height);
	}

	double fiX(int g,double s)
	{
		return g * s;
	}


	void update()
	{
		static default_random_engine engine(time(nullptr));
		static normal_distribution<double>gaussin(0, 1);
		

		for (int i = 0; i < peaks; i++)
		{
			height[i] = height[i] + h_s * gaussin(engine);
			if (height[i] > h_max)
				height[i] = h_max;
			else if (height[i] < h_min)
				height[i] = h_min;
			//cout << "height: " << height[i] << endl;

			width[i] = width[i] + w_s * gaussin(engine);
			if (width[i] > w_max)
				width[i] = w_max;
			else if (width[i] < w_min)
				width[i] = w_min;
		}
		
		
	}

	double fitnessFunction_used(vector<double>pos,int g,double s)
	{
		double result = 0.0;
		vector<double>compare;
		vector<vector<double>>parameter;

		for (int i = 0; i < pos.size(); i++)
		{
			double temp = 0.0;
			temp = (pos[i] + fiX(g,s));
			result += pow(temp, 2);
		}

		for (int i = 0; i < peaks; i++)
		{
			compare.push_back(-(height[i]) / (1 + width[i] * result));
		}
		
		int index = 0;
		for (int i = 1; i < compare.size(); i++)
		{
			if(better(compare[i],compare[index]))
				index = i;
		}
		

		return compare[index];
	}

	double fi;
	double h_max;
	double h_min;
	double w_max;
	double w_min;
	double h_s;
	double w_s;
	int peaks;
	double s;

	vector<double>height;
	vector<double>width;
};





//VRC functions

class Particle
{
public:
	vector<double> position;
	vector<double>velocity;
	vector<double>pBest;

	Particle() {}

	Particle(vector<double> position, vector<double>velocity, vector<double>best_position, double best_fitness)
	{
		this->position = position;
		this->velocity = velocity;
		this->pBest = best_position;

	}
};

class VRC
{
public:
	VRC(vector<double>fitness, vector<Particle>particles,vector<vector<double>>mean_Xid_g, 
		vector<double>mean_fi_g, int g,double s,Peakfitness &fit,double dmax,double dmin)
	{
		this->dmax = dmax;
		this->dmin = dmin;
		this->dmin;
		this->g = g;
		cout << "g:" << g << endl;
		this->s = s;
		this->fitness = fitness;
		this->particles = particles;
		this->mean_Xid_g = mean_Xid_g;
		this->mean_fi_g = mean_fi_g;
		

		relo = 1e-62;
		this->fit = fit;
		Avg_Pos_Progress();
		
		Avg_sensitivity_all();
		Avg_sensitivity();
		Calc_fitness_diff();
		
		find_fit_best();
		Relocation_radius();

		
		Relocation_radius_id();
		
		relocation();
	}

	void Avg_Pos_Progress()
	{
		vector<double>result;
		for (int i = 0; i < mean_Xid_g.size(); i++)
		{
			result.push_back(0);
		}

		for (int i = 0; i < mean_Xid_g.size(); i++)
		{ 
			for (int j = 0; j < mean_Xid_g[0].size(); j++)
			{
				result[i] += pow(mean_Xid_g[i][j], 2);
			}
			result[i] = sqrt(result[i]);
		}
		
		Xi = result;
	}

	 void Avg_sensitivity_all()
	{
		vector<double>result = Xi;
		for (int i = 0; i < Xi.size(); i++)
		{
			result[i] = mean_fi_g[i] / (Xi[i] + relo);
		}
		Si = result;
	}

	void Avg_sensitivity()
	{
		vector<vector<double>>result = mean_Xid_g;
		for (int i = 0; i < result.size(); i++)
		{
			for (int j = 0; j < result[0].size(); j++)
			{
				result[i][j] = mean_Xid_g[i][j] * Si[i] / (Xi[i] + relo);
			}
		}

		Sid = result;
	}

	void Calc_fitness_diff()
	{
		vector<double>fit_current;
		for (int i = 0; i < particles.size(); i++)
			fit_current.push_back(fit.fitnessFunction_used(particles[i].position,g,s));
		fi = fit_current;
		for (int i = 0; i < particles.size(); i++)
			fi[i] = fit_current[i] - fitness[i];
	}

	double find_fit_best()
	{
		int index = 0;
		for (int i = 0; i < particles.size(); i++)
		{
			if (better(fit.fitnessFunction_used(particles[i].position,g,s), fit.fitnessFunction_used(particles[index].position,g,s)))
				index = i;
		}
		return fit.fitnessFunction_used(particles[index].position,g,s);
	}

	void Relocation_radius()
	{
		Ri = fi;
		for (int i = 0; i < Ri.size(); i++)
		{
			if (fi[i] <= 0)
			{
				Ri[i] = -fi[i] / (Si[i] + relo);
			}
			else
			{
				Ri[i] = min((find_fit_best() - fitness[i])/(Si[i] + relo),fi[i]/(Si[i] + relo));
			}

		}
	}

	void Relocation_radius_id()
	{
		Rid = Sid;
		for (int i = 0; i < Rid.size(); i++)
		{
			for (int j = 0; j < Rid[0].size(); j++)
			{
				Rid[i][j] = Ri[i] * mean_Xid_g[i][j] / (Xi[i] + relo);
				if (Rid[i][j] > (dmax - dmin))
					Rid[i][j] = dmax - dmin;
				else if (Rid[i][j] < -(dmax - dmin))
					Rid[i][j] = -dmax + dmin;
				if (Rid[i][j] < relo && Rid[i][j] > -relo)
					Rid[i][j] = relo;
			}
		}
	}

	void relocation()
	{
		for (int i = 0; i < particles.size(); i++)
		{
			for (int j = 0; j < Rid[0].size(); j++)
			{
				double lastpos = particles[i].position[j];
				particles[i].position[j] = particles[i].position[j] + randDouble(0, 1)*Rid[i][j];
				if (particles[i].position[j] > dmax)
				{
					particles[i].position[j] = dmax + lastpos - particles[i].position[j];
				}
				else if (particles[i].position[j] < dmin)
				{
					particles[i].position[j] = dmin + lastpos - particles[i].position[j];
				}
			}
		}
	}


	vector<Particle> particles;

private:

	vector<double>fitness;
	int g;
	double s;
	double relo;
	double dmax;
	double dmin;
	
	Peakfitness fit;
	vector<vector<double>>mean_Xid_g;
	vector<double>mean_fi_g;

	vector<double>Xi;
	vector<double>Si;
	vector<vector<double>>Sid;
	vector<double>fi;
	vector<double>Ri;
	vector<vector<double>>Rid;
};
