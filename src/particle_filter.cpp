/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

#define SmallValue 0.0001

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	default_random_engine gen;	
	double std_x, std_y, std_theta; 
	// Standard deviations for x, y, and theta	
	std_x = std[0];	 
	std_y = std[1];	 
	std_theta = std[2];	 	
	// This line creates a normal (Gaussian) distribution for x	
	normal_distribution<double> dist_x(x, std_x);		
	
	normal_distribution<double> dist_y(y, std_y);	
	
	normal_distribution<double> dist_theta(theta, std_theta);		
	for (int i = 0; i < num_particles; ++i) 
		{		
		Particle Part;	
		Part.id = i;
		Part.x = dist_x(gen);		 
		Part.y = dist_y(gen);		 
		Part.theta = dist_theta(gen);	
		Part.weight = 1;
		particles.push_back(Part);

		}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;	
	double std_x, std_y, std_theta; 
	// Standard deviations for x, y, and theta	
	std_x = std_pos[0];	 
	std_y = std_pos[1];	 
	std_theta = std_pos[2];	
	for(auto it = particles.begin();it != particles.end(); ++it)
	{
		double x0 = it->x;
		double y0 = it->y;
		double theta0 = it->theta;



		double xf,yf,thetaf;

		if(fabs(yaw_rate) > SmallValue)
		{
			double temp = velocity/yaw_rate;
			double delta_theta = yaw_rate*delta_t;
			xf = x0 +  temp*(sin(theta0 + delta_theta) - sin(theta0));
			yf = y0 + temp*(cos(theta0) - cos(theta0 + delta_theta));
			thetaf = theta0 + delta_theta;
			
		}
		else
		{
			xf = x0 + velocity*delta_t*cos(theta0);
			yf = y0 + velocity*delta_t*sin(theta0);
			thetaf = theta0;
		}

		normal_distribution<double> dist_x(xf, std_x);    
		normal_distribution<double> dist_y(yf, std_y);    
		normal_distribution<double> dist_theta(thetaf, std_theta);

		it->x = dist_x(gen);
		it->y = dist_y(gen);
		it->theta = dist_theta(gen);

	}
	
	

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for(auto it = observations.begin(); it != observations.end(); ++it)
	{
		auto min_d = dist(it->x,it->y,predicted.begin()->x,predicted.begin()->y);
		for(auto l = predicted.begin();l != predicted.end(); ++l)
		{
			auto distance =dist(it->x,it->y,l->x,l->y);
			if(min_d < distance)
			{
				min_d = distance;
				it->id = l->id;
			}
			
		}
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	for(auto p = particles.begin();p != particles.end(); ++p)
	{
		vector<LandmarkObs> p_obs;
		//map observations for particle from particle coordinate to map coordinate
		for(auto o = observations.begin();o != observations.end(); ++o)
		{
			
			p_obs.x = p->x + cos(p->theta)*o->x-sin(p->theta)*o->y;
			p_obs.y = p->y + sin(p->theta)*o->x + cos(p->theta)*o->y;
		}

		//do data association based on the landmarks with sensor range
		vector<LandmarkObs> landmarks_p_obss;
		for(auto l = map_landmarks.landmark_list.begin();l != map_landmarks.landmark_list.end(); ++l)
		{
			auto distance = dist(p->x,p->y,l->x,l->y);
			if(distance <= sensor_range)
			{
				LandmarkObs landmarks_p_obs(l->id, l->x, l->y);
				landmarks_p_obss.push_back(landmarks_p_obs);
			}
		}

 		dataAssociation(landmarks_p_obss, p_obs);

		p->weight = 1;

		
		double sig_x = std_landmark[0];
		double sig_y = std_landmark[1];

		// calculate normalization term
		double gauss_norm = (1/(2 * M_PI * sig_x * sig_y));

		for(auto it = p_obs.begin();it != p_obs.end(); ++it)
		{
			double x_obs= it->x;
			double y_obs= it->y;
			for(auto l = map_landmarks.begin();l != map_landmarks.end(); ++l)
			{
				if(it.id == l.id)
				{
					double mu_x= l.x;
					double mu_y= l.y;
								// calculate exponent
					double exponent= ((x_obs - mu_x)**2)/(2 * sig_x**2) + ((y_obs - mu_y)**2)/(2 * sig_y**2)
					// calculate weight using normalization terms and exponent
					p->weight *= gauss_norm * exp(-exponent);

				}

			}



		}

		weights.push_back(p->weight);
		
	}

	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	  std::default_random_engine gen;  
	std::discrete_distribution<> d(weights.begin(), weights.end());  
	std::vector<Particle> new_particles;  
	for (size_t i = 0; i < particles.size(); ++i) 
	{    
		const Particle &src = particles[d(gen)];    
		new_particles.push_back(src);  
	}  
	particles.clear();  
	particles.insert(particles.end(), new_particles.begin(), new_particles.end());
	

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}