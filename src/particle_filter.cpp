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

	cout << "init" << is_initialized << endl;
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
	for(int i=0; i<num_particles; ++i)
	{
		double x0 = particles[i].x;
		double y0 = particles[i].y;
		double theta0 = particles[i].theta;



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

		cout << "predict" << xf << " " << yf << " " << thetaf << endl;

		normal_distribution<double> dist_x(xf, std_x);    
		normal_distribution<double> dist_y(yf, std_y);    
		normal_distribution<double> dist_theta(thetaf, std_theta);

		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);

	}
	
	cout <<"std_x" << std_x << endl;

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for(int i = 0; i<observations.size();++i )
	{
		auto min_d = dist(observations[i].x,observations[i].y,predicted[0].x,predicted[0].y);
		for(int j=0; j<predicted.size();++j)
		{
			auto distance =dist(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y);
			if(min_d < distance)
			{
				min_d = distance;
				observations[i].id = predicted[j].id;
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

	for(int i=0; i<num_particles; ++i)
	{
		vector<LandmarkObs> observationsPart;
		LandmarkObs p_obs;
		//map observations for particle from particle coordinate to map coordinate
		for(int j=0; j<observations.size();++j)
		{
			
			p_obs.x = particles[i].x + cos(particles[i].theta)*observations[j].x-sin(particles[i].theta)*observations[j].y;
			p_obs.y = particles[i].y + sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y;
			observationsPart.push_back(p_obs);
		}

		//do data association based on the landmarks with sensor range
		vector<LandmarkObs> landmarks_p_obss;
		for(int k=0; k < map_landmarks.landmark_list.size();++k)
		{
			auto distance = dist(particles[i].x,particles[i].y,map_landmarks.landmark_list[k].x_f,map_landmarks.landmark_list[k].y_f);
			if(distance <= sensor_range)
			{
				LandmarkObs landmarks_p_obs;
				landmarks_p_obs.id = map_landmarks.landmark_list[k].id_i;
				landmarks_p_obs.x = map_landmarks.landmark_list[k].x_f;
				landmarks_p_obs.y = map_landmarks.landmark_list[k].y_f;
				landmarks_p_obss.push_back(landmarks_p_obs);
			}
		}

 		dataAssociation(landmarks_p_obss, observationsPart);

		particles[i].weight = 1;

		
		double sig_x = std_landmark[0];
		double sig_y = std_landmark[1];

		// calculate normalization term
		double gauss_norm = (1/(2 * M_PI * sig_x * sig_y));

		cout << "2" << endl;

		for(int j=0; j<observationsPart.size();++j)
		{
			double x_obs= observationsPart[j].x;
			double y_obs= observationsPart[j].y;
			for(int k=0; k<landmarks_p_obss.size();++k)
			{
				if(observationsPart[j].id == landmarks_p_obss[k].id)
				{
					double mu_x= landmarks_p_obss[k].x;
					double mu_y= landmarks_p_obss[k].y;
								// calculate exponent
					double exponent= ((x_obs - mu_x)*(x_obs - mu_x))/(2 * sig_x*sig_x) + ((y_obs - mu_y)*(y_obs - mu_y))/(2 * sig_y*sig_y);
					// calculate weight using normalization terms and exponent
					particles[i].weight *= gauss_norm * exp(-exponent);

					cout << "weight " << particles[i].weight << endl;

				}

			}



		}

		weights.push_back(particles[i].weight);
		
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

	cout << "3" << endl;
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
