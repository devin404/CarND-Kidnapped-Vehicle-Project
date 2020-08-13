/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  
  std::default_random_engine gen;
  double std_x, std_y, std_theta;  // Standard deviations for x, y, and theta

  // Standard deviations for x, y, and theta
  std_x = std[0];
  std_y = std[1];
  std_theta = std[2];
  
  // Creates a normal (Gaussian) distribution for x, y and theta
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> heading_theta(theta, std_theta);

  for(int i=0; i<(num_particles); i++)
  {
    Particle p;
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = heading_theta(gen);
    p.weight = 1;
    
    particles.push_back(p);
    weights.push_back(1);
  }

  is_initialized = 1;
  
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

   
   double pred_x;
   double pred_y;
   double pred_theta;
   double particle_x;
   double particle_y;
   double particle_theta;
  
    // Adding Nosie
  std::default_random_engine gen;
  double std_x, std_y, std_theta;  // Standard deviations for x, y, and theta

  // Standard deviations for x, y, and theta
  std_x = std_pos[0];
  std_y = std_pos[1];
  std_theta = std_pos[2];
  
  // Creates a normal (Gaussian) distribution for x, y and theta
  normal_distribution<double> dist_x(0, std_x);
  normal_distribution<double> dist_y(0, std_y);
  normal_distribution<double> heading_theta(0, std_theta);  
  std::cout<<"Yaw"<<yaw_rate<<"\n";
  
    for(int i =0; i< num_particles; i++)
    {   
       particle_x = particles[i].x;
       particle_y = particles[i].y;
       particle_theta = particles[i].theta;
         
      //No change in yaw
      if (fabs(yaw_rate) < 0.0001) 
      {  
        pred_x  = particle_x + velocity * delta_t * cos(particle_theta);
        pred_y  = particle_y + velocity * delta_t * sin(particle_theta);      
      }
      else
      { 
        //std::cout<<"\nPredict"<<particles[i].x<<" "<<particles[i].y;
        pred_x = particle_x + (velocity / yaw_rate) * (sin(particle_theta + yaw_rate * delta_t) - sin(particle_theta));
        pred_y = particle_y + (velocity / yaw_rate) * (cos(particle_theta) - cos(particle_theta + yaw_rate * delta_t));
        pred_theta = particle_theta + yaw_rate * delta_t;
        //std::cout<<" "<<particles[i].x<<" "<<particles[i].y;        
      } 
  
  
  particles[i].x = pred_x + dist_x(gen);
  particles[i].y = pred_y + dist_y(gen);
  particles[i].theta = pred_theta + heading_theta(gen);
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
   
   double distance;  
   double nearest_distance; 
   double pred_x, pred_y;
   double obs_x, obs_y;
   int id;
  
   for(unsigned int i = 0; i<observations.size();i++)
   {         
     nearest_distance = std::numeric_limits<double>::max();
     obs_x = observations[i].x;
     obs_y = observations[i].y;
     
     for(unsigned int k = 0; k<predicted.size();k++)
     {
       pred_x = predicted[k].x;
       pred_y = predicted[k].y; 
       distance = dist(obs_x, obs_y, pred_x, pred_y);
       
       if(distance < nearest_distance)
       {    
         id = predicted[k].id;
         nearest_distance = distance;      
       }         
     }
     observations[i].id = id;     
   }    
}
  

double multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent); 
  
  return weight;
}  
  
  
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
   double weight_n; 
   double weight;
   int id;
   double mu_x, mu_y; 
     
  for(int i =0; i< num_particles; i++)
  {
    weight = 1;
    //normalize_weight = 0;
    vector<LandmarkObs> predicted;
    vector<LandmarkObs> transformed_observations;
    
    for(unsigned int k =0; k< map_landmarks.landmark_list.size(); k++)
    {
      if(dist(particles[i].x, particles[i].y,
              map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f)<sensor_range)
      {
        LandmarkObs L;
        L.x = map_landmarks.landmark_list[k].x_f;
        L.y = map_landmarks.landmark_list[k].y_f;
        L.id = map_landmarks.landmark_list[k].id_i;
        
        predicted.push_back(L);
        
      }
      
    }
    
       
    for(unsigned int j =0; j< observations.size(); j++)
    {
      LandmarkObs T;

      // transform to map x coordinate
      T.x = particles[i].x + (cos(particles[i].theta) * observations[j].x) 
              - (sin(particles[i].theta) * observations[j].y);

      // transform to map y coordinate
      T.y = particles[i].y + (sin(particles[i].theta) * observations[j].x)
              + (cos(particles[i].theta) * observations[j].y);
      
      T.id = observations[j].id;
      
      //std::cout<<"\nT "<<T.x<<" "<<T.y;
      transformed_observations.push_back(T);
      
     }
    
     //Landmark Association 
     dataAssociation(predicted, transformed_observations);
   

    
     for(unsigned int n =0; n< transformed_observations.size(); n++)
     {
       id = transformed_observations[n].id;

       for(unsigned int k =0; k< predicted.size(); k++)
       {
         if(predicted[k].id == id)
         {
           mu_x = predicted[k].x;
           mu_y = predicted[k].y; 
           break;
         }

       }
       
       weight_n = multiv_prob(std_landmark[0], std_landmark[1], transformed_observations[n].x, transformed_observations[n].y,
                                      mu_x, mu_y);
       weight *=weight_n; 
      }
    particles[i].weight = weight;
    weights[i] = particles[i].weight; 
  }
    
}


void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  std::default_random_engine gen;
  std::discrete_distribution<int> d_weights(weights.begin(), weights.end());
  vector<Particle> particle_temp;
  
  for (int i = 0; i < num_particles; i++) 
  {
    particle_temp.push_back(particles[d_weights(gen)]);
  }
	particles = particle_temp;
    
}



void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}