#ifndef TRAJOPT_H
#define TRAJOPT_H

#include <OpenSim/OpenSim.h>
// #include "OpenSim/Common/STOFileAdapter.h"
#include <opencv2/opencv.hpp>

class TrajOpt {

private:
	SimTK::Vector _controls;
	std::string _model_file;	
	int _actuated_dof;
	int _numParams;
	bool isLowDim;
	int _nbr_pcs;
	cv::Mat _eigvectors;
	cv::Mat _mean;
	bool _randomInitial;
	bool _randomTarget;	
	std::vector<double> _rand_coords;
	std::vector<double> _rand_end_eff = std::vector<double>(3);
	std::string _app_filename;
	int _maxIter;
	int _lambda;
	double _bestSoFar;
	double _lower_bound;
	double _upper_bound;
	int _STEP;
	int _start_opt;
	int _curr_step;
	double _NOISE;
	bool _isSpasticity;
	double _G;
	double _velocity_threshold;
	std::vector<SimTK::Vector> _control_list;
	std::vector<SimTK::State> _state_list;
	int _dim_reduction;
	double _init_stepsize;
	double _est_opt_time;

public:	
	TrajOpt() : isLowDim(false), _randomInitial(true), _randomTarget(true), _app_filename("ue"), _maxIter(20), _lambda(20), _bestSoFar(SimTK::Infinity), _lower_bound(-1), _upper_bound(1), _STEP(1), _NOISE(0.), _isSpasticity(false), _G(0), _velocity_threshold(0), _est_opt_time(0.) 
	{
		_dim_reduction = 0;
		_start_opt = 1;
		_curr_step = 1;
		_init_stepsize = 0.3;
	};	

	TrajOpt(std::string model_file) : isLowDim(false), _randomInitial(true), _randomTarget(true), _app_filename("ue"), _maxIter(20), _lambda(20), _bestSoFar(SimTK::Infinity), _lower_bound(-1), _upper_bound(1), _STEP(1), _NOISE(0.), _isSpasticity(false), _G(0), _velocity_threshold(0), _est_opt_time(0.)
	{
		_model_file = model_file;
		_dim_reduction = 0;
		_start_opt = 1;
		_curr_step = 1;
		_init_stepsize = 0.3;
	};

	void set_init_stepsize(double init_stepsize){
		_init_stepsize = init_stepsize;
	}

	void set_dim_reduction(int dim_reduction){
		_dim_reduction = dim_reduction;
	}
	
	void set_NOISE(double NOISE)
	{
		_NOISE = NOISE;
	}

	void set_STEP(int STEP)
	{
		_STEP = STEP;
	}

	void set_start_opt(int start_opt)
	{
		_start_opt = start_opt;
	}	

	void set_lower_upper_bound(double lower_bound, double upper_bound)
	{
		_lower_bound = lower_bound;
		_upper_bound = upper_bound;
	}

	double get_bestSoFar()
	{
		return _bestSoFar;
	}

	void set_opt_time(double opt_time)
	{
		_est_opt_time = opt_time;
	};

	void set_maxIter(int maxIter)
	{
		_maxIter = maxIter;
	}
	
	void set_lambda(int lambda)
	{
		_lambda = lambda;
	}

	void set_app_filename(std::string app_filename)
	{
		_app_filename = app_filename;
	};

	void set_low_dim_params(int nbr_pcs, cv::Mat eigvectors, cv::Mat mean)
	{
		isLowDim = true;
		_nbr_pcs = nbr_pcs;
		_eigvectors = eigvectors.clone();
		_mean = mean.clone();			
	};

	void set_random_initial_state_end_eff(std::vector<double> rand_coords, std::vector<double> rand_end_eff)
	{
		_randomInitial = false;
		_randomTarget = false;
		_rand_coords = rand_coords;
		_rand_end_eff = rand_end_eff;
	};	

	void set_initial_state(std::vector<double> rand_coords)
	{
		_randomInitial = false;
		_rand_coords = rand_coords;
	};
	
	void set_end_eff_target(std::vector<double> rand_end_eff)
	{
		_randomTarget = false;
		_rand_end_eff = rand_end_eff;
	};

	void set_spasticity(bool spasticity, double G, double velocity_threshold)
	{
		_isSpasticity = spasticity;
		_G = G;
		_velocity_threshold = velocity_threshold;
	};
	
	std::vector<double> get_initial_coords()
	{
		return _rand_coords;
	} 	
	std::vector<double> get_end_eff()
	{
		return _rand_end_eff;
	} 	
	
	//SimTK::State traj_opt(const SimTK::State& si, OpenSim::Model& osimModel, double initialTime, std::vector<double> rand_end_eff, std::string filename, std::string filename_states);
	SimTK::State traj_opt(const SimTK::State& si, OpenSim::Model& osimModel, double initialTime, std::string filename, std::string filename_states);
	
	std::vector<double> compute_white_noise(int dimension, double standard_deviation);

	double estimate_optimization_time(const SimTK::State& si, OpenSim::Model& osimModel);

	std::vector<double> set_random_state(const OpenSim::CoordinateSet& coords, SimTK::State& si);

	void set_state(const OpenSim::CoordinateSet& coords, SimTK::State& si, std::vector<double> random_coords);

	int actuated_degrees_of_freedom(const OpenSim::CoordinateSet& coords);

	std::vector<double> set_random_end_eff_target();

	void traj_opt_setup();
	
	int get_numControls()
	{
		OpenSim::Model osimModel(_model_file.c_str());
		osimModel.initSystem();
		return osimModel.getNumControls();
	}
	
	int get_numParams()
	{
		return _numParams;
	}

	int get_actuated_dofs()
	{
		return _actuated_dof;
	}

	SimTK::Vector get_controls()
	{
		return _controls;
	}
	
	std::vector<SimTK::Vector> get_control_list()
	{
		return _control_list;	
	};
	
	std::vector<SimTK::State> get_state_list()
	{
		return _state_list;	
	};
};

#endif
