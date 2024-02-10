#ifndef OPTSYSTEM_H
#define OPTSYSTEM_H

#include <OpenSim/OpenSim.h>
#include <string>
#include <opencv2/opencv.hpp>

class OptSystem : public SimTK::OptimizerSystem { 

private:
	int numControls;
	const SimTK::State& si;
	mutable SimTK::State _sf;
	OpenSim::Model& osimModel;
	OpenSim::Manager& _manager;
	SimTK::TimeStepper& _ts;
	double _initialTime;
	double _finalTime;
	double _window_interval;
	mutable double _bestSoFar;
	std::string _filename;
	std::string _filename_states;
	std::string _filename_excitations;
	std::string _filename_states_learning;
	std::string _filename_hand_kinematics;	
	std::string _model_file;	
	std::vector<double> _rand_end_eff;
	mutable std::vector<SimTK::Vector> _control_list;
	mutable std::vector<SimTK::State> _state_list;	
	bool _isFinal;
	int _KNOTS;
	int _numParams;
	bool _isSpasticity;
	double _G;
	double _velocity_threshold;
	cv::Mat _eigvectors;
	cv::Mat _mean;
	bool _isLowDim;
	int _nbr_pcs;
public:	
	OptSystem(int numParameters,const SimTK::State& s, OpenSim::Model& aModel, OpenSim::Manager& manager, SimTK::TimeStepper& ts, std::string filename, std::string filename_states, double finalTime, double window_interval, int KNOTS) : OptimizerSystem(numParameters), osimModel(aModel), _manager(manager), si(s), _ts(ts), _filename(filename), _filename_states(filename_states), _finalTime(finalTime), _window_interval(window_interval), _isSpasticity(false), _G(0.), _velocity_threshold(0.)  
	{	
		numControls = osimModel.getNumControls();
		_bestSoFar = std::numeric_limits<double>::infinity();
		_initialTime = 0;
		//si = s;
		 _sf = si;
		_isFinal = false;
		_numParams = numParameters;
		_KNOTS = KNOTS;
		_isLowDim = false;
		_filename_excitations = _filename + "_excitations.txt";
		_filename_states_learning = _filename + "_states_learning.txt";

		std::ofstream fo;
		_filename_hand_kinematics = _filename + "_hand_kinematics.txt";	
		fo.open(_filename_hand_kinematics.c_str(), std::ios_base::trunc);
	};

	OptSystem(const OptSystem& rhs) : osimModel(rhs.osimModel), _manager(rhs._manager), _ts(rhs._ts), si(rhs.si), _sf(rhs._sf)
	{
		numControls = rhs.numControls;
		_initialTime = rhs._initialTime;
		_finalTime =	rhs._finalTime;
		_window_interval = rhs._window_interval;
		_bestSoFar = rhs._bestSoFar;
		_filename = rhs._filename;
		_filename_states = rhs._filename_states;	
	};
	
	void set_model_file(std::string model_file)
	{
		_model_file = model_file;
	};
	
	int objectiveFunc(  const SimTK::Vector &newControls, bool new_coefficients, SimTK::Real& f ) const; 
	double evaluate_cost_function(bool newCoefficients) const;
	void storing_excitations(const SimTK::Vector& excitations)  const;
	void storing_states(const SimTK::State& s) const;
	void store_hand_kinematics(const SimTK::State& s) const;
	
	void set_spasticity(bool spasticity, double G, double velocity_threshold)
	{
		_isSpasticity = spasticity;
		_G = G;
		_velocity_threshold = velocity_threshold;
	};
	
	void set_low_dim_params(int nbr_pcs, cv::Mat eigvectors, cv::Mat mean)
	{
		_eigvectors = eigvectors.clone();
		_mean = mean.clone();
		_isLowDim = true;
		_nbr_pcs = nbr_pcs;
	};

	void set_isFinal(bool isFinal)
	{
		_isFinal = isFinal;
	};

	SimTK::State get_final_state()
	{
		return _sf;	
	};
	
	void set_bestSoFar(double bestSoFar)
	{
		_bestSoFar = bestSoFar;
	};

	void set_initial_final_time(double initialTime, double finalTime)
	{
		_initialTime = initialTime;
		_finalTime = finalTime;
	};

	void set_rand_end_eff(std::vector<double> rand_end_eff)
	{
		_rand_end_eff = rand_end_eff;
	};
			
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
