#ifndef ARMCONTROLLER_H
#define ARMCONTROLLER_H

#include <OpenSim/OpenSim.h>
#include <sstream>
#include <opencv2/opencv.hpp>
 	
#define SSTR( x ) static_cast< std::ostringstream & >(( std::ostringstream() << std::dec << x ) ).str()

namespace OpenSim {

class ArmController : public Controller { 

OpenSim_DECLARE_CONCRETE_OBJECT(ArmController, Controller);

private:
	SimTK::Vector _knotsV;
	SimTK::Vector _timeV;
	mutable int _control_index;
	bool _store_excitations;
	bool _store_states;
	std::string _filename;
	std::string _filename_excitations;
	std::string _filename_states;
	cv::Mat _eigvectors;
	cv::Mat _mean;
	bool _isLowDim;
	int _nbr_pcs;	
	mutable std::vector<double> _fiber_velocity_vector;
	mutable SimTK::Vector _last_excitations;
	mutable SimTK::Vector _last_control;
	bool _isSpasticity;
	double _G;
	double _velocity_threshold;	
	bool _isFinal;
	int _nbr_actuators;
	int _dim_reduction;

public:
	ArmController(const SimTK::Vector knotsV, const SimTK::Vector timeV) : Controller()
        {
        	_knotsV = knotsV;
                _timeV = timeV;
		_control_index = 0;
		_store_excitations = false;
		_store_states = false;
		_filename = "";
		_filename_excitations = "";
		_filename_states = "";
		_isLowDim = false;
		_isSpasticity = false;
		_G = 1.;
		_velocity_threshold = 0.;
		_isFinal = false;
        };

	void set_low_dim_params(int nbr_pcs, cv::Mat eigvectors, cv::Mat mean)
	{
		_eigvectors = eigvectors.clone();
		_mean = mean.clone();
		_isLowDim = true;
		_nbr_pcs = nbr_pcs;
		_last_control = SimTK::Vector( nbr_pcs, 0.);
	};

	void set_isFinal(bool isFinal)
	{
		_isFinal = isFinal;
	};

	void set_dim_reduction(int dim_reduction)
	{
		_dim_reduction = dim_reduction;
	};

	void computeControls(const SimTK::State& s, SimTK::Vector &controls) const override;

	void update_time(SimTK::Vector timeV)
	{
		_timeV = timeV;
		_isFinal = false;
	};
	
	void update_knots(SimTK::Vector knotsV)
	{
		_knotsV = knotsV;
		_isFinal = false;
		_nbr_actuators = getModel().getActuators().getSize();
		_last_excitations = SimTK::Vector(_nbr_actuators, 0.);
	};

	void set_control_index(int control_index);
	void reset_control_index();

	void set_store_excitations(bool store_excitations);
	void set_store_states(bool store_states);
	void set_filename(std::string filename);

	void set_spasticity(bool spasticity, double G, double velocity_threshold)
	{
		_isSpasticity = spasticity;
		_G = G;
		_velocity_threshold = velocity_threshold;
	};
	
	void set_fiber_velocity_vector(std::vector<double> fiber_velocity_vector)
	{
		_fiber_velocity_vector = fiber_velocity_vector;
	};

	SimTK::Vector& get_last_excitations()
	{
		return _last_excitations;
	};
	
	SimTK::Vector get_last_control()
	{
		return _last_control;
	};

	std::string get_filename(){
		return _filename;
	};
	
	std::string get_filename_excitations(){
		return _filename_excitations;
	};
};
};

#endif
