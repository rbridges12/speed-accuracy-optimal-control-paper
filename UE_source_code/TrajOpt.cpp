#include "TrajOpt.h"
#include "ArmController.h"
#include "CMAOptimizer.h"
#include "OptSystem.h"
#include <mpi.h>
#include <math.h>
#include <string>
#include <sstream>
#include <chrono>
#include <random>
#include <memory>

bool FIRST = true;
int CURR_STEP = 0;
double lower_bound_flag = -100; //if less than lower bound, ignore joint 

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm;
        stm << n;
        return stm.str();
    }
}

void TrajOpt::traj_opt_setup()
{		
		int numtasks = 0;
		int myrank = 0;
		MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
		MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
		
		OpenSim::Model osimModel(_model_file.c_str());	
		SimTK::State& si = osimModel.initSystem();
		const OpenSim::CoordinateSet& coords = osimModel.getCoordinateSet();

		//double est_opt_time = estimate_optimization_time(si, osimModel);
		const int KNOTS = round(_est_opt_time/0.1);
		
		SimTK::Vector timeV(KNOTS);
		double initialTime = si.getTime();
		double finalTime = initialTime + _est_opt_time;

		for (int i = 1; i <= KNOTS; i++) {
			timeV[i-1] = initialTime + (finalTime - initialTime) * i / KNOTS;
		}		

		if(!isLowDim){
			_numParams = osimModel.getNumControls()*KNOTS;	
		}else{
			_numParams = _nbr_pcs * KNOTS;
		}		
		SimTK::Vector controls(_numParams, 0.01);
		
		OpenSim::ArmController *controller = new OpenSim::ArmController(timeV, controls);
	
		if(isLowDim){
			controller->set_low_dim_params(_nbr_pcs, _eigvectors, _mean);
		}
		if(_isSpasticity){	
			controller->set_spasticity(true, _G, _velocity_threshold);
		}	
		controller->setActuators(osimModel.updActuators());
		controller->set_dim_reduction(_dim_reduction);
		osimModel.addController(controller);
		// coords.get(0).setValue(si, 0.0);
		
		_actuated_dof = actuated_degrees_of_freedom(coords);

		if(myrank == 0){
			if(_randomInitial){
			
				// _rand_coords = set_random_state(coords, si);
			
				for (int nproc = 1; nproc < numtasks ; nproc++){
					MPI_Send(&_rand_coords[0], _rand_coords.size(), MPI_DOUBLE, nproc, 0, MPI_COMM_WORLD);
				}
			
			}else{
				// set_state(coords, si, _rand_coords);
			}
		
			if(_randomTarget){
				_rand_end_eff = set_random_end_eff_target();
				
				for (int nproc = 1; nproc < numtasks ; nproc++){
					MPI_Send(&_rand_end_eff[0], _rand_end_eff.size(), MPI_DOUBLE, nproc, 0, MPI_COMM_WORLD);
				}
			}
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		if(myrank > 0){
			if(_randomInitial){
				_rand_coords.resize(_actuated_dof);
				MPI_Status status;
				MPI_Recv(&_rand_coords[0], _actuated_dof, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);	
			}	

			if(_randomTarget){
				MPI_Status status;
				MPI_Recv(&_rand_end_eff[0], _rand_end_eff.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
			}

			// set_state(coords, si, _rand_coords);
		}	
		
		const OpenSim::Set<OpenSim::Muscle> &muscleSet = osimModel.getMuscles();
     		for(int i=0; i< muscleSet.getSize(); i++ ){
            		// muscleSet[i].setActivation(si, 0.01);
		}
		osimModel.equilibrateMuscles(si);
		si.setTime(initialTime);		
		SimTK::State si_n = si;

		for(int step = 0; step < _STEP; step++){
			
			std::string filename = _app_filename + patch::to_string(step);
			std::string filename_states = _app_filename;
			
			si_n = traj_opt(si_n, osimModel, si_n.getTime(), filename, filename_states);				
	
			if(fabs(initialTime - si_n.getTime()) < 1e-6){
				break;
			}
			
			initialTime = si_n.getTime();
			_curr_step = _curr_step + 1;	
		}
	
	osimModel.removeController(controller);
}

/* MOCK FUNCTION */
std::vector<double> TrajOpt::set_random_end_eff_target()
{
	std::vector<double> end_eff_target(3);
		
	// HARD CODED.
	end_eff_target[0] = 0.55;
	end_eff_target[1] = 0.;
	end_eff_target[2] = 0.25;

	return end_eff_target;
}


std::vector<double> TrajOpt::set_random_state(const OpenSim::CoordinateSet& coords, SimTK::State& si)
{
	std::random_device rd;
    	std::mt19937 e2(rd());
    	std::uniform_real_distribution<> dist(0, 1);

	std::vector<double> random_coords;

	for(int i = 0; i < coords.getSize(); i++){
	
		double lower_bound = coords.get(i).getRangeMin(); 
		double upper_bound = coords.get(i).getRangeMax(); 
	
		if(coords.get(i).getRangeMin() > lower_bound_flag){
			double random_value = lower_bound + dist(e2)*(upper_bound - lower_bound);			
			coords.get(i).setValue(si, random_value);
			random_coords.push_back(random_value);					
		}
	}

	return random_coords;
}

void TrajOpt::set_state(const OpenSim::CoordinateSet& coords, SimTK::State& si, std::vector<double> random_coords)
{
	int j = 0;

	for(int i = 0; i < coords.getSize(); i++){
		if(coords.get(i).getRangeMin() > lower_bound_flag){
			coords.get(i).setValue(si, random_coords[j]);
			j++;
		}
	}
}

int TrajOpt::actuated_degrees_of_freedom(const OpenSim::CoordinateSet& coords)
{
	int actuated_dof = 0;
	
	for(int i = 0; i < coords.getSize(); i++){
	
		double lower_bound = coords.get(i).getRangeMin(); 
		double upper_bound = coords.get(i).getRangeMax(); 
	
		if(coords.get(i).getRangeMin() > lower_bound_flag){
			actuated_dof = actuated_dof + 1;
		}
	}
	return actuated_dof;
}

/* MOCK FUNCTION */
double TrajOpt::estimate_optimization_time(const SimTK::State& si, OpenSim::Model& osimModel)
{
	// HARD-CODED.	
	double estimated_time = 45;
	return estimated_time;	
}

SimTK::State TrajOpt::traj_opt(const SimTK::State& si, OpenSim::Model& osimModel, double initialTime, std::string filename, std::string filename_states){
	try{	
		int myrank, numtasks;
		numtasks = 0;
		myrank = 0;
		MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
		MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
		printf ("Number of processes= %d My rank= %d\n", numtasks, myrank);
				
		if(_est_opt_time < 1e-4){
			_est_opt_time = 0.45;
		}

		std::cout << "est_opt_time: " << _est_opt_time << std::endl;
				
		if(_est_opt_time < 0.05){
			return  si;
		}

		double min_feedback = SimTK::Infinity;
		if(_curr_step < _STEP){	
			min_feedback = 0.5; 
		}	
		double window_interval = std::min(min_feedback, _est_opt_time);  		

		double finalTime = initialTime + _est_opt_time;	

		std::cout << "final time: " << finalTime << std::endl;		
	
		const int KNOTS = round(_est_opt_time/0.1);

		if(!isLowDim){
			_numParams = osimModel.getNumControls()*KNOTS;	
		}else{
			_numParams = _nbr_pcs * KNOTS;
		}		

		SimTK::Vector timeV(KNOTS);

		// fill time vector with discretized time points (this is basically calling linspace)
		for (int i = 1; i <= KNOTS; i++) {
			timeV[i-1] = initialTime + (finalTime - initialTime) * i / KNOTS;
		}	
		
		double initialGuess = 0.1;
		_controls = SimTK::Vector(_numParams, initialGuess);	
	
		const OpenSim::Set<OpenSim::Muscle> &muscleSet = osimModel.getMuscles();
		int nbr_muscles = muscleSet.getSize();		

		for (int j = 0; j < _numParams; j++){
			double activation = muscleSet[j % nbr_muscles].getActivation(si);
			_controls[j] = activation;
		}

		OpenSim::ArmController* controller = static_cast<OpenSim::ArmController*>(&osimModel.getControllerSet()[0]);
		controller->update_time(timeV);
		controller->update_knots(_controls);		
		controller->reset_control_index();

		// SimTK::SemiExplicitEuler2Integrator integrator(osimModel.getMultibodySystem());
		// integrator.setAccuracy(1.0e-2);
			
		// SimTK::TimeStepper ts(osimModel.getMultibodySystem(), integrator);
        // 	ts.initialize(si);
        // 	ts.setReportAllSignificantStates(true);
        // 	integrator.setReturnEveryInternalStep(true);
		
		// OpenSim::Manager manager(osimModel, integrator);	
		// manager.setInitialTime(initialTime);
		// manager.setFinalTime(finalTime);	

		OpenSim::Manager manager(osimModel);
		manager.setIntegratorAccuracy(1.0e-2);
	
		OptSystem optSystem(_numParams, si, osimModel, manager, filename, filename_states, finalTime, window_interval, KNOTS);
		optSystem.set_initial_final_time(initialTime, finalTime);
		optSystem.set_rand_end_eff(_rand_end_eff);		
		optSystem.set_model_file(_model_file);		
	
		double f = SimTK::NaN;
		
		std::string op_res_suffix = "_optimization_result";
		std::string filename_op_res = filename + op_res_suffix; 

		if(_curr_step < _start_opt){	
			if(std::ifstream(filename_op_res.c_str())){
				std::ifstream myfile(filename_op_res.c_str());
				double a;
				for (int i = 0; i < _numParams; i++) {
		     			myfile >> a;			
		     			_controls[i] = a;
				}
				myfile.close(); 
			}				
		}		

		optSystem.objectiveFunc(_controls, true, _bestSoFar);
		optSystem.set_bestSoFar(_bestSoFar);	
		
		SimTK::Vector lower_bounds(_numParams, _lower_bound);
		SimTK::Vector upper_bounds(_numParams, _upper_bound);
	
		optSystem.setParameterLimits(lower_bounds, upper_bounds);

		SimTK::CMAOptimizer opt(optSystem);
		opt.set_filename_op_result(filename_op_res);
			
		opt.setMaxIterations(_maxIter);
		opt.setAdvancedIntOption("lambda", _lambda);
		
		opt.setAdvancedBoolOption("enableMPI", true);
		optSystem.set_isFinal(false);
		opt.setAdvancedRealOption("init_stepsize", _init_stepsize);
		
		try{
			if(_curr_step >= _start_opt){
				f = opt.optimize(_controls);
			}else{
				f = _bestSoFar;	
			}		
		}
		catch (const std::exception& ex){
			std::cout << ex.what() << std::endl;
		}
		
		if (myrank == 0){
			
			std::cout << "min result: " << f << std::endl;	
			
			if(_NOISE > 0){	
		
				bool SucInt = false;
			
				while(!SucInt){
					std::cout << "Trying to integrate after adding noise..." << std::endl;	
					std::vector<double> white_noise = compute_white_noise(_numParams, _NOISE);	
					SimTK::Vector new_controls(_numParams);		
	
					for(int i = 0; i < _numParams; i++){
						new_controls[i] = _controls[i] + white_noise[i];
					}
			
					double integrateVal = 0;
					optSystem.objectiveFunc(new_controls, true, integrateVal);
				
					if(!isinf(integrateVal)){
					
						SucInt = true;
				
						for(int i = 0; i < _numParams; i++){
							_controls[i] = new_controls[i];
						}
					}
					std::cout << "min result (after noise): " << integrateVal << std::endl;	
				}		
			}
		
			std::ofstream ofile;
			ofile.open(filename_op_res.c_str());
			
			for (int i = 0; i < _controls.size(); ++i) {
				ofile << std::fixed << std::setprecision(std::numeric_limits<double>::digits10+10) << _controls[i] << std::endl;
			}
			ofile.close();
		
			optSystem.set_isFinal(true);
			optSystem.objectiveFunc(_controls, false, _bestSoFar);	
		}

		MPI_Barrier(MPI_COMM_WORLD);	
			
		if(myrank == 0){
			for (int nproc = 1; nproc < numtasks ; nproc++){
				MPI_Send(&_controls[0], _controls.size(), MPI_DOUBLE, nproc, 0, MPI_COMM_WORLD);
			}
		}else{
			MPI_Status status;
			MPI_Recv(&_controls[0], _controls.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);	
		}
			
		optSystem.set_isFinal(true);
		optSystem.objectiveFunc(_controls, true, _bestSoFar);
		_control_list = optSystem.get_control_list();		
		_state_list = optSystem.get_state_list();

		MPI_Barrier(MPI_COMM_WORLD);	
			
		SimTK::State sf = optSystem.get_final_state();
			
		return sf;
	}
	catch (const std::exception& ex)
   	{
    	    std::cout << ex.what() << std::endl;
	}
}

std::vector<double> TrajOpt::compute_white_noise(int dimension, double standard_deviation)
{
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	std::default_random_engine generator (seed);

  	std::normal_distribution<double> distribution (0.0, standard_deviation);
	std::vector<double> white_noise(dimension); 

  	for (int i=0; i < dimension; i++){
    		white_noise[i] = distribution(generator);
	}
	
	return white_noise;	
}



