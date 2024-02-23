#include "OptSystem.h"
#include "ArmController.h"
#include "CMAOptimizer.h"
#include "TrajOpt.h"
#include <mpi.h>
#include <math.h>
#include <string>
#include <sstream>
#include <chrono>
#include <random>
#include <opencv2/opencv.hpp>

#define DIETTAG 2

using namespace OpenSim;
using namespace SimTK;

void ArmController::computeControls(const SimTK::State& s, SimTK::Vector &controls) const
{
		double t = s.getTime();
		int phase;

		for (phase = 0; phase < _timeV.size(); phase++) {
			if (t <= _timeV[phase]) {
				break;
			}
		}
			
		Vector muscleControl(1);
		// const auto& socket = getSocket<Actuator>("actuators");
		auto actuators = getActuatorSet();
		for (int i = 0; i < _nbr_actuators; i++) {
			// Muscle* muscle = static_cast<Muscle*>(&getActuators().get(i));
			// auto muscle = dynamic_cast<const Muscle*>(&socket.getConnectee(i));
			auto muscle = static_cast<const Muscle*>(&actuators.get(i));
			
			muscleControl[0] =  _knotsV[phase * _nbr_actuators + i];
					
			// We twice check the range limit of muscleControl for spasticity to have an effect regardless of the optimization. 
			if (muscleControl[0] < 0.01) {
				muscleControl[0] = 0.01;
			}
			else if (muscleControl[0] > 0.99) {
				muscleControl[0] = 0.99;
			}
									
			muscle->setControls(muscleControl, controls);
			_last_excitations[i] = muscleControl[0];
		}
	
		_control_index++;
}

void ArmController::set_control_index(int control_index)
{
	_control_index = control_index;
}

void ArmController::reset_control_index()
{
	_control_index = 0;
	_last_excitations = SimTK::Vector(_nbr_actuators, 0.);
}

void ArmController::set_store_states(bool store_states)
{
	_store_states = store_states;
}
void ArmController::set_store_excitations(bool store_excitations)
{
	_store_excitations = store_excitations;
}
void ArmController::set_filename(std::string filename)
{
	_filename = filename;
	_filename_excitations = _filename + "_excitations.txt";
}

int main(int argc, char* argv[])
{
	
	// std::string model_file = "ue_rigid.osim";
	// std::string model_file = "MOBL_ARMS_fixed_41.osim";
	std::string model_file = "/home/riley/repos/598/project/UE_source_code/MOBL_ARMS_fixed_41.osim";
	// OpenSim::Model osimModel("MOBL_ARMS_fixed_41.osim");
	OpenSim::Model osimModel(model_file);
	try{				
		int rc = MPI_Init(&argc, &argv);
		if (rc != MPI_SUCCESS) {
			printf ("Error starting MPI program. Terminating.\n");
			MPI_Abort(MPI_COMM_WORLD, rc);
		}
		int numtasks = 0;
		int myrank = 0;
		MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
		MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	
		std::cout << "num tasks: " << numtasks << std::endl;
	
		std::vector<double> initial_coords(7);
		for(int i = 0; i < initial_coords.size(); i++){
			initial_coords[i] = 0;
		}	
	
		//Example 1: Initial Pose.
		initial_coords[0] = 65 * 3.14159/180;	

		/*
 		// Example 2: Initial Pose for center-out reachinging ask
		initial_coords[0] = 40.0 * 3.14159/180;	
		initial_coords[1] = 30 * 3.14159/180;	
		initial_coords[2] = 10 * 3.14159/180;	
		initial_coords[3] = 80 * 3.14159/180;	
		initial_coords[4] = 90 * 3.14159/180;			
		*/		
		
		std::vector<double> end_eff_target(3);		
		// Example 1
		end_eff_target[0] = 0.4;
		end_eff_target[1] = -0.15; 
		end_eff_target[2] = 0.25;																						   		
		/*
 		// Example 2
		end_eff_target[0] = 0.55;
		end_eff_target[1] = 0.0; 
		end_eff_target[2] = 0.25;
		*/
																						   		
		/*
 		// Example 3
 		// Since this target is further away, better set the optimization time to 0.55s, i.e., topt_full_dim.set_opt_time(0.55).
		end_eff_target[0] = 0.3;
		end_eff_target[1] = 0.4; 
		end_eff_target[2] = 0.0;
		*/
																					   			
		int maxIter = 650; 
		int lambda = 32;
	
		TrajOpt topt_full_dim(model_file);
		topt_full_dim.set_app_filename("ue_full");
		topt_full_dim.set_maxIter(maxIter);
		topt_full_dim.set_lambda(lambda);
		topt_full_dim.set_random_initial_state_end_eff(initial_coords, end_eff_target);	
		topt_full_dim.set_opt_time(0.45);	
	
		topt_full_dim.set_STEP(1);
		topt_full_dim.set_start_opt(1);
		
		topt_full_dim.set_lower_upper_bound(-10.0,10.0);		
		topt_full_dim.traj_opt_setup();					
	
		for (int rank = 1; rank < numtasks; ++rank) {
			MPI_Send(0, 0, MPI_INT, rank, DIETTAG, MPI_COMM_WORLD);
		}
		MPI_Finalize();
				
		return 0;
	}
	catch (const std::exception& ex)
   	{
    	    std::cout << ex.what() << std::endl;
   	    return 1;
	}
}

