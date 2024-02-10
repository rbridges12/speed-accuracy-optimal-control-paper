#include "OptSystem.h"
#include "ArmController.h"
#include <string>
#include <mpi.h>
#include <list>

int OptSystem::objectiveFunc(  const SimTK::Vector &newControls, bool new_coefficients, SimTK::Real& f ) const
{
    osimModel.initSystem();
    OpenSim::ArmController* controller = static_cast<OpenSim::ArmController*>(&osimModel.getControllerSet()[0]);
    controller->update_knots(newControls);
    controller->reset_control_index();
    
    try{
        if(!new_coefficients){ // Visualization
            SimTK::State sv = si;
            int prank;
            MPI_Comm_rank(MPI_COMM_WORLD, &prank);
            if (prank == 0){
                f = evaluate_cost_function(new_coefficients);
            }
        }else{
            _manager.setWriteToStorage(false);
            _manager.setPerformAnalyses(false);
            f = evaluate_cost_function(new_coefficients);
        }
    }
    catch (const std::exception& ex)
    {
        std::cout << ex.what() << std::endl;
        std::cout << "integration errors" << std::endl;
        f = SimTK::Infinity;
    }
    
    return 0;
};

double OptSystem::evaluate_cost_function(bool new_coefficients) const
{
    int prank;
    MPI_Comm_rank(MPI_COMM_WORLD, &prank);
    
    int nbr_timesteps = 0;
    double activation_cost = 0;
    const OpenSim::Set<OpenSim::Muscle> &muscleSet = osimModel.getMuscles();
    int nbr_muscles = muscleSet.getSize();
    
    std::vector<double> fiber_velocity(nbr_muscles, 0.);
     
    double jointLimits_cost = 0;
    OpenSim::CoordinateSet& modelCoordinateSet = osimModel.updCoordinateSet();
    SimTK::State st = si;
    
    _ts.initialize(st);
    
    std::ofstream ofs;
    std::ofstream of_hand;
    if(!new_coefficients){
        if (prank == 0){
            std::string bestStates = "_bestSoFar_states.sto";
            std::string filename_bestStates = _filename_states + bestStates;
            if(_initialTime < 1e-10){
                ofs.open(filename_bestStates.c_str(), std::ofstream::out);
                ofs << "version=1" << "\n";
                ofs << "nRows=100000" << "\n";
                int nColumns = osimModel.getNumStateVariables() + 1;
                ofs << "nColumns=" << nColumns << "\n";
                ofs << "inDegrees=no" << "\n";
                ofs << "endheader" << "\n";
                ofs << "time ";
                OpenSim::Array<std::string> var_names = osimModel.getStateVariableNames();
                for(int k = 0; k < var_names.size(); k++){
                    ofs << var_names.get(k) << " ";
                }
                ofs << "\n";
            }else{
                ofs.open(filename_bestStates.c_str(), std::ofstream::app);
            }
		of_hand.open(_filename_hand_kinematics.c_str(), std::ofstream::out);
        }
    }
    
    OpenSim::ArmController* armController = static_cast<OpenSim::ArmController*>(&osimModel.getControllerSet()[0]);
     
    bool passedOverlap = false;
    double initialTime = _ts.getState().getTime();
    double td_time_reflex_delay = 0.03;
     
    std::list<SimTK::Vector> list_fiber_velocities;
    double store_time_interval = 0.001;
    int store_size = td_time_reflex_delay/store_time_interval;
    double prev_time = 0;
    
    _state_list.clear();
    _control_list.clear();
    
    SimTK::Vec3 massCenter;
    SimTK::Vec3 position, velocity;
    osimModel.getBodySet().get("hand").getMassCenter(massCenter);
     
    double task_term_accuracy = 0.0025; 
    double min_time = 10e4;

    int curr_interval = 0;
    int prev_interval = -1;

    while (_ts.getState().getTime() < _finalTime) {
        
        nbr_timesteps++;
        _ts.stepTo(_finalTime);
        st = _ts.getState();
        
        osimModel.getMultibodySystem().realize(st, SimTK::Stage::Report);
       
        for(int i=0; i < nbr_muscles; i++ ){
            activation_cost = activation_cost + muscleSet[i].getActivation(st);
	    fiber_velocity[i] = muscleSet[i].getFiberVelocity(st);	
	}
       
	armController->set_fiber_velocity_vector(fiber_velocity);
        
        double jointValue, minValue, maxValue;
        for(int i=0; i < modelCoordinateSet.getSize(); i++){
            jointValue = modelCoordinateSet[i].getValue(st);
            minValue = modelCoordinateSet[i].getRangeMin();
            maxValue = modelCoordinateSet[i].getRangeMax();
            
            if (jointValue < minValue){
                jointLimits_cost = jointLimits_cost + (minValue - jointValue);
            }
            if(jointValue > maxValue){
                jointLimits_cost = jointLimits_cost + (jointValue - maxValue);
            }
        }
        
        if (st.getTime() < (initialTime + _window_interval)){
            _sf = st;
        }else{
            passedOverlap = true;
        }
        
        if(_isFinal){
            SimTK::Vector last_control = armController->get_last_control();
            _control_list.push_back(last_control);
            _state_list.push_back(st);
        }
        
        if(!new_coefficients){
            if(prank == 0){
                if(passedOverlap){
                    break;
                }else{
                    ofs << st.getTime() << " ";
                    SimTK::Vector state_values = osimModel.getStateValues(st);
                    
                    for(int k = 0; k < state_values.size(); k++){
                        ofs << state_values.get(k) << " ";
                    }
                    ofs << "\n";

		   osimModel.getSimbodyEngine().getPosition(st, osimModel.getBodySet().get("hand"), massCenter, position);	
		    for(int k = 0; k < 3; k++){
			of_hand << position[k] << "\n";
		    }	

                }
            }
        }
    }
 
    double task_cost = 0;
    osimModel.getMultibodySystem().realize(st, SimTK::Stage::Velocity);
    
    osimModel.getMultibodySystem().realize(st, SimTK::Stage::Velocity);
    osimModel.getSimbodyEngine().getPosition(st,osimModel.getBodySet().get("hand"), massCenter, position);
    osimModel.getSimbodyEngine().getVelocity(st,osimModel.getBodySet().get("hand"), massCenter, velocity);
    
    double pro_sup_angle = modelCoordinateSet[15].getValue(st);
    double pro_sup_diff = pow(pro_sup_angle - 1.57,2);
    double wrist_diff = pow(modelCoordinateSet[16].getValue(st),2) + pow(modelCoordinateSet[17].getValue(st),2);
    
    double task_term = pow(position[0] - _rand_end_eff[0], 2) + pow(position[1] - _rand_end_eff[1], 2) + pow(position[2] - _rand_end_eff[2],2);
    double velocity_term = velocity.norm();
   
    double task_weight = 5;
    double velocity_weight = 0.5;
    double time_interval = _finalTime - initialTime;
    double activation_weight = 0.5e-2/nbr_muscles * 1.0/time_interval;
    double jointLimits_weight = 1e-5; 
    // Put weights on these terms if you want to reach the final pose with the forearm in pronation
    double pro_sup_weight = 0;
    double wrist_weight = 0;
   
    // To make time a variable in the optimization, you can determine if the task is achieved inside the loop above at each timestep,
    // and break from the loop if it does after updating the min_time variable. Then, include a cost of time objective to variable f, e.g.: 
    //double cost_of_time_weight = 100;
    //double cost_of_time = min_time;
    
    std::cout << "task_term: " << task_weight * task_term << std::endl;
    std::cout << "activation cost: " <<  activation_weight * activation_cost << std::endl;
    std::cout << "joint limits costs: " << jointLimits_weight * jointLimits_cost << std::endl;
    std::cout << "velocity term: " << velocity_weight * velocity_term << std::endl;
    std::cout << "hand position: " << position[0] << " and " << position[1] << " and " << position[2] << std::endl;
    std::cout << "target position: " << _rand_end_eff[0] << " and " << _rand_end_eff[1] << " and " << _rand_end_eff[2] << std::endl;
  	
    if(task_term > task_term_accuracy){
	activation_cost = 10e8;
        velocity_term = 10e4;
    }else{
	task_term = 0;
    }

    double f = task_weight * task_term + activation_weight * activation_cost + jointLimits_weight * jointLimits_cost
     	     + velocity_weight * velocity_term + pro_sup_weight * pro_sup_diff;
    
    if (isnan(f)){
        f = SimTK::Infinity;
    }
    
    std::cout << "f: " << f << " by rank " << prank << std::endl; 
    
    return f;
}

void OptSystem::storing_excitations(const SimTK::Vector& excitations)  const
{
    std::ofstream fo;
    fo.open(_filename_excitations.c_str(), std::ios_base::app);
    
    for(int i = 0; i < excitations.size(); i++){		
        fo <<  excitations[i] << std::endl;	
    }
}

void OptSystem::storing_states(const SimTK::State& s) const
{
    OpenSim::CoordinateSet& modelCoordinateSet = osimModel.updCoordinateSet();
    
    std::ofstream fo;
    fo.open(_filename_states_learning.c_str(), std::ios_base::app);
    
    for(int i = 0; i < modelCoordinateSet.getSize(); i++){
        double jointValue = modelCoordinateSet[i].getValue(s);
        fo <<  jointValue << std::endl;	
    }		
}
