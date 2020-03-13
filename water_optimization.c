/*
 * Epanet search algorithm for obtaining valve/pumping schedules
 * Author : Dev K Arora
 * Reference : Nidhin Koshy, EPANET, Amit Sharma, Lovelesh Patel
 * 
 * 
 * |-------|    |--------|     |--------|
 * | QUEUE |--->| ENGINE |---->| OUTPUT |
 * |-------|    |--------|     |--------|
 *     ^            |              |
 *     |            V              |
 *     |     |-----------|         |
 *     |-----| MINIMISER |<--------|
 *           |-----------|
 * 
 * Queue: ENReadOutFlow();Reads outflow demand of each tank and stores it in the tankcontrol structure array. Reads the outflow demand for each tank for times 0 to timeperiod-1.
 *        Initiate_valve();
 * 	  Queueing_Engine();Reads the jobs in the file for each tank demands and valve settings for times 0 to timeperiod-1.
 *        Job_Handler();Handles the jobs from the queue.
 * 
 * Engine: update_tank_level(); Updates the tank levels in the tank structure
 *         compute_flows();Computes the flow across valves and tanks using EPANET
 *         
 * Output: display_output();Displays the final output of Tank and Valve structures.
 * 
 * Minimiser: ENOptimiseValve();Main optimsation module. Optimises the flow across all the pipes and tanks.
 *            objective_function();Calculates the total penalty of a particular simulation
 *            compute_gradient();Computes the gradient which decides the amount of valve opening for a particular scenario
 *            update_control();Updates the tank and valve flow values according to the EPANET simulation.
 *            norm_gradient();Normalize the Tank and Valve gradient structures.norm_gradient is just the sum of individual gradient values
 *            random_direction();Defines a random vector in a particular direction.
 * 
 * Queue schedules the jobs and loads into the Priority Queue. 
 * 
 * Engine takes in these jobs to start the simulation.
 * There is an external clock that determines the simulation time.
 * After every simulation, tank and valve structures are updated.
 * 
 * Output displayes the tank and valve structures in readable format.
 * 
 * Minimiser takes these tank and valve settings from engine and output,
 * calculates the objective function and finds the gradient that minimises the gradient.
 * The nedded changes are then updated as jobs in the Priority Queue and the cycle continues.
 * 
 */


// Main function starts here
/*
Type :-
1. Search
Search NetworkFilePath Duration ExternalInputFilePath StartTime EndTime
2. Target
Target NetworkFilePath Duration ExternalInputFilePath StartTime EndTime CurrentTankLevelFilePath SolutionFilePath
3. SIM
SIM NetworkFilePath Duration ExternalInputFilePath StartTime EndTime SolutionFilePath
4. Target with advance
Target NetworkFilePath Duration ExternalInputFilePath StartTime EndTime CurrentTankLevelFilePath SolutionFilePath W12[1 - 10] W3[100 - 150] W4[300 - 500]
*/


