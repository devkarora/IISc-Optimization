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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <math.h>
#include <float.h>                                                             //(2.00.12 - LR)
#include <time.h>
#include <stdbool.h>
#include "hash.h"
#include "text.h"
#include "types.h"
#include "enumstxt.h"
#include "funcs.h"
#define  EXTERN extern 
#include "vars.h"
#include "toolkit.h"

// Macro declarations
#define MAX_VALVEID 256
#define MAX_TIMEPERIOD 72
#define MAXTANKLEVEL 100
#define myMAXITER 30
#define myOUTER_MAXITER 1000
#define DEBUG 1
#define LPS_TANKUNITS 3.6 	/* 1LPS * 3600seconds/1000 = volume in one hour in meter^3*/
#define MAX_VALVEVALUE 350
#define FILENAME_JOBS "joblist.txt"
#define SIZE 200000	/* Size of Queue */
#define MAX_VELOCITY 2.5    /* Max velocity in m/s */
#define NUM_LITERS_PER_M3 1000
#define MAX_HOURS	25

const int PRINT_INTERVAL_FUNCTIONVALUE = 500;
const int PRINT_INTERVAL_VALVEVALUES = 500;

// Struct Declarations
struct TargetSolutionStruct {
	float tankId[MAXTANKLEVEL];
	int timeVal;
};

struct ValveSolutionStruct {
	float valveId[MAX_VALVEID];
	int timeVal;
};

typedef struct PRQ
{
	char keyword[20];
	int hours;
	int minutes;
	char id[20];
	float value;
}PriorityQ;

// Priority Queue Size
PriorityQ PQ[SIZE];
#define MAX_LINKSIZE 5
// EPANET TANK Structure
struct TankStruct {
	char *TankID;	//Tank ID
	int NodeIndex;	//NodeIndex
	int TimePeriod;	//TimePeriod
	int InputLink[MAX_LINKSIZE]; 	//LinkIndex of input Link. Currently assuming there is only one input link to a tank;
	float InFlow[MAX_TIMEPERIOD];	//Array of input flows at each time.
	int OutputLink[MAX_LINKSIZE];	//LinkIndex of output Link. Currently assuming there is only one output link to a tank.
	float OutFlow[MAX_TIMEPERIOD];	//OutFlow array of size timeperiod (T).
	float OutFlow_Epanet[MAX_TIMEPERIOD];
	float TankLevels[MAX_TIMEPERIOD+1];	//All tank levels from 0 to T. Array of size(T+1).
	float MaxTankLevel;
	float MinTankLevel;
};

// EPANET VALVE Structure
struct ValveStruct {
	char ValveID[MAX_VALVEID];	//Valve ID
	int ValveLink;	//LinkIndex
	int TimePeriod; //TimePeriod
	float ValveValues[MAX_TIMEPERIOD];	//Valve setting. Array of size timeperiod (T).
	float MaxValue;  // Max permissible flow across the valve
	float Diameter; // Diameter of the valve
	float Flow[MAX_TIMEPERIOD];  // Flow across the valve. Array of size timeperiod (T).
};

// Function declarations

void ENReadOutFlow(char *f_demand_input, struct TankStruct *tankcontrol, int timeperiod);
void ENOptimiseValve(struct TankStruct *tankcontrol, struct ValveStruct *valvecontrol);
void update_tank_level(struct TankStruct *tankcontrol);
double objective_function(struct TankStruct *tankcontrol_current,struct ValveStruct *valvecontrol_current);
void compute_flows(struct TankStruct *tankstruct, struct ValveStruct *valvestruct, int simulation_time);
void compute_gradient(struct TankStruct *tankcontrol_current, struct ValveStruct *valvecontrol_current, struct TankStruct *tankcontrol_gradient,struct ValveStruct *valvecontrol_gradient);
void update_control(struct TankStruct *tankcontrol_current, struct ValveStruct *valvecontrol_current,struct TankStruct *tankcontrol_gradient,struct ValveStruct *valvecontrol_gradient,int main_iteration_count);
float norm_gradient(struct TankStruct *tankcontrol_gradient,struct ValveStruct *valvecontrol_gradient);
float max_tank_level(int tank_struct_index);
void random_direction(float *random_vector, int dimension);
float tank_level_meter3_to_Level(float tank_volume, int tank_struct_index);
void display(struct TankStruct *tankcontrol, struct TankStruct *tankcontrol_gradient, struct ValveStruct *valvecontrol, struct ValveStruct *valvecontrol_gradient);
void Display_Output(struct TankStruct *tankcontrol, struct ValveStruct *valvecontrol);
int Queueing_Engine(char *f_job_input);
PriorityQ Qpop();
int Qfull();
int Qempty();
void Qdisplay();
void Qpush(char *keyword, int hours, int minutes, char *id, float value);
int Job_Handler(struct TankStruct *tankcontrol, struct ValveStruct *valvecontrol);
void Job_Scheduler(struct TankStruct *tankcontrol, struct ValveStruct *valvecontrol);
int feasiblity_checker(struct TankStruct *tankcontrol, struct ValveStruct *valvecontrol);
float Compute_MaxValue (float Diameter);
void ENReadCurrentFile(char *f_current_input, struct TankStruct *tankcontrol);
void ENReadSolutionFile(char *f_solution_level_input, struct TankStruct *tankcontrol);
void ENReadValveSolutionFile(char *f_valve_level_input, struct TankStruct *tankcontrol);
void compute_sim(struct TankStruct *tankstruct, struct ValveStruct *valvestruct);
void exitprocedure(struct TankStruct *tankcontrol, struct ValveStruct *valvecontrol);


// Global Variable Declarations

struct TankStruct *temp_tankstruct; //temp structs for tank gradient functions
struct ValveStruct *temp_valvestruct; // temp structs for valve gradient functions
int f=0,r=-1;       	//For Priority Queue
float epsilon = 0.05;  //Arbitary value taken
int iteration_count = 1;
int elapsed_time = 0; // Starting clock is considered at 00:00 hrs
int hourly_scheduler = 1;  // Determines the job to be scheduled in hours; if 1 = hourly jobs, if 2 = jobs scheduled for 2 hours

// targetMode = FALSE i.e. search mode is on by default
int targetMode = FALSE;
struct TargetSolutionStruct solutionTankLevel[MAX_HOURS];
struct TargetSolutionStruct currentLevel;
struct ValveSolutionStruct valveSolutionLevel[MAX_HOURS];
unsigned int res_counter=1;
unsigned int duration = 6;  // Duration of the simulation
unsigned int startTime = 0;  // Duration of the simulation
unsigned int w12 = 1;	// Cost of violating the midpoint rule; DAFAULT VALUE: 1; RANGE: 1~10
unsigned int w3 = 75;	// Cost of violating the periodicity for a TANK; DEFAULT VALUE: 75; RANGE: 50~150
unsigned int w4 = 300;	// Cost for making only important changes to a VALVE; DEFAULT VALUE: 300; RANGE: 100~500
unsigned int w5 = 1;	// Cost of violating the MAX VALVE value; DEFAULT VALUE: 1; RANGE: 1~10
float w6 = 0.5;	// Cost of having negative flows in a VALVE; DAFAULT VALUE: -0.5; RANGE: -0.5 ~ -0.1
int mode = 0;          // Mode of operation

float function_delta = 0.00001;  //Arbitary value taken
int functionValueCounter = 0;
int simDuration = 0;
// Logging file pointer
FILE *fptr = NULL;
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

int main(int argc, char *argv[])
{
	// file pointer initialisation
	char * f_epanet_input,
	     * f_demand_input,
	     * f_job_input,
	     * f_epanet_Rpt = "temp.out",
	     * f_valve_init,
	     * f_blank = " ",
	     * f_solution_level_input,
	     * f_valve_level_input,
	     * f_current_input ,
	     * temp_char;

	// Local Variable initialisation
	int timeperiod = 0;
	int temp_count = 0,
	    temo_count1= 0,
	    temp_count2 = 0,
	    temp_node1 = 0,
	    temp_node2 = 0,
	    temp_nodeIndex = 0,
	    open_valve_init_file_flag = 0;
	int simulation_time = 0; 
	int run_flag = 1;
	int *tank_struct_index; //a temporary array to store the corresponding index of a tank in Tank struct of epanet.

	struct TankStruct *tankcontrol; // initialise struct for tanks
	struct ValveStruct *valvecontrol; // initialise struct for valves

	// creating the log file
	char fileName[128];
	time_t currentTime;
	currentTime = time(NULL);
	struct tm tm = *localtime(&currentTime);
	sprintf(fileName, "../../Log/log-%d-%d-%d::%d:%d:%d.log", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
	fptr = fopen(fileName, "a+");
	if(fptr == NULL) {
		perror("ERROR in creating LOG file :: ");
		exit(-1);
	}


	sscanf(argv[1],"%d", &mode); 			// Mode : 0 - Search, 1- Target
	if(mode == 0) {

		fprintf(fptr, "\n************************************* Search Mode is on *************************************\n");
		fflush(fptr);
		targetMode = FALSE;
		f_epanet_input = argv[2]; 			// 2nd argument is EPANET input file
		f_demand_input = argv[3]; 			// 3rd argument is Tank demand data file
		f_job_input = argv[4];				// 4th argument is the joblist file
		sscanf(argv[5],"%d", &startTime); 		// 5th argument is StartTime [0-23 hours].
		sscanf(argv[6],"%d", &duration); 		// 6th argument is duartion [1-24 hours].
		sscanf(argv[7],"%d", &w12);	 		// 8th argument is w12.
		sscanf(argv[8],"%d", &w3);	 		// 9th argument is w3.
		sscanf(argv[9],"%d", &w4);	 		// 10th argument is w4.
		sscanf(argv[10],"%d", &w5);	 		// 11th argument is w5.
		sscanf(argv[11],"%f", &w6);	 		// 12th argument is w6.
		timeperiod = startTime + duration + 1;		
		// @Ram changed 25 to 24
		if(timeperiod > 24) {
			timeperiod = 24;			// Duration can't be more than 25 hours.resetiing it to 25 hours
			duration = timeperiod - startTime -1;
			fprintf(fptr, "Search Mode :: StartTime[%d] duration[%d] timeperiod[%d]\n", startTime, duration, timeperiod);
			fflush(fptr);
		}
	}
	else if(mode == 1) {

		fprintf(fptr, "\n************************************* Target Mode is on *************************************\n");
		fflush(fptr);
		targetMode = TRUE;
		f_epanet_input = argv[2]; 			// 2nd argument is EPANET input file
		f_demand_input = argv[3]; 			// 3rd argument is Tank demand data file
		f_job_input = argv[4];				// 4th argument is the joblist file
		f_solution_level_input = argv[5];		// 5th argument is solution file
		f_current_input = argv[6];			// 6th argument is current level file
		sscanf(argv[7],"%d", &duration); 		// 7th argument is duartion [1-24 hours].
		sscanf(argv[8],"%d", &w12);	 		// 8th argument is w12.
		sscanf(argv[9],"%d", &w3);	 		// 9th argument is w3.
		sscanf(argv[10],"%d", &w4);	 		// 10th argument is w4.
		sscanf(argv[11],"%d", &w5);	 		// 11th argument is w5.
		sscanf(argv[12],"%f", &w6);	 		// 12th argument is w6.
		// Reading startTime Value from current level file
		int timeVal;
		char line[1024];				// Read line by line and store in line
		FILE* stream = fopen(f_current_input, "r"); 	// Open current file for reading
		char * tok;
		char * tmp;
		int solFileRow = 0;
		if(stream == NULL)
			exit(-1);
		while ((fgets(line, 1024, stream)!=NULL)) {
			if(solFileRow == 1) {
				tmp = strdup(line);
				tok = strtok(tmp, ",");
				if(tok !=NULL) {
					timeVal = strtod(tok,NULL);
					fprintf(fptr, "Target Mode StartTime[%d]\n", timeVal);
					fflush(fptr);
				}
				break;
			}
			solFileRow = 1;
		}
		fclose(stream);
		timeperiod = timeVal + duration + 1;	
		// @Ram changed 25 to 24
		if(timeperiod > 24) {
			timeperiod = 24;			// Duration can't be more than 25 hours.resetting it to 25 hours
			duration = timeperiod - timeVal -1;
			fprintf(fptr, "TargetMode :: StartTime[%d] duration[%d] timeperiod[%d]\n", timeVal, duration, timeperiod);
			fflush(fptr);
		}
		printf("\n timeperiod = %d",timeperiod);
		printf("\n Duration = %d",duration);
		fflush(stdout);
	}
	else if(mode == 2) {

		fprintf(fptr, "\n************************************* Sim Mode is on *************************************\n");
		fflush(fptr);
		targetMode = FALSE;
		f_epanet_input = argv[2]; 			// 2nd argument is EPANET input file
		f_demand_input = argv[3]; 			// 3rd argument is Tank demand data file
		printf("\n %s",f_demand_input);
		fflush(stdout);
		f_job_input = argv[4];				// 4th argument is the joblist file
		f_solution_level_input = argv[5];		// 5th argument is solution file
		f_valve_level_input = argv[6];			// 6th argument is solution file
		timeperiod = 24;
	}
	else {
		printf("Invalid Mode \n");
		fflush(stdout);
		fprintf(fptr, "Invalid Mode \n");
		fflush(fptr);
		exit(-1);
	}


	// @Ram f_epanet_input = argv[2], f_epanet_Rpt = "temp.out" , f_blank = " "
	// @Ram f_blank is used to tell EPANET that there is no need to save the binary output

	ENopen(f_epanet_input, f_epanet_Rpt,f_blank); 				// Open Epanet input file
	ENsetreport("MESSAGES NO"); //suppresses warnings and erros
	ENsolveH(); //A dummy run of Epanet to populate all internal states.

	//Tank and Valve initialisation
	// @Ram Create an array of TankStruct, ValveStruct and tank_struct_index with pointer to the beginning of the array.
	temp_tankstruct = (struct TankStruct *)calloc(Ntanks,sizeof(struct TankStruct)); // to be used in compute_gradient
	temp_valvestruct = (struct ValveStruct *)calloc(Nvalves,sizeof(struct ValveStruct));
	tank_struct_index = (int *)calloc(Ntanks,sizeof(int));
	// @Ram Create another set of TankStruct, ValveStruct array with pointers tankcontrol and valvecontrol pointing to the 
	//	beginning of the array.
	tankcontrol = (struct TankStruct *)calloc(Ntanks,sizeof(struct TankStruct)); 	//allocate memory for tank structs.
	valvecontrol = (struct ValveStruct *) calloc(Nvalves,sizeof(struct ValveStruct)); //allocate memory for valve structs.

	// @Ram ENReadOutFlow is not an Epanet function. 
	ENReadOutFlow(f_demand_input, tankcontrol, timeperiod);   //Read tank demands (outflow) from file. Need to see if it can be integrated with Epanet.

	
	// @Ram Using the initial run of epanet, update the tank control struct variables namely MaxTankLevel, TimePeriod, MinTanklevel, NodeIndex
	// OutputLink and InputLink 
	// @Ram Also, update the tank struct index 
	//Initialise tank control struct pointers
	for(temp_count = 0; temp_count < Ntanks; temp_count++) {
		// Getting NodeID
		ENgetnodeindex(tankcontrol[temp_count].TankID, &temp_nodeIndex);
		temp_count2 = 0;
		// Checking the Index of the Tanks; temp_count2 will give the Tank Index
		while(Tank[temp_count2].Node != temp_nodeIndex) {
			temp_count2++;
		}
		tank_struct_index[temp_count]=temp_count2;
		// Finding Tank MAX Level
		tankcontrol[temp_count].MaxTankLevel = max_tank_level(temp_count2);
		fprintf(fptr, "\n Tank ID = %s, Tank Max Level (meter^3) = %f",tankcontrol[temp_count].TankID,tankcontrol[temp_count].MaxTankLevel);
		fflush(fptr);
		tankcontrol[temp_count].TimePeriod = timeperiod;
		tankcontrol[temp_count].MinTankLevel = 0;
		tankcontrol[temp_count].NodeIndex = temp_nodeIndex;
		// Initialising the array for InpputLink and OutputLink
		int linkLoop;
		for(linkLoop = 0; linkLoop< MAX_LINKSIZE; linkLoop++)
		{
			tankcontrol[temp_count].OutputLink[linkLoop] = 0;
			tankcontrol[temp_count].InputLink[linkLoop] = 0;
		}
		// Getting the Index of Pipes flowing into the TANK; InputLink = Pipe whose endnode is connected to TANK i.e. N2 = TankIndex ; OutputLink = Pipe whose start node is connected to TANK i.e. N1 = TankIndex
		int linkLoopN1 = 0, linkLoopN2 = 0, index_dummy=0;
		char id_dummy[20]; 
		for (temp_count2=0; temp_count2<Nlinks; temp_count2++) {
			if(Link[temp_count2].N1 == tankcontrol[temp_count].NodeIndex)
			{
				tankcontrol[temp_count].OutputLink[linkLoopN1]=temp_count2;
				fprintf(fptr, "\n Tank ID = %s, Tank Output Link ID[%d :: %d] = %s", tankcontrol[temp_count].TankID, linkLoopN1, tankcontrol[temp_count].OutputLink[linkLoopN1], Link[temp_count2].ID );
				index_dummy=Link[temp_count2].N2;
				ENgetnodeid(Link[temp_count2].N2,id_dummy);
			//	printf("\n Tank Id = %s Link Index = %d Link id = %s", tankcontrol[temp_count].TankID,index_dummy,id_dummy);
				fflush(fptr);
				linkLoopN1++;
			}
			if(Link[temp_count2].N2 == tankcontrol[temp_count].NodeIndex)
			{
				tankcontrol[temp_count].InputLink[linkLoopN2]=temp_count2;
				fprintf(fptr, "\n Tank ID = %s, Tank Input Link ID[%d :: %d] = %s", tankcontrol[temp_count].TankID, linkLoopN2, tankcontrol[temp_count].InputLink[linkLoopN2], Link[temp_count2].ID );
				fflush(fptr);
				linkLoopN2++;
			}
		}

	}

	if(targetMode == TRUE) {
		// Reading Solution File
		ENReadSolutionFile(f_solution_level_input, tankcontrol);

		// Reading Current File
		ENReadCurrentFile(f_current_input, tankcontrol);

		int loop = 0;
		for(; loop <Ntanks; loop++) {
			tankcontrol[loop].TankLevels[timeperiod+1] = solutionTankLevel[0].tankId[loop];
			tankcontrol[loop].TankLevels[currentLevel.timeVal] = currentLevel.tankId[loop];
		}
	}

	// Initialise Valve Controls

	for(temp_count =0; temp_count < Nvalves; temp_count++) {
		// Getting the link index
		valvecontrol[temp_count].ValveLink = Valve[temp_count+1].Link;
		// Getting the Valve ID
		ENgetlinkid(Valve[temp_count+1].Link,valvecontrol[temp_count].ValveID);
		if(DEBUG) {
			fprintf(fptr, "\n Valve Link of Valve %s is %d", valvecontrol[temp_count].ValveID, valvecontrol[temp_count].ValveLink);
			fflush(fptr);
		}

		valvecontrol[temp_count].TimePeriod = timeperiod;
		fprintf(fptr, "\n Valve Id = %s \t ",valvecontrol[temp_count].ValveID);
		fflush(fptr);

		// To get the Diameter of the valves for max flow calculation across a valve
		ENgetlinkvalue(valvecontrol[temp_count].ValveLink,EN_DIAMETER, &(valvecontrol[temp_count].Diameter));
		valvecontrol[temp_count].MaxValue = Compute_MaxValue(valvecontrol[temp_count].Diameter);
		fprintf(fptr, " Initial Valve Value = %f",valvecontrol[temp_count].ValveValues[0]);
		fprintf(fptr, "\n Max Flow across Valve %s = %f \t Diameter is %f",valvecontrol[temp_count].ValveID, valvecontrol[temp_count].MaxValue, valvecontrol[temp_count].Diameter);
		fflush(fptr);
	}

	if(mode == 2) {
		// Reading Tank Solution File
		ENReadSolutionFile(f_solution_level_input, tankcontrol);

		// Reading Valve Solution File
		ENReadValveSolutionFile(f_valve_level_input, tankcontrol);

		// calling sim main function
		compute_sim(tankcontrol, valvecontrol);
		exit(0);
	}

	// Initialize random number seed to 1.
	srand(5);
	if(targetMode != TRUE)
	{
		for(temp_count=0;temp_count<Ntanks;temp_count++)
		{
			tankcontrol[temp_count].TankLevels[0]= 0.5 * tankcontrol[temp_count].MaxTankLevel; 
		}
	}

	int valveCount = 0;					// Initial Tank Count = 0;
	float avg_outflow = 0.0, avg_target_outflow=0.0;
//if(targetMode!=TRUE)
//{
		while(valveCount < Nvalves) 
		{
			if(strncmp(valvecontrol[valveCount].ValveID,"80",2)==0)
			{
				valveCount++;
			}
			else if(strncmp(valvecontrol[valveCount].ValveID,"78",2)==0)
			{
				if(!targetMode)
				{
					temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
					for(;temp_count2<timeperiod;temp_count2++)
					{

						valvecontrol[valveCount].ValveValues[temp_count2]= 28.1;
					}
				}
				else
				{
				for(temp_count=0;temp_count<Ntanks;temp_count++)
					{
						if(strncmp(tankcontrol[temp_count].TankID,"Veeranageeri",11)==0)
						{		
							printf("\n Veeranageeri");
							avg_outflow = 0.0;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
							for(;temp_count2<timeperiod;temp_count2++)
							{
								avg_outflow += tankcontrol[temp_count].OutFlow[temp_count2];
							}
						avg_target_outflow = (tankcontrol[temp_count].TankLevels[timeperiod+1]-tankcontrol[temp_count].TankLevels[currentLevel.timeVal]) + avg_outflow;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
						//	avg_outflow = avg_outflow/(timeperiod-temp_count2);
							avg_target_outflow = avg_target_outflow/(timeperiod-temp_count2);
						//	printf("\n %f",avg_outflow);
							printf("\n avg target outflow %f",avg_target_outflow);
							printf("\n Current tank level %f",tankcontrol[temp_count].TankLevels[currentLevel.timeVal]);
							printf("\n Current tank level %f",tankcontrol[temp_count].TankLevels[currentLevel.timeVal]);
							printf("\n Current level %d",currentLevel.timeVal);
							printf("\n Maximum tank level %f",tankcontrol[temp_count].MaxTankLevel);
							printf("\n demand %f",avg_outflow);
							for(;temp_count2<timeperiod;temp_count2++)
							{
								valvecontrol[valveCount].ValveValues[temp_count2]=avg_outflow ;
							}

						}
					}
				}
				valveCount++;
			}

			else if(strncmp(valvecontrol[valveCount].ValveID,"76",2)==0)
			{
				if(!targetMode)
				{
					temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
					for(;temp_count2<timeperiod;temp_count2++)
					{
						valvecontrol[valveCount].ValveValues[temp_count2]= 14.4;
					}
				}
				else
				{	
					for(temp_count=0;temp_count<Ntanks;temp_count++)
					{
						if(strncmp(tankcontrol[temp_count].TankID,"NR_mohalla2",11)==0)
						{		
							printf("\n NR_mohalla2");
							avg_outflow = 0.0;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
							for(;temp_count2<timeperiod;temp_count2++)
							{
								avg_outflow += tankcontrol[temp_count].OutFlow[temp_count2];
							}
						avg_target_outflow = (tankcontrol[temp_count].TankLevels[timeperiod+1]-tankcontrol[temp_count].TankLevels[currentLevel.timeVal]) + avg_outflow;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
						//	avg_outflow = avg_outflow/(timeperiod-temp_count2);
							avg_target_outflow = avg_target_outflow/(timeperiod-temp_count2);
							printf("\n avg target outflow %f",avg_target_outflow);
							printf("\n Current tank level %f",tankcontrol[temp_count].TankLevels[currentLevel.timeVal]);
							printf("\n Maximum tank level %f",tankcontrol[temp_count].MaxTankLevel);
							printf("\n demand %f",avg_outflow);
							for(;temp_count2<timeperiod;temp_count2++)
							{
								valvecontrol[valveCount].ValveValues[temp_count2]=avg_outflow ;
							}

						}
					}

				}
				valveCount++;
			}

			else if(strncmp(valvecontrol[valveCount].ValveID,"74",2)==0)
			{
				if(!targetMode)
				{
					temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
					for(;temp_count2<timeperiod;temp_count2++)
					{
						valvecontrol[valveCount].ValveValues[temp_count2]= 57.0;
					}
				}
				else
				{	
					for(temp_count=0;temp_count<Ntanks;temp_count++)
					{
						if(strncmp(tankcontrol[temp_count].TankID,"NR_mohalla1",11)==0)
						{		
							printf("\n NR_mohalla1");
							avg_outflow = 0.0;
							avg_target_outflow =0.0;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
							for(;temp_count2<timeperiod;temp_count2++)
							{
								avg_outflow += tankcontrol[temp_count].OutFlow[temp_count2];
							}
						avg_target_outflow = (tankcontrol[temp_count].TankLevels[timeperiod+1]-tankcontrol[temp_count].TankLevels[currentLevel.timeVal]) + avg_outflow;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
						//	avg_outflow = avg_outflow/(timeperiod-temp_count2);
							avg_target_outflow = avg_target_outflow/(timeperiod-temp_count2);
							printf("\n avg target outflow %f",avg_target_outflow);
							printf("\n Current tank level %f",tankcontrol[temp_count].TankLevels[currentLevel.timeVal]);
							printf("\n Maximum tank level %f",tankcontrol[temp_count].MaxTankLevel);
							printf("\n Demand %f",avg_outflow);
							//printf("\n %f",avg_outflow);
							for(;temp_count2<timeperiod;temp_count2++)
							{
								valvecontrol[valveCount].ValveValues[temp_count2]=avg_outflow ;
							}

						}
					}

				}
				valveCount++;

			}

			else if(strncmp(valvecontrol[valveCount].ValveID,"72",2)==0)
			{
				if(!targetMode)
				{
					temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
					for(;temp_count2<timeperiod;temp_count2++)
					{
						valvecontrol[valveCount].ValveValues[temp_count2]= 79.8;
					}
				}
				else
				{	
					for(temp_count=0;temp_count<Ntanks;temp_count++)
					{
						if(strncmp(tankcontrol[temp_count].TankID,"Rajiv_nagar_I",14)==0)
						{		
							printf("\n Rajiv_nagar_I");
							avg_outflow = 0.0;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
							for(;temp_count2<timeperiod;temp_count2++)
							{
								avg_outflow += tankcontrol[temp_count].OutFlow[temp_count2];
							}
							avg_target_outflow = (tankcontrol[temp_count].TankLevels[timeperiod+1]-tankcontrol[temp_count].TankLevels[currentLevel.timeVal]) + avg_outflow;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
						//	avg_outflow = avg_outflow/(timeperiod-temp_count2);
							avg_target_outflow = avg_target_outflow/(timeperiod-temp_count2);
							//printf("\n %f",avg_outflow);
							printf("\n avg target outflow %f",avg_target_outflow);
							printf("\n Current tank level %f",tankcontrol[temp_count].TankLevels[currentLevel.timeVal]);
							printf("\n Maximum tank level %f",tankcontrol[temp_count].MaxTankLevel);
							printf("\n flows %f",avg_outflow);
							for(;temp_count2<timeperiod;temp_count2++)
							{
								valvecontrol[valveCount].ValveValues[temp_count2]=avg_outflow ;
							}
						}
					}

				}
				valveCount++;
			}

			else if(strncmp(valvecontrol[valveCount].ValveID,"70",2)==0)
			{
				if(!targetMode)
				{
					temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
					for(;temp_count2<timeperiod;temp_count2++)
					{
						valvecontrol[valveCount].ValveValues[temp_count2]= 34.0;
					}
				}
				else
				{	
					for(temp_count=0;temp_count<Ntanks;temp_count++)
					{
						if(strncmp(tankcontrol[temp_count].TankID,"Rajendra_nagar2",15)==0)
						{		
							printf("\n Rajendra_nagar2");
							avg_outflow = 0.0;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
							for(;temp_count2<timeperiod;temp_count2++)
							{
								avg_outflow += tankcontrol[temp_count].OutFlow[temp_count2];
							}
							avg_target_outflow = (tankcontrol[temp_count].TankLevels[timeperiod+1]-tankcontrol[temp_count].TankLevels[currentLevel.timeVal]) + avg_outflow;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
						//	avg_outflow = avg_outflow/(timeperiod-temp_count2);
							avg_target_outflow = avg_target_outflow/(timeperiod-temp_count2);
							//printf("\n %f",avg_outflow);
							printf("\n avg target outflow %f",avg_target_outflow);
							printf("\n Current tank level %f",tankcontrol[temp_count].TankLevels[currentLevel.timeVal]);
							printf("\n Maximum tank level %f",tankcontrol[temp_count].MaxTankLevel);
							printf("\n flows %f",avg_outflow);
							for(;temp_count2<timeperiod;temp_count2++)
							{
								valvecontrol[valveCount].ValveValues[temp_count2]=avg_outflow ;
							}
						}
					}

				}
				valveCount++;
			}

			else if(strncmp(valvecontrol[valveCount].ValveID,"68",2)==0)
			{
				if(!targetMode)
				{
					temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
					for(;temp_count2<timeperiod;temp_count2++)
					{
						valvecontrol[valveCount].ValveValues[temp_count2]= 39.6;
					}
				}
				else
				{	
					for(temp_count=0;temp_count<Ntanks;temp_count++)
					{
						if(strncmp(tankcontrol[temp_count].TankID,"Rajendra_nagar1",15)==0)
						{		
							printf("\n Rajendra_nagar1");
							avg_outflow = 0.0;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
							for(;temp_count2<timeperiod;temp_count2++)
							{
								avg_outflow += tankcontrol[temp_count].OutFlow[temp_count2];
							}
							avg_target_outflow = (tankcontrol[temp_count].TankLevels[timeperiod+1]-tankcontrol[temp_count].TankLevels[currentLevel.timeVal]) + avg_outflow;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
						//	avg_outflow = avg_outflow/(timeperiod-temp_count2);
							avg_target_outflow = avg_target_outflow/(timeperiod-temp_count2);
							//printf("\n %f",avg_outflow);
							printf("\n avg target outflow %f",avg_target_outflow);
							printf("\n Current tank level %f",tankcontrol[temp_count].TankLevels[currentLevel.timeVal]);
							printf("\n Maximum tank level %f",tankcontrol[temp_count].MaxTankLevel);
							printf("\n flows %f",avg_outflow);
							for(;temp_count2<timeperiod;temp_count2++)
							{
								valvecontrol[valveCount].ValveValues[temp_count2]=avg_outflow ;
							}
						}
					}

				}
				valveCount++;
			}


			else if(strncmp(valvecontrol[valveCount].ValveID,"66",2)==0)
			{
				if(!targetMode)
				{
					temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
					for(;temp_count2<timeperiod;temp_count2++)
					{
						valvecontrol[valveCount].ValveValues[temp_count2]= 9.9;
					}
				}
				else
				{	
					for(temp_count=0;temp_count<Ntanks;temp_count++)
					{
						if(strncmp(tankcontrol[temp_count].TankID,"Bademakan",9)==0)
						{		
							printf("\n Bademakan");
							avg_outflow = 0.0;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
							for(;temp_count2<timeperiod;temp_count2++)
							{
								avg_outflow += tankcontrol[temp_count].OutFlow[temp_count2];
							}
							avg_target_outflow = (tankcontrol[temp_count].TankLevels[timeperiod+1]-tankcontrol[temp_count].TankLevels[currentLevel.timeVal]) + avg_outflow;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
						//	avg_outflow = avg_outflow/(timeperiod-temp_count2);
							avg_target_outflow = avg_target_outflow/(timeperiod-temp_count2);
							//printf("\n %f",avg_outflow);
							printf("\n avg target outflow %f",avg_target_outflow);
							printf("\n Current tank level %f",tankcontrol[temp_count].TankLevels[currentLevel.timeVal]);
							printf("\n Maximum tank level %f",tankcontrol[temp_count].MaxTankLevel);
							printf("\n flows %f",avg_outflow);
							for(;temp_count2<timeperiod;temp_count2++)
							{
								valvecontrol[valveCount].ValveValues[temp_count2]=avg_outflow ;
							}
						}
					}

				}
				valveCount++;
			}

			else if(strncmp(valvecontrol[valveCount].ValveID,"65",2)==0)
			{
				if(!targetMode)
				{
					temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
					for(;temp_count2<timeperiod;temp_count2++)
					{
						valvecontrol[valveCount].ValveValues[temp_count2]= 14.4;
					}
				}
				else
				{	
					for(temp_count=0;temp_count<Ntanks;temp_count++)
					{
						if(strncmp(tankcontrol[temp_count].TankID,"Kesare",6)==0)
						{		
							printf("\n Kesare");
							avg_outflow = 0.0;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
							for(;temp_count2<timeperiod;temp_count2++)
							{
								avg_outflow += tankcontrol[temp_count].OutFlow[temp_count2];
							}
							avg_target_outflow = (tankcontrol[temp_count].TankLevels[timeperiod+1]-tankcontrol[temp_count].TankLevels[currentLevel.timeVal]) + avg_outflow;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
						//	avg_outflow = avg_outflow/(timeperiod-temp_count2);
							avg_target_outflow = avg_target_outflow/(timeperiod-temp_count2);
							//printf("\n %f",avg_outflow);
							printf("\n avg target outflow %f",avg_target_outflow);
							printf("\n Current tank level %f",tankcontrol[temp_count].TankLevels[currentLevel.timeVal]);
							printf("\n Maximum tank level %f",tankcontrol[temp_count].MaxTankLevel);
							printf("\n flows %f",avg_outflow);
							for(;temp_count2<timeperiod;temp_count2++)
							{
								valvecontrol[valveCount].ValveValues[temp_count2]=avg_outflow ;
							}
						}
					}

				}
				valveCount++;
			}

			else if(strncmp(valvecontrol[valveCount].ValveID,"40",2)==0)
			{
				if(!targetMode)
				{
					temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
					for(;temp_count2<timeperiod;temp_count2++)
					{
						valvecontrol[valveCount].ValveValues[temp_count2]= 37.4;
					}
				}
				else
				{	
					for(temp_count=0;temp_count<Ntanks;temp_count++)
					{
						if(strncmp(tankcontrol[temp_count].TankID,"Tippu_park",10)==0)
						{		
							printf("\n Tippu_park");
							avg_outflow = 0.0;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
							for(;temp_count2<timeperiod;temp_count2++)
							{
								avg_outflow += tankcontrol[temp_count].OutFlow[temp_count2];
							}
							avg_target_outflow = (tankcontrol[temp_count].TankLevels[timeperiod+1]-tankcontrol[temp_count].TankLevels[currentLevel.timeVal]) + avg_outflow;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
						//	avg_outflow = avg_outflow/(timeperiod-temp_count2);
							avg_target_outflow = avg_target_outflow/(timeperiod-temp_count2);
							//printf("\n %f",avg_outflow);
							printf("\n avg target outflow %f",avg_target_outflow);
							printf("\n Current tank level %f",tankcontrol[temp_count].TankLevels[currentLevel.timeVal]);
							printf("\n Maximum tank level %f",tankcontrol[temp_count].MaxTankLevel);
							printf("\n flows %f",avg_outflow);
							for(;temp_count2<timeperiod;temp_count2++)
							{
								valvecontrol[valveCount].ValveValues[temp_count2]=avg_outflow ;
							}
						}
					}

				}
				valveCount++;
			}

			else if(strncmp(valvecontrol[valveCount].ValveID,"35",2)==0)
			{
				if(!targetMode)
				{
					temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
					for(;temp_count2<timeperiod;temp_count2++)
					{
						valvecontrol[valveCount].ValveValues[temp_count2]= 19.6;
					}
				}
				else
				{	
					for(temp_count=0;temp_count<Ntanks;temp_count++)
					{
						if(strncmp(tankcontrol[temp_count].TankID,"Rajiv_nagar_II",14)==0)
						{		
							printf("\n Rajiv_nagar_II");
							avg_outflow = 0.0;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
							for(;temp_count2<timeperiod;temp_count2++)
							{
								avg_outflow += tankcontrol[temp_count].OutFlow[temp_count2];
							}
							avg_target_outflow = (tankcontrol[temp_count].TankLevels[timeperiod+1]-tankcontrol[temp_count].TankLevels[currentLevel.timeVal]) + avg_outflow;
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
						//	avg_outflow = avg_outflow/(timeperiod-temp_count2);
							avg_target_outflow = avg_target_outflow/(timeperiod-temp_count2);
							//printf("\n %f",avg_outflow);
							printf("\n avg target outflow %f",avg_target_outflow);
							printf("\n Current tank level %f",tankcontrol[temp_count].TankLevels[currentLevel.timeVal]);
							printf("\n Maximum tank level %f",tankcontrol[temp_count].MaxTankLevel);
							printf("\n flows %f",avg_outflow);
							for(;temp_count2<timeperiod;temp_count2++)
							{
								valvecontrol[valveCount].ValveValues[temp_count2]=avg_outflow ;
							}
						}
					}

				}
				valveCount++;
			}

			else if(strncmp(valvecontrol[valveCount].ValveID,"3",1)==0)
			{
				if(!targetMode)
				{
					temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
					for(;temp_count2<timeperiod;temp_count2++)
					{
						valvecontrol[valveCount].ValveValues[temp_count2]= 136.9;
					}
				}
				else
				{
					valvecontrol[valveCount].ValveValues[temp_count2]=0.0;	
					avg_outflow = 0.0;
					avg_target_outflow=0.0;
					for(temp_count=0;temp_count<Ntanks;temp_count++)
					{
						if((strncmp(tankcontrol[temp_count].TankID,"Veeranageeri",12)==0)||(strncmp(tankcontrol[temp_count].TankID,"Tippu_park",10)==0)||(strncmp(tankcontrol[temp_count].TankID,"NR_mohalla1",11)==0)||(strncmp(tankcontrol[temp_count].TankID,"NR_mohalla2",11)==0))
						{		
							printf("\n %s",tankcontrol[temp_count].TankID);
							temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
							for(;temp_count2<timeperiod;temp_count2++)
							{
								avg_outflow += tankcontrol[temp_count].OutFlow[temp_count2];
							}
								
							avg_target_outflow = avg_target_outflow + (tankcontrol[temp_count].TankLevels[timeperiod+1]-tankcontrol[temp_count].TankLevels[currentLevel.timeVal]) + avg_outflow;
						}
					}
					temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
			//		avg_outflow = avg_outflow/(timeperiod-temp_count2);
					avg_target_outflow = avg_target_outflow/(timeperiod-temp_count2);
					printf("\n avg target outflow %f",avg_target_outflow);
					//printf("\n %f",avg_outflow);
					for(;temp_count2<timeperiod;temp_count2++)
					{
						valvecontrol[valveCount].ValveValues[temp_count2]=avg_outflow ;
					}

				}
				valveCount++;
			}
			else
			{
				printf("\n Error in valve initialization");
			}

		}

//}
	/*	temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
		printf("\n %d\n %d", timeperiod,temp_count2);
		for(;temp_count2<timeperiod;temp_count2++)
		{
		printf("\n");
		for(temp_count=0;temp_count<Nvalves;temp_count++)
		{
		printf(" %f\t",valvecontrol[temp_count].ValveValues[temp_count2]);
		}
		}
	 */ 

		temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
		for(; temp_count2 < timeperiod; temp_count2++) {
			compute_flows(tankcontrol,valvecontrol,temp_count2);
		}
		int count_flag=0;
		// Run the wrapper
		while(run_flag) {
			//if(run_flag>90)
			//{
			printf("\nOUTER ITERATION Num is %d",count_flag);
			printf("\nOUTER ITERATION COUNT is %d",run_flag);
			fflush(stdout);
			fprintf(fptr, "\nOUTER ITERATION COUNT in %d",run_flag);
			fflush(fptr);
			printf("\n timeperiod = %d",timeperiod);
			printf("\n Duration = %d",duration);
			fflush(stdout);
			//}
			ENOptimiseValve(tankcontrol, valvecontrol); 
			
			
			if(run_flag>10)
			{
			if(feasiblity_checker(tankcontrol, valvecontrol)) {
				if(DEBUG) {
					printf("\n Result %u",res_counter);
					Display_Output(tankcontrol, valvecontrol);
					res_counter++;
				}
				// @Ram Commented these lines 
				Job_Scheduler(tankcontrol, valvecontrol);
				simulation_time = Job_Handler(tankcontrol, valvecontrol);
				//printf("\n sim_time1 = %d", simulation_time);
				temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
				for(; temp_count2 < timeperiod; temp_count2++) {
					compute_flows(tankcontrol,valvecontrol,temp_count2);
				}
				//compute_flows(tankcontrol, valvecontrol, simulation_time);
				update_tank_level(tankcontrol);

				if(feasiblity_checker(tankcontrol, valvecontrol)) {
					fprintf(fptr, "final output\n");
					printf("final output\n");
					fflush(stdout);
					Display_Output(tankcontrol, valvecontrol);
					break;
				}
				else
				{
					run_flag = 1;
					count_flag++;
				}
				// @Ram Comment ends here
			}
			fflush(stdout);
			}
			run_flag++;
		}
		fflush(fptr);

		// free all pointer
		for(temp_count=0; temp_count<Ntanks; temp_count++) {
			free(tankcontrol[temp_count].TankID);
		}
		free(tankcontrol);
		free(valvecontrol);
		free(temp_tankstruct);
		free(temp_valvestruct);
		fflush(fptr);
		fclose(fptr);
		ENclose();
	}

	void ENReadOutFlow(char *f_demand_input, struct TankStruct *tankcontrol, int timeperiod)
	{
		// @Ram  Reads outflow demand of each tank and stores it in the tankcontrol structure array 
		/*
		 * Reads the outflow demand for each tank for times 0 to timeperiod-1.
		 * f_demand_input	: file name of file containing outflow demand data. Each tank data should be in one line, and each field separated by a comma. First field should be the tank ID.
		 * tankcontrol	: Pointer to the TankStruct array.
		 * timeperiod	: Optimisation time period in hours.
		 */

		FILE* stream = fopen(f_demand_input, "r"); 	// Open demand file for reading
		char line[1024];				// Read line by line and store in line
		char * temp_field;				// temporary field to extract values
		int tankcount = 0;				// Initial Tank Count = 0;
		int i = 0;
		char * tok;
		char * tmp;

		while ((fgets(line, 1024, stream)!=NULL) && (tankcount < Ntanks)) {
			fprintf(fptr, "\ntankcount = %d",tankcount);
			fflush(stdout);
			fflush(fptr);
			tmp = strdup(line);
			tok = strtok(tmp, ",");
			if(tok !=NULL) {
				tankcontrol[tankcount].TankID = tok;
			}
			fprintf(fptr, "\n Tank %d: Tank ID = %s,",tankcount+1,tok);
			fflush(fptr);

			i=0;
			while((NULL != (tok = strtok(NULL, ","))) && (i<timeperiod)) {
				tankcontrol[tankcount].OutFlow[i] = strtod(tok,NULL);
				i++;
			}
			tankcount++;
		}
		fflush(fptr);
		fclose(stream);
	}

	void ENReadValveSolutionFile(char *f_valve_level_input, struct TankStruct *tankcontrol)
	{
		/*
		 * Reads the Solution InputFile for each tank for 00-23hours.
		 * f_solution_level_input 	: file name of file containing solution input file data for 24 hours.
		 */

		FILE* stream = fopen(f_valve_level_input, "r"); 	// Open Solution file for reading
		char line[1024];					// Read line by line and store in line
		int valveCount = 0;					// Initial Tank Count = 0;
		int solFileRow = 0;
		int i = 0;
		if(stream == NULL)
			exit(-1);
		char* tok;
		char * tmp;
		while ((fgets(line, 1024, stream)!=NULL) && (valveCount < Nvalves)) {
			if(solFileRow > 0 && solFileRow <= 24) {
				tmp = strdup(line);
				tok = strtok(tmp, ",");
				if(tok !=NULL) {
					valveSolutionLevel[solFileRow-1].timeVal = strtod(tok,NULL);
				}

				while((NULL != (tok = strtok(NULL, ",")))) {
					if(valveCount < Nvalves) {
						valveSolutionLevel[solFileRow-1].valveId[valveCount] = strtod(tok,NULL);
						//printf("timeval[%d] Valve Level[%f] and id[%d]\n", valveSolutionLevel[solFileRow-1].timeVal,valveSolutionLevel[solFileRow-1].valveId[valveCount], solFileRow-1);
						valveCount++;
					}
				}
			}
			solFileRow++;
		//printf("\n %d",solFileRow);
		valveCount = 0;
		}
		//@Ram added +1 next line
		simDuration = valveSolutionLevel[solFileRow-2].timeVal - valveSolutionLevel[0].timeVal + 1;
		printf("\n %d",simDuration);
		fflush(fptr);
		fclose(stream);
	}

	void ENReadSolutionFile(char *f_solution_level_input, struct TankStruct *tankcontrol)
	{
		/*
		 * Reads the Solution InputFile for each tank for 00-23hours.
		 * f_solution_level_input 	: file name of file containing solution input file data for 24 hours.
		 */

		FILE* stream = fopen(f_solution_level_input, "r"); 	// Open Solution file for reading
		char line[1024];					// Read line by line and store in line
		int tankCount = 0;					// Initial Tank Count = 0;
		int solFileRow = 0;
		int i = 0;
		if(stream == NULL)
			exit(-1);
		char* tok;
		char * tmp;
		while ((fgets(line, 1024, stream)!=NULL) && (tankCount < Ntanks)) {
			if(solFileRow > 0 && solFileRow <= 24) {
				tmp = strdup(line);
				tok = strtok(tmp, ",");
				if(tok !=NULL) {
					solutionTankLevel[solFileRow-1].timeVal = strtod(tok,NULL);
				}

				while((NULL != (tok = strtok(NULL, ",")))) {
					if(tankCount < Ntanks) {
						float temp = strtod(tok,NULL);
						solutionTankLevel[solFileRow-1].tankId[tankCount] = (temp*tankcontrol[tankCount].MaxTankLevel)/100;
						tankCount++;
					}
				}
			}
			solFileRow++;
			tankCount = 0;
		}
		fflush(fptr);
		fclose(stream);
	}

	void ENReadCurrentFile(char *f_current_input, struct TankStruct *tankcontrol)
	{
		/*
		 * Reads the current File for each tank for given times.
		 * f_current_input 	: file name of file containing target file data for a given time period.
		 */

		FILE* stream = fopen(f_current_input, "r"); 	// Open current file for reading
		if(stream == NULL)
			exit(-1);
		char line[1024];				// Read line by line and store in line
		int tankCount = 0;				// Initial Tank Count = 0;
		int solFileRow = 0;
		int i = 0;
		char* tok;
		char * tmp;

		while ((fgets(line, 1024, stream)!=NULL) && (tankCount < Ntanks)) {
			if(solFileRow == 1) {
				tmp = strdup(line);
				tok = strtok(tmp, ",");
				if(tok !=NULL) {
					currentLevel.timeVal = strtod(tok,NULL);
					fprintf(fptr, "CurrentLevel StartTime[%d]\n", currentLevel.timeVal);
					fflush(fptr);
				}

				while((NULL != (tok = strtok(NULL, ",")))) {
					if(tankCount < Ntanks) {
						float temp = strtod(tok,NULL);
						currentLevel.tankId[tankCount] = (temp*tankcontrol[tankCount].MaxTankLevel)/100;
						tankCount++;
					}
				}
			}
			solFileRow = 1;
		}
		fflush(fptr);
		fclose(stream);
	}

	void ENOptimiseValve(struct TankStruct *tankcontrol, struct ValveStruct *valvecontrol)
	{
		/*
		 * Main optimsation module. Optimises the flow across all the pipes and tanks.
		 * Main function to call
		 * tankcontrol : Tank structure compatible with EPANET
		 * valvecontrol : Valve structure compatible with EPANET
		 */

		char * f_blank = " ";
		int 	temp_count = 0,
			temp_count2 = 0;
		double  function_value_current=0,
			function_value_previous = 0,
			function_value_previous_printiter=0;

		struct TankStruct * tankcontrol_current,
				  * tankcontrol_previous,
				  * tankcontrol_gradient;

		struct ValveStruct 	* valvecontrol_current,
					* valvecontrol_previous,
					* valvecontrol_gradient;

		iteration_count = 1;
		tankcontrol_current = (struct TankStruct *)calloc(Ntanks,sizeof(struct TankStruct));
		tankcontrol_previous = (struct TankStruct *)calloc(Ntanks,sizeof(struct TankStruct));
		tankcontrol_gradient = (struct TankStruct *)calloc(Ntanks,sizeof(struct TankStruct));


		valvecontrol_current = (struct ValveStruct *) calloc(Nvalves,sizeof(struct ValveStruct));
		valvecontrol_previous = (struct ValveStruct *) calloc(Nvalves,sizeof(struct ValveStruct));
		valvecontrol_gradient = (struct ValveStruct *) calloc(Nvalves,sizeof(struct ValveStruct));

		memcpy(tankcontrol_current,tankcontrol,Ntanks*sizeof(struct TankStruct));
		memcpy(tankcontrol_previous,tankcontrol,Ntanks*sizeof(struct TankStruct));
		memcpy(tankcontrol_gradient,tankcontrol,Ntanks*sizeof(struct TankStruct));

		memcpy(valvecontrol_current,valvecontrol,Nvalves*sizeof(struct ValveStruct));
		memcpy(valvecontrol_previous,valvecontrol,Nvalves*sizeof(struct ValveStruct));
		memcpy(valvecontrol_gradient,valvecontrol,Nvalves*sizeof(struct ValveStruct));

		function_value_previous = objective_function(tankcontrol_current,valvecontrol_current);
		function_value_previous_printiter= function_value_previous;

		// Performs 1000 iterations
		while((iteration_count <= myOUTER_MAXITER)) {

			memcpy(tankcontrol_previous,tankcontrol_current,Ntanks*sizeof(struct TankStruct));
			memcpy(valvecontrol_previous,valvecontrol_current,Nvalves*sizeof(struct ValveStruct));

			compute_gradient(tankcontrol_current, valvecontrol_current,tankcontrol_gradient,valvecontrol_gradient);
			update_control(tankcontrol_current, valvecontrol_current,tankcontrol_gradient,valvecontrol_gradient, iteration_count);
			function_value_previous = function_value_current;
			clock_t clk;
			clk = clock();
			function_value_current = objective_function(tankcontrol_current,valvecontrol_current);
			clk = clock()-clk;
	//		printf("Time taken for func compu = %f",((double)clk)*1000000/CLOCKS_PER_SEC);
			iteration_count++;

			// Runs this section every 100 iteration
			if((((iteration_count%PRINT_INTERVAL_FUNCTIONVALUE)==0))) {
				float functionValueDiff = (function_value_previous_printiter - function_value_current);
				//			float functionValueDiff = function_value_current - function_value_previous;
				// @Ram The line below was changed by @Ram
				//	if( functionValueDiff >= function_delta) 
				if(functionValueDiff >= (function_delta*function_value_previous_printiter)) 
				{
					functionValueCounter = 0;
				}
				else 
				{
					functionValueCounter++;
					//printf("\n Function value previous = %f",function_value_previous_printiter);
					//printf("\n Fucntion value current = %f",function_value_current);
				}
				if(functionValueCounter >= 6)
				{
					exitprocedure(tankcontrol, valvecontrol);
				}
				//@Ram Changed the below print statement @Ram
				//	printf("\n Iteration Count = %d\tFunction Value = %f \t Percentage Difference = %f",iteration_count,function_value_current,(function_value_current-function_value_previous)/function_value_previous);
				printf("\n Iteration Count = %d\tFunction Value = %f \t Percentage Difference = %f",iteration_count,function_value_current,(function_value_previous_printiter-function_value_current)/function_value_previous_printiter);
				fflush(stdout);
				fprintf(fptr, "\n Iteration Count = %d\tFunction Value = %f \t Percentage Difference = %f",iteration_count,function_value_current,(function_value_previous_printiter-function_value_current)/function_value_previous_printiter);
				fflush(fptr);
				function_value_previous_printiter=function_value_current;
			}

			if(DEBUG) {
				// Print intermediate valve and tank values
				if(((iteration_count%PRINT_INTERVAL_VALVEVALUES) == 1)||(iteration_count == 1)) {
					//	printf("\nIteration Count = %d\n",iteration_count);
					fflush(stdout);
					fprintf(fptr, "\nIteration Count = %d\n",iteration_count);
					fflush(fptr);
					// Print final results (Tank and Valve Values)
					display(tankcontrol_current, tankcontrol_gradient, valvecontrol_current, valvecontrol_gradient);
				}
			}
		}
		fflush(fptr);

		// Copying the structures into orignal ones
		memcpy(tankcontrol,tankcontrol_current,Ntanks*sizeof(struct TankStruct));
		memcpy(valvecontrol,valvecontrol_current,Nvalves*sizeof(struct ValveStruct));

		// free all structures
		free(tankcontrol_current);
		free(tankcontrol_previous);
		free(tankcontrol_gradient);
		free(valvecontrol_current);
		free(valvecontrol_previous);
		free(valvecontrol_gradient);
	}

	double objective_function(struct TankStruct *tankcontrol_current,struct ValveStruct *valvecontrol_current)
	{
		/*
		 * Calculates the total penalty of a particular simulation
		 * tankcontrol : Tank structure compatible with EPANET
		 * valvecontrol : Valve structure compatible with EPANET 
		 */

		double func_value = 0.0;
		int temp_count;
		int temp_count2;
		int timeperiod;
		float temp_float_var = 0.0;

		timeperiod = tankcontrol_current[0].TimePeriod;
		update_tank_level(tankcontrol_current);

		//Tank Levels at each time and Penalty Computation (2M = A + B)
		temp_count = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
		for(; temp_count <= timeperiod; temp_count++) {
			for(temp_count2 = 0; temp_count2 < Ntanks; temp_count2++) {
				//if(strcmp(tankcontrol_current[temp_count2].TankID,"34")==0) continue;
				if(!targetMode)
				{
					temp_float_var = (2.0 * tankcontrol_current[temp_count2].TankLevels[temp_count]-(tankcontrol_current[temp_count2].MaxTankLevel+tankcontrol_current[temp_count2].MinTankLevel));
					// @Ram made these changes
					func_value += 10.0 * w12 * pow(temp_float_var,2)/tankcontrol_current[temp_count2].MaxTankLevel;
				}
				else
				{
					temp_float_var = (2.5 * tankcontrol_current[temp_count2].TankLevels[temp_count]-(tankcontrol_current[temp_count2].MaxTankLevel+tankcontrol_current[temp_count2].MinTankLevel));
					// @Ram made these changes
					func_value += w12 * pow(temp_float_var,2)/tankcontrol_current[temp_count2].MaxTankLevel;
				}
			}
		}
		fflush(fptr);

		//periodicity Computation and Penalty
		if (!targetMode) {
			for(temp_count2 = 0; temp_count2 < Ntanks; temp_count2++) {
				// A hack to ignore the main tank in optimisation
				//if(strcmp(tankcontrol_current[temp_count2].TankID,"34")==0) continue;

				//penalise only if the tank level at time T is less than the initial.

				temp_float_var = tankcontrol_current[temp_count2].TankLevels[timeperiod]-tankcontrol_current[temp_count2].TankLevels[0]; // penalise non periodicity.
				//			func_value += w3*duration * pow(temp_float_var,2)/tankcontrol_current[temp_count2].MaxTankLevel;
				func_value += w3*duration * fabs(temp_float_var/tankcontrol_current[temp_count2].MaxTankLevel);

			}
		}
		else {
			// Check the difference of the 24th hour tank level to the 0th hour tank level of the solution file
			for(temp_count2 = 0; temp_count2 < Ntanks; temp_count2 ++) {
				temp_float_var = tankcontrol_current[temp_count2].TankLevels[timeperiod] - solutionTankLevel[0].tankId[temp_count2];
				//			printf("\n temp_floatvar = %f",temp_float_var);
			//	printf("\n timeperiod = %d",timeperiod);
			//	printf("\n tank value = %f ",solutionTankLevel[0].tankId[temp_count2]*100/tankcontrol_current[temp_count2].MaxTankLevel);
				//	temp_float_var = tankcontrol_current[temp_count2].TankLevels[MAX_HOURS-1] - solutionTankLevel[0].tankId[temp_count2];
				//temp_float_var = tankcontrol_current[temp_count2].TankLevels[23] - solutionTankLevel[0].tankId[temp_count2];
				//@Ram made these changes
				func_value += w3*(10.0)*duration * fabs(temp_float_var/tankcontrol_current[temp_count2].MaxTankLevel);
				//	func_value += w3 * duration * pow(temp_float_var,2)/tankcontrol_current[temp_count2].MaxTankLevel;
			}
		}

		// Penalise the valve changes so only important changes are allowed
		for(temp_count = 1; temp_count < Nvalves; temp_count++) {
			temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
			// @Ram changed timeperiod to timeperiod - 1
			for(; temp_count2 < timeperiod ; temp_count2++) {
				temp_float_var = valvecontrol_current[temp_count].ValveValues[temp_count2 + 1] - valvecontrol_current[temp_count].ValveValues[temp_count2];
				// Multiplier arbitarily choosen to match the other penalty function values
				// @Ram made the following changes
				func_value += w4 * fabs(temp_float_var/valvecontrol_current[temp_count].MaxValue); 
			}
		}

		// If valve value is more than the max valve value, penalise objective value.
		for(temp_count = 1 ; temp_count <Nvalves; temp_count++) {
			temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );

			for(; temp_count2 < timeperiod; temp_count2++) {
				temp_float_var = (valvecontrol_current[temp_count].ValveValues[temp_count2]-valvecontrol_current[temp_count].MaxValue);
				if(temp_float_var > 0) {
					func_value+=temp_float_var*temp_float_var;
				}
				func_value += w5 * fabs(temp_float_var)/valvecontrol_current[temp_count].MaxValue;
			}
		}

		// Penalise the negative flows across the valves
		for(temp_count = 1; temp_count < Nvalves; temp_count++) {
			temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );

			for(; temp_count2 < timeperiod; temp_count2++) {
				ENgetlinkvalue(valvecontrol_current[temp_count].ValveLink,EN_FLOW, &(valvecontrol_current[temp_count].Flow[temp_count2]));			
				if (valvecontrol_current[temp_count].Flow[temp_count2] < 0){
					// temp_float_var = -actual_flow * weight * max_flow_of_valve
					temp_float_var = valvecontrol_current[temp_count].Flow[temp_count2] * valvecontrol_current[temp_count].MaxValue;
					//printf("\n Negative Flow across Valve %s = %f",valvecontrol_current[temp_count].ValveID, valvecontrol_current[temp_count].Flow[temp_count2]);
					func_value += (-1 * w6) * temp_float_var;
				}
			}
		}


		// Total Penalty function
		func_value = func_value/1000;
		fflush(fptr);
		return func_value;
	}

	void update_tank_level(struct TankStruct *tankcontrol)
	{
		/*
		 * Updates the tank levels in the tank structure
		 * tankcontrol : Tank structure compatible with EPANET
		 */

		int temp_count;
		int temp_count2;
		int timeperiod;

		// Tank levels are in meter^3. Tank Capacity is in meter^3. Inflows are in LPS.
		// Currently, usage data is assumed to be hourly. Need to convert everthing to meter^3 per hour.
		timeperiod = tankcontrol[0].TimePeriod;
		// @Ram Below is a fancy way of writing if else i.e., condition ? expression-true : expression-false
		temp_count = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
		for(; temp_count < timeperiod; temp_count++) {
			for(temp_count2 = 0; temp_count2 < Ntanks; temp_count2++) {
				tankcontrol[temp_count2].TankLevels[temp_count+1] = tankcontrol[temp_count2].TankLevels[temp_count]+(tankcontrol[temp_count2].InFlow[temp_count]-tankcontrol[temp_count2].OutFlow[temp_count]-tankcontrol[temp_count2].OutFlow_Epanet[temp_count])*LPS_TANKUNITS; //TankLevel is water volume in meter^3.
			}
		}
		fflush(fptr);
	}

	void compute_sim(struct TankStruct *tankstruct, struct ValveStruct *valvestruct)
	{
		/*
		 * Computes the flow across valves and tanks using EPANET
		 * tankcontrol : Tank structure compatible with EPANET
		 * valvecontrol : Valve structure compatible with EPANET
		 * simulation_time : Time period to simulate the flows
		 */

		int timeperiod;
		int temp_count;
		long t, tstep;

		time_t currentTime;
		currentTime = time(NULL);
		struct tm tm = *localtime(&currentTime);
		char tankFileName[128];
		sprintf(tankFileName, "../../Output_Files/Tank-Simulation-Output-%d-%d-%d::%d:%d:%d.csv", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);

		char valveFileName[128];
		sprintf(valveFileName, "../../Output_Files/Valve-Simulation-Output-%d-%d-%d::%d:%d:%d.csv", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);

		// @Ram timeperiod is set to 23*3600 seconds
		timeperiod = simDuration*3600;
		//printf("\n sim Duration = %d",simDuration);
		//printf("\n timeperiod = %d",timeperiod);
		ENsettimeparam(EN_DURATION, timeperiod);
		ENsettimeparam(EN_HYDSTEP, 600);
		int startTime = 0;
		
	//	startTime = 20;
	//	printf("\n start time = %d",startTime);

		FILE *tankPtr;
		FILE *valvePtr;

		tankPtr = fopen(tankFileName, "w");
		valvePtr = fopen(valveFileName, "w");

		if(tankPtr == NULL || valvePtr == NULL) {
			perror("Error in opening the file\n");
			exit(-1);
		}

		printf("\nSimulation :: Output Tank File Path -> %s\n", tankFileName);
		printf("\nSimulation :: Output Valve File Path -> %s\n", valveFileName);
		fflush(stdout);

		fprintf(fptr, "\nSimulation :: Output Tank File Path -> %s\n", tankFileName);
		fprintf(fptr, "\nSimulation :: Output Valve File Path -> %s\n", valveFileName);
		fprintf(fptr, "\n==============================================================================================\n");
		fflush(fptr); 

		fprintf(tankPtr, "Time,");
		fprintf(valvePtr, "Time,");
		fflush(tankPtr);
		fflush(valvePtr);
		/*for(temp_count=0;temp_count<24;temp_count++)
			{
				printf("\n Demand[%d]=%f",temp_count,tankstruct[0].OutFlow[temp_count]);
			}*/

		for(temp_count = 0; temp_count< Ntanks; temp_count++) {
			fprintf(tankPtr, "%s,", tankstruct[temp_count].TankID);
			fflush(tankPtr);
			if(temp_count == Ntanks-1) {
				fprintf(tankPtr, "\n");
				fflush(tankPtr);
			}
		}
		for(temp_count = 0; temp_count< Nvalves; temp_count++) {
			fprintf(valvePtr, "%s,", valvestruct[temp_count].ValveID);
			fflush(valvePtr);
			if(temp_count == Nvalves-1) {
				fprintf(valvePtr, "\n");
				fflush(valvePtr);
			}
		}
		fflush(fptr);
		fflush(tankPtr);
		fflush(valvePtr);

		fprintf(tankPtr, "%d,", solutionTankLevel[0].timeVal);
		fprintf(valvePtr, "%d,", valveSolutionLevel[0].timeVal);
		fflush(tankPtr);
		fflush(valvePtr);
		// Call epanet for computation
		ENopenH();
		float inflow = 0, outflow = 0, valveValue = 0;
		for(temp_count=0; temp_count<Nvalves; temp_count++) {
			ENsetlinkvalue(valvestruct[temp_count].ValveLink,EN_INITSETTING, valveSolutionLevel[startTime].valveId[temp_count]);
		//	ENsetlinkvalue(valvestruct[temp_count].ValveLink,EN_INITSETTING, valveSolutionLevel[startTime].valveId[temp_count]*valvestruct[temp_count].MaxValue/100);
			fprintf(valvePtr, "%f,", valveSolutionLevel[startTime].valveId[temp_count]);
		//	printf("line[%d]::startTime[%d]:: ENSetlinkvalue -> valveid[%s] :: valve level[%f] \n", __LINE__, startTime, valvestruct[temp_count].ValveID, valveSolutionLevel[startTime].valveId[temp_count]);
			fflush(valvePtr);
			if(temp_count == Nvalves-1) {
				fprintf(valvePtr, "\n");
				fflush(valvePtr);
			}
		}

		for(temp_count=0; temp_count<Ntanks; temp_count++) {
			ENsetnodevalue(tankstruct[temp_count].TankID,EN_TANKLEVEL, solutionTankLevel[startTime].tankId[temp_count]);
			fprintf(tankPtr, "%f,", (solutionTankLevel[startTime].tankId[temp_count]*100/tankstruct[temp_count].MaxTankLevel));
			fflush(tankPtr);
			if(temp_count == Ntanks-1) {
				fprintf(tankPtr, "\n");
				fflush(tankPtr);
			}
			//printf("line[%d]::startTime[%d]:: ENSetnodevalue -> TankID[%s] :: tank level[%f] \n", __LINE__, startTime, tankstruct[temp_count].TankID, solutionTankLevel[startTime].tankId[temp_count]);
		}

		ENinitH(0);
		int count=0;
		int linkLoop;
		float tempinflow, tempoutflow;
	//	printf("\n timestep %ld",tstep);
		
		double temp[Ntanks]; 
		for(temp_count=0;temp_count<Ntanks;temp_count++)
		{
		temp[temp_count] = solutionTankLevel[startTime].tankId[temp_count]; //TankLevel is water volume in meter^3.
		}
		
		//@Ram added a new variable
		float ram_tempflow=0.0;
		int ram_index=0;
		do {
			ENrunH(&t);
			ENnextH(&tstep);
			ENsettimeparam(EN_HYDSTEP, tstep);
		//	printf("\n %d",ram_index);
			//printf("t[%ld] :: tsetp[%ld] :: timeperiod[%ld]\n", t, tstep, timeperiod);
			//printf("\n t= %d",t);
			//if((t >0 && t < timeperiod) && !(t%3600))
			if((t < timeperiod) && !(t%3600))
			{
				//printf("\n t/3600 =  %f",(t/3600.0));
			//	printf("\n time  = %d,", solutionTankLevel[0].timeVal+startTime+1);
				fprintf(tankPtr, "%d,", solutionTankLevel[0].timeVal+startTime+1);
				//printf("\n %d,", solutionTankLevel[0].timeVal+startTime+1);
				fprintf(valvePtr, "%d,", valveSolutionLevel[0].timeVal+startTime+1);
				for(temp_count=0; temp_count<Ntanks; temp_count++) {
					inflow = 0;
					outflow = 0;
					for(linkLoop = 0; linkLoop < MAX_LINKSIZE; linkLoop++) {
						tempinflow = 0;
						tempoutflow = 0;
						//printf(" \n TankID[%d] Input Link[%d] OutputLink[%d]",temp_count,tankstruct[temp_count].InputLink[linkLoop],tankstruct[temp_count].OutputLink[linkLoop]);
						ENgetlinkvalue(tankstruct[temp_count].InputLink[linkLoop],EN_FLOW, &tempinflow);
						ENgetlinkvalue(tankstruct[temp_count].OutputLink[linkLoop],EN_FLOW, &tempoutflow);
						inflow += tempinflow;
						outflow += tempoutflow;
					}
					//	ENgetnodeindex(ram_id,&ram_index);
				//		ENgetlinkindex("DVRTMP17",&ram_index);
				//		ENgetlinkvalue(ram_index,EN_FLOW, &ram_tempflow);
				//		printf("\n Flow at node 19 = %f",ram_tempflow);
				//		printf("\n %d",ram_index);
						/*ENgetlinkindex("DVRTMP15",&ram_index);
						ENgetlinkvalue(ram_index,EN_FLOW, &ram_tempflow);
						printf("\n Flow at valve DVRTMP15 = %f",ram_tempflow);
						ENgetlinkindex("DVRTMP20",&ram_index);
						ENgetlinkvalue(ram_index,EN_FLOW, &ram_tempflow);
						printf("\n Flow at valve DVRTMP20 = %f",ram_tempflow);
						ENgetlinkindex("78",&ram_index);
						ENgetlinkvalue(ram_index,EN_SETTING, &ram_tempflow);
						printf("\n Setting at valve 78 = %f",ram_tempflow);
						ENgetlinkindex("76",&ram_index);
						ENgetlinkvalue(ram_index,EN_SETTING, &ram_tempflow);
						printf("\n Setting at valve 76 = %f",ram_tempflow);
						ENgetlinkindex("74",&ram_index);
						ENgetlinkvalue(ram_index,EN_SETTING, &ram_tempflow);
						printf("\n Setting at valve 74 = %f",ram_tempflow);
						ENgetlinkindex("40",&ram_index);
						ENgetlinkvalue(ram_index,EN_SETTING, &ram_tempflow);
						printf("\n Setting at valve 40 = %f",ram_tempflow);
						ENgetlinkindex("3",&ram_index);
						ENgetlinkvalue(ram_index,EN_SETTING, &ram_tempflow);
						printf("\n Setting at valve 3 = %f",ram_tempflow);
						*/
				//	}
				//	double temp = solutionTankLevel[startTime].tankId[temp_count] + (inflow-outflow-tankstruct[temp_count].OutFlow_Epanet[solutionTankLevel[0].timeVal+startTime])*LPS_TANKUNITS; //TankLevel is water volume in meter^3.
				if(strcmp(tankstruct[temp_count].TankID,"Kesare")==0)
					{
					printf( "\n Tank value = %f \n", temp[temp_count]);
					}

					temp[temp_count] = temp[temp_count] + (inflow-outflow)*LPS_TANKUNITS;
					//temp[temp_count] = temp[temp_count] + (inflow-outflow);
					//printf("temp[%f] ", temp);
					//temp = temp*100/tankstruct[temp_count].MaxTankLevel;
					//printf("[%f]\n", temp);
					fprintf(tankPtr, "%f,", (temp[temp_count]*100.0/tankstruct[temp_count].MaxTankLevel));
					fflush(tankPtr);
				if(strcmp(tankstruct[temp_count].TankID,"Rajiv_nagar_II")==0)
					{
					printf("\n Tank No[%d] :: Inflow[%f] :: outflow[%f] :: netflow[%f]", temp_count, inflow, outflow, (inflow-outflow));
					printf("\n Tank No[%d] :: Inflow[%f] :: outflow[%f] :: netflow[%f]", temp_count, inflow*LPS_TANKUNITS, outflow*LPS_TANKUNITS, (inflow-outflow)*LPS_TANKUNITS);
					printf("\n Tank No[%s]:: Time[%d] :: outflow[%f] ",tankstruct[temp_count].TankID,solutionTankLevel[0].timeVal+startTime+1,tankstruct[temp_count].OutFlow[solutionTankLevel[0].timeVal+startTime]*LPS_TANKUNITS);

						ENgetlinkindex("35",&ram_index);
						ENgetlinkvalue(ram_index,EN_SETTING, &ram_tempflow);
						printf("\n Setting at valve 35 = %f ",ram_tempflow);
					printf( "\n Tank value = %f ", temp[temp_count]);
					printf( "\n Tank Normalised value = %f \n", (temp[temp_count]*100.0/tankstruct[temp_count].MaxTankLevel));
					fflush(stdout);
					}
					if(temp_count == Ntanks-1) {
						fprintf(tankPtr, "\n");
						fflush(tankPtr);
					}
				}
				for(temp_count=0; temp_count<Nvalves; temp_count++) {	
					//printf("valvelink[%d] ", valvestruct[temp_count].ValveLink);
					ENgetlinkvalue(valvestruct[temp_count].ValveLink,EN_FLOW, &(valveValue));
					fprintf(valvePtr, "%f,", valveValue);
					fflush(valvePtr);
					if(temp_count == Nvalves-1) {
						fprintf(valvePtr, "\n");
						fflush(valvePtr);
					}
					//printf("valveValue[%f] \n", valveValue);
				}

				//printf("\n");
				startTime++;
				for(temp_count=0; temp_count<Nvalves; temp_count++) {
					ENsetlinkvalue(valvestruct[temp_count].ValveLink,EN_SETTING, valveSolutionLevel[startTime].valveId[temp_count]);
					//ENsetlinkvalue(valvestruct[temp_count].ValveLink,EN_INITSETTING, valveSolutionLevel[startTime].valveId[temp_count]*valvestruct[temp_count].MaxValue/100);
					printf("line[%d]::startTime[%d]:: ENSetlinkvalue -> valveId[%s] :: valve level[%f] \n", __LINE__, startTime, valvestruct[temp_count].ValveID, valveSolutionLevel[startTime].valveId[temp_count]);
				}
			}
			count++;
		} while (tstep > 0);
		//printf("count[%d]\n", count);

		ENcloseH();
		fflush(valvePtr);
		fclose(valvePtr);

		fflush(tankPtr);
		fclose(tankPtr);

		// SIM Tank Output File copy
		char finalOutput[1024];
		sprintf(finalOutput, "cp %s ../../result/sim_tank.csv ", tankFileName);
		int res = system(finalOutput);	

		if(res != 0)
			perror("system command fail");

		// SIM Valve Output File copy
		char valveOutput[1024];
		sprintf(valveOutput, "cp %s ../../result/sim_valve.csv ", valveFileName);
		res = system(valveOutput);	

		if(res != 0)
			perror("system command fail");

	}

	void compute_flows(struct TankStruct *tankstruct, struct ValveStruct *valvestruct, int simulation_time)
	{
		/*
		 * Computes the flow across valves and tanks using EPANET
		 * tankcontrol : Tank structure compatible with EPANET
		 * valvecontrol : Valve structure compatible with EPANET
		 * simulation_time : Time period to simulate the flows
		 */

		int timeperiod;
		int temp_count;
		long t, tstep;

		timeperiod = tankstruct[0].TimePeriod;
		simulation_time = simulation_time%timeperiod;
		//printf("\n sim_time2 = %d\n", simulation_time);

		// Call epanet for computation
		ENopenH();
		for(temp_count=0; temp_count<Nvalves; temp_count++) {
			ENsetlinkvalue(valvestruct[temp_count].ValveLink,EN_INITSETTING,valvestruct[temp_count].ValveValues[simulation_time]);
		}
		ENinitH(10);
		do {
			ENrunH(&t);
			ENnextH(&tstep);
		} while (tstep > 0);


		int linkLoop; 
		float tempInFlow, tempOutFlow;
		for(temp_count=0; temp_count<Ntanks; temp_count++) {
			tankstruct[temp_count].InFlow[simulation_time] = 0;
			tankstruct[temp_count].OutFlow_Epanet[simulation_time] = 0;
			for(linkLoop = 0; linkLoop < MAX_LINKSIZE; linkLoop++) {
				// Initialising temp variables so that Input and output flow can be measured for every inlet and outlet
				tempInFlow = 0; 
				tempOutFlow = 0;
				// Check if there is any index in the InputLink array; Mostly all tanks will have only one input
				if(tankstruct[temp_count].InputLink[linkLoop] > 0) {
					ENgetlinkvalue(tankstruct[temp_count].InputLink[linkLoop],EN_FLOW, &tempInFlow);
					tankstruct[temp_count].InFlow[simulation_time] += tempInFlow;
				}
				// Check if there is any index in the OutputLink array; Mostly all tanks will have only one output
				if(tankstruct[temp_count].OutputLink[linkLoop] > 0) {
					ENgetlinkvalue(tankstruct[temp_count].OutputLink[linkLoop],EN_FLOW, &tempOutFlow);
					tankstruct[temp_count].OutFlow_Epanet[simulation_time] += tempOutFlow;
				}
			}
		}

		ENcloseH();
	}

	void compute_gradient(struct TankStruct *tankcontrol_current, struct ValveStruct *valvecontrol_current, struct TankStruct *tankcontrol_gradient,struct ValveStruct *valvecontrol_gradient)
	{
		/*
		 * Computes the gradient which decides the amount of valve opening for a particular scenario.
		 * tankcontrol_current : Current Tank structure compatible with EPANET.
		 * valvecontrol_current : Current Valve structure compatible with EPANET.
		 * tankcontrol_gradient : Same as Tank Structure but contains the gradient values.
		 * valvecontrol_gradient : Same as Valve Structure but contains the gradient values.
		 */

		int temp_count;
		int temp_count2;
		int timeperiod;

		float function_value_1,
		      function_value_2,
		      temp_rand=0.0;

		int **rand_valve_vector; // @Ram a temporary array to store the random vector for gradient computation
		int *rand_tank_vector; // @Ram a temporary array to store the random vector for gradient computation
		timeperiod = tankcontrol_current->TimePeriod;

		rand_valve_vector = (int **)calloc(Nvalves, sizeof(int *));	
		for(temp_count=0;temp_count<Nvalves;temp_count++)
		{
			rand_valve_vector[temp_count]= (int*) calloc(timeperiod,sizeof(int));
		}
		rand_tank_vector = (int*) calloc(Ntanks,sizeof(int));

		memcpy(temp_tankstruct,tankcontrol_current,Ntanks*sizeof(struct TankStruct));
		memcpy(temp_valvestruct,valvecontrol_current,Nvalves*sizeof(struct ValveStruct));

		//compute gradient for valves
		for(temp_count =0; temp_count <Nvalves; temp_count++) {
			temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
			for(; temp_count2 < timeperiod; temp_count2++) {
				temp_rand = (double) rand() / (RAND_MAX); 
				if(temp_rand <=0.5)
				{
					rand_valve_vector[temp_count][temp_count2]=1; 
					temp_valvestruct[temp_count].ValveValues[temp_count2] = valvecontrol_current[temp_count].ValveValues[temp_count2]+rand_valve_vector[temp_count][temp_count2]*epsilon;
					if(temp_valvestruct[temp_count].ValveValues[temp_count2] < 0) { 
						// valve value cannot be negative
						temp_valvestruct[temp_count].ValveValues[temp_count2] = 0;
					}
				}
				else if(temp_rand<=1)
				{
					rand_valve_vector[temp_count][temp_count2]=-1; 
					temp_valvestruct[temp_count].ValveValues[temp_count2] = (valvecontrol_current[temp_count].ValveValues[temp_count2]+rand_valve_vector[temp_count][temp_count2]*epsilon);
					if(temp_valvestruct[temp_count].ValveValues[temp_count2] < 0) { 
						// valve value cannot be negative
						temp_valvestruct[temp_count].ValveValues[temp_count2] = 0;
					}
				}
				else
				{
					printf("\n Error in Temp_Rand");
				}

			}
		}

		if(!targetMode)
		{
			for(temp_count=0;temp_count<Ntanks;temp_count++)
			{
				temp_rand = (double) rand() / (RAND_MAX); 
				if(temp_rand<=0.5)
				{
					rand_tank_vector[temp_count]=1;
					temp_tankstruct[temp_count].TankLevels[0] = tankcontrol_current[temp_count].TankLevels[0]+ rand_tank_vector[temp_count]*epsilon;
					if(temp_tankstruct[temp_count].TankLevels[0] < 0) { 
						// Tank value cannot be negative
						temp_tankstruct[temp_count].TankLevels[0] = 0;
					}

				}	
				else if(temp_rand<1)
				{
					rand_tank_vector[temp_count] = -1;
					temp_tankstruct[temp_count].TankLevels[0] = tankcontrol_current[temp_count].TankLevels[0]+ rand_tank_vector[temp_count]*epsilon;

					if(temp_tankstruct[temp_count].TankLevels[0] < 0) { 
						// Tank value cannot be negative
						temp_tankstruct[temp_count].TankLevels[0] = 0;
					}
				}

			}
		}
		temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
		for(; temp_count2 < timeperiod; temp_count2++) {
			compute_flows(temp_tankstruct,temp_valvestruct,temp_count2);
		}
		function_value_1 = objective_function(temp_tankstruct,temp_valvestruct);

		for(temp_count =0; temp_count <Nvalves; temp_count++) {
			temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
			for(; temp_count2 < timeperiod; temp_count2++) {
				temp_valvestruct[temp_count].ValveValues[temp_count2] =  valvecontrol_current[temp_count].ValveValues[temp_count2]-rand_valve_vector[temp_count][temp_count2]*epsilon;
				if(temp_valvestruct[temp_count].ValveValues[temp_count2] < 0) { 
					// valve value cannot be negative
					temp_valvestruct[temp_count].ValveValues[temp_count2] = 0;
				}

			}
		}
		if(!targetMode)
		{
			for(temp_count = 0; temp_count<Ntanks; temp_count++)
			{
				temp_tankstruct[temp_count].TankLevels[0] = tankcontrol_current[temp_count].TankLevels[0] - rand_tank_vector[temp_count]* epsilon;
				if(temp_tankstruct[temp_count].TankLevels[0] < 0) { 
					// Tank value cannot be negative
					temp_tankstruct[temp_count].TankLevels[0] = 0;
				}
			}
		}

		temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
		for(; temp_count2 < timeperiod; temp_count2++) {
			compute_flows(temp_tankstruct,temp_valvestruct,temp_count2);
		}
		function_value_2 = objective_function(temp_tankstruct,temp_valvestruct);	

		/*	compute_flows(temp_tankstruct,temp_valvestruct,temp_count2);
			function_value_1 = objective_function(temp_tankstruct,temp_valvestruct);

			temp_valvestruct[temp_count].ValveValues[temp_count2] = (valvecontrol_current[temp_count].ValveValues[temp_count2]-epsilon);
			if(temp_valvestruct[temp_count].ValveValues[temp_count2] < 0) { 
		// valve value cannot be negative
		temp_valvestruct[temp_count].ValveValues[temp_count2] = 0;
		}
		compute_flows(temp_tankstruct,temp_valvestruct,temp_count2);
		function_value_2 = objective_function(temp_tankstruct,temp_valvestruct);
		 */
		for(temp_count =0; temp_count <Nvalves; temp_count++) {
			temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
			for(; temp_count2 < timeperiod; temp_count2++) {
				valvecontrol_gradient[temp_count].ValveValues[temp_count2] = (function_value_1-function_value_2)/(2*rand_valve_vector[temp_count][temp_count2]*epsilon); 
			}
		}
		if(!targetMode)
		{
			for(temp_count = 0; temp_count < Ntanks; temp_count++)
			{
				tankcontrol_gradient[temp_count].TankLevels[0] = (function_value_1-function_value_2)/(2*rand_tank_vector[temp_count]*epsilon);
			}
		}
		//gradient computation

		//temp_valvestruct[temp_count].ValveValues[temp_count2] = valvecontrol_current[temp_count].ValveValues[temp_count2]; 
		//make the value back to original.
		//}
		//}

		//reinitialise temp structures. This step need not be required as nothing is changed in the above computation.
		//memcpy(temp_tankstruct,tankcontrol_current,Ntanks*sizeof(struct TankStruct));
		//memcpy(temp_valvestruct,valvecontrol_current,Nvalves*sizeof(struct ValveStruct));
		//compute gradient for tanks
/*for(temp_count =0; temp_count <Ntanks; temp_count++) {
  temp_tankstruct[temp_count].TankLevels[0] = tankcontrol_current[temp_count].TankLevels[0]+epsilon;

  function_value_1 = objective_function(temp_tankstruct,temp_valvestruct);

  temp_tankstruct[temp_count].TankLevels[0] = tankcontrol_current[temp_count].TankLevels[0]-epsilon;

  function_value_2 = objective_function(temp_tankstruct,temp_valvestruct);

  tankcontrol_gradient[temp_count].TankLevels[0] = (function_value_1-function_value_2)/(2*epsilon);
  temp_tankstruct[temp_count].TankLevels[0] = tankcontrol_current[temp_count].TankLevels[0]; //not necessary in gradient compuation as function values are separable in the tanklevels. Nevertheless to allay any fears.
  }*/
fflush(fptr);
//printf("\n Here is the error");
for(temp_count=0;temp_count<Nvalves;temp_count++)
{
	free(rand_valve_vector[temp_count]);
}
//printf("\n Not here");
free(rand_valve_vector);
free(rand_tank_vector);
}

void update_control(struct TankStruct *tankcontrol_current, struct ValveStruct *valvecontrol_current,struct TankStruct *tankcontrol_gradient,struct ValveStruct *valvecontrol_gradient, int main_iteration_count)
{
	/*
	 * Updates the tank and valve flow values according to the EPANET simulation.
	 * tankcontrol_current : Current Tank structure compatible with EPANET.
	 * valvecontrol_current : Current Valve structure compatible with EPANET.
	 * tankcontrol_gradient : Same as Tank Structure but contains the gradient values.
	 * valvecontrol_gradient : Same as Valve Structure but contains the gradient values.
	 * main_iteration_count : Main iteration count is taken into consideration for adding noise.
	 */

	int temp_count;
	int temp_count2;
	int timeperiod,
	    iteration_count = 0;
	double alpha_mult = 1,
	       rho = 0.6,
	       c = 0.001,
	       temp_valve_value = 0,
	       norm_grad_var=0.0,
	       beta = 0.6;
	float	random_number = 0.0;


	double  function_value_current = 0.0,
		function_value_next = 0.0;

	timeperiod = tankcontrol_current->TimePeriod;

	/*if(main_iteration_count==1) {
		printf("\n Simulation Parameters\n rho = %f. c = %f, beta = %f",rho,c,beta);
		fflush(stdout);
		fprintf(fptr, "\n Simulation Parameters\n rho = %f. c = %f, beta = %f",rho,c,beta);
		fflush(fptr);
	}*/
	memcpy(temp_tankstruct,tankcontrol_current,Ntanks*sizeof(struct TankStruct));
	memcpy(temp_valvestruct,valvecontrol_current,Nvalves*sizeof(struct ValveStruct));

	function_value_current = objective_function(tankcontrol_current,valvecontrol_current);
	function_value_next = function_value_current;
	/*alpha_mult = 0.2;
	  if((function_value_current < 500000)&&(function_value_current >=5000))
	  {
	  alpha_mult = 0.15;
	  }
	// TankValues are not updated.
	if((function_value_current < 5000)&&(function_value_current >=50))
	{
	alpha_mult = 0.1;
	function_delta = 0.05;
	}
	else if(function_value_current<50)
	{
	alpha_mult = 0.06;
	function_delta=0.001;
	}*/
	// @Ram norm_gradient is just the sum of individual gradient values
	norm_grad_var = norm_gradient(tankcontrol_gradient,valvecontrol_gradient);

	while(((function_value_next-function_value_current)>= -c *alpha_mult* norm_grad_var) && (iteration_count < myMAXITER)) {
		//printf("\n Line Search Iter = %f",alpha_mult);
		temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
		for(; temp_count2<timeperiod; temp_count2++) {
			for(temp_count =0; temp_count <Nvalves; temp_count++) {
				//random_direction(&random_number,1);
				temp_valve_value = valvecontrol_current[temp_count].ValveValues[temp_count2]-(valvecontrol_gradient[temp_count].ValveValues[temp_count2])*alpha_mult;
				if(temp_valve_value > 0) {
					temp_valvestruct[temp_count].ValveValues[temp_count2] = temp_valve_value;
				} 
				else {
					temp_valvestruct[temp_count].ValveValues[temp_count2] = 0;
				}
			}
			compute_flows(temp_tankstruct,temp_valvestruct,temp_count2);
		}

		for(temp_count =0; temp_count <Ntanks; temp_count++) {
			// random_direction(&random_number,1);
			temp_tankstruct[temp_count].TankLevels[0] = tankcontrol_current[temp_count].TankLevels[0]-(tankcontrol_gradient[temp_count].TankLevels[0])*alpha_mult;
		}
		//printf("\n alpha mult %f",alpha_mult);
		function_value_next = objective_function(temp_tankstruct,temp_valvestruct);

		alpha_mult*=rho;
		iteration_count++;
	}

	fflush(fptr);


	// Add noise
	temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
	for(; temp_count2<timeperiod; temp_count2++) {
		for(temp_count =0; temp_count <Nvalves; temp_count++) {
			random_direction(&random_number,1);
			temp_valve_value = valvecontrol_current[temp_count].ValveValues[temp_count2]-valvecontrol_gradient[temp_count].ValveValues[temp_count2]*(1+random_number*pow((main_iteration_count/100)+1,-0.6))*alpha_mult;

			if(temp_valve_value > 0) {
				temp_valvestruct[temp_count].ValveValues[temp_count2] = temp_valve_value;
			} 
			else {
				temp_valvestruct[temp_count].ValveValues[temp_count2] = 0;
			}
		}
		compute_flows(temp_tankstruct,temp_valvestruct,temp_count2);
	}

	for(temp_count =0; temp_count <Ntanks; temp_count++) {
		random_direction(&random_number,1);
		temp_tankstruct[temp_count].TankLevels[0] = tankcontrol_current[temp_count].TankLevels[0]-(tankcontrol_gradient[temp_count].TankLevels[0])*(1+random_number*pow((main_iteration_count/100)+1,-beta))*alpha_mult;
	}


	if(((main_iteration_count%PRINT_INTERVAL_FUNCTIONVALUE)==0)) {
		fprintf(fptr, " Gradient Norm = %f alpha mult = %f",norm_grad_var,alpha_mult);
	}
	fflush(fptr);
	memcpy(tankcontrol_current, temp_tankstruct,Ntanks*sizeof(struct TankStruct));
	memcpy(valvecontrol_current,temp_valvestruct,Nvalves*sizeof(struct ValveStruct));
}

float norm_gradient(struct TankStruct *tankcontrol_gradient,struct ValveStruct *valvecontrol_gradient)
{
	/*
	 * Normalize the Tank and Valve gradient structures.
	 * tankcontrol_gradient : Same as Tank Structure but contains the gradient values.
	 * valvecontrol_gradient : Same as Valve Structure but contains the gradient values.
	 */

	int temp_count =0,
	    temp_count2 = 0,
	    timeperiod;
	float temp_gradient_norm = 0.0;

	timeperiod = tankcontrol_gradient->TimePeriod;
	for(temp_count=0; temp_count<Nvalves; temp_count++) {
		temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
		for(; temp_count2<timeperiod; temp_count2++) {
			temp_gradient_norm += valvecontrol_gradient[temp_count].ValveValues[temp_count2]*valvecontrol_gradient[temp_count].ValveValues[temp_count2];
		}
	}
	for(temp_count2=0; temp_count2<Ntanks; temp_count2++) {
		temp_gradient_norm += tankcontrol_gradient[temp_count2].TankLevels[0]*tankcontrol_gradient[temp_count2].TankLevels[0];
	}
	fflush(fptr);
	return temp_gradient_norm;
}

float max_tank_level(int tank_struct_index)
{
	/*
	 * Calculates the max tank level of the tank
	 * tank_struct_index : Index of the Tank Structure.
	 */

	float area = 0.0,
	      tank_max_Vol,
	      tank_min_Vol,
	      tank_level_LH,
	      feet3_meter3_coversion = pow(12*2.54/100,3);


	tank_max_Vol = Tank[tank_struct_index].Vmax;
	tank_min_Vol = Tank[tank_struct_index].Vmin;
	tank_level_LH = tank_max_Vol-tank_min_Vol;	// This is in feet^3. COnvert to meter^3
	tank_level_LH *= feet3_meter3_coversion;

	return tank_level_LH;
}

float tank_level_meter3_to_Level(float tank_volume, int tank_struct_index)
{
	/*
	 * Converts the water level from feet^3 to meter^3
	 * tank_volume : Tank Volume
	 * tank_struct_index : Index of the Tank Structure.
	 */

	float tank_area = 0.0,
	      tank_level_LH,
	      feet2_meter2_coversion = pow(12*2.54/100,2);

	tank_area = Tank[tank_struct_index].A;
	tank_level_LH = tank_volume/(tank_area*feet2_meter2_coversion);	// This is in feet^3. Convert to meter^3
	return tank_level_LH;
}

void random_direction(float *random_vector, int dimension)
{
	/*
	 * Defines a random vector in a particular direction.
	 * random_vector : Random vector for adding noise.
	 * dimension : Dimension for the random vector.
	 */

	const int rand_max = 500;
	while(dimension > 0) {
		random_vector[dimension-1] = rand()%(2*rand_max);
		random_vector[dimension-1] = (random_vector[dimension-1]-rand_max)/rand_max;
		dimension--;
	}
	return ;
}

void display(struct TankStruct *tankcontrol_current, struct TankStruct *tankcontrol_gradient, struct ValveStruct *valvecontrol_current, struct ValveStruct *valvecontrol_gradient)
{
	/*
	 * Displays the Tank and Valve structures.
	 * tankcontrol_current : Current Tank structure compatible with EPANET.
	 * valvecontrol_current : Current Valve structure compatible with EPANET.
	 * tankcontrol_gradient : Same as Tank Structure but contains the gradient values.
	 * valvecontrol_gradient : Same as Valve Structure but contains the gradient values.
	 */

	int temp_count;
	int temp_count2;
	int timeperiod = tankcontrol_current[0].TimePeriod;

	// Display all the Valve Settings and Tank Levels till elapsed time
	// Elapsed time is the time till the valve or demand settings are made
	fprintf(fptr, "\n==============================================================================================\n");
	fflush(fptr);

	//Printing final Valve values
	for(temp_count = 0; temp_count< Nvalves; temp_count++){
		fprintf(fptr, "\n\n\n\n Valve ID = %s",valvecontrol_current[temp_count].ValveID);
		fflush(fptr);
		temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
		for(; temp_count2 < timeperiod; temp_count2++){
			fprintf(fptr, "\n Valve Value [%d] = %f \t Gradient = %f",
					temp_count2,
					valvecontrol_current[temp_count].ValveValues[temp_count2],
					valvecontrol_gradient[temp_count].ValveValues[temp_count2]);

			fflush(fptr);
		}
	}
	fflush(fptr);

	//Printing final tank levels
	for(temp_count = 0; temp_count< Ntanks; temp_count++){
		fprintf(fptr, "\n\n\n\n Tank ID = %s",tankcontrol_current[temp_count].TankID);
		fflush(fptr);
		temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
		for(; temp_count2 < timeperiod; temp_count2++){
			fprintf(fptr, "\n Tank Level [%d] = %f, \t Tank Inflow = %f, \t Tank Outflow = %f",
					temp_count2,
					(tankcontrol_current[temp_count].TankLevels[temp_count2])*100/tankcontrol_current[temp_count].MaxTankLevel,
					tankcontrol_current[temp_count].InFlow[temp_count2],
					(tankcontrol_current[temp_count].OutFlow[temp_count2] + tankcontrol_current[temp_count].OutFlow_Epanet[temp_count2]));
			fflush(fptr);
		}
	}
	fprintf(fptr, "\n==============================================================================================\n");
	fflush(fptr);
}

void Display_Output(struct TankStruct *tankcontrol, struct ValveStruct *valvecontrol)
{
	/*
	 * Displays the final output of Tank and Valve structures.
	 * tankcontrol_current : Current Tank structure compatible with EPANET.
	 * valvecontrol_current : Current Valve structure compatible with EPANET.
	 */

	int temp_count;
	int temp_count2;
	int timeperiod = tankcontrol[0].TimePeriod;

	time_t currentTime;
	currentTime = time(NULL);
	struct tm tm = *localtime(&currentTime);
	char tankFileName[128];
	sprintf(tankFileName, "../../Output_Files/Tank-Output-%d-%d-%d::%d:%d:%d.csv", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);

	char valveFileName[128];
	sprintf(valveFileName, "../../Output_Files/Valve-Output-%d-%d-%d::%d:%d:%d.csv", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);

	FILE *tankPtr;
	FILE *valvePtr;

	tankPtr = fopen(tankFileName, "w");
	valvePtr = fopen(valveFileName, "w");

	if(tankPtr == NULL || valvePtr == NULL) {
		perror("Error in opening the file\n");
		exit(-1);
	}

	printf("\nOutput Tank File Path -> %s\n", tankFileName);
	fflush(stdout);
	printf("\nOutput Valve File Path -> %s\n", valveFileName);
	fflush(stdout);
	fprintf(fptr, "\nOutput Tank File Path -> %s\n", tankFileName);
	fprintf(fptr, "\nOutput Valve File Path -> %s\n", valveFileName); 

	// Display all the Valve Settings and Tank Levels till elapsed time
	// Elapsed time is the time till the valve or demand settings are made
	fprintf(fptr, "\n==============================================================================================\n");
	fflush(fptr);

	// Creating Tank and valve setting final output file
	fprintf(tankPtr, "Time,");
	fprintf(valvePtr, "Time,");
	fflush(tankPtr);
	fflush(valvePtr);

	for(temp_count = 0; temp_count< Ntanks; temp_count++) {
		fprintf(tankPtr, "%s,", tankcontrol[temp_count].TankID);
		fflush(tankPtr);
		if(temp_count == Ntanks-1) {
			fprintf(tankPtr, "\n");
			fflush(tankPtr);
		}
	}
	for(temp_count = 0; temp_count< Nvalves; temp_count++) {
		fprintf(valvePtr, "%s,", valvecontrol[temp_count].ValveID);
		fflush(valvePtr);
		if(temp_count == Nvalves-1) {
			fprintf(valvePtr, "\n");
			fflush(valvePtr);
		}
	}
	fflush(fptr);
	fflush(tankPtr);
	fflush(valvePtr);

	temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
	// @Ram changed timeperiod to timeperiod +1
	if(!targetMode)
	{
	for(; temp_count2 < timeperiod+1; temp_count2 += hourly_scheduler) {
		fprintf(tankPtr, "%d,", temp_count2);
		fflush(tankPtr);

		for(temp_count = 0; temp_count< Ntanks; temp_count++) {
			fprintf(tankPtr, "%f,", (tankcontrol[temp_count].TankLevels[temp_count2])*100/tankcontrol[temp_count].MaxTankLevel);
			fflush(tankPtr);
			if(temp_count == Ntanks-1) {
				fprintf(tankPtr, "\n");
				fflush(tankPtr);
			}
		}
	}
	}
	else
	{
	for(; temp_count2 < timeperiod+2; temp_count2 += hourly_scheduler) {
		fprintf(tankPtr, "%d,", temp_count2);
		fflush(tankPtr);

		for(temp_count = 0; temp_count< Ntanks; temp_count++) {
			fprintf(tankPtr, "%f,", (tankcontrol[temp_count].TankLevels[temp_count2])*100/tankcontrol[temp_count].MaxTankLevel);
			fflush(tankPtr);
			if(temp_count == Ntanks-1) {
				fprintf(tankPtr, "\n");
				fflush(tankPtr);
			}
		}
	}
	}
	
	// @Ram added a separate loop for printing valve values
	temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
	for(; temp_count2 < timeperiod; temp_count2 += hourly_scheduler) {

		fprintf(valvePtr, "%d,", temp_count2);
		fflush(valvePtr);
		for(temp_count = 0; temp_count< Nvalves; temp_count++) {
			//fprintf(valvePtr, "%f,", valvecontrol[temp_count].ValveValues[temp_count2]*100/valvecontrol[temp_count].MaxValue);
			fprintf(valvePtr, "%f,", valvecontrol[temp_count].ValveValues[temp_count2]);
			fflush(valvePtr);
			if(temp_count == Nvalves-1) {
				fprintf(valvePtr, "\n");
				fflush(valvePtr);
			}
		}
	}

	fflush(fptr);
	fflush(tankPtr);
	fclose(tankPtr);

	fflush(valvePtr);
	fclose(valvePtr);

	// Tank Output File copy
	char finalOutput[1024];
	sprintf(finalOutput, "cp %s ../../result/tank.csv ", tankFileName);
	int res = system(finalOutput);	

	if(res != 0)
		perror("system command fail");

	// Valve Output File copy
	char valveOutput[1024];
	sprintf(valveOutput, "cp %s ../../result/valve.csv ", valveFileName);
	res = system(valveOutput);	

	if(res != 0)
		perror("system command fail");

	//Printing final Valve values
	for(temp_count = 0; temp_count< Nvalves; temp_count++) {
		fprintf(fptr, "\n\n\n\n Valve ID = %s",valvecontrol[temp_count].ValveID);
		fflush(fptr);
		temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
		for(; temp_count2 < timeperiod; temp_count2 += hourly_scheduler) {
			fprintf(fptr, "\n Valve Value [%d] = %f, \t Percentage Open = %f  ",
					temp_count2,
					valvecontrol[temp_count].ValveValues[temp_count2],
					valvecontrol[temp_count].ValveValues[temp_count2]*100/valvecontrol[temp_count].MaxValue);
			fflush(fptr);
		}
	}
	fflush(fptr);

	//Printing final tank levels
	for(temp_count = 0; temp_count< Ntanks; temp_count++){
		fprintf(fptr, "\n\n\n\n Tank ID = %s",tankcontrol[temp_count].TankID);
		fflush(fptr);
		temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
		if(!targetMode)
		{
		// @Ram changed timeperiod to timeperiod +1
		for(; temp_count2 < timeperiod+1; temp_count2 += hourly_scheduler) {
			fprintf(fptr, "\n Tank Level [%d] = %f, \t Tank Inflow = %f, \t Tank Outflow = %f",
					temp_count2,
					(tankcontrol[temp_count].TankLevels[temp_count2])*100/tankcontrol[temp_count].MaxTankLevel,
					tankcontrol[temp_count].InFlow[temp_count2],
					(tankcontrol[temp_count].OutFlow[temp_count2] + tankcontrol[temp_count].OutFlow_Epanet[temp_count2]));
			fflush(fptr);
		}
		}
		else
		{
		// @Ram changed timeperiod to timeperiod + 2
		for(; temp_count2 < timeperiod+2; temp_count2 += hourly_scheduler) {
			fprintf(fptr, "\n Tank Level [%d] = %f, \t Tank Inflow = %f, \t Tank Outflow = %f",
					temp_count2,
					(tankcontrol[temp_count].TankLevels[temp_count2])*100/tankcontrol[temp_count].MaxTankLevel,
					tankcontrol[temp_count].InFlow[temp_count2],
					(tankcontrol[temp_count].OutFlow[temp_count2] + tankcontrol[temp_count].OutFlow_Epanet[temp_count2]));
			fflush(fptr);
		}
		}
		
	}
	fprintf(fptr, "\n==============================================================================================\n");
	fflush(fptr);
}

void exitprocedure(struct TankStruct *tankcontrol, struct ValveStruct *valvecontrol)
{
	/*
	 * Displays the output of Tank and Valve structures and close the running program.
	 * tankcontrol_current : Current Tank structure compatible with EPANET.
	 * valvecontrol_current : Current Valve structure compatible with EPANET.
	 */

	int temp_count;
	int temp_count2;
	int timeperiod = tankcontrol[0].TimePeriod;

	time_t currentTime;
	currentTime = time(NULL);
	struct tm tm = *localtime(&currentTime);
	char tankFileName[128];
	sprintf(tankFileName, "../../Output_Files/Tank-Output-%d-%d-%d::%d:%d:%d.csv", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);

	char valveFileName[128];
	sprintf(valveFileName, "../../Output_Files/Valve-Output-%d-%d-%d::%d:%d:%d.csv", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);

	FILE *tankPtr;
	FILE *valvePtr;

	tankPtr = fopen(tankFileName, "w");
	valvePtr = fopen(valveFileName, "w");

	if(tankPtr == NULL || valvePtr == NULL) {
		perror("Error in opening the file\n");
		exit(-1);
	}

	// @Ram The below condition was added by @Ram
	if(feasiblity_checker(tankcontrol, valvecontrol)) {
		printf("Simulation feasible\n");
	}
	else
	{
		printf("Simulation Failure :: Not feasible\n");
	}
	fflush(stdout);
	printf("\nOutput Tank File Path -> %s\n", tankFileName);
	fflush(stdout);
	printf("\nOutput Valve File Path -> %s\n", valveFileName);
	fflush(stdout);
	fprintf(fptr, "Simulation Failure :: Not feasible\n");
	fprintf(fptr, "\nOutput Tank File Path -> %s\n", tankFileName);
	fprintf(fptr, "\nOutput Valve File Path -> %s\n", valveFileName); 

	// Display all the Valve Settings and Tank Levels till elapsed time
	// Elapsed time is the time till the valve or demand settings are made
	fprintf(fptr, "\n==============================================================================================\n");
	fflush(fptr);

	// Creating Tank and valve setting final output file
	fprintf(tankPtr, "Time,");
	fprintf(valvePtr, "Time,");
	fflush(tankPtr);
	fflush(valvePtr);

	for(temp_count = 0; temp_count< Ntanks; temp_count++) {
		fprintf(tankPtr, "%s,", tankcontrol[temp_count].TankID);
		fflush(tankPtr);
		if(temp_count == Ntanks-1) {
			fprintf(tankPtr, "\n");
			fflush(tankPtr);
		}
	}
	for(temp_count = 0; temp_count< Nvalves; temp_count++) {
		fprintf(valvePtr, "%s,", valvecontrol[temp_count].ValveID);
		fflush(valvePtr);
		if(temp_count == Nvalves-1) {
			fprintf(valvePtr, "\n");
			fflush(valvePtr);
		}
	}
	fflush(fptr);
	fflush(tankPtr);
	fflush(valvePtr);

	temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
	// @Ram changed timeperiod to timeperiod +1
	for(; temp_count2 < timeperiod+1; temp_count2 += hourly_scheduler) {
		fprintf(tankPtr, "%d,", temp_count2);
		fflush(tankPtr);

		for(temp_count = 0; temp_count< Ntanks; temp_count++) {
			fprintf(tankPtr, "%f,", (tankcontrol[temp_count].TankLevels[temp_count2])*100/tankcontrol[temp_count].MaxTankLevel);
			fflush(tankPtr);
			if(temp_count == Ntanks-1) {
				fprintf(tankPtr, "\n");
				fflush(tankPtr);
			}
		}
	}
	temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
	for(; temp_count2 < timeperiod; temp_count2 += hourly_scheduler) {

		fprintf(valvePtr, "%d,", temp_count2);
		fflush(valvePtr);
		for(temp_count = 0; temp_count< Nvalves; temp_count++) {
			fprintf(valvePtr, "%f,", valvecontrol[temp_count].ValveValues[temp_count2]*100/valvecontrol[temp_count].MaxValue);
			fflush(valvePtr);
			if(temp_count == Nvalves-1) {
				fprintf(valvePtr, "\n");
				fflush(valvePtr);
			}
		}
	}

	fflush(fptr);
	fflush(tankPtr);
	fclose(tankPtr);

	fflush(valvePtr);
	fclose(valvePtr);

	// Tank Output File copy
	char finalOutput[1024];
	sprintf(finalOutput, "cp %s ../../result/tank.csv ", tankFileName);
	int res = system(finalOutput);	

	if(res != 0)
		perror("system command fail");

	// Valve Output File copy
	char valveOutput[1024];
	sprintf(valveOutput, "cp %s ../../result/valve.csv ", valveFileName);
	res = system(valveOutput);	

	if(res != 0)
		perror("system command fail");

	fprintf(fptr, "\n==============================================================================================\n");
	fflush(fptr);
	exit(1);
}

int Queueing_Engine(char *f_job_input)
{
	/*
	 * Reads the jobs in the file for each tank demands 
	 * and valve settings for times 0 to timeperiod-1.
	 * f_job_input 	: file name of file containing job demand data. 
	 * Typical file content: demand|12|15|Bademakan|56.254
	 * 						 valve|1|35|68|118
	 * 						 tank|0|0|Bademaken|150.24
	 */

	FILE* stream = fopen(f_job_input, "rw"); // Open job file for reading
	char keyword[25];
	unsigned short int hours, minutes;
	char id[25]; 
	float value;
	char *tok;
	char line[100];
	PriorityQ p;

	fprintf(fptr, "\n==============================================================================================\n");
	fprintf(fptr, "\nChecking for jobs from the file\n");
	fflush(fptr);

	if (stream == NULL) {
		fprintf(fptr, "Cannot open the file\n");
		fprintf(fptr, "\n==============================================================================================\n");
		fflush(fptr);
		return 0;
	}

	// Reads one line at a time
	while(fgets(line,100,stream) != NULL) {
		tok = strtok(line,"|");

		// Seperate the keyword
		fprintf(fptr, "%s ",tok);
		fflush(fptr);
		strcpy(keyword,tok);

		// Separate the hours
		tok = strtok(NULL ,"|");
		fprintf(fptr, "%s ",tok);
		fflush(fptr);
		hours = atoi(tok);

		// Separate the minutes
		tok = strtok(NULL ,"|");
		fprintf(fptr, "%s ",tok);
		fflush(fptr);
		minutes = atoi(tok);

		// Separate the ID
		tok = strtok(NULL ,"|");
		fprintf(fptr, "%s ",tok);
		fflush(fptr);
		strcpy(id,tok);

		// Separate the value
		tok = strtok(NULL ,"|");
		fprintf(fptr, "%s \n",tok);
		fflush(fptr);
		value = atof(tok);

		// Pushing data into Queue
		Qpush(keyword,hours,minutes,id,value);
	}
	fclose(stream);
	//remove(f_job_input);
	fprintf(fptr, "\n==============================================================================================\n");
	fflush(fptr);

	if(DEBUG) {
		Qdisplay();
	}
	return 1;
}

void Qpush(char *keyword, int hours, int minutes, char *id, float value)
{
	/*
	 * Pushes data into Priority Queue
	 * keyword : keywords like Demand, Tank, Valve
	 * hours : scheduling hours
	 * minutes : scheduling minutes
	 * id : ID of the tank or valve
	 * value : value of tank, valve, demand
	 */

	int i;       /* Function for Insert operation */

	if( Qfull()) {fprintf(fptr, "\n\n Overflow!!!!\n\n");}
	else {
		i=r;
		++r;
		/* Find location for new elem */
		// @Ram Updates the queues according to time
		while(PQ[i].hours >= hours && i >= 0 && !strncmp(PQ[i].keyword,"tank",4)) {
			if(PQ[i].hours == hours && PQ[i].minutes < minutes) {
				break;
			}
			PQ[i+1]=PQ[i];
			i--;
		}
		strcpy(PQ[i+1].keyword,keyword);
		PQ[i+1].hours=hours;
		PQ[i+1].minutes=minutes;
		strcpy(PQ[i+1].id,id);
		PQ[i+1].value=value;
	}
	fflush(fptr);
}

PriorityQ Qpop()
{
	/*
	 * Function for Priority Queue Delete operation 
	 */

	PriorityQ p;

	if(Qempty()) { 
		fprintf(fptr, "\n\nUnderflow!!!!\n\n");
		p.hours=-1;p.minutes=-1;
		return(p); 
	}
	else {
		p=PQ[f];
		f=f+1;
		return(p);
	}
	fflush(fptr);
}

int Qfull()
{
	/*
	 * Function to Check if Queue is Full 
	 */

	if(r==SIZE-1) return 1;
	return 0;
}

int Qempty()
{                    
	/*
	 * Function to Check if Queue is Empty 
	 */

	if(f > r) return 1;
	return 0;
}

void Qdisplay()
{                  
	/*
	 * Function to display status of Queue 
	 */

	int i;

	if(Qempty()) fprintf(fptr, " \n Empty Queue\n");
	else {
		fprintf(fptr, "Front->");
		fflush(fptr);
		for(i=f;i<=r;i++)
			fprintf(fptr, "[%s,%d,%d,%s,%f] ",PQ[i].keyword,PQ[i].hours,PQ[i].minutes,PQ[i].id,PQ[i].value);
		fprintf(fptr, "<-Rear\n");
	}
	fflush(fptr);
}

int Job_Handler(struct TankStruct *tankcontrol, struct ValveStruct *valvecontrol)
{
	/*
	 * Handles the jobs from the queue.
	 * tankcontrol_current : Current Tank structure compatible with EPANET.
	 * valvecontrol_current : Current Valve structure compatible with EPANET.
	 */

	char keyword[25];
	unsigned short int hours, minutes;
	char id[25]; 
	float value;
	int temp_count2;
	int temp_count;
	int simulation_time;
	int time;
	int flag = 1;
	PriorityQ p;

	// Handling Jobs
	while(!Qempty()) {
		p=Qpop(); //Popping the first element from the queue
		/*
		   if(flag == 1){
		   time = p.hours;
		   flag = 0;
		   }

		   if(time == p.hours){
		 */
		if(!strncmp(p.keyword,"demand",6)) {
			fprintf(fptr, "[%s,%d,%d,%s,%f] \n",p.keyword,p.hours,p.minutes,p.id,p.value);
			fflush(fptr);
			for(temp_count2 = 0; temp_count2 < Ntanks; temp_count2++) {
				if(!strncmp(tankcontrol[temp_count2].TankID,p.id,15)) {
					for(temp_count = p.hours; temp_count < tankcontrol[temp_count2].TimePeriod; temp_count++) {
						tankcontrol[temp_count2].OutFlow[temp_count] = p.value;
						// try with this approach as well
						//compute_flows(tankcontrol, valvecontrol, temp_count);
					}
					break;
				}
			}
		}
		else if(!strncmp(p.keyword,"valve",5)) {
			fprintf(fptr, "[%s,%d,%d,%s,%f] \n",p.keyword,p.hours,p.minutes,p.id,p.value);
			fflush(fptr);
			for(temp_count2 = 0; temp_count2 < Nvalves; temp_count2++) {
				if(!strncmp(valvecontrol[temp_count2].ValveID,p.id,2)) {
					for(temp_count = p.hours; temp_count < valvecontrol[temp_count2].TimePeriod; temp_count++) {
						valvecontrol[temp_count2].ValveValues[temp_count] = p.value;
						// try with this approach as well
						//compute_flows(tankcontrol, valvecontrol, temp_count);
					}
					break;
				}
			}
		}
		else if(!strncmp(p.keyword,"tank",4)) {
			fprintf(fptr, "[%s,%d,%d,%s,%f] \n",p.keyword,p.hours,p.minutes,p.id,p.value);
			fflush(fptr);
			for(temp_count2 = 0; temp_count2 < Ntanks; temp_count2++) {
				if(!strncmp(tankcontrol[temp_count2].TankID,p.id,15)) {
					tankcontrol[temp_count2].TankLevels[0] = p.value;
					// try with this approach as well
					//compute_flows(tankcontrol, valvecontrol, temp_count);
					break;
				}
			}
		}
		/*
		   }
		//To be used when simulating for particular time
		else{
		Qpush(p.keyword,p.hours,p.minutes,p.id,p.value);
		break;
		}
		 */
	}
	fflush(fptr);

	if(DEBUG) {
		Qdisplay();
		//Display_Output(tankcontrol, valvecontrol);
	}
	/*
	// Calculate the simulation time. 
	// sim_time = job(seconds) - clock(seconds)
	simulation_time = ((p.hours * 60) + p.minutes) - elapsed_time; 
	simulation_time = time - elapsed_time; // If only hourly manipulations are possible

	// shift the clock to the job time
	elapsed_time = ((p.hours * 60) + p.minutes);
	elapsed_time = time;
	return(simulation_time);
	 */
	return 1;
}

void Job_Scheduler(struct TankStruct *tankcontrol, struct ValveStruct *valvecontrol)
{
	/*
	 * Schedules job to the Queue.
	 * tankcontrol_current : Current Tank structure compatible with EPANET.
	 * valvecontrol_current : Current Valve structure compatible with EPANET.
	 */

	int temp_count2;
	int temp_count;
	float *avg = calloc(Nvalves,sizeof(float));
	float valve_val = 0.0;
	float *std_dev = calloc(Nvalves,sizeof(float));
	int push_avg_flag = 1; // Flag for pushing the average value to the Queue
	int push_set_flag = 1; // Flag for pushing the Valve setting to the Queue

	PriorityQ p;
	/*
	// Scheduling all tank and valve settings as Jobs

	// Initial Tank Settings being pushed to priority queue
	for(temp_count = 0; temp_count < Ntanks; temp_count++){
	Qpush("tank",0,0,tankcontrol[temp_count].TankID,tankcontrol[temp_count].TankLevels[0]);
	printf("[%s,%d,%d,%s,%f] \n","tank",0,0,tankcontrol[temp_count].TankID,tankcontrol[temp_count].TankLevels[0]);
	}

	// Valve settings being pushed to Priority Queue
	for(temp_count = 0; temp_count < Nvalves; temp_count++){
	for(temp_count2 = 0; temp_count2 < valvecontrol[temp_count].TimePeriod; temp_count2 += hourly_scheduler){
	Qpush("valve",temp_count2,0,valvecontrol[temp_count].ValveID,valvecontrol[temp_count].ValveValues[temp_count2]);
	printf("[%s,%d,%d,%s,%f] \n","valve",temp_count2,0,valvecontrol[temp_count].ValveID,valvecontrol[temp_count].ValveValues[temp_count2]);
	}
	}
	 */
	fprintf(fptr, "\n==================================================================\n");
	fflush(fptr);
	// Valve Averaging being done
	for(temp_count = 0; temp_count < Nvalves; temp_count++) {
		for(temp_count2 = 0; temp_count2 < valvecontrol[temp_count].TimePeriod; temp_count2 += hourly_scheduler) {
			avg[temp_count] += valvecontrol[temp_count].ValveValues[temp_count2];
		}
		avg[temp_count] = avg[temp_count]/valvecontrol[temp_count].TimePeriod;
		if(DEBUG) {
			fprintf(fptr, "average valve setting for %s is %f\n",valvecontrol[temp_count].ValveID,avg[temp_count]);
			fflush(fptr);
		}
	}
	fflush(fptr);

	// Calculating std deviation
	for(temp_count = 0; temp_count < Nvalves; temp_count++) {
		for(temp_count2 = 0; temp_count2 < valvecontrol[temp_count].TimePeriod; temp_count2 += hourly_scheduler) {
			std_dev[temp_count] += pow((valvecontrol[temp_count].ValveValues[temp_count2] - avg[temp_count]),2);
		}
		std_dev[temp_count] = sqrt(std_dev[temp_count]/valvecontrol[temp_count].TimePeriod);
		if(DEBUG) {
			fprintf(fptr, "standard deviation of valve %s is %f\n",valvecontrol[temp_count].ValveID,std_dev[temp_count]);
			fflush(fptr);
		}
	}
	fflush(fptr);

	fprintf(fptr, "Jobs to be scheduled\n");
	fflush(fptr);
	// Job file creation
	char fileName[48];
	sprintf(fileName, "../../Output_Files/Job-Output.csv");
	FILE *jptr = fopen(fileName, "w");
	if(jptr == NULL) {
		perror("ERROR in creating Job file :: ");
		exit(-1);
	}

	fprintf(jptr, "Time, Valve ID, Valve Value\n");
	fflush(jptr);

	// Intelligent pushing of Jobs to the Queue
	for(temp_count = 0; temp_count < Nvalves; temp_count++) {
		//printf("\n Std Dev of %s = %f", valvecontrol[temp_count].ValveID,std_dev[temp_count]);
		//printf("\n Avg of %s = %f", valvecontrol[temp_count].ValveID,avg[temp_count]);

		for(temp_count2 = 0; temp_count2 < valvecontrol[temp_count].TimePeriod; temp_count2 += hourly_scheduler) {
			// Main pushing engine
			// @Ram changed 10 to 2
	/*	if(strncmp(valvecontrol[temp_count].ValveID,"3",1)==0)
			{
				valve_val = valvecontrol[temp_count].ValveValues[temp_count2];
				Qpush("valve",temp_count2,0,valvecontrol[temp_count].ValveID,valve_val);
				fprintf(jptr, "%d, %s, %f \n", temp_count2, valvecontrol[temp_count].ValveID,valve_val);
				fflush(jptr);
				if(DEBUG) {
					fprintf(fptr, "[%s,%d,%d,%s,%f] \n","valve",temp_count2,0,valvecontrol[temp_count].ValveID,valve_val);
					fflush(fptr);
					}
				// Setting the appropriate flags
				// @Ram changed the following 
				push_avg_flag = 1;
				push_set_flag = 1;

			}*/
		//else {
			if(std_dev[temp_count] <= 1) {
		//			printf("\n This is std");
				//if(push_avg_flag) {
		//			printf("\n This is std");
					valve_val = avg[temp_count];
					Qpush("valve",temp_count2,0,valvecontrol[temp_count].ValveID,valve_val);
					fprintf(jptr, "%d, %s, %f \n", temp_count2, valvecontrol[temp_count].ValveID,valve_val);
					fflush(jptr);
					if(DEBUG) {
						fprintf(fptr, "[%s,%d,%d,%s,%f] \n","valve",temp_count2,0,valvecontrol[temp_count].ValveID,valve_val);
						fflush(fptr);
					}
					// Setting the appropriate flags
				//	push_avg_flag = 0;
				//	push_set_flag = 1;
					break;
				//}
			}
			
			/*else if(temp_count2==0){
		//			printf("\n This is Avg");
					//valve_val = avg[temp_count];
					valve_val = valvecontrol[temp_count].ValveValues[temp_count2];
					Qpush("valve",temp_count2,0,valvecontrol[temp_count].ValveID,valve_val);
					//Qpush("valve",temp_count2,temp_count,valvecontrol[temp_count].ValveID,valve_val);
					fprintf(jptr, "%d, %s, %f \n", temp_count2, valvecontrol[temp_count].ValveID,valve_val);
					fflush(jptr);
					if(DEBUG) {
						fprintf(fptr, "[%s,%d,%d,%s,%f] \n","valve",temp_count2,0,valvecontrol[temp_count].ValveID,valve_val);
						fflush(fptr);
					}
				}

			if(fabs(valvecontrol[temp_count].ValveValues[temp_count2] - valvecontrol[temp_count].ValveValues[temp_count2+1]) <= 1) {
				// Pushing only the first value as average Valve setting
				//	printf("\n This is strange avg");
		//			printf("\n This is Avg");
				//if(push_avg_flag) {
					// Setting the appropriate flags
					//push_avg_flag = 0;
					//push_set_flag = 1;
				}				
			
			else if(fabs(valvecontrol[temp_count].ValveValues[temp_count2] - valvecontrol[temp_count].ValveValues[temp_count2+1]) > 1) {
				// Pushing the actual Valve setting
				//printf("\n This is normal");
			//		printf("\n This is None");
			//	if(push_set_flag) {
			//		printf("\n This is None");
					valve_val = valvecontrol[temp_count].ValveValues[temp_count2+1];
					Qpush("valve",temp_count2,0,valvecontrol[temp_count].ValveID,valve_val);
					fprintf(jptr, "%d, %s, %f \n", temp_count2+1, valvecontrol[temp_count].ValveID,valve_val);
					fflush(jptr);
					if(DEBUG) {
						fprintf(fptr, "[%s,%d,%d,%s,%f] \n","valve",temp_count2,0,valvecontrol[temp_count].ValveID,valve_val);
						fflush(fptr);
					}
					// Setting the appropriate flags
					// @Ram changed the following 
					//push_avg_flag = 1;
					//push_set_flag = 0;
				}*/	
			

			//@Ram changed 20 to 2
		/*	else if(fabs(valvecontrol[temp_count].ValveValues[temp_count2] - avg[temp_count]) <= 0.0075*valvecontrol[temp_count].MaxValue) {
				// Pushing only the first value as average Valve setting
				//	printf("\n This is strange avg");
		//			printf("\n This is Avg");
				if(push_avg_flag) {
		//			printf("\n This is Avg");
					valve_val = avg[temp_count];
					Qpush("valve",temp_count2,0,valvecontrol[temp_count].ValveID,valve_val);
					//Qpush("valve",temp_count2,temp_count,valvecontrol[temp_count].ValveID,valve_val);
					fprintf(jptr, "%d, %s, %f \n", temp_count2, valvecontrol[temp_count].ValveID,valve_val);
					fflush(jptr);
					if(DEBUG) {
						fprintf(fptr, "[%s,%d,%d,%s,%f] \n","valve",temp_count2,0,valvecontrol[temp_count].ValveID,valve_val);
						fflush(fptr);
					}
					// Setting the appropriate flags
					push_avg_flag = 0;
					push_set_flag = 1;
				}				
			}*/
			else if(fabs(valvecontrol[temp_count].ValveValues[temp_count2] - valve_val) >= 0.001*valvecontrol[temp_count].MaxValue) {
				// Pushing the actual Valve setting
				//printf("\n This is normal");
			//		printf("\n This is None");
				if(push_set_flag) {
			//		printf("\n This is None");
					valve_val = valvecontrol[temp_count].ValveValues[temp_count2];
					Qpush("valve",temp_count2,0,valvecontrol[temp_count].ValveID,valve_val);
					fprintf(jptr, "%d, %s, %f \n", temp_count2, valvecontrol[temp_count].ValveID,valve_val);
					fflush(jptr);
					if(DEBUG) {
						fprintf(fptr, "[%s,%d,%d,%s,%f] \n","valve",temp_count2,0,valvecontrol[temp_count].ValveID,valve_val);
						fflush(fptr);
					}
					// Setting the appropriate flags
					// @Ram changed the following 
					push_avg_flag = 1;
					push_set_flag = 1;
				}	
			}
		//}
		}

		// RESET the appropriate flags
		push_set_flag = 1;
		push_avg_flag = 1;
		valve_val = 0.0;
	}
	fflush(jptr);
	fclose(jptr);
	fprintf(fptr, "\n==================================================================\n");
	fflush(fptr);

	if(DEBUG) {
		// Check the queue
		Qdisplay();
	}

	// Job Output File
	char jobOutput[1024];
	sprintf(jobOutput, "cp %s ../../result/Job-Output.csv", fileName);
	int res = system(jobOutput);	

	if(res != 0)
		perror("system command fail");

	// shift the clock to the first job time
	//elapsed_time = ((p.hours * 60) + p.minutes);
	/*p = Qpop();
	  elapsed_time = p.hours;
	  Qpush(p.keyword, p.hours, p.minutes, p.id, p.value);*/

}

int feasiblity_checker(struct TankStruct *tankcontrol, struct ValveStruct *valvecontrol)
{
	/*
	 * checks if the solution obtained is feasible or not
	 * tankcontrol: The tank structure from the simulation
	 * valvecontrol: The valve structure from the simulation
	 * 
	 * returns 0: if any discrepancy is found
	 * returns 1: if no discrepancy is found
	 * 
	 */

	int temp_count;
	int temp_count2;
	int timeperiod;
	timeperiod = tankcontrol[0].TimePeriod;

//	printf("\n time period = %d",timeperiod);
	fprintf(fptr, "Checking feasibility \n");
	fflush(fptr);
	//Tank Levels at each time for discrepancy
	for(temp_count = 0 ; temp_count < Ntanks; temp_count++) {
		temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
		for(; temp_count2 < timeperiod; temp_count2++) {
			if(tankcontrol[temp_count].TankLevels[temp_count2] <= tankcontrol[temp_count].MinTankLevel 
					| tankcontrol[temp_count].TankLevels[temp_count2] >= tankcontrol[temp_count].MaxTankLevel){
				printf("\n TankLevel[%d::%d::%f] :: MinTankLevel[%f] :: MaxTankLevel[%f]\n", temp_count2, temp_count, tankcontrol[temp_count].TankLevels[temp_count2], tankcontrol[temp_count].MinTankLevel, tankcontrol[temp_count].MaxTankLevel);
				fflush(stdout);
				fprintf(fptr, "TankLevel[%d::%d::%f] :: MinTankLevel[%f] :: MaxTankLevel[%f]\n", temp_count2, temp_count, tankcontrol[temp_count].TankLevels[temp_count2], tankcontrol[temp_count].MinTankLevel, tankcontrol[temp_count].MaxTankLevel);
				fprintf(fptr, "tank Checking feasibility \n");
				fflush(fptr);
				return 0;
			}
		}
	}


//	printf("\n Time Period = %d",timeperiod);
	if(targetMode==TRUE)
	{
	for(temp_count = 0 ; temp_count < Ntanks; temp_count++)
	{
		if(fabs(tankcontrol[temp_count].TankLevels[timeperiod]-solutionTankLevel[0].tankId[temp_count])>(0.02*tankcontrol[temp_count].MaxTankLevel))	{
		printf("\n  %-15s \t Current = %f \t Target = %f",tankcontrol[temp_count].TankID,tankcontrol[temp_count].TankLevels[timeperiod],solutionTankLevel[0].tankId[temp_count]);
		printf("\t Difference percentage = %f",(fabs(tankcontrol[temp_count].TankLevels[timeperiod]-solutionTankLevel[0].tankId[temp_count]))/tankcontrol[temp_count].MaxTankLevel);
		fflush(stdout);
		}
	}
	}
	// @Ram Added these Tank Levels not returning to initial state discrepancy
	for(temp_count = 0 ; temp_count < Ntanks; temp_count++)
	{
		//if((tankcontrol[temp_count].TankLevels[0]> (1.1* tankcontrol[temp_count].TankLevels[timeperiod]))||(tankcontrol[temp_count].TankLevels[0] < (0.9* tankcontrol[temp_count].TankLevels[timeperiod])))
		if(targetMode!=TRUE)
		{
			if(fabs(tankcontrol[temp_count].TankLevels[0]- tankcontrol[temp_count].TankLevels[timeperiod])>(0.02*tankcontrol[temp_count].MaxTankLevel))	{
				//printf("\n Failure success");
				printf("\n TankLevel[%d] :: StartTankLevel[%f] :: EndTankLevel[%f]\n", temp_count, tankcontrol[temp_count].TankLevels[0], tankcontrol[temp_count].TankLevels[timeperiod]);
				fflush(stdout);
				fprintf(fptr, "TankLevel[%d] :: StartTankLevel[%f] :: EndTankLevel[%f]\n", temp_count, tankcontrol[temp_count].TankLevels[0], tankcontrol[temp_count].TankLevels[timeperiod]);
				fprintf(fptr, " Tank Return Checking feasibility \n");
				fflush(fptr);
				return 0;
			}
		}
		else
		{
			if(fabs(tankcontrol[temp_count].TankLevels[timeperiod]-solutionTankLevel[0].tankId[temp_count])>(0.02*tankcontrol[temp_count].MaxTankLevel))	{
//				printf("\n success");
				fprintf(fptr, "TankLevel[%d] :: StartTankLevel[%f] :: EndTankLevel[%f]\n", temp_count, tankcontrol[temp_count].TankLevels[0], tankcontrol[temp_count].TankLevels[timeperiod]);
				fprintf(fptr, " Tank Return Checking feasibility \n");
				fflush(fptr);
				return 0;
			}

		}
	}

// @Ram Added Lines End here
// Valve Values for discrepancy
for(temp_count = 0 ; temp_count <Nvalves; temp_count++) {
	temp_count2 = ((targetMode == TRUE)? currentLevel.timeVal:startTime );
	for(; temp_count2 < timeperiod; temp_count2++) {
		if(valvecontrol[temp_count].ValveValues[temp_count2] > valvecontrol[temp_count].MaxValue){
			fprintf(fptr, "ValveValue[%d::%d::%f] :: MaxValveValue[%f]\n", temp_count2, temp_count, valvecontrol[temp_count].ValveValues[temp_count2], valvecontrol[temp_count].MaxValue);
			fprintf(fptr, "valve Checking feasibility \n");
			fflush(fptr);
			return 0;
		}
	}
}

return 1;
}

float Compute_MaxValue(float d)
{	
	return(MAX_VELOCITY * PI*(d * d)/4.0 / NUM_LITERS_PER_M3);
}


