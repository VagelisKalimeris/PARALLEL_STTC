/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: cond_null_dist.hpp                                               *
*                                                                             *
*******************************************************************************
******************************************************************************/


using namespace std;

/******************************************************************************
* FUNCTION NAME: sign_trpl_limit                                              *
*                                                                             *
* ARGUMENTS: Two neuron's timelines(references to vectors), the total number  *
*             of time samples recorded(int), and a time interval(int).        *
*                                                                             *
* PURPOSE: Calculates the number of firing events of ‘reduced A’              *
*                                                                             *
* RETURNS: A number(int) of events.                                           *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
int sign_trpl_limit(const int time_line_A[], int time_line_A_size, 
                        const int time_line_C[], int time_line_C_size, int Dt);
