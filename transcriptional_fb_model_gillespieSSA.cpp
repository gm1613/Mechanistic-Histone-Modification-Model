//
//  transcriptional_fb_model_gillespieSSA.cpp
//  
//
//  Created by Govind Menon (JIC) on 27/02/2021.
//


#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <fstream>
#include <sys/time.h>
#include <stdio.h>

using namespace std;

// SET UP FUNCTION THAT USES THE SYSTEM CLOCK TO GENERATE A SEED FOR THE RANDOM NUMBER GENERATOR
unsigned long int random_seed()
{
    struct timeval tv;
    gettimeofday(&tv,0);
    return (tv.tv_sec + tv.tv_usec);
}







int gillespie_simulation(int nsims, int ncycles, string file_name)
{
    // FUNCTION TO CARRY OUT GILLESPIE SIMULATIONS OF THE STOCHASTIC HISTONE MODIFICATION MODEL: OUTPUT WRITTEN DIRECTLY TO TEXT FILE
    /* INPUTS TO THE FUNCTION:
     nsims: NUMBER OF SIMULATED TRAJECTORIES
     ncycles: DURATION OF SIMULATION EXPRESSED AS NUMBER OF CELL CYCLES
     OUTPUT FORMAT (Written to text file):
     COLUMN 1: Reaction times
     COLUMN 2: H3K27me0 coverage
     COLUMN 3: H3K27me1 coverage
     COLUMN 4: H3K27me2 coverage
     COLUMN 5: H3K27me3 coverage
     COLUMN 6: Transcription marker (1 if the reaction is transcription firing, 0 otherwise)
     COLUMN 7: OFF STATE marker (1 if the fraction of H3K27me2/3 > 3PT/4, 0 otherwise)
     COLUMN 8: ON STATE marker (1 if the fraction of H3K27me2/3 < PT/4, 0 otherwise)
     */
    
    // SET OUTPUT PRECISION
    cout.precision(10);
    // INITIALISE FLAG RETURNED BY FUNCTION - RETURNS 0 IF SIMULATION COMPLETED
    int flag=1;
    
    // SET UP RANDOM NUMBER GENERATOR
    gsl_rng *r=gsl_rng_alloc(gsl_rng_mt19937); // Using the Mersene Twister 19937 algorithm from the GNU Scientific Library for pseudo-random number generation
    gsl_rng_set(r, random_seed()); // Initialise random number generator with seed generated using system clock
    
    // SPECIFYING SYSTEM PARAMETERS
    // LOCUS SIZE
    const int N = 60; //Number of H3 histones at the locus
    const int nrxn = 2*N+1; // Total number of reactions = N methylation reactions + N demethylation reactions + transcription
    
    // REACTION PARAMETERS FOR THE MODEL: THESE PARAMETERS DETERMINE THE REACTION PROPENSITIES FOR A GIVEN STATE OF THE SYSTEM
    const float kme = 8e-6; // PRC2 mediated me2 to me3 methylation rate (per H3 histone per second)
    const float fmax = 4e-3; // Maximum transcription firing rate - with no H3K27me2/me3 coverage (per second)
    const float pdem =0.004; // Demethylation per histone during transcription (per H3 histone per transcription event at the locus) {0.004 default, 0.04 for loss of silencing example}
    const float pex =1e-3; // Histone exchange proability per histone during transcription (per H3 histone per transcription event at the locus)
    const float pinherit = 0.5; // Histone inheritance probability at replication (per H3 histone per replication event at the locus)
    const float PT = 1.f/3; // K27 me2/me3 repression threshold
    
    const float kme01 = 9*kme; // PRC2 mediated me0 to me1 methylation rate (per H3 histone per second)
    const float kme12 = 6*kme; // PRC2 mediated me1 to me2 methylation rate (per H3 histone per second)
    const float gme01 = kme01/20; // Noisy methylation rate: me0 to me1 (per H3 histone per second)
    const float gme12 = kme12/20; // Noisy methylation rate: me1 to me2 (per H3 histone per second)
    const float gme23 = kme/20;   // Noisy methylation rate: me2 to me3 (per H3 histone per second)
    const float rho_me2 = 0.1; // Relative activation of PRC2 by H3K27me2
    const float fmin = 1e-4; // Minimum transcription firing rate - with H3K27me2/me3 coverage above the threshold (per second)
    const float gdem = fmin*pdem; // Noisy demethylation rate (per H3 histone per second)
    const float cycletime = 22*3600; // Cell-cycle length in seconds (Assuming 22hr cell-cycle)
    
    const float alpha = 1; // Trans-acting gene activation parameter
    const float beta = 1; // Relative local PRC2 activity parameter
    
    
    // SPECIFYING INITIAL STATE OF THE LOCUS
    int init_K27me_state = 3; // Set init_K27me_state = 3 to start loci in the full me3 modified state, and init_K27me_state = 0 to start loci in the full me0 covered state
    
    
    
    
    
    
    // SIMULATION PARAMETERS
    double tmax = ncycles*cycletime - 1; // Time limit for simulation (Simulation set to end one second before the nth replication event)

    
    //  DECLARING SYSTEM VARIABLES
    // TIME
    double t;       // Current reaction time
    int cyclecount; //Cell-cycle counter
    // STATE VARIABLES
    int *cdomK27 = new int[N]; // Array S (named cdomK27 here) representing H3K27 methylation states across the locus: {0,1,2,3} for K27me0/1/2/3
    float PK27me23; //Proportion of K27 me2/me3 marked histones
    float Ei; // Neighbour contributions to PRC2 activity
    float PK27me0, PK27me1, PK27me2, PK27me3; //Proportion of K27 me0/me1/me2/me3 marked H3 histones
    
    
    // DECLARING SIMULATION VARIABLES
    //RANDOM NUMBERS USED IN ALGORITHM
    double Xi1,Xi2, tau; //
    int reploss,transcription,histexc;
    //REACTION PROPENSITIES
    double totalprop;
    double *prpnsty= new double[nrxn]; // Array of propensities (N methylation propensities, followed by N demethylation propensities, and transcription propensity last)
    double *cmprpnsty= new double[nrxn]; // Arrayof cumulative propensities
    
    
    
    //OPEN OUTPUT FILE
    ofstream file_;
    file_.open(file_name);
    
    // START OF SIMULATIONS
    
    if(file_.is_open())
    {
        // LOOP OVER NUMBER OF SIMULATIONS (nsims with initial uniform me0 and nsims with initial uniform me3)
        for(int j=1;j<=nsims;j++)
        {
            
            //SETTING INITIAL TIME AND SYSTEM STATE (INITIALIZE AND RECORD)
            t=0;
            cyclecount=1;
            
            
            // Initialize system state
            PK27me0 = 0;
            PK27me1 = 0;
            PK27me2 = 0;
            PK27me3 = 0;
            for(int i=0;i<N;i++)
            {
                cdomK27[i] = init_K27me_state;
                PK27me0 += (cdomK27[i]==0);
                PK27me1 += (cdomK27[i]==1);
                PK27me2 += (cdomK27[i]==2);
                PK27me3 += (cdomK27[i]==3);
                
            }
            PK27me0/= N;
            PK27me1/= N;
            PK27me2/= N;
            PK27me3/= N;
            PK27me23 = PK27me2 + PK27me3;
            // WRITE CURRENT STATE TO OUTPUT FILE (IN THIS EXAMPLE WE WRITE JUST THE FRACTIONAL COVERAGES, THE TRANSCRIPTION INDICATOR, AND OFF/ON STATE INDICATORS)
            file_ << t << " " << PK27me0 << " " << PK27me1 << " " << PK27me2 << " " << PK27me3 << " " << 0 << " "<< (PK27me23>(0.75*PT)) << " "<< (PK27me23<(0.25*PT)) << "\n"; // 8 columns
            
            
            // THE WHILE LOOP BELOW ITERATIVELY PERFORMS THE STEPS TO SIMULATE A SINGLE TRAJECTORY
            while(t<=tmax)
            {
                
                // (STEP 1 IN THE TEXT) COMPUTING ALL REACTION PROPENSITIES
                 //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                // COMPUTING K27 METHYLATION PROPENSITIES: NOTE THAT THE ARRAY cdomK27 is indexed from 0 (first histone) to N-1 (last histones)
                for(int i=0;i<N;i++)
                {   if(i%2==0)
                {   if(i==0) //Boundary histone
                   { Ei = rho_me2*(cdomK27[i+1]==2)+(cdomK27[i+1]==3)
                    +rho_me2*(cdomK27[i+2]==2)+(cdomK27[i+2]==3)
                    +rho_me2*(cdomK27[i+3]==2)+(cdomK27[i+3]==3);}
                    else if(i==N-2) //Boundary histone
                    { Ei = rho_me2*(cdomK27[i-2]==2)+(cdomK27[i-2]==3)
                    +rho_me2*(cdomK27[i-1]==2)+(cdomK27[i-1]==3)
                    +rho_me2*(cdomK27[i+1]==2)+(cdomK27[i+1]==3);}
                    else //Interior histone
                    {Ei = rho_me2*(cdomK27[i-2]==2)+(cdomK27[i-2]==3)
                    +rho_me2*(cdomK27[i-1]==2)+(cdomK27[i-1]==3)
                    +rho_me2*(cdomK27[i+1]==2)+(cdomK27[i+1]==3)
                    +rho_me2*(cdomK27[i+2]==2)+(cdomK27[i+2]==3)
                    +rho_me2*(cdomK27[i+3]==2)+(cdomK27[i+3]==3);}
                }
                else  //i odd corresponds to even numbered histones
                {   if(i==1) //Boundary histone
                    { Ei = rho_me2*(cdomK27[i-1]==2)+(cdomK27[i-1]==3)
                    +rho_me2*(cdomK27[i+1]==2)+(cdomK27[i+1]==3)
                    +rho_me2*(cdomK27[i+2]==2)+(cdomK27[i+2]==3);}
                    else if(i==N-1) //Boundary histone
                    { Ei = rho_me2*(cdomK27[i-3]==2)+(cdomK27[i-3]==3)
                    +rho_me2*(cdomK27[i-2]==2)+(cdomK27[i-2]==3)
                    +rho_me2*(cdomK27[i-1]==2)+(cdomK27[i-1]==3);}
                    else //Interior histone
                    {Ei = rho_me2*(cdomK27[i-3]==2)+(cdomK27[i-3]==3)
                    +rho_me2*(cdomK27[i-2]==2)+(cdomK27[i-2]==3)
                    +rho_me2*(cdomK27[i-1]==2)+(cdomK27[i-1]==3)
                    +rho_me2*(cdomK27[i+1]==2)+(cdomK27[i+1]==3)
                    +rho_me2*(cdomK27[i+2]==2)+(cdomK27[i+2]==3);}
                }
                    
                    prpnsty[i] = beta*((cdomK27[i]==0)*(gme01+kme01*Ei) + (cdomK27[i]==1)*(gme12+kme12*Ei) + (cdomK27[i]==2)*(gme23+kme*Ei));
                }
                
                
                
                
                
                // COMPUTING K27 DEMETHYLATION PROPENSITIES and PROPORTION OF K27me2/3 MODIFICATIONS
                PK27me23=0;
                for(int i=0;i<N;i++)
                {
                    prpnsty[i+N]=gdem*(cdomK27[i]>0);
                    PK27me23 += (cdomK27[i]>1);
                }
                
                PK27me23/=N;
                
                
                
                // COMPUTING TRANSCRIPTION PROPENSITY
                
                if(PK27me23<PT)
                {   prpnsty[2*N] = alpha*fmax*(1-(PK27me23/PT)*(1-(fmin/fmax)));}
                else
                { prpnsty[2*N] = alpha*fmin;}
               //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                
                
                
                // (STEP 2 IN THE TEXT) CALCULATING TOTAL PROPENSITY AND ARRAY OF CUMULATIVE PROPENSITIES
                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                totalprop=0;
                for(int k=0;k<nrxn;k++)
                {
                    totalprop+=prpnsty[k];
                    cmprpnsty[k]=totalprop;
                }
                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                // (STEP 3 IN THE TEXT) COMPUTE WAITING TIME UNTIL NEXT REACTION
                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Xi1=gsl_rng_uniform(r);
                tau = -log(Xi1)/totalprop;
                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                // (STEP 4 IN THE TEXT) UPDATE TIME AND DECIDE WHETHER THE NEXT STEP SHOULD BE STEP 5, STEP 6, OR STEP 7
                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                t+=tau;
                if(t>tmax)// CHECK IF THE ALGORITHM SHOULD TERMINATE
                {
                    // (STEP 5 IN THE TEXT) RECORD CURRENT STATE AND TIME (TIME LIMIT OF SIMULATION)
                    // WRITE CURRENT STATE TO OUTPUT FILE (IN THIS EXAMPLE WE WRITE JUST THE FRACTIONAL COVERAGES, THE TRANSCRIPTION INDICATOR, AND OFF/ON STATE INDICATORS)
                    file_ << tmax << " " << PK27me0 << " " << PK27me1 << " " << PK27me2 << " " << PK27me3 << " " << 0 << " "<< (PK27me23>(0.75*PT)) << " "<< (PK27me23<(0.25*PT)) << "\n"; // 8 columns
                    break; // Simulation terminates BEFORE replication for the current cycle is carried out
                }
                
                if(t>cyclecount*cycletime) //CHECK IF THE NEXT REPLICATION SHOULD OCCUR BEFORE THIS TIME
                {
                    // (STEP 6 IN THE TEXT) DNA REPLICATION AT THE LOCUS - RANDOM PARTITIONING OF NUCELOSOMES BETWEEN DAUGHTER STRANDS
                    // REPLICATION: UPDATE TIME
                    t=cyclecount*cycletime; //TIME SET TO REPLICATION TIME
                    // REPLICATION: UPDATE SYSTEM STATE
                    
                    for(int i=0;i<N;i+=2)
                    {
                        
                        reploss = (gsl_rng_uniform(r)<pinherit);//Set reploss=0 if gsl_rng_uniform(r)>pinherit: for pinherit = 0.5, this results in a 50% probability of inheritance on average
                        
                        cdomK27[i]=reploss*cdomK27[i];    // We assume ihneritance of intact nucleosomes, so both the modification state of both histones must be updated
                        cdomK27[i+1]=reploss*cdomK27[i+1];
                        
                        
                        
                    }
                    // REPLICATION: RECORD CURRENT TIME AND SYSTEM STATE AFTER REPLICATION
                    PK27me0 = 0;
                    PK27me1 = 0;
                    PK27me2 = 0;
                    PK27me3 = 0;
                    for(int i=0;i<N;i++)
                    {
                        PK27me0 += (cdomK27[i]==0);
                        PK27me1 += (cdomK27[i]==1);
                        PK27me2 += (cdomK27[i]==2);
                        PK27me3 += (cdomK27[i]==3);
                        
                    }
                    PK27me0/= N;
                    PK27me1/= N;
                    PK27me2/= N;
                    PK27me3/= N;
                    
                    
                    // WRITE CURRENT STATE TO OUTPUT FILE (IN THIS EXAMPLE WE WRITE JUST THE FRACTIONAL COVERAGES, THE TRANSCRIPTION INDICATOR, AND OFF/ON STATE INDICATORS)
                    file_ << t << " " << PK27me0 << " " << PK27me1 << " " << PK27me2 << " " << PK27me3 << " " << 0 << " "<< (PK27me23>(0.75*PT)) << " "<< (PK27me23<(0.25*PT)) << "\n"; // 8 columns
                    
                    
                    cyclecount++; // CELL CYCLE COUNTER UPDATE
                }// END OF REPLICATION BLOCK
                
                else //IF NOT REPLICATION, CHOOSE NEXT REACTION TO OCCUR BASED ON RELATIVE PROPENSITIES AND UPDATE THE SYSTEM STATE
                {
                    Xi2=gsl_rng_uniform(r);
                    transcription = 1; //INITIALIZE TRANSCRIPTION FLAG
                    for(int k=0;k<2*N;k++) //ITERATE OVER ARRAY OF CUMULATIVE PROPENSITIES UP TO BUT NOT INCLUDING TRANSCRIPTION : k=0 to 2*N-1
                    {
                        if(Xi2<(cmprpnsty[k]/totalprop)) // IF ONE OF THE METHYLATION/DEMETHYLATION REACTIONS IS CHOSEN
                        {   if(k<N)
                            {cdomK27[k]+=1;}
                            else
                            {cdomK27[k-N]-=1;}
                            transcription=0; //TRANSCRIPTION FLAG SET TO ZERO IF ONE OF THE METHYLATION/DEMETHYLATION REACTIONS IS CHOSEN
                            break;
                        }
                    }
                    if(transcription==1) // IF NONE OF THE METHYLATION/DEMETHYLATION REACTIONS ARE CHOSEN, A "TRANSCRIPTION EVENT" OCCURS
                    {
                        //HISTONE K27 DEMETHYLATION BY TRANSCRIPTION (AS SPECIFIED IN STEP 7 IN THE TEXT)
                        for(int i=0;i<N;i++)
                        {
                            cdomK27[i]-=(cdomK27[i]>0)*(gsl_rng_uniform(r)<pdem);
                            
                        }
                        //NUCLEOSOME EXCHANGE BY TRANSCRIPTION (AS SPECIFIED IN STEP 7 IN THE TEXT)
                        for(int i=0;i<N;i+=2) // ITERATE OVER PAIRS OF HISTONES CORRESPONDING TO NUCLEOSOMES
                        {
                            histexc = (gsl_rng_uniform(r)>(2*pex)); //This sets histexc=0 if gsl_rng_uniform(r) < 2*pex. We need a factor of 2 here as the iteration is over pairs of histones.
                            cdomK27[i]=histexc*cdomK27[i];
                            cdomK27[i+1]=histexc*cdomK27[i+1];
                            
                        }
                    }
                    
                    // RECORD CURRENT TIME AND SYSTEM STATE AFTER REACTION
                    PK27me0 = 0;
                    PK27me1 = 0;
                    PK27me2 = 0;
                    PK27me3 = 0;
                    for(int i=0;i<N;i++)
                    {
                        PK27me0 += (cdomK27[i]==0);
                        PK27me1 += (cdomK27[i]==1);
                        PK27me2 += (cdomK27[i]==2);
                        PK27me3 += (cdomK27[i]==3);
                        
                    }
                    PK27me0/= N;
                    PK27me1/= N;
                    PK27me2/= N;
                    PK27me3/= N;
                    // WRITE CURRENT STATE TO OUTPUT FILE (IN THIS EXAMPLE WE WRITE JUST THE FRACTIONAL COVERAGES, THE TRANSCRIPTION INDICATOR, AND OFF/ON STATE INDICATORS)
                    file_ << t << " " << PK27me0 << " " << PK27me1 << " " << PK27me2 << " " << PK27me3 << " " << (transcription==1) << " "<< (PK27me23>(0.75*PT)) << " "<< (PK27me23<(0.25*PT)) << "\n"; // 8 columns
                    
                }
                
                
                
            }//WHILE LOOP OVER SINGLE TRAJECTORY
        }//FOR LOOP OVER NUMBER OF SIMULATIONS
        
        // CLOSE OUTPUT FILE ONCE ALL SIMULATIONS ARE COMPLETE
        file_.close();
    }
    
    
    //FREE ALLOCATED MEMORY
    gsl_rng_free(r);
    delete [] cdomK27;
    delete [] prpnsty;
    delete [] cmprpnsty;
    
    flag=0;
    return flag;
    
}

// MAIN FUNCTION ACCEPTS USER SPECIFIED SIMULATION PARAMETERS AND CALLS THE "gillespie_simulation" FUNCTION TO RUN SIMULATIONS OF THE MODEL
/* USER SPECIFIED PARAMETERS
 nsims: number of simulated trajectories
 ncycles: duration of simulation expressed as number of cell cycles
 */
int main(int argc, char *argv[]) // generating
{
    // SET UP TIMER
    clock_t start,end;
    double cpu_time_used;
    
    // EXTRACT SIMULATION PARAMETERS FROM USER INPUT
    int nsims = atoi(argv[1]);
    int ncycles = atoi(argv[2]);
    
    // SPECIFY OUTPUT FILE NAME
    string file_name = "transcriptional_fb_model_gillespieSSA_output.txt";
    
    int flag;
    
    
    start=clock(); // START CLOCK TO MEASURE RUNNING TIME
    
    // CALL THE FUNCTION "gillespie_simulation" TO RUN SIMULATIONS, WITH USER SPECIFIED SIMULATION PARAMETERS AS INPUT
    flag = gillespie_simulation(nsims, ncycles, file_name);
    
    end=clock(); // STOP CLOCK
    
    //CALCULATE AND OUTPUT RUNNING TIME
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    cout << "Time used:" << cpu_time_used << "\n";
    
    
    
    
}
