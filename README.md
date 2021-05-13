# Mechanistic-Histone-Modification-Model

* The file transcriptional_fb_model_gillespieSSA.cpp contains sample C++ code for the Gillespie algorithm implementation of the transcriptional feedback type model of H3K27me3 silencing.

* The Gillespie implementation is programmed as described in Menon and Howard, Investigating histone modification dynamics by mechanistic computational modelling, Histone Methyltransferases: Methods and Protocols, Methods in Molecular Biology.

* The file transcriptional_fb_model_gillespieSSA.cpp contains 3 functions: 

    (1) random_seed(): USED TO GENERATE A SEED FOR THE RANDOM NUMBER GENERATOR

    (2) gillespie_simulation(): USED TO CARRY OUT GILLESPIE SIMULATIONS OF THE MODEL; OUTPUT WRITTEN DIRECTLY TO A TEXT FILE

    (3) main(): ACCEPTS USER SPECIFIED SIMULATION PARAMETERS (number of simulated trajectories AND AND CALLS THE "gillespie_simulation" FUNCTION TO RUN SIMULATIONS OF THE MODEL

* OUTPUT FORMAT (Written to text file): transcriptional_fb_model_gillespieSSA_output.txt

     COLUMN 1: Reaction times
     
     COLUMN 2: H3K27me0 coverage
     
     COLUMN 3: H3K27me1 coverage
     
     COLUMN 4: H3K27me2 coverage
     
     COLUMN 5: H3K27me3 coverage
     
     COLUMN 6: Transcription marker (1 if the current reaction is transcription firing, 0 otherwise)
     
     COLUMN 7: OFF STATE marker (1 if the fraction of H3K27me2/3 > 3PT/4, 0 otherwise)
     
     COLUMN 8: ON STATE marker (1 if the fraction of H3K27me2/3 < PT/4, 0 otherwise)
                                        
