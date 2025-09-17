# Dynamic-System-Imitation-Fixed-Wing-dynamics-on-a-quadcopter-with-gimbal
Short simulation GitHub repository to test the accuracy of a quadcopter with gimbal, flying using only fixed-wing controller inputs.

To run a test simulation, run the 'Rrun_fw_qc_compare.m' function with the relevant inputs for the desired flight. 
If you wish to run multiple runs with different flight patterns, use 'batch_run_fw_sequences.m', the sequence can be set using a setup like this:
seqList = {
  [3 +2; 2 +5; 4 -4; 3 +4; 2 -2; 4 +2; 4 -5; 2 -1]
};
Where the first number is the control parameter, 1=thrust, 2=elevator, 3=aileron, 4=rudder, and the second number is the strength and direction of the deflection or thrust. Currently setup to do strength*2 for deflection (+5 would give a 10 degree deflection etc.
More introduction is placed in the batch script.
For longer flights use the run_60s_seq_noThrust.m, and plan the input controls as mentioned.

If you are just here to look at pretty plots, check in the 'Simulations' folder. Specifc MAE, and MAX values are in the .txt file included in each simulation.


TODO: Add instructions for setting up simulink, as this takes a little work, and correct ordering. If in doubt use the above method.