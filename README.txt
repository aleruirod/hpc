To compile the code you can use the "Makefile" within the project directory.
The scheduler script file for Delta is also in the same directory.

1. Use the "make" command to create the executable object.

2. Run the binary file with the command "make run". (You can change the number of processes utilised in the "Makefile").

3. In order to run the application in Delta, load the modules required for running mpi.

4. Add the job to the queue with the command "qsub scheduler.script.sub"

Note: You can change the deltaT and deltaX values in the "Solution.h" header file.
Note 2: Inorder to export the results in dat format for plotting the export functions must be uncommented in the "main.cpp" file.