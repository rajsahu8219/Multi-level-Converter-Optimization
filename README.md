# Multi-level Converter Optimization

This folder stores the Modular Multilevel Converter design optimization code used in the IEEE journal https://ieeexplore.ieee.org/abstract/document/9023972. This journal paper is open sourced and free to download from IEEE site.

This MATLAB design code optimizes the MMC using a population based optimization engine with respect to their mass and loss. In this code, we use MATLAB based Genetic Algorithm Toolbox GOSET, developed by Dr. Scott Sudhoff at Purdue University, and is available for download for free at https://engineering.purdue.edu/ECE/Research/Areas/PES/Software/genetic-optimization-toolbox-2.6.

This code was written ~4 years ago, and is no longer maintained. I am uploading this code on Github for educational purpuses only.

External dependency:
GOSET for GA optimization: https://engineering.purdue.edu/ECE/Research/Areas/PES/Software/genetic-optimization-toolbox-2.6.
Electric Machine Metamodel: for electric machine design, https://ieeexplore.ieee.org/abstract/document/8255662.

Files descriptions:
MMC_Design_v3: runs thw GA optimize and calls the objective function
MMC_Fitness_v3: runs the design simulation to calculate design fitness, mass and loss

<img width="762" alt="Screenshot 2023-02-25 at 9 08 42 AM" src="https://user-images.githubusercontent.com/124555189/221370172-30db8ad8-6d2e-4cb2-9e95-f33c5c18dcb9.png">
<img width="762" alt="Screenshot 2023-02-25 at 9 10 08 AM" src="https://user-images.githubusercontent.com/124555189/221370236-89ecd6ca-cbd5-45fc-bc34-3dbf6b645f83.png">
<img width="762" alt="Screenshot 2023-02-25 at 9 10 43 AM" src="https://user-images.githubusercontent.com/124555189/221370269-effd5d87-bbe1-4c01-9ad0-d945d9efd77a.png">
