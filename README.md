# Neuronal morphology and network topology modulate weak-signal responses in single neurons and small-world networks
## 1.fig2
### Run the code in the fig2 folder to obtain the data for plotting fig2.
## 2.fig3
### Run the code in the fig3 folder to obtain the data for plotting fig3.
## 3.fig4
### Run the code in the fig4 folder: modify the parameters to 0.5, 0.6, and 0.7 respectively, and then run the code to obtain the data for plotting fig4.
```matlab
p_ = 0.7;         % Ratio of dendritic area to total cell area
```
### The following are the random seed numbers and the number of trials for the experiment.
```matlab
%% Run multiple simulations
num_runs = 10;  
seeds = [1,2,3,4,5,6,7,8,9,10];  
```
## 4.fig5
### Run the code in the fig5 folder: modify the parameter to 0.5, 0.6, and 0.7 respectively, 
```matlab
p_ = 0.6;         % Ratio of cell area occupied by the dendrite
```
### and correspondingly adjust the noise intensity to 16.75, 11.50, and 7.25.After running the code, you will obtain the data for plotting fig2 (a1)-(a3).
```matlab
D = 3.75;                                   % noise intensity 
```
### Then correspondingly adjust the noise intensity to 7.0, 4.0, and 2.0.After running the code, you will obtain the data for plotting fig2 (b1)-(b3).
## 5.fig6
### Run the code in the fig6 folder: modify the parameters to 0.4, 0.6, and 0.8 respectively, and then run the code to obtain the data for plotting fig6.
```matlab
A = 0.5;   % A cannot be too large 
```
### The following are the random seed numbers and the number of trials for the experiment.
```matlab
%% Run multiple simulations
num_runs = 10;  
seeds = [1,2,3,4,5,6,7,8,9,10];  
```
## 6.fig7
### Run the code in the fig7 folder: modify the parameters to 0.5 and 0.7 respectively, and then run the code to obtain the data for plotting fig7.
```matlab
p_ = 0.5;         % Ratio of dendritic area to total cell area
```
### The following are the random seed numbers and the number of trials for the experiment.
```matlab
% Run multiple simulations (10 different random seeds)
num_runs = 10;  % Run 10 times
seeds = [1,2,3,4,5,6,7,8,9,10];  % 10 different random seeds
```
## 7.fig8
### Run the code in the fig8 folder: modify the parameters to 0.5, 0.6, and 0.7 respectively, 
```matlab
p_ = 0.5;         % Ratio of dendritic area to total cell area
```
### and for each parameter set, run the code ten times with random seeds set sequentially to 1, 2, 3, 4, 5, 6, 7, 8, 9, and 10 respectively, to obtain the data for plotting fig8.
```matlab
 seed = 1;     % Random seed
 rng(seed);  
```
## 8.fig9
### Run the code in the fig9 folder: modify the parameters to 0.02, 0.04, and 0.06 respectively, 
```matlab
gg_ = 0.04;                  % Conductance between nodes (coupling strength)    
```
### and for each parameter set, run the code ten times with random seeds set sequentially to 1, 2, 3, 4, 5, 6, 7, 8, 9, and 10 respectively, to obtain the data for plotting fig9.
```matlab
 seed = 1;     % Random seed
 rng(seed);  
```
## 9.fig10
### Run the code in the fig10 folder: modify the parameters to 0.4, 0.6, and 0.8 respectively, 
```matlab
Per = 0.4;
```
### and for each parameter set, run the code ten times with random seeds set sequentially to 1, 2, 3, 4, 5, 6, 7, 8, 9, and 10 respectively, to obtain the data for plotting fig10.
```matlab
 seed = 1;     % Random seed
 rng(seed);  
```
## 10.fig11
### Run the code 'untitled1' in the fig11 folder: modify the parameters to 0.4, 0.6, and 0.8 respectively, and then run the code to obtain the data for plotting fig11.
```matlab
Per = 0.80;
```
## 11.fig12
### Run the code in the fig12 folder: modify the parameters to 60, 80, 100 and 200 respectively, 
```matlab
TT = 80;                      % Number of neurons per layer  
```
### and for each parameter set, the value of P (the rewiring probability) ranges from 0 to 1 with a step size of 0.5. Under each value, the code is run 10 times with random seeds set sequentially to 1, 2, 3, 4, 5, 6, 7, 8, 9, and 10 respectively, to obtain the data for plotting fig12.
```matlab
  Per = 0.4;
```
```matlab
 seed = 1;     % Random seed
 rng(seed);  
```
## 12.fig13
### Run the code 'untitled1' in the fig13 to obtain the data, then combine it with the data we obtained from fig12 to plot fig13.



