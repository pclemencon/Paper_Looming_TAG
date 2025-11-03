# Perception of mechanical loomings

![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png)
`Computational fluid dynamics`

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Data`
Raw Data: >1Go files from COMSOL simulations (Folder $${\color{green}'Point3D'}$$). The first three columns indicate the (X,Y,Z) coordinates with 0.0475< X <0.0601 , 0<Y<0.015 and 0.05<Z<0.054. The following columns indicate the flow velocities for every time step. Thus, each line gives the flow time profile at a given point of the volume.

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Scripts`
1) Computation_Momentum
- scripts to extract linear momenta from COMSOL flow velocity time series: 1-Create_Dico_Velo.ipynb, 2-Compute_Lienar_Momentum.ipynb
-> returns 12 files giving the total linear momentum signature M(t) for each combination of D and v that are copied in other folders. Example: $${\color{green}'2 Volume R7.5 V5.csv'}$$
  
2) Plots
- script to plot various linear momentum and a scheme of the physical problem: 2_Plot_and_compare_Momentums_v5.ipynb
- script to compare u and high/low Re approximations
-> returns Figure "flows" and "supplementary_flows" of the Paper

![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png)
`Electrophysiology`

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Data`
- data_final.csv with the following structure
| file | v | D | maxFR | timingFR | spikes | stimulus | repet | delta | m | S | pval | R2 | mu |
| - | - | - | - | - | - | - | - | - | - | - | - | - | - |

![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Scripts`
1) 1-influence_D_v
- The script "1-Plots_and_stats_v2.ipynb" creates a layout for the regression plots in the figure 'electrophy + behaviour' of the paper (visualisation).

2) momentum_threshold_analysis
- The script "2-Extract_Fixed_Threshold.ipynb" extracts the valves times. Then, it computes for each file the maxFR and its peak + retrieve the putative momentum threshold. It returns a .csv file data_final.csv

- The script "3-Analyse_Table_Threshold.ipynb" plots the computed threshold as a function of the delays

- A script 'Graphical_abstract.ipynb' creates the comaprison vision/airflow sensing plot

![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png)
 `Behavior (Acheta domestica on a spherical treadmill)` 
 
![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Data`

- Raw data: videos of the cricket saved as multipage .tif files with Basler.  
- A table (Treadmill_clean.xlsx, Treadmill_Repet_Size.xlsx) with the following structure

| species | video | grillon | v | safe | mvt | D | comments | keep | - | - | frame_immo | frame_reac | frame_end | tps_reaction | tps_reaction_norm |
| - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - |
| acheta | file_name | animal ID | v | touched? | moved? | D | - | 0=discard, 1=keep, 2=repet | - | - | start | reaction | end | - | < or > 1 |


![#1589F0](https://placehold.co/15x15/1589F0/1589F0.png) `Scripts`

1) 1-Influence_D_v
- The script "Analyse_Treadmill_Data.ipynb" analyses the behavioural data (videos of Acheta on the spherical treadmill)
- It analyses either the file "Sizes" (N=9) or the file "Velocities" (N=78)

2) 2-Compute_Delays
- The script "Compute_delay_behav.ipynb" computes the behavioural delay and plot the threshold momentum as a function of the delay
- It requires the "Treadmill_clean.xlsx" file and the "Volume" folder.
- It generates a .csv file (delay_behaviour.csv)

-> scripts for making the plots combining behaviour and electrophysiology are located in the "Electrophysiologcal_Data" folder


%> [!NOTE]
%>j
