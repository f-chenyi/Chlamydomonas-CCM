## Modeling Chlamydomonas CCM

### Description
This bundle includes the following packages to simulate a reaction-diffusion model of the Chlamydomonas CO<sub>2</sub>-concentrating mechnism and analyze simulation data:  

>**Engineering** contains the performances and costs of each engineering state, and the general algorithms that compute simple paths from a given starting state to the target state with the best performance. 

>**FullModel** contains the equations and code to run the reaction-diffusion model and return the desired data as output. 

>**ThylakoidStacks** contains algorithms that compute the effective diffusion constant through thyklakoid stacks with realistic geometries.

For a more detailed discription of the model, see our paper [here](https://www.nature.com/articles/s41477-022-01153-7).

### Installation

**ThylakoidStacks** and **FullModel** require installation of FEniCS. Instructions can be found [here](https://fenicsproject.org/download/). A list of packages installed in the virtual environment is stored in **fenics-envs.rtf**. Note that *meshio* and *pygmsh* are NOT installed in fenics by default, but they are required for **ThylakoidStacks**.

**ThylakoidStacks** requires installation of Gmsh, which can be downloaded [here](https://gmsh.info/#Download). Note that pygmsh 6.1.1 and Gmsh 3.0.6 are recommended. Other versions may result in syntax issues. 


### Usage

>**Engineering** 
>
>      load variables.mat
>      Type "gen_distance(BarCode)" to create the distance matrix
>      Type "engineering(dist, NcffList, ATPperCO2fixed, start)" where "start" is the desired starting state to generate engineering paths.

>**FullModel** 
>
>      Type "python ccm_main.py" to run the model for one set of parameters
>      Type "python ccm_scan.py" to run the model for a range of parameters and save the specified output

>**ThylakoidStacks** 
>
>      Type "python -u EffectiveDiffusion_v2.py" to run the algorithm


See the text files in each folder and the comments in the code files for more details.

### Support

If you encounter any issues please contact:

Chenyi Fei (<cfei@princeton.edu>)

Alexandra Wilson (<aw29@princeton.edu>)
