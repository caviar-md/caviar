
This directory has a simulation script for water-salt separation by means of polarized electrodes.
To this end, we have used the 'water-salt.casl' file which creates a 'o_xyz.xyz' output.
This simulation may take a few days with one processor.
After this simulation, we can post-process the code and extract the potential value inside the simulation box on a grid. For this part, we comment the last line of 'water-salt.casl', meaning the 'run' command.
Then we can create a 3D grid, and use it inside a newly developed class, 'postprocess::Potential_sampler'.
This file reads the snapshots of 'o_xyz.xyz' then it re-simulate one frame. 
By using the simulation result, it samples the potential value in an output file.
The format of the outputted file is:

frameNumber point_index x_index y_index z_index x y z potential_1 potential_2 ... potential_N potential_sum

The outputted file, then, can be introduced for the averaging in another code, 
'caviar/tools/mean_potential_value_on_grid'.

By using this code, one can extract the potential profile in one direction of the 3D space.

