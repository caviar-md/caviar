
Water salt examples.


Resume instructions:
To extract last frame of an xyz file, one can use the following shell command:
(in case there's 81 atoms, extract 81+2 lines of the file from the end)
$ tail -n 83 o_xyz.xyz > o_xyz_to_resume.xyz
Also. one can use the complete file by adding 'last_frame' option to 'adata add_xyz_data_file' CASL command.
In this case, the last frame will be automatically extracted.


In order to use the file, add the following command at the end of the atomdata in the CASL file:
adata add_xyz_data_file file_name "o_xyz_to_resume.xyz" read_velocity replace_data

'replace_data' is the 
