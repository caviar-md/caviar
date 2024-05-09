## Linear Distribution Calculator

### What is this?

This is a post-processing code designed to analyze the output of XYZ file formats containing atomic coordinates. It calculates the linear distribution function for a system of atoms and molecules along a specified direction and length.

### Input Parameters

The code takes the following input parameters in 'main.cpp' file:

- `xyz_filename`: Path to the XYZ file containing atomic coordinates.
- `start`: Starting position for the linear distribution calculation, specified as a 3D vector (x, y, z).
- `direction`: Direction along which the linear distribution is calculated, specified as a 3D vector (x, y, z).
- `totalLength`: Total length of the distribution.
- `numBins`: Number of bins in the length.
- `atomTypeStart`: Start number of the atom type to be averaged on.
- `atomTypeEnd`: End number of the atom type to be averaged on.
- `lim3d.x.min`, `lim3d.x.max`: Limit the atom positions in the interval in the x-direction. `lim3d.x.enable` must be set to true.
- `lim3d.y.min`, `lim3d.y.max`: Limit the atom positions in the interval in the y-direction. `lim3d.y.enable` must be set to true.
- `lim3d.z.min`, `lim3d.z.max`: Limit the atom positions in the interval in the z-direction. `lim3d.z.enable` must be set to true.

### How to Use

1. **Compile**: Run the provided compilation script to compile the code.
   ```sh
   $ ./compile.sh
   ```

2. **Run**: Run the code with the default parameters.
    ```sh
    $ ./linear_distribution
    ```

### Output Files

The code generates output files into a new output folder. This output folder is created automatically by the code. It contains the results and any other generated files.

### Note: Input XYZ Files

Before running the code, ensure that the input XYZ files are copied into the output folder. The code expects to find the input XYZ files in the same directory as the output folder.

### Run with Arguments

Alternatively, you can specify input parameters via command-line arguments. The available command-line options are:

- `-f <filename>`: Specify the path to the XYZ file containing atomic coordinates.
- `-s <x> <y> <z>`: Specify the starting position for the linear distribution calculation.
- `-d <x> <y> <z>`: Specify the direction along which the linear distribution is calculated.
- `-l <length>`: Specify the total length of the distribution.
- `-b <bins>`: Specify the number of bins in the length.
- `-start <start>`: Specify the start number of the atom type range to be averaged on.
- `-end <end>`: Specify the end number of the atom type range to be averaged on.
- `-limx <xmin> <xmax>`: Specify the interval in the x-direction within which particles must be limited.
- `-limy <ymin> <ymax>`: Specify the interval in the y-direction within which particles must be limited.
- `-limz <zmin> <zmax>`: Specify the interval in the z-direction within which particles must be limited.

Here's an example of how to run the code with custom input parameters:

```sh
$ ./linear_distribution -f o_xyz_to_mines.xyz -s 26 0 0 -d -1 0 0 -l 10 -b 50
```

