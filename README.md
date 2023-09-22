## Introduction

**CAVIAR** is a versatile molecular dynamics (MD) package meticulously crafted by the skilled hands of Morad Biagooi and Ehsan Nedaaee Oskoee. While its primary mission is the simulation of soft materials at the nano scale, it boasts a broader range of capabilities. CAVIAR excels in the realm of handling intricate geometries for MD simulations.

The heart and soul of CAVIAR are coded in C++, adhering to the C++14 standard. Compilation configurations are orchestrated seamlessly with the aid of CMake. To complement its functionality, a selection of bash shell scripts has been devised for pre and post-processing tasks.

## CAVIAR Scripting Language (CASL)

Within CAVIAR resides a built-in scripting language, known as **CASL** (CAVIAR Scripting Language). It follows a simple and intuitive one-command-per-line structure, making it accessible to both novices and experts. Comprehensive documentation and practical examples are at your disposal to grasp the intricacies of CASL.

## Getting Started

Embarking on your CAVIAR journey is as straightforward as it gets. Begin by crafting your simulation script based on a working example. Dive into the documentation, explore the provided resources, and embark on the path of scripting your simulations with confidence.

## OS and Library Compatibility

CAVIAR flourishes in the Debian-based Linux ecosystem. Its list of dependencies is concise and includes:

- deal.II (Versions 8.5.1, 9.2, 9.4) library for finite element support, with adherence to C++17 standards.
- muparser, maintaining compatibility with the corresponding version used in deal.II.

It is essential to note that deviating from these specific compiler or library versions might lead to compatibility issues. Versions lower than the specified ones may not function seamlessly, especially if they lack comprehensive support for C++17 standards. On the other hand, higher versions may or may not work, depending on potential deprecations of the utilized commands. We warmly welcome any insights or information regarding these aspects for inclusion in our documentation.

## License

CAVIAR operates under the GNU Lesser General Public License as promulgated by the Free Software Foundation. This license provides you with the freedom to choose either version 3.0 or any subsequent version at your discretion. It is vital to bear in mind that:

*The program is distributed without warranty, as permitted by applicable law. Except as expressly provided otherwise in writing, the copyright holders and other parties offer the program "as is" without any warranties, whether expressed or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose. The entire risk as to the quality and performance of the program rests with you. In the event the program proves to be defective, you are responsible for all necessary servicing, repair, or correction.*

## Contact Information

For inquiries, assistance, or collaboration, please feel free to reach out to us:

- Morad Biagooi (m.biagooi .at. gmail.com)
- Ehsan Nedaaee Oskoee (nedaaee .at. iasbs.ac.ir - nedaaee .at. gmail.com)

## Reporting Bugs

In the event that you encounter any challenges during the installation or utilization of CAVIAR, we encourage you to report your concerns to 'm.biagooi .at. gmail.com.' We are committed to swiftly resolving any issues you encounter. There are three primary categories of bugs you might encounter:

1. **Compile Time Bugs**: Ideally, no errors or bugs should surface during the compilation of any release version of CAVIAR. If such issues arise, it is likely tied to discrepancies in compiler or library versions.

2. **Run Time Bugs**: If CAVIAR's interpreter reports an 'ERROR' message, please review your script for potential errors. In the case of a segmentation fault, you can employ gdb to investigate the core file and pinpoint the issue. Once identified, don't hesitate to contact us, providing your script for assistance. Another type of bug may manifest as system freezes during an MPI run.

3. **Result Bugs**: Should your simulation results deviate from expectations in a standard test simulation, there may be parameter configuration issues in your scripts. Additionally, there's a slim possibility that the implemented algorithm might not fully support the physics of your specific problem.

## Citation

If you use CAVIAR in your research or work, please consider citing the following publication:

Biagooi, M., Samanipour, M., Ghasemi, S. A., & Nedaaee Oskoee, S. (2020). CAVIAR: A simulation package for charged particles in environments surrounded by conductive boundaries. AIP Advances, 10(3).