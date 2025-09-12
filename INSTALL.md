# CAVIAR Package Installation Guide

Building the **CAVIAR** package is straightforward and highly configurable.  
You can enable or disable optional features such as MPI, OpenMP, deal.II, Eigen, and muParser through CMake.  

---

## Quick Start (Default Build)

If you just want to build CAVIAR with default settings (Release mode, no extra libraries), run:

```bash
$ mkdir build
$ cd build
$ cmake ${PATH_TO_CAVIAR}
$ make
```

The `CAVIAR` executable will be generated in the build directory.

---

## 1. Create and Enter a Build Directory

It is best practice to build outside of the source directory:

```bash
$ mkdir build
$ cd build
```

---

## 2. Run CMake Configuration

To configure the build, point CMake to the CAVIAR source directory:

```bash
$ cmake ${PATH_TO_CAVIAR}
```

At this stage, default values will be applied (e.g., `Release` build type if none is specified).

---

## 3. Customize Configuration with `ccmake` (Optional but Recommended)

To fine-tune build options, use **ccmake** (or `cmake-gui` if you prefer a graphical interface):

```bash
$ ccmake ${PATH_TO_CAVIAR}
```

Inside `ccmake`, you can configure:

- **Build type**: `Debug`, `Release`, `RelWithDebInfo`, `MinSizeRel`  
- **Parallelization options**:  
  - `CAVIAR_WITH_MPI` → Enable MPI support  
  - `CAVIAR_WITH_OPENMP` → Enable OpenMP parallelization  
  - `CAVIAR_SINGLE_MPI_MD_DOMAIN` → Experimental MPI single-domain mode  
- **External libraries**:  
  - `CAVIAR_WITH_DEALII` / `CAVIAR_WITH_DEALII_MPI` → deal.II support  
  - `CAVIAR_WITH_EIGEN` → Eigen support (`EIGEN_DIR` can be set manually)  
  - `CAVIAR_WITH_MUPARSER` → muParser support (`MUPARSER_DIR` or `MUPARSER_INCLUDE_DIR`)  
- **Debugging**:  
  - `CAVIAR_DEBUG_VERSION` → Extra checks, logging, and assertions  

Press **`c`** to configure, **`g`** to generate, then exit.

---

## 4. Compile the Package

Once configuration is complete, build the package:

```bash
$ make
```

Upon success, an executable named **`CAVIAR`** will be generated in your build directory.

---

## 5. Run an Example

To verify the build, run one of the included examples:

```bash
$ ./CAVIAR < ${PATH_TO_CAVIAR}/examples/e1-two-atoms/e1-two-atoms.fcpp
```

---

## 6. Documentation

For more details on features and usage, consult the **`doc/`** directory.

---
