
# Spacial Transcriptomics
## Please note: This Readme is a work in progress..

This repository combines **punkst** (scalable tools for analyzing high-resolution spatial transcriptomics data) and **FICTURE** (fast image-based clustering and transcriptomics utility).  
It provides a unified framework for preprocessing, analyzing, and visualizing spatial transcriptomics data at scale.

## Documentation
- [Project Documentation](https://Yichen-Si.github.io/punkst) (currently focused on punkst; FICTURE docs are being merged)

## Repository Structure
- `src/` – Core source code for punkst  
- `examples/` – Example pipelines and usage for punkst & FICTURE  
- `docs/` – Documentation site configuration (mkdocs)  
- `script/` – Utility and helper scripts  
- `ext/` – External modules and dependencies  
- `DE_preprocessing.py` – Differential expression preprocessing (FICTURE-related)  

## Installation

You can build from source or use Docker.  
Both **punkst** and **FICTURE** are included in this repo.

### From Source

#### Prerequisites
- Git  
- CMake: 3.15 to 3.23  
- C++17 compiler (GCC ≥8, Clang ≥5, MSVC 2017+)  
- TBB, OpenCV  

#### Steps
```bash
# 1) Clone the repository
git clone --recursive https://github.com/your-org/xenium_spacialtranscriptomatics.git
cd xenium_spacialtranscriptomatics

# 2) Create and enter a build directory
mkdir build && cd build

# 3) Configure
cmake ..

# 4) Build
cmake --build . --parallel
````

If you did not clone the submodule (Eigen) initially, run:

```bash
git submodule update --init
```

If TBB is not found, install it via:

* Ubuntu/Debian: `sudo apt-get install libtbb-dev`
* Fedora: `sudo yum install tbb-devel`
* macOS: `brew install tbb`

If installed locally, you may need to specify paths:

```bash
cmake .. \
  -DOpenCV_DIR=$HOME/.local/lib/cmake/opencv4 \
  -DTBB_DIR=$HOME/user/opt/tbb/lib/cmake/tbb \
  -DCMAKE_PREFIX_PATH="$HOME/.local"
```

The `punkst` binary will be placed in `bin/` under the project root.

### Using Docker

Prerequisite: Docker

```bash
docker pull philo1984/punkst:latest
```

Verify installation:

```bash
docker run --rm philo1984/punkst:latest punkst --help
docker run --rm philo1984/punkst:latest ficture --help
```

## Usage

### Running punkst

```bash
punkst/bin/punkst --help
```

### Running FICTURE

```bash
python DE_preprocessing.py --input your_data.h5 --output results/
```

---

## License

This project is licensed under the terms of the LICENSE file included in this repository.

```
