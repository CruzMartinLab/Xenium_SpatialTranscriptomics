
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

Perfect 👍 Adding this to your **README.md** for Linux users will make it much easier for others to reproduce what you did. I’ll adapt your PR text into a polished section you can drop directly into your README under something like **“Build Instructions (Linux)”**.

---

## 🚀 Build Instructions (Linux)

This section documents the steps to successfully build and run the `punkst` binary from source on Linux (or WSL). It also covers common pitfalls and troubleshooting tips.

### 1. Configure CMake
Run CMake from a separate `build/` directory:

```bash
mkdir build && cd build
cmake ..
````

> **Note:** If you have dependencies installed in non-standard paths, you may need to specify them manually:
>
> ```bash
> cmake .. \
>   -DZLIB_LIBRARY=/path/to/zlib/libz.so \
>   -DZLIB_INCLUDE_DIR=/path/to/zlib/include \
>   -DTBB_DIR=/path/to/tbb/lib/cmake/tbb \
>   -DOpenCV_DIR=/path/to/opencv/lib/cmake/opencv4
> ```

### 2. Build the Project

Compile using:

```bash
cmake --build . --parallel
```

This produces the `punkst` binary inside:

```
./bin/punkst
```

### 3. Verify File Permissions

Ensure the binary is executable:

```bash
ls -l ./bin/punkst
# -rwxr-xr-x ...
```

If not, run:

```bash
chmod +x ./bin/punkst
```

### 4. Run the Program

Test the installation with:

```bash
./bin/punkst
```

You should see output like:

```
Licensed under the CC BY-NC 4.0
To run a specific command      : ./bin/punkst [command] [options]
For detailed instructions, run : ./bin/punkst --help
```

### 5. Available Commands

Some of the supported commands include:

* `--test` (Test run)
* `--pts2tiles` (Assign points to tiles)
* `--pts2tiles-binary`
* `--tiles2hex` (Convert points in tiles to hexagons)
* `--lda4hex` (Train LDA model)
* `--topic-model` (Train LDA/HDP model)
* `--pixel-decode` (Decode pixel-level data)
* `--draw-pixel-factors`
* `--cooccurrence`
* `--merge-mtx`
* `--coloc2markers`
* `--convert-dge`
* `--multisample-prepare`
* `--merge-units`
* `--multi-conditional-test`
* `--nmf-pois-log1p`
* `--draw-pixel-features`

### 6. Example Usage

```bash
# Run tests
./bin/punkst --test

# Get help for a specific command
./bin/punkst --pts2tiles --help
```

---

### 📝 Notes for Contributors

* Always check the `bin/` directory for the final compiled binary.
* Run `./bin/punkst --help` to explore available commands.
* Ensure dependencies (**Zlib**, **BZip2**, **TBB**, **OpenCV**) are installed and properly linked in CMake.
* Do **not** commit `build/` or `bin/` — they are machine-specific and should remain local.

### 🔮 Next Steps

* Add a dedicated `BUILDING.md` with platform-specific instructions.
* (Optional) Provide a Dockerfile or Conda environment for reproducible builds.



## License

This project is licensed under the terms of the LICENSE file included in this repository.

```
