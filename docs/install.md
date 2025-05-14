## Installation Guide for **punkst**

This guide walks you through building **punkst** on Linux and macOS, including environments without root access.

---

## Build

Prerequisites

- Git
- CMake: 3.15 to 3.23
- C++17 compiler* (GCC ≥8, Clang ≥5, MSVC 2017+)
- TBB, OpenCV

*We do assume your compiler properly supports C++17. Consider updating your compiler if you encounter issues.

```bash
# 1) Clone the repository
git clone --recursive https://github.com/your-org/punkst.git
cd punkst
# 2) Create and enter a build directory
mkdir build && cd build
# 3) Configure
cmake ..
# 4) Build
cmake --build . --parallel # or make
```

If you did not clone the submodule (Eigen) initially, you can do
```bash
git submodule update --init
```

If TBB is not found, you can either install it yourself (see below) or add a flag `cmake .. -DFETCH_TBB=ON` to let cmake fetch and build it from [oneTBB](https://github.com/uxlfoundation/oneTBB?tab=readme-ov-file) (which will take a while).

If you installed some dependencies locally, you may need to specify their paths like
```bash
cmake .. \
  -DOpenCV_DIR=$HOME/.local/lib/cmake/opencv4 \
  -DCMAKE_PREFIX_PATH="$HOME/.local"
  ```
(On mac, if CMake fails to locate OpenCV (installed with brew), pass: `-DOpenCV_DIR=$(brew --prefix opencv)/lib/cmake/opencv4` or wherever OpenCV is installed.)

The `punkst` binary will be placed in `bin/` under the project root.

Verifying the Build

```bash
punkst/bin/punkst --help
```

You should see a message starting with
```
Available Commands
The following commands are available:
```

### Required Libraries

- **TBB**

System: `sudo apt-get install libtbb-dev` or `yum install tbb-devel` on linux and `brew install tbb` on macOS.

Local: `git clone` from [oneTBB](https://github.com/uxlfoundation/oneTBB) then build locally.

- **OpenCV**

`sudo apt-get install libopencv-dev`or `sudo yum install opencv-devel` on linux, `brew install opencv` on macOS. See [OpenCV installation guide](https://docs.opencv.org/4.x/d7/d9f/tutorial_linux_install.html) for more details on how to install from source.


- Other dependencies

| Library      | Ubuntu / Debian                     | CentOS / RHEL         | macOS (Homebrew)            |
|--------------|--------------------------------------|-----------------------|-----------------------------|
| **zlib**     | `sudo apt-get install zlib1g-dev`   | `sudo yum install zlib-devel` | `brew install zlib`      |
| **BZip2**    | `sudo apt-get install libbz2-dev`   | `sudo yum install bzip2-devel` | `brew install bzip2`    |
| **LibLZMA**  | `sudo apt-get install liblzma-dev`  | `sudo yum install xz-devel` | `brew install xz`          |
