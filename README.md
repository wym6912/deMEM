# deMEM: a novel divide-and-conquer framework based on de Bruijn graph for scalable multiple sequence alignment

deMEM is a program written in C++ for aligning MSA based on de Bruijn graph with divide-and-conquer framework. It runs on Linux and Windows.

## Usage

```
        Usage: samemalign_main (64-bit) [options]

Use deBruijn graph find MEM, and make alignment by a specified aligner like WMSA/MAFFT/HAlign 3
Available options:
File process:
        --in         FILE      input sequence file name (Required)
        --out        FILE      output file name (Required)
        --removechr            remove characters in sequence
Config JSON file: (Recommended)
        --json       FILE      JSON config, contain "Conditions", "MEM arugments" and "Alignment arugments" configurations
        --jsonsample FILE      output JSON config sample to FILE and exit
        --keepjson             keep JSON file when finish
Conditions:
        --thread     N         use N threads to run program (default: 1)
        --type       X         sequence type: X = 1 is Protein, X = 2 is DNA (Required)
MEM arugments:
        --threshold  t         the minimal length of MEM (default = 10)
        --memregen             re-generate MEM (default: false)
Alignment arugments:
        --highsim              this file is high similarity, will treat all sequences into one cluster; will ignore --cluster,
                               --wmsac and --mafftc
        --path       PATH      temp file will be in path, default = "swap/" (only support first-level directory)
        --threshold2 t         the minimal length of blockaligner (default = 10)
                               not recommended for t < 10
        --align      A         alignment method: A = 0 is WMSA, A = 1 is MAFFT, A = 2 is HAlign, A = 3 is abPOA (default: 0)
        --wmsaa      w         WMSA align method: w = 0 use FFT+K-band align, w = 1 only use K-band align (default: 0)
        --maffta     m         MAFFT align method: m = 0 is NW-NS-1, m = 1 is NW-NS-2, m = 2 is NW-NS-i, m = 3 is FFT-NS-1,
                                                   m = 4 is FFT-NS-2, m = 5 is FFT-NS-i, m = 6 is G-INS-i, m = 7 is L-INS-i
                                                   (default: 3)
        --gapopena   a         use user-defined gap open/gap penalty in align for WMSA, MAFFT and abPOA
        --gapexta    a         use user-defined gap extension in align for MAFFT and abPOA
        --merge      M         merge two profiles method: M = 0 is WMSA, M = 1 is MAFFT, M = 2 is abPOA (default: 0)
        --wmsam      w         WMSA merge method: w = 0 use FFT+K-band with no gap open calculation,
                                                  w = 1 use K-band with no gap open calculation,
                                                  w = 2 use FFT+K-band with gap open calculation,
                                                  w = 3 use K-band with gap open calculation (default: 0)
        --mafftm     m         MAFFT two profiles align method: m = 0 is NW-NS, m = 1 is FFT-NS (default: 1)
        --gapopenm   m         use user-defined gap open/gap penalty in merge two profiles for WMSA, MAFFT and abPOA
        --gapextm    m         use user-defined gap extension in merge two profiles in WMSA, MAFFT and abPOA
        --cluster    C         merge clusters method: C = 0 is WMSA, C = 1 is MAFFT, C = 2 is abPOA (default: 0)
        --wmsac      w         WMSA cluster merge method: w = 0 use FFT+K-band with no gap open calculation,
                                                          w = 1 use K-band with no gap open calculation,
                                                          w = 2 use FFT+K-band with gap open calculation,
                                                          w = 3 use K-band with gap open calculation (default: 1)
        --mafftc     m         MAFFT cluster merge method: m = 0 is NW-NS-1, m = 1 is NW-NS-2, m = 2 is NW-NS-i, m = 3 is FFT-NS-1,
                                                           m = 4 is FFT-NS-2, m = 5 is FFT-NS-i, m = 6 is G-INS-i, m = 7 is L-INS-i
                                                           (default: 3)
        --abpoac     p         abPOA cluster merge method: p = 0 use center sequences make pre-align and insert clusters,
                                                           p = 1 directly align clusters (default: 0)
        --gapopenc   c         use user-defined gap open/gap penalty in merge clusters for WMSA, MAFFT and abPOA
        --gapextc    c         use user-defined gap extension in merge clusters for WMSA, MAFFT and abPOA
        --allmafft             all method use MAFFT
        --allwmsa              all method use WMSA
        --allabpoa             all method use abPOA
Others:
        --noverbose            not print the info for subprograms
        --help                 print help message
        --version              show program version
Hint: If read json file and use same arugments, program will omit the arugment settings.
Example:
        samemalign_main --in seq.fasta --out seq_out.fasta
```

We strongly recommended run the program by determined json files, which in `config` folder.

## Installation

### Linux/WSL (Windows Subsystem for Linux) - from Anaconda (Recommended)

1.Install WSL for Windows. Instructional video [1](https://www.youtube.com/watch?v=X-DHaQLrBi8&t=5s) or [2](http://lab.malab.cn/%7Etfr/1.mp4) (Copyright belongs to the original work).

2.Download and install Anaconda. Download Anaconda for different systems [here](https://www.anaconda.com/products/distribution#Downloads). Instructional video of anaconda installation [1](https://www.youtube.com/watch?v=AshsPB3KT-E) or [2](http://lab.malab.cn/%7Etfr/Install_anaconda_in_Linux.mp4) (Copyright belongs to the original work).

3.Install deMEM.

```bash
#1 Create and activate a conda environment for deMEM
conda create -n demem
conda activate demem

#2 Add channels to conda
conda config --add channels wym6912
conda config --add channels malab
conda config --add channels conda-forge

#3 Install deMEM
conda install -c malab -c conda-forge -c wym6912 demem

#4 Test deMEM
samemalign_main --help
```

### Linux/WSL (Windows Subsystem for Linux) - from the source code

1. Download and Compile the source code. (Make sure your version of gcc >= 11 or clang >= 13.0.0)

```bash
#1 Download
git clone https://github.com/malabz/demem.git --recursive
# Make sure submodules are correctly cloned
git submodule init
git submodule update
git submodule foreach git submodule init
git submodule foreach git submodule update

#2 Open the folder
cd demem/src

#3 Compile and install to path
make THREADS=16 PREFIX=/path/to/you/wish all

#4 Test
./samemalign_main --help
```

### Windows - from from Anaconda (Recommended)

1.Download and install Anaconda. Download Anaconda for different systems [here](https://www.anaconda.com/products/distribution#Downloads). Instructional video of anaconda installation [1](https://www.youtube.com/watch?v=AshsPB3KT-E) or [2](http://lab.malab.cn/%7Etfr/Install_anaconda_in_Linux.mp4) (Copyright belongs to the original work).

2.Install deMEM.

```bash
#1 Create and activate a conda environment for deMEM
conda create -n demem
conda activate demem

#2 Add channels to conda
conda config --add channels wym6912
conda config --add channels malab
conda config --add channels conda-forge

#3 Install deMEM
conda install -c malab -c conda-forge -c wym6912 demem

#4 Test deMEM
samemalign_main -h
```

### Windows - from Visual Studio 2022 (Not Recommended)

1. First, install git, and download the whole project with submodules:

```powershell
git clone https://github.com/malabz/demem.git --recursive
# Make sure submodules are correctly cloned
git submodule init
git submodule update
git submodule foreach git submodule init
git submodule foreach git submodule update
```

2. Make programs from `Visual Studio 2022`:

- First of all, download and install [`Visual Studio 2022`](https://visualstudio.microsoft.com/vs/).
- Open `memalign.sln` in `Visual Studio 2022`, then switch it to `Release` mode (Select `Properties`, choose `Configuration Properties`, then press `Configuration Manager...` button, change `Active Solution configuration` to `Release`), choose `Build` -> `Build Solution`. If failed, please try compile project one by one.
- Open `x64\Release` folder, copy all the `.exe` and `.dll` files to `x64\samemalign_env\bin`.
- Copy `src\external\submafft\windows\{profiles_wrapper.exe,src}` file/folder to `x64\samemalign_env`.
