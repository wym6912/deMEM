## Figures Genration Tool

### Introduction

This folder is designed to generate the figures in manuscript "deMEM: a novel divide-and-conquer framework based on de Bruijn graph for scalable multiple sequence alignment" based on  `gen.py` scripts located in different folders. Each folder corresponds to a figure in the manuscript, and the corresponding `gen.py` script generates the figure.

These figures are saved in SVG format by default. If you need to convert them to PDF, tools like `Adobe Illustrator` can be used for conversion.

**Note** : It is recommended to use **Windows** for generating figures, as the font used in this manuscript is  **Times New Roman** , which may not be available on all Linux systems. Additionally, ensure that you're using **Python version 3.12.9** for the best compatibility with the drawing scripts.

## Directory Structure

```plaintext
./
│
├── Figure 4/
|   ├── Figure4.svg
│   ├── gen.py
│   └── data.csv
├── Figure 5/
|   ├── Figure5.svg
│   ├── gen.py
│   └── data.csv
├── Figure 6/
│   ├── gen.py
│   └── data.csv
├── Figure 7/
|   ├── gen.py
|   └── data.csv
├── README.md
└── requirements.txt
```

* **Figure ?/** : Each folder corresponds to a figure (e.g., `Figure 4`, `Figure 5`). Each folder contains a `gen.py` script, a data file (`data.csv`), and the graph file (`Figure?.svg`).
* **gen.py** : This script generates the figure based on the data and configuration file in the folder.
* **data.csv** : Contains the data required to generate the figure, formatted according to the `gen.py` script's data loading requirements.
* **Figure?.svg** : Contains the final graph file which used in manuscript.

## Usage

### 1. Generating Figures

This project is compatible with  **Python 3.12.9** . Make sure you are using this version of Python to avoid any compatibility issues with the drawing scripts.

You can check your Python version by running the following command:

```powershell
python --version
```

Each folder's `gen.py` script is responsible for generating a specific figure. Follow these steps:

1. **Navigate to the Figure Folder** : Go to the folder of the figure you want to generate, such as `Figure 1` or `Figure 2`.
2. **Run the Script** : In the folder, run the `gen.py` script to generate the corresponding figure.

```powershell
cd "Figure 4"
python gen.py
```

This command will generate the figure based on the `data.csv` file in the `Figure 4/` folder, and save it as an SVG file.

3. **View the Output** : The generated figure will be saved as `Figure4.svg` in the current directory or a specified output directory.

### 2. Converting to PDF

If you need to convert the generated SVG figures to PDF format, you can use tools like:

* **Adobe Illustrator** :

1. Open the SVG file in Adobe Illustrator.
2. Select `File` -> `Save As` and choose PDF format.

### 3. Dependencies

Each `gen.py` script depend on the following Python libraries:

* `matplotlib`
* `pandas`
* `numpy`

You can install the required dependencies using the following command:

```powershell
pip install -r requirements.txt
```

### 4. Generating Multiple Figures

To generate multiple figures, simply run the `gen.py` script in each folder. Each script will generate a separate figure, saving them in their respective folders.

## Contributing

If you have suggestions for improvements or feature requests, feel free to submit an Issue or Pull Request.

## Contact

If you encounter any issues, please contact us via:

* Email: [wym6912@outlook.com](mailto:wym6912@outlook.com)
* GitHub Issues: [https://github.com/malabz/deMEM/issues](https://github.com/malabz/deMEM/issues)
