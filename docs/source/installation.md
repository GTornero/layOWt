# Installation

To install LayOWt it is recommended to create a separate virtual environment due to the number of package dependencies.

The following steps are recommended for a successful installation:
1. Create and activate a new virtual environment.
2. Update your `pip` version.
3. Install pre-built GDAL and Fiona wheel files with `pip`.
4. Download the latest LayOWt version from the project [GitHub](https://github.com/GTornero/layOWt).
5. Install the downloaded LayOWt package with `pip`.

## 1. Creating a Virtual Environment

First step is to create a new virtual environment using `venv` for example.

If `venv` is not installed, you can install it using `pip` with the following command:

```console
pip install venv
```

Once `venv` is installed on your machine, you can create a new environment.
Open the `cmd` at the desired destination folder of your environment.
Note that if you have mutiple versions of python installed, you can choose which version to use when creating this environment.
You can see all installed versions of python using the following command:

```console
py --list

Installed Pythons found by C:\WINDOWS\py.exe Launcher for Windows
 -3.10-64 *
 -3.8-64
```

It is recommended to use python versions 3.10 and above for LayOWt.

You can create a virtual environment using the following command:

```console
python -m venv [path/to/venv] --prompt "LayOWt-venv"
```

If a specific version of python is to be used (for example 3.10), use:

```console
py -3.10 -m venv [path/to/venv] --prompt "LayOWt-venv"
```

To activate this virtual environment, use the following command if using Windows PowerShell:

```console
[path/to/venv]/Scripts/Activate.ps1
```

If using the cmd Shell, use the following command:

```console
[path/to/venv]/Scripts/activate.bat
```

You should now see the prompt you set when creating the Virtual environment as a prefix to your Shell session.

## 2. Update `pip` in the Virtual Environment

Run the following command to update `pip`:

```console
(LayOWt-venv) python -m pip install -U pip
```

## 3. Installing GDAL, Fiona, and Rasterio Wheels

Installing the `GDAL`, `Fiona`, and `Rasterio` packages for Python using `pip` is not as easy as it may seem. These packages are actually C++ libraries with Python bindings. Therefore, in order to be installed correctly, the underlying C++ code must be  built/compiled prior to be used with Python.

Luckily, the binaries for `GDAL`, `Fiona`, and `Rasterio` can be found [here](https://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal).
Select the correct versions for your Python version and download the .whl files.

Open the terminal using the LayOWt Virtual Environment and install the binaries using `pip`.
First install `GDAL`:

```console
(LayOWt-venv) pip install [path/to/whl]/GDAL-3.4.3-cp310-cp310-win_amd64.whl
```
Then install `Fiona`:

```console
(LayOWt-venv) pip install [path/to/whl]/Fiona-1.8.21-cp310-cp310-win_amd64.whl
```

Finally install `Rasterio`:

```console
(LayOWt-venv) pip install [path/to/whl]/rasterio-1.2.10-cp310-cp310-win_amd64
```

## 4. Download LayOWt from [GitHub](https://github.com/GTornero/layOWt)

Make sure you download/clone the latest stable version from the main branch.

## 5. Install LayOWt Using `pip`

The final step is to install the LayOWt package using `pip`.
Open the terminal, activate the virtual environment, and change directory to where you have extracted the LayOWt package you have downloaded/cloned from GitHub.

Run the following command and voila!

```console
(LayOWt-venv) pip install .
```

As an optional, you can also install all of the dependencies LayOWt uses to build its documentation using the following command instead of the previous:

```console
(LayOWt-venv) pip install .[doc]
```

```{note}
If you wish to install the LayOWt package in editable mode, you can run:

```console
(LayOWt-venv) pip install -e .[doc]
```