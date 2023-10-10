# Web APP (Anchor-based Algorithm)

## Description
The Flask App serves **Anchor-based Algorithm** to direct sequence RNA.

## System Requirements

### Software requirements
#### OS requirements
This project is supported for macOS/Linux. The project has been tested on macOS (Sonoma 14.0).

#### Python dependencies
1. Python version 3.7+ is required. Go to the [Python Website](https://www.python.org/downloads/) to download python and follow the instruction for installation.
2. The rest dependencies are listed inside requirements.txt, please refer to the [Prepare the Environment](#prepare-the-environment) section on how to install them. They should install within 1 minute, depending on your network speed.

## Download
Clone this repository. Click the green button "Code" and click the "Download ZIP" button. Now a zip file named "anchor_app-main.zip" is downloaded. Unzip this file. In macOS/Linux system, usually the path of the project is ~/Downloads/anchor_app-main.

## Prepare the Environment
1. Open the Terminal.
2. Create an virtual environment by entering the following command
```Bash
python3 -m venv <your_virtual_workspace_path>
```
To make it simple, we can setup the virtual environment "vir_env" in "Downloads" folder by entering the command
```Bash
python3 -m venv ~/Downloads/vir_env
```
3. Activate the virtual environment
```Bash
source <your_virtual_workspace_path>/bin/activate
```
In our case, enter the command
```Bash
source ~/Downloads/vir_env/bin/activate
```
4. Go to the root directory of anchor_app project and install all the required libraries. Enter the command
```Bash
cd <root_directory_path>
pip install -r requirements.txt
```
In our case, enter
```Bash
cd ~/Downloads/anchor_app-main
pip install -r requirements.txt
```

## Launch the Project
Run the following command to launch the project, then project **anchor_app** will be located in your default browser.
```Bash
python run.py
```

Open your default browser, and navigate to the address
http://127.0.0.1:8999/seq

## Project Structure Description
- The folder <em>**lcmsseq**</em> includes the codes of **Anchor-based Algorithm** algorithm, which is a revised version based on [**lcmsseq**](https://github.com/szostaklab/lcmsseq).

- <em>**run.py**</em> is a simple Flask server, serves at TCP port 8999.
