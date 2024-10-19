#!/usr/bin/env python3 

import subprocess
import sys, os 
import importlib.util
pwd = os.path.abspath(os.path.dirname(__file__))

def install_package(package):
    """Install a package using pip."""
    try:
        print(f"Attempting to install {package}...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        print(f"{package} installed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while installing {package}: {e}")
        
def is_package_installed(package):
    """Check if a package is already installed."""
    spec = importlib.util.find_spec(package)
    return spec is not None

def install_requirements(requirements_file=os.path.join(pwd, "card/requirements.txt")):
    """Install all dependencies listed in requirements.txt, if not already installed."""
    try:
        with open(requirements_file, "r") as f:
            packages = f.readlines()
            packages = [pkg.strip() for pkg in packages]
            for package in packages:
                package_name = package.split("==")[0]  # Handle version constraints if present
                if is_package_installed(package_name):
                    print(f"{package_name} is already installed, skipping.")
                else:
                    install_package(package)
        print("All dependencies are installed.")
    except FileNotFoundError:
        print(f"{requirements_file} not found. Please make sure it exists.")

if __name__ == "__main__":
    install_requirements()
