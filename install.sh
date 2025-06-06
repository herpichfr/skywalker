#!/bin/bash

# input parameters: --install --check --uninstall

function PrintUsage() {
    echo "Usage: $0 [--install | --check | --uninstall]"
    echo "  --install   Install dependencies"
    echo "  --check     Check basic requirements"
    echo "  --uninstall Uninstall dependencies"
}
if [[ $# -eq 0 ]]; then
    PrintUsage
    exit 1
fi

function check_basic_requirements() {
    echo "Checking basic requirements..."
    # detect Linux distribution and determine package manager
    if command -v apt-get >/dev/null 2>&1; then
        echo "Using apt-get package manager."
        PACKAGE_MANAGER="apt-get"
    elif command -v yum >/dev/null 2>&1; then
        echo "Using yum package manager."
        PACKAGE_MANAGER="yum"
    elif command -v dnf >/dev/null 2>&1; then
        echo "Using dnf package manager."
        PACKAGE_MANAGER="dnf"
    else
        echo "No supported package manager found. Please install dependencies manually."
        exit 1
    fi

    # Check if python3 is installed
    if ! command -v python3 >/dev/null 2>&1; then
        echo "Python3 is not installed. Please install it with sudo $PACKAGE_MANAGER install python3 python3-pip python3-venv"
        exit 1
    elif ! command -v pip3 >/dev/null 2>&1 ; then
        echo "pip3 is not installed. Please install it with sudo $PACKAGE_MANAGER install python3-pip"
        exit 1
    else
        if command -v dpkg >/dev/null 2>&1; then
            if ! dpkg -s python3-venv >/dev/null 2>&1; then
                echo "Python3 venv is not installed. Please install it with sudo $PACKAGE_MANAGER install python3-venv"
                exit 1
            fi
        elif command -v rpm >/dev/null 2>&1; then
            if ! rpm -q python3-venv >/dev/null 2>&1; then
                echo "Python3 venv is not installed. Please install it with sudo $PACKAGE_MANAGER install python3-venv"
                exit 1
            else
                echo "Python3 venv is installed."
            fi
        fi
        echo "Python3 is installed."
    fi

}

function install_dependencies() {
    echo "Installing dependencies..."
    # Create a virtual environment
    python3 -m venv venv
    # Activate the virtual environment
    source venv/bin/activate
    # Install required packages
    if [ -f requirements.txt ]; then
        pip install -r requirements.txt
    else
        echo "requirements.txt not found. Please provide a valid requirements file."
        exit 1
    fi
    # install astroplan from github
    CURRENT_PATH=$(pwd)
    git clone git@github.com:herpichfr/astroplan.git "$CURRENT_PATH/astroplan"
    cd "$CURRENT_PATH/astroplan" || exit 1
    python3 setup.py install
    cd "$CURRENT_PATH" || exit 1

    echo "Dependencies installed successfully."
}

function uninstall_dependencies() {
    echo "Uninstalling dependencies..."
    # Deactivate the virtual environment if it is active
    if [[ "$VIRTUAL_ENV" != "" ]]; then
        deactivate
    fi
    # find path to this script
    SCRIPT_PATH=$(realpath "$0")

    # Remove the virtual environment directory
    if [ -d "$SCRIPT_PATH/venv" ]; then
        rm -rf "$SCRIPT_PATH/venv"
        echo "Virtual environment removed."
    else
        echo "No virtual environment found to remove."
    fi
}

case "$1" in
    --install)
        check_basic_requirements
        install_dependencies
        ;;
    --check)
        check_basic_requirements
        ;;
    --uninstall)
        uninstall_dependencies
        ;;
    *)
        PrintUsage
        exit 1
        ;;
esac
