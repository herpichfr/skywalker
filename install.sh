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
        # Debian/Ubuntu split venv's ensurepip into python3-venv; elsewhere
        # it ships with Python. Test the capability, not a package name.
        if ! python3 -c "import venv, ensurepip" >/dev/null 2>&1; then
            echo "Python3 venv is not usable. On Debian/Ubuntu install it with sudo $PACKAGE_MANAGER install python3-venv"
            exit 1
        fi
        echo "Python3 is installed."
    fi

    # pip fetches the required astroplan fork straight from GitHub.
    if ! command -v git >/dev/null 2>&1; then
        echo "git is not installed. Please install it with sudo $PACKAGE_MANAGER install git"
        exit 1
    fi
    echo "git is installed."

}

function install_dependencies() {
    echo "Installing dependencies..."
    PATH_TO_THIS_SCRIPT=$(dirname "$(readlink -f "$0")")
    # Create a virtual environment
    python3 -m venv "$PATH_TO_THIS_SCRIPT/venv"
    # Deactivate any existing virtual environment
    deactivate 2>/dev/null || true
    # Activate the virtual environment
    source "$PATH_TO_THIS_SCRIPT/venv/bin/activate"
    # Install the package and its dependencies
    if [ -f "$PATH_TO_THIS_SCRIPT/pyproject.toml" ]; then
        # [web] adds plotly and dash, so --savehtml and --web work too.
        pip install -e "${PATH_TO_THIS_SCRIPT}[web]"
    else
        echo "pyproject.toml not found. Please run this script from within the skywalker repository."
        exit 1
    fi

    echo "Dependencies installed successfully."
    echo ""
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    DEFAULT='\033[0m'
    echo -e "${RED} IMPORTANT: You now need to activate the virtual environment with the following command before running the program:${DEFAULT}"
    echo "source $PATH_TO_THIS_SCRIPT/venv/bin/activate"
    echo ""
    echo -e "${GREEN} Now you can run the script to test it (you can use an object of your choice):${DEFAULT}"
    echo "skywalker --object sirius"
}

function uninstall_dependencies() {
    PATH_TO_THIS_SCRIPT=$(dirname "$(readlink -f "$0")")
    echo "To uninstall skywalker, deactivate the virtual environment (if active) and remove it:"
    echo "deactivate && rm -rf $PATH_TO_THIS_SCRIPT/venv"
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
