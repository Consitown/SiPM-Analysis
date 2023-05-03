#!/bin/bash

# Define the version of Miniconda to install
MINICONDA_VERSION="latest"

# Determine the platform (linux or macOS)
PLATFORM=$(uname)

# Determine the processor architecture
if [[ $(uname -m) == "x86_64" ]]; then
    ARCH="x86_64"
else
    ARCH="x86"
fi

# Download the Miniconda installation script
wget https://repo.anaconda.com/miniconda/Miniconda3-$MINICONDA_VERSION-$PLATFORM-$ARCH.sh -O miniconda.sh

# Run the installation script
bash miniconda.sh -b -p $HOME/miniconda

# Add the Miniconda bin directory to the PATH environment variable
echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> $HOME/.bashrc

# Refresh the current terminal session
source $HOME/.bashrc

# Remove the installation script
rm miniconda.sh

# Initialize the base environment
$HOME/miniconda/bin/conda init bash
$HOME/miniconda/bin/conda init zsh

# Please restart shell after this step