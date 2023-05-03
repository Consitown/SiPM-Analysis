#!/bin/bash

# Set everything according to https://github.com/conda-forge/root-feedstock
conda config --add channels conda-forge
conda config --set channel_priority strict

# Install ROOT
conda install root root_base

# Refresh the current terminal session
source $HOME/.bashrc