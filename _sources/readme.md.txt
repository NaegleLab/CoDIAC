## Overview
This template project sets up the basic structure with examples for Sphinx Documentation for a Python project


Clone this template into your new repository dedicated to documentation.

This is set for an example make, but run sphinx-quickstart to set the new config.py file as you want it.

You may need to install special packages for the example config.py
* conda install -c conda-forge nbsphinx
* conda install sphinxcontrib
* conda install numpydoc
* pip install myst-parser

## Instructions for generating documentation for your Project

make html

copy the entire _build/html directory to a branch of <YOUR PROJECT_PACKAGE local repo checkout> branch gh-pages

Commit and push to gh-pages branch to release to githubio
