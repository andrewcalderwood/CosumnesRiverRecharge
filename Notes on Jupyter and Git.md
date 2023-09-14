# Notes on Jupyter and Git

Options:
1. nbstripout: Create a git pre-commit hook that strips the output from a jupyter notebook so that it is uploaded to GitHub cleanly so it is more readily tracked and diffed.  
2. jupytext: You can create a text file or Markdown file that is paired with a jupyter notebook. These paired files can be edited interchangeably. The benefit of this method is that the text files are easy to track changes in Git and are cleaner to read for diff.  
	- method used by flopy developers for notebook examples


## jupytext
Installation:
`conda install jupytext -c conda-forge`

Here are some steps for pairing a notebook with a Python script: 
1. Open your .ipynb notebook in Jupyter.
2. Pair the notebook to a .py notebook using the pair command in Jupyter Lab (Ctrl+shift+C).
3. Save the notebook after making a change or executing a cell.
4. Add the .py notebook to version control (commit via Git).

## nbstripout
*Have not figured out use with already existing notebooks*

This package must be installed once for the environment, but must be activated for each repository individually.
Installation:
`conda install -c conda-forge nbstripout` 
Repository setup:
`nbstripout --install --attributes .gitattributes`  
Removal: 
`nbstripout --uninstall --attributes .gitattributes`

Exclude nbstripout from specific folders or notebooks
`docs/** filter= diff=`
`notebooks/Analysis.ipynb filter= diff=`