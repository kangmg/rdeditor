# Using rdeditor in Jupyter Notebooks

This guide explains how to use rdeditor in Jupyter Notebooks with the new wrapper interface.

## Installation

```bash
pip install rdeditor
```

## Quick Start

### 1. Simple Editing

```python
from rdeditor import MoleculeEditor

# Create editor with a SMILES string
editor = MoleculeEditor("CCO")

# Open the GUI editor
editor.edit()

# Access the edited molecule
print(f"Edited SMILES: {editor.smiles}")
```

### 2. Display Molecules

```python
from rdeditor import display_molecule
from rdkit import Chem

mol = Chem.MolFromSmiles("c1ccccc1")
display_molecule(mol)
```

## API Overview

### RDEditorNotebook

The main class providing full control over the editor.

```python
from rdeditor import RDEditorNotebook

# Create editor instance
editor = RDEditorNotebook(loglevel="INFO")

# Set molecule from SMILES
editor.set_smiles("CCO")

# Or from RDKit Mol object
from rdkit import Chem
mol = Chem.MolFromSmiles("c1ccccc1")
editor.set_molecule(mol)

# Show the editor GUI
editor.show()

# Get the edited molecule
edited_mol = editor.get_molecule()
edited_smiles = editor.get_smiles()

# Save to file
editor.save_file("molecule.mol")

# Load from file
editor.load_file("molecule.mol")

# Clear the molecule
editor.clear()

# Hide/close the editor
editor.hide()
editor.close()
```

### MoleculeEditor

Simplified interface with property-based access.

```python
from rdeditor import MoleculeEditor

# Create with SMILES or Mol object
editor = MoleculeEditor("c1ccccc1")

# Access properties
print(editor.smiles)  # Get SMILES
editor.smiles = "CCO"  # Set SMILES

print(editor.mol)     # Get RDKit Mol
editor.mol = mol      # Set RDKit Mol

# Open GUI
editor.edit()

# File operations
editor.save("molecule.mol")
editor.load("molecule.mol")

# Clear molecule
editor.clear()
```

### Display Functions

Rich molecule visualization in notebooks.

```python
from rdeditor import (
    display_molecule,
    display_smiles,
    display_molecules,
    compare_molecules,
)

# Display single molecule
from rdkit import Chem
mol = Chem.MolFromSmiles("CCO")
display_molecule(mol, size=(300, 300))

# Display from SMILES
display_smiles("c1ccccc1")

# Display multiple molecules in a grid
mols = [Chem.MolFromSmiles(s) for s in ["CCO", "c1ccccc1", "CC(=O)O"]]
labels = ["Ethanol", "Benzene", "Acetic Acid"]
display_molecules(mols, mols_per_row=3, labels=labels)

# Compare two molecules
mol1 = Chem.MolFromSmiles("CCO")
mol2 = Chem.MolFromSmiles("CCCO")
compare_molecules(mol1, mol2, label1="Ethanol", label2="Propanol")
```

## Common Workflows

### 1. Interactive Molecule Editing

```python
from rdeditor import MoleculeEditor, display_molecule

# Start with a molecule
editor = MoleculeEditor("c1ccccc1C(=O)O")

# Display original
print("Original:")
display_molecule(editor.mol)

# Edit it
editor.edit()

# Display edited
print("Edited:")
display_molecule(editor.mol)
```

### 2. Batch Processing

```python
from rdeditor import MoleculeEditor, display_molecules
from rdkit import Chem

molecules = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
}

edited = {}

for name, smiles in molecules.items():
    print(f"Editing {name}...")
    editor = MoleculeEditor(smiles)
    # editor.edit()  # Uncomment to edit interactively
    edited[name] = editor.smiles

# Display all
mols = [Chem.MolFromSmiles(s) for s in edited.values()]
display_molecules(mols, labels=list(edited.keys()))
```

### 3. Integration with RDKit

```python
from rdeditor import MoleculeEditor
from rdkit import Chem
from rdkit.Chem import Descriptors

# Create and edit molecule
editor = MoleculeEditor("CCO")
editor.edit()

# Calculate properties
mol = editor.mol
if mol:
    print(f"Molecular Weight: {Descriptors.MolWt(mol):.2f}")
    print(f"LogP: {Descriptors.MolLogP(mol):.2f}")
    print(f"H-Bond Donors: {Descriptors.NumHDonors(mol)}")
    print(f"H-Bond Acceptors: {Descriptors.NumHAcceptors(mol)}")
```

### 4. Context Manager

```python
from rdeditor import RDEditorNotebook

# Automatic cleanup
with RDEditorNotebook() as editor:
    editor.set_smiles("c1ccccc1")
    editor.show()
    result = editor.get_smiles()
    
# Editor is automatically closed
print(f"Result: {result}")
```

### 5. File I/O

```python
from rdeditor import MoleculeEditor

# Edit and save
editor = MoleculeEditor("c1ccccc1")
editor.edit()
editor.save("benzene.mol")
editor.save("benzene.smi")

# Load and edit
editor2 = MoleculeEditor()
editor2.load("benzene.mol")
editor2.edit()
```

## Advanced Features

### Custom Display

```python
from rdeditor import MoleculeDisplay
from rdkit import Chem

mol = Chem.MolFromSmiles("CCO")

# Create custom display
display = MoleculeDisplay(
    mol=mol,
    size=(400, 400),
    show_edit_button=True
)

# Show in notebook
display
```

### Comparison Workflow

```python
from rdeditor import MoleculeEditor, compare_molecules

# Original
original = "c1ccccc1C(=O)O"
editor = MoleculeEditor(original)

# Edit
editor.edit()

# Compare
from rdkit import Chem
mol1 = Chem.MolFromSmiles(original)
mol2 = editor.mol

compare_molecules(
    mol1, mol2,
    label1=f"Original: {original}",
    label2=f"Edited: {editor.smiles}"
)
```

## Troubleshooting

### Qt Event Loop Issues

If you encounter event loop issues, try:

```python
# Option 1: Use %gui qt magic command (IPython/Jupyter)
%gui qt

# Option 2: Process events manually
from rdeditor import RDEditorNotebook

editor = RDEditorNotebook()
editor.show()
editor._process_events()
```

### Display Issues

If molecules don't display properly:

```python
# Ensure RDKit display is set up
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG = True
```

### Import Errors

If you get import errors:

```bash
# Install missing dependencies
pip install rdkit pyside6 numpy pyqtdarktheme

# For Jupyter support
pip install ipython jupyter
```

## Tips and Best Practices

1. **Close editors when done**: Always close or hide editors to free resources
   ```python
   editor.close()
   ```

2. **Use context managers**: For automatic cleanup
   ```python
   with RDEditorNotebook() as editor:
       # Your code here
       pass
   ```

3. **Check molecule validity**: Always verify molecules after editing
   ```python
   mol = editor.get_molecule()
   if mol is not None:
       # Process molecule
       pass
   ```

4. **Save frequently**: Save your work to avoid losing changes
   ```python
   editor.save("backup.mol")
   ```

5. **Display before editing**: View molecules before editing
   ```python
   from rdeditor import display_molecule
   display_molecule(mol)
   editor.edit()
   ```

## Examples

See the complete tutorial notebook:
- [examples/rdeditor_jupyter_tutorial.ipynb](examples/rdeditor_jupyter_tutorial.ipynb)

## API Reference

### RDEditorNotebook

- `__init__(loglevel="WARNING")`: Initialize editor
- `show(smiles=None, mol=None, mol_file=None)`: Show editor window
- `hide()`: Hide editor window
- `close()`: Close editor and cleanup
- `set_molecule(mol)`: Set RDKit Mol object
- `get_molecule()`: Get RDKit Mol object
- `set_smiles(smiles)`: Set molecule from SMILES
- `get_smiles(isomeric=True)`: Get SMILES string
- `load_file(filename)`: Load from file
- `save_file(filename)`: Save to file
- `clear()`: Clear molecule

### MoleculeEditor

- `__init__(molecule=None)`: Initialize with SMILES or Mol
- `mol`: Property to get/set RDKit Mol object
- `smiles`: Property to get/set SMILES string
- `edit()`: Open GUI editor
- `close()`: Close editor
- `clear()`: Clear molecule
- `save(filename)`: Save to file
- `load(filename)`: Load from file

### Display Functions

- `display_molecule(mol, size=(300,300), show_edit_button=True)`
- `display_smiles(smiles, **kwargs)`
- `display_molecules(mols, mols_per_row=3, size=(200,200), labels=None)`
- `compare_molecules(mol1, mol2, label1="Molecule 1", label2="Molecule 2", size=(300,300))`

## Support

For issues and questions:
- GitHub Issues: https://github.com/EBjerrum/rdeditor/issues
- Documentation: https://github.com/EBjerrum/rdeditor

## License

rdeditor is licensed under the LGPL license.
