#!/usr/bin/env python
"""
Headless (no-GUI) wrapper for rdeditor in Jupyter Notebooks

This module provides a lightweight interface for programmatic molecule
manipulation without requiring Qt GUI. Use this when you don't need
the interactive editor window.

Example usage:
    ```python
    from rdeditor.jupyter_headless import HeadlessMoleculeEditor
    
    editor = HeadlessMoleculeEditor("CCO")
    # Edit programmatically
    editor.smiles = "c1ccccc1"
    mol = editor.mol
    ```
"""

from typing import Optional, Union
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDepictor


class HeadlessMoleculeEditor:
    """
    Lightweight molecule editor without Qt GUI dependencies.
    
    This class provides basic molecule manipulation without requiring
    the Qt event loop or GUI components. Perfect for Jupyter notebooks
    in headless environments or when you only need programmatic access.
    
    Attributes:
        mol: The current RDKit Mol object
        smiles: SMILES representation of the molecule
    """
    
    def __init__(self, molecule: Union[str, Chem.Mol, None] = None):
        """
        Initialize the headless molecule editor.
        
        Args:
            molecule: SMILES string, RDKit Mol object, or None
        """
        self._mol = None
        
        if molecule is not None:
            if isinstance(molecule, str):
                self.smiles = molecule
            elif isinstance(molecule, Chem.Mol):
                self.mol = molecule
            else:
                raise TypeError("molecule must be a SMILES string or RDKit Mol object")
    
    @property
    def mol(self) -> Optional[Chem.Mol]:
        """Get the current molecule."""
        return self._mol
    
    @mol.setter
    def mol(self, value: Chem.Mol):
        """Set the current molecule."""
        if value is not None:
            # Ensure molecule has 2D coordinates
            if value.GetNumConformers() == 0:
                rdDepictor.Compute2DCoords(value)
        self._mol = value
    
    @property
    def smiles(self) -> Optional[str]:
        """Get the SMILES string of the current molecule."""
        if self._mol is None:
            return None
        return Chem.MolToSmiles(self._mol, isomericSmiles=True)
    
    @smiles.setter
    def smiles(self, value: str):
        """Set the molecule from a SMILES string."""
        mol = Chem.MolFromSmiles(value)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {value}")
        self.mol = mol
    
    def load(self, filename: Union[str, Path]):
        """
        Load a molecule from a file.
        
        Args:
            filename: Path to MOL or SMI file
        """
        filename = str(filename)
        
        if filename.lower().endswith('.mol'):
            mol = Chem.MolFromMolFile(filename, sanitize=False, strictParsing=False)
        elif filename.lower().endswith('.smi'):
            with open(filename, 'r') as f:
                smiles = f.readline().strip()
            mol = Chem.MolFromSmiles(smiles)
        else:
            raise ValueError(f"Unsupported file format: {filename}")
        
        if mol is None:
            raise ValueError(f"Failed to load molecule from: {filename}")
        
        self.mol = mol
    
    def save(self, filename: Union[str, Path]):
        """
        Save the current molecule to a file.
        
        Args:
            filename: Path to save the molecule (MOL or SMI format)
        """
        if self._mol is None:
            raise ValueError("No molecule to save")
        
        filename = str(filename)
        
        if filename.lower().endswith('.mol'):
            Chem.MolToMolFile(self._mol, filename)
        elif filename.lower().endswith('.smi'):
            smiles = Chem.MolToSmiles(self._mol, isomericSmiles=True)
            with open(filename, 'w') as f:
                f.write(smiles + '\n')
        else:
            raise ValueError(f"Unsupported file format: {filename}")
    
    def clear(self):
        """Clear the current molecule."""
        self._mol = None
    
    def copy(self) -> 'HeadlessMoleculeEditor':
        """
        Create a copy of this editor with the same molecule.
        
        Returns:
            New HeadlessMoleculeEditor instance
        """
        new_editor = HeadlessMoleculeEditor()
        if self._mol is not None:
            new_editor.mol = Chem.Mol(self._mol)
        return new_editor
    
    def __repr__(self):
        """String representation of the editor."""
        if self._mol is not None:
            return f"HeadlessMoleculeEditor(SMILES='{self.smiles}')"
        return "HeadlessMoleculeEditor(empty)"
    
    def __enter__(self):
        """Context manager entry."""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.clear()
        return False


# Convenience functions
def edit_smiles(smiles: str) -> str:
    """
    Quick function to validate and canonicalize a SMILES string.
    
    Args:
        smiles: Input SMILES string
        
    Returns:
        Canonical SMILES string
        
    Example:
        ```python
        canonical = edit_smiles("c1ccccc1")
        print(canonical)  # c1ccccc1
        ```
    """
    editor = HeadlessMoleculeEditor(smiles)
    return editor.smiles


def convert_mol_to_smiles(mol: Chem.Mol, isomeric: bool = True) -> str:
    """
    Convert RDKit Mol to SMILES.
    
    Args:
        mol: RDKit Mol object
        isomeric: Include stereochemistry information
        
    Returns:
        SMILES string
    """
    return Chem.MolToSmiles(mol, isomericSmiles=isomeric)


def convert_smiles_to_mol(smiles: str) -> Chem.Mol:
    """
    Convert SMILES to RDKit Mol with 2D coordinates.
    
    Args:
        smiles: SMILES string
        
    Returns:
        RDKit Mol object with 2D coordinates
    """
    editor = HeadlessMoleculeEditor(smiles)
    return editor.mol


# Export public API
__all__ = [
    'HeadlessMoleculeEditor',
    'edit_smiles',
    'convert_mol_to_smiles',
    'convert_smiles_to_mol',
]


if __name__ == "__main__":
    # Example usage
    print("HeadlessMoleculeEditor Example")
    print("-" * 50)
    
    # Example 1: Basic usage
    editor = HeadlessMoleculeEditor("CCO")
    print(f"Molecule: {editor.smiles}")
    print(f"Atoms: {editor.mol.GetNumAtoms()}")
    
    # Example 2: Conversion
    editor.smiles = "c1ccccc1"
    print(f"Benzene: {editor.smiles}")
    
    # Example 3: Context manager
    with HeadlessMoleculeEditor("CC(=O)O") as ed:
        print(f"Acetic acid: {ed.smiles}")
