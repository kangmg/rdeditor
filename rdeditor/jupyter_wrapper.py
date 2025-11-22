#!/usr/bin/env python
"""
Jupyter Notebook wrapper for rdeditor

This module provides a complete wrapper for using rdeditor in Jupyter Notebooks.
It handles the Qt event loop integration and provides both GUI and programmatic interfaces.

Example usage:
    ```python
    from rdeditor.jupyter_wrapper import RDEditorNotebook
    from rdkit import Chem
    
    # Create editor
    editor = RDEditorNotebook()
    
    # Open the GUI
    editor.show()
    
    # Or use programmatically
    mol = Chem.MolFromSmiles("CCO")
    editor.set_molecule(mol)
    edited_mol = editor.get_molecule()
    ```
"""

import sys
import logging
from typing import Optional, Union
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor
from PySide6.QtWidgets import QApplication
from PySide6.QtCore import QTimer, Qt

# Import rdeditor components
from .rdEditor import MainWindow
from .molEditWidget import MolEditWidget


class RDEditorNotebook:
    """
    A wrapper class for using rdeditor in Jupyter Notebooks.
    
    This class manages the Qt application lifecycle and provides a clean
    interface for molecular editing in notebook environments.
    
    Attributes:
        window (MainWindow): The main editor window
        editor (MolEditWidget): The molecule editor widget
    """
    
    _app = None  # Class variable to store the QApplication instance
    
    def __init__(self, loglevel: str = "WARNING"):
        """
        Initialize the rdeditor for Jupyter Notebook.
        
        Args:
            loglevel: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        """
        self.loglevel = loglevel
        self._ensure_qapp()
        self.window = None
        self._is_visible = False
        
    @classmethod
    def _ensure_qapp(cls):
        """
        Ensure a QApplication instance exists.
        
        This method creates a QApplication if one doesn't exist,
        or reuses the existing instance.
        """
        if cls._app is None:
            # Check if there's already a QApplication instance
            existing_app = QApplication.instance()
            if existing_app is None:
                # Set offscreen platform for headless environments
                import os
                if 'QT_QPA_PLATFORM' not in os.environ:
                    # Only set if not already specified
                    pass
                try:
                    cls._app = QApplication(sys.argv)
                except RuntimeError as e:
                    # If QApplication creation fails, try with minimal args
                    cls._app = QApplication([])
            else:
                cls._app = existing_app
    
    def show(self, smiles: Optional[str] = None, 
             mol: Optional[Chem.Mol] = None,
             mol_file: Optional[Union[str, Path]] = None):
        """
        Show the editor window.
        
        Args:
            smiles: Optional SMILES string to load
            mol: Optional RDKit Mol object to load
            mol_file: Optional path to MOL file to load
        """
        try:
            if self.window is None:
                self.window = MainWindow(loglevel=self.loglevel)
                self.editor = self.window.editor
            
            # Load molecule if provided
            if smiles is not None:
                self.set_smiles(smiles)
            elif mol is not None:
                self.set_molecule(mol)
            elif mol_file is not None:
                self.load_file(mol_file)
            
            self.window.show()
            self._is_visible = True
            
            # Process events to ensure the window is displayed
            self._process_events()
        except Exception as e:
            print(f"Warning: Could not show GUI window: {e}")
            print("Tip: For headless environments, use methods without .show()")
            # Keep the editor accessible even if GUI fails
            if self.window is None:
                self.window = MainWindow(loglevel=self.loglevel)
                self.editor = self.window.editor
        
    def hide(self):
        """Hide the editor window without closing it."""
        if self.window is not None:
            self.window.hide()
            self._is_visible = False
            self._process_events()
    
    def close(self):
        """Close the editor window and clean up resources."""
        if self.window is not None:
            self.window.close()
            self.window = None
            self._is_visible = False
            self._process_events()
    
    def is_visible(self) -> bool:
        """Check if the editor window is currently visible."""
        return self._is_visible and self.window is not None
    
    def set_molecule(self, mol: Chem.Mol):
        """
        Set the molecule to edit.
        
        Args:
            mol: RDKit Mol object
        """
        if self.window is None:
            self.window = MainWindow(loglevel=self.loglevel)
            self.editor = self.window.editor
        
        # Ensure the molecule has 2D coordinates
        if mol.GetNumConformers() == 0:
            rdDepictor.Compute2DCoords(mol)
        
        self.editor.mol = mol
        self._process_events()
    
    def get_molecule(self) -> Optional[Chem.Mol]:
        """
        Get the current molecule from the editor.
        
        Returns:
            RDKit Mol object or None if no molecule is present
        """
        if self.window is None or self.editor is None:
            return None
        return self.editor.mol
    
    def set_smiles(self, smiles: str):
        """
        Set the molecule from a SMILES string.
        
        Args:
            smiles: SMILES string
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        self.set_molecule(mol)
    
    def get_smiles(self, isomeric: bool = True) -> Optional[str]:
        """
        Get the SMILES string of the current molecule.
        
        Args:
            isomeric: Include stereochemistry information
            
        Returns:
            SMILES string or None if no molecule is present
        """
        mol = self.get_molecule()
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, isomericSmiles=isomeric)
    
    def load_file(self, filename: Union[str, Path]):
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
        
        self.set_molecule(mol)
    
    def save_file(self, filename: Union[str, Path]):
        """
        Save the current molecule to a file.
        
        Args:
            filename: Path to save the molecule (MOL or SMI format)
        """
        mol = self.get_molecule()
        if mol is None:
            raise ValueError("No molecule to save")
        
        filename = str(filename)
        
        if filename.lower().endswith('.mol'):
            Chem.MolToMolFile(mol, filename)
        elif filename.lower().endswith('.smi'):
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            with open(filename, 'w') as f:
                f.write(smiles + '\n')
        else:
            raise ValueError(f"Unsupported file format: {filename}")
    
    def clear(self):
        """Clear the current molecule from the editor."""
        if self.window is not None and self.editor is not None:
            self.editor.mol = None
            self._process_events()
    
    def _process_events(self, timeout_ms: int = 100):
        """
        Process Qt events to update the GUI.
        
        Args:
            timeout_ms: Maximum time to process events in milliseconds
        """
        if self._app is not None:
            # Process events for a short time to update the GUI
            QTimer.singleShot(0, lambda: None)
            self._app.processEvents()
    
    def __repr__(self):
        """String representation of the editor."""
        status = "visible" if self.is_visible() else "hidden"
        mol_info = "no molecule"
        if self.window is not None and self.editor is not None:
            mol = self.get_molecule()
            if mol is not None:
                mol_info = f"{mol.GetNumAtoms()} atoms, {mol.GetNumBonds()} bonds"
        return f"RDEditorNotebook({status}, {mol_info})"
    
    def __enter__(self):
        """Context manager entry."""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()
        return False


class MoleculeEditor:
    """
    A simplified interface for quick molecular editing in Jupyter.
    
    This provides a more streamlined API focused on programmatic use.
    
    Example:
        ```python
        editor = MoleculeEditor("CCO")
        editor.edit()  # Opens GUI for editing
        new_smiles = editor.smiles
        ```
    """
    
    def __init__(self, molecule: Union[str, Chem.Mol, None] = None):
        """
        Initialize the molecule editor.
        
        Args:
            molecule: SMILES string, RDKit Mol object, or None
        """
        self._notebook_editor = RDEditorNotebook()
        
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
        return self._notebook_editor.get_molecule()
    
    @mol.setter
    def mol(self, value: Chem.Mol):
        """Set the current molecule."""
        self._notebook_editor.set_molecule(value)
    
    @property
    def smiles(self) -> Optional[str]:
        """Get the SMILES string of the current molecule."""
        return self._notebook_editor.get_smiles()
    
    @smiles.setter
    def smiles(self, value: str):
        """Set the molecule from a SMILES string."""
        self._notebook_editor.set_smiles(value)
    
    def edit(self):
        """Open the GUI editor."""
        self._notebook_editor.show()
    
    def close(self):
        """Close the editor."""
        self._notebook_editor.close()
    
    def clear(self):
        """Clear the current molecule."""
        self._notebook_editor.clear()
    
    def save(self, filename: Union[str, Path]):
        """Save the molecule to a file."""
        self._notebook_editor.save_file(filename)
    
    def load(self, filename: Union[str, Path]):
        """Load a molecule from a file."""
        self._notebook_editor.load_file(filename)
    
    def __repr__(self):
        """String representation."""
        if self.mol is not None:
            return f"MoleculeEditor(SMILES='{self.smiles}')"
        return "MoleculeEditor(empty)"


# Convenience function for quick editing
def edit_molecule(molecule: Union[str, Chem.Mol]) -> Chem.Mol:
    """
    Quick function to edit a molecule and return the result.
    
    Args:
        molecule: SMILES string or RDKit Mol object
        
    Returns:
        Edited RDKit Mol object
        
    Example:
        ```python
        from rdkit import Chem
        mol = Chem.MolFromSmiles("CCO")
        edited_mol = edit_molecule(mol)
        ```
    """
    editor = RDEditorNotebook()
    
    if isinstance(molecule, str):
        editor.set_smiles(molecule)
    elif isinstance(molecule, Chem.Mol):
        editor.set_molecule(molecule)
    else:
        raise TypeError("molecule must be a SMILES string or RDKit Mol object")
    
    editor.show()
    print("Edit the molecule in the window, then close it to continue...")
    
    return editor.get_molecule()


def quick_edit(smiles: str = "C") -> str:
    """
    Quick SMILES editing function.
    
    Args:
        smiles: Input SMILES string
        
    Returns:
        Edited SMILES string
        
    Example:
        ```python
        new_smiles = quick_edit("CCO")
        print(new_smiles)
        ```
    """
    mol = edit_molecule(smiles)
    return Chem.MolToSmiles(mol, isomericSmiles=True) if mol else ""


# Export public API
__all__ = [
    'RDEditorNotebook',
    'MoleculeEditor',
    'edit_molecule',
    'quick_edit',
]


if __name__ == "__main__":
    # Example usage
    print("RDEditor Jupyter Notebook Wrapper")
    print("-" * 50)
    
    # Example 1: Basic usage
    print("\nExample 1: Basic usage")
    editor = RDEditorNotebook()
    editor.set_smiles("CCO")
    print(f"Loaded SMILES: {editor.get_smiles()}")
    editor.show()
    
    # Example 2: Using MoleculeEditor
    print("\nExample 2: Using MoleculeEditor")
    mol_editor = MoleculeEditor("c1ccccc1")
    print(f"Molecule: {mol_editor}")
    mol_editor.edit()
