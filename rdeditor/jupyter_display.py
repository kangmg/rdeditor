#!/usr/bin/env python
"""
IPython display integration for rdeditor

This module provides rich display integration for Jupyter notebooks,
including inline molecule rendering and interactive widgets.
"""

from typing import Optional
from io import BytesIO
import base64

from rdkit import Chem
from rdkit.Chem import Draw, AllChem, rdDepictor

try:
    from IPython.display import display, HTML, SVG, Image
    from IPython.core.display import DisplayObject
    IPYTHON_AVAILABLE = True
except ImportError:
    IPYTHON_AVAILABLE = False
    DisplayObject = object


class MoleculeDisplay(DisplayObject):
    """
    Display class for molecules in Jupyter notebooks.
    
    This class provides rich display output for molecules including
    2D structure images and interactive editing buttons.
    """
    
    def __init__(self, mol: Optional[Chem.Mol] = None, 
                 smiles: Optional[str] = None,
                 size: tuple = (300, 300),
                 show_edit_button: bool = True):
        """
        Initialize the molecule display.
        
        Args:
            mol: RDKit Mol object
            smiles: SMILES string (alternative to mol)
            size: Image size as (width, height)
            show_edit_button: Show interactive edit button
        """
        if mol is None and smiles is not None:
            mol = Chem.MolFromSmiles(smiles)
        
        self.mol = mol
        self.size = size
        self.show_edit_button = show_edit_button
    
    def _repr_html_(self):
        """Return HTML representation for Jupyter."""
        if self.mol is None:
            return "<p>No molecule to display</p>"
        
        # Generate SVG
        try:
            from rdkit.Chem.Draw import rdMolDraw2D
            
            # Ensure molecule has 2D coordinates
            if self.mol.GetNumConformers() == 0:
                rdDepictor.Compute2DCoords(self.mol)
            
            drawer = rdMolDraw2D.MolDraw2DSVG(self.size[0], self.size[1])
            drawer.DrawMolecule(self.mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            
            # Clean up SVG
            svg = svg.replace('svg:', '')
            
            # Get molecule info
            smiles = Chem.MolToSmiles(self.mol, isomericSmiles=True)
            num_atoms = self.mol.GetNumAtoms()
            num_bonds = self.mol.GetNumBonds()
            mol_wt = Chem.Descriptors.MolWt(self.mol)
            
            html = f"""
            <div style="border: 1px solid #ddd; border-radius: 5px; padding: 10px; 
                        display: inline-block; background-color: #f9f9f9;">
                <div style="text-align: center;">
                    {svg}
                </div>
                <div style="margin-top: 10px; font-family: monospace; font-size: 12px;">
                    <b>SMILES:</b> {smiles}<br>
                    <b>Atoms:</b> {num_atoms} | <b>Bonds:</b> {num_bonds} | 
                    <b>MW:</b> {mol_wt:.2f}
                </div>
            """
            
            if self.show_edit_button:
                html += """
                <div style="margin-top: 10px; text-align: center;">
                    <button onclick="alert('Use RDEditorNotebook().show() to edit this molecule');"
                            style="padding: 5px 15px; background-color: #4CAF50; color: white; 
                                   border: none; border-radius: 3px; cursor: pointer;">
                        Edit Molecule
                    </button>
                </div>
                """
            
            html += "</div>"
            
            return html
            
        except Exception as e:
            return f"<p>Error displaying molecule: {str(e)}</p>"
    
    def _repr_png_(self):
        """Return PNG representation for Jupyter."""
        if self.mol is None:
            return None
        
        try:
            # Ensure molecule has 2D coordinates
            if self.mol.GetNumConformers() == 0:
                rdDepictor.Compute2DCoords(self.mol)
            
            img = Draw.MolToImage(self.mol, size=self.size)
            
            # Convert PIL Image to PNG bytes
            buffer = BytesIO()
            img.save(buffer, format='PNG')
            return buffer.getvalue()
            
        except Exception:
            return None


def display_molecule(mol: Optional[Chem.Mol] = None,
                    smiles: Optional[str] = None,
                    size: tuple = (300, 300),
                    show_edit_button: bool = True):
    """
    Display a molecule in Jupyter notebook with rich formatting.
    
    Args:
        mol: RDKit Mol object
        smiles: SMILES string (alternative to mol)
        size: Image size as (width, height)
        show_edit_button: Show interactive edit button
        
    Example:
        ```python
        from rdkit import Chem
        mol = Chem.MolFromSmiles("CCO")
        display_molecule(mol)
        ```
    """
    if not IPYTHON_AVAILABLE:
        print("IPython not available. Cannot display molecule.")
        return
    
    display(MoleculeDisplay(mol=mol, smiles=smiles, size=size, 
                           show_edit_button=show_edit_button))


def display_smiles(smiles: str, **kwargs):
    """
    Display a molecule from SMILES string.
    
    Args:
        smiles: SMILES string
        **kwargs: Additional arguments passed to display_molecule
    """
    display_molecule(smiles=smiles, **kwargs)


def display_molecules(mols: list,
                     mols_per_row: int = 3,
                     size: tuple = (200, 200),
                     labels: Optional[list] = None):
    """
    Display multiple molecules in a grid.
    
    Args:
        mols: List of RDKit Mol objects
        mols_per_row: Number of molecules per row
        size: Size of each molecule image
        labels: Optional labels for each molecule
        
    Example:
        ```python
        from rdkit import Chem
        mols = [Chem.MolFromSmiles(s) for s in ["CCO", "c1ccccc1", "CC(=O)O"]]
        display_molecules(mols, labels=["Ethanol", "Benzene", "Acetic Acid"])
        ```
    """
    if not IPYTHON_AVAILABLE:
        print("IPython not available. Cannot display molecules.")
        return
    
    try:
        from rdkit.Chem.Draw import rdMolDraw2D
        
        # Prepare molecules
        prepared_mols = []
        for mol in mols:
            if mol is not None:
                if mol.GetNumConformers() == 0:
                    rdDepictor.Compute2DCoords(mol)
                prepared_mols.append(mol)
            else:
                prepared_mols.append(None)
        
        # Generate grid HTML
        html = '<div style="display: flex; flex-wrap: wrap; gap: 10px;">'
        
        for i, mol in enumerate(prepared_mols):
            label = labels[i] if labels and i < len(labels) else f"Molecule {i+1}"
            
            if mol is not None:
                drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText().replace('svg:', '')
                
                smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                
                html += f"""
                <div style="border: 1px solid #ddd; border-radius: 5px; padding: 10px; 
                            background-color: #f9f9f9; text-align: center; 
                            width: {size[0] + 20}px;">
                    <div style="font-weight: bold; margin-bottom: 5px;">{label}</div>
                    {svg}
                    <div style="margin-top: 5px; font-family: monospace; font-size: 10px; 
                                word-break: break-all;">
                        {smiles[:50]}{'...' if len(smiles) > 50 else ''}
                    </div>
                </div>
                """
            else:
                html += f"""
                <div style="border: 1px solid #ddd; border-radius: 5px; padding: 10px; 
                            background-color: #f9f9f9; text-align: center; 
                            width: {size[0] + 20}px; height: {size[1] + 60}px;">
                    <div style="font-weight: bold; margin-bottom: 5px;">{label}</div>
                    <div style="padding: 50px;">Invalid Molecule</div>
                </div>
                """
            
            # Add line break after mols_per_row molecules
            if (i + 1) % mols_per_row == 0:
                html += '<div style="width: 100%;"></div>'
        
        html += '</div>'
        
        display(HTML(html))
        
    except Exception as e:
        print(f"Error displaying molecules: {str(e)}")


def compare_molecules(mol1: Chem.Mol, mol2: Chem.Mol,
                     label1: str = "Molecule 1",
                     label2: str = "Molecule 2",
                     size: tuple = (300, 300)):
    """
    Display two molecules side by side for comparison.
    
    Args:
        mol1: First RDKit Mol object
        mol2: Second RDKit Mol object
        label1: Label for first molecule
        label2: Label for second molecule
        size: Size of each molecule image
    """
    display_molecules([mol1, mol2], mols_per_row=2, size=size, 
                     labels=[label1, label2])


# Export public API
__all__ = [
    'MoleculeDisplay',
    'display_molecule',
    'display_smiles',
    'display_molecules',
    'compare_molecules',
]
