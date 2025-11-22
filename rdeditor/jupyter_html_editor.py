#!/usr/bin/env python
"""
HTML-based interactive molecule editor for Jupyter Notebooks

This module provides an HTML/JavaScript-based molecule editor that works
natively in Jupyter notebooks without requiring Qt or external windows.

Example usage:
    ```python
    from rdeditor.jupyter_html_editor import HTMLMoleculeEditor
    
    editor = HTMLMoleculeEditor("CCO")
    editor.display()  # Shows interactive editor in notebook
    ```
"""

from typing import Optional, Union
from pathlib import Path
import json
import base64
from io import BytesIO

from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor, Draw, Descriptors, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

try:
    from IPython.display import display
    from IPython.display import HTML as IPythonHTML
    from IPython.display import Javascript
    from ipywidgets import widgets, Output, VBox, HBox, Button, Text, Dropdown, HTML
    IPYWIDGETS_AVAILABLE = True
except ImportError:
    IPYWIDGETS_AVAILABLE = False
    widgets = None
    Output = None
    HTML = None
    IPythonHTML = None


class HTMLMoleculeEditor:
    """
    HTML-based molecule editor for Jupyter notebooks.
    
    This editor uses ipywidgets and HTML/JavaScript to provide an interactive
    molecule editing interface directly in the notebook cell output.
    
    Features:
    - Draw molecule structure
    - Edit SMILES directly
    - Add common functional groups
    - Save/load molecules
    - Real-time preview
    """
    
    def __init__(self, molecule: Union[str, Chem.Mol, None] = None):
        """
        Initialize the HTML molecule editor.
        
        Args:
            molecule: Initial molecule (SMILES string or RDKit Mol object)
        """
        if not IPYWIDGETS_AVAILABLE:
            raise ImportError("ipywidgets is required. Install with: pip install ipywidgets")
        
        self._mol = None
        self.output = Output()
        
        # Initialize with molecule if provided
        if molecule is not None:
            if isinstance(molecule, str):
                self.smiles = molecule
            elif isinstance(molecule, Chem.Mol):
                self.mol = molecule
    
    @property
    def mol(self) -> Optional[Chem.Mol]:
        """Get the current molecule."""
        return self._mol
    
    @mol.setter
    def mol(self, value: Chem.Mol):
        """Set the current molecule."""
        if value is not None and value.GetNumConformers() == 0:
            rdDepictor.Compute2DCoords(value)
        self._mol = value
    
    @property
    def smiles(self) -> Optional[str]:
        """Get SMILES of current molecule."""
        if self._mol is None:
            return ""
        return Chem.MolToSmiles(self._mol, isomericSmiles=True)
    
    @smiles.setter
    def smiles(self, value: str):
        """Set molecule from SMILES."""
        mol = Chem.MolFromSmiles(value)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {value}")
        self.mol = mol
    
    def _get_mol_svg(self, width=400, height=300) -> str:
        """Generate SVG representation of molecule."""
        if self._mol is None:
            return '<svg width="400" height="300"><text x="200" y="150" text-anchor="middle" fill="#999">No molecule loaded</text></svg>'
        
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(self._mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg.replace('svg:', '')
    
    def _get_mol_png_base64(self, width=400, height=300) -> str:
        """Generate base64-encoded PNG of molecule."""
        if self._mol is None:
            return ""
        
        img = Draw.MolToImage(self._mol, size=(width, height))
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        return base64.b64encode(buffered.getvalue()).decode()
    
    def _create_editor_widget(self):
        """Create the interactive editor widget."""
        
        # SMILES input
        smiles_input = Text(
            value=self.smiles or "",
            description='SMILES:',
            continuous_update=False,
            layout=widgets.Layout(width='500px')
        )
        
        # Buttons
        update_btn = Button(description='Update', button_style='primary', icon='check')
        clear_btn = Button(description='Clear', button_style='warning', icon='trash')
        
        # Common functional groups dropdown
        func_groups = Dropdown(
            options=[
                ('Select to add...', ''),
                ('Methyl (-CH3)', 'C'),
                ('Ethyl (-CH2CH3)', 'CC'),
                ('Phenyl (-C6H5)', 'c1ccccc1'),
                ('Hydroxyl (-OH)', 'O'),
                ('Carboxyl (-COOH)', 'C(=O)O'),
                ('Amino (-NH2)', 'N'),
                ('Carbonyl (=O)', '=O'),
            ],
            description='Add group:',
            layout=widgets.Layout(width='300px')
        )
        
        # Molecular properties display
        props_output = Output()
        
        # Structure display
        structure_output = Output()
        
        def update_display():
            """Update the structure and properties display."""
            # Clear outputs
            structure_output.clear_output(wait=True)
            props_output.clear_output(wait=True)
            
            # Update structure
            with structure_output:
                if self._mol is not None:
                    svg = self._get_mol_svg(500, 400)
                    display(IPythonHTML(f'<div style="text-align: center; border: 1px solid #ddd; padding: 10px; border-radius: 5px;">{svg}</div>'))
                else:
                    display(IPythonHTML('<div style="text-align: center; padding: 50px; color: #999;">No molecule to display</div>'))
            
            # Update properties
            with props_output:
                if self._mol is not None:
                    mol_wt = Chem.Descriptors.MolWt(self._mol)
                    n_atoms = self._mol.GetNumAtoms()
                    n_bonds = self._mol.GetNumBonds()
                    logp = Chem.Descriptors.MolLogP(self._mol)
                    hbd = Chem.Descriptors.NumHDonors(self._mol)
                    hba = Chem.Descriptors.NumHAcceptors(self._mol)
                    
                    props_html = f"""
                    <div style="background-color: #f5f5f5; padding: 15px; border-radius: 5px; font-family: monospace;">
                        <h4 style="margin-top: 0;">Molecular Properties</h4>
                        <table style="width: 100%;">
                            <tr><td><b>Formula:</b></td><td>{Chem.rdMolDescriptors.CalcMolFormula(self._mol)}</td></tr>
                            <tr><td><b>Molecular Weight:</b></td><td>{mol_wt:.2f} g/mol</td></tr>
                            <tr><td><b>Atoms:</b></td><td>{n_atoms}</td></tr>
                            <tr><td><b>Bonds:</b></td><td>{n_bonds}</td></tr>
                            <tr><td><b>LogP:</b></td><td>{logp:.2f}</td></tr>
                            <tr><td><b>H-Bond Donors:</b></td><td>{hbd}</td></tr>
                            <tr><td><b>H-Bond Acceptors:</b></td><td>{hba}</td></tr>
                            <tr><td><b>SMILES:</b></td><td style="word-break: break-all;">{self.smiles}</td></tr>
                        </table>
                    </div>
                    """
                    display(IPythonHTML(props_html))
                else:
                    display(IPythonHTML('<div style="color: #999; padding: 10px;">No molecule loaded</div>'))
        
        def on_update_click(b):
            """Handle update button click."""
            try:
                new_smiles = smiles_input.value.strip()
                if new_smiles:
                    self.smiles = new_smiles
                    update_display()
                    with self.output:
                        print(f"✓ Updated molecule: {new_smiles}")
            except Exception as e:
                with self.output:
                    print(f"✗ Error: {e}")
        
        def on_clear_click(b):
            """Handle clear button click."""
            self._mol = None
            smiles_input.value = ""
            update_display()
            with self.output:
                print("✓ Molecule cleared")
        
        def on_func_group_change(change):
            """Handle functional group selection."""
            if change['new'] and change['new'] != '':
                current = smiles_input.value.strip()
                if current:
                    # Append to existing SMILES
                    smiles_input.value = current + change['new']
                else:
                    smiles_input.value = change['new']
                func_groups.value = ''  # Reset dropdown
        
        # Connect event handlers
        update_btn.on_click(on_update_click)
        clear_btn.on_click(on_clear_click)
        func_groups.observe(on_func_group_change, names='value')
        
        # Initial display
        update_display()
        
        # Layout
        controls = VBox([
            HBox([smiles_input]),
            HBox([update_btn, clear_btn, func_groups]),
            self.output
        ])
        
        display_area = VBox([
            structure_output,
            props_output
        ])
        
        return VBox([
            HTML('<h3>🧪 Molecule Editor</h3>'),
            controls,
            HTML('<hr style="margin: 20px 0;">'),
            display_area
        ])
    
    def display(self):
        """Display the interactive editor in the notebook."""
        widget = self._create_editor_widget()
        display(widget)
    
    def show(self):
        """Alias for display()."""
        self.display()
    
    def save(self, filename: Union[str, Path]):
        """Save molecule to file."""
        if self._mol is None:
            raise ValueError("No molecule to save")
        
        filename = str(filename)
        if filename.lower().endswith('.mol'):
            Chem.MolToMolFile(self._mol, filename)
        elif filename.lower().endswith('.smi'):
            with open(filename, 'w') as f:
                f.write(self.smiles + '\n')
        else:
            raise ValueError(f"Unsupported file format: {filename}")
    
    def load(self, filename: Union[str, Path]):
        """Load molecule from file."""
        filename = str(filename)
        if filename.lower().endswith('.mol'):
            mol = Chem.MolFromMolFile(filename, sanitize=False)
        elif filename.lower().endswith('.smi'):
            with open(filename, 'r') as f:
                smiles = f.readline().strip()
            mol = Chem.MolFromSmiles(smiles)
        else:
            raise ValueError(f"Unsupported file format: {filename}")
        
        if mol is None:
            raise ValueError(f"Failed to load molecule from: {filename}")
        
        self.mol = mol


class SimpleMoleculeEditor:
    """
    Simplified HTML editor with basic SMILES input/output.
    
    This is a lightweight alternative that doesn't require ipywidgets.
    """
    
    def __init__(self, smiles: str = ""):
        """Initialize with optional SMILES."""
        self._smiles = smiles
        self._mol = None
        if smiles:
            self._mol = Chem.MolFromSmiles(smiles)
            if self._mol:
                rdDepictor.Compute2DCoords(self._mol)
    
    @property
    def smiles(self) -> str:
        """Get SMILES."""
        return self._smiles
    
    @smiles.setter
    def smiles(self, value: str):
        """Set SMILES."""
        mol = Chem.MolFromSmiles(value)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {value}")
        self._smiles = value
        self._mol = mol
        rdDepictor.Compute2DCoords(self._mol)
    
    @property
    def mol(self) -> Optional[Chem.Mol]:
        """Get molecule."""
        return self._mol
    
    def display(self, width=500, height=400):
        """Display molecule with editable SMILES input."""
        if self._mol is None:
            html = f"""
            <div style="border: 1px solid #ddd; padding: 20px; border-radius: 5px;">
                <h3>Molecule Editor</h3>
                <p>Enter SMILES: <input type="text" id="smiles_input" value="{self._smiles}" style="width: 300px;">
                <button onclick="updateMol()">Update</button></p>
                <div id="mol_display" style="text-align: center; min-height: 200px; padding: 20px;">
                    No molecule to display
                </div>
            </div>
            """
        else:
            drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
            drawer.DrawMolecule(self._mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText().replace('svg:', '')
            
            html = f"""
            <div style="border: 1px solid #ddd; padding: 20px; border-radius: 5px;">
                <h3>Molecule Editor</h3>
                <p><b>SMILES:</b> <code>{self._smiles}</code></p>
                <div style="text-align: center;">
                    {svg}
                </div>
                <p><b>Properties:</b></p>
                <ul>
                    <li>Atoms: {self._mol.GetNumAtoms()}</li>
                    <li>Bonds: {self._mol.GetNumBonds()}</li>
                    <li>Molecular Weight: {Chem.Descriptors.MolWt(self._mol):.2f}</li>
                </ul>
            </div>
            """
        
        display(IPythonHTML(html))


# Export public API
__all__ = [
    'HTMLMoleculeEditor',
    'SimpleMoleculeEditor',
]


if __name__ == "__main__":
    print("HTML Molecule Editor for Jupyter")
    print("Install ipywidgets: pip install ipywidgets")
    print("\nUsage:")
    print("  from rdeditor.jupyter_html_editor import HTMLMoleculeEditor")
    print("  editor = HTMLMoleculeEditor('CCO')")
    print("  editor.display()")
