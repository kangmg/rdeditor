#!/usr/bin/env python
"""
JSME-based molecule editor for Jupyter Notebooks

This module integrates JSME (JavaScript Molecular Editor) to provide
a full-featured 2D structure editor directly in Jupyter notebooks.

JSME is a free molecule editor written in JavaScript:
https://jsme-editor.github.io/

Example usage:
    ```python
    from rdeditor.jupyter_jsme_editor import JSMEEditor
    
    editor = JSMEEditor()
    editor.display()
    
    # Get edited molecule
    smiles = editor.get_smiles()
    mol = editor.get_mol()
    ```
"""

from typing import Optional
import uuid

from rdkit import Chem
from rdkit.Chem import rdDepictor

try:
    from IPython.display import display, HTML, Javascript
    IPYTHON_AVAILABLE = True
except ImportError:
    IPYTHON_AVAILABLE = False


class JSMEEditor:
    """
    JSME-based molecule editor for Jupyter notebooks.
    
    This editor embeds the JSME JavaScript Molecular Editor, providing
    a full-featured drawing interface with support for:
    - Drawing bonds and atoms
    - Stereochemistry
    - Templates (rings, functional groups)
    - Copy/paste
    - Undo/redo
    """
    
    # JSME CDN URL
    JSME_URL = "https://jsme-editor.github.io/dist/jsme/jsme.nocache.js"
    
    def __init__(self, smiles: str = "", width: int = 600, height: int = 400):
        """
        Initialize JSME editor.
        
        Args:
            smiles: Initial SMILES string
            width: Editor width in pixels
            height: Editor height in pixels
        """
        if not IPYTHON_AVAILABLE:
            raise ImportError("IPython is required for JSME editor")
        
        self.smiles = smiles
        self.width = width
        self.height = height
        self.editor_id = f"jsme_{uuid.uuid4().hex[:8]}"
        self._current_smiles = smiles
    
    def _generate_html(self) -> str:
        """Generate HTML for JSME editor."""
        html = f"""
        <div style="border: 2px solid #4CAF50; border-radius: 10px; padding: 20px; background-color: #f9f9f9;">
            <h3 style="margin-top: 0; color: #4CAF50;">🧪 Molecule Editor (JSME)</h3>
            
            <div style="margin: 10px 0;">
                <button onclick="getSmiles_{self.editor_id}()" 
                        style="background-color: #4CAF50; color: white; padding: 10px 20px; 
                               border: none; border-radius: 5px; cursor: pointer; margin-right: 10px;">
                    Get SMILES
                </button>
                <button onclick="clearEditor_{self.editor_id}()" 
                        style="background-color: #f44336; color: white; padding: 10px 20px; 
                               border: none; border-radius: 5px; cursor: pointer;">
                    Clear
                </button>
            </div>
            
            <div id="{self.editor_id}" style="margin: 20px 0;"></div>
            
            <div id="smiles_output_{self.editor_id}" 
                 style="background-color: white; padding: 15px; border-radius: 5px; 
                        border: 1px solid #ddd; font-family: monospace; word-break: break-all;">
                <b>SMILES:</b> <span id="smiles_text_{self.editor_id}">{self.smiles or 'Draw a molecule...'}</span>
            </div>
            
            <div id="props_output_{self.editor_id}" 
                 style="background-color: white; padding: 15px; border-radius: 5px; 
                        border: 1px solid #ddd; margin-top: 10px; display: none;">
                <b>Properties:</b>
                <div id="props_text_{self.editor_id}"></div>
            </div>
        </div>
        
        <script>
            // Load JSME library if not already loaded
            if (typeof JSApplet === 'undefined') {{
                var script = document.createElement('script');
                script.src = "{self.JSME_URL}";
                script.onload = function() {{
                    console.log('JSME library loaded');
                    initJSME_{self.editor_id}();
                }};
                script.onerror = function() {{
                    console.error('Failed to load JSME library');
                    document.getElementById("{self.editor_id}").innerHTML = 
                        '<div style="color: red; padding: 20px;">Failed to load JSME editor. Please check your internet connection.</div>';
                }};
                document.head.appendChild(script);
            }} else {{
                // JSME already loaded
                initJSME_{self.editor_id}();
            }}
            
            var jsmeApplet_{self.editor_id};
            
            function initJSME_{self.editor_id}() {{
                try {{
                    jsmeApplet_{self.editor_id} = new JSApplet.JSME("{self.editor_id}", "{self.width}px", "{self.height}px", {{
                        "options": "oldlook,star"
                    }});
                    
                    {'jsmeApplet_' + self.editor_id + '.readGenericMolecularInput("' + self.smiles + '");' if self.smiles else ''}
                    console.log('JSME editor initialized: {self.editor_id}');
                }} catch (error) {{
                    console.error('Error initializing JSME:', error);
                    document.getElementById("{self.editor_id}").innerHTML = 
                        '<div style="color: red; padding: 20px;">Error initializing JSME editor: ' + error.message + '</div>';
                }}
            }}
            
            function getSmiles_{self.editor_id}() {{
                try {{
                    if (!jsmeApplet_{self.editor_id}) {{
                        throw new Error('JSME editor not initialized');
                    }}
                    var smiles = jsmeApplet_{self.editor_id}.smiles();
                    document.getElementById("smiles_text_{self.editor_id}").textContent = smiles;
                    
                    // Store in Python-accessible way
                    if (typeof IPython !== 'undefined' && IPython.notebook && IPython.notebook.kernel) {{
                        IPython.notebook.kernel.execute('_jsme_smiles_{self.editor_id} = "' + smiles.replace(/\\/g, '\\\\\\\\').replace(/"/g, '\\\\"') + '"');
                    }}
                    
                    // Show success message
                    var output = document.getElementById("smiles_output_{self.editor_id}");
                    output.style.borderColor = "#4CAF50";
                    setTimeout(function() {{
                        output.style.borderColor = "#ddd";
                    }}, 1000);
                }} catch (error) {{
                    console.error('Error getting SMILES:', error);
                    alert('Error: ' + error.message);
                }}
            }}
            
            function clearEditor_{self.editor_id}() {{
                try {{
                    if (!jsmeApplet_{self.editor_id}) {{
                        throw new Error('JSME editor not initialized');
                    }}
                    jsmeApplet_{self.editor_id}.clear();
                    document.getElementById("smiles_text_{self.editor_id}").textContent = "Draw a molecule...";
                }} catch (error) {{
                    console.error('Error clearing editor:', error);
                    alert('Error: ' + error.message);
                }}
            }}
        </script>
        """
        return html
    
    def display(self):
        """Display the JSME editor in the notebook."""
        html = self._generate_html()
        display(HTML(html))
    
    def show(self):
        """Alias for display()."""
        self.display()
    
    def get_smiles(self) -> str:
        """
        Get the current SMILES from the editor.
        
        Note: This retrieves the SMILES after clicking "Get SMILES" button.
        Returns the SMILES stored in the Python variable.
        """
        try:
            from IPython import get_ipython
            ipython = get_ipython()
            var_name = f'_jsme_smiles_{self.editor_id}'
            if var_name in ipython.user_ns:
                return ipython.user_ns[var_name]
        except:
            pass
        return self._current_smiles
    
    def get_mol(self) -> Optional[Chem.Mol]:
        """
        Get the current molecule as RDKit Mol object.
        
        Returns:
            RDKit Mol object with 2D coordinates, or None if invalid
        """
        smiles = self.get_smiles()
        if not smiles or smiles == "Draw a molecule...":
            return None
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            rdDepictor.Compute2DCoords(mol)
        return mol


class RDKitDrawEditor:
    """
    Simple SVG-based editor using RDKit's drawing capabilities.
    
    This provides basic molecule visualization with clickable atoms/bonds
    for simple editing operations.
    """
    
    def __init__(self, smiles: str = "", width: int = 500, height: int = 400):
        """
        Initialize RDKit draw editor.
        
        Args:
            smiles: Initial SMILES string
            width: Image width
            height: Image height
        """
        self.smiles = smiles
        self.width = width
        self.height = height
        self._mol = None
        
        if smiles:
            self._mol = Chem.MolFromSmiles(smiles)
            if self._mol:
                rdDepictor.Compute2DCoords(self._mol)
    
    def _generate_interactive_svg(self) -> str:
        """Generate interactive SVG with RDKit."""
        from rdkit.Chem.Draw import rdMolDraw2D
        
        if self._mol is None:
            return '<svg width="500" height="400"><text x="250" y="200" text-anchor="middle" fill="#999">Enter SMILES to display molecule</text></svg>'
        
        drawer = rdMolDraw2D.MolDraw2DSVG(self.width, self.height)
        drawer.DrawMolecule(self._mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg.replace('svg:', '')
    
    def display(self):
        """Display the editor."""
        editor_id = f"rdkit_{uuid.uuid4().hex[:8]}"
        svg = self._generate_interactive_svg()
        
        html = f"""
        <div style="border: 2px solid #2196F3; border-radius: 10px; padding: 20px; background-color: #f9f9f9;">
            <h3 style="margin-top: 0; color: #2196F3;">🧪 Molecule Viewer</h3>
            
            <div style="margin: 10px 0;">
                <label><b>SMILES:</b></label><br>
                <input type="text" id="smiles_input_{editor_id}" value="{self.smiles}" 
                       style="width: 80%; padding: 8px; border: 1px solid #ddd; border-radius: 4px; font-family: monospace;">
                <button onclick="updateMolecule_{editor_id}()" 
                        style="background-color: #2196F3; color: white; padding: 8px 20px; 
                               border: none; border-radius: 4px; cursor: pointer; margin-left: 10px;">
                    Update
                </button>
            </div>
            
            <div id="mol_display_{editor_id}" style="text-align: center; margin: 20px 0; 
                                                      border: 1px solid #ddd; border-radius: 5px; 
                                                      background-color: white; padding: 10px;">
                {svg}
            </div>
            
            <div id="info_{editor_id}" style="background-color: white; padding: 15px; 
                                              border-radius: 5px; border: 1px solid #ddd;">
                <b>Molecular Formula:</b> {Chem.rdMolDescriptors.CalcMolFormula(self._mol) if self._mol else 'N/A'}<br>
                <b>Molecular Weight:</b> {Chem.Descriptors.MolWt(self._mol):.2f if self._mol else 'N/A'}<br>
                <b>Atoms:</b> {self._mol.GetNumAtoms() if self._mol else 0}<br>
                <b>Bonds:</b> {self._mol.GetNumBonds() if self._mol else 0}
            </div>
        </div>
        
        <script>
            function updateMolecule_{editor_id}() {{
                var smiles = document.getElementById("smiles_input_{editor_id}").value;
                IPython.notebook.kernel.execute(
                    'from rdeditor.jupyter_jsme_editor import RDKitDrawEditor; ' +
                    '_editor_{editor_id} = RDKitDrawEditor("' + smiles + '"); ' +
                    '_editor_{editor_id}.display()'
                );
            }}
        </script>
        """
        
        display(HTML(html))


# Export public API
__all__ = [
    'JSMEEditor',
    'RDKitDrawEditor',
]


if __name__ == "__main__":
    print("JSME Molecule Editor for Jupyter")
    print("\nUsage:")
    print("  from rdeditor.jupyter_jsme_editor import JSMEEditor")
    print("  editor = JSMEEditor()")
    print("  editor.display()")
    print("  smiles = editor.get_smiles()")
