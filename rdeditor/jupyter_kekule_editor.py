#!/usr/bin/env python
"""
Kekule.js-based molecule editor for Jupyter Notebooks

This module provides a fully-featured 2D structure editor using Kekule.js,
with full mouse drawing capabilities.

Kekule.js is a powerful open-source JavaScript library for chemical structure editing:
https://partridgejiang.github.io/Kekule.js/

Example usage:
    ```python
    from rdeditor.jupyter_kekule_editor import KekuleEditor
    
    editor = KekuleEditor()
    editor.display()
    
    # Get edited molecule
    smiles = editor.get_smiles()
    mol = editor.get_mol()
    ```
"""

from typing import Optional
import uuid
import json

from rdkit import Chem
from rdkit.Chem import rdDepictor

try:
    from IPython.display import display, HTML
    IPYTHON_AVAILABLE = True
except ImportError:
    IPYTHON_AVAILABLE = False


class KekuleEditor:
    """
    Kekule.js-based molecule editor for Jupyter notebooks.
    
    This editor embeds the Kekule.js chemical structure editor, providing
    a full-featured drawing interface with:
    - Full mouse drawing support (atoms, bonds, rings)
    - Drag and drop
    - Stereochemistry
    - Templates
    - Undo/redo
    - Copy/paste
    """
    
    # Kekule.js CDN URLs
    KEKULE_CSS = "https://cdn.jsdelivr.net/npm/kekule/dist/themes/default/kekule.css"
    KEKULE_JS = "https://cdn.jsdelivr.net/npm/kekule/dist/kekule.min.js"
    
    def __init__(self, smiles: str = "", width: int = 800, height: int = 600):
        """
        Initialize Kekule editor.
        
        Args:
            smiles: Initial SMILES string
            width: Editor width in pixels
            height: Editor height in pixels
        """
        if not IPYTHON_AVAILABLE:
            raise ImportError("IPython is required for Kekule editor")
        
        self.smiles = smiles
        self.width = width
        self.height = height
        self.editor_id = f"kekule_{uuid.uuid4().hex[:8]}"
        self._current_smiles = smiles
    
    def _generate_html(self) -> str:
        """Generate HTML for Kekule editor."""
        
        # Escape SMILES for JavaScript
        escaped_smiles = self.smiles.replace('\\', '\\\\').replace('"', '\\"').replace("'", "\\'")
        
        html = f"""
        <link rel="stylesheet" type="text/css" href="{self.KEKULE_CSS}">
        <script src="{self.KEKULE_JS}"></script>
        
        <div id="container_{self.editor_id}" style="border: 2px solid #2196F3; border-radius: 10px; padding: 20px; background-color: #f9f9f9; font-family: Arial, sans-serif;">
            <h3 style="margin-top: 0; color: #2196F3;">🧪 Molecule Editor (Kekule.js) - Draw with Mouse!</h3>
            
            <div style="margin: 15px 0;">
                <button onclick="getSmiles_{self.editor_id}()" 
                        style="background-color: #4CAF50; color: white; padding: 10px 20px; 
                               border: none; border-radius: 5px; cursor: pointer; margin-right: 10px; font-size: 14px;">
                    ✓ Get SMILES
                </button>
                <button onclick="clearEditor_{self.editor_id}()" 
                        style="background-color: #f44336; color: white; padding: 10px 20px; 
                               border: none; border-radius: 5px; cursor: pointer; margin-right: 10px; font-size: 14px;">
                    🗑 Clear
                </button>
                <button onclick="loadFromSmiles_{self.editor_id}()" 
                        style="background-color: #2196F3; color: white; padding: 10px 20px; 
                               border: none; border-radius: 5px; cursor: pointer; font-size: 14px;">
                    📥 Load SMILES
                </button>
            </div>
            
            <div style="margin: 10px 0;">
                <input type="text" id="smiles_input_{self.editor_id}" 
                       placeholder="Enter SMILES here..." 
                       value="{escaped_smiles}"
                       style="width: 100%; max-width: 600px; padding: 8px; border: 1px solid #ddd; 
                              border-radius: 4px; font-family: monospace; font-size: 13px;">
            </div>
            
            <div id="{self.editor_id}" style="width: {self.width}px; height: {self.height}px; 
                                               border: 2px solid #ddd; border-radius: 5px; 
                                               background-color: white; margin: 15px 0;"></div>
            
            <div id="output_{self.editor_id}" 
                 style="background-color: white; padding: 15px; border-radius: 5px; 
                        border: 1px solid #ddd; font-family: monospace; word-break: break-all; display: none;">
                <b style="color: #4CAF50;">SMILES:</b> <span id="smiles_text_{self.editor_id}"></span>
            </div>
            
            <div id="error_{self.editor_id}" 
                 style="background-color: #ffebee; color: #c62828; padding: 15px; 
                        border-radius: 5px; border: 1px solid #ef5350; margin-top: 10px; display: none;">
            </div>
        </div>
        
        <script>
        (function() {{
            var editor_{self.editor_id};
            var isInitialized_{self.editor_id} = false;
            
            function initKekule_{self.editor_id}() {{
                try {{
                    if (typeof Kekule === 'undefined') {{
                        setTimeout(initKekule_{self.editor_id}, 100);
                        return;
                    }}
                    
                    if (isInitialized_{self.editor_id}) {{
                        return;
                    }}
                    
                    // Create editor
                    var parentElem = document.getElementById('{self.editor_id}');
                    editor_{self.editor_id} = new Kekule.Editor.Composer(parentElem);
                    
                    // Set dimensions
                    editor_{self.editor_id}.setDimension('{self.width}px', '{self.height}px');
                    
                    // Enable all features
                    editor_{self.editor_id}.setEnableOperHistory(true);
                    editor_{self.editor_id}.setEnableToolbar(true);
                    
                    // Load initial molecule if provided
                    var initialSmiles = "{escaped_smiles}";
                    if (initialSmiles) {{
                        try {{
                            var mol = Kekule.IO.loadFormatData(initialSmiles, 'smi');
                            if (mol) {{
                                editor_{self.editor_id}.setChemObj(mol);
                            }}
                        }} catch (e) {{
                            console.warn('Could not load initial SMILES:', e);
                        }}
                    }}
                    
                    isInitialized_{self.editor_id} = true;
                    console.log('Kekule editor initialized: {self.editor_id}');
                }} catch (error) {{
                    console.error('Error initializing Kekule editor:', error);
                    document.getElementById('error_{self.editor_id}').textContent = 
                        'Error: ' + error.message;
                    document.getElementById('error_{self.editor_id}').style.display = 'block';
                }}
            }}
            
            window.getSmiles_{self.editor_id} = function() {{
                try {{
                    if (!editor_{self.editor_id}) {{
                        throw new Error('Editor not initialized');
                    }}
                    
                    var chemObj = editor_{self.editor_id}.getChemObj();
                    if (!chemObj) {{
                        throw new Error('No molecule to export');
                    }}
                    
                    var smiles = Kekule.IO.saveFormatData(chemObj, 'smi');
                    
                    // Clean up SMILES (remove whitespace and newlines)
                    smiles = smiles.trim().split('\\n')[0].trim();
                    
                    document.getElementById('smiles_text_{self.editor_id}').textContent = smiles;
                    document.getElementById('output_{self.editor_id}').style.display = 'block';
                    document.getElementById('error_{self.editor_id}').style.display = 'none';
                    
                    // Update input field
                    document.getElementById('smiles_input_{self.editor_id}').value = smiles;
                    
                    // Store in Python kernel
                    if (typeof IPython !== 'undefined' && IPython.notebook && IPython.notebook.kernel) {{
                        var escapedSmiles = smiles.replace(/\\\\/g, '\\\\\\\\').replace(/'/g, "\\\\'");
                        IPython.notebook.kernel.execute('_kekule_smiles_{self.editor_id} = \\'' + escapedSmiles + '\\'');
                    }}
                    
                    // Success animation
                    var output = document.getElementById('output_{self.editor_id}');
                    output.style.borderColor = '#4CAF50';
                    output.style.backgroundColor = '#e8f5e9';
                    setTimeout(function() {{
                        output.style.borderColor = '#ddd';
                        output.style.backgroundColor = 'white';
                    }}, 1000);
                }} catch (error) {{
                    console.error('Error getting SMILES:', error);
                    document.getElementById('error_{self.editor_id}').textContent = 
                        'Error: ' + error.message;
                    document.getElementById('error_{self.editor_id}').style.display = 'block';
                }}
            }};
            
            window.clearEditor_{self.editor_id} = function() {{
                try {{
                    if (!editor_{self.editor_id}) {{
                        throw new Error('Editor not initialized');
                    }}
                    
                    editor_{self.editor_id}.setChemObj(null);
                    document.getElementById('smiles_text_{self.editor_id}').textContent = '';
                    document.getElementById('output_{self.editor_id}').style.display = 'none';
                    document.getElementById('error_{self.editor_id}').style.display = 'none';
                    document.getElementById('smiles_input_{self.editor_id}').value = '';
                }} catch (error) {{
                    console.error('Error clearing editor:', error);
                    document.getElementById('error_{self.editor_id}').textContent = 
                        'Error: ' + error.message;
                    document.getElementById('error_{self.editor_id}').style.display = 'block';
                }}
            }};
            
            window.loadFromSmiles_{self.editor_id} = function() {{
                try {{
                    if (!editor_{self.editor_id}) {{
                        throw new Error('Editor not initialized');
                    }}
                    
                    var smiles = document.getElementById('smiles_input_{self.editor_id}').value.trim();
                    if (!smiles) {{
                        throw new Error('Please enter a SMILES string');
                    }}
                    
                    var mol = Kekule.IO.loadFormatData(smiles, 'smi');
                    if (!mol) {{
                        throw new Error('Invalid SMILES: ' + smiles);
                    }}
                    
                    editor_{self.editor_id}.setChemObj(mol);
                    document.getElementById('error_{self.editor_id}').style.display = 'none';
                    
                    // Show success
                    var input = document.getElementById('smiles_input_{self.editor_id}');
                    input.style.borderColor = '#4CAF50';
                    setTimeout(function() {{
                        input.style.borderColor = '#ddd';
                    }}, 1000);
                }} catch (error) {{
                    console.error('Error loading SMILES:', error);
                    document.getElementById('error_{self.editor_id}').textContent = 
                        'Error: ' + error.message;
                    document.getElementById('error_{self.editor_id}').style.display = 'block';
                }}
            }};
            
            // Initialize when Kekule is ready
            if (document.readyState === 'loading') {{
                document.addEventListener('DOMContentLoaded', initKekule_{self.editor_id});
            }} else {{
                initKekule_{self.editor_id}();
            }}
        }})();
        </script>
        """
        return html
    
    def display(self):
        """Display the Kekule editor in the notebook."""
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
            var_name = f'_kekule_smiles_{self.editor_id}'
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
        if not smiles:
            return None
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            rdDepictor.Compute2DCoords(mol)
        return mol


# Export public API
__all__ = [
    'KekuleEditor',
]


if __name__ == "__main__":
    print("Kekule.js Molecule Editor for Jupyter")
    print("\nUsage:")
    print("  from rdeditor.jupyter_kekule_editor import KekuleEditor")
    print("  editor = KekuleEditor()")
    print("  editor.display()")
    print("  smiles = editor.get_smiles()")
