# from molViewWidget import molViewWidget
# from molEditWidget import MolEditWidget
# from ptable_widget import PTable
try:
    from rdeditor._version import __version__
except ImportError:  # pragma: no cover
    __version__ = "not-installed"

from .rdEditor import MainWindow

# Jupyter notebook support
try:
    # Import headless version (no Qt dependencies)
    from .jupyter_headless import (
        HeadlessMoleculeEditor,
        edit_smiles,
        convert_mol_to_smiles,
        convert_smiles_to_mol,
    )
    
    # Import HTML-based editors (recommended for Jupyter)
    try:
        from .jupyter_html_editor import (
            HTMLMoleculeEditor,
            SimpleMoleculeEditor,
        )
        _HAS_HTML_EDITOR = True
    except ImportError:
        _HAS_HTML_EDITOR = False
    
    # Import JSME editor (JavaScript-based)
    try:
        from .jupyter_jsme_editor import (
            JSMEEditor,
            RDKitDrawEditor,
        )
        _HAS_JSME = True
    except ImportError:
        _HAS_JSME = False
    
    # Import display utilities
    from .jupyter_display import (
        MoleculeDisplay,
        display_molecule,
        display_smiles,
        display_molecules,
        compare_molecules,
    )
    
    # Build __all__ based on available features
    __all__ = [
        'MainWindow',
        'HeadlessMoleculeEditor',
        'edit_smiles',
        'convert_mol_to_smiles',
        'convert_smiles_to_mol',
        'MoleculeDisplay',
        'display_molecule',
        'display_smiles',
        'display_molecules',
        'compare_molecules',
    ]
    
    # Add HTML editors (recommended)
    if _HAS_HTML_EDITOR:
        __all__.extend([
            'HTMLMoleculeEditor',
            'SimpleMoleculeEditor',
        ])
    
    # Add JSME editor
    if _HAS_JSME:
        __all__.extend([
            'JSMEEditor',
            'RDKitDrawEditor',
        ])
    
    # For backwards compatibility, keep Qt version but mark as deprecated
    try:
        from .jupyter_wrapper import (
            RDEditorNotebook,
            MoleculeEditor,
            edit_molecule,
            quick_edit,
        )
        __all__.extend([
            'RDEditorNotebook',  # Deprecated: use HTMLMoleculeEditor
            'MoleculeEditor',    # Deprecated: use HTMLMoleculeEditor
            'edit_molecule',     # Deprecated: use HTMLMoleculeEditor
            'quick_edit',        # Deprecated: use HTMLMoleculeEditor
        ])
    except (ImportError, RuntimeError):
        pass
    
except ImportError:
    # If IPython/Jupyter not available, only export MainWindow
    __all__ = ['MainWindow']
