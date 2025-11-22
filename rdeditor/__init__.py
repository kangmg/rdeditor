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
    from .jupyter_wrapper import (
        RDEditorNotebook,
        MoleculeEditor,
        edit_molecule,
        quick_edit,
    )
    from .jupyter_display import (
        MoleculeDisplay,
        display_molecule,
        display_smiles,
        display_molecules,
        compare_molecules,
    )
    __all__ = [
        'MainWindow',
        'RDEditorNotebook',
        'MoleculeEditor',
        'edit_molecule',
        'quick_edit',
        'MoleculeDisplay',
        'display_molecule',
        'display_smiles',
        'display_molecules',
        'compare_molecules',
    ]
except ImportError:
    # If IPython/Jupyter not available, only export MainWindow
    __all__ = ['MainWindow']
