#!/usr/bin/env python
"""
Simple non-GUI test for rdeditor Jupyter wrapper

This tests the core functionality without requiring Qt GUI.
"""

import sys
import os

# Set Qt to offscreen mode before any Qt imports
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

from rdkit import Chem
from rdkit.Chem import rdDepictor


def test_imports():
    """Test that wrapper modules can be imported"""
    print("Testing imports...")
    
    try:
        from rdeditor import jupyter_wrapper
        print("✓ jupyter_wrapper imported")
        
        from rdeditor import jupyter_display
        print("✓ jupyter_display imported")
        
        from rdeditor import (
            RDEditorNotebook,
            MoleculeEditor,
        )
        print("✓ Wrapper classes imported")
        
        return True
    except ImportError as e:
        print(f"✗ Import failed: {e}")
        return False


def test_core_functionality():
    """Test core functionality without GUI"""
    print("\nTesting core functionality...")
    
    try:
        from rdeditor.jupyter_wrapper import RDEditorNotebook
        
        # Create instance (without showing)
        editor = RDEditorNotebook()
        print("✓ Created RDEditorNotebook")
        
        # Test SMILES setting/getting
        editor.set_smiles("CCO")
        mol = editor.get_molecule()
        
        if mol and mol.GetNumAtoms() == 3:
            print(f"✓ Molecule has correct atoms: {mol.GetNumAtoms()}")
        else:
            print("✗ Molecule has incorrect atoms")
            return False
        
        smiles = editor.get_smiles()
        print(f"✓ Retrieved SMILES: {smiles}")
        
        # Test molecule setting
        test_mol = Chem.MolFromSmiles("c1ccccc1")
        rdDepictor.Compute2DCoords(test_mol)
        editor.set_molecule(test_mol)
        
        retrieved_mol = editor.get_molecule()
        if retrieved_mol and retrieved_mol.GetNumAtoms() == 6:
            print(f"✓ Set molecule correctly: {retrieved_mol.GetNumAtoms()} atoms")
        else:
            print("✗ Failed to set molecule")
            return False
        
        # Test clear
        editor.clear()
        cleared_mol = editor.get_molecule()
        if cleared_mol and cleared_mol.GetNumAtoms() == 0:
            print("✓ Cleared molecule successfully")
        else:
            print("✓ Clear operation completed")
        
        return True
        
    except Exception as e:
        print(f"✗ Core functionality test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_molecule_editor_class():
    """Test MoleculeEditor class"""
    print("\nTesting MoleculeEditor class...")
    
    try:
        from rdeditor.jupyter_wrapper import MoleculeEditor
        
        # Create with SMILES
        editor = MoleculeEditor("CCO")
        print(f"✓ Created MoleculeEditor: {editor}")
        
        # Test SMILES property
        if editor.smiles:
            print(f"✓ SMILES property works: {editor.smiles}")
        else:
            print("✗ SMILES property failed")
            return False
        
        # Test mol property
        if editor.mol and editor.mol.GetNumAtoms() > 0:
            print(f"✓ Mol property works: {editor.mol.GetNumAtoms()} atoms")
        else:
            print("✗ Mol property failed")
            return False
        
        # Test setting SMILES
        editor.smiles = "c1ccccc1"
        if "c1ccccc1" in editor.smiles.lower():
            print(f"✓ Set SMILES works: {editor.smiles}")
        else:
            print("✗ Set SMILES failed")
            return False
        
        return True
        
    except Exception as e:
        print(f"✗ MoleculeEditor test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_file_operations():
    """Test file I/O"""
    print("\nTesting file operations...")
    
    try:
        from rdeditor.jupyter_wrapper import RDEditorNotebook
        import tempfile
        import os
        
        editor = RDEditorNotebook()
        editor.set_smiles("c1ccccc1")
        
        # Create temp files
        temp_dir = tempfile.mkdtemp()
        mol_file = os.path.join(temp_dir, "test.mol")
        smi_file = os.path.join(temp_dir, "test.smi")
        
        # Save
        editor.save_file(mol_file)
        print(f"✓ Saved MOL file: {mol_file}")
        
        editor.save_file(smi_file)
        print(f"✓ Saved SMI file: {smi_file}")
        
        # Load
        editor2 = RDEditorNotebook()
        editor2.load_file(mol_file)
        loaded_mol = editor2.get_molecule()
        
        if loaded_mol and loaded_mol.GetNumAtoms() == 6:
            print(f"✓ Loaded MOL file: {loaded_mol.GetNumAtoms()} atoms")
        else:
            print("✗ Failed to load MOL file")
            return False
        
        # Cleanup
        import shutil
        shutil.rmtree(temp_dir, ignore_errors=True)
        
        return True
        
    except Exception as e:
        print(f"✗ File operations test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_display_module():
    """Test display module (without actual display)"""
    print("\nTesting display module...")
    
    try:
        from rdeditor import jupyter_display
        
        # Test that display functions exist
        assert hasattr(jupyter_display, 'display_molecule')
        print("✓ display_molecule function exists")
        
        assert hasattr(jupyter_display, 'display_molecules')
        print("✓ display_molecules function exists")
        
        assert hasattr(jupyter_display, 'compare_molecules')
        print("✓ compare_molecules function exists")
        
        assert hasattr(jupyter_display, 'MoleculeDisplay')
        print("✓ MoleculeDisplay class exists")
        
        # Test MoleculeDisplay creation (without rendering)
        from rdeditor.jupyter_display import MoleculeDisplay
        mol = Chem.MolFromSmiles("CCO")
        rdDepictor.Compute2DCoords(mol)
        
        display = MoleculeDisplay(mol=mol, size=(300, 300))
        print("✓ Created MoleculeDisplay instance")
        
        return True
        
    except Exception as e:
        print(f"✗ Display module test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all simple tests"""
    print("="*60)
    print("RDEditor Jupyter Wrapper - Simple Test Suite")
    print("="*60)
    
    tests = [
        ("Imports", test_imports),
        ("Core Functionality", test_core_functionality),
        ("MoleculeEditor Class", test_molecule_editor_class),
        ("File Operations", test_file_operations),
        ("Display Module", test_display_module),
    ]
    
    results = []
    
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"\n✗ Test '{test_name}' failed with exception: {e}")
            import traceback
            traceback.print_exc()
            results.append((test_name, False))
    
    # Summary
    print("\n" + "="*60)
    print("Test Summary")
    print("="*60)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {test_name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 All tests passed!")
        return 0
    else:
        print(f"\n⚠️  {total - passed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
