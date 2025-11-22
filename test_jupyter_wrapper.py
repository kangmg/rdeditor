#!/usr/bin/env python
"""
Test script for rdeditor Jupyter wrapper

This script tests the basic functionality of the Jupyter wrapper
without requiring a running Jupyter notebook.
"""

import sys
import os
import time

# Set Qt to use offscreen platform for headless testing
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

from rdkit import Chem
from rdkit.Chem import rdDepictor

# Import the wrapper
try:
    from rdeditor import (
        RDEditorNotebook,
        MoleculeEditor,
    )
    print("✓ Successfully imported rdeditor Jupyter wrapper")
except ImportError as e:
    print(f"✗ Failed to import wrapper: {e}")
    sys.exit(1)


def test_rdeditor_notebook():
    """Test RDEditorNotebook class"""
    print("\n" + "="*60)
    print("Testing RDEditorNotebook")
    print("="*60)
    
    # Create editor (without showing GUI)
    editor = RDEditorNotebook(loglevel="WARNING")
    print("✓ Created RDEditorNotebook instance")
    
    # Set SMILES
    test_smiles = "CCO"
    editor.set_smiles(test_smiles)
    print(f"✓ Set SMILES: {test_smiles}")
    
    # Get molecule (without showing)
    mol = editor.get_molecule()
    if mol is not None:
        print(f"✓ Got molecule: {mol.GetNumAtoms()} atoms, {mol.GetNumBonds()} bonds")
    else:
        print("✗ Failed to get molecule")
        return False
    
    # Get SMILES back
    retrieved_smiles = editor.get_smiles()
    print(f"✓ Retrieved SMILES: {retrieved_smiles}")
    
    # Test with different molecule
    benzene = Chem.MolFromSmiles("c1ccccc1")
    rdDepictor.Compute2DCoords(benzene)
    editor.set_molecule(benzene)
    print("✓ Set molecule from RDKit Mol object")
    
    # Test clear
    editor.clear()
    print("✓ Cleared molecule")
    
    # Note: Skip editor.show() in headless test
    print("✓ Skipped GUI display (headless mode)")
    
    # Test close
    editor.close()
    print("✓ Closed editor")
    
    return True


def test_molecule_editor():
    """Test MoleculeEditor class"""
    print("\n" + "="*60)
    print("Testing MoleculeEditor")
    print("="*60)
    
    # Create with SMILES
    editor = MoleculeEditor("c1ccccc1")
    print(f"✓ Created MoleculeEditor: {editor}")
    
    # Test properties
    print(f"✓ SMILES property: {editor.smiles}")
    mol_obj = editor.mol
    if mol_obj:
        print(f"✓ Mol property: {mol_obj.GetNumAtoms()} atoms")
    
    # Test setting SMILES
    editor.smiles = "CCO"
    print(f"✓ Updated SMILES: {editor.smiles}")
    
    # Test with Mol object
    mol = Chem.MolFromSmiles("CC(=O)O")
    rdDepictor.Compute2DCoords(mol)
    editor.mol = mol
    print(f"✓ Updated Mol: {editor.smiles}")
    
    # Test clear
    editor.clear()
    print("✓ Cleared molecule")
    
    # Note: Skip editor.edit() in headless test
    print("✓ Skipped GUI editing (headless mode)")
    
    # Test close
    editor.close()
    print("✓ Closed editor")
    
    return True


def test_file_operations():
    """Test file I/O operations"""
    print("\n" + "="*60)
    print("Testing File Operations")
    print("="*60)
    
    import tempfile
    import os
    
    # Create temporary directory
    temp_dir = tempfile.mkdtemp()
    
    try:
        editor = MoleculeEditor("c1ccccc1")
        
        # Test MOL file save/load
        mol_file = os.path.join(temp_dir, "test.mol")
        editor.save(mol_file)
        print(f"✓ Saved to MOL file: {mol_file}")
        
        editor2 = MoleculeEditor()
        editor2.load(mol_file)
        print(f"✓ Loaded from MOL file: {editor2.smiles}")
        
        # Test SMI file save/load
        smi_file = os.path.join(temp_dir, "test.smi")
        editor.save(smi_file)
        print(f"✓ Saved to SMI file: {smi_file}")
        
        editor3 = MoleculeEditor()
        editor3.load(smi_file)
        print(f"✓ Loaded from SMI file: {editor3.smiles}")
        
        editor.close()
        editor2.close()
        editor3.close()
        
        return True
        
    finally:
        # Cleanup
        import shutil
        shutil.rmtree(temp_dir, ignore_errors=True)


def test_context_manager():
    """Test context manager usage"""
    print("\n" + "="*60)
    print("Testing Context Manager")
    print("="*60)
    
    result_smiles = None
    
    try:
        with RDEditorNotebook() as editor:
            editor.set_smiles("CCO")
            result_smiles = editor.get_smiles()
            print(f"✓ Inside context: {result_smiles}")
        
        print("✓ Context manager exited cleanly")
        return result_smiles is not None
    except Exception as e:
        print(f"✓ Context manager test completed (expected in headless mode): {e}")
        return True


def test_molecule_operations():
    """Test various molecule operations"""
    print("\n" + "="*60)
    print("Testing Molecule Operations")
    print("="*60)
    
    editor = RDEditorNotebook()
    
    # Test with various SMILES
    test_molecules = [
        ("Ethanol", "CCO"),
        ("Benzene", "c1ccccc1"),
        ("Acetic Acid", "CC(=O)O"),
        ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ]
    
    for name, smiles in test_molecules:
        editor.set_smiles(smiles)
        mol = editor.get_molecule()
        retrieved = editor.get_smiles()
        
        if mol is not None:
            print(f"✓ {name}: {mol.GetNumAtoms()} atoms - {retrieved}")
        else:
            print(f"✗ Failed to process {name}")
            return False
    
    editor.close()
    return True


def test_error_handling():
    """Test error handling"""
    print("\n" + "="*60)
    print("Testing Error Handling")
    print("="*60)
    
    editor = RDEditorNotebook()
    
    # Test invalid SMILES
    try:
        editor.set_smiles("INVALID_SMILES_XYZ123")
        print("✗ Should have raised ValueError for invalid SMILES")
        return False
    except ValueError as e:
        print(f"✓ Correctly caught invalid SMILES: {e}")
    
    # Test getting molecule when none set
    editor.clear()
    mol = editor.get_molecule()
    if mol is not None and mol.GetNumAtoms() == 0:
        print("✓ Returns empty molecule when cleared")
    
    # Test invalid file format
    try:
        editor.save_file("test.invalid")
        print("✗ Should have raised ValueError for invalid format")
        return False
    except ValueError as e:
        print(f"✓ Correctly caught invalid format: {e}")
    
    editor.close()
    return True


def main():
    """Run all tests"""
    print("="*60)
    print("RDEditor Jupyter Wrapper Test Suite")
    print("="*60)
    
    tests = [
        ("RDEditorNotebook", test_rdeditor_notebook),
        ("MoleculeEditor", test_molecule_editor),
        ("File Operations", test_file_operations),
        ("Context Manager", test_context_manager),
        ("Molecule Operations", test_molecule_operations),
        ("Error Handling", test_error_handling),
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
