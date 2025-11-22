#!/usr/bin/env python
"""Test HTML editor fixes for Jupyter compatibility"""

import sys
import os

# Set Qt to offscreen for testing
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

def test_imports():
    """Test basic imports"""
    print("Testing imports...")
    try:
        from rdeditor.jupyter_html_editor import HTMLMoleculeEditor, SimpleMoleculeEditor
        from rdeditor.jupyter_jsme_editor import JSMEEditor, RDKitDrawEditor
        print("✓ All imports successful")
        return True
    except Exception as e:
        print(f"✗ Import failed: {e}")
        return False

def test_html_molecule_editor():
    """Test HTMLMoleculeEditor widget creation"""
    print("\nTesting HTMLMoleculeEditor...")
    try:
        from rdeditor.jupyter_html_editor import HTMLMoleculeEditor
        
        # Create editor
        editor = HTMLMoleculeEditor("CCO")
        
        # Test that widget can be created
        widget = editor._create_editor_widget()
        
        # Check that widget is an ipywidgets widget
        from ipywidgets import Widget
        if isinstance(widget, Widget):
            print("✓ Widget created successfully")
            print(f"  - Widget type: {type(widget).__name__}")
            print(f"  - Has children: {hasattr(widget, 'children')}")
            if hasattr(widget, 'children'):
                print(f"  - Number of children: {len(widget.children)}")
                for i, child in enumerate(widget.children):
                    print(f"    Child {i}: {type(child).__name__}")
            return True
        else:
            print(f"✗ Widget is not a Widget instance: {type(widget)}")
            return False
            
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_simple_editor():
    """Test SimpleMoleculeEditor"""
    print("\nTesting SimpleMoleculeEditor...")
    try:
        from rdeditor.jupyter_html_editor import SimpleMoleculeEditor
        
        editor = SimpleMoleculeEditor("CCO")
        print(f"✓ SimpleMoleculeEditor created")
        print(f"  - SMILES: {editor.smiles}")
        print(f"  - Mol: {editor.mol is not None}")
        return True
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        return False

def test_jsme_html_generation():
    """Test JSME editor HTML generation"""
    print("\nTesting JSMEEditor HTML generation...")
    try:
        from rdeditor.jupyter_jsme_editor import JSMEEditor
        
        editor = JSMEEditor("CCO")
        html = editor._generate_html()
        
        # Check for key elements
        checks = [
            ('initJSME_' in html, 'initJSME function'),
            ('getSmiles_' in html, 'getSmiles function'),
            ('clearEditor_' in html, 'clearEditor function'),
            ('jsmeApplet_' in html, 'jsmeApplet variable'),
            ('JSApplet.JSME' in html, 'JSME constructor'),
            ('try {' in html or 'catch' in html.lower(), 'error handling'),
        ]
        
        all_pass = True
        for check, name in checks:
            if check:
                print(f"  ✓ {name} found")
            else:
                print(f"  ✗ {name} missing")
                all_pass = False
        
        if all_pass:
            print("✓ JSME HTML generation successful")
            return True
        else:
            print("✗ Some checks failed")
            return False
            
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests"""
    print("=" * 60)
    print("Testing HTML Editor Fixes")
    print("=" * 60)
    
    results = []
    
    results.append(("Imports", test_imports()))
    results.append(("HTMLMoleculeEditor", test_html_molecule_editor()))
    results.append(("SimpleMoleculeEditor", test_simple_editor()))
    results.append(("JSME HTML", test_jsme_html_generation()))
    
    print("\n" + "=" * 60)
    print("Test Summary")
    print("=" * 60)
    
    for name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {name}")
    
    all_pass = all(result for _, result in results)
    
    print("=" * 60)
    if all_pass:
        print("✓ All tests passed!")
        return 0
    else:
        print("✗ Some tests failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())
