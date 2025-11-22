#!/usr/bin/env python
"""Test drawing editors for Jupyter"""

import sys
import os

# Set Qt to offscreen for testing
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

def test_imports():
    """Test basic imports"""
    print("Testing imports...")
    try:
        from rdeditor.jupyter_kekule_editor import KekuleEditor
        from rdeditor.jupyter_jsme_editor import JSMEEditor
        print("✓ All drawing editors imported successfully")
        return True
    except Exception as e:
        print(f"✗ Import failed: {e}")
        return False

def test_kekule_editor():
    """Test KekuleEditor HTML generation"""
    print("\nTesting KekuleEditor...")
    try:
        from rdeditor.jupyter_kekule_editor import KekuleEditor
        
        # Create editor
        editor = KekuleEditor("CCO")
        html = editor._generate_html()
        
        # Check for key elements
        checks = [
            ('kekule' in html.lower(), 'Kekule library reference'),
            ('getSmiles_' in html, 'getSmiles function'),
            ('clearEditor_' in html, 'clearEditor function'),
            ('loadFromSmiles_' in html, 'loadFromSmiles function'),
            ('Kekule.Editor.Composer' in html, 'Kekule Composer'),
            ('try {' in html or 'catch' in html.lower(), 'error handling'),
            ('Draw with Mouse' in html, 'drawing instructions'),
        ]
        
        all_pass = True
        for check, name in checks:
            if check:
                print(f"  ✓ {name} found")
            else:
                print(f"  ✗ {name} missing")
                all_pass = False
        
        if all_pass:
            print("✓ KekuleEditor HTML generation successful")
            return True
        else:
            print("✗ Some checks failed")
            return False
            
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_jsme_fix():
    """Test JSME editor regex fix"""
    print("\nTesting JSME editor fix...")
    try:
        from rdeditor.jupyter_jsme_editor import JSMEEditor
        
        editor = JSMEEditor("CCO")
        html = editor._generate_html()
        
        # Check that the problematic regex is fixed
        checks = [
            ("escapedSmiles" in html, 'escapedSmiles variable (proper escaping)'),
            ('initJSME_' in html, 'initJSME function'),
            ('getSmiles_' in html, 'getSmiles function'),
            ('try {' in html or 'catch' in html.lower(), 'error handling'),
        ]
        
        all_pass = True
        for check, name in checks:
            if check:
                print(f"  ✓ {name}")
            else:
                print(f"  ✗ {name}")
                all_pass = False
        
        if all_pass:
            print("✓ JSME editor regex fixed")
            return True
        else:
            print("✗ Some checks failed")
            return False
            
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_rdeditor_import():
    """Test importing from main package"""
    print("\nTesting rdeditor package import...")
    try:
        from rdeditor import KekuleEditor, JSMEEditor
        print("  ✓ KekuleEditor imported from rdeditor")
        print("  ✓ JSMEEditor imported from rdeditor")
        
        # Test instantiation
        kekule = KekuleEditor()
        jsme = JSMEEditor()
        print("  ✓ Both editors can be instantiated")
        
        return True
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests"""
    print("=" * 60)
    print("Testing Drawing Editors")
    print("=" * 60)
    
    results = []
    
    results.append(("Imports", test_imports()))
    results.append(("KekuleEditor", test_kekule_editor()))
    results.append(("JSME Fix", test_jsme_fix()))
    results.append(("Package Import", test_rdeditor_import()))
    
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
        print("\n🎉 Editors are ready for Jupyter Notebooks!")
        print("\nRecommended usage:")
        print("  from rdeditor import KekuleEditor  # Best for drawing")
        print("  editor = KekuleEditor()")
        print("  editor.display()")
        return 0
    else:
        print("✗ Some tests failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())
