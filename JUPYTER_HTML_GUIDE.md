# HTML-Based Molecule Editor for Jupyter Notebooks

## 개요

Jupyter Notebook에서 Qt GUI 대신 **HTML/JavaScript 기반** 인터랙티브 에디터를 사용합니다.

## 🎯 추천 에디터

### 1. HTMLMoleculeEditor (ipywidgets 기반) ⭐ 추천

ipywidgets를 사용한 완전한 인터랙티브 에디터

```python
from rdeditor import HTMLMoleculeEditor

editor = HTMLMoleculeEditor("CCO")
editor.display()  # 노트북에 바로 표시됨

# 편집 후 결과 가져오기
smiles = editor.smiles
mol = editor.mol
```

**특징:**
- ✅ 실시간 구조 미리보기
- ✅ SMILES 직접 편집
- ✅ 자주 사용하는 작용기 추가
- ✅ 분자 속성 자동 계산
- ✅ 버튼 클릭으로 업데이트
- ✅ 파일 저장/로드

**설치:**
```bash
pip install ipywidgets
jupyter nbextension enable --py widgetsnbextension
```

### 2. JSMEEditor (JavaScript 기반) ⭐⭐ 가장 강력

JSME (JavaScript Molecular Editor)를 통합한 완전한 구조 편집기

```python
from rdeditor import JSMEEditor

editor = JSMEEditor()
editor.display()  # 풀 기능 에디터 표시

# "Get SMILES" 버튼 클릭 후
smiles = editor.get_smiles()
mol = editor.get_mol()
```

**특징:**
- ✅ 마우스로 그리기 (원자, 결합)
- ✅ 입체화학 지원
- ✅ 템플릿 (고리, 작용기)
- ✅ 복사/붙여넣기
- ✅ 실행 취소/다시 실행
- ✅ 인터넷 연결만 필요 (CDN)

### 3. SimpleMoleculeEditor (경량)

가장 간단한 버전, 의존성 최소

```python
from rdeditor import SimpleMoleculeEditor

editor = SimpleMoleculeEditor("CCO")
editor.display()
```

**특징:**
- ✅ 의존성 없음 (RDKit만 필요)
- ✅ SMILES 입력/출력
- ✅ 구조 미리보기
- ✅ 빠른 로딩

### 4. HeadlessMoleculeEditor (프로그래밍 전용)

GUI 없이 프로그래밍 방식으로만 사용

```python
from rdeditor import HeadlessMoleculeEditor

editor = HeadlessMoleculeEditor("CCO")
editor.smiles = "c1ccccc1"
mol = editor.mol
```

## 📊 비교표

| 에디터 | GUI | 그리기 | 의존성 | 추천도 |
|--------|-----|--------|--------|---------|
| HTMLMoleculeEditor | ✅ | ⚠️ SMILES 입력 | ipywidgets | ⭐⭐⭐ |
| JSMEEditor | ✅ | ✅ 마우스 | 없음 (CDN) | ⭐⭐⭐⭐ |
| SimpleMoleculeEditor | ✅ | ⚠️ SMILES 입력 | 없음 | ⭐⭐ |
| HeadlessMoleculeEditor | ❌ | ❌ | 없음 | ⭐ (프로그래밍용) |

## 🚀 빠른 시작

### 옵션 1: ipywidgets 사용 (권장)

```python
# 설치
!pip install ipywidgets

# 사용
from rdeditor import HTMLMoleculeEditor

editor = HTMLMoleculeEditor("CCO")
editor.display()

# SMILES 입력란에서 편집하고 "Update" 클릭
# 분자 구조와 속성이 자동으로 업데이트됨
```

### 옵션 2: JSME 사용 (풀 기능)

```python
from rdeditor import JSMEEditor

editor = JSMEEditor("CCO")  # 초기 분자 (선택사항)
editor.display()

# 마우스로 그리기
# "Get SMILES" 버튼 클릭하여 Python으로 가져오기
```

### 옵션 3: 간단한 표시

```python
from rdeditor import SimpleMoleculeEditor

editor = SimpleMoleculeEditor("c1ccccc1")
editor.display()
```

## 💡 사용 예제

### 예제 1: 인터랙티브 편집

```python
from rdeditor import HTMLMoleculeEditor, display_molecule

# 에디터 생성 및 표시
editor = HTMLMoleculeEditor("CCO")
editor.display()

# 에디터에서 SMILES 수정하고 Update 클릭
# 그 다음 결과 가져오기
edited_smiles = editor.smiles
edited_mol = editor.mol

# 결과 표시
display_molecule(edited_mol)
```

### 예제 2: JSME로 그리기

```python
from rdeditor import JSMEEditor

# 빈 에디터로 시작
editor = JSMEEditor(width=700, height=500)
editor.display()

# 사용자가 그림 → "Get SMILES" 클릭
# 코드로 가져오기
smiles = editor.get_smiles()
mol = editor.get_mol()

print(f"그린 분자: {smiles}")
```

### 예제 3: 작용기 추가

```python
from rdeditor import HTMLMoleculeEditor

# 벤젠에서 시작
editor = HTMLMoleculeEditor("c1ccccc1")
editor.display()

# 드롭다운에서 작용기 선택하여 추가
# 예: Carboxyl (-COOH) 선택하면 자동으로 추가됨
```

### 예제 4: 여러 분자 편집

```python
from rdeditor import HTMLMoleculeEditor, display_molecules

molecules = []

for i in range(3):
    print(f"분자 {i+1} 편집:")
    editor = HTMLMoleculeEditor()
    editor.display()
    
    # 사용자가 편집 후
    input("편집 완료 후 Enter...")
    molecules.append(editor.mol)

# 모두 표시
display_molecules(molecules, labels=[f"Mol {i+1}" for i in range(3)])
```

### 예제 5: 파일 작업

```python
from rdeditor import HTMLMoleculeEditor

# 파일에서 로드
editor = HTMLMoleculeEditor()
editor.load("molecule.mol")
editor.display()

# 편집 후 저장
editor.save("edited_molecule.mol")
editor.save("edited_molecule.smi")
```

## 🔧 고급 기능

### HTMLMoleculeEditor 커스터마이징

```python
from rdeditor.jupyter_html_editor import HTMLMoleculeEditor

editor = HTMLMoleculeEditor("CCO")

# 속성 직접 접근
print(f"SMILES: {editor.smiles}")
print(f"Mol: {editor.mol}")

# 프로그래밍 방식 편집
editor.smiles = "c1ccccc1C(=O)O"  # 벤조산
editor.display()
```

### JSME 옵션

```python
from rdeditor.jupyter_jsme_editor import JSMEEditor

# 크기 조정
editor = JSMEEditor(smiles="CCO", width=800, height=600)
editor.display()

# 분자 가져오기
mol = editor.get_mol()
if mol:
    print(f"Atoms: {mol.GetNumAtoms()}")
```

## 📦 설치 요구사항

### 필수
```bash
pip install rdkit
```

### 선택 (HTMLMoleculeEditor용)
```bash
pip install ipywidgets
jupyter nbextension enable --py widgetsnbextension
```

### JSME는 추가 설치 불필요
- CDN에서 자동 로드
- 인터넷 연결만 필요

## ⚠️ 주의사항

### HTMLMoleculeEditor
- ipywidgets가 JupyterLab에서 작동하려면 추가 설정 필요:
  ```bash
  jupyter labextension install @jupyter-widgets/jupyterlab-manager
  ```

### JSMEEditor
- 첫 로드 시 JSME 라이브러리 다운로드 (1-2초)
- 오프라인 환경에서는 작동 안됨 (CDN 의존)
- `Get SMILES` 버튼 클릭 후에만 Python에서 접근 가능

### SimpleMoleculeEditor
- SMILES만 입력 가능 (그리기 불가)
- 인터랙티브 기능 제한적

## 🆚 Qt GUI vs HTML GUI

| 항목 | Qt GUI (구버전) | HTML GUI (신버전) |
|------|-----------------|-------------------|
| 표시 위치 | 별도 창 | 노트북 내부 |
| 이벤트 루프 | 필요 | 불필요 |
| 원격 서버 | ❌ 문제 | ✅ 작동 |
| 의존성 | PySide6 | ipywidgets (선택) |
| 추천도 | ⚠️ 비권장 | ✅ 권장 |

## 🎯 권장 사항

1. **풀 기능 필요**: `JSMEEditor` 사용
2. **속성 자동 계산 필요**: `HTMLMoleculeEditor` 사용
3. **빠르고 가벼운 것 필요**: `SimpleMoleculeEditor` 사용
4. **GUI 불필요**: `HeadlessMoleculeEditor` 사용

## 📚 추가 문서

- [JUPYTER_USAGE.md](./JUPYTER_USAGE.md) - 전체 사용 가이드
- [examples/rdeditor_jupyter_tutorial.ipynb](./examples/rdeditor_jupyter_tutorial.ipynb) - 튜토리얼

## 🐛 문제 해결

### ipywidgets가 표시 안됨
```bash
jupyter nbextension enable --py widgetsnbextension
# JupyterLab의 경우:
jupyter labextension install @jupyter-widgets/jupyterlab-manager
```

### JSME가 로드 안됨
- 인터넷 연결 확인
- 브라우저 콘솔에서 에러 확인 (F12)

### 분자가 표시 안됨
```python
# RDKit 설치 확인
import rdkit
print(rdkit.__version__)
```

---

**이제 Jupyter에서 완전히 HTML 기반으로 작동합니다!** 🎉
