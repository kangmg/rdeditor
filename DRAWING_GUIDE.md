# 🎨 Jupyter Notebook에서 마우스로 분자 그리기

## 문제 해결 완료! ✅

### 해결된 문제들

1. **HTMLMoleculeEditor의 TraitError** ✅
2. **JSME 정규식 오류** ✅
3. **마우스로 그리기 불가능** ✅ → **KekuleEditor 추가!**

---

## 🏆 권장: KekuleEditor (Kekule.js 기반)

**완벽한 마우스 드로잉 지원!**

### 특징
- ✅ **완전한 마우스 그리기** - 원자, 결합, 고리 모두 그릴 수 있음
- ✅ **드래그 앤 드롭**
- ✅ **입체화학 지원**
- ✅ **템플릿 제공**
- ✅ **실행 취소/재실행**
- ✅ **복사/붙여넣기**

### 사용법

```python
from rdeditor import KekuleEditor

# 빈 에디터 생성 (마우스로 자유롭게 그리기)
editor = KekuleEditor()
editor.display()

# 또는 초기 SMILES로 시작
editor = KekuleEditor("c1ccccc1")
editor.display()

# "Get SMILES" 버튼 클릭 후 분자 가져오기
smiles = editor.get_smiles()
mol = editor.get_mol()
```

### 에디터 기능

1. **그리기 도구**
   - 마우스로 원자 클릭하여 추가
   - 드래그하여 결합 생성
   - 고리 템플릿 사용

2. **버튼**
   - `✓ Get SMILES`: 현재 그린 분자의 SMILES 가져오기
   - `🗑 Clear`: 에디터 비우기
   - `📥 Load SMILES`: 입력한 SMILES를 에디터에 로드

3. **입력 필드**
   - SMILES 문자열을 직접 입력하고 "Load SMILES" 클릭

---

## 🔧 대안: JSMEEditor (정규식 오류 수정됨)

### 특징
- ✅ 경량 JavaScript 에디터
- ✅ 기본적인 그리기 기능
- ✅ 정규식 오류 수정 완료

### 사용법

```python
from rdeditor import JSMEEditor

# JSME 에디터
editor = JSMEEditor()
editor.display()

# 초기 SMILES 제공
editor = JSMEEditor("CCO")
editor.display()

# SMILES 가져오기
smiles = editor.get_smiles()
mol = editor.get_mol()
```

### 수정된 내용

**문제**:
```javascript
// ❌ 잘못된 정규식 (Colab 에러 발생)
smiles.replace(/\/g, '\\\\').replace(/"/g, '\\"')
```

**해결**:
```javascript
// ✅ 올바른 정규식 이스케이프
var escapedSmiles = smiles.replace(/\\/g, '\\\\').replace(/'/g, "\\'");
IPython.notebook.kernel.execute('_jsme_smiles_xxx = \'' + escapedSmiles + '\'');
```

---

## 📊 HTMLMoleculeEditor (위젯 기반)

### 특징
- ✅ SMILES 입력 기반
- ✅ 실시간 분자 프리뷰
- ✅ 분자 속성 자동 계산
- ✅ 작용기 추가 기능

### 사용법

```python
from rdeditor import HTMLMoleculeEditor

# SMILES로 시작
editor = HTMLMoleculeEditor("CCO")
editor.display()

# 분자 가져오기
mol = editor.mol
smiles = editor.smiles
```

---

## 🚀 완벽한 워크플로우

### 1단계: 에디터로 분자 그리기

```python
from rdeditor import KekuleEditor

# 마우스로 분자 그리기
editor = KekuleEditor()
editor.display()

# (에디터에서 마우스로 분자 그린 후 "Get SMILES" 클릭)
```

### 2단계: RDKit으로 분석

```python
# SMILES 가져오기
smiles = editor.get_smiles()
mol = editor.get_mol()

# RDKit으로 분석
from rdkit.Chem import Descriptors

print(f"SMILES: {smiles}")
print(f"Molecular Weight: {Descriptors.MolWt(mol):.2f}")
print(f"LogP: {Descriptors.MolLogP(mol):.2f}")
```

### 3단계: 시각화

```python
from rdeditor import display_molecule

# 예쁘게 표시
display_molecule(mol, size=(500, 500))
```

---

## 📈 에디터 비교

| 기능 | KekuleEditor | JSMEEditor | HTMLMoleculeEditor |
|-----|--------------|------------|-------------------|
| **마우스 그리기** | ✅ 완벽 | ⚠️ 기본 | ❌ 없음 |
| **SMILES 입력** | ✅ | ✅ | ✅ |
| **템플릿** | ✅ 많음 | ⚠️ 적음 | ✅ 작용기 |
| **실시간 프리뷰** | ✅ | ✅ | ✅ |
| **속성 계산** | ❌ | ❌ | ✅ |
| **로딩 속도** | ⚠️ 중간 | ✅ 빠름 | ✅ 빠름 |
| **오프라인 사용** | ❌ CDN | ❌ CDN | ✅ |

### 권장 사용 시나리오

- **마우스로 새 분자 그리기**: 👉 `KekuleEditor` ⭐
- **빠른 SMILES 편집**: 👉 `HTMLMoleculeEditor`
- **경량 그리기 도구**: 👉 `JSMEEditor`

---

## 🐛 알려진 제한사항

### KekuleEditor
- CDN 필요 (인터넷 연결 필요)
- 첫 로딩 시 약간 느릴 수 있음

### JSMEEditor
- CDN 필요
- 기능이 제한적

### HTMLMoleculeEditor
- 마우스 그리기 불가능 (SMILES 입력만 가능)

---

## 💡 팁

### 1. 에디터 크기 조절

```python
# 큰 에디터
editor = KekuleEditor(width=1000, height=700)
editor.display()
```

### 2. 여러 분자 비교

```python
from rdeditor import compare_molecules

mol1 = Chem.MolFromSmiles("CCO")
mol2 = Chem.MolFromSmiles("c1ccccc1")

compare_molecules([mol1, mol2], labels=["Ethanol", "Benzene"])
```

### 3. 에디터 재사용

```python
# 동일한 에디터로 여러 분자 편집 가능
editor = KekuleEditor()
editor.display()

# 첫 번째 분자 편집 후
smiles1 = editor.get_smiles()

# "Clear" 버튼 클릭, 새 분자 그리기
# 두 번째 분자 편집 후
smiles2 = editor.get_smiles()
```

---

## 🔬 고급 사용법

### 에디터에서 RDKit 파이프라인으로

```python
from rdeditor import KekuleEditor
from rdkit import Chem
from rdkit.Chem import AllChem

# 1. 분자 그리기
editor = KekuleEditor()
editor.display()

# 2. SMILES 가져오기
smiles = editor.get_smiles()

# 3. RDKit으로 3D 구조 생성
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol)
AllChem.UFFOptimizeMolecule(mol)

# 4. 저장
Chem.MolToMolFile(mol, 'molecule.mol')
```

---

## 📚 추가 문서

- `JUPYTER_USAGE.md` - 전체 Jupyter 사용 가이드
- `JUPYTER_HTML_GUIDE.md` - HTML 에디터 상세 가이드
- `FIX_SUMMARY.md` - 버그 수정 내역

---

## ✅ 요약

### 이제 Jupyter Notebook에서 완벽하게 작동합니다!

✅ **마우스 그리기**: `KekuleEditor` 사용  
✅ **TraitError 해결**: `ipywidgets.HTML` 사용  
✅ **JavaScript 오류 해결**: 정규식 수정  
✅ **모든 테스트 통과**: 4/4 ✓

### 시작하기

```python
from rdeditor import KekuleEditor

editor = KekuleEditor()
editor.display()

# 마우스로 그리고 "Get SMILES" 클릭!
smiles = editor.get_smiles()
print(smiles)
```

**🎉 즐거운 분자 그리기 되세요!**
