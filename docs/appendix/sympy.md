---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# 数式処理プログラミング入門

今回は，数値計算を少し離れて記号計算を扱う。

## 記号計算

計算機による算術演算は，主に数値計算のために設計されている。
主に扱うことのできる数は，限られた精度の固定長整数や浮動小数点数である。
任意精度の整数や小数を扱うためには，何らかの工夫が必要となる。
また，有理数や無理数を扱うことは，通常の方法では行われない。

数値計算には数値計算のメリットがあるが，不便な面もある。
これまでの講義で扱ったように数値計算には数値誤差が付き物であり，有理数 $\frac{1}{3}$ や無理数 $\sqrt{2}$ のように浮動小数点数では表現できない値が存在する。
正確な値を表現する代わりに，近似値によって表現する。

```{code-cell}
import math
a = 20 / 30
b = math.sqrt(8)
print(a)
print(b)
```

方程式 $x^2+5x+5=0$ の解を求めたいとき，解の公式によって得られるのはやはり近似値である。

```{code-cell}
import math
a, b, c = 1, 5, 5
x = (- b - math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
print(x)
```

もちろん，$\displaystyle\frac{1}{3}$ は分数という表現を用いて正確に表現できるし，$\sqrt{2}$ は平方根記号 $\sqrt{\phantom{2}}$ を用いて正確に表現できる。
二次方程式の解も $\displaystyle x=\frac{-5\pm\sqrt{5}}{2}$ のように正確に表現できる。

計算機において数値だけでなく数式を扱うためには **記号計算**（symbolic computation）が有用である。
本講義では数値計算を扱うが，数値計算のためのアルゴリズムを導出・解析する過程では数式を用いることがある。
行列や導関数を扱う場合があるが，その際にも数値計算でなく数式処理が必要となる。

Python では標準ライブラリによって分数を扱うことができる（有理数モジュール [`fractions`](https://docs.python.org/3/library/fractions.html) を用いる）。
より広い範囲の数値表現を扱うためには [SymPy](https://www.sympy.org/en/index.html) ライブラリを用いるとよい。

```{code-cell}
from sympy import sympify
from IPython.display import display  # IPython 5.4 あるいは 6.1 以降では import 不要
a = sympify('20/30')
b = sympify('sqrt(8)')
display(a)
display(b)
```

以下では，SymPy の使い方の基本について説明する。

## SymPy とは

記号計算のためのオープンソースライブラリである。
GitHub リポジトリ上の記録によれば，現時点（2023年5月）での最新版は 1.11.1（2022年8月31日）である。

数・式・変数・関数などを表現するだけでなく，数式処理によって方程式を解いたり微分積分の計算をおこなったりできる。

詳しくは，以下のサイトが参考になる。

- [SymPy](https://www.sympy.org/en/index.html) 公式サイト
- [SymPy documentation](https://docs.sympy.org/latest/index.html) 公式ドキュメント
- [SymPy Live](https://live.sympy.org/) ブラウザで SymPy を試すことができる 
- [SymPy Gamma](https://www.sympygamma.com/) Wolfram Alpha を意識したインタフェース

日本語であれば以下のサイト・書籍が分かりやすい。

- [SymPy による数式処理とグラフ作成](https://home.hirosaki-u.ac.jp/jupyter/sympy/) 弘前大学 教育研究事務業務支援システム HEROIC
- [Pythonで超らくらくに数学をこなす本](https://www.amazon.co.jp/dp/4274227391/) 明松真司 2021，オーム社

SymPy は標準ライブラリではないため，Python がインストールされている環境であっても SymPy が利用できるとは限らない。
SymPy をインストールするには，環境に応じて `pip` や `conda` 等の適切なパッケージマネージャを用いてインストールすればよい。

```{note}
執筆時点において，同志社大学の演習用 PC には SymPy がインストールされている。
Google Colaboratory にも SymPy がインストールされている。
これらの環境を利用している場合，追加のインストール作業は必要ない。
```

## SymPy を利用するには

他のライブラリ同様，正しくインストールされていれば，インポートして利用できる。

```{code-cell}
import sympy
x = sympy.sympify('sqrt(8)')
print(x)
```

この講義資料では，`sympy` の省略形として `sym` を用いる。

```{code-cell}
import sympy as sym
```

特定の識別子だけをインポートする場合は，`from ... import ...` 構文を用いればよい。

```{code-cell}
from sympy import sympify
x = sympify('sqrt(8)')
print(x)
```

ワイルドカード（`*`）を用いたインポートは便利であるが，他の識別子との衝突の恐れがあるため，推奨しない。

```{code-cell}
from sympy import *
x = sympify('sqrt(8)')
print(x)
```

## SymPy オブジェクトを表示する

たとえば平方根 $\sqrt{8}$ を扱うことを考えよう。

Python の標準ライブラリ `math` を用いて近似値を計算することができる。

```{code-cell}
import math
math.sqrt(8)
```

関数 `math.sqrt()` の戻り値は `float` 型のオブジェクトである。

`math` の代わりに `sympy` を用いると以下のようになる。

```{code-cell}
import sympy as sym
sym.sqrt(8)
```

関数 `sympy.sqrt()` の戻り値は SymPy オブジェクトである。

SymPy オブジェクトは単純な値でなく，内部に構造を持っている。

SymPy を表示させるためには，いくつかの異なる方法がある。

文字列として表示する:

```{code-cell}
x = sym.sqrt(8)
print(x)
```

内部構造を表示する:

```{code-cell}
x = sym.sqrt(8)
print(sym.srepr(x))
```

ASCII 文字や Unicode 文字だけで，簡易的な数式表示をおこなう:

```{code-cell}
x = sym.sqrt(8)
sym.pprint(x, use_unicode=False)  # 半角文字だけで表示する
sym.pprint(x, use_unicode=True)   # 全角文字を利用して表示する
```

LaTeX 形式で表示する:

```{code-cell}
x = sym.sqrt(8)
print(sym.latex(x))
```

MathML 形式で表示する:

```{code-cell}
x = sym.sqrt(8)

# Content MathML （デフォルト）
print('**** Content MathML ****')
sym.print_mathml(x, printer='content')

 # Presentation MathML
print('\n**** Presentation MathML ****')
sym.print_mathml(x, printer='presentation')
```

[IPython](http://ipython.org/) 環境（IPython Console，Jupyter Notebook，Google Colaboratory，等）では，MathML 形式に変換したものを適切に処理して人間に読みやすい形で表示してくれる。
表示には，[`IPython.display.display()`](https://ipython.readthedocs.io/en/stable/api/generated/IPython.display.html#IPython.display.display) 関数を用いる。
`display()` は IPython 5.4 および 6.1 以降では自動的に import されるため，（Jupyter ノートブック等を利用する場合には）明示的に import する必要はない。

```{code-cell}
from IPython.display import display  # IPython 5.4 あるいは 6.1 以降では import 不要
x = sym.sqrt(8)
display(x)
```

あるいは単に IPython セルの最後の式に与えることで `display()` と同じ結果が得られる。

```{code-cell}
x = sym.sqrt(8)
x
```

```{note}
Python の代入文は戻り値を持たない。
したがって，Code セルの最後の命令が代入文であるときはセルの戻り値が存在せず（戻り値が `None` となり），Output セルには何も表示されない。
```

## さまざまな SymPy オブジェクト

### sympify()

SymPy オブジェクトを生成するために便利なのは `sympy.sympify()` 関数である。

```{warning}
よく似た名前の `sympy.sympify()`（SymPy 化） と `sympy.simplify()` （簡単化）が存在する。
どちらも便利でよく利用するため混同しないように注意が必要。
```

`sympify()` の引数には文字列でなく数値を与えることもできる。
しかし，混乱を避けるために文字列を与えるのがよい。

```{code-cell}
import sympy as sym
from IPython.display import display
x = sym.sympify('20 / 30')  # '20 / 30' という文字列（分数）を SymPy オブジェクト化する
y = sym.sympify(20 / 30)    # 20 / 30 を float として計算した結果を SymPy オブジェクト化する
display(x)
display(y)
```

後者は Python の `int` オブジェクトどうしで `/` 演算を行ったのち，その結果である `float` オブジェクトに対して `sympy.sympify()` 関数を適用している。
`float` の時点で既に誤差が生じているため `sympify()` によって生成した SymPy オブジェクトも誤差を含む近似値となっている。

SymPy では文字を含む式を扱うことができる。
たとえば $2(x-3)$ のような式を表す SymPy オブジェクトは以下のようにして生成できる。

```{code-cell}
sym.sympify('2 * (x - 3)')
```

表示を見ると分かるように，$2(x-3)$ が $2x-6$ に自動的に変換されている。
これは `sympy.sympify()` 関数の機能によるもので，デフォルトでは簡単な計算を行って等価な式や値に変換してくれる。
この機能を一時的に無効にしたい場合は名前付き引数 `evaluate=False` を指定すればよい。

```{code-cell}
sym.sympify('2 * (x - 3)', evaluate=False)
```

```{warning}
`sympy.sympify('2 (x - 3)')` のように `*` を省略すると機能しない。
```

あるいは，以下のように文字 $x$ を表すオブジェクトから $2(x-3)$ を生成することもできる。

```{code-cell}
x = sym.sympify('x')
2 * (x - 3)
```

どちらの方法でも生成されるオブジェクトは変わらない。

```{code-cell}
e1 = sym.sympify('2 * (x - 3)')
x = sym.sympify('x')
e2 = 2 * (x - 3)
e1 == e2
```

### 文字変数

$x$ や $y$ のような文字変数（シンボル）を SymPy オブジェクトとして生成する場合，`sympify()` 以外にもいくつかの方法がある。
以下はいずれも同じ結果をもたらす。

```python
x = sym.Symbol('x')
```

```python
x = sym.symbols('x')
# 複数の文字列を同時に生成することもできる
x, y, z = sym.symbols('x y z')
```

```python
# 1文字の変数のみ（foo のような長い名前には使えない）
from sympy.abc import x
# 複数の文字列を同時に生成することもできる
from sympy.abc import x, y, z
```

```{note}
Python の変数と SymPy のシンボルオブジェクトとの違いに注意。
たとえば，変数 `a` や `x` に値を代入した状態で `sympy.sympify('a * x ** 2')` によって SymPy オブジェクトを生成しても `a` や `x` の値は無視される。
変数 `x` と SymPy シンボル `Symbol('x')` は必ずしも関係を持たない。
```

変数の定義域を実数に限定することができる。

```{code-cell}
x = sym.symbols('x')
sym.sqrt(x ** 2)
```

```{code-cell}
x = sym.symbols('x', real=True)
sym.sqrt(x ** 2)
```

複数の文字を一気に定義したいときは，以下のように書くと便利である。

```{code-cell}
sym.symbols('a:z')
```

```{code-cell}
sym.symbols('x:10')
```

```{code-cell}
sym.symbols('x1:11')
```

### 式

文字変数オブジェクトに四則演算を適用したものも SymPy オブジェクトとなる。

```{code-cell}
x = sym.Symbol('x')
(x ** 2 + x + 1) / (x - 1)
```

平方根やその他の特殊関数については，通常の `math` モジュールのものは SymPy オブジェクトに対応していないため利用できない。
代わりに `sympy` モジュールに含まれるものを使わなければならない。

```{code-cell}
:tags: [raises-exception]
import math
x = sym.Symbol('x')
math.cos(x)  # Error
```

```{code-cell}
import math
x = sym.Symbol('x')
sym.cos(x)   # OK
```

### 等式

演算子 `=` や `==` の意味については注意が必要である。

Python では，演算子 `=` は代入演算子であって等号ではない。

下の例では，1行目の代入によって変数 `x` の値は SymPy オブジェクト $x$ となり，2行目の代入によって同じく SymPy オブジェクト $1-x$ となる。
等式オブジェクト $x=1-x$ が得られるわけではない。

```{code-cell}
x = sym.sympify('x')
x = 1 - x
x
```

Python では，演算子 `==` は等価演算子であり両辺の値が等しいかどうかを検査するために用いられる。
SymPy オブジェクトがデータとして一致するかを検査することはあっても，数学的に等価であるかどうかは考慮しない。
したがって，数学的に等価な式を `==` で比較しても `False` となる場合がある。

```{code-cell}
e1 = sym.sympify('2 * (x - 3)')
e2 = sym.sympify('2 * (x - 3)', evaluate=False)
e1 == e2
```

方程式を表す SymPy オブジェクトを生成するには `Eq` を用いる。

```{code-cell}
e1 = sym.sympify('2 * (x - 3)')
e2 = sym.sympify('2 * (x - 3)', evaluate=False)
sym.Eq(e1, e2)
```

等式を `sympify()` によって生成するときにも `Eq` が必要である。

```{code-cell}
# sym.sympify('2 * (x - 3) = 0')   # Error
# sym.sympify('2 * (x - 3) == 0')  # False
sym.sympify('Eq(2 * (x - 3), 0)')
```

### ベクトル・行列

ベクトルは以下のようにして生成できる（`Vector` でなく `Matrix` とする点に注意）。

```{code-cell}
# ベクトル x の定義（NumPy のように2次元にする必要はない）
x = sym.Matrix([1, 2])
display(x)
```

行列は以下のようにして生成できる。

```{code-cell}
# 行列 A の定義
A = sym.Matrix([
        [10, 20],
        [30, 40],
])
display(A)
```

行列やベクトルの積は以下のようにして計算できる（NumPy の `@` とは異なり `*` でよい）。

```{code-cell}
# 行列 A とベクトル x の積
A * x
```

逆行列は以下のようにして計算できる。

```{code-cell}
# 行列 A の逆行列
A.inv()
```

```{code-cell}
# 行列と逆行列の積は単位行列になる
A * A.inv()
```

行列式（determinant）

```{code-cell}
# 行列式
sym.det(A)
```

行列のデータを `list` 型や `numpy.ndarray` 型で用意してもよい。

```{code-cell}
# list から sympy.Matrix を生成する
A = [
    [10, 20],
    [30, 40],
]
A = sym.Matrix(A)
A
```

```{code-cell}
# numpy.ndarray から sympy.Matrix を生成する
import numpy as np
A = np.array([
    [10, 20],
    [30, 40],
])
A = sym.Matrix(A)
A
```

`sympify()` を利用してもよい（特に文字変数を多く含む場合は便利である）:

```{code-cell}
A = sym.sympify('Matrix([[a,b,c],[d,e,f],[g,h,i]])')
display(A)
A.inv()
```

共通因子を括りだすと:

```{code-cell}
detA = A.det()
display(sym.MatMul(1 / detA, detA * A.inv(), evaluate=False))
```

### 定数

円周率 $\pi$

```{code-cell}
pi = sym.pi
display(pi)
display(pi.evalf())  # float 型の近似値に変換する
```

虚数単位 $i$

```{code-cell}
i = sym.I
display(i)       # 虚数単位
display(i ** 2)  # 虚数単位を 2 乗すると -1 になる
```

自然対数の底（ネイピア数） $e$

```{code-cell}
e = sym.E
display(e)
display(e.evalf())  # float 型の近似値に変換する
```

### 文字変数を含むベクトル・行列

定数だけでなく変数や式を要素とするベクトルや行列を扱うこともできる。

```{code-cell}
# 変数ベクトル X の定義
X = sym.sympify('Matrix([x,y])')
display(X)
```

```{code-cell}
# 行列 A の定義
A = sym.Matrix([
        [10, 20],
        [30, 40],
])
display(A)
```

2次形式 $\boldsymbol{x}^\top \boldsymbol{A} \boldsymbol{x}$

```{code-cell}
from IPython.display import Math

# 2次形式
exp1 = sym.MatMul(X.T, A, X, evaluate=False)  # 行列積
exp2 = (X.T * A * X)[0]                       # 計算途中（結果は (1,1)-行列となる）
exp3 = sym.simplify(exp2)                     # 展開して整理する
display(Math(sym.latex(exp1) + '=' + sym.latex(exp2) + '=' + sym.latex(exp3)))
```

## 式の変形

### 因数分解

人間が手計算で方程式を解く際には，因数分解によって問題を分割するのが有効である。

たとえば，2 次式 $2 x^2 - x - 6$ を因数分解してみよう。

```{code-cell}
x = sym.symbols('x')
f1 = 2 * x ** 2 - x - 6
display(f1)
f2 = sym.factor(f1)
display(f2)
```

### 展開

$(x+2)^6$ を展開すると，

```{code-cell}
x = sym.symbols('x')
f1 = (x + 2) ** 6
display(f1)
f2 = sym.expand(f1)
display(f2)
```

### 代入

文字を含む SymPy オブジェクトについて，文字に値や式を代入することができる。

たとえば，関数 $f(x) = 2 x^2 - x - 6$ について $f(3)$ を求めてみよう。

```{code-cell}
x = sym.symbols('x')
f = 2 * x ** 2 - x - 6
f.subs(x, 3)
```

```{warning}
この場合，`f(3)` のように書いても機能しない。
変数 `f` の値は Python の関数オブジェクトではなく SymPy の式オブジェクトである。
```

関数 $f(x)$ の $x$ に式 $y+1$ を代入すると，以下のようになる。

```{code-cell}
x, y = sym.symbols('x y')
f1 = 2 * x ** 2 - x - 6
f1.subs(x, y + 1)
```

代入後の式は整理されていない。
括弧を展開して整理した式を得るには `sympy.expand()` を用いるとよい。

```{code-cell}
f2 = f1.subs(x, y + 1)
sym.expand(f2)
```

SymPy 関数オブジェクト（`sympy.Function` オブジェクト）は，引数に SymPy オブジェクトや `float` をとることができる。
混同しないように注意が必要である。

```{code-cell}
display(sym.cos(sym.pi / 4))
display(sym.cos(math.pi / 4))
```

## 方程式を解く

### 1 変数 1 次方程式

まず $1$ 変数の $1$ 次方程式 $2(x-3)=0$ を解いてみよう。

```{code-cell}
x = sym.symbols('x')
eq = sym.sympify('Eq(2 * (x - 3), 0)')
sym.solve(eq, x)
```

以下の書き方でもよい。

```{code-cell}
x = sym.symbols('x')
eq = sym.Eq(2 * (x - 3), 0)
sym.solve(eq, x)
```

文字定数 $a,b$ を含む方程式 $ax+b=0$ も同様にして解くことができる。

```{code-cell}
a, b, x = sym.symbols('a b x')
eq = sym.Eq(a * x + b, 0)
sym.solve(eq, x)
```

```{note}
`sympy.solve()` の第1引数が等式でないとき，自動的に「第1引数=0」という等式であると解釈する。
したがって， `sym.solve(a * x + b, x)` でも同じ結果が得られる。
```

### 1 変数 2 次方程式

$2$ 次方程式 $a x^2 + b x + c = 0$ を解くことで，解の公式が得られる。

```{code-cell}
a, b, c, x = sym.symbols('a b c x')
eq = sym.Eq(a * x ** 2 + b * x + c, 0)
sol = sym.solve(eq, x)
sol
```

結果を数式で表示すると

```{code-cell}
x1, x2 = sol
display(x1)
display(x2)
```

$3$ 次方程式，$4$ 次方程式の解の公式についても以下のように求めることができる。
これらは紙幅に余り表示が崩れるため表示しない。
確認したい場合は自身で試してみてほしい。

```{code-cell}
# 3次方程式
x = sym.symbols('x')
eq = sym.sympify('Eq(a * x ** 3 + b * x ** 2 + c * x + d, 0)')
sol = sym.solve(eq, x)
#display(sol[0])
```

```{code-cell}
# 4次方程式
x = sym.symbols('x')
eq = sym.sympify('Eq(a * x ** 4 + b * x ** 3 + c * x ** 2 + d * x + e, 0)')
sol = sym.solve(eq, x)
#display(sol[0])
```

### 連立方程式

続いて，$2$変数の$1$次方程式を解いてみよう。

```{code-cell}
# 変数
x, y = sym.symbols('x y')
# 方程式1
eq1 = sym.sympify('Eq(a * x + b * y, e)')
display(eq1)
# 方程式2
eq2 = sym.sympify('Eq(c * x + d * y, f)')
display(eq2)
# 解を求める
sol = sym.solve([eq1, eq2], [x, y])
# 解を表示する
for var, val in sol.items():
    display(sym.Eq(var, val))
```

変数の個数が多い場合は行列を用いて記述するのが一般的である。

```{code-cell}
# 定数や変数の行列を作成する
A = sym.sympify('Matrix([[a,b,c],[d,e,f],[g,h,i]])')
X = sym.sympify('Matrix([x,y,z])')
B = sym.sympify('Matrix([p,q,r])')
# 方程式を定義する
eq = sym.Eq(sym.MatMul(A, X), B)
# 方程式を表示する
display(eq)
# 方程式の解を求める
sol = sym.solve(eq, X)
# 解を表示する
for var, val in sol.items():
    display(sym.Eq(var, val))
```

## 微分・積分

### 微分

式を微分するには `sympy.diff()` を用いればよい。

```{code-cell}
import sympy as sym
from IPython.display import Math, display

# 変数 x を定義する
x = sym.symbols('x')

# 関数 f を定義する
f = 2 * x ** 3 - x - 6
display(f)

# 関数 f の導関数を求める
df = sym.diff(f, x)
display(df)
```

右辺だけでなく等式として表示すると:

```{code-cell}
# 変数 x と関数 f を定義する
x = sym.Symbol('x')
f = 2 * x ** 3 - x - 6
df = sym.diff(f, x)

# 関数 f を等式として表示する
lhs = sym.sympify('f(x)')
display(sym.Eq(lhs, f))

# 導関数 df/dx
lhs = sym.sympify('Derivative(f(x), x)')
display(sym.Eq(lhs, df))
```

高階微分を求めることもできる。

```{code-cell}
# 変数 x と関数 f を定義する
x = sym.symbols('x')
f = 2 * x ** 3 - x - 6
lhs = sym.sympify('f(x)')
rhs = f

for k in [1, 2, 3, 4]:
    lhs = sym.Derivative(lhs, x)
    rhs = sym.diff(rhs, x)
    display(sym.Eq(lhs, rhs))
```

高階微分を一気に求めることもできる。

```{code-cell}
# 関数 f の3次導関数（3回微分）
dddf = sym.diff(f, x, x, x)
display(Math("f'''(x)=" + sym.latex(dddf)))
```

あるいは，

```{code-cell}
dddf = sym.diff(f, x, 3)
display(Math("f'''(x)=" + sym.latex(dddf)))
```

特殊関数の微分も扱うことができる。

```{code-cell}
funcs = [
    sym.sin(x),
    sym.cos(x),
    sym.tan(x),
    sym.log(x),
    sym.exp(x),
]
for f in funcs:
    display(sym.Eq(sym.Derivative(f), sym.diff(f, x)))
```

```{note}
`sympy.diff()` と `sympy.Derivative()` の違いに注意。
前者は実際に式 $x^2$ を微分して $2x$ に変換してくれるが，後者は微分そのものを表す式 $\frac{d}{dx} x^2$ を得るだけである。
なお，後者に対しては `sympy.simplify()` を適用すると微分を求めることができる。
```

### 積分

式を積分するには `sympy.integrate()` を用いればよい。

```{code-cell}
import sympy as sym
from IPython.display import Math, display

# 変数 x を定義する
x = sym.symbols('x')

# 関数 f を定義する
f = 2 * x ** 3 - x - 6
display(Math("f(x)=" + sym.latex(f)))

# 関数 f の原始関数を求める
F = sym.integrate(f, x)
display(Math("\\int f(x) dx = " + sym.latex(F)))
```

積分定数は表示されない。

特殊関数の積分も扱うことができる。

```{code-cell}
funcs = [
    sym.sin(x),
    sym.cos(x),
    sym.tan(x),
    sym.log(x),
    sym.exp(x),
]
for f in funcs:
    display(sym.Eq(sym.Integral(f), sym.integrate(f, x)))
```

### 勾配

$2$ 変数関数の勾配（偏微分ベクトル）を求めるには，ベクトルで微分すればよい。

```{code-cell}
# スカラー変数 x, y, z を定義する
x, y, z = sym.symbols('x y z')

# ベクトル変数 x を定義する
xbm = sym.symbols('xbm')
xvec = sym.Matrix([x, y, z])

# 関数 f を定義する
fxyz = x ** 2 + y ** 2 + z ** 2
f = sym.Function('f')
exprs = [f(xbm), f(x,y,z), fxyz]
exprs = map(sym.latex, exprs)
display(Math('='.join(exprs)))

# f を X で微分する
exprs = [
    r'\nabla f',
    sym.latex(sym.Derivative(f(xbm), xbm)),
    sym.latex(sym.diff(fxyz, xvec)),
]
display(Math('='.join(exprs)))
```

## テイラー展開

与えられた関数について，$x=0$ の周りで 5 次までの近似を行う。

`sym.series()` の引数は `sym.series(関数, 変数, 点, 次数)` となる。

```{code-cell}
x = sym.symbols('x')
funcs = [
    sym.sin(x),
    sym.cos(x),
    sym.tan(x),
    sym.exp(x),
]
for f in funcs:
    display(sym.Eq(f, sym.series(f, x, 0, 6)))
```

## 関数グラフ

関数 $y=\sin(x)$ のグラフを $-2 \pi \leq x \leq 2 \pi$ の範囲について描画する。

```{code-cell}
from sympy import pi, plot, sin
from sympy.abc import x

plot(sin(x), (x, -2 * pi, 2 * pi));
```

関数 $y=\sin(x), \cos(x), \sin(x)+\cos(x)$ のグラフを $-2 \pi \leq x \leq 2 \pi$ の範囲について描画する。

```{code-cell}
import matplotlib_fontja
from sympy import pi, plot, sin, cos, tan
from sympy.abc import x

plot(sin(x), cos(x), sin(x) + cos(x), (x, -2 * pi, 2 * pi),
    legend=True,
    title='三角関数',
    xlabel='x',
    ylabel='y',
    ylim=[-2,2],
);
```
