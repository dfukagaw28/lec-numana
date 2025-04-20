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

# 基本的な操作と演算

以下では，Python の基本（入出力，制御構文，算術演算，文字列，配列）を一通り理解していることを前提に，本科目で頻繁に登場するベクトル・行列の扱いについて述べる。

数学については，大学初年次の線形代数と微分積分の知識を前提とする。
必要に応じて補足を追加するが，不安な場合は復習しておいてほしい。

## ベクトル

以下，任意の正整数をとり $n$ とする。

本講義では，$n$ 次元ベクトルを以下のように列ベクトル（縦ベクトル）として表記する。

$$
\boldsymbol{x} = \begin{bmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{bmatrix}
$$

横書きで収まりが良いように，以下のように行ベクトル（横ベクトル）の転置として表記する場合もある。ここで，${}^\top$ は行列の転置を表す記号である。

$$
\boldsymbol{x} = \begin{bmatrix} x_1 & x_2 & \cdots & x_n \end{bmatrix}^\top
$$

```{note}
日本の高校数学では，ベクトルの成分表示においては $\boldsymbol{x} = (x_1,x_2,\cdots,x_n)$ のように表記することになっている。この表記も便宜的に用いる場合がある。
```

2 つの $n$ 次元ベクトル $\boldsymbol{a}$，$\boldsymbol{b}$ の内積を $\boldsymbol{a}\cdot\boldsymbol{b}$ あるいは $\langle\boldsymbol{a},\boldsymbol{b}\rangle$ のように表す。

$$
\boldsymbol{a} \cdot \boldsymbol{b} = \sum_{i=1}^{n} a_i b_i
$$

$n$ 次元ベクトル $\boldsymbol{a}$ の大きさを $|\boldsymbol{a}|$ と表す。

$$
|\boldsymbol{a}| = \sqrt{\boldsymbol{a}\cdot\boldsymbol{a}} = \sqrt{\sum_{i=1}^{n} a_i^2}
$$

## 行列

以下，任意の正整数の対をとり $m,n$ とする。

本講義では，$(m,n)$-行列すなわち $m$ 行 $n$ 列の行列を以下のように表記する。

$$
\boldsymbol{A}
=\begin{bmatrix}
  a_{11} & a_{12} & \cdots & a_{1n} \\
  a_{21} & a_{22} & \cdots & a_{2n} \\
  \vdots & \vdots & \ddots & \vdots \\
  a_{m1} & a_{m2} & \cdots & a_{mn} \\
\end{bmatrix}
$$

$(i,j)$ 成分が $a_{ij}$ であるような行列を $[a_{ij}]_{ij}$ あるいは単に $[a_{ij}]$ と書くことがある。

行列の和は要素ごとの和である。2 つの行列の形状（行数・列数）は一致する必要がある。

$$
\boldsymbol{A} + \boldsymbol{B} = [a_{ij} + b_{ij}]
$$

すなわち

$$
\begin{bmatrix}
  a_{11} & a_{12} & \cdots & a_{1n} \\
  a_{21} & a_{22} & \cdots & a_{2n} \\
  \vdots & \vdots & \ddots & \vdots \\
  a_{m1} & a_{m2} & \cdots & a_{mn} \\
\end{bmatrix}
+\begin{bmatrix}
  b_{11} & b_{12} & \cdots & b_{1n} \\
  b_{21} & b_{22} & \cdots & b_{2n} \\
  \vdots & \vdots & \ddots & \vdots \\
  b_{m1} & b_{m2} & \cdots & b_{mn} \\
\end{bmatrix}
=\begin{bmatrix}
  a_{11} + b_{11} & a_{12} + b_{12} & \cdots & a_{1n} + b_{1n} \\
  a_{21} + b_{21} & a_{22} + b_{22} & \cdots & a_{2n} + b_{2n}\\
  \vdots & \vdots & \ddots & \vdots \\
  a_{m1} + b_{m1} & a_{m2} + b_{m2} & \cdots & a_{mn} + b_{mn}\\
\end{bmatrix}
$$

行列の積は要素ごとの和である。行列 $\boldsymbol{A}$，$\boldsymbol{B}$ に対して積が定義できるのは，$\boldsymbol{A}$ の列数と $\boldsymbol{B}$ の行数が一致するときに限る。

$$
\boldsymbol{A} \boldsymbol{B} = \bigl[ \sum_{k} a_{ik} b_{kj} \bigr]
$$

すなわち

$$
\begin{bmatrix}
  a_{11} & a_{12} & \cdots & a_{1n} \\
  a_{21} & a_{22} & \cdots & a_{2n} \\
  \vdots & \vdots & \ddots & \vdots \\
  a_{m1} & a_{m2} & \cdots & a_{mn} \\
\end{bmatrix}
\begin{bmatrix}
  b_{11} & b_{12} & \cdots & b_{1l} \\
  b_{21} & b_{22} & \cdots & b_{2l} \\
  \vdots & \vdots & \ddots & \vdots \\
  b_{n1} & b_{n2} & \cdots & b_{nl} \\
\end{bmatrix}
=\begin{bmatrix}
  \sum_{k=1}^{n} a_{1k} b_{k1} & \sum_{k=1}^{n} a_{1k} b_{k2} & \cdots & \sum_{k=1}^{n} a_{1k} b_{kl} \\
  \sum_{k=1}^{n} a_{2k} b_{k1} & \sum_{k=1}^{n} a_{2k} b_{k2} & \cdots & \sum_{k=1}^{n} a_{2k} b_{kl} \\
  \vdots & \vdots & \ddots & \vdots \\
  \sum_{k=1}^{n} a_{mk} b_{k1} & \sum_{k=1}^{n} a_{mk} b_{k2} & \cdots & \sum_{k=1}^{n} a_{mk} b_{kl} \\
\end{bmatrix}
$$

## Python でベクトルを扱う

計算機プログラムでベクトルを扱うためには配列を用いるのが一般的である。

Python では標準ライブラリに組み込まれている `list` 型でベクトルを扱うこともできるが，NumPy パッケージの `numpy.ndarray` 型を用いるとよい。
`numpy.ndarray` は多次元配列（$n$ dimensional array）を扱うための型である。
たとえば，ベクトル

$$
x=[10,20,30]^\top
$$

は

```{code-cell}
import numpy as np
x = np.array([10, 20, 30])
print(x)
```

のように書けばよい。

このとき，多次元配列としての情報を以下のように取得することができる。

```{code-cell}
print('配列の形状 shape:', x.shape) 
print('各要素の型 dtype:', x.dtype)
```

`numpy.ndarray` の要素はすべて同じ型を持つ必要がある。
メソッド `numpy.array()` で初期化する場合，要素の型は引数から自動的に推測される。
もし小数（浮動小数点数）を扱いたい場合は以下のように小数を与えればよい。

```{code-cell}
import numpy as np
x = np.array([10., 20., 30.])
print(x)
print(x.dtype)
```

要素の型を明示的に指定したい場合は名前付き引数 `dtype` を用いることもできる。

```{code-cell}
import numpy as np
x = np.array([10, 20, 30], dtype=np.float64)
print(x)
print(x.dtype)
```

`numpy.ndarray` を用いると，ベクトルの和やスカラー倍のような操作を行うことが容易に可能である。

```{code-cell}
import numpy as np
a = np.array([10, 20, 30])
b = np.array([40, 50, 60])

# ベクトルの和
c = a + b
print('c:', c)

# スカラー倍
d = 2 * c
print('d:', d)
```

```{note}
`list` 型でも演算子 `+` や `*` が利用できるが，意味が大きく異なる。
`[10, 20, 30] + [40, 50, 60]` はリストの結合であり，結果は `[10, 20, 30, 40, 50, 60]` となる。`[10, 20] * 3` はリストの繰り返しであり `[10, 20, 10, 20, 10, 20]` となる。いずれもベクトルの演算には不向きである。
```

ベクトルの内積を求めるには `np.dot()` または `@` を用いる。

```{code-cell}
import numpy as np
a = np.array([10, 20, 30])
b = np.array([40, 50, 60])

x = np.dot(a, b)
print('dot:', x)

y = a @ b
print('@:', y)
```

## Python で行列を扱う

行列も NumPy を用いる方法が便利である。標準の `list` 型による二次元配列で表現する方法もあるが，行列特有の演算には向かない。

$$
\boldsymbol{A}
=\begin{bmatrix}
   10 & 20 & 30 \\
   40 & 50 & 60 \\
 \end{bmatrix}
$$

は

```{code-cell}
import numpy as np
a = np.array([[10, 20, 30], [40, 50, 60]])
print(a)
```

のように書くことができる。

```{code-cell}
print('配列の形状 shape:', a.shape) 
print('各要素の型 dtype:', a.dtype)
```

行列の和は，`numpy.ndarray` の機能を用いて，以下のように容易に計算できる。

```{code-cell}
import numpy as np
a = np.array([[10, 20, 30], [40, 50, 60]])
b = np.array([[20, 40, 60], [30, 50, 70]])
print('a:')
print(a)
print('b:')
print(b)
c = a + b
print('c:')
print(c)
```

行列の積を求めるには `np.dot()` または `@` を用いる。

```{code-cell}
import numpy as np
a = np.array([[1, 2, 3], [4, 5, 6]])
b = np.array([[1, 2], [3, 4], [5, 6]])
print('a:')
print(a)
print('b:')
print(b)

c = np.dot(a, b)
print('c:')
print(c)

d = a @ b
print('d:')
print(d)
```

## 多次元配列に対する基本操作

以下では `import numpy as np` を省略する。

### 特別な行列を生成する

すべての要素が $0$ であるような配列を生成するには，`np.zeros()` を用いる。第 1 引数は配列の形状（shape）を表す `tuple` 型の値である。

```{code-cell}
np.zeros((3, 4))
```

すべての要素が $1$ であるような配列を生成するには，`np.ones()` を用いる。
第 1 引数は配列の形状（shape）を表す `tuple` 型の値である。

```{code-cell}
np.ones((3, 4))
```

$n$ 次単位行列を生成するには，`np.eye()` を用いる。
第 1 引数は行列の次数（行数や列数のこと）である。

```{code-cell}
np.eye(4)
```

ランダムな値を持つ行列を生成するには，`np.random.rand()` などを用いる。
この場合，値は $0$ 以上 $1$ 未満の範囲で一様分布にしたがって生成される。

```{code-cell}
np.random.rand(3, 4)
```

```{note}
引数の括弧の付け方に注意。
`np.zeros((3, 4))` の場合，引数は 1 個で，引数の型は `tuple` である。
`np.random.rand(3, 4)` の場合，引数は 2 個で，引数の型は `int` である。

前者を `np.zeros(3, 4)` のように書くことはできない。
同様に，後者を `np.random.rand((3, 4))` のように書くことはできない。
```

### 次元・形状・サイズ

NumPy の多次元配列 `ndarray` では次元・形状・サイズという用語がよく用いられる。
それぞれ似た概念で紛らわしいため注意が必要である。

たとえば，ベクトル `np.array([1, 2, 3])` の次元は `1` であり，形状は `(3,)` であり，サイズは `3` である。

行列 `np.array([[1, 2], [3, 4], [5, 6]])` の次元は `2` であり，形状は `(3, 2)` であり，サイズは `6` である。

ベクトルの次元は $1$，行列の次元は $2$ である。物理や機械学習で用いられるテンソル（tensor）では $3$ 以上の次元を扱う。
一般に， **次元** （degree）とは，多次元配列の軸の本数と考えるとよい。
次元が $d$ である場合，1 つの要素を特定するためには $d$ 個の添字が必要である。

多次元配列 `ndarray` の **形状**（shape）とは，それぞれの軸の長さである。
次元が $d$ であるとき，形状は $d$ 個の整数から成るタプル（`tuple`）である。

多次元配列 `ndarray` の **サイズ**（size）とは，すべての要素の個数である。
形状が $(m_1,m_2,\cdots,m_d)$ ならばサイズ $n$ は $n=\prod_i m_i$ である。

これらの属性値は，以下のようにして確認することができる。

```{code-cell}
x = np.array([1, 2, 3])
print('次元:', x.ndim)
print('形状:', x.shape)
print('サイズ:', x.size)
```

```{code-cell}
a = np.array([[1, 2], [3, 4], [5, 6]])
print('次元:', a.ndim)
print('形状:', a.shape)
print('サイズ:', a.size)
```

### 型や形状を変更する


いったん生成した多次元配列の形状を変形することができる。

```{code-cell}
x = np.array([1, 2, 3])
x.astype(np.float64)
```

```{code-cell}
x = np.array([1.2, 2.4, 3.6])
x.astype(np.int64)
```

```{code-cell}
x = np.array(range(12))
x.reshape(2, 6)
```

```{note}
形状を変更した結果，次元が変化してもよい。
ただし，サイズが変わるような変更は許されない。

形状を表すタプルのうち 1 つの要素の値に `-1` を指定することで，自動的に形状を決定することもできる。たとえば形状が `(2, 6)` である配列 `x` に対して `x.reshape(3, -1)` とすることで形状は自動的に `(3, 4)` となる。
```

### 要素・行・列の抽出

以下のような 3 次正方行列を考える。

```{code-cell}
x = np.array(range(1, 10)).reshape(3, 3)
print(x)
```

このとき，$i$ 行 $j$ 列目の要素を取り出すには `x[i, j]` と書けばよい（添字が $0$ から始まる点に注意）。

```{code-cell}
print(x[0, 0])
print(x[0, 1])
...
print(x[2, 2])
```

第 $i$ 行の要素すべてをベクトルとして取り出すには `x[i, :]` と書けばよい。
このとき，`x[i]` と書いてもよい。

```{code-cell}
print(x[0, :])
print(x[1, :])
print(x[2])
```

第 $j$ 列の要素すべてをベクトルとして取り出すには `x[:, j]` と書けばよい。
このとき，`x[j]` と書くことはできない（第 $j$ 行ベクトルの意味になってしまうため）。

```{code-cell}
print(x[:, 0])
print(x[:, 1])
print(x[:, 2])
```

### 合計

Python の組み込み関数 `sum()` によって和を求めることができるが，精度に問題がある。
たとえば，以下のように $0.1, -0.5, 0.4$ の和は $0$ にならない [^sum-accuracy]。

```{code-cell}
x = [0.1, -0.5, 0.4]
print(x)

y = sum(x)
print(y)
```

組み込み関数 `sum()` よりも `numpy.sum()` 関数や `numpy.ndarray.sum()` メソッドの方が正確かつ高速な傾向にある。

```{code-cell}
x = np.array([0.1, -0.5, 0.4])
print(x)

y = np.sum(x)
print(y)

z = x.sum()
print(z)
```

[^sum-accuracy]: Python 3.12 から組み込み関数 [`sum()`](https://docs.python.org/3/library/functions.html#sum) のアルゴリズムが変わって精度が向上した（c.f. [gh-100425: Improve accuracy of builtin sum() for float inputs](https://github.com/python/cpython/issues/100425)）。Python 3.11 までは `sum([0.1] * 10)` の結果が `1` にならなかった。

ベクトルでなく行列に対して和をとることもできる。すべての要素の和，行ごとの和，列ごとの和などについて計算できる。

```{code-cell}
x = np.array([[1, 2, 3], [4, 5, 6]])
print(x)

y = np.sum(x)
print('すべての要素の和:', y)

y = np.sum(x, axis=1)
print('行ごとの和:', y)

y = np.sum(x, axis=0)
print('列ごとの和:', y)
```

### 行列の基本変型

線形代数における行列の（行）基本変形とは，以下のようなものである。

- ある行に，他の行を何倍かしたものを加える。
- ある行に $0$ でない数をかける。
- 2 つの行を入れ換える。

それぞれの操作を Python で書いてみよう。

以下では，$0$ 以上 $1$ 未満の乱数によって初期化した $(3,4)$-行列 `A` を考える。

```{code-cell}
A = np.array(range(12), dtype=np.float64).reshape(3, 4)
A
```

```{code-cell}
# 第 i 行に，第 j 行を p 倍したものを加える

i, j, p = 1, 2, 0.5

B = A.copy()
B[i] += B[j] * p
B
```

```{code-cell}
# 第 i 行に p をかける

i, p = 1, 0.5

B = A.copy()
B[i] *= p
B
```

```{code-cell}
# 第 i 行と第 j 行を入れ換える

i, j = 1, 2

B = A.copy()
B[(i, j), :] = B[(j, i), :]
B
```
