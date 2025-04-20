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

# 連立1次方程式の数値解法(2)

前回は連立1次方程式の直接解法の1つとして，ガウス消去法を紹介した。

連立1次方程式は重要な問題であるから，ユーザーが自分でアルゴリズムを実装せずとも，既に優れた実装が存在するだろうと容易に想像できる。
実際，Python では NumPy の実装をそのまま利用するのが便利である。
NumPy ライブラリは主に C で実装されており， Python で実装したものよりかなり高速である。
以下は，その使用例である。

```{code-cell}
"""
numpy.linalg.solve() の使用例
"""

import numpy as np

n = 500
A = np.random.rand(n, n)
b = np.random.rand(n, 1)

x = np.linalg.solve(A, b)
```

この `numpy.linalg.solve()` で使用されているアルゴリズムは，LU 分解を用いるものである。

- [numpy.linalg.solve -- NumPy v1.26 Manual](https://numpy.org/doc/1.26/reference/generated/numpy.linalg.solve.html)
- [NAG Library Routine Document F07AAF (DGESV)](https://www.nag.com/numeric/fl/nagdoc_fl26/pdf/f07/f07aaf.pdf)

NumPy とよく似ているが SciPy にも同様の関数 `scipy.linalg.solve()` がある。

- [scipy.linalg.solve -- SciPy v1.13.0 Manual](https://docs.scipy.org/doc/scipy-1.13.0/reference/generated/scipy.linalg.solve.html)

SciPy と NumPy の `linalg.solve()` は同様の使い方ができるが，SciPy の方は多くのオプション引数を持っており，状況に応じてより適切な結果を得るための工夫ができるようになっている。

このような便利なものが提供されているため，あえて自ら実装する必要はないとも言える。
しかしながら，その一方で，仕組みを理解することによって，より適切に利用することができる場合もある。
たとえば，上の関数 `scipy.linalg.solve()` のオプション引数を適切に指定すればプログラムの効率化につながる。
そのような意味で数値計算アルゴリズムのおおよその流れを理解しておくことは重要である。

今回は，LU 分解を用いるいくつかの方法を紹介する。

## LU 分解法

### 単純ガウス消去法の行列形式

次の行列 $\boldsymbol{A}$ を考える。

$$
\boldsymbol{A} =
  \begin{bmatrix}
    2  & 4  & 3  \\
    -4 & -5 & -3 \\
    6  & 6  & 7
  \end{bmatrix}
$$

この行列に対して，以下の行列 $\boldsymbol{G}_1$ を考える。

$$
\boldsymbol{G}_1 =
  \begin{bmatrix}
    1             & 0 & 0 \\
    \color{red}2  & 1 & 0 \\
    \color{red}-3 & 0 & 1
  \end{bmatrix}
$$

このとき，行列の積 $\boldsymbol{G}_1 \boldsymbol{A}$ は以下のようになる。

$$
\boldsymbol{G}_1 \boldsymbol{A}
  =
  \begin{bmatrix}
    1             & 0 & 0 \\
    \color{red}2  & 1 & 0 \\
    \color{red}-3 & 0 & 1
  \end{bmatrix}
  \begin{bmatrix}
    2  & 4  & 3  \\
    -4 & -5 & -3 \\
    6  & 6  & 7
  \end{bmatrix}
  =
  \begin{bmatrix}
    2 & 4  & 3 \\
    0 & 3  & 3 \\
    0 & -6 & -2
  \end{bmatrix}
$$

第 $1$ 列に注目すると，第 $2$ 行目以降の成分はすべて $0$ になっている。
これは，前進消去の最初のステップを適用した状況に似ている。
実際，以下の 2 つの操作は等価であると考えられる。
- 行列 $\boldsymbol{G}_1$ を左からかける
- 前進消去の最初のステップ（係数行列の 2 行目以降の第 1 列目の各成分を $0$ にする）

前進消去の残りの各ステップについても，何らかの行列を左からかける操作と等価である。
このことを確認しよう。

続いて，以下の行列 $\boldsymbol{G}_2$ を考える。

$$
\boldsymbol{G}_2 =
  \begin{bmatrix}
    1 & 0            & 0 \\
    0 & 1            & 0 \\
    0 & \color{red}2 & 1
  \end{bmatrix}
$$

すると，行列の積 $\boldsymbol{G}_2 (\boldsymbol{G}_1 \boldsymbol{A})$ は

$$
\boldsymbol{G}_2 (\boldsymbol{G}_1 \boldsymbol{A})
  =
  \begin{bmatrix}
    1 & 0            & 0 \\
    0 & 1            & 0 \\
    0 & \color{red}2 & 1
  \end{bmatrix}
  \begin{bmatrix}
    2 & 4  & 3 \\
    0 & 3  & 3 \\
    0 & -6 & -2
  \end{bmatrix}
  =
  \begin{bmatrix}
    2 & 4 & 3 \\
    0 & 3 & 3 \\
    0 & 0 & 4
  \end{bmatrix}
$$

これは，前進消去を終えた直後の状態と一致する。

つまり，方程式 $\boldsymbol{A} \boldsymbol{x} = \boldsymbol{b}$ に対する単純ガウス消去法の前進消去の操作は，行列 $\boldsymbol{G}_2 \boldsymbol{G}_1$ を方程式の両辺に左からかける操作と等価である[^matmul-is-associative]。

[^matmul-is-associative]: 行列積は一般に可換性を持たないが，結合性を持つことに注意。つまり，$\boldsymbol{G}_2 (\boldsymbol{G}_1 \boldsymbol{A}) = (\boldsymbol{G}_2 \boldsymbol{G}_1) \boldsymbol{A}$ である。

一般に，行列 $\boldsymbol{A}$ が $n$ 次正方行列であるならば，前進消去の操作は行列 $\boldsymbol{G} = \boldsymbol{G}_{n-1} \cdots \boldsymbol{G}_2 \boldsymbol{G}_1$ を方程式の両辺に左からかけることと等価である。

このような行列の列 $\boldsymbol{G}_1,\boldsymbol{G}_2,\cdots,\boldsymbol{G}_{n-1}$ は次のように定義できる。

$$
  & \boldsymbol{G}_{k}
    = \begin{bmatrix}
        1      & 0      & \cdots & 0              & 0      & 0      & \cdots & 0 \\
        0      & 1      & \cdots & 0              & 0      & 0      & \cdots & 0 \\
        \vdots & \vdots & \ddots & \vdots         & \vdots & \vdots & \ddots & \vdots \\
        0      & 0      & \cdots & 1              & 0      & 0      & \cdots & 0 \\
        0      & 0      & \cdots & -g^{(k)}_{k+1} & 1      & 0      & \cdots & 0 \\
        0      & 0      & \cdots & -g^{(k)}_{k+2} & 0      & 1      & \cdots & 0 \\
        \vdots & \vdots & \ddots & \vdots         & \vdots & \vdots & \ddots & \vdots \\
        0      & 0      & \cdots & -g^{(k)}_{n}   & 0      & 0      & \cdots & 1 \\
      \end{bmatrix},
      \quad k=1,2,\cdots,n-1,
\\& g^{(k)}_i = \frac{a^{(k)}_{ik}}{a^{(k)}_{kk}}, \quad k=1,2,\cdots,n-1, \quad i=k+1,k+2,\cdots,n
$$

ただし $a^{(k)}_{ij}$ は以下のように定義される行列 $\boldsymbol{A}^{(k)}$ の第 $(i,j)$ 成分である。

$$
  & \boldsymbol{A}^{(1)} = \boldsymbol{A}
\\& \boldsymbol{A}^{(k+1)} = \boldsymbol{G}_{k} \boldsymbol{A}^{(k)}
      \quad k=1,2,\cdots,n-1
$$

````{prf:example}
:label: lu-forward-1

以下の $3$ 次正方行列 $\boldsymbol{A}$ に対して，$\boldsymbol{G}_1,\boldsymbol{G}_2$ を求めてみよう。

```{math}
  \boldsymbol{A}
  =
  \begin{bmatrix}
    2  & 4  & 3  \\
    -4 & -5 & -3 \\
    6  & 6  & 7
  \end{bmatrix}
```

定義より，

```{math}
    \boldsymbol{A}^{(1)}
  & = \boldsymbol{A}
    =
    \begin{bmatrix}
      \color{blue}2  & 4  & 3  \\
      \color{blue}-4 & -5 & -3 \\
      \color{blue}6  & 6  & 7
    \end{bmatrix}
\\  \boldsymbol{G}_1
  & =
    \begin{bmatrix}
      1       & 0 & 0 \\
      -({\color{blue}-4})/{\color{blue}2} & 1 & 0 \\
      -{\color{blue}6}/{\color{blue}2}    & 0 & 1
    \end{bmatrix}
    =
    \begin{bmatrix}
      1             & 0 & 0 \\
      \color{red}2  & 1 & 0 \\
      \color{red}-3 & 0 & 1
    \end{bmatrix}
\\  \boldsymbol{A}^{(2)}
  & = \boldsymbol{G}_1 \boldsymbol{A}^{(1)}
    =
    \begin{bmatrix}
      1  & 0 & 0 \\
      2  & 1 & 0 \\
      -3 & 0 & 1
    \end{bmatrix}
    \begin{bmatrix}
      2  & 4  & 3  \\
      -4 & -5 & -3 \\
      6  & 6  & 7
    \end{bmatrix}
    =
    \begin{bmatrix}
      2 & 4  & 3 \\
      0 & \color{blue}3  & 3 \\
      0 & \color{blue}-6 & -2
    \end{bmatrix}
\\  \boldsymbol{G}_2
  & =
    \begin{bmatrix}
      1 & 0 & 0 \\
      0 & 1 & 0 \\
      0 & -({\color{blue}-6})/{\color{blue}3} & 1
    \end{bmatrix}
    =
    \begin{bmatrix}
      1 & 0 & 0 \\
      0 & 1 & 0 \\
      0 & \color{red}2 & 1
    \end{bmatrix}
\\  \boldsymbol{A}^{(3)}
  & = \boldsymbol{G}_2 \boldsymbol{A}^{(2)}
    =
    \begin{bmatrix}
      1 & 0 & 0 \\
      0 & 1 & 0 \\
      0 & 2 & 1
    \end{bmatrix}
    \begin{bmatrix}
      2 & 4  & 3 \\
      0 & 3  & 3 \\
      0 & -6 & -2
    \end{bmatrix}
    =
    \begin{bmatrix}
      2 & 4 & 3 \\
      0 & 3 & 3 \\
      0 & 0 & 4
    \end{bmatrix}
```
````

### LU 分解

上のように定義した行列 $\boldsymbol{G}_1,\boldsymbol{G}_2,\cdots,\boldsymbol{G}_{n-1}$ を用いて，新たな行列 $\boldsymbol{G}$ を次のように定義する。

$$
\boldsymbol{G}=\boldsymbol{G}_{n-1}\boldsymbol{G}_{n-2}\cdots\boldsymbol{G}_{2}\boldsymbol{G}_{1}
$$

この行列 $\boldsymbol{G}$ を方程式 $\boldsymbol{A} \boldsymbol{x} = \boldsymbol{b}$ の両辺に左からかけると

$$
\boldsymbol{G}\boldsymbol{A}\boldsymbol{x}=\boldsymbol{G}\boldsymbol{b}
$$

となる。これが，単純ガウス消去法における前進消去を実行して得られる方程式と一致する。
前進消去の性質から，$\boldsymbol{G}\boldsymbol{A}$ は上三角行列である。
これを $\boldsymbol{U}$ とおくと

$$
\boldsymbol{U}\boldsymbol{x}=\boldsymbol{G}\boldsymbol{b}
$$ (sle-ux-gb)

と書くことができる。

定義より行列 $\boldsymbol{G}$ は正則な下三角行列である（$\because$ 正則な下三角行列の積は，正則な下三角行列である）。
正則な下三角行列 $\boldsymbol{G}$ は逆行列 $\boldsymbol{G}^{-1}$ を持つ。
逆行列 $\boldsymbol{G}^{-1}$ もやはり下三角行列である。

ここで $\boldsymbol{L}=\boldsymbol{G}^{-1}$ とおいて式 {eq}`sle-ux-gb` の両辺に左から $\boldsymbol{L}$ をかけると以下のようになる。

$$
\boldsymbol{L}\boldsymbol{U}\boldsymbol{x}
=\boldsymbol{L}\boldsymbol{G}\boldsymbol{b}
=\boldsymbol{G}^{-1}\boldsymbol{G}\boldsymbol{b}
=\boldsymbol{b}
$$

もとの方程式 $\boldsymbol{A} \boldsymbol{x} = \boldsymbol{b}$ と比較すると，$\boldsymbol{A} = \boldsymbol{L} \boldsymbol{U}$ である。
これは行列 $\boldsymbol{A}$ を下三角行列 $\boldsymbol{L}$ と上三角行列 $\boldsymbol{U}$ の積に分解することに相当する。

このような分解

$$
\boldsymbol{A} = \boldsymbol{L} \boldsymbol{U}
$$

を，**LU分解**（LU decomposition）とよぶ。


````{prf:example}
:label: lu-g1g2

上と同じ例で，$\boldsymbol{G}_1$，$\boldsymbol{G}_2$ の逆行列はそれぞれ次のようになる（符号に注意）。

```{math}
    \boldsymbol{G}^{-1}_1
  & =
    \begin{bmatrix}
      1             & 0 & 0 \\
      \color{red}-2 & 1 & 0 \\
      \color{red}3  & 0 & 1
    \end{bmatrix}
\\  \boldsymbol{G}^{-1}_2
  & =
    \begin{bmatrix}
      1 & 0             & 0 \\
      0 & 1             & 0 \\
      0 & \color{red}-2 & 1
    \end{bmatrix}
```

LU 分解の下三角行列 $\boldsymbol{L} = \boldsymbol{G}^{-1} = (\boldsymbol{G}_2 \boldsymbol{G}_1)^{-1} = \boldsymbol{G}^{-1}_1 \boldsymbol{G}^{-1}_2$ は次のようになる。

```{math}
    \boldsymbol{L}
  & = \boldsymbol{G}^{-1}_1 \boldsymbol{G}^{-1}_2
\\
  & =
    \begin{bmatrix}
      1  & 0 & 0 \\
      -2 & 1 & 0 \\
      3  & 0 & 1
    \end{bmatrix}
    \begin{bmatrix}
      1 & 0  & 0 \\
      0 & 1  & 0 \\
      0 & -2 & 1
    \end{bmatrix}
    =
    \begin{bmatrix}
      1             & 0             & 0 \\
      \color{red}-2 & 1             & 0 \\
      \color{red}3  & \color{red}-2 & 1
    \end{bmatrix}
```

ここで，$\boldsymbol{L}$ の係数が $\boldsymbol{G}^{-1}_1$，$\boldsymbol{G}^{-1}_2$ の係数から直接得られることに注意されたい。さらに，$\boldsymbol{G}^{-1}_1$ や $\boldsymbol{G}^{-1}_2$ は前進消去の過程で各行に乗じる数の符号を反転させたものである。

この性質は，この例だけでなく一般の行列（ただし LU 分解可能なものに限る）について成り立つ。

この事実と，上三角行列 $\boldsymbol{U}$ が前進消去の結果として得られる行列 $\boldsymbol{G}\boldsymbol{A}$ であることを合わせると，LU 分解はガウス消去法の前進消去によって得られることが分かる。
````

LU 分解について以下のことが知られている。

```{prf:theorem}
正則な行列は LU 分解可能である。
```

LU 分解は一意であるとは限らないが，次のように少し条件を加えると一意に定まる。

```{prf:theorem}
正則な行列 $\boldsymbol{A}$ に対する LU 分解 $\boldsymbol{A} = \boldsymbol{L} \boldsymbol{U}$ のうち，下三角行列 $\boldsymbol{L}$ の対角成分がすべて $1$ であるようなものは一意に定まる。
```

LU 分解によって，方程式 $\boldsymbol{A} \boldsymbol{x} = \boldsymbol{b}$ は次のように2つの方程式に分解することができる。

$$
\boldsymbol{L} \boldsymbol{y} = \boldsymbol{b}
\\
\boldsymbol{U} \boldsymbol{x} = \boldsymbol{y}
$$

この2つの方程式は，係数が三角行列であるため扱いやすく，後退代入と同様の手続きによって効率良く解くことができる。後者の方程式を解く過程は，ガウス消去法の後退代入そのものである。前者の方程式を解く過程を **前進代入**（forward substitution）ともいう。

以上をふまえると，方程式 $\boldsymbol{A} \boldsymbol{x} = \boldsymbol{b}$ を LU 分解に基づいて解くためのアルゴリズムは，次のようになる。

```{prf:algorithm} LU分解法
:label: alg:lu

**Inputs** $\boldsymbol{A}$: $n$次正方行列，$\boldsymbol{b}$: $n$次元ベクトル

**Outputs** 方程式 $\boldsymbol{A} \boldsymbol{x} = \boldsymbol{b}$ の 解 $\boldsymbol{x}$

1. $\boldsymbol{A}$ の LU 分解 $\boldsymbol{A} = \boldsymbol{L} \boldsymbol{U}$ を求める。
2. 方程式 $\boldsymbol{L} \boldsymbol{y} = \boldsymbol{b}$ を $\boldsymbol{y}$ について解く。
3. 方程式 $\boldsymbol{U} \boldsymbol{x} = \boldsymbol{y}$ を $\boldsymbol{x}$ について解く。
4. 解 $\boldsymbol{x}$ を返す。
```

LU 分解法は，ガウス消去法と比べて手間が減るわけではなく，ほぼ同じである[^lu-complexity]。
しかしながら，同じ係数行列 $\boldsymbol{A}$ を持つ複数の方程式 $\boldsymbol{A} \boldsymbol{x} = \boldsymbol{b_1}$，$\boldsymbol{A} \boldsymbol{x} = \boldsymbol{b_2}$，$\boldsymbol{A} \boldsymbol{x} = \boldsymbol{b_3}$，$\cdots$ などを解く必要があるとき，LU 分解の結果を流用でき，全体の手間を大幅に減らすことができる。

[^lu-complexity]: 行列 $\boldsymbol{A}$ が $n$ 次正方行列のとき，単純ガウス消去法の計算量は $O(n^3)$ である。(a1) 前進消去の左辺の処理，(a2) 右辺の処理，(a3) 後退代入がそれぞれ (b1) LU 分解，(b2) 方程式 $\boldsymbol{L} \boldsymbol{y} = \boldsymbol{b}$ の求解，(b3) 方程式 $\boldsymbol{U} \boldsymbol{x} = \boldsymbol{y}$ の求解に概ね相当する。

```{note}
ここでは単純ガウス消去法を用いて説明したが，部分ピボット選択付きガウス消去法についても同様に行列形式で表すことができ，LU 分解に帰着することができる。具体的には，置換行列 $\boldsymbol{P}$ を用いて $\boldsymbol{A}$ を $\boldsymbol{A}' = \boldsymbol{P} \boldsymbol{A}$ に変換してから LU 分解 $\boldsymbol{P} \boldsymbol{A} = \boldsymbol{L} \boldsymbol{U}$ を行えばよい。
```

LU分解法の計算機プログラムを実装する際には，ガウス消去法をほぼそのまま利用することができる。
ただし，多くの成分が $0$ であることが分かっている三角行列 $\boldsymbol{L}$，$\boldsymbol{U}$ をそのまま格納すると無駄が大きいことから，下三角行列 $\boldsymbol{L}$ と上三角行列 $\boldsymbol{U}$ を 1 つの行列に重ねて格納する。

具体的には，

$$
\boldsymbol{L}
= \begin{bmatrix}
  1         & 0         & 0         & \cdots & 0         & 0 \\
  l_{21}    & 1         & 0         & \cdots & 0         & 0 \\
  l_{31}    & l_{32}    & 1         & \cdots & 0         & 0 \\
  \vdots    & \vdots    & \vdots    & \ddots & \vdots    & \vdots\\
  l_{n-1,1} & l_{n-1,2} & l_{n-1,3} & \cdots & 1         & 0 \\
  l_{n1}    & l_{n2}    & l_{n3}    & \cdots & l_{n,n-1} & 1 \\
  \end{bmatrix}
\\
\boldsymbol{U}
= \begin{bmatrix}
  u_{11} & u_{12} & u_{13} & \cdots & u_{1,n-1}   & u_{1n} \\
  0      & u_{22} & u_{23} & \cdots & u_{2,n-1}   & u_{2n} \\
  0      & 0      & u_{33} & \cdots & u_{3,n-1}   & u_{3n} \\
  \vdots & \vdots & \vdots & \ddots & \vdots      & \vdots \\
  0      & 0      & 0      & \cdots & u_{n-1,n-1} & u_{n-1,n} \\
  0      & 0      & 0      & \cdots & 0           & u_{nn} \\
  \end{bmatrix}
$$

という 2 つの行列を以下のように格納する。

$$
  \begin{bmatrix}
  u_{11}    & u_{12}    & u_{13}    & \cdots & u_{1,n-1}   & u_{1n} \\
  l_{21}    & u_{22}    & u_{23}    & \cdots & u_{2,n-1}   & u_{2n} \\
  l_{31}    & l_{32}    & u_{33}    & \cdots & u_{3,n-1}   & u_{3n} \\
  \vdots    & \vdots    & \vdots    & \ddots & \vdots      & \vdots\\
  l_{n-1,1} & l_{n-1,2} & l_{n-1,3} & \cdots & u_{n-1,n-1} & u_{n-1,n} \\
  l_{n1}    & l_{n2}    & l_{n3}    & \cdots & l_{n,n-1}   & u_{nn} \\
  \end{bmatrix}
$$

$\boldsymbol{L}$ の対角成分は $1$ であることが分かっているため省略して問題ない。
格納先は，新たな領域を用意せず，行列 $\boldsymbol{A}$ に上書きして格納するとさらに効率が良い。

Python による実装例を示す。

```{code-cell}
"""
LU 分解法の実装例
"""

import numpy as np

# 部分ピボット選択
def pivot(A, b, k):
    h = k + np.abs(A[k:, k]).argmax()
    if h != k:
        A[[k, h]] = A[[h, k]]
        b[[k, h]] = b[[h, k]]

# 前進消去
def forward_elim(A, b):
    n = A.shape[0]
    for k in range(n-1):
        pivot(A, b, k)
        A[k+1:, k] /= A[k, k]
        A[k+1:, k+1:] -= A[k+1:, k:k+1] @ A[k:k+1, k+1:]

# LU 分解
def lu_decompose(A, b):
    forward_elim(A, b)

# 前進代入で Lx=b を解く（L は単位下三角行列）
def forward_subst(L, b):
    n = L.shape[0]
    x = np.zeros(n)
    b = b.reshape(-1)
    for k in range(n):
        x[k] = b[k] - L[k, :k] @ x[:k]
    return x

# 後退代入で Ux=b を解く（U は上三角行列）
def backward_subst(U, b):
    n = U.shape[0]
    x = np.zeros(n)
    b = b.reshape(-1)
    for k in reversed(range(n)):
        x[k] = (b[k] - U[k, k+1:] @ x[k+1:]) / U[k, k]
    return x

# LU 分解法（部分ピボット選択付き）
def lu_solve(A, b):
    lu_decompose(A, b)
    y = forward_subst(A, b)
    x = backward_subst(A, y)
    return x
```

```{code-cell}
"""
LU 分解法の使用例
"""

A = np.array([[2., 4., 1.], [-1., -2., 2.], [4., 2., -3.]])
b = np.array([0., 10., 2.]).reshape(-1, 1)
x = lu_solve(A, b)
print(x)
```

```{note}
本来は LU 分解に $\boldsymbol{b}$ は不要であるが，上のプログラムでは関数 `lu_decompose()` の引数として $\boldsymbol{A}$ と $\boldsymbol{b}$ の両方を渡している。
そのうえで，前進消去において $\boldsymbol{A}$ と $\boldsymbol{b}$ を同時に更新している。
これは，ピボット操作の情報を格納することを避けるためである。

別の方法として，ピボット操作に相当する置換行列を記憶しておき $\boldsymbol{L} \boldsymbol{U} \boldsymbol{x} = \boldsymbol{P} \boldsymbol{b}$ を解く方法も考えられる。その際，$n$ 次正方行列である置換行列 $\boldsymbol{P}$ 全体を記憶することを避け，長さ $n$ の添字のリストとして置換の情報を保存することもできる。
たとえば Scipy の [scipy.linalg.lu_factor()](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lu_factor.html) はそのような実装の一例である。
```

## コレスキー法

コレスキー法は，方程式の係数行列 $\boldsymbol{A}$ が正定値対称行列であるときに用いることができ，LU 分解法の計算コストを半分くらいに抑えることができる。

なお，行列 $\boldsymbol{A}$ が正定値対称行列でない場合は転置行列 $\boldsymbol{A}^\top$ を方程式の両辺にかけて $\boldsymbol{A}^\top \boldsymbol{A} \boldsymbol{x} = \boldsymbol{A}^\top \boldsymbol{b}$ を考えればよい。行列 $\boldsymbol{A}$ が正則であれば，新たな係数行列 $\boldsymbol{A}^\top \boldsymbol{A}$ は正定値対称行列となる。

### 正定値行列

正方行列 $\boldsymbol{A}$ が **正定値**（positive definite）であるとは，零ベクトルでない任意のベクトル $\boldsymbol{x}$ に対して

$$
\boldsymbol{x}^\top \boldsymbol{A} \boldsymbol{x} > 0
$$

が成り立つことをいう[^real]。

[^real]: ただしベクトル $\boldsymbol{x}$ の要素はすべて実数とする。以下も同様。

````{prf:example}
:label: positive-definite-1

以下の $2$ 次正方行列 $\boldsymbol{A}$ は正定値である。

```{math}
  \boldsymbol{A}
  =
  \begin{bmatrix}
    1 & 2 \\
    2 & 5
  \end{bmatrix}
```

なぜなら，任意のベクトル $\boldsymbol{x}=[x,y]^\top$ に対して

```{math}
    \boldsymbol{x}^\top \boldsymbol{A} \boldsymbol{x}
  & =
    \begin{bmatrix} x & y \end{bmatrix}
    \begin{bmatrix}
      1 & 2 \\
      2 & 5
    \end{bmatrix}
    \begin{bmatrix} x \\ y \end{bmatrix}
    = x^2 + 4 xy + 5 y^2
    = \left( x + 2 y \right)^2 + y^2
    \geq 0
```

であり，$x \neq 0$ かつ $y \neq 0$ のときこれは正の値をとるからである。
````

````{prf:example}
:label: positive-definite-2

以下の $2$ 次正方行列 $\boldsymbol{A}$ は正定値ではない。

```{math}
  \boldsymbol{A}
  =
  \begin{bmatrix}
    1 & 2 \\
    3 & 4
  \end{bmatrix}
```

なぜなら，ベクトル $\boldsymbol{x}=[-5,2]^\top$ に対して

```{math}
    \boldsymbol{x}^\top \boldsymbol{A} \boldsymbol{x}
  & =
    \begin{bmatrix} -5 & 2 \end{bmatrix}
    \begin{bmatrix}
      1 & 2 \\
      3 & 4
    \end{bmatrix}
    \begin{bmatrix} -5 \\ 2 \end{bmatrix}
    = -9 < 0
```

が成り立つからである。
````

左辺は 2 つのベクトル $\boldsymbol{x}$ と $\boldsymbol{A} \boldsymbol{x}$ の内積であることに注意すると，行列 $\boldsymbol{A}$ が正定値であるとは，任意のベクトル $\boldsymbol{x}$ が $\boldsymbol{A} \boldsymbol{x}$ となす角が $90$ 度より小さくなる（一方のベクトルに垂直な平面で空間を2分割したとき，他方も同じ側を向いている）ということを意味する。

### コレスキー分解

コレスキー法の基本となるのは，コレスキー分解である。

```{prf:theorem} コレスキー分解（Cholesky decomposition）
正定値対称行列 $\boldsymbol{A}$ は下三角行列 $\boldsymbol{L}$ を用いて $\boldsymbol{A} = \boldsymbol{L} \boldsymbol{L}^\top$ のように分解できる。また，対角成分がすべて正であるような $\boldsymbol{L}$ は一意に定まる。
```

````{prf:example}
:label: cholesky

正定値対称行列 $\boldsymbol{A}$ を次のように定義する。

$$
\boldsymbol{A} =
  \begin{bmatrix}
    4 & 2 & 6  \\
    2 & 3 & 1 \\
    6 & 1 & 12
  \end{bmatrix}
$$

このとき，$\boldsymbol{A}$ のコレスキー分解は次のようになる。

$$
\boldsymbol{A} =
  \begin{bmatrix}
    2 & 0         & 0 \\
    1 & \sqrt{2}  & 0 \\
    3 & -\sqrt{2} & 1
  \end{bmatrix}
  \begin{bmatrix}
    2 & 1        & 3         \\
    0 & \sqrt{2} & -\sqrt{2} \\
    0 & 0        & 1
  \end{bmatrix}
$$
````

$3$ 次正方行列 $\boldsymbol{A} = [a_{ij}]$ のコレスキー分解を考える。
下三角行列 $\boldsymbol{L}$ を次のように定義する。

$$
L = \begin{bmatrix}
      l_{11} & 0      & 0 \\
      l_{21} & l_{22} & 0 \\
      l_{31} & l_{32} & l_{33}
    \end{bmatrix}
$$

コレスキー分解の式 $\boldsymbol{A} = \boldsymbol{L} \boldsymbol{L}^\top$ に代入すると，

$$
  & \boldsymbol{A}
    = \begin{bmatrix}
        l_{11} & 0      & 0 \\
        l_{21} & l_{22} & 0 \\
        l_{31} & l_{32} & l_{33}
      \end{bmatrix}
      \begin{bmatrix}
        l_{11} & l_{21} & l_{31} \\
        0      & l_{22} & l_{32} \\
        0      & 0      & l_{33}
      \end{bmatrix}
    = \begin{bmatrix}
        l_{11}^2      & l_{11} l_{21}                 & l_{11} l_{31} \\
        l_{11} l_{21} & l_{21}^2 + l_{22}^2           & l_{21} l_{31} + l_{22} l_{32} \\
        l_{11} l_{31} & l_{21} l_{31} + l_{22} l_{32} & l_{31}^2 + l_{32}^2 + l_{33}^2
      \end{bmatrix}
$$

これを解くと，

$$
  & l_{11} = \sqrt{a_{11}}, \quad
    l_{21} = \frac{a_{21}}{l_{11}}, \quad
    l_{31} = \frac{a_{31}}{l_{11}}, \\
  & l_{22} = \sqrt{a_{22} - l_{21}^2}, \quad
    l_{32} = \frac{a_{32} - l_{21} l_{31}}{l_{22}}, \\
  & l_{33} = \sqrt{a_{33} - l_{31}^2 - l_{32}^2}
$$

同様に，$3$ 次正方行列 $\boldsymbol{A} = [a_{ij}]$ のコレスキー分解は次のようになる。

$$
  & \boldsymbol{A}
    = \begin{bmatrix}
        l_{11}    & 0         & \cdots & 0           & 0 \\
        l_{21}    & l_{22}    & \cdots & 0           & 0 \\
        \vdots    & \vdots    & \ddots & \vdots      & \vdots \\
        l_{n-1,1} & l_{n-1,2} & \cdots & l_{n-1,n-1} & 0 \\
        l_{n1}    & l_{n2}    & \cdots & l_{n,n-1}   & l_{nn} \\
      \end{bmatrix}
      \begin{bmatrix}
        l_{11} & l_{21} & \cdots & l_{n-1,1}   & l_{n1} \\
        0      & l_{22} & \cdots & l_{n-1,2}   & l_{n2} \\
        \vdots & \vdots & \ddots & \vdots      & \vdots \\
        0      & 0      & \cdots & l_{n-1,n-1} & l_{n,n-1} \\
        0      & 0      & \cdots & 0           & l_{nn} \\
      \end{bmatrix}
$$

よって $i,j$ について以下が成り立つ。

$$
a_{ij} = \sum_{h=1}^{\min\{i,j\}} l_{ih} l_{jh}
$$ (aij-lih-ljh)

行列 $\boldsymbol{A}$ は対称行列であるから，$a_{ij}=a_{ji}$ が成り立つ。
したがって，$1 \leq j \leq i \leq n$ の範囲だけを考えればよい。

まず，$j=1$ のときは {eq}`aij-lih-ljh` より

$$
a_{i1} = l_{i1} l_{11}
$$

となる。$i=1,2,\cdots,n$ に対して，これを $l_{i1}$ についてそれぞれ解くと，

$$
  & l_{11} = \sqrt{a_{11}}
\\
  & l_{i1} = \frac{a_{i1}}{l_{11}} \quad i=2,3,\cdots,n
$$

同様にして $j=2,3,\cdots,n$ および $i=j,j+1,\cdots,n$ について順に $l_{ij}$ を求めることができる。

$$
  & l_{jj} = \sqrt{a_{jj} - \sum_{h=1}^{j-1} l_{jh}^2}
\\
  & l_{ij} = \frac{1}{l_{jj}} \left( a_{ij} - \sum_{h=1}^{j-1} l_{ih} l_{jh} \right) \quad i=j+1,j+2,\cdots,n
$$

Python によるコレスキー分解の実装例を示す。

```{code-cell}
"""
コレスキー分解の実装例
"""

import numpy as np

def cholesky(A):
    n = A.shape[0]
    L = np.zeros(A.shape)
    for j in range(n):
        L[j, j] = np.sqrt(A[j, j] - np.sum([L[j, h] ** 2 for h in range(j)]))
        for i in range(j + 1, n):
            L[i, j] = (A[i, j] - np.sum([L[i, h] * L[j, h] for h in range(j)])) / L[j, j]
    return L
```

```{code-cell}
"""
コレスキー分解 cholesky() の使用例
"""

# 行列 A を初期化する
A = np.array([[4, 2, 6], [2, 3, 1], [6, 1, 12]])
print('行列A:\n', A, '\n')

# コレスキー分解を実行する
L = cholesky(A)

print('行列Aのコレスキー分解 L:\n', L, '\n')

print('L と L^t の積は A に等しい:\n', L @ L.T)
```

コレスキー分解 $\boldsymbol{A} = \boldsymbol{L} \boldsymbol{L}^\top$ は LU 分解の一種であると考えることができる。したがって，LU 分解法と同様の手順で，コレスキー分解を用いて連立方程式を解くことができる。

```{prf:algorithm} コレスキー分解法
:label: alg:cholesky

**Inputs** $\boldsymbol{A}$: $n$次正方行列，$\boldsymbol{b}$: $n$次元ベクトル

**Outputs** 方程式 $\boldsymbol{A} \boldsymbol{x} = \boldsymbol{b}$ の 解 $\boldsymbol{x}$

1. $\boldsymbol{A}$ のコレスキー分解 $\boldsymbol{A} = \boldsymbol{L} \boldsymbol{L}^\top$ を求める。
2. 方程式 $\boldsymbol{L} \boldsymbol{y} = \boldsymbol{b}$ を $\boldsymbol{y}$ について解く。
3. 方程式 $\boldsymbol{L}^\top \boldsymbol{x} = \boldsymbol{y}$ を $\boldsymbol{x}$ について解く。
4. 解 $\boldsymbol{x}$ を返す。
```

LU 分解と比べると $\boldsymbol{U}$ を求める必要がなくなるため，計算コストが約半分に抑えられる。

コレスキー分解法の実装例は省略する（LU 分解法と同様である）。

### 修正コレスキー分解

コレスキー分解にはいくつかの問題が知られており，実際にはそのまま用いられることは少ない。
たとえば，コレスキー分解は平方根の計算をともなうが，この計算コストは比較的高い。

このような問題に対処するために用いられるのが次に紹介する **修正コレスキー分解**（modified Cholesky decomposition）である。

修正コレスキー分解は，行列 $\boldsymbol{A}$ を $\boldsymbol{L} \boldsymbol{D} \boldsymbol{L}^\top$ のように分解する。ただし，$\boldsymbol{L}$ は単位下三角行列[^unit-triangular]であり，$\boldsymbol{D}$ は対角行列である。

[^unit-triangular]: 単位下三角行列は，対角成分がすべて $1$ であるような下三角行列である。

````{prf:example}
:label: cholesky-mod

正定値対称行列 $\boldsymbol{A}$ を次のように定義する。

```{math}
\boldsymbol{A} =
  \begin{bmatrix}
    4 & 2 & 6  \\
    2 & 3 & 1 \\
    6 & 1 & 12
  \end{bmatrix}
```

このとき，$\boldsymbol{A}$ の修正コレスキー分解は次のようになる。

```{math}
\boldsymbol{A} =
  \begin{bmatrix}
    1   & 0  & 0 \\
    1/2 & 1  & 0 \\
    3/2 & -1 & 1
  \end{bmatrix}
  \begin{bmatrix}
    4 & 0 & 0 \\
    0 & 2 & 0 \\
    0 & 0 & 1
  \end{bmatrix}
  \begin{bmatrix}
    1 & 1/2 & 3/2 \\
    0 & 1   & -1 \\
    0 & 0   & 1
  \end{bmatrix}
```
````

コレスキー分解 $\boldsymbol{A}=\boldsymbol{L}\boldsymbol{L}^\top$ における下三角行列 $\boldsymbol{L}$ を

$$
\boldsymbol{L}=\boldsymbol{L}'\boldsymbol{D}'
$$

のように単位下三角行列と対角行列の積に分解できたとすると，

$$
\boldsymbol{A}
=\boldsymbol{L}\boldsymbol{L}^\top
=(\boldsymbol{L}'\boldsymbol{D}')(\boldsymbol{L}'\boldsymbol{D}')^\top
=\boldsymbol{L}'(\boldsymbol{D}')^2\boldsymbol{L}'^\top
=\boldsymbol{L}' \boldsymbol{D} \boldsymbol{L}'^\top
$$

となり，修正コレスキー分解が得られる。

行列 $\boldsymbol{L}$，$\boldsymbol{L}'$，$\boldsymbol{D}$，$\boldsymbol{D}'$ の $(i,j)$ 成分をそれぞれ $l_{ij}$，$l'_{ij}$，$d_{ij}$，$d'_{ij}$ とすると，

$$
  & d_{ij} = (d'_{ij})^2 = \begin{cases} (l_{jj})^2 & (i = j) \\ 0 & (i \neq j) \end{cases}
\\
  & l'_{ij} = \begin{cases} 0 & (i < j) \\ 1 & (i = j) \\ l_{ij} / d'_{jj} & (i > j) \end{cases}
$$

が成り立つ。したがって，$1 \leq j < i \leq n$ に対して次が成り立つ。

$$
  l'_{ij}
  &= \frac{1}{l_{jj} d'_{jj}} \left( a_{ij} - \sum_{h=1}^{j-1} l_{ih} l_{jh} \right)
\\
  &= \frac{1}{d_{jj}} \left( a_{ij} - \sum_{h=1}^{j-1} l'_{ih} l'_{jh} d_{hh} \right)
$$

ただし，$j=1,2,\cdots,n$ に対して

$$
  d_{jj}
  = (d'_{jj})^2
  = l_{jj}^2
  = a_{jj} - \sum_{h=1}^{j-1} l_{jh}^2
  = a_{jj} - \sum_{h=1}^{j-1} (l'_{jh})^2 d_{hh}
$$

である。

Python による修正コレスキー分解の実装例を示す。

```{code-cell}
"""
修正コレスキー分解の実装例
"""

import numpy as np

def cholesky_modified(A):
    n = A.shape[0]
    L = np.zeros(A.shape)
    D = np.zeros(n)
    for j in range(n):
        D[j] = A[j, j] - np.sum(L[j, :j] ** 2 * D[:j])
        L[j, j] = 1
        for i in range(j + 1, n):
            L[i, j] = (A[i, j] - np.sum(L[i, :j] * L[j, :j] * D[:j])) / D[j]
    return L, D
```

```{code-cell}
"""
修正コレスキー分解 cholesky_modified() の使用例
"""

# 行列 A を初期化する
A = np.array([[4, 2, 6], [2, 3, 1], [6, 1, 12]])
print('行列A:\n', A, '\n')

# 修正コレスキー分解を実行する
L, D = cholesky_modified(A)

print('行列Aの修正コレスキー分解 L, D:\n', L, '\n\n', np.diag(D), '\n')

print('L と D と L^t の積は A に等しい:\n', L @ np.diag(D) @ L.T)
```


```python
for i in range(j + 1, n):
    L[i, j] = (A[i, j] - np.sum(L[i, :j] * L[j, :j] * D[:j])) / D[j]
```

は

```python
L[j+1:n, j] = (A[j+1:n, j] - np.sum(L[j+1:n, :j] * L[j, :j] * D[:j]), axis=1) / D[j]
```

と書いてもよい。

修正コレスキー分解を用いて連立方程式 $\boldsymbol{A}\boldsymbol{x}=\boldsymbol{b}$ を解くアルゴリズムは次のようになる。

```{prf:algorithm} 修正コレスキー分解法
:label: alg:cholesky-mod

**Inputs** $\boldsymbol{A}$: $n$次正方行列，$\boldsymbol{b}$: $n$次元ベクトル

**Outputs** 方程式 $\boldsymbol{A} \boldsymbol{x} = \boldsymbol{b}$ の 解 $\boldsymbol{x}$

1. $\boldsymbol{A}$ の修正コレスキー分解 $\boldsymbol{A} = \boldsymbol{L} \boldsymbol{D} \boldsymbol{L}^\top$ を求める。
2. 方程式 $\boldsymbol{L} \boldsymbol{D} \boldsymbol{y} = \boldsymbol{b}$ を $\boldsymbol{y}$ について解く。
3. 方程式 $\boldsymbol{L}^\top \boldsymbol{x} = \boldsymbol{y}$ を $\boldsymbol{x}$ について解く。
4. 解 $\boldsymbol{x}$ を返す。
```

修正コレスキー分解法の実装例は省略する（LU 分解法やコレスキー分解法と同様である）。
