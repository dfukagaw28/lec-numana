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

# 連立1次方程式の数値解法(1)

連立1次方程式（線型方程式系; system of linear equations; SLE）は，日本の現行の学習指導要領（[中学校学習指導要領（平成 29 年告示）](https://www.mext.go.jp/content/20230120-mxt_kyoiku02-100002604_02.pdf)）によれば中学2年生で学ぶ基礎的な内容である[^youryou]。
一方で応用範囲は広く，多くの問題が連立1次方程式に帰着されるという意味において，重要な問題である。

[^youryou]: 正確には中学2年生で学ぶのは連立二元一次方程式である。

小さな連立1次方程式は加減法や代入法のような単純な方法で解くことができた。
変数の数，式の数が増えたとしても，解法プログラムを作成するのは容易であるように思える[^cramer]。
実際には，数値誤差を考えると工夫が必要である。

[^cramer]: 連立1次方程式の解を陽に表すクラメルの公式（Cramer's rule）が知られているが，実際に解を求める方法としては適さない。なぜならば，クラメルの公式は行列式の計算を $n+1$ 回繰り返す必要があり，行列式の計算はガウス消去法と同程度の計算コストを要するからである。ガウス消去法を1回だけ実行する方が効率的である。

連立1次方程式は，大きく分けて直接解法と反復解法の 2 通りである。
**直接解法**（direct method）には，ガウスの消去法やコレスキー法などがある。
**反復解法**（iterative method）には，ヤコビ法，ガウス・ザイデル法，逐次緩和法（SOR法），共役勾配法（CG法）などがある。

今回は直接解法の代表的アルゴリズムであるガウスの消去法について解説する。

## ガウスの消去法

ガウスの消去法は掃き出し法とも呼ばれる。
連立1次方程式を解くための代表的なアルゴリズムである。

まず，ガウスの消去法のアルゴリズムの基本的な流れについて述べる。
この基本のアルゴリズムは，後で述べる工夫を加えないため， **単純ガウス消去法** と呼ばれることがある。

### 例1（単純ガウス消去法）

以下の連立1次方程式を考える。

$$
\begin{alignat}{4}
 3 x & +{} &   y & +{} & 2 z & ={} & 13 \\
 5 x & +{} &   y & +{} & 3 z & ={} & 20 \\
 4 x & +{} & 2 y & +{} &   z & ={} & 13
\end{alignat}
$$ (simeq1)

一般に，連立1次方程式に以下の操作を行っても，方程式の解は変わらない。

- 2つの行を入れ替える
- ある行を $\alpha$ 倍する（$\alpha\ne 0$）
- ある行に他の行を足す

つまり，第 $i$ 行目に非零な実数 $\alpha$ をかけたものを第 $j$ 行目に足しても
方程式の解は変わらない。

以下，この性質を繰り返し利用して連立方程式の項を階段状に消去する。
この操作を **前進消去** （forward elimination）とよぶ。

上の連立方程式において，
まず $x$ の項を第2行と第3行目から消去する。
第 1 行の $-\frac{5}{3}$ 倍を
第 2 行に足し，
さらに，
第 1 行の $-\frac{4}{3}$ 倍を
第 3 行に足すと
以下のようになる。

$$
\begin{alignat}{4}
 3 x & +{} &             y & +{} &           2 z & ={} & 13 \\
     & -{} & \frac{2}{3} y & -{} & \frac{1}{3} z & ={} & -\frac{5}{3} \\
     &     & \frac{2}{3} y & -{} & \frac{5}{3} z & ={} & -\frac{13}{3}
\end{alignat}
$$

以下では，行列表示を用いる。
最初の連立方程式 {eq}`simeq1` は以下のように書ける。

$$
  \begin{bmatrix}
    3 & 1 & 2 \\
    5 & 1 & 3 \\
    4 & 2 & 1
  \end{bmatrix}
  \begin{bmatrix}
    x \\ y \\ z
  \end{bmatrix}
  =
  \begin{bmatrix}
    13 \\ 20 \\ 13
  \end{bmatrix}
$$

$x$の項を消去した連立方程式は以下のように書ける。

$$
  \begin{bmatrix}
    3 & 1             & 2 \\
    0 & -\dfrac{2}{3} & -\dfrac{1}{3} \\
    0 &  \dfrac{2}{3} & -\dfrac{5}{3}
  \end{bmatrix}
  \begin{bmatrix}
    x \\ y \\ z
  \end{bmatrix}
  =
  \begin{bmatrix}
    13 \\ -\dfrac{5}{3} \\ -\dfrac{13}{3}
  \end{bmatrix}
$$

引き続き，第3行の $y$ を消去するために，第2行の $1$ 倍を第3行に足す。

```{math}
:label: simple-gauss-elim-fwd
  \begin{bmatrix}
    3 & 1             & 2 \\
    0 & -\dfrac{2}{3} & -\dfrac{1}{3} \\
    0 & 0             & -2
  \end{bmatrix}
  \begin{bmatrix}
    x \\ y \\ z
  \end{bmatrix}
  =
  \begin{bmatrix}
    13 \\ -\dfrac{5}{3} \\ -6
  \end{bmatrix}
```

左辺の行列について，階段状に各成分が消去されて $0$ となり，右上三角行列（非零な要素が右上に位置し，左下には$0$だけが並ぶ）となった点に注目してほしい。
このような形式の行列は **階段形**（echelon form）あるいは **行階段形**（row echelon form）であるという。

行階段形において，最後の行に注目すると，これは 1 変数の 1 次方程式であり，容易に解を求めることができる。
下から2番目の行は 2 変数の 1 次方程式となるが，2 個の変数のうち 1 個にすでに解を代入すると，もう一方の変数について解を求めることができる。

このように，後ろから順にさかのぼって解を求める。
求めた解を次々に代入してゆくところがポイントである。
この操作を **後退代入** （back substitution）とよぶ。

まず式 {eq}`simple-gauss-elim-fwd` の連立方程式の第3行目は

$$
  -2 z = -6
$$

であるから，$z=3$ が得られる。

次に連立方程式 {eq}`simple-gauss-elim-fwd` の第2行目は

$$
  -\frac{2}{3} y - \frac{1}{3} z = -\frac{5}{3}
$$

であるが，これを $z=3$ を代入して $y$ について解くと以下のように $y$ の値が得られる。

```{math}
  &
  -\frac{2}{3} y - \frac{1}{3} \cdot 3 = -\frac{5}{3}
\\&
  y
  = \left( -\frac{5}{3} + \frac{1}{3} \cdot 3 \right) \cdot \left( -\frac{2}{3} \right)^{-1}
  = 1
```

同様に，第1行目は

$$
  3 x + y + 2 z = 13
$$

である。これに $y,x$ の値を代入したうえで $x$ について解けばよい。

```{math}
  &
  3 x_1 + 1 + 2 \cdot 3 = 13
\\&
  x
  = \frac{13 - 1 - 2 \cdot 3}{3}
  = 2
```

以上で連立方程式の解 $(x,y,z)=(2,1,3)$ が求められた。

### 単純ガウス消去法のアルゴリズム

上の例で用いた手順を一般化しよう。

変数の数と式の数をともに $n$ とすると，連立1次方程式は以下のように書ける。

$$
  \begin{bmatrix}
    a_{11} & a_{12} & \cdots & a_{1n} \\
    a_{21} & a_{22} & \cdots & a_{2n} \\
    \vdots & \vdots & \ddots & \vdots \\
    a_{n1} & a_{n2} & \cdots & a_{nn}
  \end{bmatrix}
  \begin{bmatrix}
    x_1 \\ x_2 \\ \vdots \\ x_n
  \end{bmatrix}
  =
  \begin{bmatrix}
    b_1 \\ b_2 \\ \vdots \\ b_n
  \end{bmatrix}
$$

この式は，行列 $\boldsymbol{A}$，ベクトル $\boldsymbol{b}$，$\boldsymbol{x}$ によって以下のように書くことができる。ただし，$\boldsymbol{A}$，$\boldsymbol{b}$ は定数であり，$\boldsymbol{x}$ は変数である。

$$
\boldsymbol{A} \boldsymbol{x} = \boldsymbol{b}
$$

この定式化のもとで，単純ガウス消去法の手順は以下のようになる。

- 前進消去
  1. 第 1 行を用いて $a_{21},a_{31},\ldots,a_{n1}$ を消去する
  2. 第 2 行を用いて $a_{32},a_{42},\ldots,a_{n2}$ を消去する
  3. (中略)
  4. 第 $n-1$ 行を用いて $a_{n,n-1}$ を消去する
- 後退代入
  1. 第 $n$ 行を $x_n$ について解く
  2. 第 $n-1$ 行に解 $x_n$ を代入したうえで $x_{n-1}$ について解く
  3. 第 $n-2$ 行に解 $x_n,x_{n-1}$ を代入したうえで $x_{n-2}$ について解く
  4. （中略）
  5. 第 $1$ 行に解 $x_n,x_{n-1},\ldots,x_2$ を代入したうえで $x_1$ について解く

この手順を擬似コードで表すと次の {prf:ref}`alg:gauss-elim-1` のようになる。

```{prf:algorithm} 単純ガウス消去法
:label: alg:gauss-elim-1

**Inputs** $A$: $n$次正方行列，$\boldsymbol{b}$: $n$次元ベクトル

**Outputs** 方程式 $A \boldsymbol{x} = \boldsymbol{b}$ の 解 $\boldsymbol{x}$

1. for $k = 1 , 2, \ldots, n-1$ do
   1. for $i = k+1 , k+2, \ldots, n$ do
      1. $\alpha \leftarrow - a_{ik} / a_{kk}$
      2. for $j = k, k+1, \ldots, n$ do
         1. $a_{ij} \leftarrow a_{ij} + \alpha \cdot a_{kj}$
      3. $b_i \leftarrow b_{i} + \alpha \cdot b_{k}$
2. for $k = n , n-1 , \ldots , 1$ do
   1. $x_{k} \leftarrow \frac{1}{a_{kk}} \left\{ b_{k}-\sum_{j=k+1}^{n}{a_{kj} \cdot x_{j}} \right\}$
3. return $\boldsymbol{x}$
```

この {prf:ref}`alg:gauss-elim-1` を Python で素直に書くと，以下のようになる。
配列の添字が $0$ から開始する点に注意。

```{code-cell}
"""
単純ガウス消去法（list版）
"""

# 前進消去
def forward_elim_1(A, b):
    n = len(A)
    for k in range(n-1):
        for i in range(k + 1, n):
            alpha = - A[i][k] / A[k][k]
            for j in range(k, n):
                A[i][j] += alpha * A[k][j]
            b[i] += alpha * b[k]

# 後退代入
def backward_subst_1(A, b):
    n = len(A)
    x = [0] * n
    for k in reversed(range(n)):
        s = b[k]
        for j in range(k + 1, n):
            s -= A[k][j] * x[j]
        x[k] = s / A[k][k]
    return x

# 単純ガウス消去法
def gauss_elim_1(A, b):
    forward_elim_1(A, b)
    x = backward_subst_1(A, b)
    return x

# 実行例
A = [[3, 1, 2], [5, 1, 3], [4, 2, 1]]
b = [13, 20, 13]
x = gauss_elim_1(A, b)
print(x)
```

NumPy を用いて書き直すと，たとえば以下のように書くことができる。

```{code-cell}
"""
単純ガウス消去法（NumPy版）
"""

import numpy as np

# 前進消去（NumPy版）
def forward_elim_2(A, b):
    n = A.shape[0]
    for k in range(n-1):
        alpha = - A[k+1:, k:k+1] / A[k, k]
        A[k+1:, k:] += alpha @ A[k:k+1, k:]
        b[k+1:] += alpha * b[k]

# 後退代入（NumPy版）
def backward_subst_2(A, b):
    n = A.shape[0]
    x = np.zeros(n)
    for k in reversed(range(n)):
        s = b[k, 0] - A[k, k+1:] @ x[k+1:]
        x[k] = s / A[k, k]
    return x

# 単純ガウス消去法
def gauss_elim_2(A, b):
    forward_elim_2(A, b)
    x = backward_subst_2(A, b)
    return x

# 実行例
A = np.array([[3., 1., 2.], [5., 1., 3.], [4., 2., 1.]])
b = np.array([13., 20., 13.]).reshape(-1, 1)
x = gauss_elim_2(A, b)
print(x)
```

### 例2（ピボット選択）

実は，上で述べた単純ガウス消去法のアルゴリズムにはバグがある。
例えば，以下のようなパラメータを持つ連立方程式を扱うとゼロ除算が発生して異常な結果が出力される。

```{code-cell}
:tags: [raises-exception]
# 単純ガウス消去法（list版）の実行例
A = [[2, 4, 1], [-1, -2, 2], [4, 2, -3]]
b = [0, 10, 2]
x = gauss_elim_1(A, b)
print(x)
```

```{code-cell}
# 単純ガウス消去法（NumPy版）の実行例
A = np.array([[2, 4, 1], [-1, -2, 2], [4, 2, -3]], dtype=float)
b = np.array([0, 10, 2], dtype=float).reshape(-1, 1)
x = gauss_elim_2(A, b)
print(x)
```

この例において起きている現象は，以下のようなものである。

まず，以下のような連立方程式がある。

$$
\begin{alignat}{4}
 2 x & +{} & 4 y & +{} &   z & ={} & 0 \\
 - x & -{} & 2 y & +{} & 2 z & ={} & 10 \\
 4 x & +{} & 2 y & -{} & 3 z & ={} & 2
\end{alignat}
$$

第1式を定数倍したものを第2，3式に加えて変数 $x$ を消去する。

$$
\begin{alignat}{4}
 2 x & +{} & 4 y & +{} & z             & ={} & 0 \\
     &     &     & +{} & \frac{5}{2} z & ={} & 10 \\
     & -{} & 6 y & -{} & 5 z           & ={} & 2
\end{alignat}
$$

次に，第2式を定数倍して第3式の $y$ を消去することを考えたい。
ところが，第2式の $y$ の項は存在しない。
前のステップで $x$ を消去した際に $y$ の係数が偶然にも $0$ になってしまったためである。

単純ガウス消去法が失敗したからと言って，この連立方程式が解けないわけではない。
たとえば，次のように第2式と第3式を入れ換えてしまえば解くことができる。

$$
\begin{alignat}{4}
 2 x & +{} & 4 y & +{} & z             & ={} & 0 \\
     & -{} & 6 y & -{} & 5 z           & ={} & 2 \\
     &     &     & +{} & \frac{5}{2} z & ={} & 10
\end{alignat}
$$

係数を確認すると，既に前進消去が終わった直後と同じ状態である。

よって，後退代入を行うと，

$$
 z &= 10 \cdot \left( \frac{5}{2} \right)^{-1} = 4 \\
 y &= \left( 2 + 5 z \right) \cdot \left( -6 \right)^{-1}
    = \frac{22}{-6} = - \frac{11}{3} \\
 x &= \left( 0 - 4 y  - z \right) \cdot \left( 2 \right)^{-1}
    = \left( 0 + 4 \cdot \frac{11}{3}  - 4 \right) \cdot \left( 2 \right)^{-1}
    = \frac{44-12}{6} = \frac{16}{3}
$$

以上により，連立方程式の解 $(x,y,z)=(\frac{16}{3},-\frac{11}{3},4)$ を得る。

このように，ガウスの消去法において，前進消去の途中で **軸**（pivot）となる項の係数が $0$ となって項が消えてしまったとき，他の式と交換することによってアルゴリズムを継続するという工夫のことを **部分ピボット選択** という。

一般に，部分ピボット選択付きガウス消去法のことを単に **ガウス消去法** とよぶ。

上に挙げたものと同じ例の行列表示バージョンを用いて，部分ピボット選択付きガウス消去法のアルゴリズムを確認しよう。

まず，初期状態は以下のようになる。

$$
  \begin{bmatrix}
    2 & 4 & 1 \\
    -1 & -2 & 2 \\
    4 & 2 & -3
  \end{bmatrix}
  \begin{bmatrix}
    x_1 \\ x_2 \\ x_3
  \end{bmatrix}
  =
  \begin{bmatrix}
    0 \\ 10 \\ 2
  \end{bmatrix}
$$

さきほどは $a_{11}$ を軸にして前進消去の最初のステップを実行したが，ここでも式の交換を行うことができる。

ここでは，第 1 列すなわち $|a_{11}|$，$|a_{21}|$，$|a_{31}|$ を比較し，絶対値が最大となる $a_{31}=4$ を選ぶ。すなわち，第 3 行と第 1 行を交換し，$a_{11}$ の位置に $4$ が来るようにする。

$$
  \begin{bmatrix}
    \fbox{4} & 2 & -3 \\
    -1 & -2 & 2 \\
    2 & 4 & 1
  \end{bmatrix}
  \begin{bmatrix}
    x_1 \\ x_2 \\ x_3
  \end{bmatrix}
  =
  \begin{bmatrix}
    2 \\ 10 \\ 0
  \end{bmatrix}
$$

交換を行ったうえで，第1式を定数倍して第2式以降に足し合わせることで係数 $a_{21}$，$a_{31}$ を消去する。

$$
  \begin{bmatrix}
    \fbox{4} & 2 & -3 \\
    0 & -3/2 & 5/4 \\
    0 & 3 & 5/2
  \end{bmatrix}
  \begin{bmatrix}
    x_1 \\ x_2 \\ x_3
  \end{bmatrix}
  =
  \begin{bmatrix}
    2 \\ 21/2 \\ -1
  \end{bmatrix}
$$

次に，2回目のピボットを選択する。
$|a_{22}|=3/2$ と $|a_{32}|=3$ を比較して絶対値が大きい $|a_{32}|=3$ を選択する。
つまり，第2式と第3式を交換する。

$$
  \begin{bmatrix}
    4 & 2 & -3 \\
    0 & \fbox{3} & 5/2 \\
    0 & -3/2 & 5/4
  \end{bmatrix}
  \begin{bmatrix}
    x_1 \\ x_2 \\ x_3
  \end{bmatrix}
  =
  \begin{bmatrix}
    2 \\ -1 \\ 21/2
  \end{bmatrix}
$$

交換を行ったうえで，第2式を定数倍して第3式に足し合わせることで係数 $a_{32}$ を消去する。

$$
  \begin{bmatrix}
    4 & 2 & -3 \\
    0 & \fbox{3} & 5/2 \\
    0 & 0 & 5/2
  \end{bmatrix}
  \begin{bmatrix}
    x_1 \\ x_2 \\ x_3
  \end{bmatrix}
  =
  \begin{bmatrix}
    2 \\ -1 \\ 10
  \end{bmatrix}
$$

あとは後退代入によって以下のように解が得られる。

$$
 x_3 &= 10 \cdot \left( \frac{5}{2} \right)^{-1} = 4 \\
 x_2 &= \left( -1 - \frac{5}{2} x_3 \right) \cdot \left( 3 \right)^{-1}
      = - \frac{11}{3} \\
 x_1 &= \left( 2 - 2 x_2  + 3 x_3 \right) \cdot \left( 4 \right)^{-1}
      = \left( 2 + \frac{22}{3}  + 12 \right) \cdot \left( 4 \right)^{-1}
    = \frac{22+42}{12} = \frac{16}{3}
$$

### 部分ピボット選択付きガウス消去法のアルゴリズム

部分ピボット選択付きガウス消去法の手順は以下のようになる。

- 前進消去
  1. 第 $1$ 列に着目し，係数 $|a_{h1}|$ が最大になるような添え字 $h=1,2,\ldots,n$ を見つけ，第 $1$ 行と第 $h$ 行とを入れ替える
  2. 第 $1$ 行を用いて $a_{21},a_{31},\ldots,a_{n1}$ を消去する
  3. 第 $2$ 列に着目し，係数 $|a_{h2}|$ が最大になるような添え字 $h=2,3,\ldots,n$ を見つけ，第 $2$ 行と第 $h$ 行とを入れ替える
  4. 第 $2$ 行を用いて $a_{32},a_{42},\ldots,a_{n2}$ を消去する
  5. （中略）
  6. 第 $n-1$ 列に着目し，係数 $|a_{h,n-1}|$ が最大になるような添え字 $h=n-1,n$ を見つけ，第 $n-1$ 行と第 $h$ 行とを入れ替える
  7. 第 $n-1$ 行を用いて $a_{n,n-1}$ を消去する
- 後退代入
  1. 第 $n$ 行を $x_n$ について解く
  2. 第 $n-1$ 行に解 $x_n$ を代入したうえで $x_{n-1}$ について解く
  3. 第 $n-2$ 行に解 $x_n,x_{n-1}$ を代入したうえで $x_{n-2}$ について解く
  4. （中略）
  5. 第 $1$ 行に解 $x_n,x_{n-1},\ldots,x_2$ を代入したうえで $x_1$ について解く

この手順を擬似コードで表すと次の {prf:ref}`alg:gauss-elim-2` のようになる。

```{prf:algorithm} 部分ピボット選択付きガウス消去法
:label: alg:gauss-elim-2

**Inputs** $A$: $n$次正方行列，$\boldsymbol{b}$: $n$次元ベクトル

**Outputs** 方程式 $A \boldsymbol{x} = \boldsymbol{b}$ の 解 $\boldsymbol{x}$

1. for $k = 1 , 2, \ldots, n-1$ do
   1. $h^{*} \leftarrow \mathop{\rm argmax}\limits_{h=k,k+1,\ldots,n} |a_{hk}|$
   2. for $j = k, k+1, \ldots, n$ do
      1.  $a_{k,j}$ と $a_{h^*,j}$ を入れ替える
   3. $b_{k}$ と $b_{h^*}$ を入れ替える
   4. for $i = k+1, k+2, \ldots, n$ do
      1. $\alpha \leftarrow - a_{ik} / a_{kk}$
      2. for $j = k, k+1, \ldots, n$ do
         1. $a_{ij} \leftarrow a_{ij} + \alpha \cdot a_{kj}$
      3. $b_i \leftarrow b_{i} + \alpha \cdot b_{k}$
2. for $k = n , n-1 , \ldots , 1$ do
   1. $x_{k} \leftarrow \frac{1}{a_{kk}} \left\{ b_{k}-\sum_{j=k+1}^{n}{a_{kj} \cdot x_{j}} \right\}$
3. return $\boldsymbol{x}$
```

この {prf:ref}`alg:gauss-elim-2` を Python で以下のように書くことができる。

```{code-cell}
"""
部分ピボット選択付きガウス消去法（NumPy版）
"""

import numpy as np

# 部分ピボット選択
def pivot(A, b, k):
    h = k + np.abs(A[k:, k]).argmax()
    if h != k:
        A[[k, h], k:] = A[[h, k], k:]
        b[[k, h]] = b[[h, k]]

# 前進消去
def forward_elim_3(A, b):
    n = A.shape[0]
    for k in range(n-1):
        pivot(A, b, k)  # この行を追加
        alpha = - A[k+1:, k:k+1] / A[k, k]
        A[k+1:, k:] += alpha @ A[k:k+1, k:]
        b[k+1:] += alpha * b[k]

# 後退代入
def backward_subst_3(A, b):
    n = A.shape[0]
    x = np.zeros(n)
    for k in reversed(range(n)):
        s = b[k, 0] - A[k, k+1:] @ x[k+1:]
        x[k] = s / A[k, k]
    return x

# ガウス消去法（部分ピボット選択付き）
def gauss_elim_3(A, b):
    forward_elim_3(A, b)
    x = backward_subst_3(A, b)
    return x
```

```{code-cell}
# 部分ピボット選択付きガウス消去法の実行例
A = np.array([[2., 4., 1.], [-1., -2., 2.], [4., 2., -3.]])
b = np.array([0., 10., 2.]).reshape(-1, 1)
x = gauss_elim_3(A, b)
print(x)
```

理論上は，この部分ピボット選択付きガウス消去法は，行列 $\boldsymbol{A}$ が正則[^regular]であれば必ず解を求めることができる。
逆に，部分ピボット選択付きガウス消去法によって $\boldsymbol{A} \boldsymbol{x} = \boldsymbol{b}$ の解を求まるのは，$\boldsymbol{A}$ が正則である場合に限る。
ただし，数値誤差の影響で問題が生じる可能性がある点には注意を要する。

[^regular]: $n$ 次正方行列 $\boldsymbol{A}$ が **正則**（regular）であるとは，行列 $\boldsymbol{A}$ が逆行列を持つことをいう。さらに，行列が正則であることは，(1) 行列式が $0$ でないこと，(2) 行列のランクが $n$ であること，(3) 行列のカーネル（核）が零ベクトル $\boldsymbol{0}$ に限ること，(4) すべての固有値が $0$ でないこと，とそれぞれ同値である。

### スケーリング

連立方程式の各式の両辺にゼロでない実数をかけても解は変わらない。
ただし数値計算においては誤差の影響で解が変わる場合がある。
このことを利用して，数値誤差を回避する方法が知られている。

たとえば，十分に大きな正の数 $M \gg 1$ について，以下の方程式を考えよう。

$$
\begin{align}
  \begin{bmatrix}
    1 & M \\
    1 & 1
  \end{bmatrix}
  \begin{bmatrix} x_1 \\ x_2 \end{bmatrix}
  =
  \begin{bmatrix} M - 1 \\ 1 \end{bmatrix}
\end{align}
$$ (sle-scaling-2)

$\delta=\frac{1}{M-1}$ とおくと，真の解は $(x_1,x_2)=(\delta,1-\delta)$ である。

この方程式を部分ピボット選択付きガウス消去法で解いてみよう。

前進消去を実行すると次のようになる。

$$
\begin{align}
  \begin{bmatrix}
    1 & M \\
    0 & -(M - 1)
  \end{bmatrix}
  \begin{bmatrix} x_1 \\ x_2 \end{bmatrix}
  =
  \begin{bmatrix} M - 1 \\ -(M - 2) \end{bmatrix}
\end{align}
$$ (sle-scaling-3)

この場合，後退代入によって以下の解を得る。

$$
\begin{align}
  x_2
    &= \frac{M-2}{M-1}
     = 1 - \frac{1}{M-1}
     = 1 - \delta
\\x_1
    &= (M - 1) - M x_2
\\
    &= (M - 1) - M (1 - \delta)
\\
    &= (M - 1) - (M - 1 - \delta)
\\
    &\approx (M - 1) - (M - 1) \quad (\because M \gg 1 \gg \delta)
\\
    &= 0
\end{align}
$$

よって，近似解は $(\tilde{x}_1,\tilde{x}_2)=(0,1-\delta)$ となる。

真の解 $(x_1,x_2)=(\delta,1-\delta)$ と比較すると，十分に近い値と言えなくもないが，微小量 $\delta$ の情報が完全に失われるのは問題である。
たとえば，何らかの確率分布の推定をするときに $0$ と $10^{-10}$ では大きな違いをもたらす場合がある。
大きな方程式の場合は，誤差が 他の解にさらなる影響を与えることも考えなければならない。

対策として，各方程式の両辺を定数で割って係数の大きさをある程度揃えることによって計算誤差を抑える方法が知られている。この操作を **スケーリング**（scaling）という。

具体的には，方程式の各行を，各係数 $a_{ij}$ の絶対値 $|a_{ij}|$ の最大値が $1$ となるように揃える。
そのためには，第 $i$ 行の両辺を $s_i=\max_j |a_{ij}|$ で割ればよい。

上の例でスケーリングを適用する。
{eq}`sle-scaling-2` において，第 $1$ 式を $s_1=\max_j |a_{1j}|=\max\{1,M\}=M$ で割る。
第 $2$ 式では，$s_2=\max_j |a_{2j}|=\max\{1,1\}=1$となり，両辺を $s_2=1$ で割っても変化しない。
スケーリングの結果，方程式は以下のように変形される。
ただし，$\varepsilon=\frac{1}{M}$ とおく。

$$
\begin{align}
  \begin{bmatrix}
    \varepsilon & 1 \\
    1 & 1
  \end{bmatrix}
  \begin{bmatrix} x_1 \\ x_2 \end{bmatrix}
  =
  \begin{bmatrix} 1 - \varepsilon \\ 1 \end{bmatrix}
\end{align}
$$ (sle-scaling-4)

ピボット選択のために第 $1$ 行と第 $2$ 行を入れ換える。

$$
\begin{align}
  \begin{bmatrix}
    \fbox{1} & 1 \\
    \varepsilon & 1
  \end{bmatrix}
  \begin{bmatrix} x_1 \\ x_2 \end{bmatrix}
  =
  \begin{bmatrix} 1 \\ 1 - \varepsilon \end{bmatrix}
\end{align}
$$

次に，前進消去を行う。第 $1$ 行を $-\varepsilon$ 倍したものを第 $2$ 行に加えることによって以下の式を得る。

$$
\begin{align}
  \begin{bmatrix}
    \fbox{1} & 1 \\
    0 & 1 - \varepsilon
  \end{bmatrix}
  \begin{bmatrix} x_1 \\ x_2 \end{bmatrix}
  =
  \begin{bmatrix} 1 \\ 1 - 2 \varepsilon \end{bmatrix}
\end{align}
$$

ここから，後退代入によって，以下の解を得る。

$$
\begin{align}
  x_2
    &= \frac{1 - 2 \varepsilon}{1 - \varepsilon}
     = \frac{M - 2}{M - 1}
     = 1 - \delta
\\x_1
    &= 1 - x_2
     = 1 - (1 - \delta)
     = \delta
\end{align}
$$

これは上で述べた真の解に等しい[^fn-scaling]。

[^fn-scaling]: 幸いにも，今回は真の解を得ることができたが，スケーリングによって必ず真の解が得られるわけではない。

実際に，スケーリングの効果を Python プログラムを用いて確かめてみよう。

```{code-cell}
"""
スケーリングの有無による結果を比較する
"""

M = 1e10
eps = 1 / M
delta = 1 / (M - 1)
print('真の解:', (delta, (1-delta)))

A = np.array([[1., M], [1., 1.]])
b = np.array([M - 1, 1]).reshape(-1, 1)
x = gauss_elim_3(A, b)
print('スケーリングなし:', x)

A = np.array([[eps, 1.], [1., 1.]])
b = np.array([1 - eps, 1]).reshape(-1, 1)
x = gauss_elim_3(A, b)
print('スケーリングあり:', x)
```

この例では，スケーリングを行うことによってたしかに精度が向上していることが確認できた。

スケーリングを組み込んだガウス消去法は以下のようになる。

```{code-cell}
"""
部分ピボット選択・スケーリング付きガウス消去法（NumPy版）
"""

import numpy as np

# スケーリング
def scaling(A, b, k):
    v = np.abs(A[k:, k:]).max(axis=1).reshape(-1, 1)
    A[k:, k:] /= v
    b[k:] /= v

# 部分ピボット選択
def pivot(A, b, k):
    h = k + np.abs(A[k:, k]).argmax()
    if h != k:
        A[[k, h], k:] = A[[h, k], k:]
        b[[k, h]] = b[[h, k]]

# 前進消去
def forward_elim_4(A, b):
    n = A.shape[0]
    for k in range(n-1):
        scaling(A, b, k)  # この行を追加
        pivot(A, b, k)
        alpha = - A[k+1:, k:k+1] / A[k, k]
        A[k+1:, k:] += alpha @ A[k:k+1, k:]
        b[k+1:] += alpha * b[k]

# 後退代入
def backward_subst_4(A, b):
    n = A.shape[0]
    x = np.zeros(n)
    for k in reversed(range(n)):
        s = b[k, 0] - A[k, k+1:] @ x[k+1:]
        x[k] = s / A[k, k]
    return x

def gauss_elim_4(A, b):
    forward_elim_4(A, b)
    x = backward_subst_4(A, b)
    return x
```

```{code-cell}
# 部分ピボット選択・スケーリング付きガウス消去法の実行例
A = np.array([[1., 1/eps], [1., 1.]])
b = np.array([1/eps - 1, 1]).reshape(-1, 1)
x = gauss_elim_4(A, b)
print(x)
```
