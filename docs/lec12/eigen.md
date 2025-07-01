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

$\newcommand\bm[1]{\boldsymbol #1}$
$\newcommand\setC{\mathbb{C}}$
$\newcommand\naiseki[2]{\left\langle #1, #2\right\rangle}$
$\DeclareMathOperator{\sgn}{sgn}$
$\DeclareMathOperator{\ker}{ker}$

# 固有値問題

線形代数，特に行列の扱いは，工学の分野でも基礎的かつ重要な意味を持つ。
この講義でも，たびたび行列を用いた式が登場した。
機械学習においても線形代数が欠かせない。
行列は方程式などのシステム（系）を表現するため，行列の構造を知ることによってそのシステムの特徴・振る舞いが理解しやすくなる。

行列の構造を調べるには，固有値と固有ベクトルが重要な役割を果たす。
たとえば，データ分析でよく用いられる主成分分析（PCA）においては，多次元データを低次元に縮約するために固有値と固有ベクトルを利用する。

今回は，与えられた行列の固有値と固有ベクトルを数値的に求める方法について学ぶ。
固有値と固有ベクトルの定義や性質については他の講義で履修済みであるものと想定しているが，念のため，おさらいから始める。

## おさらい: 固有値，固有ベクトル

$n$ 次正方行列 $\bm{A}$ およびスカラー $\lambda \in \setC$ に対して，

$$
\bm{A} \bm{x} = \lambda \bm{x} \quad \text{かつ} \quad \bm{x} \ne \bm{0}
$$

を満たすベクトル $\bm{x} \in \setC^n$ が存在するとき，$\lambda$ を $\bm{A}$ の **固有値**（eigenvalue）といい，$\bm{x}$ を固有値 $\lambda$ に対する $\bm{A}$ の **固有ベクトル**（eigenvector）という。

たとえば，

$$
\bm{A} = \begin{bmatrix} 8 & 1 \\ 4 & 5 \end{bmatrix},
\quad
\lambda = 9,
\quad
\bm{x} = \begin{bmatrix} 1 \\ 1 \end{bmatrix}
$$

のとき

$$
\bm{A}\bm{x}
=
\begin{bmatrix} 8 & 1 \\ 4 & 5 \end{bmatrix}
\begin{bmatrix} 1 \\ 1 \end{bmatrix}
=
\begin{bmatrix} 9 \\ 9 \end{bmatrix}
= \lambda \bm{x}
$$

が成り立つ。
よって $\lambda=9$ は $\bm{A}$ の固有値である。
固有値 $9$ に対応する固有ベクトル（のひとつ）は $(1,1)^\top$ である。
同時に，$4$ も $A$ の固有値であり，対応する固有ベクトル（のひとつ）は $(1,-4)^\top$ である。
$A$ の固有値は上に挙げた $4,9$ の 2 個だけである。

一般に，正整数 $n$ に対して，$n$ 次正方行列の固有値は高々 $n$ 個存在する。
このことは，固有値が固有方程式 $|A-\lambda I|=0$ の解であることから示される。

```{note}
$n$ 次正方行列 $\bm{A}=(a_{ij})$ に対して，

$$
\sum_{\sigma\in S_n}
  (\sgn \sigma)
  \cdot a_{1,\sigma(1)} a_{2,\sigma(2)} \cdots a_{n,\sigma(n)}
$$

を $\bm{A}$ の**行列式**（determinant）といい，$\det(\bm{A})$ や $|\bm{A}|$ で表す。
ただし，$S_n$ は対称群（$n$ 次の置換全体の集合）であり，置換 $\sigma \in S_n$ に対して $\sgn \sigma$ は置換 $\sigma$ の符号を表す。

行列式の定義より，$|\bm{A}-\lambda I|$ は $\lambda$ の $n$ 次多項式になる。

たとえば，$n=2$ のとき

$$
|\bm{A} - \lambda I|
= \begin{vmatrix} a_{11} - \lambda & a_{12} \\ a_{21} & a_{22} - \lambda \end{vmatrix}
= (a_{11} - \lambda) (a_{22} - \lambda) - a_{12} a_{21}
$$

であり，これは $\lambda$ の 2 次多項式である。
```

ある固有値 $\lambda$ に対応する 固有ベクトルは一つとは限らない。
$\bm{x}$ が $\lambda$ に対する $A$ の固有ベクトルであるとき，
任意のスカラー $c \ne 0$ に対して
$c\bm{x}$ も同様に $\lambda$ に対する $A$ の固有ベクトルである。

たとえば，$\bm{A} = \begin{bmatrix} 8 & 1 \\ 4 & 5 \end{bmatrix}$ の固有値 $9$ に対する固有ベクトルは，

$$
\begin{bmatrix} 1 \\ 1 \end{bmatrix}, \quad
\begin{bmatrix} 2 \\ 2 \end{bmatrix}, \quad
\begin{bmatrix} -2 \\ -2 \end{bmatrix}, \quad
\begin{bmatrix} \sqrt{2} \\ \sqrt{2} \end{bmatrix}, \quad
\cdots
$$

などがある。
これに零ベクトルを加えた集合

$$
\left\{ \bm{x} \in \setC^2 \mid A \bm{x} = 9 \bm{x} \right\}
=
\biggl\{ c \begin{bmatrix} 1 \\ 1 \end{bmatrix} \biggm\vert c \in \setC \biggr\}
$$

を，行列 $A$ の固有値 $9$ に対する **固有空間**（eigenspace）とよぶ。

```{note}
言い換えると，行列 $A$ の $\lambda$ に対する固有空間は，行列 $A - \lambda I$ の核空間（kernel） $\ker(A-\lambda I)$ である。
```

一般に，$n$ 次正方行列 $A$ の固有値 $\lambda$ に対する固有空間

$$
\left\{ \bm{x} \in \setC^n \mid A \bm{x} = \lambda \bm{x} \right\}
$$

は $\setC^n$ の部分ベクトル空間（vector subspace）である。

## おさらい: 対角化

$n$ 次正方行列 $A$ の $n$ 個の固有ベクトルを $\bm{v}_1,\bm{v}_2,\cdots,\bm{v}_n$ とし，すべての固有ベクトルを並べてできる行列を

$$
P=\begin{bmatrix} \bm{v}_1 & \bm{v}_2 & \cdots & \bm{v}_n \end{bmatrix}
$$

とすると，

$$
AP
&= \begin{bmatrix} A \bm{v}_1 & A \bm{v}_2 & \cdots & A \bm{v}_n \end{bmatrix}
\\
&= \begin{bmatrix} \lambda_1 \bm{v}_1 & \lambda_2 \bm{v}_2 & \cdots & \lambda_n \bm{v}_n \end{bmatrix}
\\
&= P D
$$

が成り立つ。ただし，行列 $D$ は $A$ の固有値を対角成分とする対角行列，すなわち，

$$
D=
\begin{bmatrix}
  \lambda_1 & 0 & \cdots & 0 \\
  0 & \lambda_2 & \cdots & 0 \\
  \vdots & \vdots & \ddots & \vdots \\
  0 & 0 & \cdots & \lambda_n
\end{bmatrix}
$$

とする。
もし $P$ が正則ならば $AP=PD$ の両辺に右から逆行列 $P^{-1}$ をかけることによって

$$
A = P D P^{-1}
$$

を得る。
これは，行列 $A$ の分解の一種である。
一般に，$n$ 次正方行列 $A$ が **対角化可能である**（diagonalizable）とは，正則行列 $P$ と対角行列 $D$ を用いて

$$
A = P D P^{-1}
$$

と表されることをいい，このような分解を，行列 $A$ の **固有値分解**（eigendecomposition）とよぶ（この条件を満たす行列 $P$ を対角化行列とよぶ）。

$A$ が対角化可能であれば，$A$ のべき乗は

$$
A^k
&= (P D P^{-1})^k
\\
&= (P D P^{-1}) (P D P^{-1}) \cdots (P D P^{-1})
\\
&= P D (P^{-1} P) D (P^{-1} P) D \cdots (P^{-1} P) D P^{-1}
\\
&= P D^k P^{-1}
$$

のように表すことができる。
対角行列 $D$ の対角成分を $\{d_i\}_{i=1}^{n}$ とすると，$D^k$ は $\{d_i^k\}_{i=1}^{n}$ を対角成分とする対角行列であり，容易に計算できる。

任意の $n$ 次正方行列 $A$ について，以下はすべて同値である。

1. $A$ は対角化可能である。
1. $A$ は $n$ 個の相異なる固有値を持つ。
1. $A$ の固有ベクトルから成る $\setC^n$ の基底が存在する。

$A$ が $n$ 個の相異なる固有値を持つとき，対角化行列 $P$ は $A$ の固有ベクトルを順に並べたものである。

```{prf:example}
:label: example-diag-1

2 次正方行列

$$
A = \begin{bmatrix} 8 & 1 \\ 4 & 5 \end{bmatrix}
$$

とのとき，$A$ の固有値 $4,9$ に対応する固有ベクトルはそれぞれ

$$
\begin{bmatrix} 1 \\ -4 \end{bmatrix}, \quad
\begin{bmatrix} 1 \\ 1 \end{bmatrix}
$$

である。
この固有ベクトルを並べて得られる行列

$$
P = \begin{bmatrix} 1 & 1 \\ -4 & 1 \end{bmatrix}
$$

は

$A$ の対角化行列である。

実際，

$$
  & P^{-1} A P
    = \begin{bmatrix} 1/5 & -1/5 \\ 4/5 & 1/5 \end{bmatrix}
      \begin{bmatrix} 8 & 1 \\ 4 & 5 \end{bmatrix}
      \begin{bmatrix} 1 & 1 \\ -4 & 1 \end{bmatrix}
    = \frac{1}{5}
      \begin{bmatrix} 4 & -4 \\ 36 & 9 \end{bmatrix}
      \begin{bmatrix} 1 & 1 \\ -4 & 1 \end{bmatrix}
    = \begin{bmatrix} 4 & 0 \\ 0 & 9 \end{bmatrix}
$$

が成り立つ。
```

## 例1（固有値）

固有値を計算機プログラムを用いて求める場合，たとえば，以下のように NumPy の [`numpy.linalg.eig()` 関数](https://numpy.org/doc/1.26/reference/generated/numpy.linalg.eig.html) を利用することができる。

```{code-cell}
import numpy as np
from numpy import linalg as LA

a = np.array([
    [8., 1.],
    [4., 5.],
])

w, v = LA.eig(a)

print('eigenvalues:', *w)
print('eigenvectors:')
print(v)
```

この関数の実装には，LAPACK の `_geev` が用いられている。

以下では，固有値を数値計算によって求めるために用いられるいくつかのアルゴリズムや，その背景となる数学的性質を紹介する。

## べき乗法

```{prf:theorem}
:label: thm-dominant-eigenvector

$A$ は対角化可能な行列で，$A$ の固有値 $\lambda_1,\lambda_2,\ldots,\lambda_n$ は

$$
|\lambda_1| > |\lambda_2| \ge |\lambda_3| \ge \cdots \ge |\lambda_n|
$$

を満たすとする。初期ベクトル $\bm{x}^{(0)}$ を適当に選ぶと，漸化式

$$
\bm{x}^{(k)} = A \bm{x}^{(k-1)}, \quad k=1,2,3,\ldots
$$

によって定まる列 $\bm{x}^{(1)},\bm{x}^{(2)},\ldots$ は行列 $A$ の絶対値最大の固有値 $\lambda_1$ に対応する固有ベクトルに収束する。
```

```{prf:proof}
$A$ が対角化可能であるから，固有ベクトル $\bm{v}_1,\bm{v}_2,\ldots,\bm{v}_n$ は基底を作る。
つまり，任意の初期ベクトル $\bm{x}^{(0)}$ に対して実数 $c_1,c_2,\ldots,c_n$ が存在して

$$
\bm{x}^{(0)} = c_1 \bm{v}_1 + c_2 \bm{v}_2 + \cdots + c_n \bm{v}_n
$$

が成り立つ。このとき

$$
  & \bm{x}^{(1)}
    = A \bm{x}^{(0)}
    = c_1 A \bm{v}_1 + c_2 A \bm{v}_2 + \cdots + c_n A \bm{v}_n
\\&\quad
    = c_1 \lambda_1 \bm{v}_1 + c_2 \lambda_2 \bm{v}_2 + \cdots + c_n \lambda_n \bm{v}_n
$$

同様に，任意の自然数 $k$ に対して

$$
\bm{x}^{(k)}
= c_1 \lambda_1^k \bm{v}_1 + c_2 \lambda_2^k \bm{v}_2 + \cdots + c_n \lambda_n^k \bm{v}_n
$$

が成り立つ。$\lambda_1 \neq 0$ であるから，

$$
\bm{x}^{(k)}
= \lambda_1^k \left\{
    c_1 \bm{v}_1
    + c_2 \biggl( \frac{\lambda_2}{\lambda_1} \biggr)^k \bm{v}_2
    +\cdots
    + c_n \biggl( \frac{\lambda_n}{\lambda_1} \biggr)^k \bm{v}_n
  \right\}
$$

$k=2,3,\cdots$ について $\Bigl|\frac{\lambda_k}{\lambda_1}\Bigr|<1$ に注目すると，$\displaystyle\lim_{k\to\infty} \frac{\bm{x}^{(k)}}{\lambda_1^k} = c_1 \bm{v}_1$ である。
$c_1 \neq 0$ となるように初期ベクトル $\bm{x}^{(0)}$ を選択すると，
$\bm{x}^{(k)}$ は $\bm{v}_1$ のスカラー倍，すなわち固有値 $\lambda_1$ に対応する $A$ の固有ベクトルに収束する。
```

この定理により，$\bm{x},A\bm{x},A^2\bm{x},\cdots$ を次々と求めることによって $A$ の固有ベクトル（のひとつ）が求まることが分かる。

この定理に基づいて固有ベクトルを求める方法を **べき乗法**（power iteration）または累乗法とよぶ。

ただし，$\bm{x}^{(k)}$ があるひとつのベクトルに収束するのではない点に注意する必要がある。なぜなら，固有ベクトルのスカラー倍も固有ベクトルであり，その結果，発散したり零ベクトルに収束したりする可能性が残るからである。

そこでベクトルの正規化（ノルムが $1$ になるようにスカラー倍する）を用いる。
すると，$\displaystyle\frac{\bm{x}^{(k)}}{|\bm{x}^{(k)}|}$ は $\displaystyle\frac{\bm{v}_1}{|\bm{v}_1|}$ に収束する。

$\bm{x}^{k}$ を正規化して更新式を

$$
\bm{x}^{(k)} = \frac{A \bm{x}^{(k-1)}}{|A \bm{x}^{(k-1)}|}
$$

のようにすれば，正規化された固有ベクトル $\displaystyle\frac{\bm{v}_1}{|\bm{v}_1|}$ に収束する。

仮に行列 $A$ の絶対値最大の固有値が複数存在する場合も同様のことが言える（収束先の固有ベクトルは初期値に応じて変わる）。

固有ベクトルが求められれば，対応する固有値は **レイリー商**（Rayleigh quotient）を用いて計算できる。

$$
R(\bm{x})
 = \frac{\naiseki{\bm{x}}{A\bm{x}}}{\naiseki{\bm{x}}{\bm{x}}}
 = \frac{\naiseki{\bm{x}}{\lambda\bm{x}}}{\naiseki{\bm{x}}{\bm{x}}}
 = \frac{\lambda\naiseki{\bm{x}}{\bm{x}}}{\naiseki{\bm{x}}{\bm{x}}}
 = \lambda
$$

ただし，$\naiseki{\cdot}{\cdot}$ はベクトルの内積である。
特に，$|\bm{x}|=1$ と正規化すれば，
$R(\bm{x})=\naiseki{\bm{x}}{A\bm{x}}$ が成立する。

このようにして，絶対値最大の固有値 $\lambda_1$ を求めることができる。

任意の $n$ 次元行列に対して，その行列の固有値を求める問題の時間計算量は $O(n^3)$ であることが知られている。
つまり，固有値を求めるために要する時間は，行列の次元 $n$ の 3 乗に比例して増加する。
一方，べき乗法は更新に $O(n^2)$ を要する。
行列のサイズと収束の速さによって，べき乗法で効率良く固有ベクトルを求めることができる場合がある。

## 例2（べき乗法）

べき乗法によって $A = \begin{bmatrix} 8 & 1 \\ 4 & 5 \end{bmatrix}$ の最大固有値を求める Python プログラムの例を以下に示す。

```{code-cell}
import numpy as np

# 小数点以下の桁数
np.set_printoptions(precision=6)

def power_iteration(A, x0, eps=1e-7, maxrepeat=100, verbose=True):
    # 初期値
    x = x0

    print(0, x)
    for k in range(1, maxrepeat + 1):
        x_old = x
        x = A @ x
        x = x / np.linalg.norm(x)
        diff = np.linalg.norm(x - x_old)
        if verbose:
            print(k, x, diff)
        if diff < eps:
            break

    return x

# 行列
A = np.array([[8., 1.], [4., 5.]])

# 初期ベクトル
x0 = np.array([1., 0.,])

# べき乗法によって固有ベクトルを求める
x = power_iteration(A, x0)

# 絶対値最大固有値を求める
eigenvalue = np.dot(x, A @ x) / np.dot(x, x)

# 結果を表示する
print()
print('eigenvector =', x)
print('eigenvalue  =', eigenvalue)
```

## QR法

べき乗法は，固有値のうち絶対値が最大のものだけを求めるアルゴリズムであった。
すべての固有値を求める数値解法のうち代表的なものが **QR 法**（QR method）である。

QR 法は，行列の QR 分解を利用する方法である。
QR 分解にはいくつかのアルゴリズム（Givens 回転，Householder 変換，Gram-Schmidt 直交化）が知られている。
ここでは，Gram-Schmidt 直交化を

行列 $A,B$ をそれぞれ $n$ 次正方行列とする。
正則行列 $M$ によって $B = M^{-1} A M$ が成り立つとき，$A$ と $B$ とは **相似**（similar）である，という。
また，$M^{-1} A M$ を $A$ の **相似変換**（similarity transformation）という。

```{prf:theorem}
:label: thm-simtrans-eigenvector

行列 $A$ の固有値を $\lambda$，対応する固有ベクトルを $\bm{x}$ とする。
このとき，$A$ の相似変換 $M^{-1} A M$ も固有値 $\lambda$ をもち，固有値 $\lambda$ に対する $M^{-1} A M$ の固有ベクトルは $M^{-1} \bm{x}$ である。

特に，$M$ が直交行列であるとき，固有ベクトルは $M^\top\bm{x}$ である。
```

たとえば，行列 $A$ を $P D P^{-1}$ のように対角化可能であるならば，行列 $A$ と行列 $D$ とは相似である。

```{prf:theorem} グラム・シュミットの直交化
:label: thm-gram-schmidt

一次独立な $n$ 個のベクトル $\bm{a}_1,\bm{a}_2,\ldots,\bm{a}_n$ を考える。
このとき，$k=1,2,\ldots$ について漸化式

$$
& \bm{u}_k
  = \bm{a}_k
  - \sum_{j=1}^{k-1} \naiseki{\bm{a}_k}{\bm{q}_j} \cdot \bm{q}_j \\
& \bm{q}_k = \frac{\bm{u}_k}{\|\bm{u}_k\|_2}
$$

によって定まるベクトルの列 $\{\bm{q}_k\}$ は正規直交基底をなす。
すなわち，以下が成り立つ。

$$
\naiseki{\bm{q}_i}{\bm{q}_j} = \delta_{ij}
$$

ただし $\delta_{ij}$ はクロネッカーのデルタである。
```

```{prf:example}
:label: example-gram-schmidt-1

3 つのベクトル $\bm{a}_1$，$\bm{a}_2$，$\bm{a}_3$ を

$$
\bm{a}_1 = (1,2,3)^\top,
\quad
\bm{a}_2 = (1,1,2)^\top,
\quad
\bm{a}_3 = (1,1,1)^\top
$$

と定義する。
このとき，

$$
\bm{u}_1
  &= \bm{a}_1
   = (1,2,3)^\top,
\\
\bm{q}_1
  &= \frac{\bm{u}_1}{\|\bm{u}_1\|_2}
   = \frac{1}{\sqrt{14}} (1,2,3)^\top,
\\
\bm{u}_2
  &= \bm{a}_2 - \naiseki{\bm{a}_2}{\bm{q}_1} \cdot \bm{q}_1
\\
  &= (1,1,2)^\top - \frac{1+2+6}{\sqrt{14}} \cdot \frac{1}{\sqrt{14}} (1,2,3)^\top
   = \frac{1}{14} ( 5, -4, 1 )^\top,
\\
\bm{q}_2
  &= \frac{1}{\sqrt{42}} ( 5, -4, 1 )^\top,
\\
\bm{u}_3
  &= \bm{a}_3
     - \naiseki{\bm{a}_3}{\bm{q}_1} \cdot \bm{q}_1
     - \naiseki{\bm{a}_3}{\bm{q}_2} \cdot \bm{q}_2
\\
  &= (1,1,1)^\top
     - \frac{1+2+3}{\sqrt{14}} \cdot \frac{1}{\sqrt{14}} (1,2,3)^\top
     - \frac{5-4+1}{\sqrt{42}} \cdot \frac{1}{\sqrt{42}} ( 5, -4, 1 )^\top
\\
  &= (1,1,1)^\top
     - \frac{3}{7} \cdot (1,2,3)^\top
     - \frac{1}{21} \cdot (5,-4,1)^\top
\\
  &= \frac{1}{21} ( 21 - 9 - 5, 21 - 18 + 4, 21 - 27 - 1 )^\top
   = \frac{1}{3} (1,1,-1)^\top,
\\
\bm{q}_3
  &= \frac{1}{\sqrt{3}} (1,1,-1)^\top
$$

以上により，正規直交基底を得る。

$$
\bm{q}_1
  &= \frac{1}{\sqrt{14}} (1,2,3)^\top,
\\
\bm{q}_2
  &= \frac{1}{\sqrt{42}} ( 5, -4, 1 )^\top,
\\
\bm{q}_3
  &= \frac{1}{\sqrt{3}} (1,1,-1)^\top
$$

これが正規直交基底であることは，$|\bm{q}_k|=1$，$k=1,2,3$ および

$$
\naiseki{\bm{q}_1}{\bm{q}_2}
  &= \frac{1}{\sqrt{14 \cdot 42}} (1 \cdot 5 + 2 \cdot (-4) + 3 \cdot 1)
   = \frac{1}{14 \sqrt{3}} (5-8+3)
   = 0,
\\
\naiseki{\bm{q}_1}{\bm{q}_3}
  &= \frac{1}{\sqrt{14 \cdot 3}} (1 \cdot 1 + 2 \cdot 1 + 3 \cdot (-1))
   = \frac{1}{\sqrt{42}} (1+2-3)
   = 0,
\\
\naiseki{\bm{q}_2}{\bm{q}_3}
  &= \frac{1}{\sqrt{42 \cdot 3}} (5 \cdot 1 + (-4) \cdot 1 + 1 \cdot (-1))
   = \frac{1}{3 \sqrt{14}} (5-4-1)
   = 0
$$

から確かめられる。
```

正規直交基底 $\{\bm{q}_1,\bm{q}_2,\ldots,\bm{q}_n\}$ に対して $Q=[\bm{q}_1 \bm{q}_2 \cdots \bm{q}_n]$ とおくと $QQ^\top=Q^\top Q=I_n$ であるから $Q$ は直交行列である。

```{prf:theorem} QR分解
:label: thm-qr-decomposition

$A$ を正則な $n$ 次実行列とする。
このとき，次のような分解が一意に可能である。

$$
A = QR
$$

ただし，$Q$ は直交行列であり，$R$ は対角成分が正であるような上三角行列である。
```

```{prf:proof}

まず，分解の存在を示す。

$A=[\bm{a}_1, \bm{a}_2, \cdots, \bm{a}_n]$ とおく。
$A$ は正則行列であるから，$\{\bm{a}_1, \bm{a}_2, \ldots, \bm{a}_n\}$ は基底をなす。
グラム・シュミットの直交化によって，正規直交基底 $\{\bm{q}_1,\bm{q}_2,\ldots,\bm{q}_n\}$
を得る。このとき，上三角行列 $R=[r_{jk}]$ を以下のように定義する。

$$
  r_{jk} =
  \begin{cases}
    \naiseki{\bm{a}_k}{\bm{q}_j} & j < k
    \\
    \| \bm{u}_k \|_2 & j = k
    \\
    0 & j > k
  \end{cases}
$$

ただし，

$$
  \bm{u}_k
  = \bm{a}_k - \sum_{j=1}^{k-1} \naiseki{\bm{a}_k}{\bm{q}_j} \bm{q}_j
  = \bm{a}_k - \sum_{j=1}^{k-1} r_{jk} \bm{q}_j
$$

このように定義した $A,Q,R$ について，$A=QR$ が成り立つ。

次に，分解の一意性を示す。

$A$ が以下のように 2 通りに QR 分解できるとしたうえで，両者が一致することを示せばよい。

$$
A = Q_1 R_1 = Q_2 R_2
$$

$A$ が正則であることから，$Q_1$，$Q_2$，$R_1$，$R_2$ も正則である。

$$
& Q_2 = Q_1 R_1 R_2^{-1}
\\
& \therefore Q_2^{-1} = R_2 R_1^{-1} Q_1^{-1}
\\
& \therefore Q_2^{-1} Q_1 = R_2 R_1^{-1}
$$

行列 $R_1$，$R_2$ は，いずれも対角成分が正であるような上三角行列である。
したがって，逆行列 $R_1^{-1}$ も対角成分が正であるような上三角行列である。2 つの対角成分が正であるような上三角行列の積は，やはり，対角成分が正であるような上三角行列であるから，$R_2 R_1^{-1}$ は対角成分が正であるような上三角行列である。

一方，$R_2 R_1^{-1}$ は直交行列でもある。
なぜならば，$R_2 R_1^{-1} = Q_2^{-1} Q_1$ であり，直交行列の積はやはり直交行列であるからである。

以上により，$R_2 R_1^{-1}$ は対角成分が正であるような上三角行列であり，かつ，直交行列である。
そのような行列は単位行列しか存在しない。
つまり，以下が成り立つ。

$$
& Q_2^{-1} Q_1 = R_2 R_1^{-1} = I
\\
& \therefore Q_1 = Q_2, \quad R_1 = R_2
$$

任意の 2 つの QR 分解が一致することから，QR 分解は一意である。
```

ここで，

$$
  A_{1} = A, \quad
  A_{k} \to Q_{k} R_{k}, \quad
  A_{k+1} = R_{k} Q_{k} \to Q_{k+1} R_{k+1} \quad (k=1,2,\ldots)
$$

のように行列の列 $A_1,A_2,\ldots$ をそれぞれ QR 分解する操作を繰り返すことを考える。ただし $\to$ は QR 分解を意味する。

すなわち以下が成り立つ。

$$
  A_{1} &= A = Q_{1} R_{1}, \\
  A_{2} &= R_{1} Q_{1} = (Q_{1}^{-1} A) Q_{1} = Q_1^t A_1 Q_1 \\
        &= Q_{2} R_{2}, \\
  A_{3} &= R_{2} Q_{2} = Q_2^t A_2 Q_2 = Q_2^t (Q_1^t A_1 Q_1) Q_2
         = (Q_1 Q_2)^t A_1 Q_1 Q_2 \\
        &= Q_{3} R_{3}, \\
$$

同様に繰り返すと以下の関係が得られる。

$$
  A_{k} = Q_{k} R_{k}, \\
  A_{k+1} = (Q_1 Q_2 \cdots Q_k)^t A_1 (Q_1 Q_2 \cdots Q_k)
$$

すなわち $A_{k+1}$ は $A$ の相似変換となっている。
したがって，$A_{k+1}$ の固有値は $A$ の固有値と一致する。

1. $A_1 \leftarrow A$ とする。
1. $k=1,2,\ldots$ について以下を繰り返す。
   1. QR 分解 $A_k = Q_k R_k$ を求める。
   1. $A_{k+1} \leftarrow R_k Q_k$ を求める。

ある条件のもとで，$A_{k}$ が上三角行列に収束し，対角成分に固有値が絶対値が大きい順に並ぶことが知られている。

このようにしてすべての固有値を求める方法が QR 法である。

## 例3（QR 分解）

たとえば $A=\begin{bmatrix} 8 & 1 \\ 4 & 5 \end{bmatrix}$ のとき，
グラム・シュミットの直交化により，

$$
A_1
  &=
   \begin{bmatrix}
     8 & 1 \\ 4 & 5
   \end{bmatrix}
   =
   \begin{bmatrix}
     \frac{2}{\sqrt{5}} & \frac{-1}{\sqrt{5}} \\
     \frac{1}{\sqrt{5}} & \frac{2}{\sqrt{5}}
   \end{bmatrix}
   \begin{bmatrix}
     4\sqrt{5} & \frac{7}{\sqrt{5}} \\
     0 & \frac{9}{\sqrt{5}}
   \end{bmatrix}
\\
A_2
  &=
   \begin{bmatrix}
     4\sqrt{5} & \frac{7}{\sqrt{5}} \\
     0 & \frac{9}{\sqrt{5}}
   \end{bmatrix}
   \begin{bmatrix}
     \frac{2}{\sqrt{5}} & \frac{-1}{\sqrt{5}} \\
     \frac{1}{\sqrt{5}} & \frac{2}{\sqrt{5}}
   \end{bmatrix}
  =
  \begin{bmatrix}
    47/5 & -6/5 \\
    9/5 & 18/5
  \end{bmatrix}
  =
  \begin{bmatrix}
    9.4 & -1.2 \\
    1.8 & 3.6
  \end{bmatrix} \\
  &=
  \begin{bmatrix}
    0.9822 & -0.1881 \\
    0.1881 &  0.9822
  \end{bmatrix}
  \begin{bmatrix}
    9.5708 & -0.5015 \\
    0.0000 &  3.7614
  \end{bmatrix}
\\
A_3
  &=
  \begin{bmatrix}
    9.5708 & -0.5015 \\
    0.0000 &  3.7614
  \end{bmatrix}
  \begin{bmatrix}
    0.9822 & -0.1881 \\
    0.1881 &  0.9822
  \end{bmatrix}
  =
  \begin{bmatrix}
    9.3061 & -2.2929 \\
    0.7075 &  3.6944
  \end{bmatrix} \\
  &=
  \begin{bmatrix}
    0.9971 & -0.0758 \\
    0.0758 &  0.9971
  \end{bmatrix}
  \begin{bmatrix}
    9.3330 & -2.0062 \\
    0.0000 &  3.8576
  \end{bmatrix}
$$

これを繰り返すと，$A_k$ は上三角行列

$$
  \begin{bmatrix}
    9 & -3 \\
    0 &  4
  \end{bmatrix}
$$

に収束する。
このときの対角成分 $9,4$ が $A$ の固有値である。

Python プログラムの例は，以下のようになる。

```{code-cell}
import sympy as sym
from sympy import Matrix
from IPython.display import display

# SymPy を用いる記号計算バージョン

def round(a, digits):
    m, n = a.shape
    for i in range(m):
        for j in range(n):
            a[i, j] = a[i, j].round(digits)

def qr_decompose_sym(a):
    n = a.shape[0]
    assert a.shape == (n, n)
    u = sym.zeros(n, n)
    q = sym.zeros(n, n)
    r = sym.zeros(n, n)
    for k in range(n):
        u[:, k] = a[:, k]
        for j in range(k):
            r[j, k] = a[:, k].dot(q[:, j])
            u[:, k] -= r[j, k] * q[:, j]
        r[k, k] = sym.sqrt(u[:, k].dot(u[:, k]))
        q[:, k] = u[:, k] / r[k, k]
    return q, r

def qr_method_sym(a, repeat=20):
    n = a.shape[0]
    assert a.shape == (n, n)

    for k in range(repeat):
        q, r = qr_decompose_sym(a)
        round(q, 4)
        round(r, 4)
        a = sym.Mul(r,q)

    return a

a = Matrix([
    [8, 1],
    [4, 5],
])
display(a)
a = qr_method_sym(a)
display(a)
```

大規模行列については，計算量を減らすために $A$ をヘッセンベルグ行列に相似変換してから QR 法を適用するとよい。詳しくは参考書（[幸谷 2021]など）を参照。
