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

# 関数近似と補間 (2) 関数近似

前回扱った関数補間は，与えられた点を通るような関数を求める問題であった。
つまり，いくつかの点 $(x_0,y_0),(x_1,y_1),\cdots,(x_n,y_n)$ が与えられたとき，各 $k$ について $f(x_k)=y_k$ が成り立つような関数 $f$ を求める問題であった。

不明な関数を決定するという目的において，点列 $\{(x_k,y_k)\}$ はいわば「正解」であり， 関数 $f(x_k)$ の値は厳密に $y_k$ に等しい必要がある。

一方，ある現象を観測し得られた結果に誤差が含まれる場合のように，必ずしも $f(x_k)=y_k$ が厳密に成り立つ必要のない状況も考えられる。
この場合も，できるだけ $f(x_k)=y_k$ に近づけることが望ましいと考えられる。
つまり，与えられたデータにできるだけ「近い」関数を求めるという問題を考えることができる。
このときに関数の「近さ」を何らかの形で定義する必要がある。
最も基本的なものは **最小二乗関数近似**（least squares function approximation）である。

## 連続関数の最小二乗近似

関数 $f$ を関数 $g$ で近似するとして，その近さを測るには，関数の差のノルム $\|f-g\|$ を用いる。

まず，関数のノルム[^norm]を定義する。

[^norm]: 関数のノルムとは，ベクトルの長さのようなものである。[Online Etymology Dictionary](https://www.etymonline.com/search?q=norm) によれば，ノルム（norm）の語源はフランス語 norme（標準，規準）やラテン語 norma（大工の直角定規）とされる。

```{prf:definition} 関数のノルム
:label: thm-fun-norm

閉区間 $[a,b]$ 上で定義された関数 $f$ および自然数 $p$ について，

$$
  \left( \int_a^b |f(x)|^p dx \right)^{1/p}
$$

を関数 $f$ の **$p$-ノルム**（$p$-norm）とよび，$\|f\|_{p,[a,b]}$ と書く。
定義域 $[a,b]$ が文脈から明らかであれば，単に $\|f\|_p$ と書くこともある。
```

閉区間 $[a,b]$ 上で定義された任意の関数 $f$，$g$ について，関数のノルム $\|\cdot\|$ は（ベクトルの長さと同様に）以下の性質を持つ。

* 独立性: $\|f\| = 0$ ならば $(\forall x \in [a,b])[f(x)=0]$

* 斉次性: 任意の実数 $\alpha$ について，$\|\alpha f\| = |\alpha| \|f\|$

* 劣加法性: $\|f + g\| \leq \|f\| + \|g\|$

```{prf:example}
:label: example-fun-norm-1

閉区間 $[0,2]$ 上で定義された関数 $f(x)=x(x-1)$ の 2-ノルム $\|f\|_2$ は以下のように求めることができる。

$$
  & ( \|f\|_2 )^2
    = \int_0^2 |x(x-1)|^2 dx
    = \int_0^2 x^2 (x-1)^2 dx
    = \int_0^2 (x^4 - 2 x^3 + x^2) dx
\\&\quad
    = \left[ \frac{1}{5} x^5 - \frac{2}{4} x^4 + \frac{1}{3} x^3 \right]_0^2
    = \frac{32}{5} - 8 + \frac{8}{3}
    = \frac{96-120+40}{15}
    = \frac{16}{15}
\\& \therefore \|f\|_2 = \sqrt{\frac{16}{15}} = \frac{4\sqrt{15}}{15}
$$
```

```{prf:example}
:label: example-fun-norm-2

閉区間 $[0,1]$ 上で定義された関数が定数 $f(x)=c$ であるならば，

$$
  \|f\|_{p}
  = \left( \int_0^1 |c|^p dx \right)^{1/p}
  = \left( \left[ |c|^p \cdot x \right]_0^1 \right)^{1/p}
  = \left( |c|^p \right)^{1/p}
  = |c|
$$
```

```{prf:example}
:label: example-fun-norm-3

閉区間 $[0,1]$ 上で定義された関数が $f(x)=x^n$ （$n$ は正の整数）であるならば，

$$
  \|f\|_2
  = \left( \int_0^1 (x^n)^2 dx \right)^{1/2}
  = \sqrt{ \left[ \frac{1}{2n+1} x^{2n+1} \right]_0^1 }
  = \sqrt{ \frac{1}{2n+1} }
$$
```

$\newcommand\naiseki[2]{\left\langle #1 , #2 \right\rangle}$
$\newcommand\bm[1]{{\boldsymbol #1}}$
$\newcommand\setR{\mathbb{R}}$

```{prf:definition} 関数の内積
:label: thm-fun-iprod

閉区間 $[a,b]$ 上で定義された 2 つの関数 $f,g$ について，

$$
  \int_a^b f(x) g(x) dx
$$

を関数 $f$ と $g$ の **内積**（inner product）とよび，$\naiseki{f}{g}$ と書く。
```

定義により，演算 $\naiseki{\cdot}{\cdot}$ は内積の基本的性質を満たす。
すなわち，閉区間 $[0,1]$ 上で定義された任意の実関数 $f$，$g$，$h$ について，以下が成り立つ。

* 線形性: 任意の実数 $p$，$q$ に対して $\naiseki{pf+qg}{h} = p \naiseki{f}{h} + q \naiseki{g}{h}$

* 対称性: $\naiseki{f}{g} = \naiseki{g}{f}$

* 正定値性: $\naiseki{f}{f} \geq 0$，等号成立は $(\forall x) [f(x)=0]$ のときに限る

つまり，関数の内積は，一般のベクトル空間における内積と同様の性質を持つ。

```{prf:example}
:label: example-fun-iprod-1

たとえば，非負整数 $m$，$n$ について，閉区間 $[0,1]$ 上で定義された関数 $f(x)=x^m$ と $g(x)=x^n$ の内積 $\left\langle f, g \right\rangle$ は以下のように求めることができる。

$$
  \naiseki{f}{g}
  = \int_0^1 x^{m+n} dx
  = \left[ \frac{1}{m+n+1} x^{m+n+1} \right]_0^1
  = \frac{1}{m+n+1}
$$
```

定義から分かる通り，$\naiseki{f}{f}$ は $\|f\|_2^2=(\|f\|_2)^2$ と等しい。

```{prf:definition} 1次従属，1次独立
:label: thm-fun-linearly-independent

閉区間 $[a,b]$ 上に関数系 $\{\phi_i(x)\}_{i=1}^{n}$ が定義されている。

このとき，係数 $c_i,\ (i=1,2,\ldots,n)$ について

$$
  & \left( \forall x \in [a,b] \right)
    \bigl[ c_1 \phi_1(x) + c_2 \phi_2(x) + \cdots + c_n \phi_n(x) = 0 \bigr]
\\&
    \quad\Longrightarrow\quad
    c_1 = c_2 = \cdots = c_n = 0
$$

が成立するならば，関数系 $\{\phi_i(x)\}_{i=1}^{n}$ は **1次独立である**（linearly independent）という。
1次独立でない関数系は **1次従属である**（linearly dependent）。
```

関数系 $\{\phi_i(x)\}_{i=1}^{n}$ が1次従属であるならば，少なくとも 1 つ以上が非零であるような係数 $c_i,\ (i=1,2,\ldots,n)$ について

$$
\left( \forall x \in [a,b] \right)
\bigl[ c_1 \phi_1(x) + c_2 \phi_2(x) + \cdots + c_n \phi_n(x) = 0 \bigr]
$$

が成り立つ。
このとき，ある $k$ について係数 $c_k$ は非零であるから，$\phi_k(x)$ を他の関数 $\{\phi_i(x)\}_{i \neq k}$ の1次結合で表せる。

$$
  \phi_k(x)
  = c'_1 \phi_1(x)
    + c'_2 \phi_2(x) + \cdots
    + c'_{k-1} \phi_{k-1}(x)
    + c'_{k+1} \phi_{k+1}(x) + \cdots
    + c'_n \phi_n(x)
$$

ただし，$c'_i = -\frac{c_i}{c_k}$ （$i=1,2,\cdots,i-1,i+1,\cdots,n$）とする。

````{prf:example}
:label: example-fun-lindep-1

たとえば，任意の閉区間 $[a,b]\subseteq\mathbb{R}$ および正の整数 $n$ に対して $1,x,x^2,x^3,\cdots,x^n$ は 1 次独立である。

このことを証明してみよう。

実数 $c_0,c_1,\cdots,c_n$ が任意の $x \in [a,b]$ に対して以下を満たすと仮定する。

```{math}
:label: eq-fun-lindep-1
\sum_{k=0}^{n} c_k x^k = c_0 x^0 + c_1 x^1 + c_2 x^2 + \cdots + c_n x^n = 0
```

ここで，$c_0=c_1=\cdots=c_n=0$ でないと仮定する。
このとき，上式 {eq}`eq-fun-lindep-1` は高々 $n$ 次の代数方程式であるから，その解は高々 $n$ 個である。
つまり，閉区間 $[a,b]$ に相異なる $n+1$ 個の点 $x_1,x_2,\cdots,x_{n+1}$ をとると，少なくとも 1 個の点 $x=x_k$ において等式 {eq}`eq-fun-lindep-1` は成り立たない。
これは仮定に反する。

したがって，背理法により，$c_0=c_1=\cdots=c_n=0$ である。
以上より，$1,x,x^2,x^3,\cdots,x^n$ は 1 次独立である。
````

以上の定義を用いて以下の定理を導入する。

```{prf:theorem}
:label: thm-minsq-con

閉区間 $[a,b]$ 上に連続関数 $f(x)$ と1次独立な関数系 $\{\phi_k(x)\}_{k=1}^{n}$ を定義する。
実数 $c_k\ (k=1,2,\ldots,n)$ に対して，関数 $g(x)$ を $g(x)=\sum_{k=1}^{n} c_k \phi_k(x)$ と定義する。
このとき，ノルム $\|f(x)-g(x)\|_2$ を最小にする係数 $\{c_k\}$ は以下の連立1次方程式の解である。

$$
\begin{bmatrix}
  \naiseki{\phi_1}{\phi_1} & \naiseki{\phi_1}{\phi_2} & \cdots & \naiseki{\phi_1}{\phi_n} \\
  \naiseki{\phi_2}{\phi_1} & \naiseki{\phi_2}{\phi_2} & \cdots & \naiseki{\phi_2}{\phi_n} \\
  \vdots & \vdots & \ddots & \vdots \\
  \naiseki{\phi_n}{\phi_1} & \naiseki{\phi_n}{\phi_2} & \cdots & \naiseki{\phi_n}{\phi_n}
\end{bmatrix}
\begin{bmatrix}
  c_1 \\ c_2 \\ \vdots \\ c_n
\end{bmatrix}
=
\begin{bmatrix}
  \naiseki{f}{\phi_1} \\ \naiseki{f}{\phi_2} \\ \vdots \\ \naiseki{f}{\phi_n}
\end{bmatrix}
$$
```

```{prf:proof}
関数 $f$ と関数 $g$ の差の 2-ノルムの 2 乗は以下のようになる。

$$
(\|f-g\|_2)^2
& = \naiseki{f-g}{f-g}
  = \naiseki{f}{f} - 2 \naiseki{f}{g} + \naiseki{g}{g}
\\
& = \naiseki{f}{f} - 2 \naiseki{f}{\sum_{k=1}^{n} c_k \phi_k} + \naiseki{\sum_{k=1}^{n} c_k \phi_k}{\sum_{k=1}^{n} c_k \phi_k}
\\
& = \naiseki{f}{f} - 2 \sum_{k=1}^{n} c_k \naiseki{f}{\phi_k} + \sum_{i=1}^{n} \sum_{j=1}^{n} c_i c_j \naiseki{\phi_i}{\phi_j}
$$

これを $F(c_1,c_2,\cdots,c_n)$ とおく。
関数 $F$ が最小となるような $c_1,c_2,\cdots,c_n$ を求めればよい。
$F$ が最小値をとるとき各 $k$ について偏微分 $\frac{\partial F}{\partial c_k}$ が $0$ になる。
$F$ を $c_k$ で偏微分すると

$$
\frac{\partial F}{\partial c_k}
& = - \frac{\partial}{\partial c_k} \left( 2 \sum_{k=1}^{n} c_k \naiseki{f}{\phi_k} \right)
    + \frac{\partial}{\partial c_k} \left( \sum_{i=1}^{n} \sum_{j=1}^{n} c_i c_j \naiseki{\phi_i}{\phi_j} \right)
\\
& = - 2 \naiseki{f}{\phi_k}
    + 2 \sum_{i=1}^{n} c_i \naiseki{\phi_i}{\phi_k}
$$

であるから，

$$
\frac{\partial F}{\partial c_k} = 0
\iff
\sum_{i=1}^{n} c_i \naiseki{\phi_i}{\phi_k} = \naiseki{f}{\phi_k}
$$

これは，この定理に示す連立1次方程式そのものである。
```

## 例（連続関数の最小二乗近似）

たとえば，閉区間 $[0,1]\subseteq\mathbb{R}$ において，関数 $f(x)=x^2$ を関数系 $1,x$ の 1 次結合 $g(x) = c_1 + c_2 x$ で近似することを考える。

ノルム $\|f(x)-g(x)\|_2$ を最小にするには，{prf:ref}`thm-minsq-con` より，以下の方程式を $c_1,c_2$ について解けばよい。

$$
\begin{align}
  \begin{bmatrix}
    \naiseki{1}{1} & \naiseki{1}{x} \\
    \naiseki{x}{1} & \naiseki{x}{x}
  \end{bmatrix}
  \begin{bmatrix}
    c_1 \\ c_2
  \end{bmatrix}
  =
  \begin{bmatrix}
    \naiseki{x^2}{1} \\ \naiseki{x^2}{x}
  \end{bmatrix}
\end{align}
$$

ここで，$m+n \geq 0$ のとき

$$
\begin{align}
  \naiseki{x^m}{x^n}
  = \frac{1}{m+n+1}
\end{align}
$$

であるから，上の方程式は以下と等価である。

$$
\begin{align}
  \begin{bmatrix}
    1 & 1/2 \\
    1/2 & 1/3
  \end{bmatrix}
  \begin{bmatrix}
    c_1 \\ c_2
  \end{bmatrix}
  =
  \begin{bmatrix}
    1/3 \\ 1/4
  \end{bmatrix}
\end{align}
$$

これを解くと $(c_1,c_2)=(-1/6,1)$ を得る。
したがって，$f(x) = x^2$ の最小2乗近似1次式は以下のようになる。

$$
\begin{align}
  g(x) = - \frac{1}{6} + x
\end{align}
$$

$y=x^2$および最小2乗近似1次関数は下図のようになる。

```{code-cell}
:tags: [remove-input, remove-stderr]

from sympy import pi, plot, sin, sympify
from sympy.abc import x

plot(x ** 2, sympify('-1/6+x'), (x, -0.1, 1.1),
    legend=True,
);
```

```{code-cell}
:tags: [remove-input, remove-stderr]

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange

xmin, xmax = -0.1, 1.1

# y = x^2
X = np.linspace(xmin, xmax, 101)
Y = X ** 2
plt.plot(X, Y, color='blue', label='$y=x^2$')

# y = -1/6 + x
X = np.linspace(0, 1, 101)
Y = -1/6 + X
plt.plot(X, Y, color='red', label=r'$y=-\frac{1}{6}+x$')

plt.scatter([0,1], [0,1], color='red')
plt.tick_params(labelsize=18)
plt.gca().set_xlim(xmin, xmax)
plt.gca().set_ylim(-.2, 1.1)
plt.grid(True)
plt.legend(fontsize=18)
plt.show()
```

## 離散データの最小二乗近似

次に，$m$ 組の与えられたデータ
$(x_1,y_1),(x_2,y_2),\ldots,(x_m,y_m)$
を通る未知の関数 $f(x)$ を
関数 $g(x) = \sum_{k=1}^{n} c_k \phi_k(x)$
で近似することを考える。

このとき，2 つのベクトル

$$
  & \bm{f} = (f(x_1),f(x_2),\ldots,f(x_m)) = (y_1,y_2,\ldots,y_m)
\\
  & \bm{g} = (g(x_1),g(x_2),\ldots,g(x_m))
    = \left( \sum_k c_k \phi_k(x_1),\sum_{k} c_k \phi_k(x_2),\ldots,\sum_k c_k \phi_k(x_m) \right)
$$

に対して，$F(c_1,c_2,\ldots,c_n)$ を以下のように定義する。
ただし，ベクトル $\bm{v}$ に対して $\|\bm{v}\|$ はベクトル $\bm{v}$ のノルムである。

$$
F(c_1,c_2,\ldots,c_n)
= \| \bm{f} - \bm{g} \|^2
= \sum_{i=1}^{m} \left( y_i - \sum_{k=1}^{n} c_k \phi_k(x_i) \right)^2
$$

以下では，この関数 $F(c_1,c_2,\ldots,c_n)$ を最小化するような $(c_1,c_2,\ldots,c_n)$ を求める問題を考える。

```{note}
与えられたデータの値 $y_k=f(x_k)$ とモデルから求めた値 $g(x_k)$ の差は，文脈によって誤差（error）や残差（residual error）とよばれるものである。
ここでは両者の違いには深く立ち入らない。
```

関数 $F$ の最小化問題を解くために，以下の定理を利用することができる。

```{prf:theorem}
:label: thm-minsq-dis

閉区間 $[a,b]$ 上に定義された1次独立な関数系 $\{\phi_k(x)\}_{k=1}^{n}$ を定義する。
また，実数 $c_k\ (k=1,2,\ldots,n)$ を用いて閉区間 $[a,b]$ 上の関数 $g(x)$ を $g(x)=\sum_{k=1}^{n} c_k \phi_k(x)$ と定義する。

2 つの $m$ 次元ベクトル $\bm{x},\bm{y}\in\setR^m$ が $i=1,2,\cdots,n$ について $x_i \in [a,b]$ を満たすとき，ベクトル $\bm{g}$ と $n$ 個のベクトル $\{ \bm{\phi}_k \}_{k=1}^{n}$ を以下のように定義する。

$$
& \bm{g} = (g(x_1),g(x_2),\ldots,g(x_m))
\\
& \bm{\phi}_k = (\phi_k(x_1),\phi_k(x_2),\ldots,\phi_k(x_m))
$$

このとき，ノルム $\|\bm{y}-\bm{g}\|_2$ を最小にする係数 $\{c_k\}$ は以下の連立1次方程式の解である。

$$
\begin{bmatrix}
  \naiseki{\bm\phi_1}{\bm\phi_1} & \naiseki{\bm\phi_1}{\bm\phi_2} & \cdots & \naiseki{\bm\phi_1}{\bm\phi_n} \\
  \naiseki{\bm\phi_2}{\bm\phi_1} & \naiseki{\bm\phi_2}{\bm\phi_2} & \cdots & \naiseki{\bm\phi_2}{\bm\phi_n} \\
  \vdots & \vdots & \ddots & \vdots \\
  \naiseki{\bm\phi_n}{\bm\phi_1} & \naiseki{\bm\phi_n}{\bm\phi_2} & \cdots & \naiseki{\bm\phi_n}{\bm\phi_n}
\end{bmatrix}
\begin{bmatrix}
  c_1 \\ c_2 \\ \vdots \\ c_n
\end{bmatrix}
=
\begin{bmatrix}
  \naiseki{\bm{y}}{\bm\phi_1} \\ \naiseki{\bm{y}}{\bm\phi_2} \\ \vdots \\ \naiseki{\bm{y}}{\bm\phi_n}
\end{bmatrix}
$$

ただしベクトル $\bm{a}$，$\bm{b}$ に対して $\naiseki{\bm{a}}{\bm{b}}$ は内積を表す。
```

証明は，関数の最小二乗近似の場合と同様であるため，省略する。

入力データが観測誤差などの誤差を含み，その誤差がたとえば正規分布のような形で発生するという前提のもとでは，最小二乗近似は最尤推定と一致する。
このことから，回帰分析のように，データに対してモデルをあてはめる際にモデルのパラメータを推定するための方法として最小二乗近似は有力な方法の一つであるといえる。

## 例（離散データの最小二乗近似）

```{list-table}
:header-rows: 1

* - $k$
  - $x_k$
  - $y_k$
* - 1
  - 0.00
  - 1.8
* - 2
  - 0.25
  - 1.4
* - 3
  - 0.50
  - 1.2
* - 4
  - 0.75
  - 2.0
* - 5
  - 1.00
  - 1.6
```

```{code-cell}
:tags: [remove-input]

import matplotlib.pyplot as plt

X = [0.00, 0.25, 0.50, 0.75, 1.00]
Y = [1.8, 1.4, 1.2, 2.0, 1.6]
plt.scatter(X, Y, color='blue')
plt.tick_params(labelsize=18)
plt.grid(True)
plt.show()
```

まず，これらの点を1次式で最小二乗近似することを考える。
そのためには，関数系 $\{1,x\}$ を考えればよい。
関数 $g(x)$ は $g(x) = c_1 + c_2 x$ となる。

このとき，ベクトル $\bm{g}$，$\bm{\phi}_1$，$\bm{\phi}_2$ は

$$
& \bm{g} = (g(x_1),g(x_2),\cdots,g(x_5))
\\
& \bm{\phi}_1 = (1,1,1,1,1)
\\
& \bm{\phi}_2 = (0.00,0.25,0.50,0.75,1.00)
$$

よって

$$
& \naiseki{\bm\phi_1}{\bm\phi_1} = 1^2 + 1^2 + 1^2 + 1^2 + 1^2 = 5
\\
& \naiseki{\bm\phi_1}{\bm\phi_2} = 0.00 + 0.25 + 0.50 + 0.75 + 1.00 = 2.50
\\
& \naiseki{\bm\phi_2}{\bm\phi_2} = (0/4)^2 + (1/4)^2 + (2/4)^2 + (3/4)^2 + (4/4)^2 = 1.875
\\
& \naiseki{\bm{y}}{\bm\phi_1} = 1.8 + 1.4 + 1.2 + 2.0 + 1.6 = 8.0
\\
& \naiseki{\bm{y}}{\bm\phi_2} = 1.8 \cdot 0.00 + 1.4 \cdot 0.25 + 1.2 \cdot 0.50 + 2.0 \cdot 0.75 + 1.6 \cdot 1.00 = 4.05
$$

したがって，連立1次方程式は

$$
&
\begin{bmatrix}
  \naiseki{\bm\phi_1}{\bm\phi_1} & \naiseki{\bm\phi_1}{\bm\phi_2} \\
  \naiseki{\bm\phi_2}{\bm\phi_1} & \naiseki{\bm\phi_2}{\bm\phi_2}
\end{bmatrix}
\begin{bmatrix}
  c_1 \\ c_2
\end{bmatrix}
=
\begin{bmatrix}
  \naiseki{\bm{y}}{\bm\phi_1} \\ \naiseki{\bm{y}}{\bm\phi_2}
\end{bmatrix}
\\&\iff
\begin{bmatrix}
  5 & 2.5 \\
  2.5 & 1.875
\end{bmatrix}
\begin{bmatrix}
  c_1 \\ c_2
\end{bmatrix}
=
\begin{bmatrix}
  8 \\ 4.05
\end{bmatrix}
$$

よって，$(c_1,c_2)=(1.56,0.08)$ となる。

入力された点および最小二乗近似によって得られた関数 $g(x)=c_1+c_2(x)$ を図示すると，以下のようになる。

```{code-cell}
:tags: [remove-input]

import matplotlib.pyplot as plt

X = [0.00, 0.25, 0.50, 0.75, 1.00]
Y = [1.8, 1.4, 1.2, 2.0, 1.6]
plt.scatter(X, Y, color='blue')
X = np.linspace(0.0, 1.0, 101)
Y = 1.56 + 0.08 * X
plt.plot(X, Y, color='red')
plt.tick_params(labelsize=18)
plt.grid(True)
plt.show()
```

係数 $\bm{c}$ を手計算でなくプログラムによって求める例は，以下のようになる（当然ながら，手計算で求めた場合と結果が一致する）。

```{code-cell}
import matplotlib.pyplot as plt
import numpy as np

def plot_least_squares(X, Y, N):
    # 入力データをプロットする
    plt.scatter(X, Y, color='blue')

    # 近似に用いる関数系（多項式 x^0, x^1, ..., x^{n}）
    vphi = [X ** k for k in range(N+1)]

    # 方程式の係数を求める
    A = np.array([[vphi[i] @ vphi[j] for j in range(N+1)] for i in range(N+1)])
    b = np.array([Y @ vphi[i] for i in range(N+1)])

    # 方程式 Ac=b を c について解く
    c = np.linalg.solve(A, b)

    # 最小二乗近似によって求めた関数 g をプロットする
    X = np.linspace(0.0, 1.0, 101)
    Y = c @ [X ** k for k in range(N+1)]
    plt.plot(X, Y, color='red')

    # グラフのパラメータを調整する
    plt.tick_params(labelsize=18)
    plt.grid(True)

    # グラフを描画する
    plt.show()

# 入力データ（x座標，y座標）
X = np.array([0.00, 0.25, 0.50, 0.75, 1.00])
Y = np.array([1.8, 1.4, 1.2, 2.0, 1.6])

# 近似に用いる多項式の最大次数
N = 1

plot_least_squares(X, Y, N)
```

入力データを見る限り，近似に用いる関数として直線（1次関数）を用いるのは適切でないように思える。

2次関数によって近似した例は以下のようになる。

```{code-cell}
plot_least_squares(X, Y, 2)
```

同様に，3次関数，4次関数によって近似した例をそれぞれ示す。

```{code-cell}
plot_least_squares(X, Y, 3)
```

```{code-cell}
plot_least_squares(X, Y, 4)
```

与えられた点が 5 点であるから，4次関数による最小二乗近似は誤差を 0 にすることができ，補間と同じ結果が得られる。

最小二乗近似で用いる関数系は，多項式関数 $x^0,x^1,x^2,\cdots$ 以外にも，指数関数・対数関数・三角関数などさまざまなものが利用できる。
かなり自由度が高いが，同様のアルゴリズムによって最小二乗近似を実行することができる。

自由度の高さが問題になる場合も考えられる。
関数系を多項式に限った場合でも，実際のデータに対して最小二乗近似によってフィッティングをおこなう際には最大次数をどのように決定するかが問題になる。
次数を決定するための指標として，たとえば，赤池情報量規準（AIC; An Information Criterion; Akaike's Information Criterion）などのモデル選択規準が用いられる。
詳しくは，データサイエンス系の講義科目（定量的データ分析など）を参照してほしい。
