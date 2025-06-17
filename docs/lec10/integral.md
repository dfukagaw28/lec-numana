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

# 数値積分

一般には，与えられた関数の積分を解析的に求めるのは難しい場合がある。
一方，数値的に近似値を求めるのであれば可能な場合も多い。
古くは，積分の概念が発見される以前に，アルキメデスが取り尽くし法を用いることによって図形の面積を数値的に求めたと言われている。

## 台形公式

閉区間 $[a,b]$ 上で定義された関数 $f(x)$ の定積分は，曲線 $y=f(x)$ と $x$ 軸，および直線 $x=a$，$x=b$ で囲まれた部分の面積である（{numref}`fig-integral-1`）。ただし，$a \leq x \leq b$ において $f(x)>0$ とする。

```{code-cell}
:tags: [remove-input, remove-output]

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from myst_nb import glue
import numpy as np
from scipy.interpolate import lagrange

def make_figure(f, a, b, approx=False, d=1, n=1, vline=False):
    # 関数 f(X) の各点を数値的に求める
    X = np.linspace(a-1, b+1, 1001)
    Y = f(X)

    # グラフを描画するための準備
    fig , ax = plt.subplots()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_position('zero')
    ax.set_yticks([])
    ax.set_xlim(min(X), max(X))
    ax.set_ylim(-1, max(Y) + 1)

    # グラフ上の点
    PX = np.linspace(a, b, d * n + 1)
    PY = f(PX)
    if d == 1 and n == 1:
        PL = ['$a$', '$b$']
    else:
        PL = [f'$x_{{{k}}}$' for k, xk in enumerate(PX)]

    # ラグランジュ補間による近似関数
    h = (b - a) / (n * d)
    if approx:
        gg = [lagrange(PX[d*k:d*k+d+1], PY[d*k:d*k+d+1]) for k in range(n)]
        
        g = lambda x: gg[np.minimum(n - 1, np.floor((x - a) / h))](x)

    # 領域を塗りつぶす
    if approx:
        for k in range(n):
            X2 = np.linspace(PX[d*k], PX[d*k+d], 1001)
            g = gg[k]
            ax.fill_between(X2, g(X2), 0, facecolor='#ffcccc')
    else:
        ax.fill_between(X, Y, 0, where=(X >= a) & (X <= b), facecolor='#ccccff')

    # 区間の区切りを表示する
    if vline:
        for k in range(n + 1):
            ax.axvline(PX[d*k], 0, 1, linestyle=':', color='gray')

    # グラフ y=f(x) を描画する
    if approx:
        ax.plot(X, Y, color='blue', alpha=0.5)
        for k in range(n):
            X2 = np.linspace(PX[d*k], PX[d*k+d], 1001)
            g = gg[k]
            ax.plot(X2, g(X2), color='red')
    else:
        ax.plot(X, Y, color='blue', alpha=1)

    # グラフ上の点を描画する
    ax.scatter(PX, PY, color='blue')
    ax.xaxis.set_major_locator(ticker.FixedLocator(PX))
    ax.set_xticks(PX)
    ax.set_xticklabels(PL, fontsize=16)

    return fig


# 積分範囲
a, b = 0, 20

# 関数 f を定義する
f = lambda x: 4 - 0.2 * x * np.sin(0.14 * x) + 0.3 * np.cos(0.9 * x)

glue("integral_1", make_figure(f, a, b), display=False)

glue("trapezoid_rule_1", make_figure(f, a, b, approx=True), display=False)

glue("trapezoid_rule_2", make_figure(f, a, b, approx=True, n=5), display=False)

glue("simpson_rule_1", make_figure(f, a, b, approx=True, d=2), display=False)

glue("simpson_rule_2", make_figure(f, a, b, approx=True, d=2, n=3, vline=True), display=False)
```

```{glue:figure} integral_1
:name: "fig-integral-1"

定積分の例
```

このとき，2 点 $(a,f(a))$，$(b,f(b))$ を通る直線 $y=g(x)$ によって関数 $f(x)$ を近似（線形補間）することによって，$f(x)$ の定積分を $g(x)$ の定積分で近似できる（{numref}`fig-trapezoid-rule-1`）。

```{glue:figure} trapezoid_rule_1
:name: "fig-trapezoid-rule-1"

台形による近似の例
```

関数 $g(x)$ の閉区間 $[a,b]$ における定積分は台形の面積に相当する（$a \leq x \leq b$ において $f(x) > 0$ より $g(x) > 0$ に注意）。
すなわち，関数 $f(x)$ の定積分を台形の面積によって以下のように近似できる。

$$
I = \int_a^b f(x) dx = \frac{(f(a) + f(b))(b - a)}{2} = \frac{h}{2}(y_0+y_1)
$$

ただし，$x_0=a$，$x_1=b$，$y_0=f(x_0)=f(a)$，$y_1=f(x_1)=f(b)$，$h=x_1-x_0=b-a$ とする。

このようにして定積分の値を近似的に求める方法を **台形公式**（trapezoidal rule）という。
台形則ともいう。

なお，関数 $f(x)$ が負の値をとる場合にも，上の式はそのまま利用できる（定積分の値は負の面積となる）。

積分区間 $[a,b]$ を $n$ 分割して区間ごとに台形公式を適用することもでき，そちらの方が近似精度が高い（{numref}`fig-trapezoid-rule-2`）。

```{glue:figure} trapezoid_rule_2
:name: "fig-trapezoid-rule-2"

複合台形公式
```

このようにして積分区間を複数の区間に分割し，区間ごとに台形公式を適用する方法を **複合台形公式**（composite trapezoidal rule）という。
一般に，台形公式といえば複合台形公式を指すことも多い。

まず，閉区間 $[a,b]$ 上に相異なる $n+1$ 個の値 $\{x_i\}_{i=0}^{n}$ を $a=x_0<x_1<x_2<\cdots<x_n=b$ を満たすようにとる。
ここでは，各区間が等分になるようにする。
すなわち，以下が成り立つように $\{x_k\}$ を定める。

$$
x_k = a + kh \quad (k=0,1,2,\cdots,n),
\quad
h = \frac{b - a}{n}
$$

このとき，定積分は以下のように近似できる。

$$
  & \int_a^b f(x) dx
    = \int_{x_{0}}^{x_{n}} f(x) dx
    = \sum_{k=1}^{n} \int_{x_{k-1}}^{x_{k}} f(x) dx
\\&\quad
    \approx \sum_{k=1}^{n} \frac{h}{2} (f(x_{k-1}) + f(x_{k}))
\\&\quad
    = \frac{h}{2} \left\{ (y_0 + y_1) + (y_1 + y_2) + \cdots + (y_{n-2} + y_{n-1}) + (y_{n-1} + y_{n}) \right\}
\\&\quad
    = \frac{h}{2} \left\{ y_0 + 2 y_1 + 2 y_2 + \cdots + 2 y_{n-1} + y_{n} \right\}
\\&\quad
    = h \left\{ \frac{y_0}{2} + y_1 + \cdots + y_{n-1} + \frac{y_n}{2} \right\}
$$

ただし，$y_k=f(x_k)$ とする。

## 例（台形公式）

たとえば，関数 $f(x)=\sqrt{1-x^2}$ を $0$ から $1$ まで定積分する場合を考える。

```{code-cell}
import matplotlib.pyplot as plt
import numpy as np

X = np.linspace(0, 1, 1001)
Y = np.sqrt(1 - X*X)

fig , ax = plt.subplots()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_position('zero')
ax.spines['bottom'].set_position('zero')
ax.set_aspect('equal')
ax.set_xlim(-0.1, 1.1)
ax.set_ylim(-0.1, 1.1)

ax.fill_between(X, Y, 0, facecolor='none', edgecolor='lightgray', hatch='//')

ax.plot(X, Y)

plt.show()
```

$y=f(x)$ のグラフは円弧であり，その定積分は中心角 90 度の扇形の面積にあたるため，面積は容易に求めることができ， $\pi \left( \frac{1}{2} \right)^2 = \frac{\pi}{4}$ となる。

あるいは，置換積分を用いて代数的に積分を求めることもできる。
たとえば，$x=\sin\theta$ とおくと $dx=\cos\theta d\theta$ であるから，

$$
\int_0^1 f(x) dx
&= \int_0^1 \sqrt{1-x^2} dx
\\
&= \int_{0}^{\pi/2} \sqrt{1-(\sin\theta)^2} \cos\theta d\theta
\\
&= \int_{0}^{\pi/2} \cos^2\theta d\theta
\\
&= \int_{0}^{\pi/2} \frac{1 + \cos 2\theta}{2} d\theta
\\
&= \left[ \frac{1}{2} \theta + \frac{1}{4} \sin 2\theta \right]_{0}^{\pi/2}
\\
&= \frac{1}{2}\frac{\pi}{2} + \frac{1}{4} \sin\pi
\\
&= \frac{\pi}{4}
$$

となる。

あるいは，SymPy を用いて解くこともできる。

```{code-cell}
import sympy as sym
x = sym.sympify('x')
f = sym.sympify('sqrt(1-x**2)')
lhs = sym.Integral(f, (x, 0, 1))
rhs = sym.integrate(f, (x, 0, 1)).factor()
sym.Eq(lhs, rhs)
```

以下では，この定積分を台形公式によって求め，値が近似的に一致することを確認する。

既に見たように，厳密な値は $\pi/4$ である。
Python では math パッケージや numpy パッケージにおいて定数 `pi` が定義されており，円周率の近似値を利用することができる[^precision-of-pi]。

```{code-cell}
import math
print(math.pi / 4)
```

```{code-cell}
import numpy as np
print(np.pi / 4)
```

[^precision-of-pi]: なお，`math.pi` と `numpy.pi` の型はいずれも float であるため，その精度は 10 進数で 15～16 桁程度である（IEEE754 倍精度浮動小数点数の仮数部は 52+1 bit であり， $\log_{10} 2^{52+1} \approx 15.95$ である）。
円周率 $\pi$ の（近似）値は，
標準パッケージ math では [ここ](https://github.com/python/cpython/blob/v3.13.5/Include/pymath.h#L14)（CPython v3.13.5），
NumPy では [ここ](https://github.com/numpy/numpy/blob/v2.3.0/numpy/_core/include/numpy/npy_math.h#L77) （NumPy v2.3.0）でそれぞれ定義されている。
桁数の違いに驚くかもしれないが，Python では現在のところ `3.141592653589793238462643383279502884 == 3.141592653589793` が `True` になることをふまえると納得できるかもしれない。

台形公式によって定積分を求める関数の実装例を以下に示す。

```{code-cell}
from typing import Callable
import numpy as np

def trapezoid(f: Callable[[float], float], a: float, b: float, n: int):
    """
    複合台形公式を用いて定積分を求める

    f: 1変数関数
    a: 積分区間の始点
    b: 積分区間の終点
    n: 積分区間の分割数（=分点数-1）
    """

    assert n > 0
    assert a < b

    x = np.linspace(a, b, n + 1)
    y = f(x)
    h = (b - a) / n

    return h * ((y[0] + y[-1]) / 2 + np.sum(y[1:-1]))
```

まず，$n=1$ すなわち全体を 1 つの台形とみなす場合を考える。

```{code-cell}
import numpy as np

## 関数 f(x) を定義する
f = lambda x: np.sqrt(1 - x ** 2)
a, b = 0.0, 1.0
n = 1

## 台形公式を用いて定積分を求める
trapezoid(f, a, b, n)
```

次に，分割数 $n$ を $1$ から $10$ まで順に増加させながら，台形公式の精度を確認する。

```{code-cell}
for n in range(1, 11):
    s = trapezoid(f, a, b, n)
    print(n, s)
```

10分割した場合でも，小数点以下 2 桁目で誤差が出てしまっており，精度は高くないことが分かる。

10万分割の例も示す。

```{code-cell}
print('近似値:', trapezoid(f, 0, 1, 100_000))
print('真の値:', np.pi / 4)
```

真の値と比較することにより，小数点以下 7 桁目まで正確な値になっていることが分かる。

## シンプソン公式

引き続き，閉区間 $[a,b]$ 上で定義された関数 $f(x)$ の定積分 $I=\int_a^b f(x) dx$ について考える。

台形公式では $f(x)$ を 1 次関数で近似（補間）した。
3点を通る2次関数によって近似することによって近似精度の向上が期待できる。

```{glue:figure} simpson_rule_1
:name: "fig-simpson-rule-1"

シンプソン公式の例
```

この手法を **シンプソン公式**（Simpson's rule）という[^simpson]。

[^simpson]: シンプソン公式の名前は英国の数学者 Thomas Simpson（1710-1761）に由来する。ただし，ドイツの天文学者ヨハネス・ケプラー（Johannes Kepler, 1571--1630）も同じ公式を先に発見していたとされており，特にドイツではケプラー公式（Keplersche Fassregel）とも呼ばれる。

関数 $f(x)$ を2次関数で近似するためには，3点を通るラグランジュ補間多項式を求めればよい。

$$
x_0 = a, \quad
x_1 = a + h, \quad
x_2 = a + 2h = b, \quad
h = \frac{b-a}{2}
$$

とおけば，ラグランジュの基本多項式 $\ell_k(x)$ は

$$
\ell_0(x) = \frac{(x-x_1)(x-x_2)}{(x_0-x_1)(x_0-x_2)}
\\
\ell_1(x) = \frac{(x-x_2)(x-x_0)}{(x_1-x_2)(x_1-x_0)}
\\
\ell_2(x) = \frac{(x-x_0)(x-x_1)}{(x_2-x_0)(x_2-x_1)}
$$

となり，$f(x)$ は以下のように近似できる。

$$
f(x) \approx \ell_0(x) \cdot y_0 + \ell_1(x) \cdot y_1 + \ell_2(x) \cdot y_2
$$

ここで，区間 $[a,b]$ における $\ell_k(x)$ の定積分 $\alpha_k = \int_a^b \ell_k(x)dx$ を求める。

$$
   \alpha_0
   &= \int_a^b \ell_0(x) dx
\\ &= \int_a^b \frac{(x-x_1)(x-x_2)}{(x_0-x_1)(x_0-x_2)} dx
\\ &= \int_a^{a+2h} \frac{(x-a-h)(x-a-2h)}{(-h)(-2h)} dx
\\ &= \int_0^{2h} \frac{(t-h)(t-2h)}{2h^2} dt \quad\text{置換積分 $(t=x-a)$}
\\ &= \int_0^{2} \frac{(s-1)(s-2)}{2} h ds \quad\text{置換積分 $(s=t/h)$}
\\ &= \frac{h}{2} \int_0^{2} \{ (s-1)(s-2) \} ds
\\ &= \frac{h}{2} \int_0^{2} \{ s^2 - 3 s + 2 \} ds
\\ &= \frac{h}{2} \left[ \frac{1}{3} s^3 - \frac{3}{2} s^2 + 2 s \right]_0^{2}
\\ &= \frac{h}{2} \left( \frac{1}{3} \cdot 2^3 - \frac{3}{2} \cdot 2^2 + 2 \cdot 2 \right)
\\ &= \frac{h}{2} \left( \frac{8}{3} - 6 + 4 \right)
\\ &= \frac{h}{2} \cdot \frac{2}{3}
\\ &= \frac{1}{3} h
$$

同様にして

$$
   \alpha_1
   &= \int_a^b \ell_1(x) dx
\\ &= \int_a^b \frac{(x-x_2)(x-x_0)}{(x_1-x_2)(x_1-x_0)} dx
\\ &= \int_a^{a+2h} \frac{(x-a-2h)(x-a)}{(-h)(+h)} dx
\\ &= \int_0^{2} \frac{(s-2)(s)}{-1} h ds \quad\text{置換積分 $(s=(x-a)/h)$}
\\ &= h \int_0^{2} \{ 2 s - s^2 \} ds
\\ &= h [ s^2 - \frac{1}{3} s^3 ]_0^2
\\ &= \frac{4}{3} h
\\
   \alpha_2
   &= \int_a^b \ell_2(x) dx
    = \int_a^b \frac{(x-x_0)(x-x_1)}{(x_2-x_0)(x_2-x_1)} dx
    = \frac{1}{3} h
$$

以上より，定積分 $I$ は以下のように近似できる。

$$
I &= \int_a^b f(x) dx
\\
  &\approx \sum_{k=0}^2 \int_a^b \ell_k(x) \cdot y_k dx
\\
  &= \sum_{k=0}^2 \left( \int_a^b \ell_k(x) dx \right) y_k
\\
  &= \sum_{k=0}^2 \alpha_k \cdot y_k
\\
  &= \alpha_0 \cdot y_0 + \alpha_1 \cdot y_1 + \alpha_2 \cdot y_2
\\
  &= \frac{h}{3} (y_0 + 4 y_1 + y_2)
$$

台形公式と同様に，区間をさらに細かく分割し，連続する3点ごとにシンプソン公式を適用することで，近似精度を向上させることができる。

```{glue:figure} simpson_rule_2
:name: "fig-simpson-rule-2"

複合シンプソン公式の例
```

具体的には，区間 $[a,b]$ を $2n$ 等分して

$$
x_k = a + kh,
\quad
y_k = f(x_k),
\quad
(0 \leq k \leq 2n)
\quad
h = \frac{b - a}{2n}
$$

とおく。このとき，定積分は以下のように近似できる。

$$
  & \int_a^b f(x) dx
    = \sum_{k=1}^{n} \int_{x_{2k-2}}^{x_{2k}} f(x) dx
\\&\quad
    \approx \sum_{k=1}^{n} \left\{ \frac{h}{3} (y_{2k-2} + 4 y_{2k-1} + y_{2k}) \right\}
\\&\quad
    = \frac{h}{3} \left\{ y_0 + 4 y_1 + 2 y_2 + 4 y_3 + 2 y_4 + \cdots + 2 y_{2n-2} + 4 y_{2n-1} + y_{2n} \right\}
\\&\quad
    = \frac{h}{3} (y_{0} + y_{2n}) + \frac{4h}{3} \sum_{k=1}^{n} y_{2k-1} + \frac{2h}{3} \sum_{k=1}^{n-1} y_{2k}
$$

これを **複合シンプソン公式** という。

一般に，単にシンプソン公式といえば複合シンプソン公式を指すことも多い。
以下では，複合シンプソン公式

## 例（シンプソン公式）

さきほどと同様に，区間 $[0,1]$ 上の関数 $f(x)=\sqrt{1-x^2}$ の定積分をシンプソン公式を用いて求めてみよう。

```{code-cell}
import numpy as np

def simpson(f: Callable[[float], float], a: float, b: float, n: int):
    """
    複合シンプソン公式を用いて定積分を求める

    f: 1変数関数
    a: 積分区間の始点
    b: 積分区間の終点
    n: 積分区間の分割数（=(分点数-1)/2）
    """

    assert n > 0
    assert a < b

    x = np.linspace(a, b, 2 * n + 1)
    y = f(x)
    h = (b - a) / (2 * n)

    # c = [1, 4, 2, 4, 2, ..., 2, 4, 1]
    c = np.tile([2., 4.], n + 1)[:-1]
    c[0] = c[-1] = 1.

    return (h / 3) * (c @ y)
```

積分区間の分割数 $n$ を $1$ から $10$ まで順に増加させながら，台形公式とシンプソン公式の精度を比較する。

```{code-cell}
import numpy as np

f = lambda x: np.sqrt(1 - x ** 2)
a, b = 0, 1

print('積分区間の分割数ごとの比較')
print('trapezoid() vs simpson()')
for n in range(1, 11):
    s1 = trapezoid(f, a, b, n)
    s2 = simpson(f, a, b, n)
    print(f'{n:3d} {s1:.16f} {s2:.16f}')
print('正解:', np.pi / 4)
```

このように，台形公式よりシンプソン公式の方が精度が高いことが分かる。

上の例では分点数（つまりサンプル点の数）の個数が異なるものを比較している（積分区間の分割数を $n$ とするとき，台形公式では $n+1$ 個に対してシンプソン公式では $2n+1$ 個である）。
そのため，フェアな比較になっていない可能性がある。
サンプル点の個数を揃えて比較してみよう。

```{code-cell}
import numpy as np

f = lambda x: np.sqrt(1 - x ** 2)
a, b = 0, 1

print('分点数ごとの比較')
print('trapezoid() vs simpson()')
for n in range(1, 11):
    s1 = trapezoid(f, a, b, 2 * n)
    s2 = simpson(f, a, b, n)
    print(f'{2*n+1:3d} {s1:.16f} {s2:.16f}')
print('正解:', np.pi / 4)
```

サンプル点の個数を揃えても，やはり，シンプソン公式の方が精度が高いことが分かる。

最後に，1万分割，10万分割の例も示す。

```{code-cell}
print(simpson(f, 0, 1, 10_000))
print(simpson(f, 0, 1, 100_000))
```

## ニュートン・コーツ公式

台形公式やシンプソン公式は，いくつかの共通点を持っている。

- 関数 $f(x)$ の定積分を求める代わりに近似関数の定積分を求める。
- 積分区間を小区間に分割し，区間ごとにラグランジュ補間多項式を求め，近似関数とする。

ラグランジュ補間多項式を求める際のサンプル点の個数が $n=2$ の場合（線形近似）が台形公式に相当し，サンプル点の個数が $n=3$ の場合（2次近似）が台形公式に相当する。

一般の $d$ について，$d+1$ 個のサンプル点からラグランジュ補間多項式を $d$ 次多項式として求め，$f(x)$ の定積分を求める方法を考えてみよう。

まず，サンプル点を以下のように定義する。

$$
x_k = a + kh, \quad
y_k = f(x_k) \quad (k = 0, 1, 2, \cdots, d), \quad
h = \frac{b-a}{d}
$$

ラグランジュの基本多項式 $\ell_k(x)$ は

$$
\ell_k(x) = \frac{(x-x_0)(x-x_1)\cdots(x-x_{k-1})(x-x_{k+1})\cdots(x-x_{d-1})(x-x_d)}
                 {(x_k-x_0)(x_k-x_1)\cdots(x_k-x_{k-1})(x_k-x_{k+1})\cdots(x_k-x_{d-1})(x_k-x_d)}
$$

となる。
定積分 $I=\int_a^b f(x)dx$ の値は，以下のように近似できる。

$$
I &= \int_a^b f(x) dx
\\ &\approx \int_a^b \left\{ \sum_{k=0}^{d} \ell_k(x) \cdot y_k \right\} dx
\\ &= \sum_{k=0}^{d} \left\{ \int_a^b \ell_k(x) dx \right\} \cdot y_k
$$

これを **$d$ 次のニュートン・コーツ公式**（Newton-Cotes formula of degree $d$）とよぶ。

台形公式は $1$ 次のニュートン・コーツ公式であり，シンプソン公式は $2$ 次のニュートン・コーツ公式である。

## NumPy, SciPy を用いた数値積分

台形公式は，NumPy の [`numpy.trapezoid()` 関数](https://numpy.org/doc/stable/reference/generated/numpy.trapezoid.html) を用いて計算することもできる[^numpy-trapezoid]。
結果は同一になる（はず）。
引数の与え方の違いに注意が必要。

[^numpy-trapezoid]: `numpy.trapezoid()` 関数は，NumPy 2.0 で追加された。NumPy 1.x 環境では `numpy.trapz()` という名前であった。

```{code-cell}
import numpy as np

f = lambda x: np.sqrt(1 - x ** 2)
a, b = 0, 1

print('trapezoid() vs numpy.trapezoid()')
for n in range(1, 11):
    s1 = trapezoid(f, a, b, n)
    x = np.linspace(a, b, n + 1)
    s2 = np.trapezoid(f(x), x)
    print(f'{n:3d} {s1:.16f} {s2:.16f}')
```

シンプソン公式は，SciPy で利用可能である。
[`scipy.integrate.simpson()` 関数](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.simpson.html) を用いて計算することもできる。
こちらも引数の与え方の違いに注意が必要。結果は変わらないはずであるが，わずかに（小数点以下16桁目くらいに）誤差が生じている。

```{code-cell}
import numpy as np
from scipy import integrate

f = lambda x: np.sqrt(1 - x ** 2)
a, b = 0, 1

print('simpson() vs scipy.integrate.simpson()')
for n in range(1, 11):
    s1 = simpson(f, a, b, n)
    x = np.linspace(a, b, 2 * n + 1)
    s2 = integrate.simpson(f(x), x=x)
    print(f'{n:3d} {s1:.16f} {s2:.16f}')
```

今回紹介した他にも，SciPy で数値積分をおこなうさまざまな関数が用意されている。
たとえば，ラグランジュ補間の代わりにエルミート補間（分割した区間の境界における微分係数を利用して補正をおこなう）を用いるロンバーグ積分（Romberg）や，係数だけでなく分点 $\{x_k\}$ を適応的に求めることで精度を向上させるガウス型積分公式などがよく知られている。
特にガウス型積分公式のような適応型の（adaptive）数値積分法は精度が高い。SciPy の [`quad()` 関数](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html#scipy.integrate.quad) のように適応型の数値積分ルーチンを呼び出す際には，引数として，分点における関数値の配列でなく関数オブジェクトとして与えるように設計されていることが多い。
詳しくは，[SciPy のドキュメント「Integration and ODEs」](https://docs.scipy.org/doc/scipy/reference/integrate.html)を参照されたい。
