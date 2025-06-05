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

# 関数近似と補間 (1) 関数の補間

ある１日の気温の変化を調べて，午前 9 時の気温が 21 °C ，正午の気温が 24 °C であったとする。

```{list-table}
:header-rows: 1

* - 時刻
  - 気温（°C）
* - 午前9時
  - 21
* - 午前10時
  - ？
* - 午前11時
  - ？
* - 午前12時
  - 24
```

このとき，午前 10 時の気温については，2つの時刻において観測した気温からおおよその値を推測できる。
たとえば，午前 10 時の気温がおよそ 22 °C であったと推測するのは妥当だろう。

```{code-cell}
:tags: [remove-input]

import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib_fontja

X = [dt.datetime(2022, 1, 1, h) for h in [9, 12]]
Y = [21, 24]
plt.plot(X, Y, linestyle=':', color='blue')
plt.scatter(X, Y, color='red')
plt.xlabel('時刻（時）', fontsize=18)
plt.ylabel('気温（°C）', fontsize=18)
plt.tick_params(labelsize=18)
plt.gca().xaxis.set_major_locator(mdates.HourLocator())
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.grid(True)
plt.show()
```

気温だけでなく一般のデータを考えるとき，時刻 $x=x_1$ における変量 $y$ の値が $y=y_1$，時刻 $x=x_2$ における変量 $y$ の値が $y=y_2$，であれば，座標平面上に与えられた $2$ 点 $(x_1,y_1)$，$(x_2,y_2)$ を結ぶ線分を描くことによって区間 $[x_1, x_2]$ における変量 $y$ の値を推定することは，やはり妥当といえよう。

$2$ 点でなく $3$ 個以上の点が与えられているときはどうだろうか。変数 $x$ の値が小さい順に並べたうえで隣り合う点どうしを線分で結ぶことによって，与えられた点をすべて通るような，区分的に線形な関数を考えることができる。
この場合，グラフは折れ線になる。

```{code-cell}
:tags: [remove-input]

import matplotlib.pyplot as plt

X = [10, 20, 30, 50]
Y = [100, 300, 200, 300]
plt.plot(X, Y, linestyle=':', color='blue')
plt.scatter(X, Y, color='red')
plt.tick_params(labelsize=18)
plt.gca().set_ylim(0, 350)
plt.grid(True)
plt.show()
```

あるいは，すべての点を通るような曲線を考えることもできる。

```{code-cell}
:tags: [remove-input, remove-stderr]

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange

X = [10, 20, 30, 50]
Y = [100, 300, 200, 300]
K = [3, 4, 5]
C = ['blue', 'green', 'orange']

for k in range(3):
    f = np.poly1d(np.polyfit(X, Y, K[k]))
    x = np.linspace(min(X), max(X), 101)
    y = f(x)
    plt.plot(x, y, linestyle=':', color=C[k])

plt.scatter(X, Y, color='red')
plt.tick_params(labelsize=18)
plt.gca().set_ylim(0, 350)
plt.grid(True)
plt.show()
```

一般に，未知の関数 $f:\mathbb{R}\to\mathbb{R}$ について，単調増加する $n$ 個の実数列 $x_0,x_1,\cdots,x_n$（$x_0 < x_1 < \cdots < x_n$）およびそれらに対する関数値 $\{y_k=f(x_k)\}$，$k=0,1,\cdots,n$ が与えられたする。
このとき，座標平面上の $n+1$ 個の点 $(x_0,y_0),(x_1,y_1),\cdots,(x_n,y_n)$ をすべて通るような曲線 $y=f(x)$ を定めることによって区間 $[x_0, x_n]$ における任意の $x$ に対する関数値 $f(x)$ を推測することを，関数の **補間**（interpolation）または **内挿** という [^extrapolation]。

[^extrapolation]: 推測すべき点が観測データの範囲（$x_0 \leq x \leq x_n$）におさまる場合を補間あるいは内挿ととよび，推測すべき点が観測データの範囲外（$x < x_0$ や $x > x_n$）にある場合を補外あるいは外挿とよんで区別する場合がある。
また，「補完」でなく「補間」である点も注意が必要。

なお，座標平面上に与えられた $(n+1)$ 個の点をすべて通る曲線 $y=f(x)$ が多項式関数 $f(x)$ によって与えられるとすると，$f(x)$ は高々 $n$ 次式である。 
また，そのような $n$ 次多項式は一意に定まる。

```{prf:theorem}
:label: thm-interp

座標平面上に $n+1$ 個の点 $(x_0,y_0),(x_1,y_1),\cdots,(x_n,y_n)$ が与えられる。
ただし， $x_0,x_1,\cdots,x_n$ はすべて異なるとする。
これら $n$ 個の点をすべて通るような $n$ 次関数は一意に定まる。
```

```{prf:proof}
$a_0,a_1,\cdots,a_n$ をパラメータとして $n$ 次関数を以下のように表すことができる。

$$
y = a_0 + a_1 x + a_2 x^2 + \cdots + a_n x^n
$$

$n+1$ 個の点列 $\{(x_k,y_k)\}$ を通ることから，以下を得る。

$$
y_0 &= a_0 + a_1 x_0 + a_2 x_0^2 + \cdots + a_n x_0^n
\\
y_1 &= a_0 + a_1 x_1 + a_2 x_1^2 + \cdots + a_n x_1^n
\\
y_2 &= a_0 + a_1 x_2 + a_2 x_2^2 + \cdots + a_n x_2^n
\\
&\vdots
\\
y_n &= a_0 + a_1 x_n + a_2 x_n^2 + \cdots + a_n x_n^n
$$

行列を用いて以下のように書くことができる。

$$
\begin{bmatrix}
1 & x_0 & x_0^2 & \cdots & x_0^n \\
1 & x_1 & x_1^2 & \cdots & x_1^n \\
1 & x_2 & x_2^2 & \cdots & x_2^n \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
1 & x_n & x_n^2 & \cdots & x_n^n \\
\end{bmatrix}
\begin{bmatrix}
a_0 \\ a_1 \\ \vdots \\ a_n
\end{bmatrix}
=
\begin{bmatrix}
y_0 \\ y_1 \\ \vdots \\ y_n
\end{bmatrix}
$$

左辺の行列 $\boldsymbol{V}$ は **ヴァンデルモンド行列**（Vandermonde matrix）とよばれる。
ヴァンデルモンド行列の行列式 $\det\boldsymbol{V}$ を **ヴァンデルモンドの行列式**（Vandermonde's determinant）とよぶ。
この場合，ヴァンデルモンドの行列式は

$$
\det \boldsymbol{V}
= \prod_{i<j} (x_j - x_i)
= \prod_{i=1}^{n-1} \prod_{j=i+1}^{n} (x_j - x_i)
$$

となる。
$x_k$ は互いに異なるため，この行列式は $0$ ではない。
つまり，ヴァンデルモンド行列は正則であり，逆行列を持つ。
したがって，上の方程式は唯一の解を持つ。
```

$n$ が小さい場合について考えてみよう。

```{prf:example} 2点を通る直線（初等的解法）
:label: example-interp-2

定数 $a, b$ によって表される直線 $\ell:y=ax+b$ が $2$ 点 $(x_1,y_1)$，$(x_2,y_2)$ を通るとき，$a, b$ の値は連立方程式

$$
  y_1 = a x_1 + b
  \\
  y_2 = a x_2 + b
$$

を解くことによって得られる（ただし $x_1 \neq x_2$ とする）。

$$
a &= \frac{y_1 - y_2}{x_1 - x_2}
\\
b &= \frac{x_1 y_2 - x_2 y_1}{x_1 - x_2}
$$

すなわち，$2$ 点を通る直線の式は

$$
y = \frac{y_1 - y_2}{x_1 - x_2} x + \frac{x_1 y_2 - x_2 y_1}{x_1 - x_2}
$$

となる。
```

## ラグランジュ補間

ラグランジュ補間は，与えられた $n+1$ 個の点を通るような $n$ 次多項式関数による補間である。
そのような多項式は，{prf:ref}`example-interp-2` のように単純に連立方程式を解くことによっても得られるが，与えられた点の座標から直接計算することもできる。
以下では与えられた点の個数が 2 個，3 個の例を確認し，最後に一般化する。

```{prf:example} 2点を通る直線（ラグランジュ補間）
:label: example-interp-2-2

{prf:ref}`example-interp-2` と同じ問題を，別の方法で解いてみよう。

まず，$1$次関数 $h_1(x)$，$h_2(x)$ を用いて，求める直線 $\ell$ の式を以下のように表す。

$$
y = h_1(x) \cdot y_1 + h_2(x) \cdot y_2
$$

直線 $\ell$ が $2$ つの点 $(x_1,y_1)$，$(x_2,y_2)$ を通るためには，関数 $h_k(x)$，$k=1,2$ が以下の条件を満たせば十分である。

$$
h_1(x_1) = 1
\\
h_1(x_2) = 0
\\
h_2(x_1) = 0
\\
h_2(x_2) = 1
$$

これらの条件を満たす $1$次関数 $h_1(x)$，$h_2(x)$ は一意に定まり，それぞれ以下のようになる。

$$
h_1(x) = \frac{x - x_2}{x_1 - x_2}
\\
h_2(x) = \frac{x - x_1}{x_2 - x_1}
$$

したがって，求める直線 $\ell$ の方程式は

$$
y = \frac{x - x_2}{x_1 - x_2} y_1 + \frac{x - x_1}{x_2 - x_1} y_2
$$

これを変形すると，{prf:ref}`example-interp-2` で求めたものと同じ式となる。

これは直線 $\ell$ が $2$ 点を通るための十分条件であるが，直線の一意性から，これが唯一の直線であることが分かる。
```

```{prf:example} 3点を通る直線（ラグランジュ補間）
:label: example-interp-2-3

同様に，$3$ 点 $(x_1,y_1),(x_2,y_2),(x_3,y_3)$ を通る曲線 $C$ を考えよう。
曲線 $C$ は，$2$ 次関数 $f(x)$ を用いて $y=f(x)$ と表すことができる。
$f(x)$ は，$2$次関数 $\ell_1(x)$，$\ell_2(x)$，$\ell_3(x)$ を用いて以下のように表すことができる。

$$
f(x) = \ell_1(x) \cdot y_1 + \ell_2(x) \cdot y_2 + \ell_3(x) \cdot y_3
$$

曲線 $C$ が $3$ つの点 $(x_1,y_1)$，$(x_2,y_2)$，$(x_3,y_3)$ を通るためには，関数 $\ell_k(x)$，$k=1,2,3$ が以下の条件を満たせば十分である。

$$
\ell_k(x_i)
= \delta_{ki}
= \begin{cases} 1 & \text{$k = i$ のとき} \\ 0 & \text{$k \neq i$ のとき} \end{cases}
$$

ただし，$\delta_{ki}$ はクロネッカーのデルタ（Kronecker delta）である。
これらの条件を満たすには，$2$次関数 $\ell_k(x)$ をそれぞれ以下のように定義すればよい。

$$
\ell_k(x) = \frac{\prod_{i \neq k}(x - x_i)}{\prod_{i \neq k}(x_k - x_i)}
$$

$k=1,2,3$ について書き下すと，

$$
\ell_1(x) = \frac{(x - x_2)(x - x_3)}{(x_1 - x_2)(x_1 - x_3)}
\\
\ell_2(x) = \frac{(x - x_1)(x - x_3)}{(x_2 - x_1)(x_2 - x_3)}
\\
\ell_3(x) = \frac{(x - x_1)(x - x_2)}{(x_3 - x_1)(x_3 - x_2)}
$$

したがって，求める曲線 $C$ の方程式は以下のようになる。

$$
f(x)
&= \sum_{k=1}^{n} \frac{\prod_{i \neq k}(x - x_i)}{\prod_{i \neq k}(x_k - x_i)} \cdot y_k
\\
&= \frac{(x - x_2)(x - x_3)}{(x_1 - x_2)(x_1 - x_3)} \cdot y_1
 + \frac{(x - x_1)(x - x_3)}{(x_2 - x_1)(x_2 - x_3)} \cdot y_2
 + \frac{(x - x_1)(x - x_2)}{(x_3 - x_1)(x_3 - x_2)} \cdot y_3
$$
```

この例のようにして得られる多項式を **ラグランジュ型の補間多項式**（interpolation polynomial in the Lagrange form）あるいは単に **ラグランジュ多項式**（Lagrange polynomial）とよぶ。
また，ラグランジュ多項式で用いる多項式 $\ell_k(x)$ を **基本多項式**（basis polynomial） とよぶ。

## 例（ラグランジュ補間）

ラグランジュ補間をおこなうプログラムの例を以下に示す。

```{code-cell}
import numpy as np
import matplotlib.pyplot as plt

X = [10, 20, 30, 50]
Y = [100, 300, 200, 300]

def lagrange(X, Y, x):
    n = len(X)
    s = []
    for k in range(n):
        t = Y[k] * np.prod([(x - X[i]) / (X[k] - X[i]) for i in range(n) if i != k])
        s.append(t)
    return np.sum(s)

xs = np.linspace(min(X), max(X), 101)
ys = [lagrange(X, Y, x) for x in xs]
plt.plot(xs, ys, linestyle=':', color='blue')
plt.scatter(X, Y, color='red')
plt.tick_params(labelsize=18)
plt.gca().set_ylim(0, 350)
plt.grid(True)
plt.show()
```

以下は，2022年6月1日の京田辺市の[気温のデータ](https://www.data.jma.go.jp/obd/stats/etrn/view/hourly_a1.php?prec_no=61&block_no=0598&year=2022&month=6&day=1&view=)を用いた例である。
実際には，1時間ごと（あるいは 10 分ごと）の観測値が利用できるが，あえてデータを抜き取ってある。

```{code-cell}
import numpy as np
import matplotlib.pyplot as plt

X_all = np.arange(1, 25)
Y_all = np.array([
    13.7, 12.4, 12.4, 12.2, 12.2, 12.7,
    14.3, 17.5, 19.2, 22.0, 24.2, 24.4,
    24.9, 25.5, 26.9, 26.1, 25.7, 24.8,
    23.4, 19.6, 18.0, 17.6, 16.5, 15.9,
])

plt.scatter(X_all, Y_all, color='green', label='hidden')

K = [1, 6, 11, 16, 21]
K = np.array([k in K for k in range(len(X_all))])

def lagrange(X, Y, x):
    n = len(X)
    s = []
    for k in range(n):
        t = Y[k] * np.prod([(x - X[i]) / (X[k] - X[i]) for i in range(n) if i != k])
        s.append(t)
    return np.sum(s)

X = X_all[K]
Y = Y_all[K]
xs = np.linspace(1, 24, 1001)
ys = [lagrange(X, Y, x) for x in xs]
plt.plot(xs, ys, linestyle=':', color='blue')
plt.scatter(X, Y, color='red', label='observed')

plt.xlabel('時刻（時）', fontsize=18)
plt.ylabel('気温（°C）', fontsize=18)
plt.tick_params(labelsize=18)
plt.legend(fontsize=18)
plt.gca().set_xlim(0, 25)
plt.gca().set_ylim(10, 30)
plt.grid(True)
plt.show()
```

24個の点のうち5点だけ用いても，それなりに高い精度で推測できている様子が見られる。

## ラグランジュ補間の問題点

上述の気温データについて，単純に考えると，観測値を増やせば補間の精度が良くなるように思える。

実際に12個の点を通るラグランジュ補間を求める例を以下に示す。

```{code-cell}
import numpy as np
import matplotlib.pyplot as plt

X_all = np.arange(1, 25)
Y_all = np.array([
    13.7, 12.4, 12.4, 12.2, 12.2, 12.7,
    14.3, 17.5, 19.2, 22.0, 24.2, 24.4,
    24.9, 25.5, 26.9, 26.1, 25.7, 24.8,
    23.4, 19.6, 18.0, 17.6, 16.5, 15.9,
])

plt.scatter(X_all, Y_all, color='green', label='hidden')

K = np.array([k % 2 == 1 for k in range(len(X_all))])

def lagrange(X, Y, x):
    n = len(X)
    s = []
    for k in range(n):
        t = Y[k] * np.prod([(x - X[i]) / (X[k] - X[i]) for i in range(n) if i != k])
        s.append(t)
    return np.sum(s)

X = X_all[K]
Y = Y_all[K]
xs = np.linspace(1, 24, 1001)
ys = [lagrange(X, Y, x) for x in xs]
plt.plot(xs, ys, linestyle=':', color='blue')
plt.scatter(X, Y, color='red', label='observed')

plt.xlabel('時刻（時）', fontsize=18)
plt.ylabel('気温（°C）', fontsize=18)
plt.tick_params(labelsize=18)
plt.legend(fontsize=18)
plt.gca().set_xlim(0, 25)
plt.gca().set_ylim(10, 30)
plt.grid(True)
plt.show()
```

たしかに多くの観測点を通ってはいるものの，点と点の間については局所的に誤差が大きい箇所が見られる。
特に端の方では極端な（実際には考えられない）ピークを示している。
このような現象は **オーバーシュート**（overshoot）や **ルンゲの現象**（Runge's phenomenon）とよばれる。

オーバーシュートを防ぐためには，関数が複雑になり過ぎないように，適度に小さな次数を設定する必要がある。
つまり，曲線が通る点の個数が多すぎてはいけない。

とはいえ，サンプル点が多い方が精度が下がるというのは直感に合わない。
精度を保つための方法として，サンプル点を減らす以外には，たとえば，区分的にラグランジュ補間を適用するという方法が考えられる。

以下は，全体を 4 つの区間に分けて区分的にラグランジュ補間を適用した例である。

```{code-cell}
:tags: [remove-input, remove-stderr]

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange

n = 24
X_all = np.arange(1, n + 1)
Y_all = np.array([
    13.7, 12.4, 12.4, 12.2, 12.2, 12.7,
    14.3, 17.5, 19.2, 22.0, 24.2, 24.4,
    24.9, 25.5, 26.9, 26.1, 25.7, 24.8,
    23.4, 19.6, 18.0, 17.6, 16.5, 15.9,
])

plt.scatter(X_all, Y_all, color='green', label='hidden')

K_all = np.array([k for k in range(n) if k % 2 == 1])
Ks = [K_all[i:i+4] for i in range(0, 12, 3)]
for h in range(0, 12, 3):
    x = X_all[K_all[h]]
    plt.axvline(x=x, linestyle=':', color='black')

for K in Ks:
    K = [k in K for k in range(n)]
    X = X_all[K]
    Y = Y_all[K]
    f = lagrange(X, Y)
    xs = np.linspace(min(X), max(X), 1001)
    ys = [f(x) for x in xs]
    plt.plot(xs, ys, linestyle=':', color='blue')
plt.scatter(X_all[K_all], Y_all[K_all], color='red', label='observed')

plt.xlabel('時刻（時）', fontsize=18)
plt.ylabel('気温（°C）', fontsize=18)
plt.tick_params(labelsize=18)
plt.legend(fontsize=18)
plt.gca().set_xlim(0, 25)
plt.gca().set_ylim(10, 30)
plt.grid(True)
plt.show()
```

極端な誤差はなくなり安定しているように見える。
しかし，区間ごとの境界において関数が微分不可能となっている。
関数が連続なだけでなく，1次導関数や2次導関数も連続であれば，より滑らかな関数になると考えられる。

## スプライン補間

低次の多項式によって区分的に補間をおこないながら，かつ，区間の境界を滑らかに接続する方法のひとつに，**スプライン補間**（spline interpolation）がある。

Wikipedia によれば，[スプライン（spline）](https://ja.wikipedia.org/wiki/%E3%82%B9%E3%83%97%E3%83%A9%E3%82%A4%E3%83%B3%E6%9B%B2%E7%B7%9A)とは，「製図などに用いられる一種の自在定規で、しなやかで弾力のある細長い板」のことである。

```{figure} ./images/wikipedia_spline.png
:alt: spline
:align: center
スプライン
```

この物理的な「スプライン」は，数学的には，3 次スプライン曲線を得るための道具である。
スプライン曲線（spline curve）は 3 次スプライン曲線以外にもさまざまなものが知られている。

まず，1 次スプライン曲線とは区分線形補間（piecewise linear interpolation）のことである。
つまり，与えられた点に対して，隣り合う点の間を 1 次関数によって補間する。
得られるのは，折れ線グラフである。

一般に，$n$ 次スプライン曲線では $n$ 次関数を用いる。
よく用いられるものは **自然３次スプライン**（natural cubic spline） である[^nurbs]。
以下では自然３次スプラインについて説明する。

[^nurbs]: コンピュータグラフィックスで用いられるものは NURBS（Non-Uniform Rational B-Spline; 非一様有理Bスプライン）とよばれるもので，制御点を必ずしも通らないという性質を持つ，拡張版スプライン補間である。

```{prf:example} 3 点に対するスプライン補間
:label: example-spline-1

曲線 $y=S(x)$ が 3 点 $(-1,16),(0,4),(1,12)$ を通り，かつ，以下の条件を満たすように関数 $S(x)$ を求める。

- 区間 $[-1,0]$ において，$3$ 次関数 $S_1(x)$ が定義されて，$S(x)=S_1(x)$ が成り立つ。
- 区間 $[0,1]$ において，$3$ 次関数 $S_2(x)$ が定義されて，$S(x)=S_2(x)$ が成り立つ。
- $S(-1)=S_1(-1)=16$，$S(0)=S_1(0)=S_2(0)=4$，$S(1)=S_2(1)=12$ である。
- $S(x)$ の 1 次導関数 $S'(x)$ と 2 次導関数 $S''(x)$ は $x=0$ において連続である。
- $S''(-1)=S''(1)=0$ である。

3 次関数 $S_1(x)$，$S_2(x)$ を以下のようにおく。

$$
S_1(x) = a x^3 + b x^2 + c x + d
\\
S_2(x) = e x^3 + f x^2 + g x + h
$$

8 個のパラメータ $a,b,c,d,e,f,g,h$ を定めると関数 $S_1(x)$，$S_2(x)$ が定まり，さらに関数 $S(x)$ も定まる。
満たすべき条件も 8 個でありパラメータの個数と同数であるから，典型的な場合においては，パラメータが一意に定まる。

まず，微分を求める。

$$
S'_1(x) &= 3 a x^2 + 2 b x + c
\\
S'_2(x) &= 3 e x^2 + 2 f x + g
\\
S''_1(x) &= 6 a x + 2 b
\\
S''_2(x) &= 6 e x + 2 f
$$

条件より，

$$
& S_1(-1) = - a + b - c + d = 16
\\
& S_1(0) = d = 4
\\
& S_2(0) = h = 4
\\
& S_2(1) = e + f + g + h = 12
\\
& S'_1(0) = S'_2(0) \quad \text{つまり} \quad c = g
\\
& S''_1(0) = S''_2(0) \quad \text{つまり} \quad b = f
\\
& S''_1(-1) = -6a+2b = 0
\\
& S''_2(1) = 6e+2f = 0
$$

これを解いて以下を得る。

$$
(a,b,c,d)&=(5,15,-2,4)
\\
(e,f,g,h)&=(-5,15,-2,4)
$$

すなわち，

$$
S_1(x) &= 5 x^3 + 15 x^2 - 2 x + 4
\\
S_2(x) &= -5 x^3 + 15 x^2 - 2 x + 4
$$
```

```{code-cell}
:tags: [remove-input]

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange

s1 = lambda x: 5 * x ** 3 + 15 * x ** 2 - 2 * x + 4
s2 = lambda x: -5 * x ** 3 + 15 * x ** 2 - 2 * x + 4

X = np.linspace(-1,1,101)
plt.plot(X, s1(X), linestyle=':', color='blue')
plt.plot(X, s2(X), linestyle=':', color='green')
X = np.linspace(-1,0,101)
plt.plot(X, s1(X), color='blue', label='S1')
X = np.linspace(0,1,101)
plt.plot(X, s2(X), color='green', label='S2')

X = [-1,0,1]
Y = [16,4,12]
plt.scatter(X, Y, color='red')

#X = [-1,0,1]
#Y = [16,4,12]
#f = lagrange(X, Y)
#X = np.linspace(-1,1,101)
#plt.plot(X, f(X), color='orange', label='Lagrange')

plt.tick_params(labelsize=18)
plt.grid(True)
plt.legend(fontsize=18)
plt.show()
```

```{prf:example} 3 点に対するスプライン補間
:label: example-spline-3-1

点 $(x_1,y_1),(x_2,y_2),(x_3,y_3)$ が与えられたとき，以下の条件を満たす関数 $S(x)$ を求める。

- 区間 $[x_1,x_2]$ において，$3$ 次関数 $S_1(x)$ が定義されて，$S(x)=S_1(x)$ が成り立つ。
- 区間 $[x_2,x_3]$ において，$3$ 次関数 $S_2(x)$ が定義されて，$S(x)=S_2(x)$ が成り立つ。
- $k=1,2,3$ について $S(x_k)=y_k$ が成り立つ。すなわち，$S_1(x_1)=y_1$，$S_1(x_2)=S_2(x_2)=y_2$，$S_2(x_3)=y_3$ が成り立つ。
- $S(x)$ の 1 次導関数 $S'(x)$ と 2 次導関数 $S''(x)$ は $x=x_2$ において連続である。
- $S''(x_1)=0$ かつ $S''(x_3)=0$ である。

$k=1,2$ について，関数 $S_k(x)$ は 3 次関数であるから，以下のように表すことができる。

$$
S_k(x) = a_k x^3 + b_k x^2 + c_k x + d_k
$$

つまり，この問題には求めるべきパラメータが 8 個存在する。
満たすべき条件も 8 個であるから，典型的な場合においてはパラメータが一意に定まる。

ここで， $S_1(x_2)=S_2(x_2)=y_2$ を満たすことから，新たなパラメータ $p_k,q_k$ を用いて以下のように書きかえることができる。

$$
S_k(x) = a_k (x - x_2)^3 + p_k (x - x_2)^2 + q_k (x - x_2) + y_2
$$

$S_k(x)$ を微分すると

$$
S'_k(x) &= 3 a_k (x - x_2)^2 + 2 p_k (x - x_2) + q_k
\\
S''_k(x) &= 6 a_k (x - x_2) + 2 p_k
$$

$S'(x)$ と $S''(x)$ が連続であるということは，つまり，$S'_1(x_2)=S'_2(x_2)$ および $S''_1(x_2)=S''_2(x_2)$ が成り立つということである。

ここから，ある定数 $p,q$ が存在して，

$$
p_1 = p_2 = p, \quad q_1 = q_2 = q
$$

が得られる。

つまり，あとは以下を解けばよい。

$$
& S_1(x) = a_1 (x - x_2)^3 + p (x - x_2)^2 + q (x - x_2) + y_2
\\
& S_2(x) = a_2 (x - x_2)^3 + p (x - x_2)^2 + q (x - x_2) + y_2
\\
& S_1(x_1) = y_1
\\
& S_2(x_3) = y_3
\\
& S''_1(x_1) = 0
\\
& S''_2(x_3) = 0
$$

これらを解くと，以下が得られる。

$$
a_1 &= \frac{p}{3 h_1}
\\
a_2 &= - \frac{p}{3 h_2}
\\
p &= \frac{3}{2} \frac{- g_1 + g_2}{h_1 + h_2}
\\
q &= \frac{g_1 h_2 + g_2 h_1}{h_1 + h_2}
$$

ただし，$h_1=x_2-x_1$，$h_2=x_3-x_2$，$g_1=\frac{y_2-y_1}{x_2-x_1}$，$g_2=\frac{y_3-y_2}{x_3-x_2}$ である。
```

一般に，座標平面上の $n+1$ 個の点 $(x_0,y_0),(x_1,y_1),\cdots,(x_n,y_n)$ が与えられたとき，自然 3 次スプラインは以下の条件を満たす関数 $S(x)$ として定義される。

- $k=0,1,\cdots,n-1$ について，区間 $[x_k,x_{k+1}]$ において，$3$ 次関数 $S_k(x)$ が定義されて，$S(x)=S_k(x)$ が成り立つ。
- $k=0,1,\cdots,n-1$ について $S_k(x_k)=y_k$ かつ $S_k(x_{k+1})=y_{k+1}$ が成り立つ。
- $S(x)$ の 1 次導関数 $S'(x)$ と 2 次導関数 $S''(x)$ は $x=x_1,x_2,\cdots,x_{n-1}$ において連続である。
- $S''_0(x_0)=0$ かつ $S''_{n-1}(x_n)=0$ である。

求めるべきパラメータは $4n$ 個（$n$ 個の区間ごとに，3次関数の係数が $4$ 個ずつ），条件は $2n+2(n-1)+2=4n$ 個であるから，（典型的な場合においては）条件を満たすパラメータが一意に定まる。

これらの条件を満たす関数 $\{S_k(x)\}$ を求めよう。

各区間 $[x_k,x_{k+1}]$（$k=0,1,\cdots,n-1$）において，$S_k(x)$ が 3 次関数であり $S_k(x_k)=y_k$ を満たすことから，$S_k(x)$ を以下のように表すことができる。

$$
S_k(x) = a_k (x - x_k)^3 + b_k (x - x_k)^2 + c_k (x - x_k) + y_k
$$

ただし，$a_k,b_k,c_k$ は定数である。これを繰り返し微分することにより，以下を得る。

$$
S'_k(x) &= 3 a_k (x - x_k)^2 + 2 b_k (x - x_k) + c_k
\\
S''_k(x) &= 6 a_k (x - x_k) + 2 b_k
$$

$S''(x)$ が連続であることから，未知のパラメータ $\sigma_k$（$k=0,\cdots,n$）を $S''(x_k) = \sigma_k$ が成り立つように定めることができる。

このとき，$S''_k(x_k) = 2 b_k = \sigma_k$ であるから以下を得る。

$$
b_k = \frac{\sigma_k}{2}
$$

同様に，$S''_k(x_{k+1}) = 6 a_k (x_{k+1} - x_k) + 2 b_k = \sigma_{k+1}$ より，以下を得る。

$$
a_k = \frac{\sigma_{k+1} - \sigma_{k}}{6 (x_{k+1} - x_{k})}
$$

さらに，$S(x)$ が $x=x_{k+1}$ で連続，すなわち $S_k(x_{k+1}) = y_{k+1}$ が成り立つから

$$
c_k = \frac{y_{k+1} - y_{k}}{x_{k+1} - x_{k}} - \frac{\sigma_{k+1} + 2\sigma_{k}}{6} (x_{k+1} - x_k)
$$

最後に，$S'(x)$ が $x=x_1,x_2,\cdots,x_{n-1}$ において連続であるから $S'_k(x_{k+1}) = S'_{k+1}(x_{k+1})$ より

$$
3 a_k (x_{k+1} - x_{k})^2 + 2 b_k (x_{k+1} - x_{k}) + c_k = c_{k+1}
$$

パラメータ $a_k,b_k,c_k$ の式を代入すると以下の式を得る。

$$
h_{k} \sigma_{k} + 2(h_{k} + h_{k+1}) \sigma_{k+1} + h_{k+1} \sigma_{k+2} = 6 (g_{k+1} - g_{k})
$$

ただし，$h_k=x_{k+1}-x_{k}$，$g_k=\frac{y_{k+1}-y_{k}}{x_{k+1}-x_{k}}$ とする。

上の方程式を $k=0,1,\cdots,n-2$ について連立すると，以下の連立方程式を得る。

$$
\begin{bmatrix}
h_0 & 2(h_0 + h_1) & h_1 & 0 & 0 & \cdots & 0 & 0 & 0 \\
0 & h_1 & 2(h_1 + h_2) & h_2 & 0 & \cdots & 0 & 0 & 0 \\
0 & 0 & h_2 & 2(h_2 + h_3) & h_3 & \cdots & 0 & 0 & 0 \\
\vdots & \vdots & \vdots & \vdots & \vdots & \ddots & \ddots & \vdots & \vdots \\
0 & 0 & 0 & 0 & 0 & \cdots & h_{n-2} & 2(h_{n-2} + h_{n-1}) & h_{n-1}
\end{bmatrix}
\begin{bmatrix}
\sigma_0 \\ \sigma_1 \\ \sigma_2 \\ \vdots \\ \sigma_n
\end{bmatrix}
=
\begin{bmatrix}
6 (g_1 - g_0) \\
6 (g_2 - g_1) \\
6 (g_3 - g_2) \\
\vdots \\
6 (g_{n-1} - g_{n-2})
\end{bmatrix}
$$

この連立1次方程式を境界条件 $\sigma_{0}=\sigma_{n}=0$ のもとで解き，解 $\{\sigma_k\}$ から $\{a_k\}$，$\{b_k\}$，$\{c_k\}$ を求めることによって，関数 $S_k(x)$ を定めることができる。

## 例（スプライン補間）

上述と同じ気温のデータについて，スプライン補間を適用した例を以下に示す。

ここでは，上で導出した方程式を解く方法ではなく，SciPy ライブラリの関数 [`scipy.interpolate.CubicSpline()`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html) を用いた。オプション引数 `bc_type='natural'` を与えることによって自然スプライン補間を適用できる。

```{code-cell}
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

X_all = np.arange(1, 25)
Y_all = np.array([
    13.7, 12.4, 12.4, 12.2, 12.2, 12.7,
    14.3, 17.5, 19.2, 22.0, 24.2, 24.4,
    24.9, 25.5, 26.9, 26.1, 25.7, 24.8,
    23.4, 19.6, 18.0, 17.6, 16.5, 15.9,
])

plt.scatter(X_all, Y_all, color='green', label='hidden')

K = np.array([k % 2 == 1 for k in range(len(X_all))])

X = X_all[K]
Y = Y_all[K]
f = CubicSpline(X, Y, bc_type='natural')
xs = np.linspace(1, 24, 1001)
ys = f(xs)
plt.plot(xs, ys, linestyle=':', color='blue')
plt.scatter(X, Y, color='red', label='observed')

plt.xlabel('時刻（時）', fontsize=18)
plt.ylabel('気温（°C）', fontsize=18)
plt.tick_params(labelsize=18)
plt.legend(fontsize=18)
plt.gca().set_xlim(0, 25)
plt.gca().set_ylim(10, 30)
plt.grid(True)
plt.show()
```

グラフより，かなり高い精度で自然な補間を実現できていることが分かる。
