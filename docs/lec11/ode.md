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

(sec:ode)=
# 常微分方程式

$\newcommand\setR{\mathbb{R}}$

## 微分方程式とは

1次方程式や2次方程式などの，いわゆる方程式（equation）とは，

$$
x+10=15 \quad \text{や} \quad x^2=e^x
$$

のような等式が成立するような数値 $x$ を求める問題であった。
この場合，解は実数や整数などの値，すなわち数である。

これに対して，微分方程式（differential equation）とは

$$
\frac{dy}{dx}(x) = x^2 \quad \text{や} \quad \frac{dy}{dx}(x)+3y(x)=0
$$

のような等式が成立するような関数 $y(x)$ を求める問題である。
微分方程式の解は，数ではなく，関数である。

たとえば，ある直線状の道路を走る自動車について，時刻 $t$ における位置 $x(t)$，速度を $v(t)$ とおくと，

$$
\frac{dx}{dt}(t) = v(t)
$$

となり，微分方程式となる。
速度が一定の値 $V$ をとるとき，この方程式を解くと関数 $x(t) = V t + C$（$C$ は定数）が得られる。定数 $C$ は，たとえば，位置 $x(t)$ の時刻 $t=0$ における値 $x(0)$ を初期値として定めることによって得られる。

数値を解とする方程式は「点」のように静的な対象を求める問題であるのに対して，関数を解とする方程式は「関係」や「動き」のように動的な対象を求める問題とも言える。
日常生活における動的な現象を数学モデルで扱うとき，微分方程式によって表現されることが多い。

数学の講義や参考書で扱われる演習問題では解析的に解ける問題を扱うが，自然現象の数理モデルなど，多くの微分方程式は解析的に解くことが難しい[^ode]。

[^ode]: たとえばNewtonの運動方程式は単純な微分方程式であるが，物体（質点）が3つあるだけで解析的に解けなくなる（三体問題）。

そこで，計算機による数値計算を用いて微分方程式の解を求めることが行われる。

2つの変数 $x,y$ とその導関数によって表される微分方程式において，一方の変数 $y$ の値が他方の変数 $x$ の値によって定まるとき（すなわち，$y$ が $x$ の関数であるとき），$y$ を **従属変数**（dependent variable），$x$ を **独立変数**（independent variable）とよぶ。

3 個以上の変数を含む微分方程式についても同様に，従属変数 $y_1,y_2,\cdots$ の値が独立変数 $x_1,x_2,\cdots$ の値によって決まるものとする。

独立変数の個数が $1$ 個である微分方程式を **常微分方程式**（ordinary differential equation; ODE）とよぶ。
これに対して，$2$ 個以上の独立変数をもつ微分方程式を **偏微分方程式**（partial differential equation）とよぶ。
今回は，常微分方程式だけを扱う。

常微分方程式においては，導関数を

$$
\frac{dy}{dx}(x), \quad
\frac{d^2y}{dx^2}(x), \quad
\frac{d^3y}{dx^3}(x), \quad
\cdots, \quad
\frac{d^ny}{dx^n}(x), \quad
\cdots
$$

などと書く代わりに

$$
y', \quad y'', \quad y''', \quad \cdots, \quad y^{(n)}, \quad \cdots
$$

と書くことがある[^diff]。
常微分方程式では独立変数を1個であり混同するおそれがないことから，微分方程式を簡潔に表現するためによく用いられる。

[^diff]: 導関数を $\displaystyle \frac{dy}{dx}, \frac{d^2y}{dx^2}, \frac{d^3y}{dx^3}$ のように表すのはライプニッツの記法（Leibniz's notation）である。
これに対して，$y',y'',y'''$ のように表すのをラグランジュの記法（Lagrange's notation），$\dot{y},\ddot{y},\dddot{y}$ のように表すのをニュートンの記法（Newton's notation）とよぶ。

微分方程式に含まれる導関数の最高階数をその微分方程式の **階数**（order）とよぶ。
階数 $n$ の方程式を **$n$ 階微分方程式**（$n$-th order differential equation）とよぶ。
たとえば，常微分方程式 $y'+xy=0$ の階数は $1$ （$y'$ が最高階数の導関数）であり，常微分方程式 $y''=3y$ の階数は $2$ （$y''$ が最高階数の導関数）である。

$1$ 階の常微分方程式は $3$ つの変数 $x,y,y'$ についての方程式として表すことができる。特に，標準形では

$$
  y' = f(x,y)
$$

のように書くことができる。

````{prf:example}
:label: ode-standard-1

たとえば，

$$
y' = 3 x^2 y
$$

は標準形の $1$ 階常微分方程式である。

これは，以下のように変形できる。

$$
& \frac{1}{y} dy = 3 x^2 dx
\\
& \log y = x^3 + C \quad (\text{$C$は積分定数})
\\
& y = K \cdot \exp(x^3) \quad (\text{$K$は定数})
$$
````

標準形の $1$ 階常微分方程式は，一般には解の存在や一意性が保証されない。
ただし，一定の仮定のもとでは，解の存在や一意性が保証されるという事実が知られている[^ode-existence-of-solution]。

[^ode-existence-of-solution]: たとえば，ピカール＝リンデレーフの定理，ペアノの存在定理，カラテオドリの存在定理など。詳細については，たとえば，細矢（2019）{cite}`hosoya_2019` などを参照されたい。

微分方程式の解の一意性については，変数 $x,y$ の初期値 $x_0,y_0$ を定めることで解決できる場合がある。

```{prf:definition} 初期値問題
:label: ode-initial-value-problem

標準形の1階常微分方程式 $y'=f(x,y)$ について，独立変数 $x$ の値が $x_0$ であるときにおける従属変数 $y$ の値 $y(x_0) = y_0$ が与えられているとする。

このとき，初期条件のもとで微分方程式を解く問題，すなわち，

$$
  & y'(x) = f(x,y(x)) \\
  & y(x_0) = y_0
$$

を（1階常微分方程式の）**初期値問題**（initial value problem）とよぶ。
```

## オイラー法

区間は $[a,b]$ 上の $1$階常微分方程式の初期値問題を考える。

```{math}
:label: ivp-1
  & y' = f(x,y) \\
  & y(a) = y_0
```

区間 $[a,b]$ を刻み幅 $\displaystyle h=\frac{b-a}{n}$ で $n$ 等分し，$n+1$ 個の分点を $x_i=a+ih$（$i=0,1,\cdots,n$） とおく。

分点 $x_i$ に対応する従属変数 $y$ 真の値を $y_i=y(x_i)$ と書き，数値解を $Y_i$ と書くことにする。

初期値より，分点 $x=x_0$ に対する数値解 $Y_0$ は厳密な値 $y_0$ が分かっている。
次の分点 $x=x_1$ までの区間 $[x_0,x_1]$ における関数 $y(x)$ を直線で近似する。
点 $(x_0,Y_0)$ を通り，傾きが $y'(x_0)$ であるような直線の式は

$$
y = y_0 + (x - x_0) \cdot y'(x_0)
  = y_0 + (x - x_0) \cdot f(x_0,y_0)
$$

である。
これは，関数 $y(x)$ の $x=x_0$ における Taylor 展開に相当する。

この直線の式から，$x=x_1$ における数値解 $Y_1$ を求める。

$$
Y_1 = Y_0 + (x_1 - x_0) \cdot f(x_0,Y_0)
$$

次の区間 $[x_1,x_2]$ についても同様に，点 $(x_1,Y_1)$ を通り傾きが $f(x_1,Y_1)$ であるような直線によって関数 $y(x)$ を近似し，$x=x_2$ における値 $y(x_2)=y_2$ の近似値 $Y_2$ を求める。

$$
Y_2 = Y_1 + (x_2 - x_1) \cdot f(x_1,Y_1)
$$

以降も同様にして，非負整数 $k$ について，分点 $x=x_k$ における数値解 $Y_k$ を逐次的に求める。

$$
Y_{k+1} = Y_{k} + (x_{k+1} - x_{k}) \cdot f(x_{k},Y_{k})
$$

$\setR^2$ 上の関数 $\Phi(x,y)$ を

$$
\Phi(x,y) = y + h \cdot f(x,y)
$$

と定義すると，数値解 $\{Y_k\}$ は漸化式

$$
Y_0 = y_0, \quad Y_{k+1} = \Phi(x_k, Y_k)
$$

によって得られる。

## 例1（オイラー法）

区間 $[0,1]$ において，以下の 1 階常微分方程式の初期値問題を考える。

```{math}
:label: euler-ex-1

  & y' = 2xy \\
  & y(0) = 3
```

以下では，初期値問題を数値的に解く例としてこの問題を考える。

```{note}
ちなみに，この初期値問題は解析的に解くことができる。
解は

$$
y = 3 \cdot \exp(x^2) = 3 \cdot e^{x^2}
$$

である（これが解になることを確かめよ）。
一般の初期値問題について，このような解が得られるとは限らない。
以下では，解析的な解をいったん忘れて，あくまで数値的に解く状況を考える。
```

微分方程式を数値的に解いた結果として得られるのは，標本点 $x_0,x_1,\cdots,x_n$ における関数 $y(x)$ の値 $y_0,y_1,\cdots,y_0$ の数値解 $Y_0,Y_1,\cdots,Y_n$ である。
関数が得られるわけではないことに注意されたい。
とはいえ，標本点を十分に細かくとれば関数のグラフを十分に滑らかに描画することもできるし，補間や近似を用いれば多項式などの形で近似関数を求めることもできる。
多くの場合は，それで十分に実用的である。

ここでは，区間 $[0,1]$ を刻み幅 $\displaystyle h=\frac{1}{n}$ で等分し，$n+1$ 個の標本点 $\displaystyle 0,\frac{1}{n},\frac{2}{n},\cdots,1$ をとる。

Taylor の定理より，各 $i=0,1,\cdots,n-1$ について，

$$
  y(x_{i+1}) = y(x_{i}) + h \cdot y'(x_{i}) + \frac{h^2}{2} \cdot y''(x_{i}+\theta_i h)
$$

を満たす $0<\theta_i<1$ が存在する。$h \ll 1$ のとき，

$$
  y(x_{i+1}) \approx y(x_{i}) + h \cdot y'(x_{i})
$$

である。
式 {eq}`euler-ex-1` を代入すると

$$
  y(x_{i+1}) \approx y(x_{i}) + h \cdot 2 x_{i} y(x_{i})
$$

を得る。
関数 $\Phi(x,y)$ を $\Phi(x,y)=y+2hxy$ とおくと，

$$
  y_{i+1} \approx \Phi(x_{i},y_{i})
$$

を得る。
関数値を $y_i=y(x_i)$ ，数値解を $Y_i$ と書くことにすると，

$$
  Y_{i+1} = \Phi(x_{i},Y_{i})
$$

によって数値解を次々に求めることができる。

例えば，$n=5$ として $x = 0, 0.2, 0.4, \cdots, 1$ （つまり $x_0=0,x_1=0.2,x_2=0.4,\cdots,x_5=1$）における数値解 $Y_i$（$i=0,1,2,\cdots,5$）を求めると，{numref}`euler-values-1` のようになる。
参考として真の値 $y_i=y(x_i)=3\cdot\exp(x_i^2)$（$i=0,1,2,\cdots,5$）も併記する。

```{list-table} オイラー法の適用例（$y'=2xy$，$y(0)=3$）
:header-rows: 1
:name: euler-values-1

* - $x$
  - $0$
  - $0.2$
  - $0.4$
  - $0.6$
  - $0.8$
  - $1$
* - 数値解 $Y_i$
  - $3$
  - $3$
  - $3.24$
  - $3.7584$
  - $4.6604$
  - $6.1517$
* - 真の値 $y_i$
  - $3$
  - $3.1224$
  - $3.5205$
  - $4.3000$
  - $5.6894$
  - $8.1548$
```

これらの値を求める Python プログラムは，たとえば以下のようになる。

```{code-cell}
import numpy as np
import pandas as pd

def ivp_euler(f, a, b, y0, n):
    # X の値
    X = np.linspace(a, b, n + 1)

    # 標本点の間隔
    h = (b - a) / n

    # 数値解を求める
    Y = np.zeros_like(X)
    Y[0] = y0
    for k in range(n):
        Y[k + 1] = Y[k] + h * f(X[k], Y[k])

    return X, Y

# 分割数
n = 5

# 区間
a, b = 0, 1

# 初期値 y(a)
y0 = 3

# 関数 f(x,y)
f = lambda x, y: 2 * x * y

# オイラー法で初期値問題を解く
X, Y = ivp_euler(f, a, b, y0, n)

# 真の値を求める
y = 3 * np.exp(X ** 2)

# データフレームを作成し，表示させる
pd.options.display.precision = 4
pd.DataFrame({'数値解 Y': Y, '真の値 y': y}, index=X)
```

このような解法を **オイラー法**（Euler method）とよぶ。
{numref}`euler-values-1` からも分かる通り，オイラー法は誤差が大きく，実用にあたっては注意が必要である。

## ホイン法

以下では，オイラー法による近似値を $\bar{Y}_{i+1}$ と書くことにする。すなわち，

$$
  \bar{Y}_{i+1} = Y_i + h f(x_i,Y_i)
$$

オイラー法では区間 $(x_i,x_{i+1})$ における関数 $y(x)$ を直線で近似する。
多くの場合，関数は直線ではなく曲線であり，傾きは一定ではない。
$h$ が十分に小さいとき，一段法の更新式

$$
Y_{i+1} = Y_{i} + h \cdot \Phi(x_{i},x_{i+1},Y_{i};h)
$$

における $\Phi(x_i,x_{i+1},Y_i;h)$ の値を
$y'(x_i)$ と $y'(x_{i+1})$ の中間に設定するのがよい。
$y'(x_{i+1})=f(x_{i+1},y_{i+1})$ の真の値は不明であるため代わりに$f(x_{i+1},\bar{Y}_{i+1})$を用いたときの $y(x_{i+1})$ の予測値を $Y^*_{i+1}$ と書くことにする。
すなわち，

$$
  Y^*_{i+1} = Y_i + h f(x_{i+1},\bar{Y}_{i+1})
$$

この $\bar{Y}_{i+1}$ と $Y^*_{i+1}$ の算術平均を $Y_{i+1}$ とするのが **ホイン法**（Heun's method）である。
すなわち，

$$
  Y_{i+1} = \frac{1}{2} \left\{ \bar{Y}_{i+1} + Y^*_{i+1} \right\}
$$

## 例2（ホイン法）

{numref}`euler-values-1` と同じ初期値問題に対してホイン法を適用する。
結果を {numref}`heun-values-1` に示す。

```{list-table} ホイン法の適用例（$y'=2xy$，$y(0)=3$）
:header-rows: 1
:name: heun-values-1

* - $x$
  - $0$
  - $0.2$
  - $0.4$
  - $0.6$
  - $0.8$
  - $1$
* - 数値解 $Y_i$
  - $3$
  - $3.12$
  - $3.5144$
  - $4.2847$
  - $5.6490$
  - $8.0441$
* - 真の値 $y_i$
  - $3$
  - $3.1224$
  - $3.5205$
  - $4.3000$
  - $5.6894$
  - $8.1548$
```

これらの値を求める Python プログラムは，たとえば以下のようになる。

```{code-cell}
import numpy as np
import pandas as pd

def ivp_heun(f, a, b, y0, n):
    # X の値
    X = np.linspace(a, b, n + 1)

    # 標本点の間隔
    h = (b - a) / n

    # 数値解を求める
    Y = np.zeros_like(X)
    Y[0] = y0
    for k in range(n):
        k1 = f(X[k],     Y[k])
        k2 = f(X[k] + h, Y[k] + h * k1)
        Y[k + 1] = Y[k] + h * (k1 + k2) / 2

    return X, Y

# 分割数
n = 5

# 区間
a, b = 0, 1

# 初期値 y(a)
y0 = 3

# 関数 f(x,y)
f = lambda x, y: 2 * x * y

# ホイン法で初期値問題を解く
X, Y = ivp_heun(f, a, b, y0, n)

# 真の値を求める
y = 3 * np.exp(X ** 2)

# データフレームを作成し，表示させる
pd.options.display.precision = 4
pd.DataFrame({'数値解 Y': Y, '真の値 y': y}, index=X)
```

(sec:rungekutta)=
## ルンゲ・クッタ法

微分方程式の数値解法としてよく利用されるのが次数$4$の **ルンゲ・クッタ法**（Runge-Kutta method）である。
次数$4$のルンゲ・クッタ法は，古典的ルンゲ・クッタ法や RK4 ともよばれる。

まず，次数 $4$ のルンゲ・クッタ法のアルゴリズムを天下り式に示す。
その導出は，次回以降におこなう。

オイラー法と同じように，数値解 $Y_i$ を逐次的に求めてゆく。

$$
  & Y_0 = y_0 \\
  & Y_{i+1} = Y_{i} + h \cdot \Phi(x_i,Y_i) \\
  & \Phi(x_i,Y_i) = \frac{1}{6} (\alpha_i + 2 \beta_i + 2 \gamma_i + \delta_i)
$$

ただし，

$$
  & \alpha_i
    = f(x_i,Y_i)
\\& \beta_i
    = f\left( x_i+\frac{h}{2}, Y_i+\frac{h}{2} \alpha_i \right)
\\& \gamma_i
    = f\left( x_i + \frac{h}{2}, Y_i+\frac{h}{2} \beta_i \right)
\\& \delta_i
    = f(x_i + h, Y_i + h \gamma_i)
$$

とする。

## 例3（ルンゲ・クッタ法）

{numref}`euler-values-1` と同じ初期値問題に対してルンゲ・クッタ法を適用すると，{numref}`runge-kutta-values-1` のようになる。

```{list-table} ルンゲ・クッタ法の適用例（$y'=2xy$，$y(0)=3$）
:header-rows: 1
:name: runge-kutta-values-1

* - $x$
  - $0$
  - $0.2$
  - $0.4$
  - $0.6$
  - $0.8$
  - $1$
* - 数値解 $Y_i$
  - $3$
  - $3.1224$
  - $3.5205$
  - $4.3000$
  - $5.6893$
  - $8.1543$
* - 真の値 $y_i$
  - $3$
  - $3.1224$
  - $3.5205$
  - $4.3000$
  - $5.6894$
  - $8.1548$
```

```{code-cell}
import numpy as np
import pandas as pd

def ivp_runge_kutta_4(f, a, b, y0, n):
    # X の値
    X = np.linspace(a, b, n + 1)

    # 標本点の間隔
    h = (b - a) / n

    # 数値解を求める
    Y = np.zeros_like(X)
    Y[0] = y0
    for k in range(n):
        k1 = f(X[k],         Y[k])
        k2 = f(X[k] + h / 2, Y[k] + h / 2 * k1)
        k3 = f(X[k] + h / 2, Y[k] + h / 2 * k2)
        k4 = f(X[k] + h,     Y[k] + h * k3)
        Y[k + 1] = Y[k] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6

    return X, Y

# 分割数
n = 5

# 区間
a, b = 0, 1

# 初期値 y(a)
y0 = 3

# 関数 f(x,y)
f = lambda x, y: 2 * x * y

# ルンゲ・クッタ法で初期値問題を解く
X, Y = ivp_runge_kutta_4(f, a, b, y0, n)

# 真の値を求める
y = 3 * np.exp(X ** 2)

# データフレームを作成し，表示させる
pd.options.display.precision = 4
pd.DataFrame({'数値解 Y': Y, '真の値 y': y}, index=X)
```

オイラー法やホイン法と比較して，かなり精度が高いことが分かる。

SciPy で常微分方程式の初期値問題（IVP）を解く際には，[`scipy.integrate.solve_ivp()` 関数](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp)を利用できる。
ドキュメントによれば， `method='DOP853'` を指定することによって次数 8 の陽的ルンゲ・クッタ法を利用できる（これが推奨されている）。
デフォルトは `'RK45'` つまり次数 4 の陽的ルンゲ・クッタ法が用いられるようになっている。
関数によっては，より安定な陰的ルンゲ・クッタ法を使うべき場合もあるが，計算コストが高い。

```{code-cell}
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

# 分割数
n = 5

# 区間
a, b = 0, 1

# 初期値 y(a)
y0 = 3

# 関数 f(x,y)
f = lambda x, y: 2 * x * y

# ルンゲ・クッタ法で初期値問題を解く
X = np.linspace(a, b, n + 1)
Y = solve_ivp(f, (a, b), (y0,), method='DOP853', t_eval=X).y.reshape(-1)

# 真の値を求める
y = 3 * np.exp(X ** 2)

# データフレームを作成し，表示させる
pd.options.display.precision = 4
pd.DataFrame({'数値解 Y': Y, '真の値 y': y}, index=X)
```
