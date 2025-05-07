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

# 非線形方程式の数値解法 (2) 反復法

前回に引き続き，今回は非線型方程式の解法としての反復法を紹介する。

## 反復法

二分法は，解を含む区間を絞り込むことによって解を探索するという直接法のアルゴリズムであった。
これに対して，解 $\alpha$ に収束するような数列 $\{x_n\}$ の値を逐次的に求め，数列が収束したときの $x_n$ を近似解とする方法を **反復法**（iterative method）とよぶ。

反復法では，方程式 $f(x)=0$ を同値変形によって

$$
x = g(x)
$$

に変形し，適当な初期値 $x=x_0$ から逐次的に $x_1=g(x_0), x_2=g(x_1), \cdots$ を求める。
真の解 $x=\alpha$ では $f(\alpha)=0$ すなわち $\alpha=g(\alpha)$ が成り立つため，$x=\alpha$ は関数 $g$ の **不動点**（fixed point）である。

真の解 $\alpha$ は不明であるから，一般の場合に解の絶対誤差 $|x_n-\alpha|$ を求めることはできない。
したがって，数列の収束判定には以下のいずれかの条件を用いるのが一般的である。

1. $|x_{n+1}-x_{n}| < \delta$
2. $\displaystyle\frac{|x_{n+1}-x_{n}|}{|x_n|} < \delta$
3. $|f(x_{n})| < \varepsilon$

数列 $\{x_n\}$ がつねに収束するとは限らない。
次の定理は収束の十分条件を与える [^complete-space]。

[^complete-space]: ある空間が完備（complete）であることの説明は省略する。
実数全体の集合や実数上の閉区間は完備である。

```{prf:definition} リプシッツ連続

完備な区間 $I$ 上で定義された関数 $g:I \to I$ に対してある非負定数 $L \geq 0$ が存在し任意の $x,y \in I$ に対して

$$
|g(x) - g(y)| \leq L |x - y|
$$

が成り立つとき，$g$ は **リプシッツ連続**（Lipschitz continuous）であるという。

このとき，定数 $L$ をリプシッツ定数（Lipschitz constant），不等式をリプシッツ条件とよぶ。
```

```{prf:definition} 縮小写像

完備な区間 $I$ 上で定義された関数 $g:I \to I$ がリプシッツ連続であり，かつ，リプシッツ定数 $L$ が $1$ より小さいとき，関数 $g$ を **縮小写像**（contraction mapping）とよぶ。
```

```{prf:theorem} バナッハの不動点定理（Banach fixed-point theorem）
:label: thm-fixedpoint

完備な区間 $I$ 上で定義された関数 $g:I \to I$ が縮小写像であるならば，$g$ は $I$ において唯一の不動点を持ち，それは漸化式 $x_{n+1}=g(x_{n})$ によって得られる。
```

定理の証明は省略するが，その概要だけを述べる。

もし $g$ が縮小写像であるならば，ある $L<1$ が存在して

$$
|x_{n+1} - x_{n}|
= |g(x_{n}) - g(x_{n-1})|
\leq L |x_{n} - x_{n-1}|
\leq \cdots
\leq L^{n} |x_{1} - x_{0}|
$$

が成り立つ。
$0\leq L<1$ より $|x_{n+1}-x_{n}|$ は $0$ に収束する。
したがって，点列 $x_0,x_1,x_2,\cdots$ は収束する。

次に，どのような場合に関数 $g$ が縮小写像になり得るかについて確認する。

平均値の定理より，開区間 $(a,b)$ 上で定義された関数 $g$ について

$$
g'(c) = \frac{g(b) - g(a)}{b - a}
$$

かつ $a < c < b$ を満たす $c$ が存在する。
つまり，そのような $c$ は

$$
|g(b) - g(a)| = |g'(c)| \cdot |b - a|
$$

を満たす。

区間 $(a,b)$ において関数 $g$ が $C^1$ 級 [^class-ck] であるならば $|g'(c)| \leq L$ を満たす定数 $L$ が存在する。たとえば，$L=\sup_{x \in (a,b)} |g'(x)|$ はその例である。

[^class-ck]: 開区間 $I$ 上で定義された関数 $f$ と非負整数 $k$ について，関数 $f$ が $k$ 回微分可能であり，導関数 $f',f'',\cdots,f^{(k)}$ が連続であるとき，$f$ は **$C^k$ 級**（class $C^k$）であるという。

したがって，

$$
|g(b) - g(a)| = |g'(c)| \cdot |b - a| \leq L \cdot |b - a|
$$

が成り立つ。つまり，関数 $g$ はリプシッツ連続である。

このことから以下の系が得られる。

```{prf:corollary}
:label: cor:banach

完備な区間 $I$ 上で定義された関数 $g:I \to I$ が $|g'(x)| < 1$ を満たすとき，$g$ は $I$ において唯一の不動点を持ち，それは漸化式 $x_{n+1}=g(x_{n})$ によって得られる。
```

なお，$|g'(x)| < 1$ は十分条件であり，必ずしもこの条件を満たす必要はない。

## ニュートン法

非線形方程式を解くための反復法として最も基本的で有名なものが **ニュートン法**（Newton's method）である。
ニュートン法はニュートン・ラフソン法（Newton-Raphson method）ともよばれる[^newton-raphson]。

[^newton-raphson]: ニュートン（Isaac Newton, 1642-1727）が発見したニュートン法は，現在の形とはかなり異なる。ニュートンの方法をラフソン（Joseph Raphson, 1668-1715）が改良することによって，一般の方程式に適用可能な，現在知られる方法になった。歴史的な経緯についてはたとえば {cite}`osada_2010` や {cite}`osada_2018` 等を参照。

開区間 $(a,b)$ に解を持つ方程式 $f(x)=0$ を考える。
ただし，関数 $f(x)$ は開区間 $(a,b)$ で $C^2$ 級であり，
かつ，$f'(x) \neq 0$ とする。

このとき，テイラーの定理（Taylor's theorem）より，
任意の $x_0,x \in (a,b)$ に対して以下を満たすような実数 $\xi$ が存在する。

$$
f(x) = f(x_0) + f'(x_0) (x-x_0) + \frac{f''(\xi)}{2} (x-x_0)^{2}
$$

ちなみに $x_0 \le \xi \le x$ または $x \le \xi \le x_0$ である。

$|x-x_0|$ が十分に小さければ

$$
f(x) \approx f(x_0) + f'(x_0) (x - x_0)
$$

という近似式が得られる。
右辺をグラフで表すと，点 $(x_0,f(x_0))$ における関数 $f(x)$ の接線となる。
つまり，上の近似式は，関数 $f(x)$ を 1 次関数で近似することに相当する。
直感的には，多くの関数は（例外はあるが）ある微小区間だけを考えれば線形関数（グラフが直線や平面などになる）とみなすことができる。
これが，ニュートン法の考え方の基本である。

この近似式を用いて $f(x)=0$ を解くと以下の近似解を得る（仮定より $f'(x_0) \neq 0$ である点に注意）。

$$
x \approx x_0 - \frac{f(x_0)}{f'(x_0)}
$$

ここで，関数 $g$ を

$$
g(x) = x - \frac{f(x)}{f'(x)}
$$

と定義する。
この方程式を変形することによって，$f(x)=0$ と $x=g(x)$ は同値であることが確認できる。
もし関数 $g$ が不動点を持つならば，その不動点は方程式 $f(x)=0$ の解である。

漸化式 $x_{n+1}=g(x_n)$ を用いて数列 $x_0,x_1,x_2,\ldots$ を定義する。
この数列が収束するかどうかは関数 $f(x)$ や初期値 $x_0$ に依存する。
初期値が解 $\alpha$ に十分に近ければ数列 $\{x_n\}$ は解 $\alpha$ に収束することが知られている。

ニュートン・ラフソン法の擬似コードを {prf:ref}`alg:newton-raphson` に示す。

```{prf:algorithm} ニュートン・ラフソン法
:label: alg:newton-raphson

**Inputs** 関数 $f:\mathbb{R}\to\mathbb{R}$，初期値 $x_0$，許容誤差 $\varepsilon>0$，最大反復数 $N$

**Outputs** 方程式 $f(x)=0$ の 解 $x\in\mathbb{R}$

1. $x \leftarrow x_0$
2. for $n = 1 , 2, \ldots, N$ do
   1. $\displaystyle d \leftarrow -\frac{f(x)}{f'(x)}$
   2. $x \leftarrow x + d$
   3. if $|d| < \varepsilon$ then
      1. return $x$
3. "失敗" と表示する
```

## 例2（ニュートン法）

以下は方程式 $x^6+5x-4=0$ の解をニュートン法を用いて求めるプログラムの例である。
初期値を $x_0=1$，許容誤差を $10^{-6}$ としている。

```{code-cell}
"""
ニュートン法の実行例
"""

# ニュートン法によって方程式 f(x)=0 の解を求める。
# 引数 df は f の導関数である。
# 初期値を x=x0 とし，誤差が eps 未満になるか，もしくは
# 繰り返し回数が上限 N 回に達したら終了する。
def newton(f, df, x0, eps, N):
    x = x0
    for n in range(N):
        print(f'{n:3d} x={x:.6f} f(x)={f(x):.6f}')
        d = - f(x) / df(x)
        x += d
        if abs(d) < eps:
            return x
    raise Exception('失敗')

# 定数や関数を定義する
f = lambda x: x ** 6 + 5 * x - 4
df = lambda x: 6 * x ** 5 + 5
x0 = 1.0
eps = 1e-6
N = 30

# ニュートン法を実行する
x = newton(f, df, x0, eps, N)

# 解を表示する
print('x =', x)
```

```{code-cell}
:tags: [remove-input]

import numpy as np
import matplotlib.pyplot as plt

def plot_line_segment(x1, y1, x2, y2):
    X = np.linspace(x1, x2, 100)
    Y = np.linspace(y1, y2, 100)
    plt.plot(X, Y, linestyle=':', color='blue')

def get_span(X, eps=1e-10):
    xmin, xmax = np.min(X), np.max(X)
    w = np.max([eps, xmax - xmin])
    xmin, xmax = xmin - w * 0.1, xmax + w * 0.1
    return xmin, xmax

def plot_newton(f, df, x0, eps, N):
    X = []
    Y = []
    x = x0
    failure = True
    for n in range(N):
        P = (x, f(x))
        #print(f'{n:3d} x={x:.6f} f(x)={f(x):.6f}')
        X.append(P[0])
        Y.append(P[1])
        plot_line_segment(x, 0, *P)
        d = - f(x) / df(x)
        x += d
        plot_line_segment(*P, x, 0)
        plt.scatter(*P,  color='blue')
        if abs(d) < eps:
            failure = False
            break
    if failure:
        print('失敗')
    xmin, xmax = get_span(X)
    ymin, ymax = get_span(Y)
    XX = np.linspace(xmin, xmax, 100)
    # x 軸
    plt.axhline(0, linestyle=':', color='red')
    # y = f(x)
    plt.plot(XX, f(XX), color='black')
    #
    for k in range(len(X) - 1):
        plt.plot([X[k], X[k], X[k+1]], [0, Y[k], 0], linestyle=':', color='blue')
        plt.scatter(X[k], Y[k], color='blue')
    #
    plt.ylim(ymin, ymax)
    return x

f = lambda x: x ** 6 + 5 * x - 4
df = lambda x: 6 * x ** 5 + 5
x0 = 1.0
eps = 1e-6
N = 30
X = plot_newton(f, df, x0, eps, N)

plt.show()
```

4回の反復によって近似解 $x=0.7611184552121928$ を得ることができた。

許容誤差は $\varepsilon=10^{-6}$ としたが，これは必ずしも誤差限界が $\varepsilon$ であることを意味しない。
ニュートン法の誤差の解析については後述する（c.f. {numref}`sec:convergence-condition-newton`）。

## 例3（ニュートン法）

次に，方程式 $f(x) = x - \cos(x)=0$ をニュートン法で解いてみよう。
この方程式に対する代数的な解は（おそらく）存在しない。
しかし，グラフを描いてみると $0 < x < 1$ の範囲に解が存在することが分かる。
たしかに $f(0)=-1<0$ かつ $f(1) \approx 0.4597 > 0$ であり，中間値の定理より解の存在が保証される。

```{code-cell}
:tags: [remove-input]

import numpy as np
import matplotlib.pyplot as plt

f = lambda x: x - np.cos(x)
X = np.linspace(-5, 5, 101)
Y = f(X)
plt.plot(X, Y, color='black')
plt.axhline(y=0, linestyle=':', color='red')
plt.show()
```

初期値を $x=1$ とすると，$x$ の値が単調に減少して近似解に収束する。

```{code-cell}
# 定数や関数を定義する
f = lambda x: x - np.cos(x)
df = lambda x: 1 + np.sin(x)
eps = 1e-6
N = 30

# x の初期値を設定する
x0 = 1.0
print('初期値:', x0)

# ニュートン法を実行する
x = newton(f, df, x0, eps, N)

# 解を表示する
print('x =', x)
```

```{code-cell}
:tags: [remove-input]

plot_newton(f, df, x0, eps, N)
plt.show()
```

初期値を解から少し離れた $x=3$ としてやり直すと，繰り返しのたびに $x$ の値が振動するが，数ステップで近似解に収束する。

```{code-cell}
# x の初期値を設定する
x0 = 3.0
print('初期値:', x0)

# ニュートン法を実行する
x = newton(f, df, x0, eps, N)

# 解を表示する
print('x =', x)
```

```{code-cell}
:tags: [remove-input]

plot_newton(f, df, x0, eps, N)
plt.show()
```

初期値を解からさらに遠い $x=4$ としてやり直すと，$x$ の値は大きく振動し，ニュートン法は収束しない。

```{code-cell}
:tags: [raises-exception]

# x の初期値を設定する
x0 = 4.0
print('初期値:', x0)

# ニュートン法を実行する
x = newton(f, df, x0, eps, N)

# 解を表示する
print('x =', x)
```

```{code-cell}
:tags: [remove-input]

plot_newton(f, df, x0, eps, N)
plt.show()
```

このように，ニュートン法を実行するときは解に十分近い初期値を選ぶ必要がある。
ニュートン法が収束しない場合は，二分法を用いて初期値を解に十分に近づけてからやり直すとよい。

(sec:convergence-condition-newton)=
## 収束条件，収束の速さ

近似解の誤差について詳しく確認してみよう。

方程式 $f(x)=0$ の真の解を $x=\alpha$ とし，ニュートン法の $n$ 回の繰り返しによって求めた近似解を $x_n$ とする。

ニュートン法の反復式より

$$
x_{n+1} = g(x_n) = x_{n} - \frac{f(x_n)}{f'(x_n)}
\label{newton-1}
$$ (newton-1)

テイラーの定理により以下を満たす $\xi_n$ が存在し，かつ，$x_n < \xi_n < \alpha$ または $\alpha < \xi_n < x_n$ が成り立つ。

$$
f(\alpha)
= f(x_{n})
  + f'(x_{n}) (\alpha - x_{n})
  + \frac{f''(\xi_n)}{2} (\alpha - x_{n})^2
$$ (newton-2)

式 {eq}`newton-1` および式 {eq}`newton-2` および $f(\alpha)=0$ より，

$$
  x_{n+1} - \alpha
  &= x_{n} - \frac{f(x_n)}{f'(x_n)} - \alpha
\\
  &= -\frac{f(x_n) + (\alpha - x_{n}) f'(x_{n})}{f'(x_n)}
\\
  &= \frac{(\alpha - x_{n})^2 f''(\xi_n)}{2 f'(x_n)}
$$

が成り立つ。
適当な区間 $I$ において $x \in I \Rightarrow g(x) \in I$ が成り立ち，かつ， $0 < A \leq |f'(x)|$ かつ $|f''(x)| < B$ を満たすような定数 $A,B$ について，

$$
|x_{n+1} - \alpha|
\leq \frac{B}{2A} |x_{n} - \alpha|^2
$$

が成り立つ。
$\displaystyle |x_{n} - \alpha| < \frac{2A}{B}$ のとき $|x_{n+1} - \alpha| < |x_{n} - \alpha|$ が成り立ち数列 $\{x_n\}$ は収束する。
つまり，初期値を解の十分近くにとれば（$A$ や $B$ は定数と考えられるため）ニュートン法は必ず収束する。

上式より，ニュートン法の（繰返し1回あたりの）収束次数は $2$ である。
ただし，繰返し1回につき関数評価を2回（$f(x_n)$ と $f'(x_n)$）おこなうため，関数評価1回あたりの収束次数は $\sqrt{2} \approx 1.414$ である。

## 連立非線形方程式に対するニュートン法

以下のような，$n$ 個の変数にかんする $m$ 個の非線形方程式を考える。

$$
& f_1(x_1,x_2,\cdots,x_n) = 0
\\
& f_2(x_1,x_2,\cdots,x_n) = 0
\\
& \quad \vdots
\\
& f_m(x_1,x_2,\cdots,x_n) = 0
$$

ベクトル変数 $\boldsymbol{x}=[x_1,x_2,\cdots,x_n]^\top \in \mathbb{R}^n$ および
関数 $\boldsymbol{f}: \mathbb{R}^n \to \mathbb{R}^m$ を用いると，以下のように書き直すことができる。

$$
\boldsymbol{f}(\boldsymbol{x}) = \boldsymbol{0}
$$

多変数ベクトル値関数 $\boldsymbol{f}$ の **ヤコビ行列**（Jacobian matrix）とは，以下のように定義される行列値関数である。

$$
\boldsymbol{J}(\boldsymbol{x})
= \left[ \frac{\partial f_i}{\partial x_j}(\boldsymbol{x}) \right]_{ij}
= \begin{bmatrix}
    \displaystyle\frac{\partial f_1}{\partial x_1}(\boldsymbol{x}) & \cdots & \displaystyle\frac{\partial f_1}{\partial x_n}(\boldsymbol{x}) \\
    \vdots & \ddots & \vdots \\
    \displaystyle\frac{\partial f_m}{\partial x_1}(\boldsymbol{x}) & \cdots & \displaystyle\frac{\partial f_m}{\partial x_n}(\boldsymbol{x})
  \end{bmatrix}
$$

ただし，多変数関数 $f(x_1,x_2,\cdots,x_n)$ に対して $\displaystyle\frac{\partial f}{\partial x_j}$ は $f$ の $x_j$ による偏導関数であり，次のように定義される。

$$
\frac{\partial f}{\partial x_j}(x_1,x_2,\cdots,x_n)
= \lim_{h \to 0} \frac{f(x_1,x_2,\cdots,x_j+h,\cdots,x_n) - f(x_1,x_2,\cdots,x_j,\cdots,x_n)}{h}
$$

解 $\boldsymbol{x}=\boldsymbol{\alpha}$ に十分近い $\boldsymbol{x}$，$\boldsymbol{x}^{(0)}$ について，テイラーの定理より以下が成り立つ。

$$
\boldsymbol{f}(\boldsymbol{x})
\approx
\boldsymbol{f}(\boldsymbol{x}^{(0)})
+ \boldsymbol{J}(\boldsymbol{x}^{(0)}) (\boldsymbol{x} - \boldsymbol{x}^{(0)})
$$

よって，方程式 $\boldsymbol{f}(\boldsymbol{x}) = \boldsymbol{0}$ は次の線形方程式によって近似できる。

$$
  & \boldsymbol{A}_0 \cdot \boldsymbol{x}
    = \boldsymbol{b}_0
\\& \text{ただし}\quad
    \boldsymbol{A}_0 = \boldsymbol{J}(\boldsymbol{x}^{(0)}),
    \quad
    \boldsymbol{b}_0 = \boldsymbol{J}(\boldsymbol{x}^{(0)}) \cdot \boldsymbol{x}^{(0)}
                     - \boldsymbol{f}(\boldsymbol{x}^{(0)})
$$

求めた近似解によって，新たな近似線形方程式が得られる。
近似解を繰り返し求めて，十分に収束したところで停止すると，もとの方程式の近似解が得られる。

もし $m=n$（未知数の個数と方程式の本数が等しい）であり，かつ，ヤコビ行列 $\boldsymbol{J}(\boldsymbol{x})$ が逆行列を持つならば，

$$
\boldsymbol{g}(\boldsymbol{x})
=
\boldsymbol{x}
- \Bigl[ \boldsymbol{J}(\boldsymbol{x}) \Bigr]^{-1}
  \boldsymbol{f}(\boldsymbol{x})
$$

のように定義した関数 $\boldsymbol{g}$ について $\boldsymbol{x}=\boldsymbol{g}(\boldsymbol{x})$ と $\boldsymbol{f}(\boldsymbol{x})=\boldsymbol{0}$ は同値である。

```{note} 計算上の工夫
関数 $\boldsymbol{g}$ の計算を定義通りにおこなう場合，ヤコビ行列の逆行列が必要である。

$$
\boldsymbol{x}_{n+1}
=
\boldsymbol{g}(\boldsymbol{x}_{n})
=
\boldsymbol{x}_{n}
- \Bigl[ \boldsymbol{J}(\boldsymbol{x}_{n}) \Bigr]^{-1}
  \boldsymbol{f}(\boldsymbol{x}_{n})
$$

一般的に逆行列の計算は計算コストが高く，数値誤差の安定性においても不利であるため，なるべく避けることが望ましい。

実は，ベクトル $\boldsymbol{x}_{n+1}$ の値を求める際に逆行列を避けることが可能である。

以下では標記を単純にするため，行列 $\boldsymbol{J}(\boldsymbol{x}_{n})$ を $\boldsymbol{J}_{n}$ と書き，
ベクトル $\boldsymbol{f}(\boldsymbol{x}_{n})$ を $\boldsymbol{f}_{n}$ と書くことにする。
ベクトル $\boldsymbol{x}_{n}$ の値が与えられたとき，行列 $\boldsymbol{J}_{n}=\boldsymbol{J}(\boldsymbol{x}_{n})$ とベクトル $\boldsymbol{f}_{n}=\boldsymbol{f}(\boldsymbol{x}_{n})$ の値も計算することができる。

$$
\boldsymbol{x}_{n+1}
=
\boldsymbol{x}_{n}
- \Bigl( \boldsymbol{J}_{n} \Bigr)^{-1}
  \boldsymbol{f}_{n}
$$

次に，両辺に左から $\boldsymbol{J}_{n}$ をかける。

$$
\boldsymbol{J}_{n}
\boldsymbol{x}_{n+1}
=
\boldsymbol{J}_{n}
\boldsymbol{x}_{n}
- \boldsymbol{f}_{n}
$$

右辺のベクトルを $\boldsymbol{b}_{n}$ とおく。
ベクトル $\boldsymbol{x}_{n}$ の値が与えられたとき，ベクトル $\boldsymbol{b}_{n}$ の値も求めることができる。

すると，上の等式は次のように書き換えることができる。

$$
\boldsymbol{J}_{n}
\boldsymbol{x}_{n+1}
=
\boldsymbol{b}_{n}
$$

行列 $\boldsymbol{J}_{n}$ と ベクトル $\boldsymbol{b}_{n}$ の値は与えられており，ベクトル変数 $\boldsymbol{x}_{n+1}$ の値を未知である。

変数 $\boldsymbol{x}_{n+1}$ の値を求めるには，上の線形方程式を解けばよい。

次の例に示すプログラムでは，この工夫を用いて $\boldsymbol{x}$ の値を更新している。
```

## 例4（連立非線形方程式）

以下の非線形連立方程式を解いてみよう。

$$
  & f(x,y,z) = x^2 + y^2 + xz - x - y - 1 = 0
\\
  & g(x,y,z) = x^3 + z^3 + 3 x^2 - z^2 + 2yz - 2z - 4 = 0
\\
  & h(x,y,z) = 3xy + 4yz + 2xz - 2x - 3y - 4z = 0
$$

ヤコビ行列は次の通り。

$$
\boldsymbol{J}(x,y,z)
& = \begin{bmatrix}
      \displaystyle\frac{\partial f}{\partial x}(x,y,z)
    & \displaystyle\frac{\partial f}{\partial y}(x,y,z)
    & \displaystyle\frac{\partial f}{\partial z}(x,y,z)
    \\
      \displaystyle\frac{\partial g}{\partial x}(x,y,z)
    & \displaystyle\frac{\partial g}{\partial y}(x,y,z)
    & \displaystyle\frac{\partial g}{\partial z}(x,y,z)
    \\
      \displaystyle\frac{\partial h}{\partial x}(x,y,z)
    & \displaystyle\frac{\partial h}{\partial y}(x,y,z)
    & \displaystyle\frac{\partial h}{\partial z}(x,y,z)
    \end{bmatrix}
\\
& = \begin{bmatrix}
      2 x + z - 1
    & 2 y - 1
    & x
    \\
      3 x^2 + 6 x
    & 2 z
    & 3 z^2 + 2 y - 2 z - 2
    \\
      3 y + 2 z - 2
    & 3 x + 4 z - 3
    & 2 x + 4 y - 4
    \end{bmatrix}
$$

ヤコビ行列は，Python の記号計算・数式処理ライブラリである SymPy を用いて以下のように自動的に求めることができる。

```{code-cell}
"""
SymPy を用いてヤコビ行列を求める
"""

from IPython.display import display
from sympy import Matrix
from sympy.abc import x, y, z

# 関数 f を定義する
f = Matrix([
    x * x + y * y + x * z - x - y - 1,
    x * x * x + z * z * z + 3 * x * x - z * z + 2 * y * z - 2 * z - 4,
    3 * x * y + 4 * y * z + 2 * x * z - 2 * x - 3 * y - 4 * z,
])

# 変数 X を定義する
X = Matrix([x, y, z])

# 関数 f を表示する
display(f)

# f のヤコビ行列を求める
J = f.jacobian(X)

# ヤコビ行列を表示する
display(J)
```

連立非線形方程式に対するニュートン法の実装例を次に示す。
線形方程式を解く箇所は `numpy.linalg.solve()` を用いた。

```{code-cell}
"""
連立非線形方程式をニュートン法で解く
"""

from inspect import signature
import numpy as np
from numpy import linalg as LA
from sympy import Matrix, symbols
from sympy.utilities.lambdify import lambdify

# 乱数の初期化
rng = np.random.default_rng(1234567)

# 表示オプション
np.set_printoptions(precision=5, sign=' ', floatmode='fixed')

# ヤコビ行列を求める
def make_jacobian(f):
    arity = len(signature(f).parameters)
    X = symbols(f'x0:{arity}')

    # SymPy 関数
    f_sym = Matrix(f(*X))
    J_sym = f_sym.jacobian(X)

    # NumPy ユニバーサル関数
    f = lambdify(X, f_sym, 'numpy')
    J = lambdify(X, J_sym, 'numpy')

    return f, J

# ニュートン法（多変数版）
def newton(f, x0, eps, N):
    # f のヤコビ行列を求める
    f, J = make_jacobian(f)

    # 初期解 x = x0
    x = np.array(x0).reshape(-1)

    # 最大で N 回まで繰り返す
    for n in range(N):
        # 途中経過を表示する
        print(f'{n:3d} x={x} f(x)={f(*x).reshape(-1)}')

        # 行列 A を求める
        A = J(*x)

        # ベクトル b を求める
        b = A @ x - f(*x).reshape(-1)

        # 古い解を記録する
        x_old = x

        # 線形方程式 A x = b を解く
        x = LA.solve(A, b)

        # 解の変化が eps 未満であれば終了する
        if LA.norm(x_old - x, 1) < eps:
            return x
    raise Exception('失敗')

# 関数 f を定義する
f = lambda x, y, z: [
    x * x + y * y + x * z - x - y - 1,
    x * x * x + z * z * z + 3 * x * x - z * z + 2 * y * z - 2 * z - 4,
    3 * x * y + 4 * y * z + 2 * x * z - 2 * x - 3 * y - 4 * z,
]

# ニュートン法の初期値などを設定する
X0 = np.array(rng.random(3))
eps = 1e-6
N = 30

# ニュートン法を実行する
X = newton(f, X0, eps, N)

# 結果を表示する
print('X =', X)
```

解の一つが見つかる。この連立方程式は（少なくとも）5つの実数解を持つ。

より単純な別の例として，以下の非線形連立方程式を解いてみよう。

$$
  & f(x,y) = x + y - 3 = 0
\\
  & g(x,y) = x y - 2 = 0
$$

```{code-cell}
"""
連立非線形方程式をニュートン法で解く
"""

# 表示オプション
np.set_printoptions(precision=10, sign=' ', floatmode='fixed')

# 関数 f を定義する
f = lambda x, y: [
    x + y - 3,
    x * y - 2,
]

# 定数を定義する
eps = 1e-10
N = 30

# 初期値を設定する
X0 = np.array([10., -10.])

# ニュートン法を実行する
X = newton(f, X0, eps, N)

# 解を表示する
print('X =', X)
```

$f(2,1)=2+1-3=0$ かつ $g(x,y)=2\cdot 1-2=0$ より，$(x,y)=(2,1)$ は確かに解である。
この方程式は $x,y$ について対称であるから，値を入れ替えた $(x,y)=(2,1)$ と $(1,2)$ がともに解である。
数値解法によってどちらの解（の近似解）が見つかるかは初期値に依存する。

## SciPy を用いて非線形方程式を解く

SciPy でニュートン法を利用するには，[`scipy.optimize.root_scalar()`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root_scalar.html) を利用することができる。

```{code-cell}
"""
SciPy でニュートン法を用いる
"""

from scipy import optimize

# 関数とその導関数を定義する
f = lambda x: x ** 6 + 5 * x - 4
df = lambda x: 6 * x ** 5 + 5

# ニュートン法を実行する
sol = optimize.root_scalar(f, x0=1.0, fprime=df, method='newton', xtol=1e-6)

# 解を表示する
print(sol)
print('x =', sol.root)
```

二分法とニュートン法以外にも，多くのアルゴリズムが実装されている。
二分法は収束が保証されており確実ではあるが，いくつかの方法のなかで最も収束が遅い。
ニュートン法は解の近くを初期解とした場合の収束速度が非常に速い。

`scipy.optimize.root_scalar()` の引数 `method` を省略した場合，状況に応じて適切なアルゴリズムを推定して実行することになっている。
たとえば `bracket` が引数として与えられたときは `brentq`，初期解 `x0` と一階導関数 `fprime` が与えられたときは `newton`，二階導関数 `fprime2` も与えられたときは `halley` が選択される。
とはいえ，これらの選択が「最適」であることの保証があるわけではない点に注意が必要である。

`scipy.optimize.root_scalar()` はその名の通りスカラー関数（1変数関数）に適用するためのものである。
ベクトル関数（多変数関数）に対しては [`scipy.optimize.root()`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html) を用いることができる。
関数 `f` や初期値 `x0` の形式に注意が必要である。

```{code-cell}
"""
連立非線形方程式を scipy.optimize.root() で解く
"""

from scipy import optimize

# 表示オプション
np.set_printoptions(precision=20, sign=' ', floatmode='fixed')

# 関数を定義する
f = lambda x: [
    x[0] + x[1] - 3,
    x[0] * x[1] - 2,
]

# 初期値を設定する
x0 = np.array([10., -10.])

# 方程式を解く
sol = optimize.root(f, x0)

# 解を表示する
print(sol)
print('X =', sol.x)
```

[SciPy のマニュアル](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html)によれば，関数 `scipy.optimize.root()` がデフォルトで用いるアルゴリズムはパウエルのハイブリッド法と呼ばれるものである。
ヤコビ行列は f から数値的に計算されたものが用いられる。
計算を効率化するために，オプション引数としてヤコビ行列を与えることもできる。
