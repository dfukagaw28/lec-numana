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

# 多変数関数

## おさらい: 1変数の微分

````{prf:definition}
:label: def-1d-diff

1 変数関数 $f(x)$ が $x=a$ およびその近傍で定義されているとする。
このとき，$f(x)$ が $x=a$ において微分可能であるとは，以下の極限 $\alpha$ が存在することである。

$$
\lim_{h\to 0} \frac{f(a + h) - f(a)}{h} = \alpha
$$
````

別の表現で言えば，以下を満たすような定数 $\alpha$ が存在することである。

$$
f(x) = f(a) + \alpha (x-a) + g(x)
\quad\text{かつ}\quad
\lim_{x\to a} \frac{g(x)}{x-a} = 0
$$

上の条件を満たす定数 $\alpha$ を $x=a$ における $f(x)$ の微分係数とよび，$f'(a)$ と書く。

微分係数 $f'(a)$ は $a$ の関数である。$f'$ を $f$ の（1次）導関数とよぶ。

```{image} ./images/bibun1.gif
:alt: 1階微分
:width: 400px
:align: center
```

## 2変数の微分

````{prf:definition}
:label: def-2d-diff

2 変数関数 $f(x,y)$ が $(x,y)=(a,b)$ およびその近傍で定義されているとする。
このとき，$f(x,y)$ が $(x,y)=(a,b)$ において微分可能（全微分可能）であるとは，以下を満たすような定数 $\alpha$，$\beta$ が存在することである。

$$
\begin{align*}
& f(x,y) = f(a,b) + \alpha (x-a) + \beta (y-b) + g(x,y)
\\
& \text{かつ}\quad
\lim_{(x,y)\to (a,b)} \frac{g(x,y)}{|(x,y)-(a,b)|} = 0
\end{align*}
$$
````

ベクトル $(\alpha,\beta)\in\mathbb{R}^2$ を $(x,y)=(a,b)$ における $f(x,y)$ の全微分係数とよび，$f'(a,b)$ と書く。

定数 $\alpha$ を $(x,y)=(a,b)$ における $f(x,y)$ の変数 $x$ に関する偏微分係数とよび，$\displaystyle\frac{\partial f(a,b)}{\partial x}$ や $f_x(a,b)$ などと書く。
定数 $\beta$ を $(x,y)=(a,b)$ における $f(x,y)$ の変数 $y$ に関する偏微分係数とよび，$\displaystyle\frac{\partial f(a,b)}{\partial y}$ や $f_y(a,b)$ などと書く。

このとき，

$$
\alpha = \frac{\partial f(a,b)}{\partial x} = \lim_{x\to a}\frac{f(x,b)-f(a,b)}{x-a}
$$

および

$$
\beta = \frac{\partial f(a,b)}{\partial y} = \lim_{y\to b}\frac{f(a,y)-f(a,b)}{y-b}
$$

が成り立つ。

以上より，関数 $f$ が全微分可能かつ各変数に関して偏微分可能であるとき，全微分係数と偏微分係数との間には以下の関係が成り立つ。

$$
f'(a,b) = \biggl( \frac{\partial f(a,b)}{\partial x}, \frac{\partial f(a,b)}{\partial y} \biggr)
$$

ただし，各変数に関して偏微分可能であるからといって，全微分可能であるとは限らない。

```{code-cell}
:tags: [remove-input, remove-output]

import matplotlib.pyplot as plt
from matplotlib import rcParams
from myst_nb import glue
import numpy as np

plt.rcParams.update({
  'mathtext.fontset': 'cm',
  'font.size': 18,
})

def make_graph(f, x, y):

    # メッシュグリッドを作成する
    X, Y = np.meshgrid(x, y)
    Z = f(X, Y)

    # Figureを作成する
    fig = plt.figure(figsize=(14, 6))

    # 2D の Contour グラフを作成する
    ax = fig.add_subplot(1, 2, 1)
    contour = ax.contourf(X, Y, Z, 20, cmap='viridis')
    fig.colorbar(contour, ax=ax)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_aspect('equal')

    # 3D グラフを作成する
    subsample_step = 5
    X_sub, Y_sub, Z_sub = (A[::subsample_step, ::subsample_step] for A in (X, Y, Z))
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    ax.plot_wireframe(X_sub, Y_sub, Z_sub, color='blue')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_zlabel('$f(x, y)$')

    #plt.tight_layout()
    plt.subplots_adjust(wspace=0.01)
    plt.show()

    return fig

f = lambda x, y: 2 * x ** 2 - 3 * y
x = np.linspace(-2, 2, 100)
y = np.linspace(-2, 2, 100)
fig = make_graph(f, x, y)
glue('fun-2d-diff-1', fig, display=False)

f = lambda x, y: np.minimum(np.abs(x), np.abs(y))
x = np.linspace(-2, 2, 500)
y = np.linspace(-2, 2, 500)
fig = make_graph(f, x, y)
glue('fun-2d-diff-2', fig, display=False)
```

````{prf:example}
:label: example-2d-diff-1

関数 $f(x,y)=2x^2-3y$ とおくと，$f(x,y)$ は任意の点 $(x,y)$ において全微分可能である。

たとえば， $(x,y)=(1,1)$ における偏微分係数は

$$
\frac{\partial f(1,1)}{\partial x} = 4, \quad
\frac{\partial f(1,1)}{\partial y} = -3
$$

となる。 $\displaystyle g(x,y) = f(x,y) - f(1,1) - \frac{\partial f(1,1)}{\partial x} (x-1) - \frac{\partial f(1,1)}{\partial y} (y-1)$ とおくと

$$
g(x,y)
= (2x^2-3y) - (-1) - 4 (x-1) + 3 (y-1)
= 2 (x - 1)^2
$$

となり，

$$
\lim_{(x,y)\to (1,1)} \frac{g(x,y)}{|(x,y)-(1,1)|}
=
\lim_{(x,y)\to (1,1)} \frac{2 (x - 1)^2}{\sqrt{(x-1)^2+(y-1)^2}}
=
\lim_{x\to 1} 2 |x - 1|^2
=
0
$$

を満たす。

```{glue:figure} fun-2d-diff-1
:name: "fig-fun-2d-diff-1"

2変数関数 $f(x,y)=2x^2-3y$ のグラフ
```
````

````{prf:example}
:label: example-2d-diff-2

関数 $f(x,y)=\min\{|x|,|y|\}$ とおくと，
$x$ 軸上，$y$ 軸上では $f(x,y)=0$ である。
したがって，$f(x,y)$ は $(x,y)=(0,0)$ において偏微分可能であり，偏微分係数は

$$
\frac{\partial f(0,0)}{\partial x} = \frac{\partial f(0,0)}{\partial y} = 0
$$

である。しかし， $\displaystyle g(x,y) = f(x,y) - f(0,0) - \frac{\partial f(0,0)}{\partial x} (x-0) - \frac{\partial f(0,0)}{\partial y} (y-0)$ とおくと

$$
g(x,y) = \min\{|x|,|y|\}
$$

となり，

$$
\lim_{(x,y)\to (0,0)} \frac{g(x,y)}{|(x,y)-(0,0)|}
=
\lim_{(x,y)\to (0,0)} \frac{\min\{|x|,|y|\}}{\sqrt{x^2+y^2}}
=
\lim_{x\to 0} \frac{|x|}{\sqrt{2x^2}}
=
\frac{1}{\sqrt{2}}
\neq 0
$$

となるから，$f(x,y)$ は $(x,y)=(0,0)$ において全微分可能でない。

```{glue:figure} fun-2d-diff-2
:name: "fig-fun-2d-diff-2"

2変数関数 $f(x,y)=\min\{|x|,|y|\}$ のグラフ
```
````

## 多変数の微分

````{prf:definition}
:label: def-multi-diff

$n$ 変数関数 $f: \mathbb{R}^n\to\mathbb{R}$ が $\boldsymbol{x}=\boldsymbol{a}\in\mathbb{R}^n$ およびその近傍で定義されているとする。
このとき，$f(\boldsymbol{x})$ が $\boldsymbol{x}=\boldsymbol{a}$ において微分可能（全微分可能）であるとは，以下を満たすような定数ベクトル $\boldsymbol{\alpha}$ が存在することである。

$$
f(\boldsymbol{x}) = f(\boldsymbol{a}) + \boldsymbol{\alpha}^\top (\boldsymbol{x}-\boldsymbol{a}) + g(\boldsymbol{x})
\quad\text{かつ}\quad
\lim_{\boldsymbol{x}\to \boldsymbol{a}} \frac{g(\boldsymbol{x})}{|\boldsymbol{x}-\boldsymbol{a}|} = 0
$$
````

上の条件を満たすベクトル $\boldsymbol{\alpha}\in\mathbb{R}^n$ を $\boldsymbol{x}=\boldsymbol{a}$ における全微分係数とよび，$f'(\boldsymbol{a})$ あるいは $\displaystyle\frac{d f(\boldsymbol{a})}{d\boldsymbol{x}}$ などと書く。

全微分係数は勾配ともよばれる。関数 $f(\boldsymbol{x})$ 

全微分係数ベクトルは，偏微分係数を並べたものである。つまり，以下が成り立つ。

$$
f'(\boldsymbol{a}) = \biggl( \frac{\partial f(\boldsymbol{a})}{\partial x_1}, \frac{\partial f(\boldsymbol{a})}{\partial x_2}, \cdots, \frac{\partial f(\boldsymbol{a})}{\partial x_n} \biggr)
$$

## 1変数のテイラー展開

$1$ 変数関数 $f(x)$ に対して， 点 $x=a$ における $m$ 次の **テイラー多項式** $P_m(x)$ は

$$
P_m(x) = f(a) + f'(a) (x-a) + \frac{f''(a)}{2!} (x-a)^2 + \frac{f'''(a)}{3!} (x-a)^3 + \cdots + \frac{f^{(m)}(a)}{m!} (x-a)^m
$$

である。$\sum$ を用いて書き直すと，

$$
P_m(x) = \sum_{i=0}^{m} \frac{f^{(i)}(a)}{i!} (x - a)^i
$$

テイラーの定理は，$f(x)$ がテイラー多項式で近似できる，ということを主張する。
すなわち，$x=a$ の近傍において

$$
f(x) \approx P_m(x)
$$

関数 $f(x)$ と $m$ 次のテイラー多項式のズレ（差）を **剰余項** とよび $R_m(x) = f(x) - P_m(x)$ と書く。
このとき，

$$
\lim_{x\to a} \frac{R_m(x)}{(x-a)^{m}} = 0
$$

が成り立つ。

テイラー展開を以下のように書くこともできる。

$$
f(x+h) = \sum_{i=0}^{m} \frac{f^{(i)}(x)}{i!} h^i + O(h^{m+1})
$$

ただし， $O(h^{m+1})$ は $h\to 0$ において $0$ に収束する関数である。
$m$ が大きいほど急速に $0$ に収束する。

````{prf:example}
:label: example-1d-taylor

関数 $f(x)=\sin(x)$ とおくと，点 $x=-\frac{\pi}{3}$ における関数 $f(x)$ の $m$ 次のテイラー多項式 $P_m(x)$ は

$$
\begin{align*}
P_0(x)
&= f\left(-\frac{\pi}{3}\right) = \sin\left(-\frac{\pi}{3}\right)
\\
&= - \frac{\sqrt{3}}{2}
\\
P_1(x)
&= f\left(-\frac{\pi}{3}\right) + f'\left(-\frac{\pi}{3}\right) (x+1)
\\
&= - \frac{\sqrt{3}}{2} + \frac{1}{2} (x+1)
\\
P_2(x)
&= f\left(-\frac{\pi}{3}\right) + f'\left(-\frac{\pi}{3}\right) (x+1) + \frac{1}{2} f''\left(-\frac{\pi}{3}\right) (x+1)^2
\\
&= - \frac{\sqrt{3}}{2} + \frac{1}{2} (x+1) + \frac{\sqrt{3}}{4} (x+1)^2
\\
\cdots
\end{align*}
$$

```{figure} ./images/taylor1.gif
:alt: $y=\sin x$ の $x=-1$ におけるテイラー多項式
:width: 400px
:align: center
$y=\sin x$ の $x=-1$ におけるテイラー多項式
```

````
















(sec-2d-taylor)=
## 2変数のテイラー展開(1)

$2$ 変数関数 $f(x,y)$ に対して， 点 $(x,y)=(a,b)$ における $1$ 次のテイラー多項式 $P_1(x,y)$ は

$$
\begin{align*}
P_1(x,y)
&= f(a,b) + f_x(a,b) (x-a) + f_y(a,b) (y-b) \\
&= f(a,b) + \Bigl( (x-a) \frac{\partial}{\partial x} + (y-b) \frac{\partial}{\partial y} \Bigr) f(a,b)
\end{align*}
$$

演算子 $\displaystyle \Bigl( (x-a) \frac{\partial}{\partial x} + (y-b) \frac{\partial}{\partial y} \Bigr)$ を $D$ と書くと以下のように簡潔に書き替えることができる。

$$
P_1(x,y) = f(a,b) + D f(a,b)
$$

同様に，点 $(x,y)=(a,b)$ における関数 $f(x,y)$ の $2$ 次のテイラー多項式 $P_2(x,y)$ は

$$
\begin{align*}
P_2(x,y)
&= P_1(x,y) + \frac{1}{2} f_{xx}(a,b) (x-a)^2 + f_{xy}(a,b) (x-a)(y-b) + \frac{1}{2} f_{yy}(a,b) (y-b)^2 \\
&= P_1(x,y) + \frac{1}{2} \Bigl( (x-a) \frac{\partial}{\partial x} + (y-b) \frac{\partial}{\partial y} \Bigr)^2 f(a,b) \\
&= P_1(x,y) + \frac{1}{2} D^2 f(a,b) \\
&= f(a,b) + D f(a,b) + \frac{1}{2} D^2 f(a,b) \\
&= \sum_{i=0}^{2} \frac{D^i f(a,b)}{i!}
\end{align*}
$$

同様に，点 $(x,y)=(a,b)$ における関数 $f(x,y)$ の $m$ 次のテイラー多項式 $P_m(x,y)$ は

$$
P_m(x,y) = \sum_{i=0}^{m} \frac{D^i f(a,b)}{i!}
$$

$1$ 変数のときと同様に，剰余項 $R_m(x,y) = f(x,y) - P_m(x,y)$ について

$$
\lim_{(x,y)\to(a,b)} \frac{R_m(x,y)}{|(x,y)-(a,b)|^{m}} = 0
$$

が成り立つ。したがって，テイラー多項式 $P_m(x,y)$ は関数 $f(x,y)$ の近似である。

$$
f(x,y) \approx P_m(x,y)
$$

テイラー展開を以下のように書くこともできる。

$$
f(x+h,y+k) = \sum_{i=0}^{m} \frac{f^{(i)}(x,y)}{i!} h^i + O(|\sqrt{h^2+k^2}|^{m+1})
$$

## 2変数のテイラー展開(2)

関数 $f(x,y)$ について，$y$ が $x$ の関数である場合を考える。

このとき，新たな関数 $g(x)$ を用いて

$$
g(x) = f(x, y(x))
$$

と書くことができる。

$g(x)$ を $x$ で微分すると，

$$
\begin{align*}
g'(x)
&= \frac{d}{dx} f(x, y(x))
\\
&= \frac{\partial f(x,y(x))}{\partial x}
   + \frac{\partial f(x,y(x))}{\partial y} y'(x)
\\
&= \Bigl( \frac{\partial}{\partial x}
          + y'(x) \frac{\partial}{\partial y} \Bigr) f(x,y(x))
\end{align*}
$$

これを利用すると， $f(x, y(x))$ の点 $x=a$ におけるテイラー展開は以下のようになる。

$$
\begin{align*}
f(x, y(x))
&= \biggl\{ \sum_{i=0}^{m} \frac{g^{(i)}(a)}{i!} (x-a)^i \biggr\} + R_m(x,y(x))
\\
&= \biggl\{ \sum_{i=0}^{m} \frac{(x-a)^i}{i!} \Bigl( \frac{\partial}{\partial x} + y'(x) \frac{\partial}{\partial y} \Bigr)^i f(a,y(a)) \biggr\} + R_m(x,y(x))
\end{align*}
$$

ただし $R_m(x,y(x))$ は $m$ 次の剰余項である。
