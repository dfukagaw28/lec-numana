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

# ルンゲ・クッタ法の導出

## はじめに

本節では， [](../lec11/ode.md) で扱った [ルンゲ・クッタ法（Runge-Kutta method）](sec:rungekutta) の導出を行う。
本節の内容は，主に {cite:t}`yamamoto_2003` を参考にした。

次の初期値問題を考える。

$$
& y' = f(x,y)
\\& y(a) = y_0
$$

ただし，$\displaystyle y'=\frac{dy}{dx}$ である。

この初期値問題を「数値的に解く」ことが目的である。
このことについて，おさらいをしておく。

変数 $x$ が区間 $[a,b]$ を動くとき，$x$ の関数 $y(x)$ を求めたい。
ただし，解析的に解くことが難しい場合もあるため，ここでは数値的に解くことを考える。
初期値問題を数値的に解く，というのは，関数を数値的に定めることになる。
具体的には，$x$ が動く区間を $n$ 等分したうえで $x_0=a$，$x_n=b$ とおき，分点

$$
x_0 = a, \quad
x_1 = a + h, \quad
x_2 = a + 2h, \quad
\cdots, \quad
x_n = a + nh = b
$$

における $y$ の値 $y_i = y(x_i)$ を求めればよい。
区間の幅 $b-a$ を $n$ 等分した $\displaystyle h=\frac{b-a}{n}$ を **刻み幅**（step width） とよぶ。

````{prf:example}
:label: example-ode-1

たとえば区間 $[0,1]$ において $y'=3x^2$，$y(0)=1$ を解いてみよう。

解析的に解くと $y(x)=x^3+1$ である。

数値的に解く場合は，たとえば刻み幅を $h=0.1$ として $x=0,0.1,0.2,\cdots,1.0$ のときの関数値 $y=1,1.01,1.08,\ldots,2$ を求めることができればよい。
ただし，一定の数値誤差も許容する。
````

真の値 $y_i=y(x_i)$ に対して，数値解を $Y_i$ と書くことにする。

初期値は与えられているため，$Y_0=y_0$ がつねに成立すると考えてよい。
数値解を求める際，典型的には，$x=x_0$ に近い方から順に $Y_1,Y_2,\cdots,Y_n$ を求めることになる。

ある非負整数 $i$ について，$Y_0,Y_1,\cdots,Y_i$ までの数値解が与えられたとき，$Y_{i+1}$ を求める状況を考えよう。
このとき，数値解 $Y_0,Y_1,\cdots,Y_i$ のうち直前の $Y_i$ だけを用いる方法を **1段法**（single-step method）とよぶ。
1段法のうち，$f$ の導関数を必要としないものを **Runge-Kutta 型公式** とよぶ。
Runge-Kutta 型公式は，一般に，以下の形を持つ。

$$
& Y_0 = y_0,
\\
& Y_{i+1} = Y_{i} + h \cdot \Phi(x_{i},Y_{i}) \quad (i \geq 0)
$$ (rk-formula)

## Runge-Kutta 型公式の次数

数値解 $Y_i$ は，真の解 $y_i=y(x_i)$ の近似値である。

近似解 $Y_i$ の **大域離散化誤差** とは，

$$
\max_{i}\{|y(x_{i}) - Y_{i}|\}
$$

のことをいう。

また，$x \in [a,b-h]$ に対して，

$$
\tau(x,h) = \frac{y(x+h) - y(x)}{h} - \Phi(x,y(x))
$$

を **局所打ち切り誤差** といい，

$$
\tau(h) = \max_{x \in [a,b-h]} |\tau(x,h)|
$$

を **大域打ち切り誤差** という。

大域打ち切り誤差 $\tau(h)$ が正の整数 $m$ について

$$
\tau(h) = O(h^m)
$$

を満たすとき，そのような $m$ のうち最大の $m$ について，式 {eq}`rk-formula` は $m$ 次の公式であるという。

## Euler の公式

Euler の公式は，1 次の Runge-Kutta 型公式である。
このことを確認しよう。

まず，Euler の公式は以下の形で与えられる。

$$
& Y_0 = y_0,
\\
& Y_{i+1} = Y_{i} + h \cdot f(x_{i},Y_{i}) \quad (i \geq 0)
$$ (euler-formula)

すなわち，

$$
\Phi(x,y) = f(x,y)
$$

である。

関数 $y(x)$ が $C^2$ 級ならば Taylor の定理より $\xi \in (x,x+h)$ が存在して

$$
y(x+h)
  &= y(x) + h y'(x) + \frac{1}{2} h^2 y''(\xi)
\\
  &= y(x) + h f(x,y(x)) + \frac{1}{2} h^2 \frac{d}{dx} f(\xi,y(\xi))
$$

$$
\therefore
  \frac{y(x+h) - y(x)}{h} - f(x,y(x))
  =  \frac{h}{2} \cdot \frac{d}{dx} f(\xi,y(\xi))
$$

したがって，

$$
& \tau(x,h) = \frac{h}{2} \cdot \frac{d}{dx} f(\xi,y(\xi))
\\
& \therefore \tau(h)
  = \max_{x\in[a,b-h]} |\tau(x,h)|
  \leq \frac{h}{2} \cdot \max_{x\in[a,b]} \left| \frac{d}{dx} f(x,y(x)) \right|
  = O(h)
$$

以上より，Euler の公式は，1 次の Runge-Kutta 型公式である。

## Heun 法

以下では，Euler 公式による近似値を $\bar{Y}_{i+1}$ と書くことにする。すなわち，

$$
  \bar{Y}_{i+1} = Y_{i} + h f(x_{i},Y_{i})
$$ (y-bar-formula)

{numref}`ode-euler` に見るように，Euler 公式では，区間 $(x_{i},x_{i+1})$ における関数 $y(x)$ を直線で近似する。

```{figure} ./images/euler.png
---
name: ode-euler
---
Euler 法（1次の Runge Kutta 法）
```

しかし，多くの場合，関数は直線ではなく曲線であり，傾きは一定ではない。
$h$ が十分に小さいとき，更新式 {eq}`rk-formula` における $\Phi(x_{i},Y_{i})$ の値を $y'(x_{i})$ と $y'(x_{i+1})$ の中間に設定するのがよい。
$y'(x_{i+1})=f(x_{i+1},y_{i+1})$ の真の値は不明であるため，代わりに $f(x_{i+1},\bar{Y}_{i+1})$ を用いたときの $y(x_{i+1})$ の予測値を $Y^*_{i+1}$ と書くことにする。すなわち，

$$
  Y^*_{i+1} = Y_i + h f(x_{i+1},\bar{Y}_{i+1})
$$ (y-ast-formula)

この $\bar{Y}_{i+1}$ と $Y^*_{i+1}$ の算術平均を $Y_{i+1}$ とするのが Heun 法（ホイン法）である。

すなわち，

$$
  Y_{i+1} = \frac{1}{2} \left\{ \bar{Y}_{i+1} + Y^*_{i+1} \right\}
$$ (heun-formula)

{numref}`ode-heun` は， Heun 法の考え方を示す図である。

```{figure} ./images/heun.png
---
name: ode-heun
---
Heun 法（2次の Runge Kutta 法）
```

### Heun 法の誤差

Heun 公式は，2 次の Runge-Kutta 型公式であり，$\tau(x) = O(h^2)$ を満たす。
以下ではこれを確かめる。

まず $\Phi(x,y)$ を求める。Runge-Kutta 型公式の定義によれば

$$
\Phi(x_i,Y_i)
= \frac{Y_{i+1} - Y_i}{h}
$$

式 {eq}`heun-formula` より

$$
\Phi(x_i,Y_i)
= \frac{1}{2} \biggl\{ \frac{\bar{Y}_{i+1} - Y_i}{h} + \frac{Y^*_{i+1} - Y_i}{h} \biggr\}
$$

さらに式 {eq}`y-bar-formula` および式 {eq}`y-ast-formula` より

$$
\Phi(x_i,Y_i)
= \frac{1}{2} \bigl\{ f(x_{i},Y_{i}) + f(x_{i+1},\bar{Y}_{i+1}) \bigr\}
$$

ここで， [2 変数の Taylor の定理](sec-2d-taylor) により

$$
  & f(x_{i+1},\bar{Y}_{i+1})
    = f(x_i+h,Y_i + h f(x_i,Y_i))
\\&\quad
    = f(x_i,Y_i)
    + \sum_{i=1}^{2} \frac{h^i}{i!} \biggl( \frac{\partial}{\partial x} + f(x_i,Y_i) \frac{\partial}{\partial y} \biggr)^i f(x_i,Y_i)
    + O(h^{3})
$$

よって

$$
\Phi(x_i,Y_i)
    = f(x_i,Y_i)
    + \frac{1}{2} \sum_{i=1}^{2} \frac{h^i}{i!} \biggl( \frac{\partial}{\partial x} + f(x_i,Y_i) \frac{\partial}{\partial y} \biggr)^i f(x_i,Y_i)
    + O(h^{3})
$$

つまり

$$
\Phi(x,y)
&= f(x,y)
 + \frac{1}{2} \sum_{i=1}^{2} \frac{h^i}{i!} \biggl( \frac{\partial}{\partial x} + f(x,y) \frac{\partial}{\partial y} \biggr)^i f(x,y)
 + O(h^{3})
\\
&= f(x,y)
 + \frac{1}{2} \sum_{i=1}^{2} \frac{h^i}{i!} \biggl( \frac{d}{dx} \biggr)^i f(x,y)
 + O(h^{3})
\\
&= f(x,y)
 + \frac{h}{2} \frac{d}{dx} f(x,y)
 + \frac{h^2}{4} \frac{d^2}{dx^2} f(x,y)
 + O(h^{3})
$$ (heun-phi)

次に，$\displaystyle\frac{Y_{i+1} - Y_i}{h}$ を求める。
Taylor の定理により，

$$
\frac{y(x+h) - y(x)}{h}
&= \frac{1}{h} \biggl\{ \sum_{i=1}^{3} \frac{h^i}{i!} y^{(i)}(x) + O(h^{4}) \biggr\}
\\
&= \sum_{i=1}^{3} \frac{h^{i-1}}{i!} y^{(i)}(x) + O(h^{3})
\\
&= \sum_{i=1}^{3} \frac{h^{i-1}}{i!} \Bigl( \frac{d}{dx} \Bigr)^{i-1} f(x,y) + O(h^{3})
\\
&= f(x,y)
 + \frac{h}{2} \frac{d}{dx} f(x,y)
 + \frac{h^{2}}{6} \frac{d^2}{dx^2} f(x,y)
 + O(h^{3})
$$ (heun-y-slope)

式 {eq}`heun-phi` と式 {eq}`heun-y-slope` より，局所打切り誤差は

$$
\tau(x,h)
&= \frac{y(x+h) - y(x)}{h} - \Phi(x,y)
\\
&= -\frac{h^2}{12} \frac{d^2}{dx^2} f(x,y) + O(h^{3})
$$

$$
\tau(h)
&= \max_x |\tau(x,h)|
\\
&\le h^2 \biggl| \max_x \frac{d^2}{dx^2} f(x,y) \biggr| + O(h^3)
\\
&=
O(h^2)
$$

以上より，Heun の公式は，2 次の Runge-Kutta 型公式である。

## 修正 Euler 法

2次の Runge-Kutta 公式は，Heun 法だけではない。

定数 $\alpha,\beta,p,q$ に対して，関数 $\Phi(x_i,Y_i)$ を以下のように定義する。

$$
& \Phi(x_i,Y_i) = \alpha k_1 + \beta k_2
\\
& k_1 = f(x_i,Y_i)
\\
& k_2 = f(x_i+ph,Y_i+qhk_1)
$$ (rk2-formula)

Taylor の定理より，

$$
k_2 = f(x_i,Y_i) + ph f_x(x_i,Y_i) + qh f(x_i,Y_i) f_y(x_i,Y_i) + O(h^2)
$$

よって

$$
Y_{i+1}
&= Y_{i} + h \cdot \Phi(x_i,Y_i)
\\
&= Y_{i} + h \cdot (\alpha k_1 + \beta k_2)
\\
&= Y_{i} + \alpha h f + \beta h \left\{ f + ph f_x + qh f f_y + O(h^2) \right\}
\\
&= Y_{i} + h (\alpha + \beta) f + h^2 \beta \left\{ p f_x + q f f_y \right\} + O(h^3)
\\
&= Y_{i} + h \left\{ (\alpha + \beta) f + h \beta \left\{ p f_x + q f f_y \right\} + O(h^2) \right\}
$$

つまり，

$$
\Phi(x,y) = (\alpha + \beta) f(x,y) + h \beta \left\{ p f_x(x,y) + q f(x,y) f_y(x,y) \right\} + O(h^2)
$$

である。

したがって，局所打ち切り誤差は

$$
\tau(x,h)
&= \frac{y(x+h) - y(x)}{h} - \Phi(x,y(x))
\\
&= f(x,y(x))
   + \frac{h}{2} \cdot f'(x,y(x))
   + O(h^2)
   - \Phi(x,y(x))
\\
&= f(x,y(x))
   + \frac{h}{2} \cdot \{ f_x + f f_y \}
   - \Phi(x,y(x))
   + O(h^2)
\\
&= (1 - \alpha - \beta) f
   + h \left\{ \frac{1}{2} - \beta p \right\} f_x
   + h \left\{ \frac{1}{2} - \beta q \right\} f f_y
   + O(h^2)
$$

ここで，

$$
\alpha + \beta = 1, \quad
\beta p = \beta q = \frac{1}{2}
$$

を満たすように $\alpha,\beta,p,q$ を選べば，

$$
\tau(x)
  = \max_{x \in [a,b-h]}|\tau(x,h)|
  = O(h^2)
$$

となり，このとき，式 {eq}`rk2-formula` は 2次の Runge-Kutta 型公式となる。

特に $\alpha=\beta=\frac{1}{2}$，$p=q=1$ の場合は Heun 法そのものである。

また，$\alpha=0,\beta=1,p=q=\frac{1}{2}$ の場合は **修正 Euler 法** とよばれる。

## 3 次の Runge-Kutta 公式

3 次の Runge-Kutta 公式は一般に以下の形式を持つ。

$$
  & Y_{i+1} = Y_i + h \Phi(x_i,Y_i),\quad Y_0=y_0
\\
  & \Phi(x_i,Y_i) = \alpha k_1 + \beta k_2 + \gamma k_3
\\
  & k_1 = f(x_i,Y_i)
\\
  & k_2 = f(x_i+ph,Y_i+qhk_1)
\\
  & k_3 = f(x_i+rh,Y_i+shk_2+thk_1)
$$

たとえば，以下のような値を選べば $\tau(x,h)=O(h^3)$ となる。

$$
  \alpha = \gamma = \frac{1}{6},\quad
  \beta = \frac{2}{3},\quad
  p = q = \frac{1}{2},\quad
  r = 1,\quad
  s = 2,\quad
  t = -1
$$

## 4 次の Runge-Kutta 公式

4 次の Runge-Kutta 公式は一般に以下の形式を持つ。

$$
  & Y_{i+1} = Y_i + h \Phi(x_i,Y_i)
\\
  & \Phi(x_i, Y_i) = \lambda_1 k_1 + \lambda_2 k_2 + \lambda_3 k_3 + \lambda_4 k_4
\\
  & k_1 = f(x_i, Y_i)
\\
  & k_2 = f(x_i + \alpha_2 h, Y_i + \beta_2 k_1 h)
\\
  & k_3 = f(x_i + \alpha_3 h, Y_i + (\beta_3 k_1 + \gamma_3 k_2) h)
\\
  & k_4 = f(x_i + \alpha_4 h, Y_i + (\beta_4 k_1 + \gamma_4 k_2 + \delta_4 k_3) h)
$$

$\tau(x,h)=O(h^4)$ となるように $\lambda_1,\lambda_2,\lambda_3,\lambda_4,\alpha_2,\alpha_3,\alpha_4,\beta_2,\beta_3,\beta_4,\gamma_3,\gamma_4,\delta_4$ を選べばよい。
Runge は以下のように選んだ。

$$
  & (\lambda_1,\lambda_2,\lambda_3,\lambda_4)
    = \left( \frac{1}{6}, \frac{1}{3}, \frac{1}{3}, \frac{1}{6} \right),
    \quad
    (\alpha_2,\beta_2)
    = \left( \frac{1}{2}, \frac{1}{2} \right),
\\& (\alpha_3,\beta_3,\gamma_3)
    = \left( \frac{1}{2}, 0, \frac{1}{2} \right),
    \quad
    (\alpha_4,\beta_4,\gamma_4,\delta_4)
    = \left( 1, 0, 0, 1 \right)
$$

一方，Kutta は以下のように選んだ。

$$
  & (\lambda_1,\lambda_2,\lambda_3,\lambda_4)
    = \left( \frac{1}{8}, \frac{3}{8}, \frac{3}{8}, \frac{1}{8} \right),
    \quad
    (\alpha_2,\beta_2)
    = \left( \frac{1}{3}, \frac{1}{3} \right),
\\& (\alpha_3,\beta_3,\gamma_3)
    = \left( \frac{2}{3}, -\frac{1}{3}, 1 \right),
    \quad
    (\alpha_4,\beta_4,\gamma_4,\delta_4)
    = \left( 1, 1, -1, 1 \right)
$$

いずれの場合も $\tau(x,h)=O(h^4)$ となることが確認できる。
