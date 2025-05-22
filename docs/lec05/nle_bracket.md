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

# 非線形方程式の数値解法 (1) 囲い込み法

1次方程式は，線形方程式ともよばれる。
線形方程式は，代数的な解（四則演算と冪根を求める操作を有限回組み合わせて得られるもの）が存在し，ガウスの消去法などの方法で解くことができる。
変数の個数と等式の個数，およびそれらの内容によっては解が存在しないこともあるが，解の存在判定も求解と同様の手順で行うことができる。

一方，非線形方程式については複雑であり数値解法に頼らざるを得ない。
例外的に，代数方程式（多項式関数 $f(x)$ について $f(x)=0$ と表される方程式）のように解の存在や個数が保証されるものについて特別な解法が知られている場合もあるが，一般の非線形方程式について応用できるわけではない。

与えられた非線形方程式 $f(x)=0$ について，もし，解の存在する範囲がある程度分かっている場合には，線型方程式と同様に，次のような解法が考えられる。

- 囲い込み法（bracketing method）: 関数 $f$ の値を利用して，解の存在範囲を繰り返し絞り込む。たとえば，二分法，ブレント法など。
- 反復法（iterative method）: 方程式 $f(x)=0$ と同値な方程式 $x = g(x)$ で，関数 $g$ が（解の近傍において）縮小写像となるものを見つける。数列 $\{x_{n}\}$ を漸化式 $x_{n+1}=g(x_{n})$ によって定め，近似解を求める。たとえば，ニュートン・ラフソン法など。

まず今回は囲い込み法を紹介し，次回は反復法を扱う。

## 二分法

二分法の基礎となるのは中間値の定理である。

```{prf:theorem} 中間値の定理
:label: thm-intermediate

関数 $f(x)$ は閉区間 $[a,b]$ で連続であるとする。
また，$\min\{f(a),f(b)\}=m$，$\max\{f(a),f(b)\}=M$ とおく。
このとき，任意の数 $y\in(m,M)$ に対してある数 $c\in(a,b)$ が存在して $f(c)=y$ を満たす。
```

特に $y=0$ としたとき以下の系を得る。

```{prf:corollary}
:label: cor-intermediate

関数 $f(x)$ は閉区間 $[a,b]$ で連続であるとする。
$f(a)f(b)<0$ のとき，$f(x)=0$ の解（のうち少なくとも1つ）が区間 $(a,b)$ に存在する。
```

```{code-cell}
:tags: [remove-input]

import numpy as np
import matplotlib.pyplot as plt

a,b,c=[-1.5,0.2,1.3]
f = lambda x: (x - a) * (x - b) * (x - c)
X = np.linspace(-2.5,2.5,101)
Y = f(X)
plt.plot(X, Y, color='black')
plt.axhline(y=0, linestyle=':', color='red')
plt.axvline(x=-2, linestyle=':', color='blue')
plt.axvline(x=+2, linestyle=':', color='blue')
plt.scatter([-2,2], [f(-2),f(2)], color='blue')
plt.scatter([a,b,c], [f(a),f(b),f(c)], color='red')
plt.text(-2.05,-0.25,'a',ha='right',va='top',color='blue')
plt.text(+2.05,-0.25,'b',ha='left',va='top',color='blue')
plt.show()
```

ある閉区間 $[a,b]$ 上で定義された連続関数 $f$ について $f(a)f(b)<0$ が成り立つとき， $a<c<b$ を満たす任意の $c$ について以下のいずれかが成り立つ。

1. $f(a)f(c)<0$
   - このとき， {prf:ref}`cor-intermediate` により，区間 $(a,c)$ に方程式 $f(x)=0$ の解が存在する。
2. $f(a)f(c)>0$
   - $f(a)f(b)<0$ より $f(c)f(b)<0$ が成り立つ。
   - このとき， {prf:ref}`cor-intermediate` により，区間 $(c,b)$ に方程式 $f(x)=0$ の解が存在する。
3. $f(a)f(c)=0$
   - このとき，$f(c)=0$ である。つまり，$x=c$ は方程式 $f(x)=0$ の解である。

いずれの場合でも，解の存在する区間の幅を半分に絞り込むか，もしくは，解そのものを見つけることができる。
この操作は繰り返すことができ，区間の幅を繰り返し半減させることによって解を急速に絞り込むことができる。

このような解の探索方法を **二分法**（bisection method）とよぶ。

```{code-cell}
:tags: [remove-input]

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib_fontja
import numpy as np

cmap = plt.get_cmap('tab10')
colors = cmap.colors

def draw_bracket(k, a, b, y, color):
    plt.gca().add_patch(patches.FancyArrowPatch(
        (a, y),
        (b, y),
        arrowstyle='<->',
        shrinkA=0,
        shrinkB=0,
        mutation_scale=10,
        color=color,
    ))
    plt.text(max(a, b) + 0.1, y, f'iter={k}', ha='left', va='center', color=color)

def draw_secant(x1, x2, y1, y2):
    if y1 == y2:
        plt.axhline(y=y1, linestyle=':', linewidth=1, color='black')
    else:
        x3 = x2 - y2 * (x2 - x1) / (y2 - y1)
        y3 = 0
        xmin = min(x1, x2, x3)
        xmax = max(x1, x2, x3)
        ymin = min(y1, y2, y3)
        ymax = max(y1, y2, y3)
        plt.plot([xmin, xmax], [ymin, ymax], linestyle=':', linewidth=1, color='black')

def show_bracketing(a, b, f, repeat=5, margin=0.5, yoffset=-5, ystep=0.6, method='bisect'):
    X = np.linspace(min(a, b) - margin, max(a, b) + margin, 101)
    Y = f(X)

    # y = 0
    plt.axhline(y=0, linestyle='-', linewidth=.7, color='black')

    for k in range(1, 1 + repeat):
        y = yoffset - k * ystep
        color = colors[k]
        draw_bracket(k, a, b, y, color)

        #plt.scatter([a, b], [f(a), f(b)], color='red')
        plt.plot([a, a], [y, f(a)], linestyle=':', color=color)
        plt.plot([b, b], [y, f(b)], linestyle=':', color=color)

        if method == 'bisect':
            c = (a + b) / 2
            a, b = (a, c) if f(a) * f(c) < 0 else (c, b)
        elif method == 'secant':
            draw_secant(a, b, f(a), f(b))
            c = b - f(b) * (b - a) / (f(b) - f(a))
            a, b = b, c
        elif method == 'regula-falsi':
            draw_secant(a, b, f(a), f(b))
            c = b - f(b) * (b - a) / (f(b) - f(a))
            a, b = (a, c) if f(a) * f(c) < 0 else (c, b)
        else:
            raise Exception

    plt.plot(X, Y, color='black')
    plt.title(method)
    plt.show()
```

```{code-cell}
:tags: [remove-input]

f = lambda x: 0.02 * (x + 3) ** 4 - 2
a, b = -2, 1

show_bracketing(a, b, f, repeat=4, yoffset=-2.2, ystep=0.35, method='bisect')
```

二分法の疑似コードを以下に示す。

```{prf:algorithm} 二分法
:label: alg:bisection

**Inputs** 区間 $[a,b]$，精度 $\varepsilon>0$，関数 $f:\mathbb{R}\to\mathbb{R}$

**Outputs** 方程式 $f(x)=0$ の 解 $x\in\mathbb{R}$

1. while $|a - b| \geq \varepsilon$
   1. $\displaystyle c \leftarrow \frac{a+b}{2}$
   2. if $f(a) f(c) < 0$ then
      1. $b \leftarrow c$
   3. else
      1. $a \leftarrow c$
2. $\displaystyle c \leftarrow \frac{a+b}{2}$
3. return $c$
```

二分法によって求めることのできる解の数値誤差について考えてみよう。

探索すべき区間の初期値を $[a_0,b_0]=[a,b]$ とする。
$1$ 回目の反復によって区間は $[a_0,b_0]$ から $\displaystyle \Bigl[ a_0, \frac{a_0+b_0}{2} \Bigr]$ または $\displaystyle \Bigl[ \frac{a_0+b_0}{2}, b_0 \Bigr]$ に更新される。これを $[a_1,b_1]$ と書くことにする。
同様に，$n$ 回目の反復によって得られる区間を $[a_n,b_n]$ とする。
区間幅 $b_n-a_n$ は $\displaystyle \frac{b-a}{2^n}$ と等しい。
真の解を $x=\alpha$，近似解を $x=c_n=\frac{a_n+b_n}{2}$ とすると，

$$
a_n \leq \alpha \leq b_n
$$

したがって

$$
\begin{align*}
& c_n - \alpha \leq c_n - a_n = \frac{a_n+b_n}{2} - a_n = \frac{b_n-a_n}{2}
\\
& c_n - \alpha \geq c_n - b_n = \frac{a_n+b_n}{2} - b_n = - \frac{b_n-a_n}{2}
\end{align*}
$$

すなわち，真の解 $\alpha$ と近似解 $c_n$ との絶対誤差について以下が成り立つ。

$$
|c_n-\alpha| \leq \frac{b_n-a_n}{2} = \frac{b-a}{2^{n+1}}
$$

すなわち誤差限界は $\displaystyle \frac{b-a}{2^{n+1}}$ となる。
誤差限界を $\varepsilon$ 以下とするためには，

$$
\log_2 \frac{b-a}{\varepsilon} - 1
$$

回の反復で十分である。
たとえば初期区間幅を $b-a=10$，誤差限界を $\varepsilon=10^{-8}$ とすると， $\log_2{10^9}-1\approx 29$ 回の反復で十分である（関数によって，より少ない場合もあり得る）。

二分法を適用するためには，まず最初に $f(a)f(b)<0$ を満たすような区間 $[a,b]$ を見つけ出す必要がある。
もし $f(x)=x^3+bx^2+cx+d$ のように $\displaystyle\lim_{x\to\infty}f(x)$ と $\displaystyle\lim_{x\to-\infty}f(x)$ とで符号が異なる場合は，十分に広い区間に標本となる点 $x_0,x_1,x_2,\cdots,x_N$ をとって，$f(x_k)f(x_{k+1})<0$ を満たすような区間 $[x_k,x_{k+1}]$ を調べればよい。もし $f$ が多変数関数である場合は，探索範囲が広大となるが，符号が異なる 2 点さえ見つかれば，その 2 点を結ぶ直線探索によって解を見つけることができる。

## 例1（二分法）

たとえば関数 $f(x)=x^6+5x-4$ について方程式 $f(x)=0$ の解を求めてみよう。
解を 64 ビット浮動小数点数で近似するとして，真の解を $x=\alpha$，近似解を $x=\tilde{\alpha}$ と書くことにする。

計算機プログラムでは，しらみつぶしに探すことも有力な選択肢に成り得る。
すべての 64 ビット浮動小数点数を試せば $x=\tilde{\alpha}$ が見つかるはずである。
ただし 64 ビット浮動小数点数は $2^{64} \approx 16 \times 10^{18}$ 通り存在し，1個試すのに 1 ナノ秒（10億分の1秒）しかかからないとしても，500 年以上かかることになる。これは現実的でない。

そこで，適当な範囲 $[a,b]$ および適当な刻み幅 $h$ で $x$ を変化させて $f(x)$ が $0$ に近いものを探す方法を考える。
正整数 $N$ について $\displaystyle h=\frac{b-a}{N}$ が成り立つように選べば，候補は $N+1$ 個の点となる。

たとえば区間 $[-9,9]$ を刻み幅 $h=1$ で分割し，各点 $x_0,x_1,x_2,\cdots,x_N$ における関数値 $f(x_i)$ を計算する。各点 $(x_0,f(x_0)),(x_1,f(x_1)),\cdots,(x_N,f(x_N))$ のグラフを $xy$ 座標平面上にプロットして繋いだものを以下に示す。

```{code-cell}
"""
サンプル点における関数 f の値を求め，プロットする
"""

import numpy as np
import matplotlib.pyplot as plt

# 区間 [-9, 9] 上に 19 個の点を生成する
X = np.linspace(-9, 9, 19)

# 関数 f を定義する
f = lambda x: x ** 6 + 5 * x - 4

# 関数値を求める
Y = f(X)

# 19個の点について関数値を表示する
for i, x in enumerate(X):
    print(f'{x:+4.1f} {Y[i]:+.6e}')

# 関数 y=f(x) のグラフを描画する
plt.plot(X, Y, marker='o')
plt.axhline(0, linestyle=':', color='red')
plt.show()
```

これだけでは解は分からない。
サンプル点における $f(x)$ の値を調べると，$f(x)$ が最も $0$ に近いのは $x=1$ のときである。
このことに注目すると，解は $x=1$ の付近に絞り込まれたと考えてよいだろう。

周辺を拡大すると，次の図のようになる。

```{code-cell}
:tags: [remove-input]

import numpy as np
import matplotlib.pyplot as plt

X = np.linspace(-2,2,5)
f = lambda x: x ** 6 + 5 * x - 4
Y = f(X)
plt.plot(X, Y, marker='o')
plt.axhline(0, linestyle=':', color='red')
plt.show()
```

符号の変化から，この範囲には，少なくとも 2 個の解が存在することが分かる。

ここでは区間 $[0,1]$ に注目して，この区間における解を二分法によって求めてみよう。

```{code-cell}
"""
二分法の実行例
"""

# 二分法を繰り返し適用することによって方程式 f(x)=0 の解を求める。
# 初期区間を [a,b] とし，誤差が tol 未満になったら終了する。
def bisection(f, a, b, tol, maxiter=100):
    for iter in range(1, maxiter + 1):
        # 途中経過を表示する
        c = (a + b) / 2
        print(f'{iter:3d} 区間[a={a:.6f},c={c:.6f},b={b:.6f}]'
              + f' 幅{b-a:.6f} f(a)={f(a):.6f} f(b)={f(b):.6f}')

        if abs(a - b) < tol * 2:
            break

        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return c

# 定数や関数を定義する
a, b, tol = 0.0, 1.0, 1e-6
f = lambda x: x ** 6 + 5 * x - 4

# 二分法を実行する
x = bisection(f, a, b, tol)

# 解を表示する
print('x =', x)
```

20回の反復によって近似解 $x=0.7611188888549805$ を得ることができた。
許容誤差を `tol = 1e-6` と設定したことをふまえると，近似解は約 $0.76112$ といえる[^bisection-example]。

[^bisection-example]: 解が $0.7611188888549805$ の上下 $\pm 10^{-6}$ の区間に含まれることが保証される。すなわち $(0.7611178\cdots, 0.7611198\cdots)$ に含まれるから小数第6位までで丸めると $0.761118$ または $0.761119$ または $0.761120$ になり，定まらない。
小数第5位までで丸めると $0.76112$ となる。

## 割線法・挟み撃ち法

割線法（secant method）と挟み撃ち法（regula falsi, false position method）は二分法と同様に，区間 $[a,b]$ を二分割してどちらか一方を選ぶという操作を繰り返すが，分割点 $c$（$a<c<b$）の選び方が異なる。

区間 $[a,b]$ に解が存在すると分かっているとき，二分法では新しい点 $c$ を次のように求める。

$$
c = \frac{a + b}{2}
$$

一方，割線法や挟み撃ち法では $f(a),f(b)$ の値を用いて $c$ を次のように求める。

$$
c = b - f(b) \frac{b - a}{f(b) - f(a)}
$$

この $c$ は，$xy$ 平面において点 $(a,f(a))$ と点 $(b,f(b))$ を結ぶ直線と $x$ 軸（直線 $y=0$）との交点である。

割線法では $f(c)$ を使わず，最近求めた2点に挟まれた区間を選ぶ。
このとき，新たに選ばれた区間に解が含まれる保証はないという点に注意が必要である（したがって，割線法は「囲い込み法」ではない）。

```{code-cell}
:tags: [remove-input]

f = lambda x: 0.02 * (x + 3) ** 4 - 2
a, b = -2, 1

show_bracketing(a, b, f, repeat=4, yoffset=-2.2, ystep=0.35, method='secant')
```

挟み撃ち法では $f(c)$ の符号を用いて $[a,c]$ と $[c,b]$ のうち解を含む方の区間を選ぶ。

```{code-cell}
:tags: [remove-input]

f = lambda x: 0.02 * (x + 3) ** 4 - 2
a, b = -2, 1

show_bracketing(a, b, f, repeat=4, yoffset=-2.2, ystep=0.35, method='regula-falsi')
```

(sec:convergence-condition-bisectionn)=
## 収束条件，収束の速さ

近似解の誤差について詳しく確認してみよう。

方程式 $f(x)=0$ の真の解を $x=\alpha$ とし， $n$ 回の繰り返しによって求めたブラケットを $[a_n,b_n]$，近似解を $x_n$ とする。
また，$n$ 回目の近似解の誤差を $e_n$，誤差限界を $\varepsilon_{n}$ とおく。
つまり，$x_n = \alpha + e_n$ かつ $|e_{n}| \leq \varepsilon_{n}$ が成り立つ。

二分法では

$$
x_n = \frac{a_n+b_n}{2}
\quad \text{かつ} \quad
a_n \leq \alpha \leq b_n
\quad \text{かつ} \quad
|b_{n+1} - a_{n+1}| = \frac{1}{2} |b_{n} - a_{n}|
\quad \text{かつ} \quad
\varepsilon_{n} = \frac{1}{2} |b_{n} - a_{n}|
$$

が成り立つため，

$$
\varepsilon_{n+1} = \frac{1}{2} \varepsilon_{n}
$$ (bisect-eps)

である。

実数 $\alpha$ に収束する数列 $\{x_n\}$ に対して，ある定数 $\kappa$（$0 < \kappa < 1$）に収束する数列 $\{\kappa_n\}$ が存在して

$$
x_{n+1} - \alpha = \kappa_n (x_n - \alpha)
$$

を満たすとき，$\{ x_n \}$ は **線形収束する**（converges linearly）という。
また，ある定数 $M$（$0 < M < \infty$）が存在して

$$
|x_{n+1} - \alpha| = M |x_n - \alpha|^p
$$

を満たすとき，$\{x_n\}$ は **$p$ 次収束する**（converges with order $p$）という。

式 {eq}`bisect-eps` より，二分法は線形収束（$1$ 次収束）する。

一方，割線法は黄金比を $\phi = \frac{1+\sqrt{5}}{2} \approx 1.618$ とするとき $\phi$ 次収束することが知られている（証明は省略する）。
ただし，割線法は二分法と違って必ずしもブラケットに真の解が含まれるとは限らない。
たとえ収束が速くても収束したときの値と真の解のズレが大きければ実用には適さないかもしれない。

挟み撃ち法は割線法を改良して必ずブラケットが真の解を含むようにしたものである。
挟み撃ち法の収束の速さは二分法と割線法の中間くらい（$1$と$1.618$との間）である。

収束の速さについて考えると，二分法より挟み撃ち法の方が望ましいと言えそうである。
では，挟み撃ち法よりも良い方法があるだろうか。

実は，割線法や挟み撃ち法よりも二分法の方が収束が速くなるような病的な関数（pathological function）の存在が知られている。
そのような欠点を改良する方法として，たとえば以下のようなものが知られている。

- Ridders 法 [^Ridders1979] は $2$ 次収束することが知られている。
  ただし，Ridders 法は繰返しのたびに関数値を2回求める必要があるため，関数評価1回あたりの収束次数は $\sqrt{2} \approx 1.414$ となる。

- Brent 法（van Wijngaarden-Dekker-Brent 法） [^Brent1973] やその変種である Bus-Dekker 法 [^BusAndDekker1975] はいずれも $1$ と $1.618$ の間の収束次数を持つことが知られている。Brent 法はブラケット・二分法・逆二乗外挿を用いる手法である。Bus-Dekker 法は逆二乗外挿の代わりに双曲線外挿を用いる。実装はやや複雑であるためここでは省略する。

- TOMS アルゴリズム 748 法 [^APS1995] は Brent 法をさらに改良したものである。パラメータ $k$ を持ち，収束次数は $k=1$（収束次数は約 $1.65$）または $k=2$（収束次数は約 $1.66$）が最適とされている。

[^Ridders1979]: Ridders, C. (1979) "A new algorithm for computing a single root of a real continuous function". IEEE Transactions on Circuits and Systems. Volume 26, Issue 11, pp. 979–980. [doi:10.1109/TCS.1979.1084580](https://doi.org/10.1109/TCS.1979.1084580).

[^Brent1973]: Brent, R. P. (1973) "Algorithms for Minimization Without Derivatives". Englewood Cliffs, NJ: Prentice-Hall, 1973. Ch. 3-4.

[^BusAndDekker1975]: Bus, J. C. P., Dekker, T. J. (1975) "Two Efficient Algorithms with Guaranteed Convergence for Finding a Zero of a Function", ACM Transactions on Mathematical Software, Volume 1, Issue 4, pp. 330-345. Section 3: "Algorithm M". [doi:10.1145/355656.355659](https://doi.org/10.1145/355656.355659)

[^APS1995]: Alefeld, G. E. and Potra, F. A. and Shi, Yixun (1995) "Algorithm 748: Enclosing Zeros of Continuous Functions". ACM Transactions on Mathematical Software Volume 221, Issue 3, pp. 327-344 [doi:10.1145/210089.210111](https://doi.org/10.1145/210089.210111)

## SciPy を用いた解法

上で紹介した解法のうちいくつかは SciPy で実装されており，手軽に利用できる。
SciPy で利用可能な手法とそれらの特徴は，SciPy のドキュメント SciPy API ＞ Optimization and root finding ＞ Root finding ＞ [Scalar functions](https://docs.scipy.org/doc/scipy/reference/optimize.html#scalar-functions)にまとめられている。

たとえば，二分法は以下のようにして関数 [`scipy.optimize.bisect()`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.bisect.html) を呼び出すことによって実行できる。

```{code-cell}
"""
二分法の実行例
"""

from scipy import optimize

# 定数や関数を定義する
a, b, tol = 0.0, 1.0, 1e-6
f = lambda x: x ** 6 + 5 * x - 4

# 二分法を実行する
root, results = optimize.bisect(f, a, b, xtol=tol, full_output=True)

# 実行結果を表示する
print(results)
```

あるいは関数 [`scipy.optimize.root_scalar()`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root_scalar.html) を用いて以下のように実行できる。

```{code-cell}
"""
二分法の実行例
"""

from scipy import optimize

# 定数や関数を定義する
a, b, tol = 0.0, 1.0, 1e-6
f = lambda x: x ** 6 + 5 * x - 4

# 二分法を実行する
results = optimize.root_scalar(f, method='bisect', bracket=[a, b], xtol=tol)

# 実行結果を表示する
print(results)
```

いくつかの関数について，繰返し回数を比較してみよう。

```{code-cell}
import pandas as pd
from scipy import optimize
from sympy import Eq, sympify

methods = [
    'bisect',
    'brentq',
    'brenth',
    'ridder',
    'toms748',
]

# 定数や関数を定義する
a, b, tol = 0.0, 1.0, 1e-16
f = lambda x: x ** 6 + 5 * x - 4

# bracketing によって解を求める
display(Eq(sympify('f(x)'), f(sympify('x'))))
print(f'bracket={[a,b]} {tol=}')
log = []
for method in methods:
    results = optimize.root_scalar(f, method=method, bracket=[a, b], xtol=tol)
    log.append({
        'method': method,
        'iterations': results.iterations,
        'root': results.root,
    })
display(pd.DataFrame(log))
```

```{code-cell}
import pandas as pd
from scipy import optimize
from sympy import Eq, sympify

methods = [
    'bisect',
    'brentq',
    'brenth',
    'ridder',
    'toms748',
]

# 定数や関数を定義する
a, b, tol = -2, 1.99, 1e-16
f = lambda x: -1 / (x - 2) - 1

# bracketing によって解を求める
display(Eq(sympify('f(x)'), f(sympify('x'))))
print(f'bracket={[a,b]} {tol=}')
log = []
for method in methods:
    results = optimize.root_scalar(f, method=method, bracket=[a, b], xtol=tol)
    log.append({
        'method': method,
        'iterations': results.iterations,
        'root': results.root,
    })
display(pd.DataFrame(log))
```

関数 `f` やブラケット `bracket` の初期値によって，アルゴリズム `method` の優劣が入れ替わる場合があるものの，
一般的に `brentq` や `toms748` は `bisect` よりも繰返しが少なく済む傾向が見られる。
