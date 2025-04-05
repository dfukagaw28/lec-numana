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

# 数値計算と誤差

## 浮動小数点数の四則演算

2 つの数 $x, y$ が $\beta$ 進 $t$ 桁の浮動小数点数で与えられているとする。

```{math}
  x &= D_x \times \beta^{e_x} \\
  y &= D_y \times \beta^{e_y}
```

ただし，$D_x$ と $D_y$ はそれぞれ仮数部であり，いずれも 1 の位が最上位桁であるように正規化されているものとする。

このとき，$x$ と $y$ の和（あるいは差）は以下のように計算できる。

1. 指数部を揃える。
1. 仮数部の和（あるいは差）を求める。
1. 正規化する。

また，$x$ と $y$ の積（あるいは商）は以下のように計算できる。

1. 仮数部の積（あるいは商）を求める。
1. 指数部の和（あるいは差）を求める。
1. 正規化する。

たとえば $\beta=10$，$t=5$，$x=1.2345 \times 10^{1}$，$y=4.5678 \times 10^{-2}$ のとき，

```{math}
  x + y
  = (1.2345 + 0.0045678) \times 10^1
  = 1.2390678 \times 10^1
  \approx 1.2391 \times 10^1
```

指数部の小さい方の数について，下位桁が計算結果に与える影響が小さいことに注意されたい。
たとえば，
```{math}
  12.345 + 0.045678 \approx 12.391 \\
  12.345 + 0.045600 \approx 12.391
```
である。
このように，絶対値が大きく異なる値の加減算において，絶対値の小さい値の下位桁（あるいは全部）の情報が失われる現象を **情報落ち**（loss of trailing digits）という。

情報落ちは，例えば，多くの小さな値の和をとるときに生じやすい。
情報落ちを避けるためには，計算の順序を工夫して，たとえば，小さい値から順に足し合わせるとよい。

もう一つの例を見てみよう。
たとえば $x=\sqrt{101}$，$y=\sqrt{100}$ とする。
$x$ と $y$ をそれぞれ 10 進 5 桁の浮動小数点数で近似したものを $\tilde{x}=1.0050 \times 10^{1}$，$\tilde{y}=1.0000 \times 10^{1}$ とする。
このとき，2 数の差は

```{math}
  x - y
  \approx \tilde{x} - \tilde{y}
  = (1.0050 - 1.0000) \times 10^1
  = 0.0050 \times 10^1
  \approx 5.0000 \times 10^{-2}
```

となる。
真の値は $x-y=\sqrt{101}-\sqrt{100}=4.98756\cdots \times 10^{-2}$ であるから，計算結果の仮数部 $5.0000$ はかなり不正確である。
特に下位桁は，$x$ や $y$ では丸められた桁に由来するため意味を持たない。
つまり，実際の有効桁数（意味を持つ桁の幅）は 5 桁ではないと考えるべきである。
このように，近い値の差をとることで有効桁数が失われる現象のことを **桁落ち**（loss of significance）という。

## アルゴリズムと数値誤差

桁落ちは，数値計算で頻繁に生じる。
さきほどの計算 $\sqrt{101} - \sqrt{100}$ において桁落ちを避けるためには，たとえば，以下のように計算を工夫するとよい。
```{math}
  \sqrt{101} - \sqrt{100}
  &= \frac{1}{\sqrt{101} + \sqrt{100}} \\
  &\approx \frac{1}{(1.0050 + 1.0000) \times 10^1} \\
  &= \frac{1}{2.0050 \times 10^1} \\
  &= (1.0000 \times 10^0) \div (2.0050 \times 10^1) \\
  &= (1.0000 \div 2.0050) \times 10^{0-1} \\
  &= 0.4987531172\cdots \times 10^{-1} \\
  &\approx 4.9875 \times 10^{-2} \\
```
これは真の値にかなり近い。

この例は，以下の事実を示している。
ある二つの計算式があり，それらが数学的には互いに等価な場合であっても，計算式が示す計算方法（すなわちアルゴリズム）の違いによって数値計算の結果が大きく異なる場合がある。

もう一つ例を挙げる。
実数 $a,b,c$ に対して，二次方程式 $ax^2+bx+c=0$ の解を求めたい。
ただし $b^2-4ac>0$ とする。
解の公式を用いて

$$
x = \frac{-b \pm \sqrt{b^2-4ac}}{2a}
$$

のまま計算をすると，精度が大幅に下がることがある。
たとえば $b > 0$ かつ $|b^2|$ が $|4ac|$ よりもずっと大きいと $b^2-4ac$ の値において下位桁の情報落ちが生じ，$b^2 - 4ac \approx b$ となる。このとき解の一方はほぼゼロとなり桁落ちが生じる。

そのようなときは，まず解のうち絶対値の大きな $\displaystyle\alpha=\frac{-b-\sqrt{b^2-4ac}}{2a}$ を求める。次に，解と係数の関係 $\alpha\beta=\frac{c}{a}$ を用いてもう一つの解 $\beta$ を求める。これにより，精度の低下を防ぐことができる。

## Python を用いて誤差を確認する

上に挙げた例は Python プログラムで以下のように書くことができる。
ここでは [decimal モジュール](https://docs.python.org/ja/3/library/decimal.html)を利用している。このモジュールは標準ライブラリに含まれる。

```{code-cell}
from decimal import BasicContext, Decimal, getcontext, setcontext
setcontext(BasicContext)  # 丸め
getcontext().prec = 5     # 精度は10進5桁

x = Decimal(101).sqrt()
y = Decimal(100).sqrt()

print('x:', x)
print('y:', y)
print('x - y       =', x - y)
print('1 / (x + y) =', 1 / (x + y))
```

以下の例は，情報落ちによって加算の結合則が成り立たない例である。

```{code-cell}
from decimal import BasicContext, Decimal, getcontext, setcontext
setcontext(BasicContext)  # 丸め
getcontext().prec = 3     # 精度は10進3桁

x = Decimal('4.56')
y = Decimal('6.78')
z = Decimal('1.23')
print('x:', x)
print('y:', y)
print('z:', z)
print('(x + y) + z =', (x + y) + z)
print('x + (y + z) =', x + (y + z))
```

3 個の数の和ではわずかな違いに過ぎないが，多くの数の和をとる場合には大きな違いをもたらす可能性がある。

最後に，二次方程式の例を見てみよう。

```{code-cell}
from decimal import BasicContext, Decimal, getcontext, setcontext
setcontext(BasicContext)  # 丸め
getcontext().prec = 3     # 精度は10進3桁

a = Decimal(0.1)
b = Decimal(10)
c = Decimal(0.1)
D = b * b - 4 * a * c
x1 = (- b - D.sqrt()) / (2 * a)
x2 = (- b + D.sqrt()) / (2 * a)
print(f'x1: {x1:=+10.4f}')
print(f'x2: {x2:=+10.4f}')
x2 = c / a / x1
print(f'x2: {x2:=+10.4f}')
```

## 近似値と誤差

観測や数値計算によって得られた近似値 $\hat{x}$ の真の値を $x$ とする。

このとき，

$$
e(\hat{x})=\hat{x}-x
$$

を値 $\hat{x}$ の **誤差**（error） といい，

$$
|e(\hat{x})|=|\hat{x}-x|
$$

を $\hat{x}$ の **絶対誤差**（absolute error）という。

通常は真の値 $x$ を知ることはできないため，誤差や絶対誤差の正確な値を知ることもできない。しかし，絶対誤差の上限を知ることができる場合がある。
そこで，

$$
|e(\hat{x})| \leq \varepsilon(\hat{x})
$$

を満たすような $\varepsilon(\hat{x})$ を $\hat{x}$ の **誤差限界**（limit of error）という。

誤差の大きさを評価するときは，真の値の大きさも考慮する必要がある。
たとえば，靴のサイズで 1cm のズレがあったらほぼ間違いなく支障が出るが，100m 走やマラソンのコースを測るときに 1cm のズレがあってもそれほど支障はないことが多い。

そこで，真の値に対する誤差の大きさ

$$
e_{\mathrm{rel}}(\hat{x})=\frac{e(\hat{x})}{x}
$$

を **相対誤差**（relative error）とよぶ。
さらに，

$$
|e_{\mathrm{rel}}(\hat{x})| \leq \varepsilon_{\mathrm{rel}}(\hat{x})
$$

を満たすような $\varepsilon_{\mathrm{rel}}(\hat{x})$ を $\hat{x}$ の **相対誤差限界**（limit of error）という。

誤差について考えるとき，符号は重要でないことが多い。
したがって，絶対誤差 $|e(\hat{x})|$ を単に「誤差」とよんだり，
相対誤差の絶対値 $|e_{\mathrm{rel}}(\hat{x})|$ を単に「相対誤差」とよぶことも多い。

一般に，近似値 $\hat{x}$ が $\beta$ 進 $t$ 桁の浮動小数点数

$$
  (d_0 . d_1 d_2 \cdots d_{t-1})_{(\beta)} \times \beta^e
  \quad
  \text{ただし $d_0 \neq 0$}
$$

によって表されるとき，$\hat{x}$ の相対誤差限界は

$$
\varepsilon(\hat{x})_{\mathrm{rel}} = \frac{1}{2} \beta^{1-t}
$$

である。

````{prf:example}
前述の例で，10進数 $123.456$ を 2 進 5 桁の浮動小数点数に変換すると $1.1111_{(2)} \times 2^6$ となる。

このときの絶対誤差は $|e(\hat{x})|=|124-123.456|=0.544$ であり，
相対絶対誤差は $|e_{\mathrm{rel}}(\hat{x})|=|e(\hat{x})/x|=0.544/123.456=0.00440\cdots$ である。

一般に，ある数 $x$ を 2 進 5 桁の浮動小数点数に変換したときの相対誤差限界は $0.00001_{(2)}=2^{-5}$ 以下，すなわち，3.2% 以下である。これが相対誤差限界となる。
````

````{prf:example}
近似値 $\hat{x}=1.23$ について，真の値 $x$ の範囲は $1.225$ 以上 $1.235$ 未満であると考えられる（四捨五入によって丸めた場合）。したがって，誤差 $e(\hat{x})=\hat{x}-x$ は

$$
1.23 - 1.235 < e(\hat{x}) \leq 1.23 - 1.225
$$

となる。すなわち，絶対誤差 $|e(\hat{x})|$ は $0.005$ 以下である。
よって $\hat{x}$ の誤差限界は $0.005$ である。

近似値を $\hat{y}=1.2300$ と書いた場合は，誤差限界は $0.00005$ となる。
つまり，末尾の $0$ の有無は意味を持つ。

近似値 $\hat{z}=0.0123$ の誤差限界もやはり $0.00005$ であり，$\varepsilon(\hat{y})$ と $\varepsilon(\hat{z})$ は等しい。
````

````{prf:example}
$\hat{y}=1.2300$ と $\hat{z}=0.0123$ の相対誤差限界は大きく異なる。

実際，$\hat{y}$ の相対誤差は

$$
\frac{1.2300 - 1.23005}{1.23005} < e_{\mathrm{rel}}(\hat{y}) \leq \frac{1.2300 - 1.22995}{1.22995}
$$

より $\varepsilon_{\mathrm{rel}}(\hat{y}) \approx 0.00004$ である。

$\hat{z}$ の相対誤差は

$$
\frac{0.0123 - 0.01235}{0.01235} < e_{\mathrm{rel}}(\hat{z}) \leq \frac{0.0123 - 0.01225}{0.01225}
$$

より $\varepsilon_{\mathrm{rel}}(\hat{z}) \approx 0.004$ である。

以上より，$\varepsilon_{\mathrm{rel}}(\hat{y})$ と $\varepsilon_{\mathrm{rel}}(\hat{y})$ はおよそ 100 倍の違いがある。
これは，両者の有効桁数の違いによるものである。

$\hat{x}=1.23$ と $\hat{z}=0.0123$ は有効桁数が同じであるから，相対誤差限界は等しい。
````
