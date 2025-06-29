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

# 数値シミュレーション

## シミュレーションとは

現実世界における現象，たとえば人口増加・感染症の拡大・株価の変動など，を分析・予測・推定することによって，人々の行動の指針や政策を作成することができる可能性がある。

ただし，人口増加のような大規模な現象や生体実験のように簡単に実行することができない現象については，計算機上で **シミュレート**（simulate）する方が望ましい。

**シミュレーション**（simulation）は必ずしも計算機上でおこなわれるものではない。たとえば，自動車や飛行機の運転を体感できる装置や，流体工学で用いられる風洞実験施設のようなものも，一般には，シミュレーションに含まれる。

本稿では計算機上でおこなわれる計算機シミュレーションについて扱う。
計算機シミュレーションの典型的な手順は以下のとおりである {cite}`hashimoto_makino_2021` 。

1. モデルの作成
1. 計算アルゴリズムを考える
1. プログラムの作成と実行
1. 結果の表示と分析

**モデル**（model）とは，一般的には，ある種の模型・見本のようなものを指す（たとえば，プラモデル，ファッションモデル，モデルケースなど）。
モデルは，対象を簡易的に表現したもので，モデルを利用する目的は，対象の本質的な性質を抽出することである。

何らかの物理的な装置を用いる物理モデル（風洞実験，フライトシミュレータ）に対して，数式や計算機によって表現されたモデルを論理モデルという。
論理モデルには，たとえば，人間の脳（というより神経細胞）を模して考えられたニューラルネットワークや，物理現象を数式によって記述する運動方程式のようなものも含まれる。

モデルを作成する際の原則は，たとえば，以下のとおりである {cite}`hashimoto_makino_2021`。

- 対象のどのような性質や現象に注目するかなど，モデル作成の目的を明確にする。
- 初めは単純なモデルを考えてから段階的に複雑なモデルを考える。
- 技術的・科学的・社会学的な法則や原理を極力利用すること。しかし，ときとして現象を直感的に数学的に表現するほうがわかりやすい場合がある。

## モンテカルロ法

特定の値（必ずしも確率的な現象でなくともよい）を求めるために確率的シミュレーションを繰返し実施する方法のことを **モンテカルロ法**（Monte Carlo method）とよぶ。
1つの値を求めるだけでなく，自然現象や社会現象など複雑・大規模な現象を限られた標本から近似して再現するような方法もモンテカルロ法に含まれる。
ある種の確率分布から疑似乱数を生成する方法（von Neumann の棄却法）や，ベイズ推定の MCMC（Markov Chain Monte Carlo; マルコフ連鎖モンテカルロ法）のように，多くの応用をもつ重要な手法のひとつである。

ここでは，円周率（$\pi=3.1415\cdots$）という有名な数学的定数を単純な確率的シミュレーションにもとづいて求める方法を紹介する。

円周率の近似値は，古代バビロニアから用いられていたとされる。
その後，実用的な範囲を超え，数学や計算機科学における興味の対象として円周率の近似値を求める試みが続けられている。
2022年現在，半年近く費やして100兆桁まで求めたというのが世界記録である。
一方，NASA（アメリカ航空宇宙局）のように高精度な軌道計算が求められる組織でも[小数点以下15桁程度しか用いていない](https://www.jpl.nasa.gov/edu/news/2016/3/16/how-many-decimals-of-pi-do-we-really-need/)と言われている。

円周率の値を，高速に，高精度に求めるためには数学的な工夫が重要である。
人間の手で円周率を求める際と同様に，級数展開などの公式を用いることになる。

一方，乱数を用いたシミュレーションによって円周率などの定数を求める簡便な方法が知られている。
そのひとつがモンテカルロ法である。

まず，一様乱数を用いて区間 $[0,1)$ の乱数を 2 つ生成し，$x$，$y$ とする。

```{code-cell}
import numpy as np
x, y = np.random.rand(2)
```

これにより，座標平面上の 4 点 $\mathrm{O}(0,0)$，$\mathrm{A}(1,0)$，$\mathrm{C}(1,1)$，$\mathrm{B}(0,1)$ を頂点とする正方形領域内の 1 点をランダムに選ぶことができる。
厳密にいえば， [`numpy.random.rand()`](https://numpy.org/doc/stable/reference/random/generated/numpy.random.rand.html) の戻り値は半開区間であり，$x=0$ や $y=0$ となることはあるが $x=1$ や $y=1$ とならないという点に注意が必要である。
ただし，この違いが今回の計算結果に影響を与えることはないため，無視してよい。

```{code-cell}
:tags: [remove-input]

import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

plt.scatter([x], [y])

ax.tick_params(labelsize=18)
ax.axis('scaled')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.grid(True)

plt.show()
```

上記の方法でランダムに選択した 1 点 $(x,y)$ が，原点 $\mathrm{O}(0,0)$ を中心とする円の内部に位置する確率を考える。
正方形の面積が $1$，円のうち正方形と重なる扇型の面積が $\displaystyle\frac{\pi}{4}$ であるから，求める確率は

$$
\frac{\pi}{4}
$$

である。

したがって，ランダムな点を大量に生成すれば，大数の法則により，ランダムに生成した点が円内に収まる割合は $\displaystyle\frac{\pi}{4}$ に近づくことが分かる。

以下は Python プログラムによるシミュレーションの例である。

```{code-cell}
import numpy as np
import matplotlib.pyplot as plt

# 乱数生成器の初期化（再現性のために種を指定する）
rng = np.random.default_rng(12345678)

# 標本数
N = 50

X = rng.random(N)
Y = rng.random(N)

# 描画
fig = plt.figure()
ax = fig.add_subplot(111)

I = X ** 2 + Y ** 2 < 1
plt.scatter(X[I], Y[I])
plt.scatter(X[~I], Y[~I])

ax.tick_params(labelsize=18)
ax.axis('scaled')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.grid(True)

plt.show()

pi = sum(I) / len(X) * 4
print(pi)
```

標本数を増やすと以下のようになる。

```{code-cell}
:tags: [remove-input]

import numpy as np
import matplotlib.pyplot as plt

# 乱数生成器の初期化（再現性のために種を指定する）
rng = np.random.default_rng(12345678)

# 標本数
N = 10000

X = rng.random(N)
Y = rng.random(N)

# 描画
fig = plt.figure()
ax = fig.add_subplot(111)

I = X ** 2 + Y ** 2 < 1
plt.scatter(X[I], Y[I], s=1)
plt.scatter(X[~I], Y[~I], s=1)

ax.tick_params(labelsize=18)
ax.axis('scaled')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

plt.show()

pi = sum(I) / len(X) * 4
print(pi)
```

## 複雑系

数値計算によるシミュレーションにおいては数値誤差・計算誤差が問題となる場合がある。
運が良ければ，サンプル数や計算回数を増やすことによって誤差が急速に小さくなる場合がある。
というより，数値計算においてはそのような性質を満たすようにモデルやアルゴリズムを構築することが多い。

これに対して，複雑系科学やカオス理論において扱われる **カオス**（chaos）という現象は，強い**初期値鋭敏性**（sensitive dependence on the initial conditions）を持っている。
カオスとは本来「混沌」を指す言葉であり，決定性モデルにもとづくシミュレーションにおいてカオス現象が現れることは直感的には考えにくいが，この初期値鋭敏性により予測できない「混沌」が観測されることからこの名が利用される。

カオス理論の発端になったローレンツ [^Lorenz] は，以下のモデルを用いて大気の対流に関するシミュレーションをおこない，その結果，ローレンツ・アトラクタとよばれるパターンを発見した。

[^Lorenz]: Edward Norton Lorenz (1917-2008) はマサチューセッツ工科大学の気象学者。

$$
\frac{dx}{dt} &= s(y-x)
\\
\frac{dy}{dt} &= rx-y-xz
\\
\frac{dz}{dt} &= xy - bz
$$

以下に，ローレンツ・アトラクタを図示するための Python プログラムの例を示す。

```{code-cell}
def lorenz(x, y, z, s=10, r=28, b=8/3):
    x_dot = s * (y - x)
    y_dot = r * x - y - x * z
    z_dot = x * y - b * z
    return x_dot, y_dot, z_dot

dt = 0.01
xyz_ini = (0.0, 1.0, 1.05)
num_steps = 10000

xs = np.empty(num_steps)
ys = np.empty(num_steps)
zs = np.empty(num_steps)
xs[0], ys[0], zs[0] = xyz_ini

for i in range(1, num_steps):
    x_dot, y_dot, z_dot = lorenz(xs[i-1], ys[i-1], zs[i-1])
    xs[i] = xs[i - 1] + x_dot * dt
    ys[i] = ys[i - 1] + y_dot * dt
    zs[i] = zs[i - 1] + z_dot * dt

fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": "3d"})
ax.plot(xs, ys, zs, lw=0.5)
plt.show()
```

ただし，このプログラムはあくまで一例であり，ローレンツの計算によるローレンツ・アトラクタと必ずしも一致しない。
また，微分方程式の解法は，単純なオイラー法によるプログラムである。

- https://matplotlib.org/stable/gallery/mplot3d/lorenz_attractor.html
- https://en.wikipedia.org/wiki/Lorenz_system
