# Vector_Soliton

## 自身の研究との関連

光の特性を決めるものとして振幅,位相,偏光があるが自分は偏光について研究している.偏光についてより深い理解をするために,授業の内容を生かして偏光に関係したソリトン,ベクトルソリトンというものについて調べる.

## ベクトルソリトン

光の分野ででてくる通常のソリトンは1つの偏光成分しか持たないがベクトルソリトンは2つの偏光成分を持つ.その特異な性質から~
ベクトルソリトンとして

- bright-bright(BB) soliton
- bright-dark(BD) soliton
- dark-dark(DD) soliton

がある.もっと言えばN個のソリトンのうちs個はdark solitonでm個はbright solitonのようなソリトンを作ることができる.

## 2-ベクトルソリトンの導出

以下に示す2次元の非線形シュレディンガー方程式を考える.

$$
i\frac{\partial q _ j}{\partial t}+\frac{\partial^2 q _ j}{\partial x^2}+2Q(q _ 1,q _ 2)q _ j=0 \quad (j=1,2) \quad (1)
$$

ただし $q _ 1,q _ 2$ はenvelop waveで

- $Q(q _ 1, q _ 2)=a\vert q _ 1 \vert ^2+c\vert q _ 2 \vert ^2+bq _ 1 q^\ast _ 2+b^\ast q^\ast _ 1 q _ 2$
  
また $a,b,c$ はそれぞれ

- $a:$ self-phase modulation(SPM)
- $b:$ four-wave mixing(FWM)
- $c:$ cross-phase modulation(XPM)

の強さを決めるパラメーターである.このパラメーターの違いによってソリトン同士が反射したり透過したりする.

ここで

$$
\begin{aligned}
q _ 1 &= \psi _ 1-b^{\ast}\psi _ 2 \newline
q _ 2 &= a \psi _ 2
\end{aligned}
$$

と表してみる.これを(1)式に代入すると

$$
\begin{aligned}
i\frac{\partial \psi _ 1}{\partial t}+\frac{\partial^2 \psi _ 1}{\partial x^2}+2a(\vert \psi _ 1 \vert ^2+ \sigma \vert \psi _ 2 \vert ^2 )\psi _ 1 &=0 \quad (2) \newline
i\frac{\partial \psi _ 2}{\partial t}+\frac{\partial^2 \psi _ 2}{\partial x^2}+2a(\vert \psi _ 1 \vert ^2+ \sigma \vert \psi _ 2 \vert ^2 )\psi _ 2 &=0 \quad(3)
\end{aligned} 
$$

を得る.さらにGを複素関数,Fを実関数とし

$$
\psi _ j = \frac{G _ j}{F}
$$

とすると(2),(3)式は

$$
\begin{aligned}
(i D _ t+D _ x^2 -\lambda)G _ j \cdot F &= 0 \newline
(D _ x^2-\lambda)F \cdot F &= 2a(\vert G _ 1 \vert ^2+\sigma \vert G _ 2 \vert ^2)
\end{aligned}
$$

とまとめられる.

## 

<!-- ## 非線形シュレディンガー方程式

非線形カー効果が発生するとき伝播する光は以下の非線形シュレディンガー方程式を満たす.

$$
i\frac{\partial \psi}{\partial z}=\frac{1}{2k _ 0}\Bigl(\frac{\partial^2 \psi}{\partial r^2} + \frac{1}{r}\frac{\partial \psi}{\partial r}\Bigr)+\kappa \vert \psi^2 \vert \psi
$$

- $k _ 0:$真空中の波数
- $\kappa:$非線形パラメーター

ただし円筒座標系で表している.

ここで

$$
\begin{aligned}
\~{r}&=\frac{r}{a _ 0} \newline
\~{z}&=\frac{z}{L _ D} \newline
L _ D&\coloneqq2k _ 0a _ 0^2
\end{aligned}
$$

とすると

$$
i\frac{\partial \psi}{\partial \~{z}}=\frac{\partial^2 \psi}{\partial \~{r}^2} + \frac{1}{\~{r}}\frac{\partial \psi}{\partial \~{r}}+\~{\kappa} \vert \psi^2 \vert \psi
$$

と書ける.改めてこれを

$$
i\frac{\partial \psi}{\partial z}=\frac{\partial^2 \psi}{\partial r} + \frac{1}{r}\frac{\partial \psi}{\partial r}+\kappa \vert \psi^2 \vert \psi
$$

と書くことにする.ここで定常解

$$
\psi=A(r)e^{i\delta z}
$$

を仮定する.代入すると

$$
\frac{d^2 A}{d r^2}+\frac{dA}{dr}+\delta A+\kappa A^3=0
$$ -->

```python
def a(b):
  return b**2
```

## 参考文献
[?]Agalar Agalarov, Vladimir Zhulego, and Telman Gadzhimuradov Phys. Rev. E 91, 042909 – Published 17 April 2015