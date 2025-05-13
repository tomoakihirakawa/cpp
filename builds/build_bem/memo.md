

三角要素の頂点情報を使って，３次元の補間をする方法．

１次多項式で次のように補間する．

$$
f(x,y,z)=ax+by+cz+d
$$

しかし，
３頂点の情報は，よりも未知変数が１つ多く，解くことができない．
そこで，以下の形で解くことにする．

$$
f(x,y)=ax+by+f_0
$$

$$
\begin{bmatrix}
a\\
b\\
f_0
\end{bmatrix}
=
\begin{bmatrix}
-(1/x_1) & 1/x_1 & 0\\
(-x_1 + x_2)/(x_1 y_2) & -(x_2/(x_1 y_2)) & 1/y_2\\
1 & 0 & 0
\end{bmatrix}
\begin{bmatrix}
f_0\\
f_1\\
f_2
\end{bmatrix}
$$

$$
\left({\begin{gathered}
{{f}_{1}}\\
{{f}_{2}}\\
{{f}_{3}}
\end{gathered}}\right)=\left({\begin{gathered}
{\left({{x}_{1},{y}_{1},{z}_{1}}\right)}\\
{\left({{x}_{2},{y}_{2},{z}_{2}}\right)}\\
{\left({{x}_{3},{y}_{3},{z}_{3}}\right)}
\end{gathered}}\right)\left({\begin{gathered}
{a}\\
{b}\\
{c}
\end{gathered}}\right)
$$


$R^{-1}$が$R^{\top}$であることと，$I^{-1}$が対角成分のみの行列であることを利用すれば，次のように書ける．

\begin{align*}
{\left({{I}_{G}}\right)}_{il}&={\left({{R}^{-{1}}}\right)}_{ij}{I}_{jk}{R}_{kl}={R}_{ji}{I}_{jj}{R}_{jl}={R}_{1i}{R}_{1l}{I}_{x}+{R}_{2i}{R}_{2l}{I}_{y}+{R}_{3i}{R}_{3l}{I}_{z}\\
{\left({{I}_{G}^{-{1}}}\right)}_{il}&={\left({{R}^{-{1}}}\right)}_{ij}{\left({{I}^{-{1}}}\right)}_{jk}{R}_{kl}={R}_{ji}{\left({{I}^{-{1}}}\right)}_{jj}{R}_{jl}=\frac{{R}_{1i}{R}_{1l}}{{I}_{x}}+\frac{{R}_{2i}{R}_{2l}}{{I}_{y}}+\frac{{R}_{3i}{R}_{3l}}{{I}_{z}}
\end{align*}


\begin{align*}
0=\sum\limits_{k_{\vartriangle}}^{}{\sum\limits_{\xi_{1}w_{1}}^{}{\sum\limits_{\xi_{0}w_{0}}^{}{\left\{{{{1}\over{\left\|{x_{k_{\vartriangle}}\left({{\rm \xi}}\right)-x_{i_{\circ}}}\right\|}}\left({\sum\limits_{j=0}^{2}{{\left({\phi_{n}}\right)}_{k_{\vartriangle}j}N_{j}\left({\xi}\right)}}\right)W_{k_{\vartriangle}}-\left\{{\alpha_{i_{\circ}}\phi_{i_{\circ}}+\left({-{{x_{k_{\vartriangle}}\left({\xi}\right)-x_{i_{\circ}}}\over{{\left\|{x_{k_{\vartriangle}}\left({\xi}\right)-x_{i_{\circ}}}\right\|}^{3}}}\cdot{{{{\partial x_{k_{\vartriangle}}}\over{\partial\xi_{0}}}\times{{\partial x_{k_{\vartriangle}}}\over{\partial\xi_{1}}}}\over{\left\|{{{\partial x_{k_{\vartriangle}}}\over{\partial\xi_{0}}}\times{{\partial x_{k_{\vartriangle}}}\over{\partial\xi_{1}}}}\right\|}}}\right)\left({\sum\limits_{j=0}^{2}{{\left({\phi}\right)}_{k_{\vartriangle}j}N_{j}\left({\xi}\right)}}\right)W_{k_{\vartriangle}}}\right\}}\right\}}}}
\end{align*}