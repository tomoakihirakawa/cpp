
参考：

https://www.youtube.com/watch?v=E9_kyXjtRHc&list=PLnJ8lIgfDbkp5DtCPtP2rcqEEUJk-PM8N


保存量$\boldsymbol \Lambda$，
保存量$\boldsymbol \Lambda$の流束$\boldsymbol \Pi$，
湧き出し$\boldsymbol \Gamma$，を用いて，
流れの基礎方程式となる保存則は，

```math
\int_V \frac{\partial {\boldsymbol \Lambda}}{\partial t} dV
+ \int_S {\boldsymbol \Pi} \cdot {\boldsymbol n} dS
= \int_V {\boldsymbol \Gamma} dS
```

```math
{\boldsymbol \Lambda} = \left( 
\begin{array}{c}
\rho\\
\rho {\boldsymbol u}\\
\rho E
\end{array} \right)
,
{\boldsymbol \Pi} = \left(
\begin{array}{c}
\rho {\boldsymbol u}\\
\rho {\boldsymbol u}{\boldsymbol u} - {\boldsymbol T}\\
\rho {\boldsymbol u} E - {\boldsymbol T} \cdot {\boldsymbol u} + {\boldsymbol q}
\end{array} \right)
,
{\boldsymbol \Gamma} = \left(
\begin{array}{c}
0\\
\rho {\boldsymbol f}\\
\rho {\boldsymbol f} \cdot {\boldsymbol u}
\end{array} \right)
```

ここで，${\boldsymbol u}{\boldsymbol u}$は，${\boldsymbol u}\otimes{\boldsymbol u}$であり，または$u_i u_j$である．${\boldsymbol T}$は，応力テンソルであり，非圧縮性流体の場合，
${\boldsymbol T} = -p{\boldsymbol I} + \mu \left( \nabla {\boldsymbol u} + (\nabla {\boldsymbol u})^\top \right)$
である．
$\nabla {\boldsymbol u}$は，速度勾配テンソルであり，
$\nabla\otimes {\boldsymbol u} = \nabla_i {\boldsymbol u}_j$や
$\frac{\partial u_i}{\partial x_j}$である．

## 確認

運動量保存の式が，$\rho {\boldsymbol u}$, $\rho {\boldsymbol u} {\boldsymbol u} + p{\boldsymbol I} - \mu (\nabla {\boldsymbol u} + (\nabla {\boldsymbol u})^\top)$, $\rho {\boldsymbol f}$であることをNS方程式を変形して確認する．

```math
\begin{align}
\frac{D {\rho \boldsymbol u}}{D t} = \frac{\partial {\rho \boldsymbol u}}{\partial t} + {\rho \boldsymbol u} \cdot \nabla {\boldsymbol u} = -\nabla p + \nabla \cdot (\mu \nabla {\boldsymbol u}) + \rho {\boldsymbol g}\\
\frac{\partial \rho u_i}{\partial t} + \rho u_j \frac{\partial u_i}{\partial x_j} = -\frac{\partial p}{\partial x_i} + \frac{\partial}{\partial x_j} \left( \mu \frac{\partial u_i}{\partial x_j} \right) + \rho g_i\\
\frac{\partial \rho u_i}{\partial t} + \frac{\partial (\rho u_i) u_j }{\partial x_j} - u_i \frac{\partial \rho u_j }{\partial x_j} = -\frac{\partial p}{\partial x_i} + \frac{\partial}{\partial x_j} \left( \mu \frac{\partial u_i}{\partial x_j} \right) + \rho g_i\\
\int_V\frac{\partial \rho u_i}{\partial t} dV 
+ \int_S{\rho u_i u_j }n_j dS 
- \int_V u_i \frac{\partial \rho u_j }{\partial x_j} dV 
= -\int_V \frac{\partial p}{\partial x_i} dV
+ \int_V \frac{\partial}{\partial x_j} \left( \mu \frac{\partial u_i}{\partial x_j} \right) dV
+ \int_V \rho g_i dV\\
\int_V\frac{\partial \rho u_i}{\partial t} dV + \int_S{\rho u_i u_j }n_j dS 
- \int_V \left( u_i \frac{\partial \rho u_j }{\partial x_j} dV -  \frac{\partial p}{\partial x_i} dV
+ \frac{\partial}{\partial x_j} \left( \mu \frac{\partial u_i}{\partial x_j} \right) 
\right)dV = \int_V \rho g_i dV
\end{align}
```



```math
\begin{align}
\frac{D {\boldsymbol u}}{D t} 
=
\rho \frac{\partial {\boldsymbol u}}{\partial t} 
+ \rho {\boldsymbol u} \cdot \nabla {\boldsymbol u} = -\nabla p + \nabla \cdot (\mu \nabla {\boldsymbol u}) + \rho {\boldsymbol f}
\end{align}
```

