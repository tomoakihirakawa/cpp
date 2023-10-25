

# ケーブルについて

[鋼構造シリーズ１１　ケーブル・スペース構造の基礎と応用](http://library.jsce.or.jp/Image_DB/committee/steel_structure/bklist/47254.html)にわかりやすい解説があった．
この[4 ケーブル構造の動的問題](http://library.jsce.or.jp/Image_DB/committee/steel_structure/book/47254/47254-0057.pdf)によると，ケーブルの運動方程式は以下のように表される．

```math
\frac{\partial}{\partial s}\left(T\frac{\partial \boldsymbol{r}}{\partial s}\right) = \rho \frac{\partial^2 \boldsymbol{r}}{\partial t^2} + \rho g \boldsymbol{e}_z + \boldsymbol{F}
```

ここで，

* 曲線座標を$`s`$
* ケーブル上位置ベクトルを$`\boldsymbol{r}(x,t)`$
* ケーブルの張力を$`T(s,t)`$
* 単位長さあたりの質量を$`\rho`$
* 重力加速度を$`g`$
* 外力を$`\boldsymbol{F}(s,t)`$（単位長さあたりの力となるだろう）