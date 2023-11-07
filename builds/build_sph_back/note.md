asobi=0.1
を入れた

遊びがあるとズレが上部まで広がっていくようだ．
いや，圧力勾配の計算があわないようだ．

Rungekutta を修正することにし
RK = 1
\label{SPH:gradP1}だと下に下がってきた
\label{SPH:gradP0}だとどうだ？ 綺麗に下がってきた．

密度を更新しないことにすると？
そこまで変化はない

後退オイラーのように変更してみる．


stokesで
// \label{SPH:rho_next}
double rho_next(auto p) {
   // return rho_next_(p);
   // return rho_next_(p);
   return _WATER_DENSITY_;
};
とだけ変更して計算．
