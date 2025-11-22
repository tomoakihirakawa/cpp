ケース：case5
63291節点
方法:Fourier M2L double

初回（1step目のRKの1回目）は収束が早いので，2回目以降の計算を比較すること．

------- LOG -------


net->velocityTranslational() = {updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0,0,0}
0.0047955,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
use tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.998468,100.179} s
RK_step = 4/4, RK_time = 0.01, simulation_time = 0
setNeumannVelocity: tank
move_name_velocity = fixed
setting acceleration
setNeumannVelocity: wavemaker
move_name_velocity = piston
(RigidBodyObject) velocity is explicityly given as piston
T6d velocity(const std::string &name, const std::vector<std::string> strings, double t)piston
A = 0.136, w = 3.14159, k = 1.10953, h = 1.36, {T, L} = {2, 5.66291}
T6d velocity(const std::string &name, const std::vector<std::string> strings, double t)piston
A = 0.136, w = 3.14159, k = 1.10953, h = 1.36, {T, L} = {2, 5.66291}
T6d velocity(const std::string &name, const std::vector<std::string> strings, double t)piston
A = 0.136, w = 3.14159, k = 1.10953, h = 1.36, {T, L} = {2, 5.66291}
net->velocity = {0.00958981,0,0,0,0,0}
setting acceleration
指定がないので加速度はdefault_acceleration
setNeumannVelocity
Elapsed time: {5e-05,100.179} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 71076
Creating sources on surfaces...
バケットの作成．極の追加
added :126578
バケットの作成．極の追加, Elapsed time : {6.48731,6.48731}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 1136 deepest level buckets.
After shrink: 1136 deepest level buckets.
After grow: 1136 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.731082,7.2184}
initializeFMM
SimpleM2L
[FMM:init] Step0: update poles start, #poles=126578
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=1505, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=1136, total leaf sources=126578
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.034603,0.034603}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.001426,0.001426}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000149,0.001575}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.005675,0.00725}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.01045,0.0177}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.021757,0.039457}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {0.719478,0.758935}
setM2L, level=6, buckets_at_a_level.size()=992, Elapsed time : {0.76163,1.52056}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000207,1.52077}
setM2L, Elapsed time : {1.52078,1.52078}
[FMM:init] Step3: setM2L done
level=0
mean A->buckets_for_M2L=0
level=1
mean A->buckets_for_M2L=0
level=2
mean A->buckets_for_M2L=6
level=3
mean A->buckets_for_M2L=9
level=4
mean A->buckets_for_M2L=10.5
level=5
mean A->buckets_for_M2L=58.5918
level=6
mean A->buckets_for_M2L=28.1532
level=7
[FMM:init] Debug: coupling stats printed
[FMM:init] Step4: setL2L start
[FMM:init] Step4: setL2L done
各レベルの各セルのL2Pの相手を保存する (#targets=63291)
setL2P, Elapsed time : {0.016482,0.016482} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=63291)
setDirectIntegration, Elapsed time : {1.58714,1.60362} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=1118.71 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.001774,0.001774}
reset Moments, Elapsed time : {0.001097,0.002871}
increment_M_reuse, Elapsed time : {0.002331,0.005202}
M2M, Elapsed time : {0.009319,0.014521}
M2L, Elapsed time : {0.047604,0.062125}
L2L, Elapsed time : {0.010022,0.072147}
updateFMM, Elapsed time : {7e-06,0.072154}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00134646s 1.51528 %
reset Moments               0.00112098s 1.26153 %
increment_M_reuse           0.0021662s 2.4378 %
M2M                         0.0092392s 10.3976 %
M2L                         0.0473607s 53.2987 %
L2L                         0.00995342s 11.2014 %
others (direct integration) 0.017672s 19.8877 %
total                       0.088859 s
--------------------------------------------------
Elapsed time for GMRES: {4.68149,4.68149} [s] for size = 50
GMRES expected error = 0.00142949 True error = 0.00142949
destructing gmres
Restart count :1
----- Average time for 50 iterations -----
update poles                0.00133732s 1.51437 %
reset Moments               0.0011237s 1.27246 %
increment_M_reuse           0.00219958s 2.49078 %
M2M                         0.00924222s 10.4658 %
M2L                         0.0469651s 53.1828 %
L2L                         0.00995702s 11.2752 %
others (direct integration) 0.017484s 19.7986 %
total                       0.0883089 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00124314s 1.40649 %
reset Moments               0.00112744s 1.27559 %
increment_M_reuse           0.00220296s 2.49244 %
M2M                         0.00923152s 10.4446 %
M2L                         0.0472985s 53.5137 %
L2L                         0.00995994s 11.2687 %
others (direct integration) 0.0173223s 19.5985 %
total                       0.0883858 s
--------------------------------------------------
Elapsed time for GMRES: {9.45908,9.45908} [s] for size = 100
GMRES expected error = 1.48758e-06 True error = 1.48758e-06
destructing gmres
Restart count :2
----- Average time for 50 iterations -----
update poles                0.00134452s 1.50744 %
reset Moments               0.00113224s 1.26944 %
increment_M_reuse           0.00222016s 2.48919 %
M2M                         0.0092399s 10.3596 %
M2L                         0.0472191s 52.9409 %
L2L                         0.00998894s 11.1994 %
others (direct integration) 0.0180472s 20.2341 %
total                       0.089192 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00129104s 1.42806 %
reset Moments               0.00114236s 1.2636 %
increment_M_reuse           0.0022649s 2.50528 %
M2M                         0.00926226s 10.2453 %
M2L                         0.0482889s 53.414 %
L2L                         0.0100378s 11.1032 %
others (direct integration) 0.0181177s 20.0406 %
total                       0.090405 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00116498s 1.3248 %
reset Moments               0.0011255s 1.2799 %
increment_M_reuse           0.00219974s 2.50151 %
M2M                         0.00925522s 10.5249 %
M2L                         0.0466455s 53.0446 %
L2L                         0.00996042s 11.3268 %
others (direct integration) 0.0175851s 19.9975 %
total                       0.0879365 s
--------------------------------------------------
Elapsed time for GMRES: {14.609,14.609} [s] for size = 150
GMRES expected error = 3.41927e-11 True error = 3.41927e-11
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {0.00142949,1.48758e-06,3.41927e-11}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: {41.1942,41.1942} s