ケース：case4
26158節点
方法:SimpleM2L

初回（1step目のRKの1回目）は収束が早いので，2回目以降の計算を比較すること．

------------ LOG -----------


Last login: Thu Sep 18 13:19:20 on ttys040
/Users/tomoaki/.zshrc is loaded
build_bem(main)$ ./SimpleM2L /Users/tomoaki/Ruehl20160d075_SimpleM2L
input_directory : "/Users/tomoaki/Ruehl20160d075_SimpleM2L"
ALE: {pseudo_quad}
ALEPERIOD: {1}
GRAVITY: {9.81}
WATER_DENSITY: {1000.0}
element: {linear}
end_time: {10.0}
end_time_step: {100000}
input_files: {tank.json,wavemaker.json,water.json,wg3.json,wg6.json}
max_dt: {0.1}
output_directory: {/Users/tomoaki/Ruehl20160d075_SimpleM2L/output}
LINEAR_ELEMENT
ALE_ON_PSEUDO_QUADRATIC_ELEMENT

WATER_DENSITY: 1000
GRAVITY: 9.81
"/Users/tomoaki/Ruehl20160d075_SimpleM2L/tank.json"
isFixed: {true}
name: {tank}
objfile: {/Users/tomoaki/Dropbox/code//cpp/obj/Ruehl2016//tank.obj}
type: {RigidBody}
velocity: {fixed}
/Users/tomoaki/Dropbox/code//cpp/obj/Ruehl2016//tank.obj is opened
/Users/tomoaki/Dropbox/code//cpp/obj/Ruehl2016//tank.obj is closed
Load3DFile : lines of file are splited
Load3DFile : file is loaded
Load3DFile : complex is generated
Load3DFile : max and min are calculated
/* -------------------------------------------------------------------------- */
        this->getName(): tank, 0x140008000
           Lines.size(): 76
          Points.size(): 28
           Faces.size(): 52
           this->bounds: {{-7,22},{-3,3},{-1,3.5}}
                              0     1     2     3     4     5     6     7     8     9    10    11    12    13    14
       Lines of points :      0     0     0     0    12     5     4     3     3     0     1     0     0     0     0
   Connection of faces :      0     0    74 全ての線が2つの面と接続しているため，格子は閉じた面を形成している
/* -------------------------------------------------------------------------- */
setOutputInfo
setTypes
RigidBody
set velocity
/Users/tomoaki/Ruehl20160d075_SimpleM2L/output/tank_init.vtu  VV_points.size() : 52 |||||||||
"/Users/tomoaki/Ruehl20160d075_SimpleM2L/wavemaker.json"
isFixed: {false}
name: {wavemaker}
objfile: {/Users/tomoaki/Dropbox/code//cpp/obj/Ruehl2016//wavemaker.obj}
type: {RigidBody}
velocity: {piston,0,0.136,2,1.36,1,0,0}
/Users/tomoaki/Dropbox/code//cpp/obj/Ruehl2016//wavemaker.obj is opened
/Users/tomoaki/Dropbox/code//cpp/obj/Ruehl2016//wavemaker.obj is closed
Load3DFile : lines of file are splited
Load3DFile : file is loaded
Load3DFile : complex is generated
Load3DFile : max and min are calculated
/* -------------------------------------------------------------------------- */
        this->getName(): wavemaker, 0x128008000
           Lines.size(): 744
          Points.size(): 250
           Faces.size(): 496
           this->bounds: {{-5,0},{-1,1},{0,3}}
                              0     1     2     3     4     5     6     7     8     9    10    11    12    13    14
       Lines of points :      0     0     0     0     4     4   242     0     0     0     0     0     0     0     0
   Connection of faces :      0     0   744 全ての線が2つの面と接続しているため，格子は閉じた面を形成している
/* -------------------------------------------------------------------------- */
setOutputInfo
setTypes
RigidBody
set velocity
/Users/tomoaki/Ruehl20160d075_SimpleM2L/output/wavemaker_init.vtu  VV_points.size() : 496 |||||||||
"/Users/tomoaki/Ruehl20160d075_SimpleM2L/water.json"
name: {water}
objfile: {/Users/tomoaki/Dropbox/code//cpp/obj/Ruehl2016//water0d075.obj}
reverseNormal: {false}
triangles: {12}
type: {Fluid}
vertices: {8}
/Users/tomoaki/Dropbox/code//cpp/obj/Ruehl2016//water0d075.obj is opened
/Users/tomoaki/Dropbox/code//cpp/obj/Ruehl2016//water0d075.obj is closed
Load3DFile : lines of file are splited
Load3DFile : file is loaded
Load3DFile : complex is generated
Load3DFile : max and min are calculated
/* -------------------------------------------------------------------------- */
        this->getName(): water, 0x14ffe0000
           Lines.size(): 78468
          Points.size(): 26158
           Faces.size(): 52312
           this->bounds: {{0,20},{-1,1},{0,1.36}}
                              0     1     2     3     4     5     6     7     8     9    10    11    12    13    14
       Lines of points :      0     0     0     0     0  2850 20523  2732    53     0     0     0     0     0     0
   Connection of faces :      0     0 78468 全ての線が2つの面と接続しているため，格子は閉じた面を形成している
/* -------------------------------------------------------------------------- */
setOutputInfo
setTypes
Fluid
set velocity
/Users/tomoaki/Ruehl20160d075_SimpleM2L/output/water_init.vtu  VV_points.size() : 52312 |||||||||
"/Users/tomoaki/Ruehl20160d075_SimpleM2L/wg3.json"
name: {wg3}
position: {9.48,0.0,1.86,9.48,0.0,0.8600000000000001}
type: {wavegauge}
type = wavegauge
skipped
"/Users/tomoaki/Ruehl20160d075_SimpleM2L/wg6.json"
name: {wg6}
position: {16.778,-0.031,2.0,16.778,-0.031,0.5000000000000001}
type: {wavegauge}
type = wavegauge
skipped
setting done
flipIf
flipIf
flipIf
net.getPoints() = 26158
Total : 26158
Total variables: 0
Total case double-node : 26158
node reduction : 0
CORNER : 0
Total CORNER faces : 0
Neumann : 0
Dirichlet : 0
===========================================================================
       dt :0.01
time_step :0
real time :0
---------------------------------------------------------------------------
RK_step = 1/4, RK_time = 0, simulation_time = 0
water makeBucketPoints(tank makeBucketPoints(2.01457)wavemaker makeBucketPoints(
0.616441)
2.99541)
BucketPoints.data1D.size() = 28
      bounds : {{-14.25,29.25},{-4.5,4.5},{-2.125,4.625}}
          dL : 2.99541
bucket sizes : {15, 4, 3}
tank makeBucketPoints done
tank makeBucketFaces(2.99541)
BucketPoints.data1D.size() = 250
      bounds : {{-6.25,1.25},{-1.5,1.5},{-0.75,3.75}}
          dL : 0.616441
bucket sizes : {13, 5, 8}
wavemaker makeBucketPoints done
wavemaker makeBucketFaces(0.616441)
      bounds : {{-14.25,29.25},{-4.5,4.5},{-2.125,4.625}}
          dL : 2.99541
bucket sizes : {15, 4, 3}
      bounds : {{-14.25,29.25},{-4.5,4.5},{-2.125,4.625}}
          dL : 2.99541
bucket sizes : {15, 4, 3}
tank makeBucketFaces done
tank makeBucketTetras(2.99541)
      bounds : {{-14.25,29.25},{-4.5,4.5},{-2.125,4.625}}
          dL : 2.99541
bucket sizes : {15, 4, 3}
tank makeBucketTetras done
      bounds : {{-6.25,1.25},{-1.5,1.5},{-0.75,3.75}}
          dL : 0.616441
bucket sizes : {13, 5, 8}
      bounds : {{-6.25,1.25},{-1.5,1.5},{-0.75,3.75}}
          dL : 0.616441
bucket sizes : {13, 5, 8}
wavemaker makeBucketFaces done
wavemaker makeBucketTetras(0.616441)
      bounds : {{-6.25,1.25},{-1.5,1.5},{-0.75,3.75}}
          dL : 0.616441
bucket sizes : {13, 5, 8}
wavemaker makeBucketTetras done
BucketPoints.data1D.size() = 26158
      bounds : {{-5,25},{-1.5,1.5},{-0.34,1.7}}
          dL : 2.01457
bucket sizes : {15, 2, 2}
water makeBucketPoints done
water makeBucketFaces(2.01457)
      bounds : {{-5,25},{-1.5,1.5},{-0.34,1.7}}
          dL : 2.01457
bucket sizes : {15, 2, 2}
      bounds : {{-5,25},{-1.5,1.5},{-0.34,1.7}}
          dL : 2.01457
bucket sizes : {15, 2, 2}
water makeBucketFaces done
water makeBucketTetras(2.01457)
      bounds : {{-5,25},{-1.5,1.5},{-0.34,1.7}}
          dL : 2.01457
bucket sizes : {15, 2, 2}
water makeBucketTetras done
makeBuckets
Elapsed time: {0.998158,0.998158} s
waterの境界条件を決定 setBoundaryTypes
water setContactFaces()
step2 面の境界条件を判定
step3 線の境界条件を決定
step4 点の境界条件を決定
setBoundaryTypes終了
setBoundaryTypes
Elapsed time: {0.802935,1.80109} s
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
net->velocity = {0,0,0,0,0,0}
setting acceleration
指定がないので加速度はdefault_acceleration
setNeumannVelocity
Elapsed time: {5.7e-05,1.80115} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 31039
Setting up buckets...
Buckets set up complete.
Creating sources on surfaces...
バケットの作成．極の追加
added :52312
バケットの作成．極の追加, Elapsed time : {0.013599,0.013599}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 1 deepest level buckets.
After shrink: 1 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {1.63053,1.64413}
initializeFMM
SimpleM2L
[FMM:init] Step0: update poles start, #poles=52312
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=52312
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.010132,0.010132}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.001389,0.001389}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.00019,0.001579}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.006381,0.00796}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.013588,0.021548}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.02635,0.047898}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {0.788377,0.836275}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000206,0.836481}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000208,0.836689}
setM2L, Elapsed time : {0.836698,0.836698}
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
level=7
[FMM:init] Debug: coupling stats printed
[FMM:init] Step4: setL2L start
[FMM:init] Step4: setL2L done
各レベルの各セルのL2Pの相手を保存する (#targets=26158)
setL2P, Elapsed time : {0.007941,0.007941} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=26158)
setDirectIntegration, Elapsed time : {0.61172,0.619661} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=997.913 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000943,0.000943}
reset Moments, Elapsed time : {0.000403,0.001346}
increment_M_reuse, Elapsed time : {0.002477,0.003823}
M2M, Elapsed time : {0.003688,0.007511}
M2L, Elapsed time : {0.067277,0.074788}
L2L, Elapsed time : {0.005447,0.080235}
updateFMM, Elapsed time : {4e-06,0.080239}
integration (direct integration is dominant), Elapsed time : {0.043119,0.043119}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00057246s 0.75731 %
reset Moments               0.00035656s 0.471695 %
increment_M_reuse           0.00099966s 1.32246 %
M2M                         0.00426104s 5.63695 %
M2L                         0.0563033s 74.4839 %
L2L                         0.00599498s 7.93079 %
others (direct integration) 0.0071032s 9.39686 %
total                       0.0755912 s
--------------------------------------------------
Elapsed time for GMRES: {3.96285,3.96285} [s] for size = 50
GMRES expected error = 0 True error = 0
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {0}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 7.7998 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {7.79982,9.60097} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.022776,9.62375} s
name = tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 0, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 1, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 2, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 3, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
BVP.solveForPhiPhin_t-> {Φt,Φtn}とnet->accelerationが決まる
Elapsed time: {0.651036,10.2748} s
name = tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.382412,10.6572} s
RK_step = 2/4, RK_time = 0.005, simulation_time = 0
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
net->velocity = {0.0047955,0,0,0,0,0}
setting acceleration
指定がないので加速度はdefault_acceleration
setNeumannVelocity
Elapsed time: {4.6e-05,10.6572} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 31039
Creating sources on surfaces...
バケットの作成．極の追加
added :52312
バケットの作成．極の追加, Elapsed time : {0.622361,0.622361}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 392 deepest level buckets.
After shrink: 392 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.122911,0.745272}
initializeFMM
SimpleM2L
[FMM:init] Step0: update poles start, #poles=52312
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=52312
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.012257,0.012257}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.000704,0.000704}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000183,0.000887}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.005581,0.006468}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.011032,0.0175}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.02149,0.03899}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {0.740833,0.779823}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000196,0.780019}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000153,0.780172}
setM2L, Elapsed time : {0.780183,0.780183}
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
level=7
[FMM:init] Debug: coupling stats printed
[FMM:init] Step4: setL2L start
[FMM:init] Step4: setL2L done
各レベルの各セルのL2Pの相手を保存する (#targets=26158)
setL2P, Elapsed time : {0.00718,0.00718} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=26158)
setDirectIntegration, Elapsed time : {0.574994,0.582174} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=997.913 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000943,0.000943}
reset Moments, Elapsed time : {0.000483,0.001426}
increment_M_reuse, Elapsed time : {0.001045,0.002471}
M2M, Elapsed time : {0.003384,0.005855}
M2L, Elapsed time : {0.023186,0.029041}
L2L, Elapsed time : {0.00381,0.032851}
updateFMM, Elapsed time : {1.2e-05,0.032863}
integration (direct integration is dominant), Elapsed time : {0.006857,0.006857}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00058018s 1.50285 %
reset Moments               0.00035748s 0.925989 %
increment_M_reuse           0.00100504s 2.60338 %
M2M                         0.00316226s 8.19128 %
M2L                         0.0233959s 60.603 %
L2L                         0.00340214s 8.81265 %
others (direct integration) 0.0067022s 17.3609 %
total                       0.0386052 s
--------------------------------------------------
Elapsed time for GMRES: {2.04452,2.04452} [s] for size = 50
GMRES expected error = 0.000623408 True error = 0.000623408
destructing gmres
Restart count :1
----- Average time for 50 iterations -----
update poles                0.00058946s 1.43344 %
reset Moments               0.00037822s 0.919751 %
increment_M_reuse           0.00104906s 2.55109 %
M2M                         0.00315862s 7.68109 %
M2L                         0.0254409s 61.8669 %
L2L                         0.00347552s 8.45172 %
others (direct integration) 0.00703024s 17.096 %
total                       0.041122 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00059034s 1.51722 %
reset Moments               0.00038398s 0.986861 %
increment_M_reuse           0.00101954s 2.6203 %
M2M                         0.00319158s 8.20263 %
M2L                         0.0236337s 60.7406 %
L2L                         0.00344456s 8.85281 %
others (direct integration) 0.00664552s 17.0795 %
total                       0.0389092 s
--------------------------------------------------
Elapsed time for GMRES: {4.2809,4.2809} [s] for size = 100
GMRES expected error = 1.58529e-07 True error = 1.58529e-07
destructing gmres
Restart count :2
----- Average time for 50 iterations -----
update poles                0.00059106s 1.50672 %
reset Moments               0.00037346s 0.95202 %
increment_M_reuse           0.00101488s 2.58712 %
M2M                         0.00316758s 8.07476 %
M2L                         0.0239476s 61.0469 %
L2L                         0.00346368s 8.82957 %
others (direct integration) 0.00666992s 17.0029 %
total                       0.0392282 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00059152s 1.49222 %
reset Moments               0.0003836s 0.967706 %
increment_M_reuse           0.00104098s 2.62608 %
M2M                         0.0031758s 8.01158 %
M2L                         0.0241608s 60.9503 %
L2L                         0.00345384s 8.71299 %
others (direct integration) 0.0068336s 17.2391 %
total                       0.0396401 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00059024s 1.48289 %
reset Moments               0.0003757s 0.94389 %
increment_M_reuse           0.0010554s 2.65153 %
M2M                         0.00315938s 7.93747 %
M2L                         0.0242733s 60.9831 %
L2L                         0.00344354s 8.65138 %
others (direct integration) 0.00690578s 17.3497 %
total                       0.0398034 s
--------------------------------------------------
Elapsed time for GMRES: {6.48865,6.48865} [s] for size = 150
GMRES expected error = 1.07025e-13 True error = 1.07337e-13
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {0.000623408,1.58529e-07,1.07025e-13}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 15.6988 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {15.6988,26.356} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.022786,26.3788} s
name = tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0047955,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 0, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 1, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 2, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 3, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
BVP.solveForPhiPhin_t-> {Φt,Φtn}とnet->accelerationが決まる
Elapsed time: {0.805291,27.1841} s
name = tank
net->velocityTranslational() = {updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0047955,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.393893,27.578} s
RK_step = 3/4, RK_time = 0.005, simulation_time = 0
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
net->velocity = {0.0047955,0,0,0,0,0}
setting acceleration
指定がないので加速度はdefault_acceleration
setNeumannVelocity
Elapsed time: {5e-05,27.578} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 31039
Creating sources on surfaces...
バケットの作成．極の追加
added :52312
バケットの作成．極の追加, Elapsed time : {0.637994,0.637994}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 392 deepest level buckets.
After shrink: 392 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.12171,0.759704}
initializeFMM
SimpleM2L
[FMM:init] Step0: update poles start, #poles=52312
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=52312
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.012593,0.012593}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.000725,0.000725}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000155,0.00088}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.005622,0.006502}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.010664,0.017166}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.021981,0.039147}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {0.752964,0.792111}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000187,0.792298}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000152,0.79245}
setM2L, Elapsed time : {0.792458,0.792458}
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
level=7
[FMM:init] Debug: coupling stats printed
[FMM:init] Step4: setL2L start
[FMM:init] Step4: setL2L done
各レベルの各セルのL2Pの相手を保存する (#targets=26158)
setL2P, Elapsed time : {0.006865,0.006865} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=26158)
setDirectIntegration, Elapsed time : {0.561421,0.568286} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=997.913 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000916,0.000916}
reset Moments, Elapsed time : {0.000395,0.001311}
increment_M_reuse, Elapsed time : {0.001016,0.002327}
M2M, Elapsed time : {0.003446,0.005773}
M2L, Elapsed time : {0.023193,0.028966}
L2L, Elapsed time : {0.003803,0.032769}
updateFMM, Elapsed time : {4e-06,0.032773}
integration (direct integration is dominant), Elapsed time : {0.006542,0.006542}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00057834s 1.45613 %
reset Moments               0.00036868s 0.92825 %
increment_M_reuse           0.00103442s 2.60443 %
M2M                         0.00317476s 7.9933 %
M2L                         0.0241952s 60.918 %
L2L                         0.00345374s 8.69571 %
others (direct integration) 0.00691256s 17.4042 %
total                       0.0397177 s
--------------------------------------------------
Elapsed time for GMRES: {2.09582,2.09582} [s] for size = 50
GMRES expected error = 5.91976e-09 True error = 5.91976e-09
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {5.91976e-09}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 4.96171 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {4.96176,32.5398} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.023158,32.563} s
name = tank
net->velocityTranslational() = {updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0047955,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 0, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 1, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 2, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 3, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
BVP.solveForPhiPhin_t-> {Φt,Φtn}とnet->accelerationが決まる
Elapsed time: {0.790058,33.353} s
name = updating wavemaker's (RigidBodyObject) velocitytank
name = wavemaker
net->velocityTranslational() = {0.0047955,0,0}

net->velocityTranslational() = {use wavemaker's (RigidBodyObject) predetermiend velocity
0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.409702,33.7627} s
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
Elapsed time: {5.8e-05,33.7628} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 31039
Creating sources on surfaces...
バケットの作成．極の追加
added :52312
バケットの作成．極の追加, Elapsed time : {0.589473,0.589473}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 392 deepest level buckets.
After shrink: 392 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.117261,0.706734}
initializeFMM
SimpleM2L
[FMM:init] Step0: update poles start, #poles=52312
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=52312
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.011753,0.011753}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.000653,0.000653}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.00014,0.000793}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.005591,0.006384}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.010671,0.017055}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.024018,0.041073}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {0.727032,0.768105}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000218,0.768323}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000182,0.768505}
setM2L, Elapsed time : {0.768515,0.768515}
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
level=7
[FMM:init] Debug: coupling stats printed
[FMM:init] Step4: setL2L start
[FMM:init] Step4: setL2L done
各レベルの各セルのL2Pの相手を保存する (#targets=26158)
setL2P, Elapsed time : {0.007005,0.007005} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=26158)
setDirectIntegration, Elapsed time : {0.565666,0.572671} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=997.913 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000856,0.000856}
reset Moments, Elapsed time : {0.000379,0.001235}
increment_M_reuse, Elapsed time : {0.000989,0.002224}
M2M, Elapsed time : {0.003147,0.005371}
M2L, Elapsed time : {0.023247,0.028618}
L2L, Elapsed time : {0.00345,0.032068}
updateFMM, Elapsed time : {1.1e-05,0.032079}
integration (direct integration is dominant), Elapsed time : {0.00677,0.00677}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00058724s 1.49692 %
reset Moments               0.0003619s 0.922509 %
increment_M_reuse           0.00102224s 2.60576 %
M2M                         0.0031377s 7.99822 %
M2L                         0.0240721s 61.3614 %
L2L                         0.00346082s 8.82188 %
others (direct integration) 0.006588s 16.7933 %
total                       0.03923 s
--------------------------------------------------
Elapsed time for GMRES: {2.07163,2.07163} [s] for size = 50
GMRES expected error = 0.000623266 True error = 0.000623266
destructing gmres
Restart count :1
----- Average time for 50 iterations -----
update poles                0.0006029s 1.50743 %
reset Moments               0.00038284s 0.957215 %
increment_M_reuse           0.0010831s 2.70807 %
M2M                         0.00319968s 8.00016 %
M2L                         0.024341s 60.8599 %
L2L                         0.00350896s 8.77345 %
others (direct integration) 0.00687668s 17.1938 %
total                       0.0399952 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00060316s 1.51477 %
reset Moments               0.00038688s 0.971607 %
increment_M_reuse           0.00109204s 2.74254 %
M2M                         0.00319518s 8.02434 %
M2L                         0.0243386s 61.1237 %
L2L                         0.00345754s 8.68323 %
others (direct integration) 0.00674518s 16.9398 %
total                       0.0398186 s
--------------------------------------------------
Elapsed time for GMRES: {4.26445,4.26445} [s] for size = 100
GMRES expected error = 1.58634e-07 True error = 1.58634e-07
destructing gmres
Restart count :2
----- Average time for 50 iterations -----
update poles                0.00058864s 1.49416 %
reset Moments               0.00036198s 0.918823 %
increment_M_reuse           0.00101398s 2.57381 %
M2M                         0.00319462s 8.10899 %
M2L                         0.0241131s 61.207 %
L2L                         0.00342498s 8.69372 %
others (direct integration) 0.00669872s 17.0035 %
total                       0.039396 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00059106s 1.51168 %
reset Moments               0.0003617s 0.925072 %
increment_M_reuse           0.00098904s 2.52954 %
M2M                         0.003188s 8.15352 %
M2L                         0.0238574s 61.0169 %
L2L                         0.00342784s 8.76693 %
others (direct integration) 0.00668462s 17.0964 %
total                       0.0390997 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.0005947s 1.51741 %
reset Moments               0.00037354s 0.95311 %
increment_M_reuse           0.00099898s 2.54896 %
M2M                         0.0031913s 8.1428 %
M2L                         0.0239273s 61.0521 %
L2L                         0.00342346s 8.73517 %
others (direct integration) 0.00668238s 17.0505 %
total                       0.0391917 s
--------------------------------------------------
Elapsed time for GMRES: {6.42996,6.42996} [s] for size = 150
GMRES expected error = 1.07666e-13 True error = 1.0954e-13
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {0.000623266,1.58634e-07,1.07666e-13}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 15.6742 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {15.6743,49.437} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.023046,49.4601} s
ALEのU_update_BEMを計算
Elapsed time: 1.81498 s
name = tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.00958981,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 0, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 1, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 2, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 3, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
BVP.solveForPhiPhin_t-> {Φt,Φtn}とnet->accelerationが決まる
Elapsed time: {0.833844,52.1089} s
name = tankupdating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.00958981,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity

net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.400087,52.509} s
Total elapsed time: 52.509 s
============================================================================
ElapsedTimeBIEDiscretization={7.75767,15.6589,4.91771,15.6328}
ElapsedTimeSolve={0,0,0,0}
ElapsedTimeALE=1.81498
ElapsedTimeTotal=52.509
============================================================================
simulation_timeを取得
/Users/tomoaki/Ruehl20160d075_SimpleM2L/output/water_0.vtu  VV_points.size() : 52312 |||||||||
/Users/tomoaki/Ruehl20160d075_SimpleM2L/output/water.pvd
Creating /Users/tomoaki/Ruehl20160d075_SimpleM2L/output/water.pvd ...
Done.
/Users/tomoaki/Ruehl20160d075_SimpleM2L/output/tank_0.vtu  VV_points.size() : 52 |||||||||
/Users/tomoaki/Ruehl20160d075_SimpleM2L/output/tank.pvd
Creating /Users/tomoaki/Ruehl20160d075_SimpleM2L/output/tank.pvd ...
Done.
/Users/tomoaki/Ruehl20160d075_SimpleM2L/output/wavemaker_0.vtu  VV_points.size() : 496 |||||||||
/Users/tomoaki/Ruehl20160d075_SimpleM2L/output/wavemaker.pvd
Creating /Users/tomoaki/Ruehl20160d075_SimpleM2L/output/wavemaker.pvd ...
Done.
flipIf
flipIf
flipIf
net.getPoints() = 26158
Total : 26158
Total variables: 29016
Total case double-node : 29016
node reduction : 0.0984974
CORNER : 561
Total CORNER faces : 3419
Neumann : 18554
Dirichlet : 7043
===========================================================================
       dt :0.01
time_step :1
real time :0.01
---------------------------------------------------------------------------
RK_step = 1/4, RK_time = 0.01, simulation_time = 0.01
water makeBucketPoints(tank makeBucketPoints(2.01457)
wavemaker makeBucketPoints(2.99541)
0.616441)
BucketPoints.data1D.size() = 28
      bounds : {{-14.25,29.25},{-4.5,4.5},{-2.125,4.625}}
          dL : 2.99541
bucket sizes : {15, 4, 3}
tank makeBucketPoints done
tank makeBucketFaces(2.99541)
BucketPoints.data1D.size() = 250
      bounds : {{-6.24995,1.25005},{-1.5,1.5},{-0.75,3.75}}
          dL : 0.616441
bucket sizes : {13, 5, 8}
wavemaker makeBucketPoints done
wavemaker makeBucketFaces(0.616441)
      bounds : {{-14.25,29.25},{-4.5,4.5},{-2.125,4.625}}
          dL : 2.99541
bucket sizes : {15, 4, 3}
      bounds : {{-14.25,29.25},{-4.5,4.5},{-2.125,4.625}}
          dL : 2.99541
bucket sizes : {15, 4, 3}
tank makeBucketFaces done
tank makeBucketTetras(2.99541)
      bounds : {{-14.25,29.25},{-4.5,4.5},{-2.125,4.625}}
          dL : 2.99541
bucket sizes : {15, 4, 3}
tank makeBucketTetras done
      bounds : {{-6.24995,1.25005},{-1.5,1.5},{-0.75,3.75}}
          dL : 0.616441
bucket sizes : {13, 5, 8}
      bounds : {{-6.24995,1.25005},{-1.5,1.5},{-0.75,3.75}}
          dL : 0.616441
bucket sizes : {13, 5, 8}
wavemaker makeBucketFaces done
wavemaker makeBucketTetras(0.616441)
      bounds : {{-6.24995,1.25005},{-1.5,1.5},{-0.75,3.75}}
          dL : 0.616441
bucket sizes : {13, 5, 8}
wavemaker makeBucketTetras done
BucketPoints.data1D.size() = 26158
      bounds : {{-4.99994,25},{-1.5,1.5},{-0.340041,1.70021}}
          dL : 2.01457
bucket sizes : {15, 2, 2}
water makeBucketPoints done
water makeBucketFaces(2.01457)
      bounds : {{-4.99994,25},{-1.5,1.5},{-0.340041,1.70021}}
          dL : 2.01457
bucket sizes : {15, 2, 2}
      bounds : {{-4.99994,25},{-1.5,1.5},{-0.340041,1.70021}}
          dL : 2.01457
bucket sizes : {15, 2, 2}
water makeBucketFaces done
water makeBucketTetras(2.01457)
      bounds : {{-4.99994,25},{-1.5,1.5},{-0.340041,1.70021}}
          dL : 2.01457
bucket sizes : {15, 2, 2}
water makeBucketTetras done
makeBuckets
Elapsed time: {0.974463,0.974463} s
waterの境界条件を決定 setBoundaryTypes
water setContactFaces()
step2 面の境界条件を判定
step3 線の境界条件を決定
step4 点の境界条件を決定
setBoundaryTypes終了
setBoundaryTypes
Elapsed time: {0.686161,1.66062} s
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
Elapsed time: {4.4e-05,1.66067} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 31039
Setting up buckets...
Buckets set up complete.
Creating sources on surfaces...
バケットの作成．極の追加
added :52312
バケットの作成．極の追加, Elapsed time : {0.014302,0.014302}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 1 deepest level buckets.
After shrink: 1 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {1.6431,1.6574}
initializeFMM
SimpleM2L
[FMM:init] Step0: update poles start, #poles=52312
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=52312
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.010481,0.010481}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.001183,0.001183}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000185,0.001368}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.005557,0.006925}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.011796,0.018721}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.024339,0.04306}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {0.779756,0.822816}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000217,0.823033}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000192,0.823225}
setM2L, Elapsed time : {0.823237,0.823237}
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
level=7
[FMM:init] Debug: coupling stats printed
[FMM:init] Step4: setL2L start
[FMM:init] Step4: setL2L done
各レベルの各セルのL2Pの相手を保存する (#targets=26158)
setL2P, Elapsed time : {0.007547,0.007547} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=26158)
setDirectIntegration, Elapsed time : {0.623929,0.631476} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=997.236 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000815,0.000815}
reset Moments, Elapsed time : {0.000382,0.001197}
increment_M_reuse, Elapsed time : {0.001392,0.002589}
M2M, Elapsed time : {0.003614,0.006203}
M2L, Elapsed time : {0.024707,0.03091}
L2L, Elapsed time : {0.004566,0.035476}
updateFMM, Elapsed time : {4e-06,0.03548}
integration (direct integration is dominant), Elapsed time : {0.007458,0.007458}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00053226s 1.35503 %
reset Moments               0.00036456s 0.928098 %
increment_M_reuse           0.00101508s 2.58419 %
M2M                         0.0031483s 8.01495 %
M2L                         0.0237194s 60.3849 %
L2L                         0.0034218s 8.71123 %
others (direct integration) 0.00707894s 18.0216 %
total                       0.0392803 s
--------------------------------------------------
Elapsed time for GMRES: {2.0741,2.0741} [s] for size = 50
GMRES expected error = 6.53317e-06 True error = 6.53317e-06
destructing gmres
Restart count :1
----- Average time for 50 iterations -----
update poles                0.00055356s 1.3611 %
reset Moments               0.00038632s 0.949888 %
increment_M_reuse           0.00112904s 2.7761 %
M2M                         0.00322096s 7.91973 %
M2L                         0.024556s 60.3786 %
L2L                         0.00349304s 8.58873 %
others (direct integration) 0.00733112s 18.0258 %
total                       0.0406701 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00052356s 1.33013 %
reset Moments               0.0003619s 0.919427 %
increment_M_reuse           0.00114714s 2.91437 %
M2M                         0.00312322s 7.93472 %
M2L                         0.0236867s 60.1774 %
L2L                         0.00339736s 8.63118 %
others (direct integration) 0.00712158s 18.0928 %
total                       0.0393615 s
--------------------------------------------------
Elapsed time for GMRES: {4.27696,4.27696} [s] for size = 100
GMRES expected error = 1.4078e-09 True error = 1.4078e-09
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {6.53317e-06,1.4078e-09}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 10.1859 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {10.1859,11.8465} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.022365,11.8689} s
name = tank
net->velocityTranslational() = {updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.00958981,0,0}
0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
use wavemaker's (RigidBodyObject) predetermiend velocity
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 0, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 1, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 2, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える
setPhiPhin_t終了
j = 3, alpha = 1, Norm(func) = 0, Norm(BM.dX) = 0
BVP.solveForPhiPhin_t-> {Φt,Φtn}とnet->accelerationが決まる
Elapsed time: {0.741893,12.6108} s
name = tank
net->velocityTranslational() = {updating wavemaker's (RigidBodyObject) velocity0,0,0}
use tank's (RigidBodyObject) predetermiend velocity

name = wavemaker
net->velocityTranslational() = {0.00958981,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.380153,12.991} s
RK_step = 2/4, RK_time = 0.015, simulation_time = 0.01
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
net->velocity = {0.0143818,0,0,0,0,0}
setting acceleration
指定がないので加速度はdefault_acceleration
setNeumannVelocity
Elapsed time: {4.2e-05,12.991} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 31039
Creating sources on surfaces...
バケットの作成．極の追加
added :52312
バケットの作成．極の追加, Elapsed time : {0.570752,0.570752}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 392 deepest level buckets.
After shrink: 392 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.10966,0.680412}
initializeFMM
SimpleM2L
[FMM:init] Step0: update poles start, #poles=52312
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=52312
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.010515,0.010515}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.000713,0.000713}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000168,0.000881}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.00553,0.006411}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.011708,0.018119}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.022904,0.041023}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {0.739557,0.78058}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.0002,0.78078}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000169,0.780949}
setM2L, Elapsed time : {0.780956,0.780956}
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
level=7
[FMM:init] Debug: coupling stats printed
[FMM:init] Step4: setL2L start
[FMM:init] Step4: setL2L done
各レベルの各セルのL2Pの相手を保存する (#targets=26158)
setL2P, Elapsed time : {0.007046,0.007046} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=26158)
setDirectIntegration, Elapsed time : {0.606037,0.613083} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=997.236 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000948,0.000948}
reset Moments, Elapsed time : {0.000412,0.00136}
increment_M_reuse, Elapsed time : {0.000989,0.002349}
M2M, Elapsed time : {0.003278,0.005627}
M2L, Elapsed time : {0.025808,0.031435}
L2L, Elapsed time : {0.003624,0.035059}
updateFMM, Elapsed time : {6e-06,0.035065}
integration (direct integration is dominant), Elapsed time : {0.007042,0.007042}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00058374s 1.48191 %
reset Moments               0.00035912s 0.911678 %
increment_M_reuse           0.0010068s 2.55591 %
M2M                         0.00314722s 7.98967 %
M2L                         0.0237739s 60.3535 %
L2L                         0.00341876s 8.67901 %
others (direct integration) 0.00710158s 18.0284 %
total                       0.0393911 s
--------------------------------------------------
Elapsed time for GMRES: {2.07613,2.07613} [s] for size = 50
GMRES expected error = 0.000585097 True error = 0.000585097
destructing gmres
Restart count :1
