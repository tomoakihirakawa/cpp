ケース：case5
63291節点
方法:SimpleM2L

初回は収束が早いので2回目からの計測結果を比較すること．

------------ LOG -----------


Last login: Thu Sep 18 15:04:56 on ttys019
/Users/tomoaki/.zshrc is loaded
build_bem(main)$ ./SimpleM2L /Users/tomoaki/Ruehl20160d05_SimpleM2L
input_directory : "/Users/tomoaki/Ruehl20160d05_SimpleM2L"
ALE: {pseudo_quad}
ALEPERIOD: {1}
GRAVITY: {9.81}
WATER_DENSITY: {1000.0}
element: {linear}
end_time: {10.0}
end_time_step: {100000}
input_files: {tank.json,wavemaker.json,water.json,wg3.json,wg6.json}
max_dt: {0.1}
output_directory: {/Users/tomoaki/Ruehl20160d05_SimpleM2L/output}
LINEAR_ELEMENT
ALE_ON_PSEUDO_QUADRATIC_ELEMENT

WATER_DENSITY: 1000
GRAVITY: 9.81
"/Users/tomoaki/Ruehl20160d05_SimpleM2L/tank.json"
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
        this->getName(): tank, 0x150008000
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
/Users/tomoaki/Ruehl20160d05_SimpleM2L/output/tank_init.vtu  VV_points.size() : 52 |||||||||
"/Users/tomoaki/Ruehl20160d05_SimpleM2L/wavemaker.json"
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
        this->getName(): wavemaker, 0x152bd8000
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
/Users/tomoaki/Ruehl20160d05_SimpleM2L/output/wavemaker_init.vtu  VV_points.size() : 496 |||||||||
"/Users/tomoaki/Ruehl20160d05_SimpleM2L/water.json"
name: {water}
objfile: {/Users/tomoaki/Dropbox/code//cpp/obj/Ruehl2016//water0d05.obj}
reverseNormal: {false}
triangles: {12}
type: {Fluid}
vertices: {8}
/Users/tomoaki/Dropbox/code//cpp/obj/Ruehl2016//water0d05.obj is opened
/Users/tomoaki/Dropbox/code//cpp/obj/Ruehl2016//water0d05.obj is closed
Load3DFile : lines of file are splited
Load3DFile : file is loaded
Load3DFile : complex is generated
Load3DFile : max and min are calculated
/* -------------------------------------------------------------------------- */
        this->getName(): water, 0x152c58000
           Lines.size(): 189867
          Points.size(): 63291
           Faces.size(): 126578
           this->bounds: {{0,20},{-1,1},{0,1.36}}
                              0     1     2     3     4     5     6     7     8     9    10    11    12    13    14
       Lines of points :      0     0     0     0     1  7123 49219  6783   165     0     0     0     0     0     0
   Connection of faces :      0     0189867 全ての線が2つの面と接続しているため，格子は閉じた面を形成している
/* -------------------------------------------------------------------------- */
setOutputInfo
setTypes
Fluid
set velocity
/Users/tomoaki/Ruehl20160d05_SimpleM2L/output/water_init.vtu  VV_points.size() : 126578 |||||||||
"/Users/tomoaki/Ruehl20160d05_SimpleM2L/wg3.json"
name: {wg3}
position: {9.48,0.0,1.86,9.48,0.0,0.8600000000000001}
type: {wavegauge}
type = wavegauge
skipped
"/Users/tomoaki/Ruehl20160d05_SimpleM2L/wg6.json"
name: {wg6}
position: {16.778,-0.031,2.0,16.778,-0.031,0.5000000000000001}
type: {wavegauge}
type = wavegauge
skipped
setting done
flipIf
flipIf
flipIf
net.getPoints() = 63291
Total : 63291
Total variables: 0
Total case double-node : 63291
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
water makeBucketPoints(tank makeBucketPoints(wavemaker makeBucketPoints(0.616441)
2.01457)
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
BucketPoints.data1D.size() = 63291
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
Elapsed time: {2.08517,2.08517} s
waterの境界条件を決定 setBoundaryTypes
water setContactFaces()
step2 面の境界条件を判定
step3 線の境界条件を決定
step4 点の境界条件を決定
setBoundaryTypes終了
setBoundaryTypes
Elapsed time: {1.73202,3.81719} s
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
Elapsed time: {4.9e-05,3.81724} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 71077
Setting up buckets...
Buckets set up complete.
Creating sources on surfaces...
バケットの作成．極の追加
added :126578
バケットの作成．極の追加, Elapsed time : {0.034656,0.034656}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 1 deepest level buckets.
After shrink: 1 deepest level buckets.
After grow: 1136 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {5.47625,5.5109}
initializeFMM
SimpleM2L
[FMM:init] Step0: update poles start, #poles=126578
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=1505, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=1136, total leaf sources=126578
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.045124,0.045124}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.003165,0.003165}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000183,0.003348}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.006242,0.00959}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.013413,0.023003}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.024798,0.047801}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {0.776302,0.824103}
setM2L, level=6, buckets_at_a_level.size()=992, Elapsed time : {0.798251,1.62235}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000249,1.6226}
setM2L, Elapsed time : {1.62261,1.62261}
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
setL2P, Elapsed time : {0.018974,0.018974} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=63291)
setDirectIntegration, Elapsed time : {1.73501,1.75398} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=1118.63 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.001913,0.001913}
reset Moments, Elapsed time : {0.001212,0.003125}
increment_M_reuse, Elapsed time : {0.005556,0.008681}
M2M, Elapsed time : {0.010973,0.019654}
M2L, Elapsed time : {0.074818,0.094472}
L2L, Elapsed time : {0.014235,0.108707}
updateFMM, Elapsed time : {6e-06,0.108713}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00131858s 0.861728 %
reset Moments               0.00111582s 0.729219 %
increment_M_reuse           0.00221694s 1.44883 %
M2M                         0.0116152s 7.59085 %
M2L                         0.0996848s 65.1467 %
L2L                         0.0177893s 11.6258 %
others (direct integration) 0.0192752s 12.5969 %
total                       0.153016 s
--------------------------------------------------
Elapsed time for GMRES: {8.03465,8.03465} [s] for size = 50
GMRES expected error = 0 True error = 0
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {0}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: {18.6162,18.6162} s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {18.6162,22.4335} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.056992,22.4904} s
name = tankupdating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {
net->velocityTranslational() = {0,0,0,0}
0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
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
Elapsed time: {1.56476,24.0552} s
name = tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.968162,25.0234} s
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
Elapsed time: {4.5e-05,25.0234} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 71077
Creating sources on surfaces...
バケットの作成．極の追加
added :126578
バケットの作成．極の追加, Elapsed time : {6.17972,6.17972}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 1136 deepest level buckets.
After shrink: 1136 deepest level buckets.
After grow: 1136 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.724272,6.90399}
initializeFMM
SimpleM2L
[FMM:init] Step0: update poles start, #poles=126578
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=1505, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=1136, total leaf sources=126578
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.05096,0.05096}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.001534,0.001534}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000162,0.001696}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.005498,0.007194}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.010232,0.017426}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.022056,0.039482}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {0.730334,0.769816}
setM2L, level=6, buckets_at_a_level.size()=992, Elapsed time : {0.824485,1.5943}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000215,1.59452}
setM2L, Elapsed time : {1.59452,1.59452}
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
setL2P, Elapsed time : {0.016438,0.016438} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=63291)
setDirectIntegration, Elapsed time : {1.63128,1.64772} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=1118.63 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.001812,0.001812}
reset Moments, Elapsed time : {0.001216,0.003028}
increment_M_reuse, Elapsed time : {0.002326,0.005354}
M2M, Elapsed time : {0.009026,0.01438}
M2L, Elapsed time : {0.051899,0.066279}
L2L, Elapsed time : {0.010438,0.076717}
updateFMM, Elapsed time : {9e-06,0.076726}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00133414s 1.50158 %
reset Moments               0.0011165s 1.25662 %
increment_M_reuse           0.00223562s 2.51619 %
M2M                         0.00921154s 10.3676 %
M2L                         0.0471735s 53.0938 %
L2L                         0.00990378s 11.1467 %
others (direct integration) 0.0178742s 20.1175 %
total                       0.0888493 s
--------------------------------------------------
Elapsed time for GMRES: {4.68352,4.68352} [s] for size = 50
GMRES expected error = 0.00143084 True error = 0.00143084
destructing gmres
Restart count :1
----- Average time for 50 iterations -----
update poles                0.00137088s 1.53408 %
reset Moments               0.0011268s 1.26094 %
increment_M_reuse           0.00223868s 2.50518 %
M2M                         0.00923028s 10.3291 %
M2L                         0.04717s 52.7853 %
L2L                         0.00993912s 11.1223 %
others (direct integration) 0.0182862s 20.463 %
total                       0.089362 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00136884s 1.50768 %
reset Moments               0.0011749s 1.29407 %
increment_M_reuse           0.00233718s 2.57424 %
M2M                         0.00932304s 10.2687 %
M2L                         0.0482792s 53.1761 %
L2L                         0.0101285s 11.1558 %
others (direct integration) 0.0181795s 20.0235 %
total                       0.0907912 s
--------------------------------------------------
Elapsed time for GMRES: {9.63381,9.63381} [s] for size = 100
GMRES expected error = 2.09159e-06 True error = 2.09159e-06
destructing gmres
Restart count :2
----- Average time for 50 iterations -----
update poles                0.00134196s 1.49834 %
reset Moments               0.0011176s 1.24784 %
increment_M_reuse           0.00225684s 2.51984 %
M2M                         0.00919526s 10.2668 %
M2L                         0.0476982s 53.2566 %
L2L                         0.00988416s 11.036 %
others (direct integration) 0.0180689s 20.1745 %
total                       0.0895628 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.0012544s 1.39108 %
reset Moments               0.00112358s 1.24601 %
increment_M_reuse           0.00228442s 2.53333 %
M2M                         0.009199s 10.2013 %
M2L                         0.0478684s 53.0842 %
L2L                         0.00990742s 10.9869 %
others (direct integration) 0.0185373s 20.5572 %
total                       0.0901746 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.0011899s 1.34116 %
reset Moments               0.00111828s 1.26044 %
increment_M_reuse           0.00222336s 2.50599 %
M2M                         0.0092277s 10.4007 %
M2L                         0.0472277s 53.2312 %
L2L                         0.0098891s 11.1462 %
others (direct integration) 0.0178457s 20.1143 %
total                       0.0887217 s
--------------------------------------------------
Elapsed time for GMRES: {14.6397,14.6397} [s] for size = 150
GMRES expected error = 6.19521e-11 True error = 6.19522e-11
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {0.00143084,2.09159e-06,6.19521e-11}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: {41.0472,41.0472} s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {41.0473,66.0707} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.052299,66.123} s
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
Elapsed time: {2.3431,68.4661} s
name = updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0047955,0,0}
tank
net->velocityTranslational() = {use wavemaker's (RigidBodyObject) predetermiend velocity
0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.993095,69.4592} s
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
Elapsed time: {5.5e-05,69.4592} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 71077
Creating sources on surfaces...
バケットの作成．極の追加
added :126578
バケットの作成．極の追加, Elapsed time : {6.36508,6.36508}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 1136 deepest level buckets.
After shrink: 1136 deepest level buckets.
After grow: 1136 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.733343,7.09842}
initializeFMM
SimpleM2L
[FMM:init] Step0: update poles start, #poles=126578
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=1505, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=1136, total leaf sources=126578
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.050574,0.050574}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.001421,0.001421}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000139,0.00156}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.005517,0.007077}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.010656,0.017733}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.021646,0.039379}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {0.726636,0.766015}
setM2L, level=6, buckets_at_a_level.size()=992, Elapsed time : {0.784269,1.55028}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000236,1.55052}
setM2L, Elapsed time : {1.55053,1.55053}
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
setL2P, Elapsed time : {0.016698,0.016698} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=63291)
setDirectIntegration, Elapsed time : {1.6317,1.6484} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=1118.63 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.001827,0.001827}
reset Moments, Elapsed time : {0.001189,0.003016}
increment_M_reuse, Elapsed time : {0.002386,0.005402}
M2M, Elapsed time : {0.009388,0.01479}
M2L, Elapsed time : {0.046054,0.060844}
L2L, Elapsed time : {0.01041,0.071254}
updateFMM, Elapsed time : {6e-06,0.07126}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00138806s 1.5128 %
reset Moments               0.00112716s 1.22845 %
increment_M_reuse           0.0023213s 2.5299 %
M2M                         0.00915316s 9.97571 %
M2L                         0.0490021s 53.4057 %
L2L                         0.00990976s 10.8003 %
others (direct integration) 0.0188529s 20.5471 %
total                       0.0917544 s
--------------------------------------------------
Elapsed time for GMRES: {4.83081,4.83081} [s] for size = 50
GMRES expected error = 2.63959e-08 True error = 2.63959e-08
destructing gmres
Restart count :1
----- Average time for 50 iterations -----
update poles                0.0013644s 1.53621 %
reset Moments               0.00111724s 1.25793 %
increment_M_reuse           0.00225612s 2.54022 %
M2M                         0.00924476s 10.4089 %
M2L                         0.0470645s 52.9911 %
L2L                         0.00987538s 11.1189 %
others (direct integration) 0.0178935s 20.1468 %
total                       0.088816 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.001264s 1.40326 %
reset Moments               0.0011199s 1.24328 %
increment_M_reuse           0.00224756s 2.49518 %
M2M                         0.00924174s 10.2599 %
M2L                         0.0478919s 53.1684 %
L2L                         0.00987998s 10.9685 %
others (direct integration) 0.0184308s 20.4614 %
total                       0.090076 s
--------------------------------------------------
Elapsed time for GMRES: {9.56903,9.56903} [s] for size = 100
GMRES expected error = 3.08123e-11 True error = 3.08122e-11
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {2.63959e-08,3.08123e-11}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: {26.6438,26.6438} s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {26.6438,96.103} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.053478,96.1565} s
name = tankupdating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {
net->velocityTranslational() = {0,0,0}
0.0047955,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
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
Elapsed time: {1.86866,98.0252} s
name = tank
net->velocityTranslational() = {updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0047955,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.979232,99.0044} s
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
Elapsed time: {4.3e-05,99.0044} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 71077
Creating sources on surfaces...
バケットの作成．極の追加
added :126578
バケットの作成．極の追加, Elapsed time : {6.47088,6.47088}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 1136 deepest level buckets.
After shrink: 1136 deepest level buckets.
After grow: 1136 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.741063,7.21195}
initializeFMM
SimpleM2L
[FMM:init] Step0: update poles start, #poles=126578
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=1505, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=1136, total leaf sources=126578
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.049809,0.049809}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.001404,0.001404}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000183,0.001587}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.005566,0.007153}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.010833,0.017986}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.02299,0.040976}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {0.747606,0.788582}
setM2L, level=6, buckets_at_a_level.size()=992, Elapsed time : {0.777405,1.56599}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000192,1.56618}
setM2L, Elapsed time : {1.56618,1.56618}
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
setL2P, Elapsed time : {0.01672,0.01672} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=63291)
setDirectIntegration, Elapsed time : {1.6707,1.68742} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=1118.63 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.001762,0.001762}
reset Moments, Elapsed time : {0.001199,0.002961}
increment_M_reuse, Elapsed time : {0.002466,0.005427}
M2M, Elapsed time : {0.00932,0.014747}
M2L, Elapsed time : {0.048567,0.063314}
L2L, Elapsed time : {0.010373,0.073687}
updateFMM, Elapsed time : {5e-06,0.073692}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00133544s 1.49283 %
reset Moments               0.00112472s 1.25728 %
increment_M_reuse           0.00228356s 2.5527 %
M2M                         0.00925048s 10.3407 %
M2L                         0.0476518s 53.268 %
L2L                         0.00994544s 11.1176 %
others (direct integration) 0.0178652s 19.9708 %
total                       0.0894567 s
--------------------------------------------------
Elapsed time for GMRES: {4.70292,4.70292} [s] for size = 50
GMRES expected error = 0.00142965 True error = 0.00142965
destructing gmres
Restart count :1
----- Average time for 50 iterations -----
update poles                0.00133362s 1.50842 %
reset Moments               0.00111614s 1.26243 %
increment_M_reuse           0.00223196s 2.52451 %
M2M                         0.00924772s 10.4598 %
M2L                         0.0471057s 53.2799 %
L2L                         0.0098734s 11.1675 %
others (direct integration) 0.0175032s 19.7974 %
total                       0.0884118 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00123242s 1.4172 %
reset Moments               0.00111732s 1.28484 %
increment_M_reuse           0.00220768s 2.53868 %
M2M                         0.00923142s 10.6155 %
M2L                         0.0460582s 52.9637 %
L2L                         0.00985174s 11.3288 %
others (direct integration) 0.017263s 19.8513 %
total                       0.0869618 s
--------------------------------------------------
Elapsed time for GMRES: {9.37934,9.37934} [s] for size = 100
GMRES expected error = 2.08782e-06 True error = 2.08782e-06
destructing gmres
Restart count :2
----- Average time for 50 iterations -----
update poles                0.00131978s 1.48228 %
reset Moments               0.00113424s 1.2739 %
increment_M_reuse           0.00223354s 2.50855 %
M2M                         0.00930358s 10.4491 %
M2L                         0.0475395s 53.3929 %
L2L                         0.00990518s 11.1248 %
others (direct integration) 0.0176014s 19.7686 %
total                       0.0890372 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.0012474s 1.41377 %
reset Moments               0.0011326s 1.28366 %
increment_M_reuse           0.00222968s 2.52706 %
M2M                         0.0092738s 10.5107 %
M2L                         0.0471248s 53.4101 %
L2L                         0.00989038s 11.2095 %
others (direct integration) 0.0173334s 19.6452 %
total                       0.088232 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00122064s 1.36923 %
reset Moments               0.00113406s 1.27211 %
increment_M_reuse           0.00223472s 2.50676 %
M2M                         0.00927516s 10.4043 %
M2L                         0.0476922s 53.4979 %
L2L                         0.0099155s 11.1225 %
others (direct integration) 0.0176755s 19.8272 %
total                       0.0891478 s
--------------------------------------------------
Elapsed time for GMRES: {14.5371,14.5371} [s] for size = 150
GMRES expected error = 6.1734e-11 True error = 6.17338e-11
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {0.00142965,2.08782e-06,6.1734e-11}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: {41.1571,41.1571} s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {41.1571,140.162} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.04743,140.209} s
ALEのU_update_BEMを計算
Elapsed time: 4.30137 s
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
Elapsed time: {2.29761,146.808} s
name = tank
net->velocityTranslational() = {0,0,0}
updating use tank's (RigidBodyObject) predetermiend velocity
wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.00958981,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {1.12513,147.933} s
Total elapsed time: 147.933 s
============================================================================
ElapsedTimeBIEDiscretization={18.5072,40.9477,26.5397,41.0551}
ElapsedTimeSolve={0,0,0,0}
ElapsedTimeALE=4.30137
ElapsedTimeTotal=147.933
============================================================================
simulation_timeを取得
/Users/tomoaki/Ruehl20160d05_SimpleM2L/output/water_0.vtu  VV_points.size() : 126578 |||||||||
/Users/tomoaki/Ruehl20160d05_SimpleM2L/output/water.pvd
Creating /Users/tomoaki/Ruehl20160d05_SimpleM2L/output/water.pvd ...
Done.
/Users/tomoaki/Ruehl20160d05_SimpleM2L/output/tank_0.vtu  VV_points.size() : 52 |||||||||
/Users/tomoaki/Ruehl20160d05_SimpleM2L/output/tank.pvd
Creating /Users/tomoaki/Ruehl20160d05_SimpleM2L/output/tank.pvd ...
Done.
/Users/tomoaki/Ruehl20160d05_SimpleM2L/output/wavemaker_0.vtu  VV_points.size() : 496 |||||||||
/Users/tomoaki/Ruehl20160d05_SimpleM2L/output/wavemaker.pvd
Creating /Users/tomoaki/Ruehl20160d05_SimpleM2L/output/wavemaker.pvd ...
Done.
flipIf
