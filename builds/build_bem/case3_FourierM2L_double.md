ケース：case3
14840節点
方法:Fourier M2L double


初回（1step目のRKの1回目）は収束が早いので，2回目以降の計算を比較すること．

------- LOG -------

Last login: Thu Sep 18 13:27:52 on ttys039
/Users/tomoaki/.zshrc is loaded
build_bem(main)$ ./FourierM2L_double /Users/tomoaki/Ruehl20160d1_FourierM2L_double
input_directory : "/Users/tomoaki/Ruehl20160d1_FourierM2L_double"
ALE: {pseudo_quad}
ALEPERIOD: {1}
GRAVITY: {9.81}
WATER_DENSITY: {1000.0}
element: {linear}
end_time: {10.0}
end_time_step: {100000}
input_files: {tank.json,wavemaker.json,water.json,wg3.json,wg6.json}
max_dt: {0.1}
output_directory: {/Users/tomoaki/Ruehl20160d1_FourierM2L_double/output}
LINEAR_ELEMENT
ALE_ON_PSEUDO_QUADRATIC_ELEMENT

WATER_DENSITY: 1000
GRAVITY: 9.81
"/Users/tomoaki/Ruehl20160d1_FourierM2L_double/tank.json"
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
        this->getName(): tank, 0x150020000
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
/Users/tomoaki/Ruehl20160d1_FourierM2L_double/output/tank_init.vtu  VV_points.size() : 52 |||||||||
"/Users/tomoaki/Ruehl20160d1_FourierM2L_double/wavemaker.json"
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
        this->getName(): wavemaker, 0x152bf8000
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
/Users/tomoaki/Ruehl20160d1_FourierM2L_double/output/wavemaker_init.vtu  VV_points.size() : 496 |||||||||
"/Users/tomoaki/Ruehl20160d1_FourierM2L_double/water.json"
name: {water}
objfile: {/Users/tomoaki/Dropbox/code//cpp/obj/Ruehl2016//water0d1.obj}
reverseNormal: {false}
triangles: {12}
type: {Fluid}
vertices: {8}
/Users/tomoaki/Dropbox/code//cpp/obj/Ruehl2016//water0d1.obj is opened
/Users/tomoaki/Dropbox/code//cpp/obj/Ruehl2016//water0d1.obj is closed
Load3DFile : lines of file are splited
Load3DFile : file is loaded
Load3DFile : complex is generated
Load3DFile : max and min are calculated
/* -------------------------------------------------------------------------- */
        this->getName(): water, 0x157fd8000
           Lines.size(): 44514
          Points.size(): 14840
           Faces.size(): 29676
           this->bounds: {{0,20},{-1,1},{0,1.36}}
                              0     1     2     3     4     5     6     7     8     9    10    11    12    13    14
       Lines of points :      0     0     0     0     2  2110 10684  1986    58     0     0     0     0     0     0
   Connection of faces :      0     0 44514 全ての線が2つの面と接続しているため，格子は閉じた面を形成している
/* -------------------------------------------------------------------------- */
setOutputInfo
setTypes
Fluid
set velocity
/Users/tomoaki/Ruehl20160d1_FourierM2L_double/output/water_init.vtu  VV_points.size() : 29676 |||||||||
"/Users/tomoaki/Ruehl20160d1_FourierM2L_double/wg3.json"
name: {wg3}
position: {9.48,0.0,1.86,9.48,0.0,0.8600000000000001}
type: {wavegauge}
type = wavegauge
skipped
"/Users/tomoaki/Ruehl20160d1_FourierM2L_double/wg6.json"
name: {wg6}
position: {16.778,-0.031,2.0,16.778,-0.031,0.5000000000000001}
type: {wavegauge}
type = wavegauge
skipped
setting done
flipIf
flipIf
flipIf
net.getPoints() = 14840
Total : 14840
Total variables: 0
Total case double-node : 14840
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
water makeBucketPoints(tank makeBucketPoints(wavemaker makeBucketPoints(2.99541)
0.616441)2.01457)

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
BucketPoints.data1D.size() = 14840
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
Elapsed time: {0.488384,0.488384} s
waterの境界条件を決定 setBoundaryTypes
water setContactFaces()
step2 面の境界条件を判定
step3 線の境界条件を決定
step4 点の境界条件を決定
setBoundaryTypes終了
setBoundaryTypes
Elapsed time: {0.323331,0.811715} s
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
Elapsed time: {6.6e-05,0.811781} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18373
Setting up buckets...
Buckets set up complete.
Creating sources on surfaces...
バケットの作成．極の追加
added :29676
バケットの作成．極の追加, Elapsed time : {0.005325,0.005325}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 1 deepest level buckets.
After shrink: 1 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {3.03343,3.03875}
initializeFMM
FourierM2L_double
[FMM:init] Step0: update poles start, #poles=29676
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=29676
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.010364,0.010364}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.001449,0.001449}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.00019,0.001639}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.030942,0.032581}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.052015,0.084596}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.104402,0.188998}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {3.58349,3.77249}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000359,3.77284}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000297,3.77314}
setM2L, Elapsed time : {3.77315,3.77315}
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
各レベルの各セルのL2Pの相手を保存する (#targets=14840)
setL2P, Elapsed time : {0.0046,0.0046} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=14840)
setDirectIntegration, Elapsed time : {0.174399,0.178999} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=602.81 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000393,0.000393}
reset Moments, Elapsed time : {0.000381,0.000774}
increment_M_reuse, Elapsed time : {0.002675,0.003449}
M2M, Elapsed time : {0.003644,0.007093}
M2L, Elapsed time : {0.018709,0.025802}
L2L, Elapsed time : {0.004814,0.030616}
updateFMM, Elapsed time : {5e-06,0.030621}
integration (direct integration is dominant), Elapsed time : {0.002506,0.002506}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00024858s 0.728991 %
reset Moments               0.00032256s 0.945946 %
increment_M_reuse           0.00062624s 1.83652 %
M2M                         0.00412244s 12.0896 %
M2L                         0.0203485s 59.6745 %
L2L                         0.00591924s 17.3589 %
others (direct integration) 0.0025116s 7.36557 %
total                       0.0340992 s
--------------------------------------------------
Elapsed time for GMRES: {1.79707,1.79707} [s] for size = 50
GMRES expected error = 0 True error = 0
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {0}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 9.19518 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {9.19521,10.007} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.014155,10.0211} s
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name = tank
net->velocityTranslational() = {0,0,0}
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
Elapsed time: {0.394573,10.4157} s
name = tank
net->velocityTranslational() = {0updating wavemaker's (RigidBodyObject) velocity
,0,0}
use tank's (RigidBodyObject) predetermiend velocity
name = wavemaker
net->velocityTranslational() = {0,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.216037,10.6318} s
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
Elapsed time: {4.4e-05,10.6318} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18373
Creating sources on surfaces...
バケットの作成．極の追加
added :29676
バケットの作成．極の追加, Elapsed time : {0.338723,0.338723}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 392 deepest level buckets.
After shrink: 392 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.109396,0.448119}
initializeFMM
FourierM2L_double
[FMM:init] Step0: update poles start, #poles=29676
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=29676
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.011973,0.011973}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.000655,0.000655}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000166,0.000821}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.029812,0.030633}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.058005,0.088638}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.109736,0.198374}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {3.62554,3.82391}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000248,3.82416}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000183,3.82435}
setM2L, Elapsed time : {3.82435,3.82435}
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
各レベルの各セルのL2Pの相手を保存する (#targets=14840)
setL2P, Elapsed time : {0.003931,0.003931} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=14840)
setDirectIntegration, Elapsed time : {0.158973,0.162904} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=602.81 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000388,0.000388}
reset Moments, Elapsed time : {0.000366,0.000754}
increment_M_reuse, Elapsed time : {0.000573,0.001327}
M2M, Elapsed time : {0.003371,0.004698}
M2L, Elapsed time : {0.008205,0.012903}
L2L, Elapsed time : {0.003627,0.01653}
updateFMM, Elapsed time : {5e-06,0.016535}
integration (direct integration is dominant), Elapsed time : {0.002166,0.002166}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00025058s 1.33312 %
reset Moments               0.00030682s 1.63233 %
increment_M_reuse           0.00060388s 3.21273 %
M2M                         0.00310112s 16.4984 %
M2L                         0.00906524s 48.2284 %
L2L                         0.0033473s 17.8081 %
others (direct integration) 0.00212152s 11.2868 %
total                       0.0187965 s
--------------------------------------------------
Elapsed time for GMRES: {0.997799,0.997799} [s] for size = 50
GMRES expected error = 0.00032091 True error = 0.00032091
destructing gmres
Restart count :1
----- Average time for 50 iterations -----
update poles                0.00024884s 1.32062 %
reset Moments               0.00030254s 1.60561 %
increment_M_reuse           0.00060406s 3.2058 %
M2M                         0.00308722s 16.3842 %
M2L                         0.00915538s 48.5884 %
L2L                         0.0033381s 17.7156 %
others (direct integration) 0.00210658s 11.1798 %
total                       0.0188427 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00024684s 1.31229 %
reset Moments               0.0003028s 1.60979 %
increment_M_reuse           0.00059136s 3.14388 %
M2M                         0.00307282s 16.3362 %
M2L                         0.00913652s 48.573 %
L2L                         0.00334392s 17.7775 %
others (direct integration) 0.00211562s 11.2474 %
total                       0.0188099 s
--------------------------------------------------
Elapsed time for GMRES: {2.0315,2.0315} [s] for size = 100
GMRES expected error = 1.43793e-08 True error = 1.43793e-08
destructing gmres
Restart count :2
----- Average time for 50 iterations -----
update poles                0.00024742s 1.31062 %
reset Moments               0.0003023s 1.60133 %
increment_M_reuse           0.00060014s 3.17903 %
M2M                         0.00308742s 16.3545 %
M2L                         0.00914218s 48.4275 %
L2L                         0.00335326s 17.7627 %
others (direct integration) 0.00214534s 11.3642 %
total                       0.0188781 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00024614s 1.30921 %
reset Moments               0.00030076s 1.59973 %
increment_M_reuse           0.0006015s 3.19936 %
M2M                         0.0030782s 16.3729 %
M2L                         0.00912718s 48.5472 %
L2L                         0.00333916s 17.7609 %
others (direct integration) 0.00210768s 11.2107 %
total                       0.0188006 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00024724s 1.30962 %
reset Moments               0.00030468s 1.61388 %
increment_M_reuse           0.00060056s 3.18114 %
M2M                         0.0030819s 16.3247 %
M2L                         0.00917032s 48.5748 %
L2L                         0.00333582s 17.6697 %
others (direct integration) 0.00213824s 11.3262 %
total                       0.0188788 s
--------------------------------------------------
Elapsed time for GMRES: {3.12332,3.12332} [s] for size = 150
GMRES expected error = 4.65961e-16 True error = 1.54447e-12
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {0.00032091,1.43793e-08,4.65961e-16}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 11.015 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {11.015,21.6468} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.011958,21.6588} s
name = tank
net->velocityTranslational() = {0,0,0}
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {use tank's (RigidBodyObject) predetermiend velocity
0.0047955,0,0}
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
Elapsed time: {0.49031,22.1491} s
name = tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0047955,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.213709,22.3628} s
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
Elapsed time: {4.2e-05,22.3629} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18373
Creating sources on surfaces...
バケットの作成．極の追加
added :29676
バケットの作成．極の追加, Elapsed time : {0.31557,0.31557}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 392 deepest level buckets.
After shrink: 392 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.104201,0.419771}
initializeFMM
FourierM2L_double
[FMM:init] Step0: update poles start, #poles=29676
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=29676
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.011941,0.011941}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.000679,0.000679}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.00022,0.000899}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.029757,0.030656}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.055743,0.086399}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.104554,0.190953}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {3.56874,3.75969}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000349,3.76004}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000219,3.76026}
setM2L, Elapsed time : {3.76027,3.76027}
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
各レベルの各セルのL2Pの相手を保存する (#targets=14840)
setL2P, Elapsed time : {0.003934,0.003934} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=14840)
setDirectIntegration, Elapsed time : {0.164714,0.168648} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=602.81 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000379,0.000379}
reset Moments, Elapsed time : {0.000353,0.000732}
increment_M_reuse, Elapsed time : {0.00056,0.001292}
M2M, Elapsed time : {0.003245,0.004537}
M2L, Elapsed time : {0.008176,0.012713}
L2L, Elapsed time : {0.003406,0.016119}
updateFMM, Elapsed time : {6e-06,0.016125}
integration (direct integration is dominant), Elapsed time : {0.00214,0.00214}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00025262s 1.33877 %
reset Moments               0.00031262s 1.65675 %
increment_M_reuse           0.00060916s 3.22828 %
M2M                         0.00310388s 16.4492 %
M2L                         0.009061s 48.0192 %
L2L                         0.00338112s 17.9184 %
others (direct integration) 0.00214912s 11.3894 %
total                       0.0188695 s
--------------------------------------------------
Elapsed time for GMRES: {1.00067,1.00067} [s] for size = 50
GMRES expected error = 3.24012e-09 True error = 3.24014e-09
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {3.24012e-09}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 5.76734 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {5.76736,28.1302} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.012121,28.1423} s
name = tank
net->velocityTranslational() = {updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0047955,0,0}
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
Elapsed time: {0.436216,28.5786} s
name = tank
net->velocityTranslational() = {updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0047955,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.215219,28.7938} s
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
Elapsed time: {4.6e-05,28.7938} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18373
Creating sources on surfaces...
バケットの作成．極の追加
added :29676
バケットの作成．極の追加, Elapsed time : {0.368708,0.368708}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 392 deepest level buckets.
After shrink: 392 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.119869,0.488577}
initializeFMM
FourierM2L_double
[FMM:init] Step0: update poles start, #poles=29676
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=29676
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.010346,0.010346}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.000732,0.000732}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000138,0.00087}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.029579,0.030449}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.049709,0.080158}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.107889,0.188047}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {3.55349,3.74154}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000262,3.7418}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000215,3.74201}
setM2L, Elapsed time : {3.74202,3.74202}
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
各レベルの各セルのL2Pの相手を保存する (#targets=14840)
setL2P, Elapsed time : {0.004019,0.004019} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=14840)
setDirectIntegration, Elapsed time : {0.169789,0.173808} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=602.81 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000401,0.000401}
reset Moments, Elapsed time : {0.000374,0.000775}
increment_M_reuse, Elapsed time : {0.000589,0.001364}
M2M, Elapsed time : {0.003173,0.004537}
M2L, Elapsed time : {0.007865,0.012402}
L2L, Elapsed time : {0.00355,0.015952}
updateFMM, Elapsed time : {5e-06,0.015957}
integration (direct integration is dominant), Elapsed time : {0.002169,0.002169}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00025678s 1.36454 %
reset Moments               0.00030668s 1.62971 %
increment_M_reuse           0.0005983s 3.17939 %
M2M                         0.0031026s 16.4873 %
M2L                         0.00905236s 48.1045 %
L2L                         0.0033799s 17.9609 %
others (direct integration) 0.00212148s 11.2736 %
total                       0.0188181 s
--------------------------------------------------
Elapsed time for GMRES: {0.998237,0.998237} [s] for size = 50
GMRES expected error = 0.000320864 True error = 0.000320864
destructing gmres
Restart count :1
----- Average time for 50 iterations -----
update poles                0.00025274s 1.34253 %
reset Moments               0.00030258s 1.60727 %
increment_M_reuse           0.0005935s 3.15261 %
M2M                         0.00310912s 16.5153 %
M2L                         0.00905302s 48.0887 %
L2L                         0.00338682s 17.9904 %
others (direct integration) 0.0021279s 11.3032 %
total                       0.0188257 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.0002522s 1.33869 %
reset Moments               0.00030358s 1.61142 %
increment_M_reuse           0.00059548s 3.16083 %
M2M                         0.00309512s 16.429 %
M2L                         0.00908242s 48.2099 %
L2L                         0.00338498s 17.9676 %
others (direct integration) 0.00212556s 11.2826 %
total                       0.0188393 s
--------------------------------------------------
Elapsed time for GMRES: {2.03102,2.03102} [s] for size = 100
GMRES expected error = 1.43794e-08 True error = 1.43794e-08
destructing gmres
Restart count :2
----- Average time for 50 iterations -----
update poles                0.0002532s 1.34625 %
reset Moments               0.00030306s 1.61135 %
increment_M_reuse           0.00059794s 3.17922 %
M2M                         0.0031022s 16.4942 %
M2L                         0.00906694s 48.2085 %
L2L                         0.00336938s 17.9148 %
others (direct integration) 0.00211506s 11.2457 %
total                       0.0188078 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00025276s 1.34781 %
reset Moments               0.00030334s 1.61752 %
increment_M_reuse           0.0005922s 3.15783 %
M2M                         0.00309642s 16.5112 %
M2L                         0.00904034s 48.2064 %
L2L                         0.0033679s 17.9589 %
others (direct integration) 0.00210044s 11.2003 %
total                       0.0187534 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00025336s 1.34936 %
reset Moments               0.00030672s 1.63355 %
increment_M_reuse           0.00059424s 3.16484 %
M2M                         0.00308484s 16.4295 %
M2L                         0.00907664s 48.341 %
L2L                         0.00337172s 17.9573 %
others (direct integration) 0.00208876s 11.1245 %
total                       0.0187763 s
--------------------------------------------------
Elapsed time for GMRES: {3.11197,3.11197} [s] for size = 150
GMRES expected error = 4.65114e-16 True error = 3.16118e-12
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {0.000320864,1.43794e-08,4.65114e-16}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 10.9499 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {10.95,39.7438} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.010756,39.7545} s
ALEのU_update_BEMを計算
Elapsed time: 1.04511 s
name = tank
net->velocityTranslational() = {updating 0,wavemaker0,0}
use tank's (RigidBodyObject) predetermiend velocity
's (RigidBodyObject) velocity
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
Elapsed time: {0.435966,41.2356} s
name = tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.00958981,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.190248,41.4259} s
Total elapsed time: 41.4259 s
============================================================================
ElapsedTimeBIEDiscretization={9.17455,10.9936,5.74593,10.9283}
ElapsedTimeSolve={0,0,0,0}
ElapsedTimeALE=1.04511
ElapsedTimeTotal=41.4259
============================================================================
simulation_timeを取得
/Users/tomoaki/Ruehl20160d1_FourierM2L_double/output/water_0.vtu  VV_points.size() : 29676 |||||||||
/Users/tomoaki/Ruehl20160d1_FourierM2L_double/output/water.pvd
Creating /Users/tomoaki/Ruehl20160d1_FourierM2L_double/output/water.pvd ...
Done.
/Users/tomoaki/Ruehl20160d1_FourierM2L_double/output/tank_0.vtu  VV_points.size() : 52 |||||||||
/Users/tomoaki/Ruehl20160d1_FourierM2L_double/output/tank.pvd
Creating /Users/tomoaki/Ruehl20160d1_FourierM2L_double/output/tank.pvd ...
Done.
/Users/tomoaki/Ruehl20160d1_FourierM2L_double/output/wavemaker_0.vtu  VV_points.size() : 496 |||||||||
/Users/tomoaki/Ruehl20160d1_FourierM2L_double/output/wavemaker.pvd
Creating /Users/tomoaki/Ruehl20160d1_FourierM2L_double/output/wavemaker.pvd ...
Done.
flipIf
flipIf
flipIf
net.getPoints() = 14840
Total : 14840
Total variables: 16914
Total case double-node : 16914
node reduction : 0.12262
CORNER : 402
Total CORNER faces : 2476
Neumann : 10391
Dirichlet : 4047
===========================================================================
       dt :0.01
time_step :1
real time :0.01
---------------------------------------------------------------------------
RK_step = 1/4, RK_time = 0.01, simulation_time = 0.01
water makeBucketPoints(tank2.01457)wavemaker makeBucketPoints(
 makeBucketPoints(2.99541)
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
BucketPoints.data1D.size() = 14840
      bounds : {{-4.99994,25},{-1.5,1.5},{-0.340039,1.70019}}
          dL : 2.01457
bucket sizes : {15, 2, 2}
water makeBucketPoints done
water makeBucketFaces(2.01457)
      bounds : {{-4.99994,25},{-1.5,1.5},{-0.340039,1.70019}}
          dL : 2.01457
bucket sizes : {15, 2, 2}
      bounds : {{-4.99994,25},{-1.5,1.5},{-0.340039,1.70019}}
          dL : 2.01457
bucket sizes : {15, 2, 2}
water makeBucketFaces done
water makeBucketTetras(2.01457)
      bounds : {{-4.99994,25},{-1.5,1.5},{-0.340039,1.70019}}
          dL : 2.01457
bucket sizes : {15, 2, 2}
water makeBucketTetras done
makeBuckets
Elapsed time: {0.543236,0.543236} s
waterの境界条件を決定 setBoundaryTypes
water setContactFaces()
step2 面の境界条件を判定
step3 線の境界条件を決定
step4 点の境界条件を決定
setBoundaryTypes終了
setBoundaryTypes
Elapsed time: {0.30246,0.845696} s
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
Elapsed time: {4.3e-05,0.845739} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18373
Setting up buckets...
Buckets set up complete.
Creating sources on surfaces...
バケットの作成．極の追加
added :29676
バケットの作成．極の追加, Elapsed time : {0.005645,0.005645}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 1 deepest level buckets.
After shrink: 1 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {3.02767,3.03331}
initializeFMM
FourierM2L_double
[FMM:init] Step0: update poles start, #poles=29676
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=29676
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.010351,0.010351}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.001069,0.001069}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000198,0.001267}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.029846,0.031113}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.051239,0.082352}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.106524,0.188876}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {3.56303,3.7519}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000224,3.75213}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.00017,3.7523}
setM2L, Elapsed time : {3.75231,3.75231}
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
各レベルの各セルのL2Pの相手を保存する (#targets=14840)
setL2P, Elapsed time : {0.004076,0.004076} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=14840)
setDirectIntegration, Elapsed time : {0.172628,0.176704} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=602.823 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000407,0.000407}
reset Moments, Elapsed time : {0.000389,0.000796}
increment_M_reuse, Elapsed time : {0.000827,0.001623}
M2M, Elapsed time : {0.003298,0.004921}
M2L, Elapsed time : {0.008906,0.013827}
L2L, Elapsed time : {0.003733,0.01756}
updateFMM, Elapsed time : {5e-06,0.017565}
integration (direct integration is dominant), Elapsed time : {0.002786,0.002786}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00026454s 1.36608 %
reset Moments               0.00031818s 1.64308 %
increment_M_reuse           0.00060474s 3.12287 %
M2M                         0.0031288s 16.1571 %
M2L                         0.00903528s 46.658 %
L2L                         0.00343472s 17.7368 %
others (direct integration) 0.00257864s 13.3161 %
total                       0.0193649 s
--------------------------------------------------
Elapsed time for GMRES: {1.02857,1.02857} [s] for size = 50
GMRES expected error = 2.50557e-06 True error = 2.50557e-06
destructing gmres
Restart count :1
----- Average time for 50 iterations -----
update poles                0.00026038s 1.34651 %
reset Moments               0.00030846s 1.59515 %
increment_M_reuse           0.00060294s 3.118 %
M2M                         0.00310896s 16.0774 %
M2L                         0.00905628s 46.833 %
L2L                         0.003413s 17.6497 %
others (direct integration) 0.00258738s 13.3802 %
total                       0.0193374 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00025762s 1.3389 %
reset Moments               0.00030824s 1.60198 %
increment_M_reuse           0.00059658s 3.10054 %
M2M                         0.00309838s 16.1029 %
M2L                         0.00902598s 46.9097 %
L2L                         0.00340242s 17.683 %
others (direct integration) 0.00255194s 13.2629 %
total                       0.0192412 s
--------------------------------------------------
Elapsed time for GMRES: {2.08086,2.08086} [s] for size = 100
GMRES expected error = 7.30481e-11 True error = 7.30486e-11
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {2.50557e-06,7.30481e-11}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 10.4294 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {10.4294,11.2751} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.011421,11.2865} s
name = tank
net->velocityTranslational() = {updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.00958981,0,0}
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
Elapsed time: {0.430374,11.7169} s
name = tank
net->velocityTranslational() = {0,0,0}
updating use wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.00958981,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.213464,11.9304} s
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
Elapsed time: {3.9e-05,11.9304} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18373
Creating sources on surfaces...
バケットの作成．極の追加
added :29676
バケットの作成．極の追加, Elapsed time : {0.377037,0.377037}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 392 deepest level buckets.
After shrink: 392 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.11196,0.488997}
initializeFMM
FourierM2L_double
[FMM:init] Step0: update poles start, #poles=29676
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=29676
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.011965,0.011965}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.000647,0.000647}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000175,0.000822}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.030035,0.030857}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.051663,0.08252}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.102466,0.184986}
