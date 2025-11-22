ケース：case4
26158節点
方法:Fourier M2L double

初回（1step目のRKの1回目）は収束が早いので，2回目以降の計算を比較すること．

------- LOG -------

Last login: Thu Sep 18 13:30:27 on ttys041
/Users/tomoaki/.zshrc is loaded
build_bem(main)$ ./FourierM2L_double /Users/tomoaki/Ruehl20160d075_FourierM2L_double
input_directory : "/Users/tomoaki/Ruehl20160d075_FourierM2L_double"
ALE: {pseudo_quad}
ALEPERIOD: {1}
GRAVITY: {9.81}
WATER_DENSITY: {1000.0}
element: {linear}
end_time: {10.0}
end_time_step: {100000}
input_files: {tank.json,wavemaker.json,water.json,wg3.json,wg6.json}
max_dt: {0.1}
output_directory: {/Users/tomoaki/Ruehl20160d075_FourierM2L_double/output}
LINEAR_ELEMENT
ALE_ON_PSEUDO_QUADRATIC_ELEMENT

WATER_DENSITY: 1000
GRAVITY: 9.81
"/Users/tomoaki/Ruehl20160d075_FourierM2L_double/tank.json"
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
        this->getName(): tank, 0x108008000
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
/Users/tomoaki/Ruehl20160d075_FourierM2L_double/output/tank_init.vtu  VV_points.size() : 52 |||||||||
"/Users/tomoaki/Ruehl20160d075_FourierM2L_double/wavemaker.json"
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
        this->getName(): wavemaker, 0x10abe0000
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
/Users/tomoaki/Ruehl20160d075_FourierM2L_double/output/wavemaker_init.vtu  VV_points.size() : 496 |||||||||
"/Users/tomoaki/Ruehl20160d075_FourierM2L_double/water.json"
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
        this->getName(): water, 0x13ffd0000
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
/Users/tomoaki/Ruehl20160d075_FourierM2L_double/output/water_init.vtu  VV_points.size() : 52312 |||||||||
"/Users/tomoaki/Ruehl20160d075_FourierM2L_double/wg3.json"
name: {wg3}
position: {9.48,0.0,1.86,9.48,0.0,0.8600000000000001}
type: {wavegauge}
type = wavegauge
skipped
"/Users/tomoaki/Ruehl20160d075_FourierM2L_double/wg6.json"
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
waterwavemaker makeBucketPoints(0.616441)
tank makeBucketPoints( makeBucketPoints(2.01457)
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
Elapsed time: {0.866672,0.866672} s
waterの境界条件を決定 setBoundaryTypes
water setContactFaces()
step2 面の境界条件を判定
step3 線の境界条件を決定
step4 点の境界条件を決定
setBoundaryTypes終了
setBoundaryTypes
Elapsed time: {0.661512,1.52818} s
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
Elapsed time: {5.1e-05,1.52824} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 31039
Setting up buckets...
Buckets set up complete.
Creating sources on surfaces...
バケットの作成．極の追加
added :52312
バケットの作成．極の追加, Elapsed time : {0.014294,0.014294}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 1 deepest level buckets.
After shrink: 1 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {3.04854,3.06284}
initializeFMM
FourierM2L_double
[FMM:init] Step0: update poles start, #poles=52312
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=52312
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.010258,0.010258}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.001333,0.001333}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000162,0.001495}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.030314,0.031809}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.051532,0.083341}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.107044,0.190385}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {3.56425,3.75463}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000186,3.75482}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000283,3.7551}
setM2L, Elapsed time : {3.75511,3.75511}
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
setL2P, Elapsed time : {0.007652,0.007652} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=26158)
setDirectIntegration, Elapsed time : {0.592705,0.600357} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=997.897 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000748,0.000748}
reset Moments, Elapsed time : {0.000395,0.001143}
increment_M_reuse, Elapsed time : {0.004622,0.005765}
M2M, Elapsed time : {0.004103,0.009868}
M2L, Elapsed time : {0.018096,0.027964}
L2L, Elapsed time : {0.005362,0.033326}
updateFMM, Elapsed time : {5e-06,0.033331}
integration (direct integration is dominant), Elapsed time : {0.016433,0.016433}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00047024s 1.27389 %
reset Moments               0.00035122s 0.951459 %
increment_M_reuse           0.0009732s 2.63641 %
M2M                         0.00425354s 11.5229 %
M2L                         0.0179561s 48.6432 %
L2L                         0.00596258s 16.1527 %
others (direct integration) 0.00694698s 18.8195 %
total                       0.0369138 s
--------------------------------------------------
Elapsed time for GMRES: {1.96812,1.96812} [s] for size = 50
GMRES expected error = 0 True error = 0
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {0}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 10.0141 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {10.0142,11.5424} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.022944,11.5654} s
name = tank
net->velocityTranslational() = {0,0,0}
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0,0,0}
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
Elapsed time: {0.632613,12.198} s
name = tank
net->velocityTranslational() = updating wavemaker's (RigidBodyObject) velocity
{name = wavemaker
net->velocityTranslational() = {0,0,00}
use tank's (RigidBodyObject) predetermiend velocity
,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.379873,12.5778} s
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
Elapsed time: {4.9e-05,12.5779} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 31039
Creating sources on surfaces...
バケットの作成．極の追加
added :52312
バケットの作成．極の追加, Elapsed time : {0.566269,0.566269}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 392 deepest level buckets.
After shrink: 392 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.110928,0.677197}
initializeFMM
FourierM2L_double
[FMM:init] Step0: update poles start, #poles=52312
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=52312
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.011857,0.011857}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.000672,0.000672}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000225,0.000897}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.030236,0.031133}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.051505,0.082638}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.106155,0.188793}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {3.53888,3.72767}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000353,3.72803}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000318,3.72834}
setM2L, Elapsed time : {3.72837,3.72837}
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
setL2P, Elapsed time : {0.006855,0.006855} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=26158)
setDirectIntegration, Elapsed time : {0.554409,0.561264} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=997.897 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000736,0.000736}
reset Moments, Elapsed time : {0.000386,0.001122}
increment_M_reuse, Elapsed time : {0.000937,0.002059}
M2M, Elapsed time : {0.003337,0.005396}
M2L, Elapsed time : {0.008383,0.013779}
L2L, Elapsed time : {0.003612,0.017391}
updateFMM, Elapsed time : {5e-06,0.017396}
integration (direct integration is dominant), Elapsed time : {0.0063,0.0063}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.0004723s 1.99872 %
reset Moments               0.00034334s 1.45298 %
increment_M_reuse           0.00094306s 3.99093 %
M2M                         0.00312898s 13.2415 %
M2L                         0.0090001s 38.0875 %
L2L                         0.00338962s 14.3445 %
others (direct integration) 0.00635268s 26.8839 %
total                       0.0236301 s
--------------------------------------------------
Elapsed time for GMRES: {1.27186,1.27186} [s] for size = 50
GMRES expected error = 0.000623322 True error = 0.000623322
destructing gmres
Restart count :1
----- Average time for 50 iterations -----
update poles                0.00046454s 1.96649 %
reset Moments               0.00034488s 1.45994 %
increment_M_reuse           0.00095666s 4.04972 %
M2M                         0.00311728s 13.196 %
M2L                         0.00903522s 38.2478 %
L2L                         0.00338334s 14.3223 %
others (direct integration) 0.00632094s 26.7577 %
total                       0.0236229 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00046608s 1.96809 %
reset Moments               0.00034922s 1.47463 %
increment_M_reuse           0.00095172s 4.01878 %
M2M                         0.00311584s 13.1571 %
M2L                         0.00898328s 37.9332 %
L2L                         0.00339084s 14.3183 %
others (direct integration) 0.00642484s 27.1298 %
total                       0.0236818 s
--------------------------------------------------
Elapsed time for GMRES: {2.61369,2.61369} [s] for size = 100
GMRES expected error = 1.60452e-07 True error = 1.60452e-07
destructing gmres
Restart count :2
----- Average time for 50 iterations -----
update poles                0.00046278s 1.96222 %
reset Moments               0.0003399s 1.4412 %
increment_M_reuse           0.00095314s 4.04137 %
M2M                         0.00310956s 13.1847 %
M2L                         0.00899804s 38.1522 %
L2L                         0.00338396s 14.3482 %
others (direct integration) 0.00633718s 26.87 %
total                       0.0235846 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00046278s 1.97096 %
reset Moments               0.00034516s 1.47002 %
increment_M_reuse           0.0009555s 4.06943 %
M2M                         0.00311826s 13.2805 %
M2L                         0.00893078s 38.0358 %
L2L                         0.00338778s 14.4284 %
others (direct integration) 0.00627968s 26.7449 %
total                       0.0234799 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00046264s 1.97259 %
reset Moments               0.00035072s 1.49539 %
increment_M_reuse           0.00095942s 4.09074 %
M2M                         0.00310688s 13.247 %
M2L                         0.00888798s 37.8962 %
L2L                         0.00338102s 14.4159 %
others (direct integration) 0.00630482s 26.8822 %
total                       0.0234535 s
--------------------------------------------------
Elapsed time for GMRES: {4.01224,4.01224} [s] for size = 150
GMRES expected error = 1.10035e-13 True error = 2.16946e-12
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {0.000623322,1.60452e-07,1.10035e-13}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 13.5607 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {13.5607,26.1386} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.022131,26.1607} s
name = tank
net->velocityTranslational() = {updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0,0.0047955,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
0,0}
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
Elapsed time: {0.740557,26.9013} s
name = tankupdating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0047955,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity

net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.377296,27.2786} s
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
Elapsed time: {3.8e-05,27.2786} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 31039
Creating sources on surfaces...
バケットの作成．極の追加
added :52312
バケットの作成．極の追加, Elapsed time : {0.54367,0.54367}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 392 deepest level buckets.
After shrink: 392 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.109003,0.652673}
initializeFMM
FourierM2L_double
[FMM:init] Step0: update poles start, #poles=52312
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=52312
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.011775,0.011775}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.000766,0.000766}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000166,0.000932}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.029959,0.030891}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.050429,0.08132}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.102438,0.183758}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {3.53027,3.71402}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000337,3.71436}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000319,3.71468}
setM2L, Elapsed time : {3.7147,3.7147}
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
setL2P, Elapsed time : {0.006742,0.006742} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=26158)
setDirectIntegration, Elapsed time : {0.551611,0.558353} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=997.897 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.00075,0.00075}
reset Moments, Elapsed time : {0.00039,0.00114}
increment_M_reuse, Elapsed time : {0.000945,0.002085}
M2M, Elapsed time : {0.003206,0.005291}
M2L, Elapsed time : {0.008554,0.013845}
L2L, Elapsed time : {0.003601,0.017446}
updateFMM, Elapsed time : {1e-05,0.017456}
integration (direct integration is dominant), Elapsed time : {0.006225,0.006225}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00047436s 2.01798 %
reset Moments               0.0003446s 1.46597 %
increment_M_reuse           0.0009512s 4.04651 %
M2M                         0.00313568s 13.3395 %
M2L                         0.00892104s 37.9511 %
L2L                         0.00338662s 14.4071 %
others (direct integration) 0.00629316s 26.7718 %
total                       0.0235067 s
--------------------------------------------------
Elapsed time for GMRES: {1.2667,1.2667} [s] for size = 50
GMRES expected error = 5.79693e-09 True error = 5.79693e-09
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {5.79693e-09}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 6.85035 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {6.85037,34.129} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.021974,34.151} s
name = tank
net->velocityTranslational() = {updating wavemaker's (RigidBodyObject) velocity0,0,0}
use tank's (RigidBodyObject) predetermiend velocity

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
Elapsed time: {0.837765,34.9887} s
name = tank
net->velocityTranslational() = {updating 0,0,0}
wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0047955,0,0}
use tank's (RigidBodyObject) predetermiend velocityuse wavemaker's (RigidBodyObject) predetermiend velocity

name:water: setBounds
Elapsed time: {0.379519,35.3683} s
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
Elapsed time: {4.1e-05,35.3683} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 31039
Creating sources on surfaces...
バケットの作成．極の追加
added :52312
バケットの作成．極の追加, Elapsed time : {0.57263,0.57263}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 392 deepest level buckets.
After shrink: 392 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.110975,0.683605}
initializeFMM
FourierM2L_double
[FMM:init] Step0: update poles start, #poles=52312
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=52312
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.011735,0.011735}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.000747,0.000747}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000149,0.000896}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.029729,0.030625}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.052375,0.083}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.106496,0.189496}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {3.58542,3.77492}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000296,3.77522}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000306,3.77552}
setM2L, Elapsed time : {3.77553,3.77553}
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
setL2P, Elapsed time : {0.007058,0.007058} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=26158)
setDirectIntegration, Elapsed time : {0.557642,0.5647} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=997.897 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000752,0.000752}
reset Moments, Elapsed time : {0.000404,0.001156}
increment_M_reuse, Elapsed time : {0.000966,0.002122}
M2M, Elapsed time : {0.003235,0.005357}
M2L, Elapsed time : {0.008258,0.013615}
L2L, Elapsed time : {0.003462,0.017077}
updateFMM, Elapsed time : {3e-06,0.01708}
integration (direct integration is dominant), Elapsed time : {0.006518,0.006518}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00048088s 2.01453 %
reset Moments               0.00034578s 1.44856 %
increment_M_reuse           0.0009749s 4.08411 %
M2M                         0.00314426s 13.1721 %
M2L                         0.00904778s 37.9035 %
L2L                         0.00340864s 14.2797 %
others (direct integration) 0.00646834s 27.0975 %
total                       0.0238706 s
--------------------------------------------------
Elapsed time for GMRES: {1.28584,1.28584} [s] for size = 50
GMRES expected error = 0.000623178 True error = 0.000623178
destructing gmres
Restart count :1
----- Average time for 50 iterations -----
update poles                0.00047524s 2.01523 %
reset Moments               0.00034784s 1.475 %
increment_M_reuse           0.0009437s 4.00171 %
M2M                         0.00313982s 13.3143 %
M2L                         0.00893438s 37.8858 %
L2L                         0.00339596s 14.4004 %
others (direct integration) 0.00634546s 26.9076 %
total                       0.0235824 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00047154s 1.99571 %
reset Moments               0.00034794s 1.4726 %
increment_M_reuse           0.00094824s 4.01327 %
M2M                         0.00313084s 13.2508 %
M2L                         0.0089262s 37.7786 %
L2L                         0.003388s 14.3391 %
others (direct integration) 0.00641488s 27.1499 %
total                       0.0236276 s
--------------------------------------------------
Elapsed time for GMRES: {2.61124,2.61124} [s] for size = 100
GMRES expected error = 1.60483e-07 True error = 1.60483e-07
destructing gmres
Restart count :2
----- Average time for 50 iterations -----
update poles                0.00047308s 2.00169 %
reset Moments               0.00034566s 1.46255 %
increment_M_reuse           0.00094428s 3.99542 %
M2M                         0.00313446s 13.2625 %
M2L                         0.00895252s 37.8797 %
L2L                         0.0033891s 14.3399 %
others (direct integration) 0.00639498s 27.0583 %
total                       0.0236341 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00047184s 1.97644 %
reset Moments               0.00034706s 1.45377 %
increment_M_reuse           0.0009545s 3.99821 %
M2M                         0.0031235s 13.0837 %
M2L                         0.0090334s 37.8391 %
L2L                         0.0033882s 14.1925 %
others (direct integration) 0.00655468s 27.4563 %
total                       0.0238732 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00046958s 1.985 %
reset Moments               0.00035368s 1.49507 %
increment_M_reuse           0.00095592s 4.04084 %
M2M                         0.00311152s 13.153 %
M2L                         0.00894458s 37.8103 %
L2L                         0.00338216s 14.297 %
others (direct integration) 0.006439s 27.2188 %
total                       0.0236564 s
--------------------------------------------------
Elapsed time for GMRES: {4.05852,4.05852} [s] for size = 150
GMRES expected error = 1.09735e-13 True error = 3.88981e-12
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {0.000623178,1.60483e-07,1.09735e-13}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 13.7057 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {13.7057,49.074} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.022,49.096} s
ALEのU_update_BEMを計算
Elapsed time: 1.87098 s
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
Elapsed time: {0.837097,51.8041} s
name = tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.00958981,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.341647,52.1458} s
Total elapsed time: 52.1458 s
============================================================================
ElapsedTimeBIEDiscretization={9.97456,13.5227,6.81259,13.6673}
ElapsedTimeSolve={0,0,0,0}
ElapsedTimeALE=1.87098
ElapsedTimeTotal=52.1458
============================================================================
simulation_timeを取得
/Users/tomoaki/Ruehl20160d075_FourierM2L_double/output/water_0.vtu  VV_points.size() : 52312 |||||||||
/Users/tomoaki/Ruehl20160d075_FourierM2L_double/output/water.pvd
Creating /Users/tomoaki/Ruehl20160d075_FourierM2L_double/output/water.pvd ...
Done.
/Users/tomoaki/Ruehl20160d075_FourierM2L_double/output/tank_0.vtu  VV_points.size() : 52 |||||||||
/Users/tomoaki/Ruehl20160d075_FourierM2L_double/output/tank.pvd
Creating /Users/tomoaki/Ruehl20160d075_FourierM2L_double/output/tank.pvd ...
Done.
/Users/tomoaki/Ruehl20160d075_FourierM2L_double/output/wavemaker_0.vtu  VV_points.size() : 496 |||||||||
/Users/tomoaki/Ruehl20160d075_FourierM2L_double/output/wavemaker.pvd
Creating /Users/tomoaki/Ruehl20160d075_FourierM2L_double/output/wavemaker.pvd ...
Done.
flipIf
flipIf
flipIf
net.getPoints() = 26158
Total : 26158
Total variables: 29013
Total case double-node : 29013
node reduction : 0.0984042
CORNER : 561
Total CORNER faces : 3416
Neumann : 18554
Dirichlet : 7043
===========================================================================
       dt :0.01
time_step :1
real time :0.01
---------------------------------------------------------------------------
RK_step = 1/4, RK_time = 0.01, simulation_time = 0.01
water makeBucketPoints(2.01457)
tank makeBucketPoints(wavemaker makeBucketPoints(2.99541)
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
Elapsed time: {0.851863,0.851863} s
waterの境界条件を決定 setBoundaryTypes
water setContactFaces()
step2 面の境界条件を判定
step3 線の境界条件を決定
step4 点の境界条件を決定
setBoundaryTypes終了
setBoundaryTypes
Elapsed time: {0.6151,1.46696} s
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
Elapsed time: {4.4e-05,1.46701} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 31039
Setting up buckets...
Buckets set up complete.
Creating sources on surfaces...
バケットの作成．極の追加
added :52312
バケットの作成．極の追加, Elapsed time : {0.011327,0.011327}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 1 deepest level buckets.
After shrink: 1 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {3.04405,3.05537}
initializeFMM
FourierM2L_double
[FMM:init] Step0: update poles start, #poles=52312
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=52312
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.010092,0.010092}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.001126,0.001126}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.00019,0.001316}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.029806,0.031122}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.055169,0.086291}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.106896,0.193187}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {3.56854,3.76173}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000189,3.76192}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000234,3.76215}
setM2L, Elapsed time : {3.76216,3.76216}
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
setL2P, Elapsed time : {0.006942,0.006942} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=26158)
setDirectIntegration, Elapsed time : {0.591805,0.598747} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=997.298 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000757,0.000757}
reset Moments, Elapsed time : {0.000461,0.001218}
increment_M_reuse, Elapsed time : {0.001604,0.002822}
M2M, Elapsed time : {0.004202,0.007024}
M2L, Elapsed time : {0.010456,0.01748}
L2L, Elapsed time : {0.00519,0.02267}
updateFMM, Elapsed time : {1.5e-05,0.022685}
integration (direct integration is dominant), Elapsed time : {0.007036,0.007036}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00049338s 1.99035 %
reset Moments               0.0003701s 1.49303 %
increment_M_reuse           0.00099928s 4.03121 %
M2M                         0.00315774s 12.7387 %
M2L                         0.00929262s 37.4875 %
L2L                         0.00345756s 13.9482 %
others (direct integration) 0.0070179s 28.311 %
total                       0.0247886 s
--------------------------------------------------
Elapsed time for GMRES: {1.3335,1.3335} [s] for size = 50
GMRES expected error = 6.86257e-06 True error = 6.86257e-06
destructing gmres
Restart count :1
----- Average time for 50 iterations -----
update poles                0.00048738s 1.84165 %
reset Moments               0.0003702s 1.39887 %
increment_M_reuse           0.00125432s 4.73968 %
M2M                         0.00317972s 12.0151 %
M2L                         0.0103336s 39.0475 %
L2L                         0.00344932s 13.0339 %
others (direct integration) 0.00738968s 27.9232 %
total                       0.0264643 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00047596s 1.92895 %
reset Moments               0.00036486s 1.47869 %
increment_M_reuse           0.0009874s 4.00169 %
M2M                         0.00314408s 12.7422 %
M2L                         0.00925438s 37.5057 %
L2L                         0.00343962s 13.9399 %
others (direct integration) 0.00700828s 28.4028 %
total                       0.0246746 s
--------------------------------------------------
Elapsed time for GMRES: {2.82971,2.82971} [s] for size = 100
GMRES expected error = 1.49832e-09 True error = 1.49832e-09
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {6.86257e-06,1.49832e-09}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 12.2298 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {12.2298,13.6968} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.021647,13.7184} s
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {name = tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
0.00958981,0,0}
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
Elapsed time: {0.830262,14.5487} s
name = tank
net->velocityTranslational() = {updating wavemaker's (RigidBodyObject) velocity0,0,0}
name = wavemaker
net->velocityTranslational() = {0.00958981,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity

use tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.388284,14.937} s
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
Elapsed time: {4.6e-05,14.937} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 31039
Creating sources on surfaces...
バケットの作成．極の追加
added :52312
バケットの作成．極の追加, Elapsed time : {0.582583,0.582583}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 392 deepest level buckets.
After shrink: 392 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.112383,0.694966}
initializeFMM
FourierM2L_double
[FMM:init] Step0: update poles start, #poles=52312
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=52312
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.011993,0.011993}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.000747,0.000747}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000166,0.000913}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.029948,0.030861}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.055161,0.086022}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.108585,0.194607}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {3.55186,3.74647}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000196,3.74667}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000209,3.74688}
setM2L, Elapsed time : {3.74689,3.74689}
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
setL2P, Elapsed time : {0.006706,0.006706} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=26158)
setDirectIntegration, Elapsed time : {0.546937,0.553643} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=997.298 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000762,0.000762}
reset Moments, Elapsed time : {0.000425,0.001187}
increment_M_reuse, Elapsed time : {0.000937,0.002124}
M2M, Elapsed time : {0.003294,0.005418}
M2L, Elapsed time : {0.008511,0.013929}
L2L, Elapsed time : {0.003524,0.017453}
updateFMM, Elapsed time : {5e-06,0.017458}
integration (direct integration is dominant), Elapsed time : {0.006836,0.006836}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00048554s 1.98037 %
reset Moments               0.00036528s 1.48986 %
increment_M_reuse           0.00098586s 4.02101 %
M2M                         0.0031623s 12.898 %
M2L                         0.0092831s 37.8629 %
L2L                         0.00345102s 14.0756 %
others (direct integration) 0.0067846s 27.6723 %
total                       0.0245177 s
--------------------------------------------------
Elapsed time for GMRES: {1.31947,1.31947} [s] for size = 50
GMRES expected error = 0.000583678 True error = 0.000583678
destructing gmres
Restart count :1
----- Average time for 50 iterations -----
update poles                0.00048874s 1.98671 %
reset Moments               0.00036954s 1.50217 %
increment_M_reuse           0.00098184s 3.99115 %
M2M                         0.00316224s 12.8544 %
M2L                         0.00931364s 37.8596 %
L2L                         0.00344496s 14.0037 %
others (direct integration) 0.00683948s 27.8023 %
total                       0.0246004 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.0004812s 1.95854 %
reset Moments               0.0003655s 1.48763 %
increment_M_reuse           0.00099236s 4.03902 %
M2M                         0.00315092s 12.8246 %
M2L                         0.00930928s 37.8898 %
L2L                         0.00342482s 13.9394 %
others (direct integration) 0.00684526s 27.861 %
total                       0.0245693 s
--------------------------------------------------
Elapsed time for GMRES: {2.71556,2.71556} [s] for size = 100
GMRES expected error = 3.62242e-07 True error = 3.62242e-07
destructing gmres
Restart count :2
