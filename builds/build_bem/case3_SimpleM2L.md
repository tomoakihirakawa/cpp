ケース：case3
14840節点
方法:SimpleM2L

初回（1step目のRKの1回目）は収束が早いので，2回目以降の計算を比較すること．

------------ LOG -----------

Last login: Thu Sep 18 13:25:00 on ttys041
/Users/tomoaki/.zshrc is loaded
build_bem(main)$ ./SimpleM2L /Users/tomoaki/Ruehl20160d1_SimpleM2L
input_directory : "/Users/tomoaki/Ruehl20160d1_SimpleM2L"
ALE: {pseudo_quad}
ALEPERIOD: {1}
GRAVITY: {9.81}
WATER_DENSITY: {1000.0}
element: {linear}
end_time: {10.0}
end_time_step: {100000}
input_files: {tank.json,wavemaker.json,water.json,wg3.json,wg6.json}
max_dt: {0.1}
output_directory: {/Users/tomoaki/Ruehl20160d1_SimpleM2L/output}
LINEAR_ELEMENT
ALE_ON_PSEUDO_QUADRATIC_ELEMENT

WATER_DENSITY: 1000
GRAVITY: 9.81
"/Users/tomoaki/Ruehl20160d1_SimpleM2L/tank.json"
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
/Users/tomoaki/Ruehl20160d1_SimpleM2L/output/tank_init.vtu  VV_points.size() : 52 |||||||||
"/Users/tomoaki/Ruehl20160d1_SimpleM2L/wavemaker.json"
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
        this->getName(): wavemaker, 0x150030000
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
/Users/tomoaki/Ruehl20160d1_SimpleM2L/output/wavemaker_init.vtu  VV_points.size() : 496 |||||||||
"/Users/tomoaki/Ruehl20160d1_SimpleM2L/water.json"
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
        this->getName(): water, 0x140038000
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
/Users/tomoaki/Ruehl20160d1_SimpleM2L/output/water_init.vtu  VV_points.size() : 29676 |||||||||
"/Users/tomoaki/Ruehl20160d1_SimpleM2L/wg3.json"
name: {wg3}
position: {9.48,0.0,1.86,9.48,0.0,0.8600000000000001}
type: {wavegauge}
type = wavegauge
skipped
"/Users/tomoaki/Ruehl20160d1_SimpleM2L/wg6.json"
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
water makeBucketPoints(2.01457)
wavemaker makeBucketPoints(0.616441)
tank makeBucketPoints(2.99541)
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
Elapsed time: {0.504846,0.504846} s
waterの境界条件を決定 setBoundaryTypes
water setContactFaces()
step2 面の境界条件を判定
step3 線の境界条件を決定
step4 点の境界条件を決定
setBoundaryTypes終了
setBoundaryTypes
Elapsed time: {0.311369,0.816215} s
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
Elapsed time: {4.9e-05,0.816264} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18371
Setting up buckets...
Buckets set up complete.
Creating sources on surfaces...
バケットの作成．極の追加
added :29676
バケットの作成．極の追加, Elapsed time : {0.005588,0.005588}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 1 deepest level buckets.
After shrink: 1 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {1.59277,1.59836}
initializeFMM
SimpleM2L
[FMM:init] Step0: update poles start, #poles=29676
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=29676
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.010344,0.010344}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.001344,0.001344}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000151,0.001495}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.005878,0.007373}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.01369,0.021063}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.024578,0.045641}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {0.783712,0.829353}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.0002,0.829553}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000212,0.829765}
setM2L, Elapsed time : {0.829773,0.829773}
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
setL2P, Elapsed time : {0.004513,0.004513} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=14840)
setDirectIntegration, Elapsed time : {0.174552,0.179065} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=602.773 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000385,0.000385}
reset Moments, Elapsed time : {0.000357,0.000742}
increment_M_reuse, Elapsed time : {0.001203,0.001945}
M2M, Elapsed time : {0.003639,0.005584}
M2L, Elapsed time : {0.04557,0.051154}
L2L, Elapsed time : {0.004836,0.05599}
updateFMM, Elapsed time : {5e-06,0.055995}
integration (direct integration is dominant), Elapsed time : {0.002647,0.002647}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00025342s 0.365586 %
reset Moments               0.0003209s 0.462933 %
increment_M_reuse           0.00062518s 0.90189 %
M2M                         0.00411128s 5.93097 %
M2L                         0.0553483s 79.8459 %
L2L                         0.00600544s 8.6635 %
others (direct integration) 0.00265436s 3.8292 %
total                       0.0693188 s
--------------------------------------------------
Elapsed time for GMRES: {3.6134,3.6134} [s] for size = 50
GMRES expected error = 0 True error = 0
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {0}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 6.65412 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {6.65414,7.4704} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.013492,7.4839} s
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
Elapsed time: {0.374042,7.85794} s
name = tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.209051,8.06699} s
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
Elapsed time: {4.3e-05,8.06703} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18371
Creating sources on surfaces...
バケットの作成．極の追加
added :29676
バケットの作成．極の追加, Elapsed time : {0.325915,0.325915}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 392 deepest level buckets.
After shrink: 392 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.107473,0.433388}
initializeFMM
SimpleM2L
[FMM:init] Step0: update poles start, #poles=29676
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=29676
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.011783,0.011783}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.000746,0.000746}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000174,0.00092}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.005669,0.006589}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.010975,0.017564}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.022551,0.040115}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {0.735056,0.775171}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000261,0.775432}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000261,0.775693}
setM2L, Elapsed time : {0.775704,0.775704}
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
setL2P, Elapsed time : {0.004472,0.004472} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=14840)
setDirectIntegration, Elapsed time : {0.179727,0.184199} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=602.773 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000416,0.000416}
reset Moments, Elapsed time : {0.000436,0.000852}
increment_M_reuse, Elapsed time : {0.00061,0.001462}
M2M, Elapsed time : {0.003222,0.004684}
M2L, Elapsed time : {0.023784,0.028468}
L2L, Elapsed time : {0.003747,0.032215}
updateFMM, Elapsed time : {5e-06,0.03222}
integration (direct integration is dominant), Elapsed time : {0.002171,0.002171}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00025502s 0.780154 %
reset Moments               0.0003041s 0.930299 %
increment_M_reuse           0.00061144s 1.87051 %
M2M                         0.00312276s 9.55311 %
M2L                         0.0228608s 69.9355 %
L2L                         0.00337498s 10.3247 %
others (direct integration) 0.00215932s 6.60577 %
total                       0.0326884 s
--------------------------------------------------
Elapsed time for GMRES: {1.70966,1.70966} [s] for size = 50
GMRES expected error = 0.00032394 True error = 0.00032394
destructing gmres
Restart count :1
----- Average time for 50 iterations -----
update poles                0.00025328s 0.773675 %
reset Moments               0.00030374s 0.927812 %
increment_M_reuse           0.00061762s 1.8866 %
M2M                         0.0031444s 9.60496 %
M2L                         0.0228601s 69.829 %
L2L                         0.00337532s 10.3103 %
others (direct integration) 0.00218278s 6.66757 %
total                       0.0327372 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00025266s 0.770852 %
reset Moments               0.00030384s 0.926999 %
increment_M_reuse           0.00062216s 1.89818 %
M2M                         0.00314934s 9.60847 %
M2L                         0.0228928s 69.8447 %
L2L                         0.00338056s 10.3139 %
others (direct integration) 0.00217536s 6.63691 %
total                       0.0327767 s
--------------------------------------------------
Elapsed time for GMRES: {3.44412,3.44412} [s] for size = 100
GMRES expected error = 1.55877e-08 True error = 1.55877e-08
destructing gmres
Restart count :2
----- Average time for 50 iterations -----
update poles                0.00025394s 0.77888 %
reset Moments               0.00030508s 0.935736 %
increment_M_reuse           0.00061634s 1.89043 %
M2M                         0.00314424s 9.64396 %
M2L                         0.0227323s 69.7242 %
L2L                         0.00338282s 10.3757 %
others (direct integration) 0.00216848s 6.65112 %
total                       0.0326032 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.0002528s 0.765081 %
reset Moments               0.00030668s 0.928145 %
increment_M_reuse           0.0006208s 1.87881 %
M2M                         0.00313674s 9.49312 %
M2L                         0.0231801s 70.1531 %
L2L                         0.00337458s 10.2129 %
others (direct integration) 0.0021705s 6.56886 %
total                       0.0330422 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00025344s 0.778994 %
reset Moments               0.00030688s 0.943252 %
increment_M_reuse           0.00062118s 1.90931 %
M2M                         0.00314312s 9.66095 %
M2L                         0.0226323s 69.5645 %
L2L                         0.00337954s 10.3876 %
others (direct integration) 0.00219782s 6.7554 %
total                       0.0325343 s
--------------------------------------------------
Elapsed time for GMRES: {5.23132,5.23132} [s] for size = 150
GMRES expected error = 4.67857e-16 True error = 7.16066e-15
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {0.00032394,1.55877e-08,4.67857e-16}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 12.2516 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {12.2517,20.3187} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.013046,20.3317} s
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
Elapsed time: {0.461582,20.7933} s
name = tank
net->velocityTranslational() = {updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0047955,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.208148,21.0015} s
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
Elapsed time: {4.9e-05,21.0015} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18371
Creating sources on surfaces...
バケットの作成．極の追加
added :29676
バケットの作成．極の追加, Elapsed time : {0.333561,0.333561}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 392 deepest level buckets.
After shrink: 392 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.10481,0.438371}
initializeFMM
SimpleM2L
[FMM:init] Step0: update poles start, #poles=29676
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=29676
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.012168,0.012168}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.000711,0.000711}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000195,0.000906}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.005568,0.006474}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.011851,0.018325}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.021323,0.039648}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {0.728805,0.768453}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.00024,0.768693}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000161,0.768854}
setM2L, Elapsed time : {0.768868,0.768868}
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
setL2P, Elapsed time : {0.003868,0.003868} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=14840)
setDirectIntegration, Elapsed time : {0.168333,0.172201} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=602.773 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000394,0.000394}
reset Moments, Elapsed time : {0.000361,0.000755}
increment_M_reuse, Elapsed time : {0.000585,0.00134}
M2M, Elapsed time : {0.003155,0.004495}
M2L, Elapsed time : {0.022096,0.026591}
L2L, Elapsed time : {0.00361,0.030201}
updateFMM, Elapsed time : {5e-06,0.030206}
integration (direct integration is dominant), Elapsed time : {0.002173,0.002173}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00025664s 0.780544 %
reset Moments               0.00030308s 0.921786 %
increment_M_reuse           0.00061718s 1.87709 %
M2M                         0.00314248s 9.55753 %
M2L                         0.0230255s 70.0298 %
L2L                         0.00338266s 10.288 %
others (direct integration) 0.00215206s 6.54527 %
total                       0.0328796 s
--------------------------------------------------
Elapsed time for GMRES: {1.71812,1.71812} [s] for size = 50
GMRES expected error = 3.26797e-09 True error = 3.26797e-09
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {3.26797e-09}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 3.53663 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {3.53665,24.5382} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.011852,24.55} s
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
Elapsed time: {0.475169,25.0252} s
name = tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0047955,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.211929,25.2371} s
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
Elapsed time: {4e-05,25.2371} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18371
Creating sources on surfaces...
バケットの作成．極の追加
added :29676
バケットの作成．極の追加, Elapsed time : {0.343261,0.343261}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 392 deepest level buckets.
After shrink: 392 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.108233,0.451494}
initializeFMM
SimpleM2L
[FMM:init] Step0: update poles start, #poles=29676
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=29676
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.011969,0.011969}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.000705,0.000705}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.0002,0.000905}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.00553,0.006435}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.010555,0.01699}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.023258,0.040248}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {0.728957,0.769205}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000216,0.769421}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000174,0.769595}
setM2L, Elapsed time : {0.769601,0.769601}
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
setL2P, Elapsed time : {0.003917,0.003917} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=14840)
setDirectIntegration, Elapsed time : {0.163421,0.167338} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=602.773 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000392,0.000392}
reset Moments, Elapsed time : {0.000349,0.000741}
increment_M_reuse, Elapsed time : {0.00058,0.001321}
M2M, Elapsed time : {0.00319,0.004511}
M2L, Elapsed time : {0.022631,0.027142}
L2L, Elapsed time : {0.00343,0.030572}
updateFMM, Elapsed time : {5e-06,0.030577}
integration (direct integration is dominant), Elapsed time : {0.002041,0.002041}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00025586s 0.78771 %
reset Moments               0.00030414s 0.936348 %
increment_M_reuse           0.00061538s 1.89456 %
M2M                         0.00314986s 9.6974 %
M2L                         0.0226368s 69.6913 %
L2L                         0.00337768s 10.3988 %
others (direct integration) 0.0021418s 6.59391 %
total                       0.0324815 s
--------------------------------------------------
Elapsed time for GMRES: {1.69532,1.69532} [s] for size = 50
GMRES expected error = 0.000323894 True error = 0.000323894
destructing gmres
Restart count :1
----- Average time for 50 iterations -----
update poles                0.00025422s 0.775529 %
reset Moments               0.00030266s 0.923301 %
increment_M_reuse           0.00061914s 1.88876 %
M2M                         0.00313888s 9.57554 %
M2L                         0.022952s 70.0178 %
L2L                         0.00337056s 10.2823 %
others (direct integration) 0.00214276s 6.53675 %
total                       0.0327802 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.0002544s 0.767186 %
reset Moments               0.00030594s 0.922613 %
increment_M_reuse           0.0006146s 1.85343 %
M2M                         0.00313112s 9.44242 %
M2L                         0.0233488s 70.4121 %
L2L                         0.00336952s 10.1614 %
others (direct integration) 0.0021358s 6.44086 %
total                       0.0331602 s
--------------------------------------------------
Elapsed time for GMRES: {3.46833,3.46833} [s] for size = 100
GMRES expected error = 1.55862e-08 True error = 1.55862e-08
destructing gmres
Restart count :2
----- Average time for 50 iterations -----
update poles                0.00025516s 0.775431 %
reset Moments               0.00030868s 0.938078 %
increment_M_reuse           0.00061282s 1.86236 %
M2M                         0.00314772s 9.56592 %
M2L                         0.0230176s 69.9504 %
L2L                         0.00340392s 10.3445 %
others (direct integration) 0.00215968s 6.56326 %
total                       0.0329056 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.0002558s 0.777869 %
reset Moments               0.00031256s 0.950472 %
increment_M_reuse           0.00061966s 1.88434 %
M2M                         0.00314956s 9.57758 %
M2L                         0.0229924s 69.9182 %
L2L                         0.00339014s 10.3092 %
others (direct integration) 0.00216458s 6.58233 %
total                       0.0328847 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.00025414s 0.780934 %
reset Moments               0.00030718s 0.943917 %
increment_M_reuse           0.00062258s 1.91309 %
M2M                         0.003144s 9.66103 %
M2L                         0.0226985s 69.7489 %
L2L                         0.00337662s 10.3758 %
others (direct integration) 0.00214012s 6.57626 %
total                       0.0325431 s
--------------------------------------------------
Elapsed time for GMRES: {5.23894,5.23894} [s] for size = 150
GMRES expected error = 4.6771e-16 True error = 1.41288e-14
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {0.000323894,1.55862e-08,4.6771e-16}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 12.3103 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {12.3103,37.5474} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.012591,37.56} s
ALEのU_update_BEMを計算
Elapsed time: 1.0592 s
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
Elapsed time: {0.415131,39.0344} s
name = tank
net->velocityTranslational() = {updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.00958981,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.202525,39.2369} s
Total elapsed time: 39.2369 s
============================================================================
ElapsedTimeBIEDiscretization={6.63385,12.2309,3.51644,12.29}
ElapsedTimeSolve={0,0,0,0}
ElapsedTimeALE=1.0592
ElapsedTimeTotal=39.2369
============================================================================
simulation_timeを取得
/Users/tomoaki/Ruehl20160d1_SimpleM2L/output/water_0.vtu  VV_points.size() : 29676 |||||||||
/Users/tomoaki/Ruehl20160d1_SimpleM2L/output/water.pvd
Creating /Users/tomoaki/Ruehl20160d1_SimpleM2L/output/water.pvd ...
Done.
/Users/tomoaki/Ruehl20160d1_SimpleM2L/output/tank_0.vtu  VV_points.size() : 52 |||||||||
/Users/tomoaki/Ruehl20160d1_SimpleM2L/output/tank.pvd
Creating /Users/tomoaki/Ruehl20160d1_SimpleM2L/output/tank.pvd ...
Done.
/Users/tomoaki/Ruehl20160d1_SimpleM2L/output/wavemaker_0.vtu  VV_points.size() : 496 |||||||||
/Users/tomoaki/Ruehl20160d1_SimpleM2L/output/wavemaker.pvd
Creating /Users/tomoaki/Ruehl20160d1_SimpleM2L/output/wavemaker.pvd ...
Done.
flipIf
flipIf
flipIf
net.getPoints() = 14840
Total : 14840
Total variables: 16913
Total case double-node : 16913
node reduction : 0.122568
CORNER : 402
Total CORNER faces : 2475
Neumann : 10391
Dirichlet : 4047
===========================================================================
       dt :0.01
time_step :1
real time :0.01
---------------------------------------------------------------------------
RK_step = 1/4, RK_time = 0.01, simulation_time = 0.01
water makeBucketPoints(tank makeBucketPoints(wavemaker makeBucketPoints(2.99541)
0.616441)2.01457)

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
Elapsed time: {0.495122,0.495122} s
waterの境界条件を決定 setBoundaryTypes
water setContactFaces()
step2 面の境界条件を判定
step3 線の境界条件を決定
step4 点の境界条件を決定
setBoundaryTypes終了
setBoundaryTypes
Elapsed time: {0.308611,0.803733} s
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
Elapsed time: {3.8e-05,0.803771} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18371
Setting up buckets...
Buckets set up complete.
Creating sources on surfaces...
バケットの作成．極の追加
added :29676
バケットの作成．極の追加, Elapsed time : {0.005452,0.005452}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 1 deepest level buckets.
After shrink: 1 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {1.57056,1.57601}
initializeFMM
SimpleM2L
[FMM:init] Step0: update poles start, #poles=29676
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=29676
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.010193,0.010193}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.00099,0.00099}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000206,0.001196}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.005654,0.00685}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.010645,0.017495}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.02166,0.039155}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {0.757975,0.79713}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000263,0.797393}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.000128,0.797521}
setM2L, Elapsed time : {0.797524,0.797524}
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
setL2P, Elapsed time : {0.003965,0.003965} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=14840)
setDirectIntegration, Elapsed time : {0.163754,0.167719} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=602.737 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000389,0.000389}
reset Moments, Elapsed time : {0.000376,0.000765}
increment_M_reuse, Elapsed time : {0.000583,0.001348}
M2M, Elapsed time : {0.00322,0.004568}
M2L, Elapsed time : {0.023167,0.027735}
L2L, Elapsed time : {0.003668,0.031403}
updateFMM, Elapsed time : {3e-06,0.031406}
integration (direct integration is dominant), Elapsed time : {0.002376,0.002376}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00026046s 0.790138 %
reset Moments               0.00030482s 0.92471 %
increment_M_reuse           0.00061124s 1.85427 %
M2M                         0.00314344s 9.53602 %
M2L                         0.0227212s 68.9277 %
L2L                         0.0033756s 10.2403 %
others (direct integration) 0.00254706s 7.72683 %
total                       0.0329639 s
--------------------------------------------------
Elapsed time for GMRES: {1.72406,1.72406} [s] for size = 50
GMRES expected error = 2.27281e-06 True error = 2.27281e-06
destructing gmres
Restart count :1
----- Average time for 50 iterations -----
update poles                0.00025778s 0.783564 %
reset Moments               0.0003067s 0.932264 %
increment_M_reuse           0.00061738s 1.87663 %
M2M                         0.00315474s 9.58934 %
M2L                         0.0226309s 68.7904 %
L2L                         0.0033737s 10.2549 %
others (direct integration) 0.00255716s 7.7729 %
total                       0.0328984 s
--------------------------------------------------
----- Average time for 50 iterations -----
update poles                0.0002574s 0.774727 %
reset Moments               0.00031002s 0.933103 %
increment_M_reuse           0.00061992s 1.86585 %
M2M                         0.00315312s 9.49031 %
M2L                         0.0229893s 69.1936 %
L2L                         0.00338252s 10.1808 %
others (direct integration) 0.00251232s 7.56162 %
total                       0.0332246 s
--------------------------------------------------
Elapsed time for GMRES: {3.47481,3.47481} [s] for size = 100
GMRES expected error = 7.21819e-11 True error = 7.21819e-11
destructing gmres
gmres size list = {50,100,150,200}
gmres error = {2.27281e-06,7.21819e-11}
update p->phiphin and p->phinOnFace for Dirichlet boundary
Elapsed time for solving BIE: 8.18031 s
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {8.18033,8.9841} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.012759,8.99686} s
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
Elapsed time: {0.439728,9.43659} s
name = tank
net->velocityTranslational() = {0,updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.00958981,0,0}
0,0}
use tank's (RigidBodyObject) predetermiend velocity
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.228468,9.66506} s
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
Elapsed time: {4e-05,9.6651} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18371
Creating sources on surfaces...
バケットの作成．極の追加
added :29676
バケットの作成．極の追加, Elapsed time : {0.318609,0.318609}
ツリー構造を生成，またはrebin
Escaped objects after rebin: 0
Success objects after rebin: 0
Remaining objects after rebin: 0
Before shrink: 392 deepest level buckets.
After shrink: 392 deepest level buckets.
After grow: 392 deepest level buckets.
ツリー構造を生成，またはrebin, Elapsed time : {0.107282,0.425891}
initializeFMM
SimpleM2L
[FMM:init] Step0: update poles start, #poles=29676
[FMM:init] Step0: update poles done
[FMM:init] Step1: reset moments for total nodes=513, levels=8
[FMM:init] Step1: reset moments done
[FMM:init] Step2: P2M begin, #leaves=392, total leaf sources=29676
[FMM:init] Step2: P2M (increment_M) done
[FMM:init] Step3: setM2M start
setM2M, Elapsed time : {0.01192,0.01192}
[FMM:init] Step3: setM2M done, setM2L start
setM2L, level=0, buckets_at_a_level.size()=1, Elapsed time : {0.000695,0.000695}
setM2L, level=1, buckets_at_a_level.size()=8, Elapsed time : {0.000155,0.00085}
setM2L, level=2, buckets_at_a_level.size()=16, Elapsed time : {0.005618,0.006468}
setM2L, level=3, buckets_at_a_level.size()=32, Elapsed time : {0.010925,0.017393}
setM2L, level=4, buckets_at_a_level.size()=64, Elapsed time : {0.021875,0.039268}
setM2L, level=5, buckets_at_a_level.size()=392, Elapsed time : {0.722375,0.761643}
setM2L, level=6, buckets_at_a_level.size()=0, Elapsed time : {0.000287,0.76193}
setM2L, level=7, buckets_at_a_level.size()=0, Elapsed time : {0.00028,0.76221}
setM2L, Elapsed time : {0.762219,0.762219}
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
setL2P, Elapsed time : {0.004004,0.004004} [FMM:init Step5 done]
各ターゲットの近傍ソースを保存する (#targets=14840)
setDirectIntegration, Elapsed time : {0.17349,0.177494} [FMM:init Step6 done]
mean vec_phiphin_WGNWGnN=602.737 [FMM:init Step7 done]
copy phi phin
calculate diagonal coefficient
update FMM
update poles, Elapsed time : {0.000401,0.000401}
reset Moments, Elapsed time : {0.000376,0.000777}
increment_M_reuse, Elapsed time : {0.000637,0.001414}
M2M, Elapsed time : {0.003176,0.00459}
M2L, Elapsed time : {0.022849,0.027439}
L2L, Elapsed time : {0.00344,0.030879}
updateFMM, Elapsed time : {4e-06,0.030883}
integration (direct integration is dominant), Elapsed time : {0.002601,0.002601}
Restart count :0
----- Average time for 50 iterations -----
update poles                0.00025772s 0.773141 %
reset Moments               0.000305s 0.914977 %
increment_M_reuse           0.0006153s 1.84585 %
M2M                         0.00314598s 9.43771 %
M2L                         0.0230844s 69.2515 %
L2L                         0.00338496s 10.1546 %
others (direct integration) 0.0025408s 7.62221 %
total                       0.0333342 s
--------------------------------------------------
Elapsed time for GMRES: {1.73874,1.73874} [s] for size = 50