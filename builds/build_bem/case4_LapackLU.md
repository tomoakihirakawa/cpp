ケース：case4
26158節点
方法:LapackLU

------- LOG -------

/Users/tomoaki/.zshrc is loaded
build_bem(main)$ ./main /Users/tomoaki/Ruehl20160d075                 
input_directory : "/Users/tomoaki/Ruehl20160d075"
ALE: {pseudo_quad}
ALEPERIOD: {1}
GRAVITY: {9.81}
WATER_DENSITY: {1000.0}
element: {linear}
end_time: {10.0}
end_time_step: {100000}
input_files: {tank.json,wavemaker.json,water.json,wg3.json,wg6.json}
max_dt: {0.1}
output_directory: {/Users/tomoaki/Ruehl20160d075/output}
LINEAR_ELEMENT
ALE_ON_PSEUDO_QUADRATIC_ELEMENT

WATER_DENSITY: 1000
GRAVITY: 9.81
"/Users/tomoaki/Ruehl20160d075/tank.json"
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
        this->getName(): tank, 0x138008000
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
/Users/tomoaki/Ruehl20160d075/output/tank_init.vtu  VV_points.size() : 52 |||||||||
"/Users/tomoaki/Ruehl20160d075/wavemaker.json"
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
        this->getName(): wavemaker, 0x13abd8000
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
/Users/tomoaki/Ruehl20160d075/output/wavemaker_init.vtu  VV_points.size() : 496 |||||||||
"/Users/tomoaki/Ruehl20160d075/water.json"
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
        this->getName(): water, 0x150008000
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
/Users/tomoaki/Ruehl20160d075/output/water_init.vtu  VV_points.size() : 52312 |||||||||
"/Users/tomoaki/Ruehl20160d075/wg3.json"
name: {wg3}
position: {9.48,0.0,1.86,9.48,0.0,0.8600000000000001}
type: {wavegauge}
type = wavegauge
skipped
"/Users/tomoaki/Ruehl20160d075/wg6.json"
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
water makeBucketPoints(tank makeBucketPoints(wavemaker makeBucketPoints(2.014572.99541))

0.616441)
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
Elapsed time: {1.50894,1.50894} s
waterの境界条件を決定 setBoundaryTypes
water setContactFaces()
step2 面の境界条件を判定
step3 線の境界条件を決定
step4 点の境界条件を決定
setBoundaryTypes終了
setBoundaryTypes
Elapsed time: {0.816071,2.32501} s
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
Elapsed time: {8.4e-05,2.32509} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 31040
原点を節点にとり，方程式を作成．並列化
線形要素の面の数：52312 persecent: 100
擬似二次要素の面の数：0 persecent: 0
線形要素を使ってBIEを離散化
generateBIEMatrix()
BIEの離散化に要した時間: 19.5933
LU分解を使った解法に要した時間: 14.9461
update p->phiphin and p->phinOnFace for Dirichlet boundary
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {34.6683,36.9934} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.024579,37.0179} s
name = tankupdating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity

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
Elapsed time: {5.22552,42.2435} s
name = tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.539063,42.7825} s
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
Elapsed time: {6.1e-05,42.7826} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 31040
原点を節点にとり，方程式を作成．並列化
線形要素の面の数：52312 persecent: 100
擬似二次要素の面の数：0 persecent: 0
線形要素を使ってBIEを離散化
generateBIEMatrix()
BIEの離散化に要した時間: 18.6116
LU分解を使った解法に要した時間: 13.9537
update p->phiphin and p->phinOnFace for Dirichlet boundary
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {32.7021,75.4847} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.022067,75.5067} s
name = tankupdating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0047955,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity

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
Elapsed time: {5.3816,80.8883} s
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {name = 0.0047955tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.507364,81.3957} s
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
Elapsed time: {6.4e-05,81.3958} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 31040
原点を節点にとり，方程式を作成．並列化
線形要素の面の数：52312 persecent: 100
擬似二次要素の面の数：0 persecent: 0
線形要素を使ってBIEを離散化
generateBIEMatrix()
BIEの離散化に要した時間: 18.4762
LU分解を使った解法に要した時間: 14.0018
update p->phiphin and p->phinOnFace for Dirichlet boundary
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {32.6106,114.006} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.021751,114.028} s
name = tank
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {net->velocityTranslational() = {0.0047955,0,0}
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
Elapsed time: {5.36118,119.389} s
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0047955,0,0}
name = tank
net->velocityTranslational() = {0,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
use tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.505627,119.895} s
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
Elapsed time: {6.4e-05,119.895} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 31040
原点を節点にとり，方程式を作成．並列化
線形要素の面の数：52312 persecent: 100
擬似二次要素の面の数：0 persecent: 0
線形要素を使ってBIEを離散化
generateBIEMatrix()
BIEの離散化に要した時間: 19.1142
LU分解を使った解法に要した時間: 14.2229
update p->phiphin and p->phinOnFace for Dirichlet boundary
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {33.471,153.366} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.021965,153.388} s
ALEのU_update_BEMを計算
Elapsed time: 1.86267 s
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
Elapsed time: {5.39613,160.647} s
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {name = 0.00958981,tank
0net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.504064,161.151} s
Total elapsed time: 161.151 s
============================================================================
ElapsedTimeBIEDiscretization={19.5933,18.6116,18.4762,19.1142}
ElapsedTimeSolve={14.9461,13.9537,14.0018,14.2229}
ElapsedTimeALE=1.86267
ElapsedTimeTotal=161.151
============================================================================
simulation_timeを取得
/Users/tomoaki/Ruehl20160d075/output/water_0.vtu  VV_points.size() : 52312 |||||||||
/Users/tomoaki/Ruehl20160d075/output/water.pvd
Creating /Users/tomoaki/Ruehl20160d075/output/water.pvd ...
Done.
/Users/tomoaki/Ruehl20160d075/output/tank_0.vtu  VV_points.size() : 52 |||||||||
/Users/tomoaki/Ruehl20160d075/output/tank.pvd
Creating /Users/tomoaki/Ruehl20160d075/output/tank.pvd ...
Done.
/Users/tomoaki/Ruehl20160d075/output/wavemaker_0.vtu  VV_points.size() : 496 |||||||||
/Users/tomoaki/Ruehl20160d075/output/wavemaker.pvd
Creating /Users/tomoaki/Ruehl20160d075/output/wavemaker.pvd ...
Done.
flipIf
flipIf
flipIf
net.getPoints() = 26158
Total : 26158
Total variables: 29014
Total case double-node : 29014
node reduction : 0.0984352
CORNER : 561
Total CORNER faces : 3417
Neumann : 18554
Dirichlet : 7043
===========================================================================
       dt :0.01
time_step :1
real time :0.01
---------------------------------------------------------------------------
RK_step = 1/4, RK_time = 0.01, simulation_time = 0.01
water makeBucketPoints(tank makeBucketPoints(2.01457)
2.99541)
wavemaker makeBucketPoints(0.616441)
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
      bounds : {{-4.99994,25},{-1.5,1.5},{-0.340041,1.7002}}
          dL : 2.01457
bucket sizes : {15, 2, 2}
water makeBucketPoints done
water makeBucketFaces(2.01457)
