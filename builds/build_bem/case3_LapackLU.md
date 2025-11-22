ケース：case3
14840節点
方法:Lapack LUで解く

$ ./main ~/Ruehl20160d1
input_directory : "/Users/tomoaki/Ruehl20160d1"
ALE: {pseudo_quad}
ALEPERIOD: {1}
GRAVITY: {9.81}
WATER_DENSITY: {1000.0}
element: {linear}
end_time: {10.0}
end_time_step: {100000}
input_files: {tank.json,wavemaker.json,water.json,wg3.json,wg6.json}
max_dt: {0.1}
output_directory: {/Users/tomoaki/Ruehl20160d1/output}
LINEAR_ELEMENT
ALE_ON_PSEUDO_QUADRATIC_ELEMENT

WATER_DENSITY: 1000
GRAVITY: 9.81
"/Users/tomoaki/Ruehl20160d1/tank.json"
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
        this->getName(): tank, 0x120008000
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
/Users/tomoaki/Ruehl20160d1/output/tank_init.vtu  VV_points.size() : 52 |||||||||
"/Users/tomoaki/Ruehl20160d1/wavemaker.json"
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
        this->getName(): wavemaker, 0x148008000
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
/Users/tomoaki/Ruehl20160d1/output/wavemaker_init.vtu  VV_points.size() : 496 |||||||||
"/Users/tomoaki/Ruehl20160d1/water.json"
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
        this->getName(): water, 0x126738000
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
/Users/tomoaki/Ruehl20160d1/output/water_init.vtu  VV_points.size() : 29676 |||||||||
"/Users/tomoaki/Ruehl20160d1/wg3.json"
name: {wg3}
position: {9.48,0.0,1.86,9.48,0.0,0.8600000000000001}
type: {wavegauge}
type = wavegauge
skipped
"/Users/tomoaki/Ruehl20160d1/wg6.json"
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
tankwavemaker makeBucketPoints( makeBucketPoints(2.99541)0.616441
)
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
Elapsed time: {0.944915,0.944915} s
waterの境界条件を決定 setBoundaryTypes
water setContactFaces()
step2 面の境界条件を判定
step3 線の境界条件を決定
step4 点の境界条件を決定
setBoundaryTypes終了
setBoundaryTypes
Elapsed time: {0.411861,1.35678} s
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
Elapsed time: {9e-05,1.35687} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18369
原点を節点にとり，方程式を作成．並列化
線形要素の面の数：29676 persecent: 100
擬似二次要素の面の数：0 persecent: 0
線形要素を使ってBIEを離散化
generateBIEMatrix()
BIEの離散化に要した時間: 5.90346
LU分解を使った解法に要した時間: 3.19394
update p->phiphin and p->phinOnFace for Dirichlet boundary
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {9.16507,10.5219} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.013724,10.5357} s
name = tank
net->velocityTranslational() = {updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0,0,0}
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
Elapsed time: {2.02873,12.5644} s
name = tankupdating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0,0,0}

net->velocityTranslational() = {use wavemaker's (RigidBodyObject) predetermiend velocity
0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.274928,12.8393} s
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
Elapsed time: {5.7e-05,12.8394} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18369
原点を節点にとり，方程式を作成．並列化
線形要素の面の数：29676 persecent: 100
擬似二次要素の面の数：0 persecent: 0
線形要素を使ってBIEを離散化
generateBIEMatrix()
BIEの離散化に要した時間: 4.93688
LU分解を使った解法に要した時間: 2.93553
update p->phiphin and p->phinOnFace for Dirichlet boundary
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {7.94593,20.7853} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.012482,20.7978} s
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
Elapsed time: {2.1188,22.9166} s
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {name = tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
0.0047955,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.276693,23.1933} s
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
Elapsed time: {7e-05,23.1934} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18369
原点を節点にとり，方程式を作成．並列化
線形要素の面の数：29676 persecent: 100
擬似二次要素の面の数：0 persecent: 0
線形要素を使ってBIEを離散化
generateBIEMatrix()
BIEの離散化に要した時間: 5.25348
LU分解を使った解法に要した時間: 2.95744
update p->phiphin and p->phinOnFace for Dirichlet boundary
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {8.28475,31.4781} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.012356,31.4905} s
name = updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0047955,0,tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
0}
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
Elapsed time: {2.13005,33.6205} s
name = tank
net->velocityTranslational() = {updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0047955,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.301537,33.922} s
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
Elapsed time: {6.9e-05,33.9221} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18369
原点を節点にとり，方程式を作成．並列化
線形要素の面の数：29676 persecent: 100
擬似二次要素の面の数：0 persecent: 0
線形要素を使ってBIEを離散化
generateBIEMatrix()
BIEの離散化に要した時間: 5.26831
LU分解を使った解法に要した時間: 2.95037
update p->phiphin and p->phinOnFace for Dirichlet boundary
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {8.2914,42.2135} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.011786,42.2253} s
ALEのU_update_BEMを計算
Elapsed time: 1.0765 s
name = tankupdating 
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
wavemaker's (RigidBodyObject) velocity
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
Elapsed time: {2.08969,45.3915} s
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {name = tank
net->velocityTranslational() = {0,0,0}
0.00958981,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
use tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.274113,45.6656} s
Total elapsed time: 45.6656 s
============================================================================
ElapsedTimeBIEDiscretization={5.90346,4.93688,5.25348,5.26831}
ElapsedTimeSolve={3.19394,2.93553,2.95744,2.95037}
ElapsedTimeALE=1.0765
ElapsedTimeTotal=45.6656
============================================================================
simulation_timeを取得
/Users/tomoaki/Ruehl20160d1/output/water_0.vtu  VV_points.size() : 29676 |||||||||
/Users/tomoaki/Ruehl20160d1/output/water.pvd
Creating /Users/tomoaki/Ruehl20160d1/output/water.pvd ...
Done.
/Users/tomoaki/Ruehl20160d1/output/tank_0.vtu  VV_points.size() : 52 |||||||||
/Users/tomoaki/Ruehl20160d1/output/tank.pvd
Creating /Users/tomoaki/Ruehl20160d1/output/tank.pvd ...
Done.
/Users/tomoaki/Ruehl20160d1/output/wavemaker_0.vtu  VV_points.size() : 496 |||||||||
/Users/tomoaki/Ruehl20160d1/output/wavemaker.pvd
Creating /Users/tomoaki/Ruehl20160d1/output/wavemaker.pvd ...
Done.
flipIf
flipIf
flipIf
net.getPoints() = 14840
Total : 14840
Total variables: 16912
Total case double-node : 16912
node reduction : 0.122517
CORNER : 402
Total CORNER faces : 2474
Neumann : 10391
Dirichlet : 4047
===========================================================================
       dt :0.01
time_step :1
real time :0.01
---------------------------------------------------------------------------
RK_step = 1/4, RK_time = 0.01, simulation_time = 0.01
waterwavemaker makeBucketPoints( makeBucketPoints(tank makeBucketPoints(2.014572.99541)
)
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
Elapsed time: {0.883221,0.883221} s
waterの境界条件を決定 setBoundaryTypes
water setContactFaces()
step2 面の境界条件を判定
step3 線の境界条件を決定
step4 点の境界条件を決定
setBoundaryTypes終了
setBoundaryTypes
Elapsed time: {0.41764,1.30086} s
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
Elapsed time: {7.1e-05,1.30093} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18369
原点を節点にとり，方程式を作成．並列化
線形要素の面の数：29676 persecent: 100
擬似二次要素の面の数：0 persecent: 0
線形要素を使ってBIEを離散化
generateBIEMatrix()
BIEの離散化に要した時間: 6.56581
LU分解を使った解法に要した時間: 3.12933
update p->phiphin and p->phinOnFace for Dirichlet boundary
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {9.77057,11.0715} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.012568,11.0841} s
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
Elapsed time: {2.09905,13.1831} s
name = tank
net->velocityTranslational() = {0,0,0}
updating wavemaker's (RigidBodyObject) velocityuse tank
name = wavemaker
net->velocityTranslational() = {0.00958981,0,0}
's (RigidBodyObject) predetermiend velocityuse wavemaker's (RigidBodyObject) predetermiend velocity

name:water: setBounds
Elapsed time: {0.273446,13.4566} s
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
Elapsed time: {5.7e-05,13.4566} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18369
原点を節点にとり，方程式を作成．並列化
線形要素の面の数：29676 persecent: 100
擬似二次要素の面の数：0 persecent: 0
線形要素を使ってBIEを離散化
generateBIEMatrix()
BIEの離散化に要した時間: 5.57095
LU分解を使った解法に要した時間: 2.96667
update p->phiphin and p->phinOnFace for Dirichlet boundary
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {8.60725,22.0639} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.010931,22.0748} s
name = tank
net->velocityTranslational() = {0,0,0}
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0143818,0,0}
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
Elapsed time: {2.19288,24.2677} s
name = tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0143818,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.277522,24.5452} s
RK_step = 3/4, RK_time = 0.015, simulation_time = 0.01
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
Elapsed time: {5.7e-05,24.5453} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18369
原点を節点にとり，方程式を作成．並列化
線形要素の面の数：29676 persecent: 100
擬似二次要素の面の数：0 persecent: 0
線形要素を使ってBIEを離散化
generateBIEMatrix()
BIEの離散化に要した時間: 5.25054
LU分解を使った解法に要した時間: 2.89297
update p->phiphin and p->phinOnFace for Dirichlet boundary
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {8.21501,32.7603} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.011792,32.7721} s
name = tankupdating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0143818,0,0}

net->velocityTranslational() = {use wavemaker's (RigidBodyObject) predetermiend velocity
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
Elapsed time: {2.09534,34.8674} s
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0143818,0,0}
name = tank
net->velocityTranslational() = {use wavemaker's (RigidBodyObject) predetermiend velocity
0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.274531,35.1419} s
RK_step = 4/4, RK_time = 0.02, simulation_time = 0.01
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
net->velocity = {0.0191702,0,0,0,0,0}
setting acceleration
指定がないので加速度はdefault_acceleration
setNeumannVelocity
Elapsed time: {5.6e-05,35.142} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18369
原点を節点にとり，方程式を作成．並列化
線形要素の面の数：29676 persecent: 100
擬似二次要素の面の数：0 persecent: 0
線形要素を使ってBIEを離散化
generateBIEMatrix()
BIEの離散化に要した時間: 5.20587
LU分解を使った解法に要した時間: 2.89254
update p->phiphin and p->phinOnFace for Dirichlet boundary
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {8.16884,43.3108} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.012351,43.3232} s
ALEのU_update_BEMを計算
Elapsed time: 1.11367 s
name = tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0191702,0,0}
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
Elapsed time: {2.17384,46.6107} s
name = tank
net->velocityTranslational() = updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {{0.0191702,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.294752,46.9054} s
Total elapsed time: 46.9055 s
============================================================================
ElapsedTimeBIEDiscretization={6.56581,5.57095,5.25054,5.20587}
ElapsedTimeSolve={3.12933,2.96667,2.89297,2.89254}
ElapsedTimeALE=1.11367
ElapsedTimeTotal=46.9055
============================================================================
simulation_timeを取得
/Users/tomoaki/Ruehl20160d1/output/water_1.vtu  VV_points.size() : 29676 |||||||||
/Users/tomoaki/Ruehl20160d1/output/water.pvd
Creating /Users/tomoaki/Ruehl20160d1/output/water.pvd ...
Done.
/Users/tomoaki/Ruehl20160d1/output/tank_1.vtu  VV_points.size() : 52 |||||||||
/Users/tomoaki/Ruehl20160d1/output/tank.pvd
Creating /Users/tomoaki/Ruehl20160d1/output/tank.pvd ...
Done.
/Users/tomoaki/Ruehl20160d1/output/wavemaker_1.vtu  VV_points.size() : 496 |||||||||
/Users/tomoaki/Ruehl20160d1/output/wavemaker.pvd
Creating /Users/tomoaki/Ruehl20160d1/output/wavemaker.pvd ...
Done.
flipIf
flipIf
flipIf
net.getPoints() = 14840
Total : 14840
Total variables: 16912
Total case double-node : 16912
node reduction : 0.122517
CORNER : 402
Total CORNER faces : 2474
Neumann : 10391
Dirichlet : 4047
===========================================================================
       dt :0.01
time_step :2
real time :0.02
---------------------------------------------------------------------------
RK_step = 1/4, RK_time = 0.02, simulation_time = 0.02
waterwavemaker makeBucketPoints(0.616441)
tank makeBucketPoints(2.99541)
 makeBucketPoints(2.01456)
BucketPoints.data1D.size() = 28
      bounds : {{-14.25,29.25},{-4.5,4.5},{-2.125,4.625}}
          dL : 2.99541
bucket sizes : {15, 4, 3}
tank makeBucketPoints done
tank makeBucketFaces(2.99541)
BucketPoints.data1D.size() = 250
      bounds : {{-6.24981,1.25019},{-1.5,1.5},{-0.75,3.75}}
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
      bounds : {{-6.24981,1.25019},{-1.5,1.5},{-0.75,3.75}}
          dL : 0.616441
bucket sizes : {13, 5, 8}
      bounds : {{-6.24981,1.25019},{-1.5,1.5},{-0.75,3.75}}
          dL : 0.616441
bucket sizes : {13, 5, 8}
wavemaker makeBucketFaces done
wavemaker makeBucketTetras(0.616441)
      bounds : {{-6.24981,1.25019},{-1.5,1.5},{-0.75,3.75}}
          dL : 0.616441
bucket sizes : {13, 5, 8}
wavemaker makeBucketTetras done
BucketPoints.data1D.size() = 14840
      bounds : {{-4.99976,25},{-1.5,1.5},{-0.340157,1.70078}}
          dL : 2.01456
bucket sizes : {15, 2, 2}
water makeBucketPoints done
water makeBucketFaces(2.01456)
      bounds : {{-4.99976,25},{-1.5,1.5},{-0.340157,1.70078}}
          dL : 2.01456
bucket sizes : {15, 2, 2}
      bounds : {{-4.99976,25},{-1.5,1.5},{-0.340157,1.70078}}
          dL : 2.01456
bucket sizes : {15, 2, 2}
water makeBucketFaces done
water makeBucketTetras(2.01456)
      bounds : {{-4.99976,25},{-1.5,1.5},{-0.340157,1.70078}}
          dL : 2.01456
bucket sizes : {15, 2, 2}
water makeBucketTetras done
makeBuckets
Elapsed time: {0.949448,0.949448} s
waterの境界条件を決定 setBoundaryTypes
water setContactFaces()
step2 面の境界条件を判定
step3 線の境界条件を決定
step4 点の境界条件を決定
setBoundaryTypes終了
setBoundaryTypes
Elapsed time: {0.424955,1.3744} s
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
net->velocity = {0.0191702,0,0,0,0,0}
setting acceleration
指定がないので加速度はdefault_acceleration
setNeumannVelocity
Elapsed time: {7.5e-05,1.37448} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18369
原点を節点にとり，方程式を作成．並列化
線形要素の面の数：29676 persecent: 100
擬似二次要素の面の数：0 persecent: 0
線形要素を使ってBIEを離散化
generateBIEMatrix()
BIEの離散化に要した時間: 6.47082
LU分解を使った解法に要した時間: 3.21344
update p->phiphin and p->phinOnFace for Dirichlet boundary
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {9.75609,11.1306} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.012787,11.1434} s
name = updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0191702,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
tank
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
Elapsed time: {2.0901,13.2335} s
name = tank
net->velocityTranslational() = {updating wavemaker's (RigidBodyObject) velocity
0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
name = wavemaker
net->velocityTranslational() = {0.0191702,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.283372,13.5168} s
RK_step = 2/4, RK_time = 0.025, simulation_time = 0.02
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
net->velocity = {0.0239538,0,0,0,0,0}
setting acceleration
指定がないので加速度はdefault_acceleration
setNeumannVelocity
Elapsed time: {5.8e-05,13.5169} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18369
原点を節点にとり，方程式を作成．並列化
線形要素の面の数：29676 persecent: 100
擬似二次要素の面の数：0 persecent: 0
線形要素を使ってBIEを離散化
generateBIEMatrix()
BIEの離散化に要した時間: 5.45565
LU分解を使った解法に要した時間: 2.98131
update p->phiphin and p->phinOnFace for Dirichlet boundary
BVP.solve -> {Φ,Φn}が決まる
Elapsed time: {8.50438,22.0213} s
U_BEMとU_update_BEMを計算
Elapsed time: {0.012072,22.0333} s
name = tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0239538,0,0}
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
Elapsed time: {2.20318,24.2365} s
name = tank
net->velocityTranslational() = {0,0,0}
use tank's (RigidBodyObject) predetermiend velocity
updating wavemaker's (RigidBodyObject) velocity
name = wavemaker
net->velocityTranslational() = {0.0239538,0,0}
use wavemaker's (RigidBodyObject) predetermiend velocity
name:water: setBounds
Elapsed time: {0.301165,24.5377} s
RK_step = 3/4, RK_time = 0.025, simulation_time = 0.02
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
net->velocity = {0.0239538,0,0,0,0,0}
setting acceleration
指定がないので加速度はdefault_acceleration
setNeumannVelocity
Elapsed time: {7.4e-05,24.5378} s
RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える
RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．
   unknown size : 18369
原点を節点にとり，方程式を作成．並列化
線形要素の面の数：29676 persecent: 100
擬似二次要素の面の数：0 persecent: 0
線形要素を使ってBIEを離散化
generateBIEMatrix()
BIEの離散化に要した時間: 5.551