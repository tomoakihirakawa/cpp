#ifndef lib_multipole_expansion_H
#define lib_multipole_expansion_H

#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include "basic_arithmetic_array_operations.hpp"
#include "basic_linear_systems.hpp"
#include "lib_multipole_expansion_constants.hpp"

const std::array<std::array<double, 31>, 31> sqrt_nm_nm = {{{1., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.7071067811865475, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.4082482904638631, 0.20412414523193154, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.2886751345948129, 0.09128709291752768, 0.03726779962499649, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.22360679774997896, 0.05270462766947299, 0.014085904245475275, 0.004980119205559973, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.18257418583505536, 0.03450327796711771, 0.007042952122737638, 0.0016600397351866577, 0.00052495065695726, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.1543033499620919, 0.024397501823713332, 0.0040662503039522215, 0.0007423923386456234, 0.00015827857841616382, 0.00004569108992776174, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.1336306209562122, 0.018184824186332698, 0.0025717224993681985, 0.0003877017543326962, 0.00006461695905544937, 0.000012672428274337013, 3.3868489186454308e-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.1178511301977579, 0.014085904245475275, 0.0017338549553676645, 0.00022383971222927231, 0.00003104098307414404, 4.7897276744570195e-6, 8.744806305356235e-7, 2.1862015763390587e-7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.10540925533894598, 0.011236664374387369, 0.0012260205965343742, 0.00013881949648441295, 0.000016592103373156565, 2.1420313347595755e-6, 3.0917559193401366e-7, 5.3023176577281006e-8, 1.2497682572615703e-8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.09534625892455924, 0.009174698042719672, 0.0008996531605990218, 0.00009087869293935404, 9.579455348914039e-6, 1.0710156673797877e-6, 1.2987972715583107e-7, 1.7674392192427002e-8, 2.867165019077962e-9, 6.411175885367803e-10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.08703882797784893, 0.00763381020690574, 0.0006800738654737899, 0.00006208196614828809, 5.866194404653597e-6, 5.808397976391274e-7, 6.12258905403645e-8, 7.023091305098204e-9, 9.066771887846482e-10, 1.3990332756368328e-10, 2.982748965711089e-11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.08006407690254357, 0.0064517471761563645, 0.0005267829510341783, 0.000043898579252848194, 3.764272115814879e-6, 3.3534801352299794e-7, 3.14082191407746e-8, 3.14082191407746e-9, 3.4269176584847724e-10, 4.218244040462945e-11, 6.2194615286459366e-12, 1.2695422683377336e-12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.07412493166611012, 0.005524946201098299, 0.0004164584894532388, 0.000031940908071889604, 2.509514743876586e-6, 2.0354852405013068e-7, 1.720299011446777e-8, 1.532564167533272e-9, 1.4612425993612895e-10, 1.5234507220052806e-11, 1.7954038938891261e-12, 2.539084536675467e-13, 4.9795544612191934e-14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.06900655593423542, 0.0047847437646934545, 0.00033499907004387234, 0.0000238073472372723, 1.7271664998763255e-6, 1.2873539007279484e-7, 9.93215097345444e-9, 8.003555336576621e-10, 6.813078749549884e-11, 6.2194615286459366e-12, 6.219461528645936e-13, 7.042153453631633e-14, 9.583157028764344e-15, 1.8110464480707933e-15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.06454972243679027, 0.004184137043778615, 0.00027352559530478956, 0.000018114675073333608, 1.221291144300784e-6, 8.427709566117844e-8, 5.989312389111183e-9, 4.4153796733665665e-10, 3.406539374774942e-11, 2.7814277522994136e-12, 2.4394735152773264e-13, 2.3473844845438778e-14, 2.5612064489493374e-15, 3.363028826288546e-16, 6.140022498994532e-17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.06063390625908324, 0.0036900620230837303, 0.00022625221914921153, 0.00001403156697619594, 8.839056362752623e-7, 5.6819606303568e-8, 3.746573890260891e-9, 2.549220642992589e-10, 1.802571203400791e-11, 1.3361546727538149e-12, 1.049782255353683e-13, 8.872279396506633e-15, 8.237704614677945e-16, 8.683303091354011e-17, 1.1027805953831061e-17, 1.9494590928908317e-18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.057166195047502935, 0.0032787061473526402, 0.00018929618767677276, 0.00001103997785862415, 6.528075759249468e-7, 3.929439846428888e-8, 2.4184030470526603e-9, 1.5295323857955533e-10, 9.99886600236366e-12, 6.803366586737007e-13, 4.859547561955005e-14, 3.684013499396719e-15, 3.0079844263490085e-16, 2.7012497569312447e-17, 2.7569514884577653e-18, 3.39357269271957e-19, 5.819929153796601e-20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.05407380704358751, 0.0029325639294907786, 0.0001599844784189123, 8.806855686093134e-6, 4.907868795155064e-7, 2.7785335616744933e-8, 1.6041870997851784e-9, 9.485757175887514e-11, 5.772847978055655e-12, 3.6365524063136487e-13, 2.3875136219024277e-14, 1.647540922935589e-15, 1.2080356161406092e-16, 9.55036010402298e-18, 8.312521502205675e-19, 8.230622741349724e-20, 9.837475773099464e-21, 1.6395792955165772e-21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.051298917604257706, 0.002638531611922098, 0.00013643522161876318, 7.112177573145887e-6, 3.7484467091183e-7, 2.0036290454429938e-8, 1.0898305269346695e-9, 6.054614038525942e-11, 3.4499365303551642e-12, 2.0258724866669537e-13, 1.2329067328628607e-14, 7.828965582648662e-16, 5.230947661270771e-17, 3.717472628672112e-18, 2.8511713531899113e-19, 2.4096796001085158e-20, 2.3187152763056684e-21, 2.6954517578975646e-22, 4.37260014755555e-23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.048795003647426664, 0.0023866416420646975, 0.00011729703744891826, 5.807068671424468e-6, 2.903534335712234e-7, 1.4702617627274514e-8, 7.562210510333867e-10, 3.963675305041057e-11, 2.1247529892807687e-12, 1.1696380255210221e-13, 6.64309758380035e-15, 3.914482791324331e-16, 2.4091976159852608e-17, 1.561650865454643e-18, 1.07764147796743e-19, 8.032265333695052e-21, 6.602481433136971e-22, 6.183790431507657e-23, 7.001763889559219e-24, 1.1070760764863385e-24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.046524210519923545, 0.0021692025004706648, 0.00010158221421941786, 4.7886315014996605e-6, 2.277719728643592e-7, 1.095868415392416e-8, 5.347290332617261e-10, 2.6538156330298823e-11, 1.3438117822757118e-12, 6.967339525147562e-14, 3.713604444635547e-15, 2.0442719657779699e-16, 1.1686324992580568e-17, 6.98391498455592e-19, 4.39945291113154e-20, 2.9527194609348998e-21, 2.142127842145906e-22, 1.7150748828864954e-23, 1.5656420019457741e-24, 1.728962355617497e-25, 2.667846834301295e-26, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.044455422447438706, 0.001980201902173989, 0.00008855732124870858, 3.9843822502982624e-6, 1.8073524913764726e-7, 8.283986562225028e-9, 3.845744215316325e-10, 1.8129012089060738e-11, 8.702204277741714e-13, 4.2666066753341936e-14, 2.144050525774124e-15, 1.108662133605148e-16, 5.9260483737718005e-18, 3.2922490965398895e-19, 1.9135809181461728e-20, 1.1732917402031786e-21, 7.670048049273388e-23, 5.423542987667864e-24, 4.2350755557433e-25, 3.7729051754030186e-26, 4.068428280973127e-27, 6.133386406572036e-28, 0, 0, 0, 0, 0, 0, 0, 0}, {1., 0.04256282653793743, 0.0018148850216015693, 0.00007766990876952864, 3.3423807018425985e-6, 1.449105814256371e-7, 6.342563843542926e-9, 2.808534476165007e-10, 1.2610692071105653e-11, 5.755967094246757e-13, 2.677918248384888e-14, 1.2737558160678573e-15, 6.215291969096202e-17, 3.123301730909286e-18, 1.623727252284406e-19, 8.78011141314415e-21, 4.970759255166987e-22, 2.9705968359446517e-23, 1.893983366497944e-24, 1.3069726911875188e-25, 9.965573343492674e-27, 8.673918239448955e-28, 9.143112624912157e-29, 1.3480784107615843e-29, 0, 0, 0, 0, 0, 0, 0}, {1., 0.04082482904638631, 0.0016694514082354448, 0.0000684984209857492, 2.8248272710847575e-6, 1.1729456811558224e-7, 4.912928827659812e-9, 2.0798086168523098e-10, 8.917105878934236e-12, 3.8806722232841114e-13, 1.7183905434211925e-14, 7.762897181293485e-16, 3.588400491449811e-17, 1.7029795092108926e-18, 8.32954505267255e-20, 4.217829092369243e-21, 2.222991118867921e-22, 1.2274415085370757e-23, 7.158584250066375e-25, 4.4567398861618955e-26, 3.0047334182847955e-27, 2.2395960591833703e-28, 1.906470771641401e-29, 1.966374459242337e-30, 2.83821705842792e-31, 0, 0, 0, 0, 0, 0}, {1., 0.03922322702763681, 0.0015408338784034144, 0.0000607173606927742, 2.403823329189731e-6, 9.577061382776717e-8, 3.8462404083659755e-9, 1.5598564626392325e-10, 6.400168590000365e-12, 2.6621191827802367e-13, 1.1249506769605595e-14, 4.841016930192229e-16, 2.127020725267102e-17, 9.569918674447674e-19, 4.423696472146341e-20, 2.1089145461846216e-21, 1.0415186570356942e-22, 5.3569939276208894e-24, 2.888297556246701e-25, 1.645760272484188e-26, 1.0015778060949317e-27, 6.604208479255898e-29, 4.816614068384924e-30, 4.013845056987437e-31, 4.054595797754172e-32, 5.734064367124909e-33, 0, 0, 0, 0, 0}, {1., 0.037742567804819854, 0.0014265349750363762, 0.00005407266868359784, 2.058511316097547e-6, 7.8824472197141e-8, 3.0407200297531575e-9, 1.183598666326814e-10, 4.656806667940934e-12, 1.8553161858873915e-13, 7.499671179737064e-15, 3.082345867679868e-16, 1.2910526133844702e-17, 5.525195122815191e-19, 2.4229583453505744e-20, 1.0923539830333294e-21, 5.082090666891961e-23, 2.4508017457278914e-24, 1.2315742160799112e-25, 6.49096605058148e-27, 3.6172739585553436e-28, 2.1540552956581018e-29, 1.3904367144822869e-30, 9.93169081773062e-32, 8.109191595508343e-33, 8.029296133819685e-34, 1.113463035245063e-34, 0, 0, 0, 0}, {1., 0.036369648372665396, 0.0013245044730538606, 0.000048364065160539534, 1.7731126566244665e-6, 6.535779964392945e-8, 2.4256554521012607e-9, 9.077787379191073e-11, 3.431081122865768e-12, 1.3119066562861356e-13, 5.083535034852086e-15, 2.0000901082862465e-16, 8.006768412092335e-18, 3.268749516376755e-19, 1.3643496875171323e-20, 5.838877642269419e-22, 2.5704225560302658e-23, 1.1683738891046663e-24, 5.507767332981392e-26, 2.706919965402724e-27, 1.3959873821192436e-28, 7.615735533053188e-30, 4.441587159991229e-31, 2.809106370346196e-32, 1.9667678521559944e-33, 1.5746745256446799e-34, 1.5294591045728086e-35, 2.0813302159207996e-36, 0, 0, 0}, {1., 0.03509312031717982, 0.0012330465600513202, 0.00004343221265200662, 1.5355606044085026e-6, 5.456377014455781e-8, 1.9511964368108084e-9, 7.03162386800378e-11, 2.5573768756813884e-12, 9.401104366027673e-14, 3.498728762031909e-15, 1.3205100753166277e-16, 5.063924975940288e-18, 1.9771305335356593e-19, 7.87707659356798e-21, 3.210456743489421e-22, 1.3423594001082063e-23, 5.776595112493402e-25, 2.5680097603370388e-26, 1.184534573216848e-27, 5.699094622593048e-29, 2.8784774673280942e-30, 1.5386109254700724e-31, 8.795653226764121e-33, 5.454832567602014e-34, 3.7463943886574394e-35, 2.943445419132117e-36, 2.8064650906287853e-37, 3.75029672752657e-38, 0, 0}, {1., 0.033903175181040524, 0.0011507487481140638, 0.00003914926743091642, 1.3365339529667335e-6, 4.584273638574386e-8, 1.5817246762906214e-9, 5.4968732346059175e-11, 1.926653071686668e-12, 6.820277942358616e-14, 2.442051239803031e-15, 8.858250880144039e-17, 3.2607672101605292e-18, 1.2203114585475613e-19, 4.652393819724967e-21, 1.8109418382397407e-22, 7.214964983015853e-24, 2.950418589035446e-25, 1.2423503413532843e-26, 5.406635882463521e-28, 2.4424691239684485e-29, 1.1513909869312377e-30, 5.700234783577265e-32, 2.9877348446894465e-33, 1.67543850471236e-34, 1.0196394030485391e-35, 6.874407452974273e-37, 5.303720494991511e-38, 4.9673893642040876e-39, 6.52250178038964e-40, 0}, {1., 0.032791291789197645, 0.001076426888425463, 0.000035411844786990176, 1.1687650369369893e-6, 3.874418370480919e-8, 1.2914727901603063e-9, 4.333898916301203e-11, 1.4659631761914363e-12, 5.0047157644424506e-14, 1.7267909916697387e-15, 6.030215715220695e-17, 2.134673223455801e-18, 7.672924452278737e-20, 2.8054990321811047e-21, 1.0455477577944701e-22, 3.9803322886358565e-24, 1.5516950790121656e-25, 6.2117517067664216e-27, 2.5616832286005083e-28, 1.0923053988275624e-29, 4.836809603787499e-31, 2.2358160302708322e-32, 1.0858082505169e-33, 5.584795015707868e-35, 3.074328473976332e-36, 1.8372626732939303e-37, 1.2167569296028167e-38, 9.224210478429689e-40, 8.49157403659412e-41, 1.0962574942272106e-41}}};

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT Multipole_Expansion

# 多重極展開

この実装は，\cite{Greengard1997a}に基づいている．

## Green関数の多重極展開

次のGreen関数を考える．

```math
G({\bf x},{\bf a}) = \frac{1}{\|{\bf x}-{\bf a}\|},
\quad \nabla G({\bf x},{\bf a}) = -\frac{{\bf x}-{\bf a}}{\|{\bf x}-{\bf a}\|^3}
```

グリーン関数は，球面調和関数を使って近似できる．
近似を$`G_{\rm apx}({\bf x},{\bf a},{\bf c})`$とする．

```math
G_{\rm apx}(n, {\bf x},{\bf a},{\bf c}) = \sum_{k=0}^n \sum_{m=-k}^k \left( \frac{r_{\rm near}}{r_{\rm far}} \right)^k \frac{1}{r_{\rm far}} Y(k, -m, a_{\rm near}, b_{\rm near}) Y(k, m, a_{\rm far}, b_{\rm far})=
{\bf Y}^\ast({\bf x},{\bf c})\cdot{\bf Y}({\bf a},{\bf c})
```

```math
{\bf Y}^\ast({\bf x},{\bf c}) = r_{\rm near}^k Y(k, -m, a_{\rm near},b_{\rm near}), \quad {\bf Y}({\bf a},{\bf c}) = r_{\rm far}^{-k-1} Y(k, m, a_{\rm far}, b_{\rm far})
```

ここで，$`(r_{\rm near},a_{\rm near},b_{\rm near})`$は，球面座標系に$`{\bf x}-{\bf c}`$を変換したものであり，
$`(r_{\rm far},a_{\rm far},b_{\rm far})`$は，球面座標系に$`{\bf a}-{\bf c}`$を変換したもの．$`Y(k, m, a, b)`$は球面調和関数：

```math
Y(k, m, a, b) = \sqrt{\frac{(k - |m|)!}{(k + |m|)!}} P_k^{|m|}(\cos(a)) e^{i mb}
```

$`P_k^m(x)`$はルジャンドル陪関数：

```math
P_k^m(x) = \frac{(-1)^m}{2^k k!} (1-x^2)^{m/2} \frac{d^{k+m}}{dx^{k+m}}(x^2-1)^k
```

### 球面座標系への変換

$`{\bf x}=(x,y,z)`$から球面座標$`(r,a,b)`$への変換は次のように行う．

```math
r = \|{\bf x}\|, \quad a = \arctan \frac{\sqrt{x^2 + y^2}}{z}, \quad b = \arctan \frac{y}{x}
```

$`r_\parallel=\sqrt{x^2+y^2}`$とする．$`\frac{\partial}{\partial t}(\arctan(f(t))) = \frac{f'(t)}{1 + f(t)^2}`$なので，
$`(r,a,b)`$の$`(x,y,z)`$に関する勾配は次のようになる．

```math
\nabla r = \frac{\bf x}{r},\quad
\nabla a = \frac{1}{r^2r_\parallel} \left(xz,yz,-r_\parallel^2\right),\quad
\nabla b = \frac{1}{r_\parallel^2} \left(-y,x,0\right)
```

*/

constexpr double
factorial(const int n) {
   if (n < 0)
      throw std::runtime_error("factorial of negative number");
   else if (n < factorial_list_as_double.size())
      return factorial_list_as_double[n];
   else
      return n * factorial(n - 1);
}

/* -------------------------------------------------------------------------- */

/*

\cite{Liu_2009}p.70には，体球調和関数を使った高速多重極展開のアルゴリズムが書かれている．
以下は，その定義通りの体球調和関数を計算する関数である．

*/

constexpr std::complex<double> complex_zero(0., 0.);
constexpr std::array<std::complex<double>, 3> complex_zero3D{complex_zero, complex_zero, complex_zero};

std::complex<double> Dot(const std::array<double, 3>& normal, const std::array<std::complex<double>, 3>& abc) {
   return std::get<0>(abc) * std::get<0>(normal) + std::get<1>(abc) * std::get<1>(normal) + std::get<2>(abc) * std::get<2>(normal);
};

std::complex<double> Dot(const std::array<std::complex<double>, 3>& abc, const std::array<double, 3>& normal) {
   return Dot(normal, abc);
};

std::array<std::complex<double>, 3> operator*(std::array<std::complex<double>, 3> a, const std::array<std::complex<double>, 3>& A) {
   std::get<0>(a) *= std::get<0>(A);
   std::get<1>(a) *= std::get<1>(A);
   std::get<2>(a) *= std::get<2>(A);
   return a;
};

std::array<std::complex<double>, 3> operator*(std::complex<double> a, std::array<std::complex<double>, 3> A) {
   std::get<0>(A) *= a;
   std::get<1>(A) *= a;
   std::get<2>(A) *= a;
   return A;
};

std::array<std::complex<double>, 3> operator*(std::complex<double> a, std::array<double, 3> A) {
   std::array<std::complex<double>, 3> B = {a, a, a};
   std::get<0>(B) *= std::get<0>(A);
   std::get<1>(B) *= std::get<1>(A);
   std::get<2>(B) *= std::get<2>(A);
   return B;
};

std::array<std::complex<double>, 3> operator+(std::array<std::complex<double>, 3> a, const std::array<double, 3>& A) {
   std::get<0>(a) += std::get<0>(A);
   std::get<1>(a) += std::get<1>(A);
   std::get<2>(a) += std::get<2>(A);
   return a;
};

std::array<std::complex<double>, 3> operator+(const std::array<double, 3>& A, const std::array<std::complex<double>, 3>& a) {
   return a + A;
};

struct SphericalCoordinates {
   std::array<double, 3> X;
   double r2D;
   double rho;
   double div_r2D;
   double theta, phi;
   double x_r, y_r, z_r;
   SphericalCoordinates(const std::array<double, 3>& X)
       : X(X),
         r2D(std::hypot(std::get<0>(X), std::get<1>(X))),
         rho(std::hypot(std::get<0>(X), std::get<1>(X), std::get<2>(X))),
         div_r2D(1. / r2D),
         theta(std::atan2(r2D, std::get<2>(X))),
         phi(std::atan2(std::get<1>(X), std::get<0>(X))),
         x_r(std::get<0>(X) / rho),
         y_r(std::get<1>(X) / rho),
         z_r(std::get<2>(X) / rho)
   // Jacobian({{/*grad_sph x*/ {x_r, y_r, z_r},
   //            /*grad_sph y*/ {x_r * z_r * div_r2D, y_r * z_r * div_r2D, -r2D * div_rho * div_rho},
   //            /*grad_sph zr*/ {-std::get<1>(X) * div_r2D * div_r2D, std::get<0>(X) * div_r2D * div_r2D, 0.}}}),
   {
      // compute_precomputed_values();
   }

   // int pre_max = -1;  // 0,1,2,...,pre_max
   int pre_max = -1;  // 0,1,2,...,pre_max
   // std::array<std::array<std::complex<double>, 2 * pre_max + 1>, pre_max + 1> pre_compute_sph_harmonics, pre_compute_sph_harmonics_rho_n;
   std::vector<std::vector<std::complex<double>>> pre_compute_sph_harmonics, pre_compute_sph_harmonics_rho_n;
   void precompute_sph(const int pre_max_IN) {
      this->pre_max = pre_max_IN;
      pre_compute_sph_harmonics.resize(pre_max + 1, std::vector<std::complex<double>>(2 * pre_max + 1, 0.));
      pre_compute_sph_harmonics_rho_n.resize(pre_max + 1, std::vector<std::complex<double>>(2 * pre_max + 1, 0.));
      double rho_n;
      for (int n = 0; n <= pre_max; n++) {
         rho_n = std::pow(rho, n);  // runtime value, cannot be constexpr
         for (int m = -pre_max; m <= pre_max; m++) {
            // sph = this->sph_harmonics_(n, m);  // runtime value, cannot be constexpr
            // this->pre_compute_sph_harmonics[n][m + pre_max] = sph;
            // this->pre_compute_sph_harmonics_rho_n[n][m + pre_max] = sph * rho_n;

            this->pre_compute_sph_harmonics_rho_n[n][m + pre_max] = (this->pre_compute_sph_harmonics[n][m + pre_max] = this->sph_harmonics_(n, m)) * rho_n;
         }
      }
   }

   std::array<std::array<double, 3>, 3> Jacobian_inv() const {
      return {{/*grad r*/ {x_r, y_r, z_r},
               /*grad theta*/ {x_r * z_r * div_r2D, y_r * z_r * div_r2D, -r2D / rho / rho},
               /*grad phi*/ {-std::get<1>(X) * div_r2D * div_r2D, std::get<0>(X) * div_r2D * div_r2D, 0.}}};
   }

   // void initialize(const std::array<double, 3>& X) {
   //    this->X = X;
   //    div_r2D = 1. / (r2D = std::hypot(std::get<0>(X), std::get<1>(X)));
   //    div_rho = 1. / (rho = std::hypot(std::get<0>(X), std::get<1>(X), std::get<2>(X)));
   //    theta = std::atan2(r2D, std::get<2>(X));
   //    phi = std::atan2(std::get<1>(X), std::get<0>(X));
   //    x_r = std::get<0>(X) * div_rho;
   //    y_r = std::get<1>(X) * div_rho;
   //    z_r = std::get<2>(X) * div_rho;
   //    Jacobian_inv = {{{x_r, y_r, z_r},
   //                     {x_r * z_r * div_r2D, y_r * z_r * div_r2D, -r2D * div_rho * div_rho},
   //                     {-std::get<1>(X) * div_r2D * div_r2D, std::get<0>(X) * div_r2D * div_r2D, 0.}}};
   //    compute_precomputed_values();
   // }

   inline int negOnePower(const int p) const {
      return 1 - ((p & 1) << 1);
   }

   inline std::complex<double> sph_harmonics_(const int l, const int m) const {
      constexpr double sqrt_MPI2 = std::sqrt(2.0 * M_PI);
      if (std::abs(m) <= l) {
         double s = sqrt_MPI2 / std::sqrt(l + 0.5) * negOnePower(m);
         if (m < 0)
            return std::polar(std::sph_legendre(l, -m, theta) * s, m * phi);
         else
            return std::polar(std::sph_legendre(l, m, theta) * s, m * phi);
      } else
         return 0.0;
   }

   inline std::complex<double> sph_harmonics(const int l, const int m) const {
      if (l <= pre_max && std::abs(m) <= pre_max)
         return this->pre_compute_sph_harmonics[l][m + pre_max];
      else
         return this->sph_harmonics_(l, m);
   }

   inline std::complex<double> sph_harmonics_rho(const int l, const int m) const {

      if (l > pre_max || m > pre_max || m < -pre_max)
         return this->sph_harmonics_(l, m) * std::pow(rho, l);
      else
         return this->pre_compute_sph_harmonics_rho_n[l][m + pre_max];

      // if (l <= pre_max && std::abs(m) <= pre_max)
      //    return this->pre_compute_sph_harmonics_rho_n[l][m + pre_max];
      // else
      //    return this->sph_harmonics_(l, m) * std::pow(rho, l);
   }

   inline std::complex<double> SolidHarmonicR(const int n, const int m) const { return std::pow(rho, n) * sph_harmonics(n, -m); }

   inline std::complex<double> SolidHarmonicS(const int n, const int m) const { return std::pow(rho, -n - 1) * sph_harmonics(n, m); }

   inline std::complex<double> grad_R_theta(const int n, const int m) const {
      double abs_m = std::abs(m);
      double Pnm = std::pow(-1, abs_m) * std::assoc_legendre(n, abs_m, std::cos(theta));
      double Pn1m = 0;
      if (n >= abs_m + 1)
         Pn1m = std::pow(-1, std::abs(m + 1)) * std::assoc_legendre(n, abs_m + 1, std::cos(theta));
      return std::polar(std::pow(rho, n) * sqrt_nm_nm[n][abs_m] * (-abs_m / 2. * (Pnm / std::sin(theta)) + Pn1m), m * phi);
   };

   inline std::array<std::complex<double>, 2> SolidHarmonicR_Grad_SolidHarmonicR_normal(const int n, const int m, const std::array<double, 3>& normal) const {

      // std::complex<double> Rnm = std::pow(rho, n) * sph_harmonics(n, -m);
      std::complex<double> Rnm = sph_harmonics_rho(n, -m);

      auto J = Jacobian_inv();

      // double c = std::sqrt((n - std::abs(m)) * (n + std::abs(m) + 1));
      std::complex<double> grad_Rnm1_dot_normal = (/*rについての微分*/ (double)n * (Rnm / rho)) * Dot(normal, std::get<0>(J));
      // grad_Rnm1_dot_normal += (/*thetaについての微分*/ (double)m / std::tan(theta) * Rnm + std::polar(c, phi) * Rnm1 / (-std::sin(theta))) * Dot(normal, J[1]);
      // grad_Rnm1_dot_normal += (/*thetaについての微分*/ (double)m * (std::get<2>(X) / r2D) * Rnm + std::polar(c, phi) * Rnm1) * Dot(normal, J[1]);
      grad_Rnm1_dot_normal += grad_R_theta(n, -m) * Dot(normal, std::get<1>(J));
      grad_Rnm1_dot_normal += (/*phiについての微分*/ -std::complex<double>(0., m) * Rnm) * Dot(normal, std::get<2>(J));
      return {Rnm, grad_Rnm1_dot_normal};
   }

   inline std::array<std::complex<double>, 3> Grad_SolidHarmonicR(const int n, const int m) const {

      // std::complex<double> Rnm = std::pow(rho, n) * sph_harmonics(n, -m);
      std::complex<double> Rnm = sph_harmonics_rho(n, -m);
      auto J = Jacobian_inv();
      // std::complex<double> Rnm1 = std::pow(rho, n) * sph_harmonics(n, -m - 1);
      //
      // double c = std::sqrt((n - std::abs(m)) * (n + std::abs(m) + 1));
      auto grad_Rnm1_dot_normal = (/*rについての微分*/ (double)n / rho * Rnm) * J[0];
      // grad_Rnm1_dot_normal += (/*thetaについての微分*/ (double)m / std::tan(theta) * Rnm + std::polar(c, phi) * Rnm1) * J[1];
      // grad_Rnm1_dot_normal += (/*thetaについての微分*/ (double)m * (std::get<2>(X) / r2D) * Rnm + std::polar(c, phi) * Rnm1) * J[1];
      grad_Rnm1_dot_normal += grad_R_theta(n, -m) * J[1];
      grad_Rnm1_dot_normal += (/*phiについての微分*/ -std::complex<double>(0., m) * Rnm) * J[2];
      //
      return grad_Rnm1_dot_normal;
   }
};

template <typename T>
void complex_fused_multiply_increment(const std::complex<T>& a, const std::complex<T>& b, std::complex<T>& c) {
   c.real(c.real() + a.real() * b.real() - a.imag() * b.imag());
   c.imag(c.imag() + a.real() * b.imag() + a.imag() * b.real());
}

/* -------------------------------------------------------------------------- */

#include <functional>
#include "basic_exception.hpp"

struct pole4FMM {
   /*DOC_EXTRACT
   境界要素法において，pole4FMMは，ある節点における変数値を表すわけではなく，
   数値的な面積分における和のある一項に相当する．
   この一項は，複数の節点における変数値の積に重みをかけたものである．
   */
   Tddd X;
   Tdd weights;
   Tddd normal;
   std::function<Tdd()> get_values;
   pole4FMM(const Tddd& X, const Tdd& weights, const Tddd& normal, std::function<Tdd()> get_values)
       : X(X), weights(weights), normal(normal), get_values(get_values) {}
};

template <int N>
struct ExpCoeffs {

   /*

   N = 0 [ 0,  ,  ,  ,  ] n=0   total 1 = (N+1)^2             center : 0, in 1D 0
   N = 1 [-1, 0, 1,  ,  ] n=1   total 1 + 3 = 4 = (N+1)^2     center : 1, in 1D 2
   N = 2 [-2,-1, 0, 1, 2] n=2=N total 4 + 5 = 9 = (N+1)^2     center : 2, in 1D 6
         <-- 2*N+1=5 -->
   */

   /* -------------------------------------------------------------------------- */

   std::array<double, 3> X;
   std::array<std::complex<double>, (N + 1) * (N + 1)> coeffs;
   std::array<std::complex<double>, (N + 1) * (N + 1)> coeffs_;
   // std::array<std::array<int, 2>, (N + 1) * (N + 1)> index_map;
   std::array<std::array<int, 2>, (N + 1) * (N + 1)> nm_set;
   const std::complex<double> zero = {0., 0.};
   int index(int n, int m) const { return n * (n + 1) + m; }
   const std::complex<double>& get_coeffs(const int n, const int m) const { return coeffs[index(n, m)]; }
   const std::complex<double>& get_coeffs_(const int n, const int m) const { return coeffs_[index(n, m)]; }

   // nのサイズはN+1
   // mのサイズは2N+1
   // 正しく離散畳み込みを行うには，二つの配列のサイズの合計プラス1が必要なのd，2*N+1+2*N+1+1となる
   // それを，各nに対して，n+1個の配列を用意する
   std::array<std::array<std::complex<double>, (2 * N + 1) + (2 * N + 1) + 1>, N + 1> DFT_Cn;

   void set_nm_set() {
      int ind = 0;
      for (int n = 0; n <= N; ++n)
         for (int m = -n; m <= n; ++m)
            this->nm_set[ind++] = {n, m};
   }

   /* -------------------------------------------------------------------------- */

   void initialize(const std::array<double, 3>& XIN) {
      this->X = XIN;
      std::fill(coeffs.begin(), coeffs.end(), std::complex<double>(0.0, 0.0));
      std::fill(coeffs_.begin(), coeffs_.end(), std::complex<double>(0.0, 0.0));
      set_nm_set();
   }

   // Default constructor
   ExpCoeffs() : X{0.0, 0.0, 0.0} {
      std::fill(coeffs.begin(), coeffs.end(), std::complex<double>(0.0, 0.0));
      std::fill(coeffs_.begin(), coeffs_.end(), std::complex<double>(0.0, 0.0));
      set_nm_set();
   }

   ExpCoeffs(const std::array<double, 3>& XIN) : X(XIN) {
      std::fill(coeffs.begin(), coeffs.end(), std::complex<double>(0.0, 0.0));
      std::fill(coeffs_.begin(), coeffs_.end(), std::complex<double>(0.0, 0.0));
      set_nm_set();
   }

   /* -------------------------------------------------------------------------- */

   // void increment(const Tddd& XIN, const std::array<double, 2> weights, const Tddd& normal) {
   //    SphericalCoordinates R(XIN - this->X);
   //    auto set_coeffs = [&](int n, int m) -> std::array<std::complex<double>, 2> {
   //       return {R.SolidHarmonicR(n, m) * std::get<0>(weights), Dot(normal, R.Grad_SolidHarmonicR(n, m)) * std::get<1>(weights)};
   //    };
   //    this->increment_coeffs(set_coeffs);
   // }

   // template <typename T>
   // void increment_moments(const T& poles) {
   //    for (const auto& pole : poles) {
   //       auto weights = pole->weights;
   //       // SphericalCoordinates R(this->X - pole->X);
   //       // SphericalCoordinates R(pole->X - this->X);
   //       SphericalCoordinates R(pole->X - this->X);
   //       if constexpr (std::is_pointer_v<typename T::value_type> || std::is_same_v<typename T::value_type, std::shared_ptr<pole4FMM>>) {
   //          this->increment_coeffs([&](int n, int m) -> std::array<std::complex<double>, 2> {
   //             // return {R.SolidHarmonicR(n, m) * std::get<0>(pole->weights),
   //             //         R.Grad_SolidHarmonicR_normal(n, m, pole->normal) * std::get<1>(pole->weights)};
   //             return R.SolidHarmonicR_Grad_SolidHarmonicR_normal(n, m, pole->normal) * weights;
   //          });
   //       } else {
   //          this->increment_coeffs([&](int n, int m) -> std::array<std::complex<double>, 2> {
   //             // return {R.SolidHarmonicR(n, m) * std::get<0>(pole.weights),
   //             //         R.Grad_SolidHarmonicR_normal(n, m, pole.normal) * std::get<1>(pole.weights)};
   //             return R.SolidHarmonicR_Grad_SolidHarmonicR_normal(n, m, pole->normal) * weights;
   //          });
   //       }
   //    }
   // }
   std::vector<std::tuple<std::function<Tdd()>, std::vector<std::tuple<int, int, int, std::array<std::complex<double>, 2>>>>> saved_getValue_SolidHarmonicR_Grad_SolidHarmonicR_normal;

   void increment_moments_reuse() {
      coeffs.fill(0.0);
      coeffs_.fill(0.0);

      for (const auto& [get_values_func, tuple_vector] : saved_getValue_SolidHarmonicR_Grad_SolidHarmonicR_normal) {
         auto coefficients = get_values_func();

         for (const auto& [ind, n, m, values] : tuple_vector) {
            coeffs[ind] += std::get<0>(coefficients) * std::get<0>(values);
            coeffs_[ind] += std::get<1>(coefficients) * std::get<1>(values);
         }
      }
   }

   template <typename T>
   void increment_moments(const T& poles) {
      for (const auto& pole : poles) {
         auto weights = pole->get_values() * pole->weights;  // Call the function to get the values
         SphericalCoordinates R(pole->X - this->X);

         this->increment_coeffs3([&](int ind, int n, int m) -> std::array<std::complex<double>, 2> {
            auto tmp = R.SolidHarmonicR_Grad_SolidHarmonicR_normal(n, m, pole->normal) * weights;

            saved_getValue_SolidHarmonicR_Grad_SolidHarmonicR_normal.emplace_back(
                pole->get_values,
                std::vector<std::tuple<int, int, int, std::array<std::complex<double>, 2>>>{std::make_tuple(ind, n, m, tmp)});

            return tmp;
         });
      }
   }

   void increment_coeffs(const auto& func) {
      for (size_t ind = 0; ind < nm_set.size(); ++ind) {
         auto [n, m] = nm_set[ind];
         auto coefficients = func(n, m);
         coeffs[ind] += std::get<0>(coefficients);
         coeffs_[ind] += std::get<1>(coefficients);
      }
   }

   void increment_coeffs3(const auto& func) {
      for (size_t ind = 0; ind < nm_set.size(); ++ind) {
         auto [n, m] = nm_set[ind];
         auto coefficients = func(ind, n, m);
         coeffs[ind] += std::get<0>(coefficients);
         coeffs_[ind] += std::get<1>(coefficients);
      }
   }

   std::array<double, 2> IGIGn_using_L(const std::array<double, 3>& a) const {
      // SphericalCoordinates Xc2O(this->X - a);
      SphericalCoordinates P(a - this->X);
      std::array<std::complex<double>, 2> ret = {0, 0};
      std::complex<double> R;
      for (const auto& [n, m] : this->nm_set) {
         // auto R = P.SolidHarmonicR(n, -m);
         // auto R = std::pow(P.rho, n) * P.sph_harmonics(n, m);
         R = P.sph_harmonics_rho(n, m);
         complex_fused_multiply_increment(R, this->get_coeffs(n, m) /*L*/, std::get<0>(ret));
         complex_fused_multiply_increment(R, this->get_coeffs_(n, m) /*L*/, std::get<1>(ret));
      }
      return {std::get<0>(ret).real(), std::get<1>(ret).real()};
   }

   /* -------------------------------------------------------------------------- */

   /*DOC_EXTRACT coeffs

   ### 定数の読み込み

   `AAA_M2M_FMM`，`AAA_M2L_FMM`，`AAA_L2L_FMM`は．`定数.nb`あらかじめ計算しておいたものを読み込む．

   */

   void M2M(const ExpCoeffs<N>& M) {
      // SphericalCoordinates r(M_shifted.X - M.X);
      SphericalCoordinates rab(M.X - this->X);
      rab.precompute_sph(2 * N);
      std::complex<double> R, c = {1., 0.}, AAA;
      double rho;
      this->increment_coeffs([&](int j, int k) -> std::array<std::complex<double>, 2> {
         std::array<std::complex<double>, 2> ret = {0, 0};
         auto& AAA_M2M_FMM_j_k = AAA_M2M_FMM[j][k + N_AAA_M2M_FMM];
         for (int n = 0; n <= j; ++n) {
            // rho = std::pow(rab.rho, n);
            for (int m = -n; m <= n; ++m) {
               AAA = AAA_M2M_FMM_j_k[n][m + N_AAA_M2M_FMM];
               if (AAA.real() != 0.0 || AAA.imag() != 0.0) {
                  // R = AAA * rho * rab.sph_harmonics(n, -m);
                  R = AAA * rab.sph_harmonics_rho(n, -m);
                  complex_fused_multiply_increment(R, M.get_coeffs(j - n, k - m), std::get<0>(ret));
                  complex_fused_multiply_increment(R, M.get_coeffs_(j - n, k - m), std::get<1>(ret));
               }
            }
         }
         return ret;
      });
   }

   void M2L(const ExpCoeffs<N>& M) {
      std::complex<double> S_, AAA;
      // SphericalCoordinates r(LocalExp.X - M.X);
      SphericalCoordinates rab(M.X - this->X);
      rab.precompute_sph(2 * N);

      this->increment_coeffs([&](int j, int k) -> std::array<std::complex<double>, 2> {
         std::array<std::complex<double>, 2> ret = {0, 0};
         auto& AAA_M2L_FMM_j_k = AAA_M2L_FMM[j][k + N_AAA_M2L_FMM];
         double rho;
         for (int n = 0; n <= N; ++n) {
            rho = std::pow(rab.rho, j + n + 1);
            for (int m = -n; m <= n; ++m) {
               AAA = AAA_M2L_FMM_j_k[n][m + N_AAA_M2L_FMM];
               if (AAA.real() != 0.0 || AAA.imag() != 0.0) {
                  S_ = AAA * rab.sph_harmonics(j + n, m - k) / rho;
                  complex_fused_multiply_increment(S_, M.get_coeffs(n, m), std::get<0>(ret));
                  complex_fused_multiply_increment(S_, M.get_coeffs_(n, m), std::get<1>(ret));
               }
            }
         }
         return ret;
      });

      // for (const auto& [j, k, n, m, AAA] : this->jknmAAA_set_for_M2L) {
      //    S_ = AAA * rab.sph_harmonics(j + n, m - k) / std::pow(rab.rho, j + n + 1);
      //    complex_fused_multiply_increment(S_, M.get_coeffs(n, m), this->coeffs[index(j, k)]);
      //    complex_fused_multiply_increment(S_, M.get_coeffs_(n, m), this->coeffs_[index(j, k)]);
      // }
   }

   void L2L(const ExpCoeffs<N>& L) {
      std::complex<double> R, AAA;
      double rho;
      // SphericalCoordinates r(LocalExp.X - L.X);
      SphericalCoordinates rab(L.X - this->X);
      rab.precompute_sph(2 * N);
      this->increment_coeffs([&](int j, int k) -> std::array<std::complex<double>, 2> {
         std::array<std::complex<double>, 2> ret = {0, 0};
         auto& AAA_L2L_FMM_j_k = AAA_L2L_FMM[j][k + N_AAA_L2L_FMM];
         for (int n = j; n <= N; ++n) {
            // rho = std::pow(rab.rho, n - j);
            for (int m = -n; m <= n; ++m) {
               AAA = AAA_L2L_FMM_j_k[n][m + N_AAA_L2L_FMM];
               if (AAA.real() != 0.0 || AAA.imag() != 0.0) {
                  // R = AAA * rab.sph_harmonics(n - j, m - k) * rho;
                  R = AAA * rab.sph_harmonics_rho(n - j, m - k);
                  complex_fused_multiply_increment(R, L.get_coeffs(n, m), std::get<0>(ret));
                  complex_fused_multiply_increment(R, L.get_coeffs_(n, m), std::get<1>(ret));
               }
            }
         }
         return ret;
      });
   }
};

/* -------------------------------------------------------------------------- */

double G(const std::array<double, 3>& X_near, const std::array<double, 3>& X_far) { return 1 / Norm(X_near - X_far); }

std::array<double, 3> gradG(const std::array<double, 3>& X_near, const std::array<double, 3>& X_far) {
   return -(X_near - X_far) / std::pow(Norm(X_near - X_far), 3.);
}

double Gapx(unsigned p,
            std::array<double, 3> X_near,
            std::array<double, 3> X_far,
            const std::array<double, 3>& center) {
   SphericalCoordinates sph_near(X_near - center);
   SphericalCoordinates sph_far(X_far - center);

   std::complex<double> accum = 0;
   int k, m;
   for (k = 0; k <= p; ++k)
      for (m = -k; m <= k; ++m)
         accum += sph_near.SolidHarmonicR(k, m) * sph_far.SolidHarmonicS(k, m);
   return accum.real();
}

std::array<double, 3> gradGapx(unsigned p,
                               const std::array<double, 3>& X_near_IN,
                               const std::array<double, 3>& X_far_IN,
                               const std::array<double, 3>& center) {
   SphericalCoordinates sph_near(X_near_IN - center);
   SphericalCoordinates sph_far(X_far_IN - center);
   int k, m;

   std::array<std::complex<double>, 3> grad = {0., 0., 0.};

   for (k = 0; k <= p; ++k)
      for (m = -k; m <= k; ++m)
         grad += sph_near.Grad_SolidHarmonicR(k, m) * sph_far.SolidHarmonicS(k, m);

   return {grad[0].real(), grad[1].real(), grad[2].real()};
}

std::array<double, 3> ToSphericalCoordinates(const double x, const double y, const double z) {
   return {std::hypot(x, y, z), std::atan2(std::hypot(x, y), z), std::atan2(y, x)};
};

// \label{ToSphericalCoordinates}
std::array<double, 3> ToSphericalCoordinates(const std::array<double, 3>& X) {
   return ToSphericalCoordinates(std::get<0>(X), std::get<1>(X), std::get<2>(X));
};

#endif