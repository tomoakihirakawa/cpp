#ifndef networkLine_H
#define networkLine_H
#pragma once

#include "Network.hpp"

inline T2Tddd networkLine::getLocationsTuple() const {
   return {this->Point_A->getXtuple(), this->Point_B->getXtuple()};
};
// inline bool networkLine::setBounds()
// {
// 	CoordinateBounds::setBounds(getLocationsTuple());
// 	for (auto &f : this->getFaces())
// 		if (f)
// 		{
// 			f->setPointsFromLines();
// 			f->setBounds();
// 		}
// 	return true;
// };
inline void networkLine::setBoundsSingle() {
   CoordinateBounds::setBounds(getLocationsTuple());
};
inline bool networkLine::Replace(netP *oldP, netP *newP) {
   auto bool1 = this->Switch(oldP, newP);  // 1
   auto bool2 = oldP->Erase(this);         // 2
   auto bool3 = newP->Add(this);           // 3
   // このステップがdouble replace
   //  switchでないと，順番に意味のあるFaceではおかしくなるので注意
   if (bool1 && bool2 & bool3)
      return true;
   else
      return false;
};
// inline bool networkLine::Replace(netF *oldF, netF *newF, netL *newL)
// {
// 	auto bool1 = this->Switch(oldF, newF); // 1
// 	auto bool2 = oldF->Erase(this);		   // 2
// 	auto bool3 = newF->Add(this);		   // 3
// 										   //このステップがdouble replace
// 										   // switchでないと，順番に意味のあるFaceではおかしくなるので注意

// 	if (newL != nullptr)
// 	{
// 		auto bool4 = oldF->Add(newL);
// 		auto bool5 = newL->Add(oldF); //許されない，Pointの場合

// 		if (bool1 && bool2 & bool3 && bool4 & bool5)
// 			return true;
// 		else
// 			return false;
// 	}

// 	if (bool1 && bool2 & bool3)
// 		return true;
// 	else
// 		return false;
// 	// throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
// 	// return true;
// 	// return network::doubleReplace(this, oldL, newL, newF);
// 	// return network::doubleReplace(this, oldF, newF, newL);
// };

inline networkLine::networkLine(Network *network_IN,
                                netP *sPoint_IN,
                                netP *ePoint_IN)
    : CoordinateBounds({0, 0, 0}),
      status(false),
      Faces(0),
      network(network_IN),
      Point_A(nullptr),
      Point_B(nullptr),
      Neumann(false),
      Dirichlet(false),
      CORNER(false) {
#ifdef DEM
   this->tension = 0.;
#endif
   network->Lines.emplace(this);
   set(sPoint_IN, ePoint_IN);
   sPoint_IN->Add(this);
   ePoint_IN->Add(this);
   setBoundsSingle();
   // setBounds();
};

inline bool networkLine::isIntxn() {
   //  return intxn; /*face-face intersection*/
   auto fs = this->Faces;
   if (fs.size() < 2)
      return false;
   else {
      for (auto i = 0; i < fs.size(); i++)
         for (auto j = i + 1; j < fs.size(); j++)
            if (fs[i]->getNetwork() != fs[j]->getNetwork())
               return true;
      return false;
   }
};

// networkLine
inline double networkLine::length() const {
   try {
      return Norm(this->Point_A->getXtuple() - this->Point_B->getXtuple());
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
//
// inline V_d networkLine::getNormal() const
// {
// 	VV_d normals;
// 	for (const auto &f : this->Faces)
// 		normals.emplace_back(ToVector(f->normal));
// 	return Mean(normals);
// };

inline Tddd networkLine::getNormal() const {
   Tddd ret = {0, 0, 0};
   for (const auto &f : this->Faces)
      ret += f->normal;
   return ret / (double)(this->Faces.size());
};

inline V_netFp networkLine::getFacesPenetrating() const {
   V_netFp ret({});
   netFp f;
   for (const auto &p : this->XPoints)
      if ((f = p->getXFace()))
         ret.emplace_back(f);
   return ret;
};

// class boundsSetter
// {
// public:
// 	V_netPp Points;
// 	V_netFp Faces;
// 	V_netLp Lines;
// 	boundsSetter() : Points({}), Lines({}), Faces({}){};
// 	//面や点がまだ存在するかチェックする昨日が必要だろう
// 	void add(const V_netPp &ps)
// 	{
// 		network::add(this->Points, ps);
// 		for (const auto &p : ps)
// 		{
// 			this->add(p->getFaces());
// 			this->add(p->getLines()); //その時点のlineを保存しておくことが大事
// 		}
// 	};
// 	void add(const V_netFp &fs) { network::add(this->Faces, fs); };
// 	void add(const V_netLp &ls) { network::add(this->Lines, ls); };
// 	void setBounds()
// 	{
// 		for (const auto &p : this->Points)
// 			p->setBounds();
// 		for (const auto &f : this->Faces)
// 		{
// 			f->setPointsFromLines();
// 			f->setBounds();
// 		}
// 		for (const auto &l : this->Lines)
// 			l->setBounds();

// 		// Print("setBounds done");
// 	};
// };

bool isLinkedDoubly(const netLp l, const netPp p) {
   try {
      auto [q0, q1] = l->getPoints();
      if (q0 && q0 == p) {
         for (const auto &m : p->getLines())
            if (m && m == l)
               return true;
      }
      if (q1 && q1 == p) {
         for (const auto &m : p->getLines())
            if (m && m == l)
               return true;
      }
      return false;
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};

bool isLinkedDoubly(const netLp l, const netFp f) {
   try {
      for (const auto &q : l->getFaces())
         if (q && q == f) {
            auto [l0, l1, l2] = f->getLines();
            if ((l0 && l0 == l) || (l1 && l1 == l) || (l2 && l2 == l))
               return true;
         }
      return false;
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};

bool isLinkedDoubly(const netPp p, const netLp l) { return isLinkedDoubly(l, p); };
bool isLinkedDoubly(const netFp f, const netLp l) { return isLinkedDoubly(l, f); };

// bool isConsistent(const networkLine *const l)
// {
// 	std::cout << Blue << "isConsistent";

// 	//２点を結んでいるかどうかのチェック
// 	if (l->getPoints().size() != 2)
// 	{
// 		std::stringstream ss;
// 		ss << "l->getPoints() = " << l->getPoints();
// 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
// 	};

// 	for (const auto &f : l->getFaces())
// 	{
// 		auto fps = f->getPoints(l);

// 		if (fps.size() != 3)
// 			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "fps.size() != 3");
// 		// line->face <=== consistent ===> line
// 		//このlineを参照できるか
// 		if (!MemberQ(fps[0]->getLines(), l))
// 		{
// 			std::cout << " l = " << l << std::endl;
// 			std::cout << " fps[0]->getLines() = " << fps[0]->getLines() << std::endl;
// 			std::cout << " l->getPoints() = " << l->getPoints() << std::endl;
// 			for (const auto &L : fps[0]->getLines())
// 				std::cout << " fps[]->getLines()->getPoints() = " << L->getPoints() << std::endl;
// 			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "0");
// 		}
// 		if (!MemberQ(fps[1]->getLines(), l))
// 		{
// 			std::cout << " l = " << l << std::endl;
// 			std::cout << " fps[1]->getLines() = " << fps[1]->getLines() << std::endl;
// 			std::cout << " l->getPoints() = " << l->getPoints() << std::endl;
// 			for (const auto &L : fps[1]->getLines())
// 				std::cout << " fps[]->getLines()->getPoints() = " << L->getPoints() << std::endl;
// 			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "1");
// 			return false;
// 		}
// 		if (MemberQ(fps[2]->getLines(), l))
// 		{
// 			std::cout << " l = " << l << std::endl;
// 			std::cout << " fps[2]->getLines() = " << fps[2]->getLines() << std::endl;
// 			std::cout << " l->getPoints() = " << l->getPoints() << std::endl;
// 			for (const auto &L : fps[2]->getLines())
// 				std::cout << " fps[]->getLines()->getPoints() = " << L->getPoints() << std::endl;

// 			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "2");
// 			return false;
// 		}

// 		// line->face->getPoints <=== consistent ===> line->getPoints
// 		for (const auto &p : l->getPoints())
// 			if (!MemberQ(f->getPoints(), p))
// 			{
// 				Print(p);
// 				Print(f->getPoints());
// 				throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "line->face->getPoints <=== consistent ===> line->getPoints");
// 				return false;
// 			}
// 	}
// 	std::cout << Green << " done" << colorOff << std::endl;
// 	return true;
// };
//% ------------------------------------------------------ */
//%                          辺の分割                        */
//% ------------------------------------------------------ */
// #define debug_divide

inline netPp networkLine::divide(const Tddd &midX) {
   // isConsistent(this);
   // 2面の場合も対応できるように，ポインター予め準備しておく
   netFp oldF = nullptr;
   V_netLp oldFLines;
   netLp bL = nullptr, fL = nullptr;
   //
   netFp newF = nullptr;
   V_netPp oldF_ps;
   netPp fP = nullptr, oP = nullptr, bP = nullptr;
   //
   netFp oldF1 = nullptr;
   V_netPp oldF1_ps;
   netPp fP1 = nullptr, oP1 = nullptr, bP1 = nullptr;
   //
   netFp newF1 = nullptr;
   V_netLp oldFLines1;
   netLp bL1 = nullptr, fL1 = nullptr;

   int c = 0;  // カラーバー
   try {
      if (this->Faces.size() == 1 || this->Faces.size() == 2) {
         // std::cout << "this->Faces = " << this->Faces << std::endl;
         if (this->Faces.size() > 0) {
            // 点・線・面のベクトルのそれぞれの設定を忘れずに
            oldF = this->Faces[0];
            bL = oldF->getLine(this, -1);
            fL = oldF->getLine(this, 1);
            auto [a, b, c] = oldF->getPoints(this);
            fP = b;
            oP = c;
            bP = a;
            newF = new networkFace(oldF);
         }

         if (this->Faces.size() == 2) {
            oldF1 = this->Faces[1];
            bL1 = oldF1->getLine(this, -1);
            fL1 = oldF1->getLine(this, 1);

            // oldF1_ps = oldF1->getPoints(this);
            // fP1 = oldF1_ps[1];
            // oP1 = oldF1_ps[2];
            // bP1 = oldF1_ps[0];

            auto [a, b, c] = oldF1->getPoints(this);
            fP1 = b;
            oP1 = c;
            bP1 = a;

            newF1 = new networkFace(oldF1);
         }
         ///////////////////////////////////////////////////////
         // boundsSetter bSetter;
         // bSetter.add({fP, oP, bP});
         ///////////////////////////////////////////////////////
         netPp newP = nullptr;
         netLp newDivL, newMidL;
         try {
            // std::cout << ColorFunction(c++) << "|" << colorOff;
            newP = new networkPoint(this->getNetwork() /*属性*/, isFinite(midX) ? midX : (fP->getXtuple() + bP->getXtuple()) / 2.);
            newDivL = new networkLine(this->getNetwork(), newP, fP);
#ifdef debug_divide
            std::cout << "debug divide, " << __FILE__ << ", " << __PRETTY_FUNCTION__ << ", " << __LINE__ << std::endl;
#endif
            /*           /     \
             *          /       \
             *         /         \
             *        /    / \    \              /   / \   \
             *     bL/-><-/   \-><-\fL        bL/ <-/   \-> \fL
             *      /    /oldF \    \          /   /newF \   \
             *     /     ---V---     \             ---|---
             *    /         A         \               V
             *  bP-><------this-----><-fP        ----this----
             *                              newP-><-------newDivL-----><-fP
             */
            if (!(this->Switch(fP, newP)))
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
            newP->Add(this);
#ifdef debug_divide
            std::cout << "debug divide, " << __FILE__ << ", " << __PRETTY_FUNCTION__ << ", " << __LINE__ << std::endl;
#endif
            /*
             *           /     \
             *          /       \
             *         /         \
             *        /    / \    \              /   / \   \
             *     bL/-><-/   \-><-\fL        bL/ <-/   \-> \fL
             *      /    /oldF \    \          /   /newF \   \
             *     /     ---V---     \             ---|---
             *    /         A         \               V
             *  bP-><------this      <-fP        ----this----
             *                \-------><-newP-><-------newDivL-----><-fP
             */
            if (!(fP->Erase(this)))
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
#ifdef debug_divide
            std::cout << "debug divide, " << __FILE__ << ", " << __PRETTY_FUNCTION__ << ", " << __LINE__ << std::endl;
#endif
            /*              oP                               oP
             *           V     V                              V
             *          A       A                              A
             *         /         \                              \
             *        /    / \    \              /   / \   \     \
             *     bL/-><-/   \-><-\fL        bL/ <-/   \-> \fL   \
             *      /    /oldF \    \          /   /newF \   \     \
             *     V     ---V---     V             ---|---          \
             *    A         A         A               V              V
             *  bP-><------this        fP        ----this----         A
             *                \-------><-newP-><-------newDivL-----><-fP
             */
            // std::cout << ColorFunction(c++) << "|" << colorOff;
            newMidL = new networkLine(this->getNetwork(), newP, oP);
#ifdef debug_divide
            std::cout << "debug divide, " << __FILE__ << ", " << __PRETTY_FUNCTION__ << ", " << __LINE__ << std::endl;
#endif
            /*             oP ------><--------                     oP
             *           V     V             |                      V
             *          A       A            |                       A
             *         /         \           |                        \
             *        /    / \    \          |         /   / \   \     \
             *     bL/-><-/   \-><-\fL    newMidL   bL/ <-/   \-> \fL   \
             *      /    /oldF \    \        |       /   /newF \   \     \
             *     V     ---V---     V       V           ---|---          \
             *    A         A        A       A              V             V
             *  bP-><------this      fP      |         ----this----       A
             *                 \--------><--newP-><-------newDivL-----><-fP*/
            newMidL->set(oldF, newF);
#ifdef debug_divide
            std::cout << "debug divide, " << __FILE__ << ", " << __PRETTY_FUNCTION__ << ", " << __LINE__ << std::endl;
#endif
            /*             oP ------><--------                     oP
             *           V     V             |                      V
             *          A       A            |                       A
             *         /         \           |                        \
             *        /    / \    \          |         /   / \   \     \
             *     bL/-><-/   \-><-\fL       |      bL/ <-/   \-> \fL   \
             *      /    /oldF \ <--------newMidL------->/newF \   \     \
             *     V     ---V---     V       V           ---|---          \
             *    A         A        A       A              V             V
             *  bP-><------this      fP      |         ----this----       A
             *                 \--------><--newP-><-------newDivL-----><-fP
             */
            if (!(oldF->Switch(fL, newMidL))) {
               std::stringstream ss;
               // ss << "oldF->getLines() = " << oldF->getLines();
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
            }
#ifdef debug_divide
            std::cout << "debug divide, " << __FILE__ << ", " << __PRETTY_FUNCTION__ << ", " << __LINE__ << std::endl;
#endif
            /*             oP ------><--------                     oP
             *           V                   |                      V
             *          A                    |                       A
             *         /                     |                        \
             *        /    / \               |         /   / \   \     \
             *     bL/-><-/   \------><---newMidL   bL/ <-/   \-> \fL   \
             *      /    /oldF \<--fL        |   \------>/newF \   \     \
             *     V     ---V---             V           ---|---          \
             *    A         A                A              V             V
             *  bP-><------this              |         ----this----       A
             *                 \--------><--newP-><-------newDivL-----><-fP
             * std::cout << ColorFunction(c++) << "|" << colorOff;
             */
            if (!(fL->Switch(oldF, newF))) {
               std::stringstream ss;
               ss << "fL->getFaces() = " << fL->getFaces();
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
            }
#ifdef debug_divide
            std::cout << "debug divide, " << __FILE__ << ", " << __PRETTY_FUNCTION__ << ", " << __LINE__ << std::endl;
#endif
            /*
             *             oP ------><--------                      oP
             *           V                   |                       V
             *          A                    |                        A
             *         /                     |                         \
             *        /    / \               |          /  / \    \     \
             *     bL/-><-/   \------><---newMidL    bL/<-/   \-><-\fL   \
             *      /    /oldF \<--fL        |   \------>/newF \   \      \
             *     V     ---V---             V           ---|---           \
             *    A         A                A              V              V
             *  bP-><------this              |         ----this----        A
             *                 \--------><--newP-><-------newDivL-----><-fP*/
            if (!(newF->Switch(bL, newMidL)))
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
#ifdef debug_divide
            std::cout << "debug divide, " << __FILE__ << ", " << __PRETTY_FUNCTION__ << ", " << __LINE__ << std::endl;
#endif
            /*             oP ------><--------                      oP
             *           V                   |                       V
             *          A                    |                        A
             *         /                     |                         \
             *        /    / \               |             / \    \     \
             *     bL/-><-/   \------><---newMidL---><----/   \-><-\fL   \
             *      /    /oldF \　　　        |           /newF \    \     \
             *     V     ---V---             V           ---|---           \
             *    A         A                A              V              V
             *  bP-><------this              |         ----this----        A
             *                 \--------><--newP-><-------newDivL-----><-fP
             *
             */
            // std::cout << ColorFunction(c++) << "|" << colorOff;

            if (!(newF->Switch(this, newDivL))) {
               std::stringstream ss;
               // ss << "newF->getLines() = " << newF->getLines() << ", this = " << this;
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
            }
#ifdef debug_divide
            std::cout << "debug divide, " << __FILE__ << ", " << __PRETTY_FUNCTION__ << ", " << __LINE__ << std::endl;
#endif

            newDivL->Add(newF);

            /*            oP ------><--------                      oP
             *           V                   |                       V
             *          A                    |                        A
             *         /                     |                         \
             *        /    / \               |             / \    \     \
             *     bL/-><-/   \------><-----newMidL--><---/   \-><-\fL   \
             *      /    /oldF \　　　        |           /newF \    \     \
             *     V     ---V---             V           ---|---           \
             *    A         A                A              V              V
             *  bP-><------this              |              A              A
             *                 \--------><--newP-><-------newDivL-----><-fP
             */
            // check
         } catch (std::exception &e) {
            std::cerr << e.what() << colorOff << std::endl;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         };
#ifdef debug_divide
         std::cout << "debug divide, " << __FILE__ << ", " << __PRETTY_FUNCTION__ << ", " << __LINE__ << std::endl;
#endif

         // もし2面ある場合
         if (newF1 != nullptr) {
#if defined(debug_divide)
            std::cout << Blue << "|";
#endif
            // std::cout << ColorFunction(c++) << "|" << colorOff;
            auto newMidL1 = new networkLine(this->getNetwork(), newP, oP1);
            // 逆側も
            /*             oP ------><--------                      oP
             *           V                   |                       V
             *          A                    |                        A
             *         /                     |                         \
             *        /    / \               |             / \    \     \
             *     bL/-><-/   \------><-----newMidL--><---/   \-><-\fL   \
             *      /    /oldF \　　　        |           /newF \    \     \
             *     V     ---V---             V           ---|---           \
             *    A         A                A              V              V
             *  bP-><------this      fP      |              A              A
             *   V          |  \-----+--><--newP-><-------newDivL-----><-fP
             *    A         V        V       V                this
             *     \     ---A---     A       A              ---A---
             *   fL1\-><-\oldF1/    /        |         fL1<-\newF1/
             *       \    \   /-><-/bL1      |               \   /->bL1
             *        \    \ /    /        newMidL1           \ /
             *          V       V            |
             *           A     A             |
             *             oP1 ------><-------
             */
            // std::cout << ColorFunction(c++) << "|" << colorOff;
#if defined(debug_divide)
            std::cout << red << "|";
#endif
            if (!(oldF1->Switch(bL1, newMidL1)))
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
               /*             oP ------><--------                      oP
                *           V                   |                       V
                *          A                    |                        A
                *         /                     |                         \
                *        /    / \               |             / \    \     \
                *     bL/-><-/   \------><-----newMidL--><---/   \-><-\fL   \
                *      /    /oldF \　　　        |           /newF \    \     \
                *     V     ---V---             V           ---|---           \
                *    A         A                A              V              V
                *  bP-><------this      fP      |              A              A
                *   V          |  \-----+--><--newP-><-------newDivL-----><-fP
                *    A         V        V       V                this
                *     \     ---A---     A       A              ---A---
                *   fL1\-><-\oldF1/    /        |         fL1<-\newF1/
                *       \    \   /  <-/bL1      |               \   /->bL1
                *        \    \ / ---------->newMidL1            \ /
                *          V       V            |
                *           A     A             |
                *             oP1 ------><-------
                */
               // std::cout << ColorFunction(c++) << "|" << colorOff;
#if defined(debug_divide)
            std::cout << red << "|";
#endif
            if (!(bL1->Switch(oldF1, newF1)))
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
               /*             oP ------><--------                      oP
                *           V                   |                       V
                *          A                    |                        A
                *         /                     |                         \
                *        /    / \               |             / \    \     \
                *     bL/-><-/   \------><-----newMidL--><---/   \-><-\fL   \
                *      /    /oldF \　　　        |           /newF \    \     \
                *     V     ---V---             V           ---|---           \
                *    A         A                A              V              V
                *  bP-><------this              |              A              A
                *   V          |  \--------><--newP-><-------newDivL-----><-fP
                *    A         V                V                this       V
                *     \     ---A---             A              ---A---      A
                *   fL1\-><-\oldF1/             |         fL1<-\newF1/     /
                *       \    \   /              |               \   /-><-/bL1
                *        \    \ / ---------->newMidL1            \ /    /
                *          V                    |                      /
                *           A                   |                     /
                *             oP1 ------><-------
                */
               // std::cout << ColorFunction(c++) << "|" << colorOff;
#if defined(debug_divide)
            std::cout << red << "|";
#endif
            if (!(newF1->Switch(fL1, newMidL1)))
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
            newDivL->Add(newF1);
            /*           oP ------><--------                      oP
             *           V                   |                       V
             *          A                    |                        A
             *         /                     |                         \
             *        /    / \               |             / \    \     \
             *     bL/-><-/   \------><-----newMidL--><---/   \-><-\fL   \
             *      /    /oldF \　　　        |           /newF \    \     \
             *     V     ---V---             V           ---|---           \
             *    A         A                A              V              V
             *  bP-><------this              |              A              A
             *   V          |  \--------><--newP-><-------newDivL-----><-fP
             *    A         V                V               | this      V
             *     \     ---A---             A              -V--A--      A
             *   fL1\-><-\oldF1/---------->newMidL1<--------\newF1/    /
             *       \    \   /              |               \   /-><-/bL1
             *        \    \ /               |                \ /    /
             *          V                    |                      /
             *           A                   |                     /
             *             oP1 ------><-------
             */
            // std::cout << ColorFunction(c++) << "|" << colorOff;
#if defined(debug_divide)
            std::cout << red << "|";
#endif
            newMidL1->set(oldF1, newF1);
            if (!(newF1->Switch(this, newDivL))) {
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
            }
            /*           oP ------><--------                      oP
             *           V                   |                       V
             *          A                    |                        A
             *         /                     |                         \
             *        /    / \               |             / \    \     \
             *     bL/-><-/   \------><-----newMidL--><---/   \-><-\fL   \
             *      /    /oldF \　　　        |           /newF \    \     \
             *     V     ---V---             V           ---|---           \
             *    A         A                A              V              V
             *  bP-><------this              |              A              A
             *   V          |  \--------><--newP-><-------newDivL-----><-fP
             *    A         V                V                 V         V
             *     \     ---A---             A              ---A---      A
             *   fL1\-><-\oldF1/-----><---newMidL1----><----\newF1/    /
             *       \    \   /              |               \   /-><-/bL1
             *        \    \ /               |                \ /    /
             *          V                    |                      /
             *           A                   |                     /
             *             oP1 ------><-------
             */
            // std::cout << ColorFunction(c++) << "|" << colorOff;
#if defined(debug_divide)
            std::cout << Blue << "|" << std::endl;
#endif
            newP->setBoundsSingle();
            bP->setBoundsSingle();
            fP->setBoundsSingle();
            //
            bL->setBoundsSingle();
            newMidL->setBoundsSingle();
            fL->setBoundsSingle();
            fL1->setBoundsSingle();
            newMidL1->setBoundsSingle();
            bL1->setBoundsSingle();
            this->setBoundsSingle();
            newDivL->setBoundsSingle();
            //
            newP->setBoundsSingle();
            newP->setFaces(newP->Lines);
            for (const auto &f : newP->getFaces())
               f->setGeometricProperties(ToX(f->setPoints(f->Lines)));
            bP->setBoundsSingle();
            bP->setFaces(bP->Lines);
            for (const auto &f : bP->getFaces())
               f->setGeometricProperties(ToX(f->setPoints(f->Lines)));
            fP->setBoundsSingle();
            fP->setFaces(fP->Lines);
            for (const auto &f : fP->getFaces())
               f->setGeometricProperties(ToX(f->setPoints(f->Lines)));
            //
            newP->Dirichlet = this->Dirichlet;
            newP->Neumann = this->Neumann;
            newP->CORNER = this->CORNER;
         }
#ifdef debug_divide
         std::cout << Blue << "|" << colorOff << std::endl;
         std::cout << "debug divide, " << __FILE__ << ", " << __PRETTY_FUNCTION__ << ", " << __LINE__ << std::endl;
#endif
         return newP;
      } else {
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
                             "この線は面を持たないか何かの理由で，divideできません\n this->Faces.size()=" + std::to_string(this->Faces.size()));
      }
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };

   // Print("divide done");
};
//@ ------------------------------------------------------ */
//@                          辺の削除                        */
//@ ------------------------------------------------------ */
/*networkLine::divide_code*/
/*networkLine::divide_code*/
inline bool networkLine::isMergeable() const {
   // isConsistent(this);
   // auto oldPoints = this->getPoints(); //最後にsetBoundsする必要がる
   //* ------------------------------------------------------ */
   auto [p0, p1] = this->getPoints();
   if (p0->getNeighbors().size() < 4 || p1->getNeighbors().size() < 4) {
      Print("今の所，mergeできない条件");
      return false;
   }

   auto AB = this->getFaces();
   if (AB.size() != 2) {
      Print("AB.size() != 2");
      return false;
   }
   auto A = AB[0];
   auto B = AB[1];
   auto Aps = A->getPoints(this);
   auto Bps = B->getPoints(this);
   //* ------------------------------------------------------ */
   int c = 0;
   for_each(Join(Bps, Aps), [&](const auto &p) {		if (p->getLines().empty())
		{
			throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p->getLines().empty()");
		} });

   /*ここが間違っていたconcaveの場合，拒否する*/
   if (isConcavePolygon({std::get<1>(Aps)->getXtuple(),
                         std::get<2>(Aps)->getXtuple(),
                         std::get<1>(Bps)->getXtuple(),
                         std::get<2>(Bps)->getXtuple()})) {
      Print("Concave Polygonです．mergeしません．");
      return false;
   }
   //* ------------------------------------------------------ */
   auto Abl = A->getLineBack(this);
   auto Afl = A->getLineFront(this);
   auto Bbl = B->getLineBack(this);
   auto Bfl = B->getLineFront(this);

   auto a = (*Abl)(A);
   auto b = (*Bfl)(B);
   if (MemberQ(std::get<1>(Aps)->getFaces(), a) || MemberQ(std::get<1>(Aps)->getFaces(), b)) {
      std::cout << "面が潰れる条件なので，mergeしない" << std::endl;
      return false;
   }
   //* ------------------------------------------------------ */
   auto Abl_F_Except_A = TakeExcept(Abl->getFaces(), A);
   auto Bfl_F_Except_B = TakeExcept(Bfl->getFaces(), B);
   if (Abl_F_Except_A.empty()) {
      Print("Abl_F_Except_A.empty()");
      if (Bfl_F_Except_B.empty()) {
         Print("Bfl_F_Except_B.empty()");
         return false;
      }
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   }
   return true;
};
//@ ------------------------------------------------------ */
//@ ------------------------------------------------------ */
// #define debug_merge
inline netPp networkLine::merge() {
   auto AB = this->getFaces();
   auto A = AB[0];
   auto B = AB[1];
   auto Aps = A->getPoints(this);
   auto Bps = B->getPoints(this);
   //
   auto Abl = A->getLineBack(this);
   auto Afl = A->getLineFront(this);
   auto Bbl = B->getLineBack(this);
   auto Bfl = B->getLineFront(this);

   auto a = (*Abl)(A);
   auto b = (*Bfl)(B);

   auto aa = (*Afl)(A);
   auto bb = (*Bbl)(B);

   // aaとbbが同じ！？divideのミスでしょう．

   std::cout << "A " << A << std::endl;
   std::cout << "B " << B << std::endl;
   // std::cout << "a points " << a->getPoints(Aps[0]) << std::endl;
   // std::cout << "b points " << b->getPoints(Aps[0]) << std::endl;
   std::cout << "a " << a << std::endl;
   std::cout << "b " << b << std::endl;
   std::cout << "aa " << aa << std::endl;
   std::cout << "bb " << bb << std::endl;
   std::cout << "(*Afl)(A) " << (*Afl)(A) << std::endl;
   std::cout << "(*Bbl)(B) " << (*Bbl)(B) << std::endl;

   auto Aps0lines = TakeExcept(std::get<0>(Aps)->getLines(), {this, Abl, Bfl});
   try {
#if defined(debug_merge)
      std::cout << green << "|" << colorOff;
#endif

      /*                          Aps[2]
       *                          V    V
       *                         A      A
       *                        /         \
       *      (will be deleted)/    / \    \
       *                   Abl/-><-/   \-><-\Afl
       *                     /    /  FA \ 　  \　　
       *                    V     ---V---      V
       *  (will be deleted)A         A          A
       *        Bps[1],Aps[0]---><--this---><---Aps[1],Bps[0]
       *                  V          |          V
       *                   A         V         A
       *   (will be deleted)\     ---A---     /
       *                  Bfl\-><-\     /-><-/Bbl
       *                      \    \ FB/    /
       *                       \    \ /    /
       *                        V        V
       *                         A      A
       *                          Bps[2]
       */
#if defined(debug_merge)
      std::cout << green << "|" << colorOff;
#endif
      auto mid_X = (std::get<0>(Aps)->getXtuple() + std::get<1>(Aps)->getXtuple()) / 2.;

#if defined(debug_merge)
      std::cout << green << "|" << colorOff;
#endif

      /*                      Aps[2]
       *                     V    V
       *                     A      A
       *                   /         \          ______
       *                  /    / \    \         \    /
       *         a-><-Abl/-><-/   \-><-\Afl-><---\  /
       *                /    /  FA \ 　  \　　     \/
       *               V     ---V---      V
       *              A         A          A
       *  Bps[1],Aps[0]---><--this---><---Aps[1],Bps[0]
       *             V          |          V
       *              A         V         A
       *               \     ---A---     /
       *        b-><-Bfl\-><-\     /-><-/Bbl
       *                 \    \ FB/    /
       *                  \    \ /     /
       *                    V        V
       *                     A      A
       *                      Bps[2]
       */
#if defined(debug_merge)
      std::cout << green << "|" << colorOff;
#endif

      if (!(a->Switch(Abl, Afl)))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      if (!(Afl->Switch(A, a)))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      if (!(Abl->Erase(a)))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         // if (!(Abl->Erase(A)))
         // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         // if (!(A->Erase(Abl)))
         // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         // if (!(A->Erase(Afl)))
         // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         // if (!(A->Erase(this)))
         // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         // if (!(this->Erase(A)))
         // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         //
         // if (!(Abl->Erase(Aps[2])))
         // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         // if (!(Abl->Erase(Aps[0])))
         // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         // if (!(Aps[2]->Erase(Abl)))
         // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         // if (!(Aps[0]->Erase(Abl)))
         // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

         /*                 Aps[2]
          *                           V
          *     ______                 A
          *     \ a / ----><----------->\Afl       ______
          *      \ /             / \     \         \    /
          *               Abl   /   \     \----><---\  /
          *                    /  A  \  　  \　　     \/
          *                    -------      V
          *                                   A
          *  Bps[1],Aps[0]---><--this---><---Aps[1],Bps[0]
          *             V          |          V
          *              A         V         A
          *               \     ---A---     /
          *(*Bfl)(B)-><-Bfl\-><-\     /-><-/Bbl
          *                 \    \ B /    /
          *                  \    \ /    /
          *                   V         V
          *                    A       A
          *                      Bps[2]
          */

#if defined(debug_merge)
      std::cout << green << "|" << colorOff;
#endif

      if (!(b->Switch(Bfl, Bbl)))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      if (!(Bbl->Switch(B, b)))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      if (!(Bfl->Erase(b)))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // if (!(Bfl->Erase(B)))
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // if (!(B->Erase(Bfl)))
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // if (!(B->Erase(Bbl)))
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // if (!(B->Erase(this)))
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // if (!(this->Erase(B)))
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // //
      // if (!(Bfl->Erase(Bps[2])))
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // if (!(Bfl->Erase(Aps[0])))
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // if (!(Bps[2]->Erase(Bfl)))
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // if (!(Aps[0]->Erase(Bfl)))
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

      // Aps[0]->Erase(this);
      /*                   Aps[2]
       *                           V
       *     ______                 A
       *     \ a / ----><----------->\Afl       ______
       *      \ /             / \     \         \    /
       *               Abl   /   \     \----><---\  /
       *                    /  FA \  　  \　　     \/
       *                    -------      V
       *                                   A
       *  Bps[1],Aps[0]    <--this---><---Aps[1],Bps[0]
       *                                   V
       *                                  A
       *                     -------     /
       *              Bfl    \     /    /Bbl
       *                      \ B /    /
       *       /b \            \ /    /
       *      /____\----><---------->/Bbl
       *                            A
       *                      Bps[2]
       */

      // for (const auto &l : Aps[0]->getLines())
      // auto tmp = TakeExcept(Aps[0]->getLines(), {this, Abl, Bfl});
      for (auto &l : Aps0lines) {
         if (!(l->Switch(std::get<0>(Aps), std::get<1>(Aps))))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         // if (!(Aps[0]->Erase(l)))
         // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         std::get<1>(Aps)->Add(l);
         // Aps[1]が3lineしか持っていない場合もあり得るので，Addできない場合もある
         // throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      }
      // if (!Aps[0]->getLines().empty())
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // if (!this->Erase(Aps[0]))
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // if (!this->Erase(Aps[1]))
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // if (!Aps[1]->Erase(this))
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

      /*                   Aps[2]
       *                           V
       *     ______                 A
       *     \ a / ----><----------->\Afl       ______
       *      \ /             / \     \         \    /
       *               Abl   /   \     \----><---\  /
       *                    /  A  \  　  \　　     \/
       *                   --------       V           |
       *                                   A          V
       *  Bps[1],Aps[0]       this          Aps[1],Bps[0]<----
       *                                   V          A
       *                                  A           |
       *                    --------     /
       *              Bfl    \     /    /Bbl
       *                      \ B /    /
       *       /b \            \ /    /
       *      /____\----><---------->/Bbl
       *                            A
       *                      Bps[2]
       */

#if defined(debug_merge)
      if (!A->getLines().empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      if (!B->getLines().empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // if (!Abl->getPoints().empty())
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // if (!Bfl->getPoints().empty())
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      if (!Abl->getFaces().empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      if (!Bfl->getFaces().empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // if (!this->getPoints().empty())
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // if (!this->getFaces().empty())
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      if (!Aps[0]->getLines().empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

      /* ------------------------------------------------------ */
      if (Aps[2]->getLines().empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      if (Aps[1]->getLines().empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      if (Bps[2]->getLines().empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
#endif

      delete A;
      delete B;
      delete Abl;
      delete Bfl;
      delete this;
      delete std::get<0>(Aps);  // thisもきえます

      /*                   Aps[2]
       *                           V
       *     ______                 A
       *     \ a / ----><----------->\Afl       ______
       *      \ /                     \         \    /
       *                               \----><---\  /
       *                             　  \　　     \/
       *                                  V           |
       *                                   A          V
       *                                  Aps[1],Bps[0]<----
       *                                   V          A
       *                                   A           |
       *                                  /
       *                                 /Bbl
       *                                /
       *       /b \                    /
       *      /____\----><---------->/Bbl
       *                            A
       *                      Bps[2]
       */

#if defined(debug_merge)
      // lineを基準としてpointを取得できるか
      //  for (const auto &f : Afl->getFaces())
      //  	std::cout << "f->getPoints(Afl) = " << f->getPoints(Afl) << std::endl;

      // for (const auto &f : Bbl->getFaces())
      // 	std::cout << "f->getPoints(Bbl) = " << f->getPoints(Bbl) << std::endl;

      for (const auto &l : V_netLp{Afl, Bbl})
         for (const auto &f : l->getFaces())
            if (!isLinkedDoubly(f, l))
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

      for (const auto &p : V_netPp{std::get<1>(Aps), Aps[2], Bps[2]})
         for (const auto &l : p->getLines())
            if (!isLinkedDoubly(p, l))
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
#endif
      //
      // std::cout << std::get<1>(Aps)->getNetwork()->isMember(Aps[0]) << std::endl;
      // std::cout << std::get<1>(Aps)->getNetwork() << std::endl;
      // std::cout << std::get<1>(Aps)->getStorage() << std::endl;
      // std::cout << "check in merge 1" << std::endl;
      // std::get<1>(Aps)->getNetwork()->displayStates();
      // std::cout << Aps[0] << std::endl;
      // std::cout << "network : " << std::get<1>(Aps)->getNetwork() << std::endl;
      // std::cout << std::get<1>(Aps)->getNetwork()->isMember(Aps[0]) << std::endl;
      // std::cout << std::get<1>(Aps)->getStorage()->isMember(Aps[0]) << std::endl;
      // std::cout << std::get<1>(Aps)->getLines() << std::endl;
      // std::cout << Aps[2]->getLines() << std::endl;
      // std::cout << Bps[2]->getLines() << std::endl;
      // std::cout << "a->setBounds()" << std::endl;
      // a->setBounds();
      // std::cout << "b->setBounds()" << std::endl;
      // b->setBounds();
      // std::cout << "aa->setBounds()" << std::endl;
      // aa->setBounds();
      // std::cout << "bb->setBounds()" << std::endl;
      // bb->setBounds();
      // std::cout << "Aps[2]->setBounds();" << std::endl;
      // Aps[2]->setBounds();
      // std::cout << "Bps[2]->setBounds();" << std::endl;
      // Bps[2]->setBounds();
      // std::cout << "std::get<1>(Aps)->setX(mid_X);" << std::endl;
      std::get<1>(Aps)->setX(mid_X);
      return std::get<1>(Aps);
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      std::stringstream ss;
      ss << "Aps " << Aps << std::endl;
      ss << "Bps " << Bps << std::endl;
      ss << "Abl " << Abl << std::endl;
      ss << "Afl " << Afl << std::endl;
      ss << "Bbl " << Bbl << std::endl;
      ss << "Bfl " << Bfl << std::endl;
      ss << "a " << a << std::endl;
      ss << "b " << b << std::endl;
      ss << "aa " << aa << std::endl;
      ss << "bb " << bb << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
   };
};

inline netPp networkLine::mergeIfMergeable() {
   if (isMergeable())
      return this->merge();
   else
      return nullptr;
};

/*networkLine::divide_code*/

// #define debug_flip
inline bool networkLine::flip() {
#if defined(debug_flip)
   std::cout << Green << "flip |" << colorOff;
#endif
   try {
      // isConsistent(this);

      ///////////
      auto AB = this->getFaces();

      if (AB.size() != 2)
         return false;

      // std::cout << ColorFunction(c++) << "|" << colorOff;
      // auto oldPoints = this->getPoints(); //最後にsetBoundsする必要がる

      auto A = AB[0];
      auto B = AB[1];
      auto Aps = A->getPoints(this);
      auto Bps = B->getPoints(this);

      // std::cout << ColorFunction(c++) << "|" << colorOff;
      // for (const auto &p : Join(Bps, Aps))
      // 	if (p->getLines().empty())
      // 	{
      // 		// mk_vtu("./vtu/p->getLines().empty().vtu", {{p}});
      // 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p->getLines().empty()");
      // 	}

#if defined(debug_flip)
      std::cout << green << "|" << colorOff;
#endif

      /*ここが間違っていたconcaveの場合，拒否する*/
      // 凹凸包かどうかのチェックは，チェックできなかった場合falseを返す．チェックできなかった場合もflipさせないためには，!isConvexPolygonとしなければならない
      //  if (!isConvexPolygon({std::get<1>(Aps)->getX(), Aps[2]->getX(), Bps[1]->getX(), Bps[2]->getX()}))
      //  {
      //  	return false;
      //  }
      //  if (!isConvexPolygon({std::get<1>(Aps)->getX(), Aps[2]->getX(), Bps[1]->getX(), Bps[2]->getX()}))
      //  {
      //  	return false;
      //  }

      //  std::cout << ColorFunction(c++) << "|" << colorOff;
      //                     Aps[2]
      //                   V    V
      //                  A      A
      //                 /         \
		//                /    / \    \
		//            Abl/-><-/   \-><-\Afl
      //              /    /  FA \ 　  \　　
      //             V     ---V---      V
      //            A         A          A
      // Bps[1],Aps[0]---><--this---><---Aps[1],Bps[0]
      //           V          |          V
      //            A         V         A
      //             \     ---A---     /
      //           Bfl\-><-\     /-><-/Bbl
      //               \    \ FB/    /
      //                \    \ /    /
      //                  V        V
      //                   A      A
      //                    Bps[2]

#if defined(debug_flip)
      std::cout << green << "|" << colorOff;
#endif

      auto Abl = A->getLine(this, -1);
      auto Afl = A->getLine(this, 1);
      auto Bbl = B->getLine(this, -1);
      auto Bfl = B->getLine(this, 1);

      auto [p0, p1] = this->getPoints();

      // std::cout << ColorFunction(c++) << "|" << colorOff;
      if (!(std::get<2>(Aps)->Add(this)))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

      if (!(std::get<2>(Bps)->Add(this)))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

      if (!(p0->Erase(this)))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

      if (!(p1->Erase(this)))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

      this->set(std::get<2>(Aps), std::get<2>(Bps));  // setBoundsされないので注意

#if defined(debug_flip)
      std::cout << green << "|" << colorOff;
#endif

      //                --><-std::get<2>(Aps)
      //                |  V    V
      //                | A      A
      //                |/         \
        //                |    / \    \
        //            Abl/+><-/   \-><-\Afl
      //              / |  /  FA \ 　  \　　
      //             V  |  ---V---      V
      //            A   |      A          A
      // Bps[1],Aps[0]  ----this------    Aps[1],Bps[0]
      //           V          |      |   V
      //            A         V      |   A
      //             \     ---A---   |  /
      //           Bfl\-><-\     /-><|-/Bbl
      //               \    \ FB/    | /
      //                \    \ /    /|
      //                  V        V |
      //                   A      AV |
      //                  Bps[2]-> <-

      // std::cout << ColorFunction(c++) << "|" << colorOff;
#if defined(debug_flip)
      std::cout << green << "|" << colorOff;
#endif

      if (!(Abl->Switch(A, B)))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

      if (!(Bbl->Switch(B, A)))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

#if defined(debug_flip)
      std::cout << green << "|" << colorOff;
#endif

      //                --><-Aps[2]
      //                |  V    V
      //                | A      A
      //                |/         \
        //                |    / \    \
        //            Abl/  <-/   \-><-\Afl
      //              /||  /  FA \ 　  \　　
      //             V ||  ---V---      V
      //            A  ||     A   A----+  A
      // Bps[1],Aps[0] |----this------ |  Aps[1],Bps[0]
      //           V   |      |      | | V
      //            A   \     V      | | A
      //             \   V ---A---   | |/
      //           Bfl\-><-\     /-> | /Bbl
      //               \    \ FB/    |/
      //                \    \ /    /|
      //                  V        V |
      //                   A      AV |
      //                  Bps[2]-> <--

      // std::cout << ColorFunction(c++) << "|" << colorOff;
      A->setLines({Afl, this, Bbl});
      B->setLines({Bfl, this, Abl});
#if defined(debug_flip)
      std::cout << green << "|" << colorOff;
#endif

      //                --><-Aps[2]
      //                |  V    V
      //                | A      A
      //                |/         \
        //                |    / \    \
        //            Abl/|   /   \-><-\Afl
      //              /||  /  FA \ 　  \　　
      //             V ||  ---V---V      V
      //            A  ||     A   A----+  A
      // Bps[1],Aps[0] |----this------ |  Aps[1],Bps[0]
      //           V   |      |      | | V
      //            A   V     V      | | A
      //             \   A ---A---   | |/
      //           Bfl\-><-\     /   | /Bbl
      //               \    \ FB/    |/
      //                \    \ /    /|
      //                  V        V |
      //                   A      AV |
      //                  Bps[2]-> <--

#if defined(debug_flip)
      std::cout << "オブジェクトが互いに参照できる状態にあるか（リストにお互いに保存されているか）チェック" << std::endl;
      // line--><--face
      for (const auto &f : Afl->getFaces())
         if (!isLinkedDoubly(f, Afl))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      for (const auto &f : this->getFaces())
         if (!isLinkedDoubly(f, this))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      for (const auto &f : Bbl->getFaces())
         if (!isLinkedDoubly(f, Bbl))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      for (const auto &f : Abl->getFaces())
         if (!isLinkedDoubly(f, Abl))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      for (const auto &f : Bfl->getFaces())
         if (!isLinkedDoubly(f, Bfl))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

      // line--><--point
      // for (const auto &f : Afl->getPoints())
      // 	if (!isLinkedDoubly(f, Afl))
      // 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // for (const auto &f : this->getPoints())
      // 	if (!isLinkedDoubly(f, this))
      // 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // for (const auto &f : Bbl->getPoints())
      // 	if (!isLinkedDoubly(f, Bbl))
      // 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // for (const auto &f : Abl->getPoints())
      // 	if (!isLinkedDoubly(f, Abl))
      // 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      // for (const auto &f : Bfl->getPoints())
      // 	if (!isLinkedDoubly(f, Bfl))
      // 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

      // point--><--line
      for (const auto &l : Aps[0]->getLines())
         if (!isLinkedDoubly(l, Aps[0]))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      for (const auto &l : Aps[1]->getLines())
         if (!isLinkedDoubly(l, Aps[1]))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      for (const auto &l : Bps[2]->getLines())
         if (!isLinkedDoubly(l, Bps[2]))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      for (const auto &l : Aps[2]->getLines())
         if (!isLinkedDoubly(l, Aps[2]))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

      // face--><--line
      for (const auto &l : A->getLines())
         if (!isLinkedDoubly(l, A))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      for (const auto &l : B->getLines())
         if (!isLinkedDoubly(l, B))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

            // std::cout << "Abl->getPoints() = " << Abl->getPoints() << std::endl;
            // std::cout << "Afl->getPoints() = " << Afl->getPoints() << std::endl;
            // std::cout << "Bfl->getPoints() = " << Bfl->getPoints() << std::endl;
            // std::cout << "Bbl->getPoints() = " << Bbl->getPoints() << std::endl;
            // std::cout << "this->getPoints() = " << this->getPoints() << std::endl;
#endif

      Afl->setBoundsSingle();  // lineのセット
      Abl->setBoundsSingle();  // lineのセット
      Bfl->setBoundsSingle();  // lineのセット
      Bbl->setBoundsSingle();  // lineのセット
      //
      A->setGeometricProperties(ToX(A->setPoints(A->Lines)));
      B->setGeometricProperties(ToX(B->setPoints(B->Lines)));
      //
      p0->setBoundsSingle();
      p0->setFaces(p0->Lines);
      for (const auto &f : p0->getFaces())
         f->setGeometricProperties(ToX(f->setPoints(f->Lines)));
      p1->setBoundsSingle();
      p1->setFaces(p1->Lines);
      //
      for (const auto &f : p1->getFaces())
         f->setGeometricProperties(ToX(f->setPoints(f->Lines)));
      return true;
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
//
inline bool networkLine::islegal() const {
   try {
      auto fs = this->getFaces();
      if (fs.size() < 2)
         return true;

      auto [Ap0, Ap1, Ap2] = ToX(fs[0]->getPoints(this));
      auto [Bp0, Bp1, Bp2] = ToX(fs[1]->getPoints(this));
      // double sumangleA = M_PI - std::abs(MyVectorAngle(Aps[2] - Aps[1], Aps[0] - Aps[2]));
      // double sumangleB = M_PI - std::abs(MyVectorAngle(Bps[2] - Bps[1], Bps[0] - Bps[2]));

      double sumangle = std::abs(MyVectorAngle(Ap1 - Ap2, Ap0 - Ap2));
      sumangle += std::abs(MyVectorAngle(Bp1 - Bp2, Bp0 - Bp2));
      if (sumangle > 2. * M_PI || 0. > sumangle)
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      else if (M_PI /*180.0.....1 若干大きい場合はOKとする*/ > sumangle)
         return true /*正*/;
      else
         return false /*不正*/;
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};
inline bool networkLine::flipIfIllegal() {
   if (!islegal() && !isIntxn())
      return (this->flip());  // flipはflipが成功した場合trueを返す．convexでない場合flipされない場合がある
   else
      return false;
};

inline bool networkLine::isFlat(const double minangle = M_PI / 180.) const {
   auto fs = this->getFaces();
   if (fs.size() != 2)
      return false;
   else if (Dot(fs[0]->normal, fs[1]->normal) > cos(minangle))
      return true;
   else
      return false;
};

inline bool networkLine::canflip(const double min_inner_angle = M_PI / 180.) const {
   /*
              f2 *------* f1,F0
                     |this/ |
           F1,f0 *------* F2
   */
   /*
   現在の最小角度よりも大きかったらtrue．
   そうでなくとも，指定されたmin_inner_angleよりも角度が大きくなればそれでもtrueを返す
   */
   auto f0f1 = this->getFaces();
   auto [f0, f1, f2] = f0f1[0]->getPoints(this);
   auto [F0, F1, F2] = f0f1[1]->getPoints(this);
   //
   auto nextTrig0 = T3Tddd{f0->getXtuple(), F2->getXtuple(), f2->getXtuple()};
   auto nextTrig1 = T3Tddd{F0->getXtuple(), f2->getXtuple(), F2->getXtuple()};
   //
   Tddd angles0 = TriangleAngles(nextTrig0);  // flip後の三角形
   Tddd angles1 = TriangleAngles(nextTrig1);  // flip後の三角形
   bool isfiniteangles = (isFinite(angles0) && isFinite(angles1));
   bool isfiniteareas = (isFinite(TriangleArea(nextTrig0)) && isFinite(TriangleArea(nextTrig1)));
   auto n0 = TriangleNormal(nextTrig0);
   auto n1 = TriangleNormal(nextTrig1);
   bool isfinitenormal = (isFinite(n0) && isFinite(n1));
   bool isPositive = (Dot(f0f1[0]->normal, n0) >= 0. && Dot(f0f1[0]->normal, n1) >= 0.) && (Dot(f0f1[1]->normal, n0) >= 0. && Dot(f0f1[1]->normal, n1) >= 0.);
   bool currentMin = Min(Tdd{Min(f0f1[0]->getAngles()), Min(f0f1[1]->getAngles())});
   double min0 = Min(angles0);
   double min1 = Min(angles1);
   bool notTooSmall0 = (min0 > min_inner_angle) || (min0 > currentMin);
   bool notTooSmall1 = (min1 > min_inner_angle) || (min1 > currentMin);
   return (isfiniteangles && isfiniteareas && isfinitenormal && notTooSmall1 && notTooSmall0 && isPositive);
};

inline bool networkLine::flipIfBetter(const double min_degree_to_flat,
                                      const double min_inner_angle) {
   //@ flipを実行するには面の法線方向が成す全ての角度かこれよりも小さくなければならない
   //@ フリップ前後の両方で不正な辺と判定された場合，
   //@ 線の数と面の面積の差をチェックし，差が少ない方を選択する．
   if (!canflip(min_inner_angle))
      return false;
   else if (this->isFlat(min_degree_to_flat /*ここで引っかかってしまいフリップされないことがよくある*/) && !islegal() && !isIntxn()) {
      {
         auto [p0, p1] = this->getPoints();
         auto f0f1 = this->getFaces();
         int s0 = p0->getLines().size();
         int s1 = p1->getLines().size();
         auto p2 = f0f1[0]->getPointOpposite(this);
         auto p3 = f0f1[1]->getPointOpposite(this);
         int s2 = p2->getLines().size();
         int s3 = p3->getLines().size();
         // if (s0 > 3 || s1 > 3 || s2 > 3 || s3 > 3)
         //    if (s0 - 1 < 4 || s1 - 1 < 4 || s2 + 1 < 4 || s3 + 1 < 4)
         //       return false;  // 3以下はつくらない
      }
      /* -------------------------------------------------------------------------- */
      // auto [p0, p2] = this->getPoints();
      auto f0f1 = this->getFaces();
      auto [f0pb, f0pf, f0po] = f0f1[0]->getPoints(this);
      auto [f1pb, f1pf, f1po] = f0f1[1]->getPoints(this);
      // 面積の減少
      // double diffAinit = std::abs(TriangleArea(T3Tddd{f0pb->getXtuple(), f0pf->getXtuple(), f0po->getXtuple()}) - TriangleArea(T3Tddd{f1pb->getXtuple(), f1pf->getXtuple(), f1po->getXtuple()}));
      // double diffAnext = std::abs(TriangleArea(T3Tddd{f0pb->getXtuple(), f1po->getXtuple(), f0po->getXtuple()}) - TriangleArea(T3Tddd{f0pf->getXtuple(), f0po->getXtuple(), f1po->getXtuple()}));
      // double diff = (diffAnext - diffAinit) / (f0f1[0]->getArea() + f0f1[1]->getArea());
      /*
              f0po *------* f0pf,f1pb
                        |   /  |
      f0pb,f1pf *------* f1po
      */

      double min_init = Min(TriangleAngles(T3Tddd{f0pb->getXtuple(), f0pf->getXtuple(), f0po->getXtuple()}));
      auto tmp_angle = Min(TriangleAngles(T3Tddd{f1pb->getXtuple(), f1pf->getXtuple(), f1po->getXtuple()}));
      if (min_init > tmp_angle)
         min_init = tmp_angle;

      double min_later = Min(TriangleAngles(T3Tddd{f0pb->getXtuple(), f1po->getXtuple(), f0po->getXtuple()}));
      tmp_angle = Min(TriangleAngles(T3Tddd{f1pb->getXtuple(), f0po->getXtuple(), f1po->getXtuple()}));
      if (min_later > tmp_angle)
         min_later = tmp_angle;
      //
      // int s0 = f0pb->getLines().size();
      // int s1 = f0pf->getLines().size();
      // int s2 = f0po->getLines().size();
      // int s3 = f1po->getLines().size();
      // double s_mean = 6; //(s0 + s1 + s2 + s3) / 4.;
      // double v_init = std::pow(s0 - s_mean, 2) + std::pow(s1 - s_mean, 2) + std::pow(s2 - s_mean, 2) + std::pow(s3 - s_mean, 2);
      // double v_next = std::pow(s0 - 1 - s_mean, 2) + std::pow(s1 - 1 - s_mean, 2) + std::pow(s2 + 1 - s_mean, 2) + std::pow(s3 + 1 - s_mean, 2);
      // double c = 0.;
      int s0 = f0pb->getLines().size();
      int s1 = f0pf->getLines().size();
      // int s2 = f0po->getLines().size();
      // int s3 = f1po->getLines().size();
      // double s_mean = 6; //(s0 + s1 + s2 + s3) / 4.;
      // double v_init = std::pow(s0 - s_mean, 2) + std::pow(s1 - s_mean, 2) + std::pow(s2 - s_mean, 2) + std::pow(s3 - s_mean, 2);
      // double v_next = std::pow(s0 - 1 - s_mean, 2) + std::pow(s1 - 1 - s_mean, 2) + std::pow(s2 + 1 - s_mean, 2) + std::pow(s3 + 1 - s_mean, 2);

      int next_s0 = s0 - 1;
      int next_s1 = s1 - 1;
      // int next_s2 = s2 + 1;
      // int next_s3 = s3 + 1;
      /*
      @ <-------------(-c)------------ (c) --------------
      @ <---- flip -----|--- topology --|----- none -----
      */
      // if (min_init < min_later)
      // if (std::abs(min_init - min_later) < M_PI / 180. * 5)
      // {
      // 	this->flipIfTopologicalyBetter(min_degree_to_flat);
      // 	return true;
      // }
      // else
      if (min_init <= min_later && !(next_s0 <= 4 && next_s1 <= 4)) {
         this->flip();
         return true;
      } else
         return false;
   } else
      return false;
};

inline bool networkLine::flipIfTopologicalyBetter(const double min_degree_of_line,
                                                  const double min_degree_of_face,
                                                  const int s_meanIN) {
   //@ flipを実行するには面の法線方向が成す全ての角度かこれよりも小さくなければならない
   //@ フリップ前後の両方で不正な辺と判定された場合，
   //@ 線の数と面の面積の差をチェックし，差が少ない方を選択する．
   try {
      if (!canflip(min_degree_of_face) /*flipした後の三角形の最小角度*/)
         return false;
      //@ flipを実行するには面の法線方向が成す全ての角度かこれよりも小さくなければならない
      //@ フリップ前後の両方で不正な辺と判定された場合，
      //@ 線の数と面の面積の差をチェックし，差が少ない方を選択する．
      if (this->isFlat(min_degree_of_line) && !isIntxn()) {
         auto [p0, p1] = this->getPoints();
         auto f0f1 = this->getFaces();
         int s0 = p0->getLines().size();
         int s1 = p1->getLines().size();
         auto p2 = f0f1[0]->getPointOpposite(this);
         auto p3 = f0f1[1]->getPointOpposite(this);
         int s2 = p2->getLines().size();
         int s3 = p3->getLines().size();
         // if (s0 > 3 || s1 > 3 || s2 > 3 || s3 > 3)
         //    if (s0 - 1 < 4 || s1 - 1 < 4 || s2 + 1 < 4 || s3 + 1 < 4)
         //       return false;  // 3以下はつくらない

         double s_mean = s_meanIN;  //(s0 + s1 + s2 + s3) / 4.;
         double v_init = std::pow(s0 - s_mean, 2) + std::pow(s1 - s_mean, 2) + std::pow(s2 - s_mean, 2) + std::pow(s3 - s_mean, 2);
         double v_next = std::pow(s0 - s_mean - 1, 2) + std::pow(s1 - s_mean - 1, 2) + std::pow(s2 - s_mean + 1, 2) + std::pow(s3 - s_mean + 1, 2);
         // flipはflipが成功した場合trueを返す．convexでない場合flipされない場合がある

         // if (p2->CORNER && !p3->CORNER)
         // {
         // 	if (s2 <= 5)
         // 	{
         // 		this->flip();
         // 		return true;
         // 	}
         // }
         // else if (!p2->CORNER && p3->CORNER)
         // {
         // 	if (s3 <= 5)
         // 	{
         // 		this->flip();
         // 		return true;
         // 	}
         // }

         if (v_init > v_next || (s2 <= 4 || s3 <= 4) || (s0 >= 8 || s1 >= 8)) {
            this->flip();
            return true;
         } else
            return false;
      } else
         return false;
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};

// inline bool networkLine::flipIfBetter()
// {
// 	if (!this->isFlat(1.) && !islegal() && !isIntxn())
// 	{
// 		auto fs = this->getFaces();
// 		auto diff_area = std::abs(fs[0]->getArea() - fs[1]->getArea());
// 		this->flip(); // flipはflipが成功した場合trueを返す．convexでない場合flipされない場合がある
// 		// まだ不正な線のままで，点の数が違いすぎる場合戻す．
// 		if (!this->isFlat(1.) && !islegal() && !isIntxn())
// 		{
// 			auto fs = this->getFaces();
// 			if (std::abs(std::abs(fs[0]->getArea() - fs[1]->getArea())) > diff_area)
// 			{
// 				this->flip();
// 				return false;
// 			}
// 		}
// 		return false;
// 	}
// 	else
// 		return false;
// };

inline void networkLine::divideIfIllegal() {
   if (!islegal())
      this->divide();
};

   // inline bool networkLine::isXLine() const
   // {
   //     auto fs = this->Faces;
   //     if (this->Points.size() == 2)
   //     {
   //         if (fs.size() < 2)
   //             return false;
   //         else
   //         {
   //             for (auto i = 0; i < fs.size(); i++)
   //                 for (auto j = i + 1; j < fs.size(); j++)
   //                     if (fs[i]->getNetwork() != fs[j]->getNetwork())
   //                         return true;
   //             return false;
   //         }
   //     }
   //     else
   //         return false;
   //     /*face-face intersection*/
   // };

#endif