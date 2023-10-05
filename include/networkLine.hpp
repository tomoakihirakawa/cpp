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
   auto bool1 = this->replace(oldP, newP);  // 1
   auto bool2 = oldP->Erase(this);          // 2
   auto bool3 = newP->Add(this);            // 3
   // このステップがdouble replace
   //  switchでないと，順番に意味のあるFaceではおかしくなるので注意
   if (bool1 && bool2 & bool3)
      return true;
   else
      return false;
};
// inline bool networkLine::Replace(netF *oldF, netF *newF, netL *newL)
// {
// 	auto bool1 = this->replace(oldF, newF); // 1
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
    : CoordinateBounds(Tddd{{0., 0., 0.}}),
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

   // divideによって境界面が大きく変わる可能性がある．
   // もし境界面が変形すると困る場合は，注意する．

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

   std::unordered_set<networkFace *> related_faces;

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
            for (auto &f : fP->getFaces())
               if (f != nullptr)
                  related_faces.emplace(f);
            for (auto &f : oP->getFaces())
               if (f != nullptr)
                  related_faces.emplace(f);
            for (auto &f : bP->getFaces())
               if (f != nullptr)
                  related_faces.emplace(f);
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

            for (auto &f : fP1->getFaces())
               if (f != nullptr)
                  related_faces.emplace(f);
            for (auto &f : oP1->getFaces())
               if (f != nullptr)
                  related_faces.emplace(f);
            for (auto &f : bP1->getFaces())
               if (f != nullptr)
                  related_faces.emplace(f);
         }
         ///////////////////////////////////////////////////////
         // boundsSetter bSetter;
         // bSetter.add({fP, oP, bP});
         ///////////////////////////////////////////////////////
         netPp newP = nullptr;
         netLp newDivL, newMidL;
         try {
            // std::cout << ColorFunction(c++) << "|" << colorOff;
            newP = new networkPoint(this->getNetwork() /*属性*/, isFinite(midX) ? midX : 0.5 * (fP->X + bP->X));
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
            if (!(this->replace(fP, newP)))
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
            if (!(oldF->replace(fL, newMidL))) {
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
            if (!(fL->replace(oldF, newF))) {
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
            if (!(newF->replace(bL, newMidL)))
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

            if (!(newF->replace(this, newDivL))) {
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
            if (!(oldF1->replace(bL1, newMidL1)))
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
            if (!(bL1->replace(oldF1, newF1)))
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
            if (!(newF1->replace(fL1, newMidL1)))
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
            if (!(newF1->replace(this, newDivL))) {
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
            //
            /* ---------------------------------------------------------- */
            // auto set = [&](networkPoint *newP) {
            //    newP->setFaces();
            //    for (const auto &p : newP->getNeighbors())
            //       p->setFaces();
            //    for (const auto &l : newP->getLines())
            //       l->setBoundsSingle();
            //    for (const auto &f : newP->getFacesFromLines())
            //       f->setGeometricProperties(ToX(f->setPoints()));
            // };
            // set(newP);
            // set(bP);
            // set(fP);
            // set(oP);
            // set(oP1);
            /* ---------------------------------------------------------- */
            // fP->getNetwork()->setGeometricProperties();
            for (auto &f : related_faces)
               f->setGeometricProperties(ToX(f->setPoints()));
            /* ---------------------------------------------------------- */
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
   std::ranges::for_each(Join(Bps, Aps), [&](const auto &p) {		if (p->getLines().empty())
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

inline netPp networkLine::merge() {
   auto net = this->getNetwork();
   // net->setGeometricProperties();
   std::cout << "choose a deleting point, lines and faces" << std::endl;
   const auto Fs = this->getFaces();
   const auto f_del0 = Fs[0];
   const auto f_del1 = Fs[1];
   const auto p_del = f_del0->getPointBack(this);
   const auto p_rem = f_del0->getPointFront(this);
   //
   const auto l_del0 = f_del0->getLineBack(this);
   const auto l_remain0 = f_del0->getLineFront(this);
   const auto l_del1 = f_del1->getLineFront(this);
   const auto l_remain1 = f_del1->getLineBack(this);
   //
   const auto faces = p_del->getFaces();
   const auto lines = p_del->getLines();
   const auto points = Join(p_del->getNeighbors(), p_rem->getNeighbors());
   // is merge able ?
   for (const auto &f : faces) {
      if (f != f_del0 && f != f_del1) {
         if (f->MemberQ(l_remain0) || f->MemberQ(l_remain1) || f->MemberQ(this) || f->MemberQ(p_rem))
            return nullptr;  // not mergeable
         if (f->MemberQ(l_del0) && f->MemberQ(l_del1))
            return nullptr;  // not mergeable
         if (l_del0 == l_del1 || p_del == p_rem || l_remain0 == l_remain1)
            return nullptr;  // not mergeable
      }
   }
   //
   // std::cout << Red << "delete f_del0 " << colorOff << f_del0 << std::endl;
   // std::cout << Red << "delete f_del1 " << colorOff << f_del1 << std::endl;
   // std::cout << Red << "delete l_del0 " << colorOff << l_del0 << std::endl;
   // std::cout << Red << "delete l_del1 " << colorOff << l_del1 << std::endl;
   // std::cout << Red << "delete p_del " << colorOff << p_del << std::endl;
   // std::cout << Red << "delete this " << colorOff << this << std::endl;
   //
   // std::cout << "set" << std::endl;
   p_rem->setX((p_rem->X + p_del->X) / 2.);
   //
   // std::cout << "replace" << std::endl;
   for (const auto &f : faces) {
      if (f != f_del0 && f != f_del1) {
         if (f->replace(l_del0, l_remain0)) {
            l_remain0->replace(f_del0, f);
            p_rem->replace(f_del0, f);
         } else if (f->replace(l_del1, l_remain1)) {
            l_remain1->replace(f_del1, f);
            p_rem->replace(f_del1, f);
         }
         f->replace(p_del, p_rem);
      }
   }
   // std::cout << "replace and add" << std::endl;
   p_rem->Erase(this);
   for (const auto &l : lines) {
      if (l != l_del0 && l != l_del1 && l != this) {
         l->replace(p_del, p_rem);
         p_rem->Add(l);
      }
      p_del->Erase(l);
   }
   for (const auto &p : points) {
      p->Erase(this);
      p->Erase(l_del0);
      p->Erase(l_del1);
   }

   std::cout << "delete f_del0 " << f_del0 << std::endl;
   delete f_del0;
   std::cout << "delete f_del1 " << f_del1 << std::endl;
   delete f_del1;
   std::cout << "delete l_del0 " << l_del0 << std::endl;
   delete l_del0;
   std::cout << "delete l_del1 " << l_del1 << std::endl;
   delete l_del1;
   std::cout << "delete p_del " << p_del << std::endl;
   delete p_del;
   std::cout << "delete this " << this << std::endl;
   delete this;
   //
   /* ----------------------------------------- */
   // p_rem->setFaces();
   // for (const auto &p : p_rem->getNeighbors())
   //    p->setFaces();
   // for (const auto &l : p_rem->getLines())
   //    l->setBoundsSingle();
   // for (const auto &f : p_rem->getFacesFromLines())
   //    f->setGeometricProperties(ToX(f->setPoints()));
   // /* ----------------------------------------- */
   net->setGeometricProperties();
   /* ----------------------------------------- */
   // std::cout << "done" << std::endl;
   return p_rem;
};

/* -------------------------------------------------------------------------- */

// #define debug_merge
// inline netPp networkLine::merge() {
//    auto AB = this->getFaces();
//    auto A = AB[0];
//    auto B = AB[1];
//    auto Aps = A->getPoints(this);
//    auto A0 = std::get<0>(Aps);
//    auto A1 = std::get<1>(Aps);
//    auto Bps = B->getPoints(this);
//    //
//    auto Abl = A->getLineBack(this);
//    auto Afl = A->getLineFront(this);
//    auto Bbl = B->getLineBack(this);
//    auto Bfl = B->getLineFront(this);

//    auto a = (*Abl)(A);
//    auto b = (*Bfl)(B);

//    auto aa = (*Afl)(A);
//    auto bb = (*Bbl)(B);

//    // aaとbbが同じ！？divideのミスでしょう．

//    std::cout << "A " << A << std::endl;
//    std::cout << "B " << B << std::endl;
//    // std::cout << "a points " << a->getPoints(Aps[0]) << std::endl;
//    // std::cout << "b points " << b->getPoints(Aps[0]) << std::endl;
//    std::cout << "a " << a << std::endl;
//    std::cout << "b " << b << std::endl;
//    std::cout << "aa " << aa << std::endl;
//    std::cout << "bb " << bb << std::endl;
//    std::cout << "(*Afl)(A) " << (*Afl)(A) << std::endl;
//    std::cout << "(*Bbl)(B) " << (*Bbl)(B) << std::endl;

//    auto Aps0lines = TakeExcept(std::get<0>(Aps)->getLines(), {this, Abl, Bfl});
//    try {
// #if defined(debug_merge)
//       std::cout << green << "|" << colorOff;
// #endif

//       /*                          Aps[2]
//        *                          V    V
//        *                         A      A
//        *                        /         \
//        *      (will be deleted)/    / \    \
//        *                   Abl/-><-/   \-><-\Afl
//        *                     /    /  A  \ 　  \　　
//        *                    V     ---V---      V
//        *  (will be deleted)A         A          A
//        *        Bps[1],Aps[0]---><--this---><---Aps[1],Bps[0]
//        *                  V          |          V
//        *                   A         V         A
//        *   (will be deleted)\     ---A---     /
//        *                  Bfl\-><-\     /-><-/Bbl
//        *                      \    \ B /    /
//        *                       \    \ /    /
//        *                        V        V
//        *                         A      A
//        *                          Bps[2]
//        */

// #if defined(debug_merge)
//       std::cout << green << "|" << colorOff;
// #endif
//       auto mid_X = (std::get<0>(Aps)->getXtuple() + std::get<1>(Aps)->getXtuple()) / 2.;

// #if defined(debug_merge)
//       std::cout << green << "|" << colorOff;
// #endif

//       /*                      Aps[2]
//        *                     V    V
//        *                    A      A
//        *                   /         \          ______
//        *       ___        /    / \    \         \    /
//        *       \a/-><-Abl/-><-/   \-><-\Afl-><---\  /
//        *        V       /    /  A  \ 　  \　　     \/
//        *               V     ---V---      V
//        *              A         A          A
//        *  Bps[1],Aps[0]---><--this---><---Aps[1],Bps[0]
//        *             V          |          V
//        *              A         V         A
//        *               \     ---A---     /
//        *        b-><-Bfl\-><-\     /-><-/Bbl
//        *                 \    \ B /    /
//        *                  \    \ /     /
//        *                    V        V
//        *                     A      A
//        *                      Bps[2]
//        */
// #if defined(debug_merge)
//       std::cout << green << "|" << colorOff;
// #endif
//       if (!(this->Erase(A)))
//          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//       if (!(a->replace(Abl, Afl)))
//          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//       if (!(a->replace(A0, A1)))
//          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//       if (!(Afl->replace(A, a)))
//          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//       if (!(Abl->Erase(a)))
//          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//          /*                     Aps[2]
//           *                    V       V
//           *        _____      A         A
//           *         \a/ ----><-----------\Afl       ______
//           *          V      /     / \     \         \    /
//           *               Abl    / A \     \----><---\  /
//           *        \      /     /     \  　  \　　     \/
//           *         \    V      ---V---       V
//           *          V  A          A           A
//           *-->Bps[1],Aps[0]---><--this---><---Aps[1],Bps[0]
//           *         A   V          V          V
//           *        /     A         A         A
//           *               \     -------     /
//           *(*Bfl)(B)-><-Bfl\-><-\     /-><-/Bbl
//           *                 \    \ B /    /
//           *                  \    \ /    /
//           *                   V         V
//           *                    A       A
//           *                      Bps[2]
//           */

// #if defined(debug_merge)
//       std::cout << green << "|" << colorOff;
// #endif
//       if (!(this->Erase(B)))
//          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//       if (!(b->replace(Bfl, Bbl)))
//          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//       if (!(b->replace(A0, A1)))
//          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//       if (!(Bbl->replace(B, b)))
//          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//       if (!(Bfl->Erase(b)))
//          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

//       /*                   Aps[2]
//        *                   V       V
//        *       ____       A         A
//        *       \a/ ----></---------->\Afl       ______
//        *        V       /     / \     \         \    /
//        *              Abl    /   \     \----><---\  /
//        *       \      /     /  A  \  　  \　　     \/
//        *        \    V      -------       V
//        *         V  A                      A
//        *-->Bps[1],Aps[0]-><--this---><----Aps[1],Bps[0]
//        *        A   V                      V
//        *       /    A                      A
//        *      /       \      -------     /
//        *              Bfl    \     /    /Bbl
//        *                 \    \ B /    /
//        *         /\       \    \ /    /
//        *        /b_\----><-\-------->/Bbl
//        *                    A       A
//        *                      Bps[2]
//        */

//       // for (const auto &l : Aps[0]->getLines())
//       // auto tmp = TakeExcept(Aps[0]->getLines(), {this, Abl, Bfl});
//       // auto tmp = std::get<0>(Aps)->getLines();
//       auto tmp = TakeExcept(std::get<0>(Aps)->getLines(), {this, Abl, Bfl});
//       for (auto &l : tmp) {
//          {
//             if (!(l->replace(std::get<0>(Aps), std::get<1>(Aps))))
//                throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//             if (!(std::get<0>(Aps)->Erase(l)))
//                throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//             std::get<1>(Aps)->Add(l);
//          }
//          // Aps[1]が3lineしか持っていない場合もあり得るので，Addできない場合もある
//          // throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//       }

//       for (auto &l : std::get<0>(Aps)->getLines()) {
//          std::get<0>(Aps)->Erase(l);
//       }
//       this->Erase(std::get<0>(Aps));
//       Abl->Erase(std::get<0>(Aps));
//       Bfl->Erase(std::get<0>(Aps));
//       Abl->Erase(A);
//       Bfl->Erase(B);
//       this->Erase(A);
//       this->Erase(B);
//       /*                   Aps[2]
//        *                   V       V
//        *       ____       A         A
//        *       \a/ ----></---------->\Afl       ______
//        *        V       /     / \     \         \    /
//        *    (points)  Abl    /   \     \----><---\  /
//        *     (lines)        /  A  \  　  \　　     \/     /
//        *                    -------       V             /
//        *                                   A           V
//        *   Bps[1],Aps[0]     this---><----Aps[1],Bps[0]<----
//        *                                   V        A
//        *                                   A         \
//        *                     -------     /             \
//        *     (lines)  Bfl    \     /    /Bbl
//        *    (points)     \    \ B /    /
//        *       /b \       \    \ /    /
//        *      /____\----><-\-------->/Bbl
//        *                    A       A
//        *                      Bps[2]
//        */

// #if defined(debug_merge)
//       if (!Abl->getFaces().empty())
//          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//       if (!Bfl->getFaces().empty())
//          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

// #endif

//       for (const auto &f : A0->getFaces()) {
//          f->replace(A0, A1);
//          for_each(f->getLines(), [&](const auto &l) {
//             l->replace(A0, A1);
//             l->replace(A0, A1);
//          });
//       }

//       std::get<2>(Aps)->erase(A);
//       std::get<0>(Aps)->erase(A);
//       std::get<0>(Aps)->erase(B);
//       std::get<0>(Aps)->erase(a);
//       std::get<0>(Aps)->erase(b);
//       std::get<2>(Bps)->erase(B);
//       //
//       std::get<1>(Aps)->erase(A);
//       std::get<1>(Aps)->erase(B);

//       // std::cout << " deleting A = " << A << std::endl;
//       delete A;
//       // std::cout << " deleting B = " << B << std::endl;
//       delete B;
//       // std::cout << " deleting Abl = " << Abl << std::endl;
//       delete Abl;
//       // std::cout << " deleting Bfl = " << Bfl << std::endl;
//       delete Bfl;
//       // std::cout << " deleting this = " << this << std::endl;
//       delete this;
//       delete std::get<0>(Aps);  // thisもきえます

//       /*                   Aps[2]
//        *                           V
//        *     ______                 A
//        *     \ a / ----><----------->\Afl       ______
//        *      \ /                     \         \    /
//        *                               \----><---\  /
//        *                             　  \　　     \/     /
//        *                                  V            /
//        *                                   A          V
//        *                                  Aps[1],Bps[0]<----
//        *                                   V          A
//        *                                   A           \
//        *                                  /             \
//        *                                 /Bbl
//        *                                /
//        *       /b \                    /
//        *      /____\----><---------->/Bbl
//        *                            A
//        *                      Bps[2]
//        */

// #if defined(debug_merge)
//       std::cout << " l = " << std::get<0>(Aps)->getLines() << std::endl;
//       std::cout << " l = " << std::get<1>(Aps)->getLines() << std::endl;
//       //
//       std::cout << " std::get<1>(Aps)->getFaces() = " << std::get<1>(Aps)->getFaces() << std::endl;
//       std::cout << " std::get<1>(Aps)->getNeighbors() = " << std::get<1>(Aps)->getNeighbors() << std::endl;
//       //
//       std::cout << " Afl->points = " << Afl->getPoints() << std::endl;
//       std::cout << " Bbl->points = " << Bbl->getPoints() << std::endl;
//       //
//       std::cout << " Afl->faces = " << Afl->getFaces() << std::endl;
//       std::cout << " Bbl->faces = " << Bbl->getFaces() << std::endl;
//       //
//       if (ConnectedQ(std::get<1>(Aps), std::get<2>(Aps)))
//          std::cout << "Connected Q = true" << std::endl;
//       if (ConnectedQ(std::get<1>(Aps), std::get<2>(Bps)))
//          std::cout << "Connected Q = true" << std::endl;

//       for (const auto &l : V_netLp{Afl, Bbl}) {
//          if (l) {
//             for (const auto &f : l->getFaces()) {
//                if (f) {
//                   if (!isLinkedDoubly(f, l))
//                      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//                } else {
//                   std::cout << "f = nullptr" << std::endl;
//                   throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//                }
//             }
//          } else {
//             std::cout << "l = nullptr" << std::endl;
//             throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//          }
//       }
//       std::cout << "Clear " << std::endl;
// #endif
//       std::get<1>(Aps)->setX(mid_X);
//       std::cout << "mid_X = " << mid_X << std::endl;
//       std::cout << "std::get<1>(Aps)->X = " << std::get<1>(Aps)->X << std::endl;
//       return std::get<1>(Aps);
//    } catch (std::exception &e) {
//       std::cerr << e.what() << colorOff << std::endl;
//       std::stringstream ss;
//       ss << "Aps " << Aps << std::endl;
//       ss << "Bps " << Bps << std::endl;
//       ss << "Abl " << Abl << std::endl;
//       ss << "Afl " << Afl << std::endl;
//       ss << "Bbl " << Bbl << std::endl;
//       ss << "Bfl " << Bfl << std::endl;
//       ss << "a " << a << std::endl;
//       ss << "b " << b << std::endl;
//       ss << "aa " << aa << std::endl;
//       ss << "bb " << bb << std::endl;
//       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
//    };
// };

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

      if (!(Abl->replace(A, B)))
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

      if (!(Bbl->replace(B, A)))
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
      double sumangle = std::abs(VectorAngle(Bp1 - Bp2, Bp0 - Bp2)) + std::abs(VectorAngle(Ap1 - Ap2, Ap0 - Ap2));
      if (M_PI /*180.0.....1 若干大きい場合はOKとする*/ > sumangle)
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

inline bool networkLine::isAdjacentFacesFlat(const double minangle = M_PI / 180.) const {
   auto fs = this->getFaces();
   auto [p0, p1, p2] = fs[0]->getPoints();
   auto [P0, P1, P2] = fs[1]->getPoints();
   if (fs.size() != 2)
      return false;
   else if (isFlat(Cross(p1->X - p0->X, p2->X - p0->X), Cross(P1->X - P0->X, P2->X - P0->X), minangle))
      return true;
   else
      return false;
};

// Define actual types for T3Tddd, Tddd, and others if not defined

// const double DEFAULT_MIN_INNER_ANGLE = M_PI / 180.0;

// // Helper function to check if a triangle should flip based on angles
// bool shouldFlipBasedOnAngle(const double min_angle, const double currentMin, const Tddd &angles) {
//    double min_angle_calculated = Min(angles);
//    return (min_angle_calculated > min_angle) || (min_angle_calculated > currentMin);
// }

// inline bool networkLine::canFlip(const double min_inner_angle = DEFAULT_MIN_INNER_ANGLE) const {
//    try {
//       // Fetch relevant points and faces
//       auto faces = this->getFaces();
//       auto [point_f0, point_f1, point_f2] = faces[0]->getPoints(this);
//       auto [point_F0, point_F1, point_F2] = faces[1]->getPoints(this);

//       // Compute new triangles if a flip happens
//       auto nextTriangle0 = T3Tddd{point_f0->getXtuple(), point_F2->getXtuple(), point_f2->getXtuple()};
//       auto nextTriangle1 = T3Tddd{point_F0->getXtuple(), point_f2->getXtuple(), point_F2->getXtuple()};

//       if (TriangleArea(nextTriangle0) == 0.0 || TriangleArea(nextTriangle1) == 0.0)
//          return false;

//       // Compute angles and normals for new triangles
//       Tddd angles0 = TriangleAngles(nextTriangle0);
//       Tddd angles1 = TriangleAngles(nextTriangle1);
//       auto normal0 = TriangleNormal(nextTriangle0);
//       auto normal1 = TriangleNormal(nextTriangle1);

//       // Condition checks
//       bool isFiniteAngles = isFinite(angles0) && isFinite(angles1);
//       bool isFiniteAreas = isFinite(TriangleArea(nextTriangle0)) && isFinite(TriangleArea(nextTriangle1));
//       bool isFiniteNormals = isFinite(normal0) && isFinite(normal1);
//       bool isNormalsPositive = (Dot(faces[0]->normal, normal0) >= 0.0 && Dot(faces[0]->normal, normal1) >= 0.0) &&
//                                (Dot(faces[1]->normal, normal0) >= 0.0 && Dot(faces[1]->normal, normal1) >= 0.0);
//       double currentMinimumAngle = Min(Tdd{Min(faces[0]->getAngles()), Min(faces[1]->getAngles())});

//       return isFiniteAngles && isFiniteAreas && isFiniteNormals &&
//              shouldFlipBasedOnAngle(min_inner_angle, currentMinimumAngle, angles0) &&
//              shouldFlipBasedOnAngle(min_inner_angle, currentMinimumAngle, angles1) &&
//              isNormalsPositive;
//    } catch (std::exception &e) {
//       std::cerr << e.what() << std::endl;  // Removed colorOff for clarity
//       throw;                               // Re-throw the caught exception
//    }
// }

/*DOC_EXTRACT flip

### `flip`可能かどうかの判定

\ref{canFlip}{`canFlip`}でフリップ可能かどうかを判定する．直感的に次のような条件の場合，境界面が崩れるため，フリップさせたくない．

* フリップ前後で，辺に隣接する面の面積の和が大きく変化する場合，フリップさせない
* フリップ前後で，辺に隣接する面の法線ベクトルが大きく変換する場合，フリップさせない

しかし，これの判定において必要となる計算：三角形の内角や法線方向，ベクトルの成す角度の計算は，精確に判定できない領域があるようだ．
なので，その領域をおおよそ実験的に調べて，まずはその領域に入らせない条件を設ける（信頼できる三角形）．
次のような三角形は信頼しない：

* 三角形の内角が小さすぎる，または大きすぎる場合
* 内角の和が$`\pi`$にならない場合

信頼できる三角形の判定には，\ref{isValidTriangle}{`isValidTriangle`}を用いる．

*/

// \label{canFlip}
inline bool networkLine::canFlip(const double acceptable_n_diff_before_after = M_PI / 180.) const {
   try {
      auto f_and_F = this->getFaces();
      auto [f0, f1, f2] = f_and_F[0]->getPoints(this);
      auto [F0, F1, F2] = f_and_F[1]->getPoints(this);
      auto tri0_now = T3Tddd{f0->X, f1->X, f2->X};
      auto tri1_now = T3Tddd{F0->X, F1->X, F2->X};

      auto tri0 = T3Tddd{f0->X, F2->X, f2->X};
      auto tri1 = T3Tddd{F0->X, f2->X, F2->X};

      if (!isValidTriangle(tri0, 5222 * M_PI / 180.))
         return false;
      if (!isValidTriangle(tri1, 5 * M_PI / 180.))
         return false;

      //$ large difference of normal vector after and before flip
      if (!isFlat(Cross(tri0[1] - tri0[0], tri0[2] - tri0[0]), Cross(tri0_now[1] - tri0_now[0], tri0_now[2] - tri0_now[0]), acceptable_n_diff_before_after) ||
          !isFlat(Cross(tri0[1] - tri0[0], tri0[2] - tri0[0]), Cross(tri1_now[1] - tri1_now[0], tri1_now[2] - tri1_now[0]), acceptable_n_diff_before_after) ||
          !isFlat(Cross(tri1[1] - tri1[0], tri1[2] - tri1[0]), Cross(tri1_now[1] - tri1_now[0], tri1_now[2] - tri1_now[0]), acceptable_n_diff_before_after) ||
          !isFlat(Cross(tri1[1] - tri1[0], tri1[2] - tri1[0]), Cross(tri0_now[1] - tri0_now[0], tri0_now[2] - tri0_now[0]), acceptable_n_diff_before_after))
         return false;

      //$ area conservation
      return true;  // TriangleArea(tri0) + TriangleArea(tri1) == TriangleArea(tri0_now) + TriangleArea(tri1_now);

   } catch (const std::exception &e) {
      std::cerr << e.what() << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
}

inline bool networkLine::flipIfBetter(const double n_diff_tagert_face,
                                      const double acceptable_n_diff_before_after,
                                      const int min_n) {
   try {
      //@ flipを実行するには面の法線方向が成す全ての角度かこれよりも小さくなければならない
      //@ フリップ前後の両方で不正な辺と判定された場合，
      //@ 線の数と面の面積の差をチェックし，差が少ない方を選択する．
      if (!canFlip(acceptable_n_diff_before_after))
         return false;
      // else if (this->isAdjacentFacesFlat(n_diff_tagert_face /*ここで引っかかってしまいフリップされないことがよくある*/) && !islegal() && !isIntxn()) {
      else if (this->isAdjacentFacesFlat(n_diff_tagert_face /*ここで引っかかってしまいフリップされないことがよくある*/) && !islegal()) {
         auto [p0, p1] = this->getPoints();
         auto f0f1_ = this->getFaces();
         int s0 = p0->getLines().size();
         int s1 = p1->getLines().size();
         auto p2 = f0f1_[0]->getPointOpposite(this);
         auto p3 = f0f1_[1]->getPointOpposite(this);
         int s2 = p2->getLines().size();
         int s3 = p3->getLines().size();
         // if (s0 > 3 || s1 > 3 || s2 > 3 || s3 > 3)
         //    if (s0 - 1 < 4 || s1 - 1 < 4 || s2 + 1 < 4 || s3 + 1 < 4)
         //       return false;  // 3以下はつくらない

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
         //     f0po *------* f0pf,f1pb
         //          |   /  |
         //f0pb,f1pf *------* f1po
         */

         double min_init = std::min(Min(TriangleAngles(T3Tddd{f0pb->X, f0pf->X, f0po->X})), Min(TriangleAngles(T3Tddd{f1pb->X, f1pf->X, f1po->X})));
         double min_later = std::min(Min(TriangleAngles(T3Tddd{f0pb->X, f1po->X, f0po->X})), Min(TriangleAngles(T3Tddd{f1pb->X, f0po->X, f1po->X})));

         //
         // int s0 = f0pb->getLines().size();
         // int s1 = f0pf->getLines().size();
         // int s2 = f0po->getLines().size();
         // int s3 = f1po->getLines().size();
         // double s_mean = 6; //(s0 + s1 + s2 + s3) / 4.;
         // double v_init = std::pow(s0 - s_mean, 2) + std::pow(s1 - s_mean, 2) + std::pow(s2 - s_mean, 2) + std::pow(s3 - s_mean, 2);
         // double v_next = std::pow(s0 - 1 - s_mean, 2) + std::pow(s1 - 1 - s_mean, 2) + std::pow(s2 + 1 - s_mean, 2) + std::pow(s3 + 1 - s_mean, 2);
         // double c = 0.;
         // int s0 = f0pb->getLines().size();
         // int s1 = f0pf->getLines().size();
         // int s2 = f0po->getLines().size();
         // int s3 = f1po->getLines().size();
         // double s_mean = 6; //(s0 + s1 + s2 + s3) / 4.;
         // double v_init = std::pow(s0 - s_mean, 2) + std::pow(s1 - s_mean, 2) + std::pow(s2 - s_mean, 2) + std::pow(s3 - s_mean, 2);
         // double v_next = std::pow(s0 - 1 - s_mean, 2) + std::pow(s1 - 1 - s_mean, 2) + std::pow(s2 + 1 - s_mean, 2) + std::pow(s3 + 1 - s_mean, 2);

         int next_s0 = s0 - 1;
         int next_s1 = s1 - 1;
         int next_s2 = s2 + 1;
         int next_s3 = s3 + 1;
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
         if (min_init <= min_later
             // &&
             //     (next_s0 >= min_n || p0->CORNER) &&
             //     (next_s1 >= min_n || p1->CORNER) &&
             //     (next_s2 >= min_n || p2->CORNER) &&
             //     (next_s3 >= min_n || p3->CORNER)
         ) {
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

inline bool networkLine::flipIfTopologicallyBetter(const double n_diff_tagert_face,
                                                   const double acceptable_n_diff_before_after,
                                                   const int s_meanIN) {
   try {
      // Check if the flip is allowed based on the minimum angle of the triangle after flipping
      if (!canFlip(acceptable_n_diff_before_after))
         return false;
      if (this->isAdjacentFacesFlat(n_diff_tagert_face) && !isIntxn()) {
         auto [p0, p1] = this->getPoints();
         auto f0f1 = this->getFaces();
         int s0 = p0->getLines().size();
         int s1 = p1->getLines().size();
         auto p2 = f0f1[0]->getPointOpposite(this);
         auto p3 = f0f1[1]->getPointOpposite(this);
         int s2 = p2->getLines().size();
         int s3 = p3->getLines().size();
         double s_mean = s_meanIN;
         double v_init = Norm(std::array<int, 4>{s0, s1, s2, s3} - s_mean);
         double v_next = Norm(std::array<int, 4>{s0 - 1, s1 - 1, s2 + 1, s3 + 1} - s_mean);
         // ただし，例えばs0が角であった場合，角が多くの線を持つことは問題ないため，考慮に入れない．つまり辺の数は変わっても変わらないものとして，v_nextを考える．
         // v_next = Norm(std::array<int, 4>{p0->CORNER || p0->isMultipleNode ? s0 : s0 - 1, p1->CORNER || p1->isMultipleNode ? s1 : s1 - 1,
         //                                  p2->CORNER || p2->isMultipleNode ? s2 : s2 + 1, p3->CORNER || p3->isMultipleNode ? s3 : s3 + 1} -
         //               s_mean);
         if (v_init > v_next || s0 >= 8 || s1 >= 8 || s2 <= 4 || s3 <= 4) {
            this->flip();
            return true;
         } else
            return false;
      } else
         return false;
   } catch (std::exception &e) {
      std::cerr << e.what() << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   }
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