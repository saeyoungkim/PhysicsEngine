/*
made by KIM
*/


#include <Physics/PHContactEngine.h>
#include <Physics/PHConstraintEngine.h>
#include <algorithm>
#include <math.h>

using namespace PTM;
using namespace std;

namespace Spr {
	;
	
	const double EPSILON_2 = 1.0e-2;
	const double EPSILON_4 = 1.0e-4;
	const double EPSILON_8 = 1.0e-8;
	const double EPSILON_10 = 1.0e-10;
	const double EPSILON_16 = 1.0e-16;
	const double EPSILON_20 = 1.0e-20;

	// 内積を計算するための関数
	inline double Dot6(const double* v1, const double* v2) {
		return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3] + v1[4] * v2[4] + v1[5] * v2[5];
	}

	// tmp_1,tmp_2から構成される線分にold_ZMPから射影したときに一番近い点を見つけ出し、その点を返す
	Vec3d ClosestPtPointLine(const Vec3d& old_ZMP, const Vec3d& tmp_1, const Vec3d& tmp_2, double& t)
	{
		// 線分のベクトル
		Vec3d tmp_12 = tmp_2 - tmp_1;
		t = (old_ZMP - tmp_1).dot(tmp_12);

		// 線分の外にある(tmp_12の線分から垂直でtmp_1を通る直線を引いたときにその線より下にある)からこれは 
		if (t <= 0.0) {
			t = 0.0;
			return tmp_1;
		}
		else {
			double nom = tmp_12.dot(tmp_12);
			if (t >= nom) {
				t = 1.0;
				return tmp_2;
			}
			else {
				t = t / nom;
				return tmp_1 + t * tmp_12;
			}
		}
	}


	// Collision Detectionを参考にしました
	// p141～p142
	// 点三つ(tmp_1, tmp_2, tmp_3)に対してある点(old_ZMP)から一番近い点を見つけ出し、その点を返す
	Vec3d ClosestPtPointTriangle(const Vec3d& old_ZMP, const Vec3d& tmp_1, const Vec3d& tmp_2, const Vec3d& tmp_3) {

		Vec3d tmp2_tmp1 = tmp_2 - tmp_1; // ab
		Vec3d tmp3_tmp1 = tmp_3 - tmp_1; // ac
		Vec3d zmp_tmp1 = old_ZMP - tmp_1; // ap

		float d1 = tmp2_tmp1 * zmp_tmp1;
		float d2 = tmp3_tmp1 * zmp_tmp1;
		// Voronoi領域により一番近い点となる
		if (d1 <= 0.0f && d2 <= 0.0f) return tmp_1;

		Vec3d tmp2_zmp = old_ZMP - tmp_2; // bp
		float d3 = tmp2_tmp1 * tmp2_zmp;
		float d4 = tmp3_tmp1 * tmp2_zmp;
		if (d3 >= 0.0f && d4 <= d3) return tmp_2;

		float vc = d1*d4 - d3*d2;
		if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
			float v = d1 / (d1 - d3);
			return tmp_1 + v * tmp2_tmp1;
		}

		Vec3d tmp3_zmp = old_ZMP - tmp_3;
		float d5 = tmp2_tmp1 * tmp3_zmp;
		float d6 = tmp3_tmp1 * tmp3_zmp;
		if (d6 >= 0.0f && d5 <= d6) return tmp_3;

		float vb = d5*d2 - d1*d6;
		if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
			float w = d2 / (d2 - d6);
			return tmp_1 + w * tmp3_tmp1;
		}

		float va = d3*d6 - d5*d4;
		if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
			float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
			return tmp_2 + w * (tmp_3 - tmp_2);
		}

		float denom = 1.0f / (va + vb + vc);
		float v = vb * denom;
		float w = vc * denom;
		return tmp_1 + tmp2_tmp1 * v + tmp3_tmp1 * w;

	}

	// 二点(tmp_1, tmp_2)に対してある点(old_ZMP)から一番近い点を見つけ出し、その点からold_ZMPまでのベクトルを返す
	Vec3d NormalPt2ClosestPt(const Vec3d & old_ZMP, const Vec3d & tmp_1, const Vec3d & tmp_2)
	{
		Vec3d tmp_12 = tmp_2 - tmp_1;
		double t = (old_ZMP - tmp_1).dot(tmp_12);

		if (t <= 0.0) {
			t = 0.0;
			return old_ZMP - tmp_1;
		}
		else {
			double nom = tmp_12.dot(tmp_12);
			if (t >= nom) {
				t = 1.0;
				return old_ZMP - tmp_2;
			}
			else {
				t = t / nom;
				return old_ZMP - (tmp_1 + t * tmp_12);
			}
		}
	}

#define YZ()	sub_vector( PTM::TSubVectorDim<1,2>() )

	// 2DConvex hullを作るために、3DのVerticesを入れる
	// 入れたVerticesはx方向に対して切って、2DのConvex hullを作る
	// ベクトル(normal)を基準に、2DConvex hullからサポートポイントに一番近い点を取り出し、その座標を返す。
	ll PointSupport(const std::vector<Vec3d>& vec_group, const Vec3d& normal) {
		ll n = vec_group.size();
		double d1 = 0.0;
		ll curPos = 0;
		Vec2d normal2D = normal.YZ();
		double maxValue = normal2D * vec_group[0].YZ();

		for (ll i = 1; i < n; i++) {
			d1 = normal.YZ() * vec_group[i].YZ();
			if (d1 > maxValue) {
				maxValue = d1;
				curPos = i;
			}
		}

		return curPos;
	}

	// 3D上の二点(tmp_1, tmp_2)に対してposeを掛けて(各々の物体から位置を計算して)2D上座標に変換する（YZ平面を作る）
	// へっかんさせてに転換の距離を返す(YZ平面上での距離)
	inline double Norm2D(const Posed& pose_1, const Vec3d& tmp_1, const Posed& pose_2, const Vec3d& tmp_2)
	{
		return ((pose_1 * tmp_1).YZ() - (pose_2 * tmp_2).YZ()).norm();
	}

	inline double Norm3D(const Vec3d& pt1, const Vec3d& pt2) {
		return (pt1 - pt2).norm();
	}

	// 3Dの座標群に対して、x方向に対して切って2DのConvex Hullを作る。
	// 2DConvex Hullの中にある点(pt)があるかどうか内包判定する
	bool InOut2D(const std::vector<Vec3d>& vec_group, const Vec3d& pt) {

		size_t size = vec_group.size();

		for (size_t i = 0; i < size; i++) {
			Vec3d t0 = vec_group[i % size] - pt;
			Vec3d t1 = vec_group[(i + 1) % size] - pt;
			Vec3d t2 = vec_group[(i + 2) % size] - pt;
			double d_tmp0 = t0.Y() * t1.Z() - t0.Z() * t1.Y();
			double d_tmp1 = t1.Y() * t2.Z() - t1.Z() * t2.Y();

			if (d_tmp0 * d_tmp1 < 0 && std::abs(d_tmp0 - d_tmp1) >= 0.5e-7) {
				return false;
			}
		}
		// 中にある
		return true;
	}

	// 点群（3DConvexHullから生成された点群であるvec_group）からある点（old_ZMP）に一番近い点（pp）を見つけ出し、その点を返す関数
	Vec3d ClosestPtFromGroup(const Vec3d& old_ZMP, const std::vector<Vec3d>& vec_group) {
		size_t size = vec_group.size();
		double d1 = 0.0, d2 = INFINITY;
		double t = 0.0;
		Vec3d pp;
		for (ui i = 0; i < size; i++) {
			Vec3d tmp = ClosestPtPointLine(old_ZMP, vec_group[i%size], vec_group[(i + 1) % size], t);
			d1 = (old_ZMP - tmp).norm();
			if (d2 > d1) {
				d2 = d1;
				pp = tmp;
			}
		}
		return pp;
	}

	// -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  ----- 
	// PHContactEngine


	// ContactPointをアップデートしていなかったため


	PHContactEngine::PHContactEngine(const Matrix3d& local, PHShapePairForLCP* sp, Vec3d p, PHSolid* s0, PHSolid* s1) :
		PHContactPoint(local, sp, p, s0, s1), 
		ContLocal(local), PHS0(s0), PHS1(s1), 
		CDC0(DCAST(CDConvex, PHS0->GetShape(0))),
		CDC1(DCAST(CDConvex, PHS1->GetShape(0))){

		// 軸設定
		Vec3d u = local.Ex();	//	u: 物体ではなく原点の速度の向きなので - がつく．
		if (u.X() < -1 + EPSILON_16) {
			w2x = Quaterniond::Rot(Rad(180), 'z');
		}
		else if (u.X() < 1 - EPSILON_16) {
			Matrix3d matW2x = Matrix3d::Rot(u, Vec3f(1, 0, 0), 'x');
			w2x.FromMatrix(matW2x);
			w2x = w2x.Inv();
		}

		
		Pose0.Ori() = w2x * PHS0->GetOrientation();
		Pose0.Pos() = w2x * PHS0->GetCenterPosition();
		Pose1.Ori() = w2x * PHS1->GetOrientation();
		Pose1.Pos() = w2x * PHS1->GetCenterPosition();

		Vec3d cp2ctw[2]; // CenterPositionToConTactWorld
		for (int i = 0; i < 2; i++) {
			cp2ctw[i] = pose.Pos() - solid[i]->GetCenterPosition();	//剛体の中心から接触点までのベクトル(world座標系)
		}

		// local: 接触点の関節フレーム は，x軸を法線, y,z軸を接線とする

		Quaterniond qlocal;
		qlocal.FromMatrix(local);
		for (int i = 0; i < 2; i++) {
			(i == 0 ? poseSocket : posePlug).Ori() = Xj[i].q = solid[i]->GetOrientation().Conjugated() * qlocal;
			(i == 0 ? poseSocket : posePlug).Pos() = Xj[i].r = solid[i]->GetOrientation().Conjugated() * cp2ctw[i];
		}

		// 6自由度拘束必要
		movableAxes.Clear();

	}

	bool PHContactEngine::Iterate() {

		bool updated = false;

		for (int i = 0; i < 6; i++) {

			// Gauss-Seidel Update
			dv[i] = Dot6((const double*)J[0].row(i), (const double*)solid[0]->dv)
				+ Dot6((const double*)J[1].row(i), (const double*)solid[1]->dv);
			res[i] = b[i] + dA[i] * f[i] + dv[i];
			fnew[i] = f[i] - Ainv[i] * res[i];

			/*
			// Comp Response & Update f
			df[i] = fnew[i] - f[i];
			f[i] = fnew[i];

			if (std::abs(df[i]) > engine->dfEps) {
			CompResponse(df[i], i);

			}
			*/
		}

		// 法線での力
		Projection(fnew, 0, updated);
		// 接線でのトルク
		Projection(fnew, 4, updated);
		Projection(fnew, 5, updated);
		// 法線での力
		Projection(fnew, 1, updated);
		Projection(fnew, 2, updated);
		// 接線でのトルク
		Projection(fnew, 3, updated);

		return updated;
	}

	void PHContactEngine::Projection(SpatialVector & fnew, const int i, bool& updat)
	{
		switch (i) {

		case 0: // 垂直方向（接触法線方向）
		{
			if (fnew[0] < EPSILON_10) {
				fnew[0] = 0.0;
			}

			df[i] = fnew[i] - f[i];
			f[i] = fnew[i];

			if (std::abs(df[i]) > engine->dfEps) {
				updat = true;
				CompResponse(df[i], i);
			}

			break;

		}

		case 1: // 接触面方向（normalと垂直)
		{
			double vmax2cont = sqrt(vjrel[1] * vjrel[1] + vjrel[2] * vjrel[2]);
			double vmax2ideal = GetScene()->GetFrictionThreshold();

			double ftmax = (vmax2cont < vmax2ideal ? mu0 : mu) * fnew[0];
			double ft = sqrt(fnew[1] * fnew[1] + fnew[2] * fnew[2]);

			frictionMargin = ftmax - ft;

			if (vmax2cont < vmax2ideal) {
				if (ft > ftmax) {
					ft = ftmax;
					frictionMargin = 0.0;
				}
				else if (ft < -ftmax) {
					ft = -ftmax;
					frictionMargin = 0.0;
				}
			}
			else {
				if (ft > ftmax) {
					ft = ftmax;
					frictionMargin = 0.0;
				}
				else if (ft < -ftmax) {
					ft = -ftmax;
					frictionMargin = 0.0;
				}
			}

			/*
			if (frictionMargin < 0.0) {
			double k = ftmax / ft;
			fnew[1] *= k;
			fnew[2] *= k;
			frictionMargin = 0.0;
			}
			*/

			df[i] = fnew[i] - f[i];
			f[i] = fnew[i];

			if (std::abs(df[i]) > engine->dfEps) {
				updat = true;
				CompResponse(df[i], i);
			}

			break;
		}
		case 4://接触法線に作用するトルク
		{

			copMargin = EPSILON_10;

			if (fnew[0] < EPSILON_10) {

				fnew[4] = 0, fnew[5] = 0;

				df[4] = fnew[4] - f[4];
				f[4] = fnew[4];

				if (std::abs(df[4]) > engine->dfEps) {
					updat = true;
					CompResponse(df[4], 4);
				}

				break;
			}

			Vec3d cop2cont = Vec3d(0.0, -fnew[5] / fnew[0], fnew[4] / fnew[0]);

			std::cout << "Contact Point : " << ContPoint << '\n';

			OldCoP = ContPoint + ContLocal * cop2cont;

			std::cout << "Old CoP : " << OldCoP << '\n';

			if (HowToContact(CDC0, CDC1) == CONTACT::POINT) {

				NewCoP = Vec3d();
				copMargin = 0.0;

				fnew[4] = NewCoP.Z() * fnew[0];
				fnew[5] = -NewCoP.Y() * fnew[0];

			}
			else {
				// 点が中にあれば処理が必要ない
				if (InOut(OldCoP, CDC0, CDC1)) {

					copMargin = 0.0;

				}

				else {

					Vec3f w = Vec3f();

					FindClosestPoint(w, OldCoP, HowToContact(CDC0, CDC1));

					NewCoP = ContLocal.inv() * (w - ContPoint);
					copMargin = 0.0;

					fnew[4] = NewCoP.Z() * fnew[0];
					fnew[5] = -NewCoP.Y() * fnew[0];

				}
			}

			if (copMargin != 0.0) {
				Vec3d c = OldCoP - ContPoint;
				copMargin = std::min(copMargin, -c * ContLocal.Ex());
			}

			df[i] = fnew[i] - f[i];
			f[i] = fnew[i];

			if (std::abs(df[i]) > engine->dfEps) {
				updat = true;
				CompResponse(df[i], i);
			}

			break;
		}

		case 3:	//接触法線と垂直な方向のトルク
		{
			if (fnew[0] < EPSILON_10) {
				fnew[3] = 0.0;

				df[i] = fnew[i] - f[i];
				f[i] = fnew[i];

				if (std::abs(df[i]) > engine->dfEps) {
					updat = true;
					CompResponse(df[i], i);
				}

				break;
			}

			// 接触法線垂直向きのfnew1, fnew2から生じるモーメント
			// 接触法線垂直向きから生まれる力静止摩擦力のように扱う
			// 拘束している力であると考えている

			double ncop = NewCoP.y*fnew[2] - NewCoP.z*fnew[1];

			// ここの2/3の意味とあとMarginらが必要な理由というのがまだ理解できない
			double margin = (2.0 / 3.0)*copMargin*frictionMargin;
			//double margin = (2.0/3.0)*contactRadius*frictionMargin;

			if (fnew[3] - ncop > margin) fnew[3] = ncop + margin;
			if (fnew[3] - ncop < -margin) fnew[3] = ncop - margin;

			df[i] = fnew[i] - f[i];
			f[i] = fnew[i];

			if (std::abs(df[i]) > engine->dfEps) {
				updat = true;
				CompResponse(df[i], i);
			}

			break;
		}

		default:
		{
			df[i] = fnew[i] - f[i];
			f[i] = fnew[i];

			if (std::abs(df[i]) > engine->dfEps) {
				updat = true;
				CompResponse(df[i], i);
			}

			break;
		}
		}
	}

	// 近い点を見つけるとどうように接触状態を把握する必要がある
	// 面と面の接触状態であれば→ZMPの位置を算出する
	// 支持多角形が存在しない→点と何かの接触状態である
	// edge接触であるならば?→edge接触

	// 簡単なアルゴリズムの流れ
	/*
	if(接触状態 = Face, Edge)
	calculate constrainted COP
	if(COPが物体の中に属していれば)
	そのまま計算する
	else ←COPが物体の外に属していれば
	COPが一番近い点に近づけるように計算する
	else if(接触状態 = Pointたちであれば)
	ZMPが接触ポイントに近づけるように計算する

	重要関数　
	1. 接触状態関数
	2. COP内外判定関数
	3. 一番近い点を見つけてくれる
	*/

	/* 3に該当する関数 */


	// OldCoPが外に出ているときに呼ばれる関数
	bool PHContactEngine::FindClosestPoint(Vec3f & cp, const Vec3d & OldCoP, CONTACT cts) {

		double t = 0.0;

		if (cts == CONTACT::LINE) {

			if (SupportTag0 == CONTACT::LINE && SupportTag1 == CONTACT::LINE) {

				Vec3d zp_0 = ClosestPtPointLine(OldCoP, VecGroup0[0], VecGroup0[1], t);
				Vec3d zp_1 = ClosestPtPointLine(OldCoP, VecGroup1[0], VecGroup1[1], t);

				if (IsInside(CDC1, PHS1, zp_0)) {
					cp = zp_0;
					return true;
				}
				else if (IsInside(CDC0, PHS0, zp_1)) {
					cp = zp_1;
					return true;
				}
				else {
					//EuclideanProjection(SupportTag0, SupportTag1, zp_1, cp);
					GJKIteration(SupportTag0, SupportTag1, zp_0, zp_1, cp);
					return true;
				}
			}
			else if (SupportTag0 == CONTACT::FACE && SupportTag1 == CONTACT::LINE) {

				Vec3f zp_0;
				GJK2D(VecGroup0, CDC0, PHS0, zp_0, OldCoP);
				Vec3d zp_1 = ClosestPtPointLine(OldCoP, VecGroup1[0], VecGroup1[1], t);

				if (IsInside(CDC1, PHS1, zp_0)) {
					cp = zp_0;
					return true;
				}

				else if (IsInside(CDC0, PHS0, zp_1)) {
					cp = zp_1;
					return true;
				}
				else {
					//EuclideanProjection(SupportTag0, SupportTag1, zp_1, cp);
					GJKIteration(SupportTag0, SupportTag1, zp_0, zp_1, cp);
					return true;
				}

			}
			else if (SupportTag0 == CONTACT::LINE && SupportTag1 == CONTACT::FACE) {

				Vec3d zp_0 = ClosestPtPointLine(OldCoP, VecGroup0[0], VecGroup0[1], t);
				Vec3f zp_1;
				GJK2D(VecGroup1, CDC1, PHS1, zp_1, OldCoP);

				if (IsInside(CDC1, PHS1, zp_0)) {
					cp = zp_0;
					return true;
				}

				else if (IsInside(CDC0, PHS0, zp_1)) {
					cp = zp_1;
					return true;
				}
				else {
					//EuclideanProjection(SupportTag0, SupportTag1, zp_0, cp);
					GJKIteration(SupportTag0, SupportTag1, zp_0, zp_1, cp);
					return true;
				}
			}
		}
		// 面と面の接触であれば
		else {

			Vec3f zp_0;
			GJK2D(VecGroup0, CDC0, PHS0, zp_0, OldCoP);
			Vec3f zp_1;
			GJK2D(VecGroup1, CDC1, PHS1, zp_1, OldCoP);

			if (IsInside2D(VecGroup1, CDC1, PHS1, zp_0)) {
				cp = zp_0;
				return true;
			}

			else if (IsInside2D(VecGroup0, CDC0, PHS0, zp_1)) {
				cp = zp_1;
				return true;
			}
			// Dynamical運動を考えると必ずしもZMPがうまい位置にとれることはないのでこの関数を設けている
			else {
				//EuclideanProjection(SupportTag0, SupportTag1, zp_1, cp);
				GJKIteration(SupportTag0, SupportTag1, zp_0, zp_1, cp);
				return true;
			}
		}

		// 失敗したらfalse
		cp = Vec3d();
		return false;

	}

	bool PHContactEngine::GJK2D(const std::vector<Vec3d>& c_vec, CDConvex* c, PHSolid* s, Vec3f& tmp, const Vec3d& OldCoP) {

		size_t vec_size = c_vec.size();

		if (vec_size == 1) {
			c->Support(tmp, s->GetOrientation().Conjugated() * (OldCoP - s->GetCenterPosition()));
			tmp = s->GetOrientation() * tmp + s->GetCenterPosition();
			return true;
		}

		Posed pose_;
		pose_.Ori() = w2x * s->GetOrientation();
		pose_.Pos() = w2x * s->GetCenterPosition();

		Vec3d ct2x = pose_.Ori() * OldCoP;

		std::vector<Vec3d> tmp_vec(vec_size);
		for (unsigned int i = 0; i < vec_size; i++) {
			tmp_vec[i] = pose_.Ori() * c_vec[i];
		}

		Vec3d center = ct2x - tmp_vec[0];

		if (InOut2D(tmp_vec, ct2x) == true) {
			tmp = ClosestPtFromGroup(OldCoP, c_vec);
			return true;
		}

		else {
			// 面接触のときにだけ現れるので、もしc_vecの数が1個であるということは曲線の形を持つ2次元のものである

			std::vector<Vec3d> simplex;
			Vec3d tmp_0, tmp_1, tmp_2, tmp_3;
			while (1) {
				switch (simplex.size()) {
				case 0:
				{
					tmp_0 = tmp_vec[PointSupport(tmp_vec, center)];
					simplex.push_back(tmp_0);
					break;
				}
				case 1:
				{
					tmp_1 = ct2x - tmp_0;
					tmp_1 = tmp_vec[PointSupport(tmp_vec, tmp_1)];
					if (simplex[0] == tmp_1) {
						tmp = pose_.Ori().Conjugated() * tmp_1;
						return true;
					}
					simplex.push_back(tmp_1);
					break;
				}
				case 2:
				{
					double t = 0.0;
					Vec3d tmp__ = ClosestPtPointLine(ct2x, simplex[0], simplex[1], t);

					tmp_2 = ct2x - tmp__;
					tmp_2 = tmp_vec[PointSupport(tmp_vec, tmp_2)];

					if (tmp_2 == simplex[0] || tmp_2 == simplex[1]) {
						tmp = pose_.Ori().Conjugated() * tmp__;
						return true;
					}

					simplex.push_back(tmp_2);
					break;
				}
				default:
				{
					Vec3d tmp__ = ClosestPtPointTriangle(ct2x, simplex[0], simplex[1], simplex[2]);

					tmp_3 = ct2x - tmp__;
					tmp_3 = tmp_vec[PointSupport(tmp_vec, tmp_3)];

					if (tmp_3 == simplex[0] || tmp_3 == simplex[1] || tmp_3 == simplex[2]) {
						tmp = pose_.Ori().Conjugated() * tmp__;
						return true;
					}

					simplex.erase(simplex.begin());
					simplex.push_back(tmp_3);
					break;
				}
				}
			}
		}
	}

	bool PHContactEngine::GJK3D(CDConvex * cdc0, PHSolid * phs0, Vec3d & pt0, CDConvex * cdc1, PHSolid * phs1, Vec3d & pt1)
	{
		// support point
		// 初期化
		Vec3d init = phs1->GetCenterPosition() - phs0->GetCenterPosition();
		cdc0->Support((Vec3f)pt0, init);
		cdc1->Support((Vec3f)pt1, -init);
		return false;
	}

	// この関数は改良も必要だしそもそも精度が非常に重要
	// return 1 -> 点接触
	// return 2 -> 線接触
	// reteurn 3 -> 面接触
	/* 1に該当する関数 */
	CONTACT PHContactEngine::HowToContact(CDConvex * s0, CDConvex * s1)
	{
		VecGroup0.clear();
		VecGroup1.clear();

		Vec3d tmp = ContPoint - PHS1->GetCenterPosition();
		if (tmp.dot(ContLocal.Ex()) < 0) {
			Vec3d tmp_0 = PHS0->GetOrientation().Conjugated() * ContLocal.Ex();
			SupportTag0 = static_cast<CONTACT>(s0->SupportTag(tmp_0, VecGroup0));
			for (unsigned int i = 0; i < VecGroup0.size(); i++)
				VecGroup0[i] = PHS0->GetOrientation() * VecGroup0[i] + PHS0->GetCenterPosition();

			Vec3d tmp_1 = PHS1->GetOrientation().Conjugated() * -ContLocal.Ex();
			SupportTag1 = static_cast<CONTACT>(s1->SupportTag(tmp_1, VecGroup1));
			for (unsigned int i = 0; i < VecGroup1.size(); i++)
				VecGroup1[i] = PHS1->GetOrientation() * VecGroup1[i] + PHS1->GetCenterPosition();

		}
		else {
			Vec3d tmp_0 = PHS0->GetOrientation().Conjugated() * -ContLocal.Ex();
			SupportTag0 = static_cast<CONTACT>(s0->SupportTag(tmp_0, VecGroup0));
			for (unsigned int i = 0; i < VecGroup0.size(); i++)
				VecGroup0[i] = PHS0->GetOrientation() * VecGroup0[i] + PHS0->GetCenterPosition();

			Vec3d tmp_1 = PHS1->GetOrientation().Conjugated() * ContLocal.Ex();
			SupportTag1 = static_cast<CONTACT>(s1->SupportTag(tmp_1, VecGroup1));
			for (unsigned int i = 0; i < VecGroup1.size(); i++)
				VecGroup1[i] = PHS1->GetOrientation() * VecGroup1[i] + PHS1->GetCenterPosition();
		}

		//点接触
		if (SupportTag0 == CONTACT::POINT || SupportTag1 == CONTACT::POINT) {
			return CONTACT::POINT;
		}

		// 面とEDGEであるかEdgeとEdgeであるか
		else if (SupportTag0 == CONTACT::LINE || SupportTag1 == CONTACT::POINT) {
			return CONTACT::LINE;
		}

		// 面と面の接触だけ
		else if (SupportTag0 == CONTACT::FACE && SupportTag1 == CONTACT::FACE) {
			return CONTACT::FACE;
		}

		return CONTACT::UNKNOWN;
	}

	// ZMPが物体の中にあるか外にあるかを判断する
	// ただの内外判定だけではない
	// 内外判定以外にどのように接触しているかの情報
	// 中にあれば true
	// 外にあるとき false
	/* 2に該当する関数 */
	inline bool PHContactEngine::InOut(const Vec3d & pt, CDConvex* CDC0, CDConvex* CDC1)
	{
		if (CDC0->IsInside(PHS0->GetOrientation().Conjugated() * (pt - PHS0->GetCenterPosition())) && CDC1->IsInside(PHS1->GetOrientation().Conjugated() * (pt - PHS1->GetCenterPosition())))
			return true;
		return false;
	}

	inline bool	PHContactEngine::IsInside(CDConvex* c, PHSolid* s, const Vec3d& pt) {
		if (c->IsInside(s->GetOrientation().Conjugated() * (pt - s->GetCenterPosition())))
			return true;
		return false;
	}

	// 中にあればTrue
	// 外にあればFalse
	bool PHContactEngine::IsInside2D(const std::vector<Vec3d>& c_vec, CDConvex * c, PHSolid * s, const Vec3d & pt)
	{
		Posed pose_;
		pose_.Ori() = w2x * s->GetOrientation();
		pose_.Pos() = w2x * s->GetCenterPosition();

		Vec3d ct2x = pose_.Ori() * pt;

		std::vector<Vec3d> tmp_vec(c_vec.size());
		for (unsigned int i = 0; i < tmp_vec.size(); i++) {
			tmp_vec[i] = pose_.Ori() * c_vec[i];
		}

		if (tmp_vec.size() != 1) {
			return InOut2D(tmp_vec, ct2x);
		}
		else {
			Vec3f w;
			c->Support(w, s->GetOrientation().Conjugated() * (pt - s->GetCenterPosition()));
			w = s->GetOrientation() * w + s->GetCenterPosition();
			if (((pose_.Ori() * pt).YZ()).norm() - ((pose_.Ori() * w).YZ()).norm() >= 0)
				return false;
			else
				return true;
		}
	}

	bool PHContactEngine::GJKIteration(int SupportTag0, int SupportTag1, const Vec3f& z_0, const Vec3f& z_1, Vec3f& tmp) {

		Vec3d pt = OldCoP;
		Vec3d tmp_0 = z_0;
		Vec3d tmp_1 = z_1;
		Vec3f tmp_tmp, tmp_tmp_0, tmp_tmp_1;

		// 切り替えながら判定するための変数
		int cross = 0;
		double t = 0.0;

		double stop = Norm2D(Pose0, tmp_0, Pose1, tmp_1);
		if (SupportTag0 == 2 && SupportTag1 == 2) {

			tmp_tmp_0 = z_0;
			tmp_tmp_1 = z_1;

			while (stop > EPSILON_2) {
				tmp_tmp_0 = ClosestPtPointLine(tmp_1, VecGroup0[0], VecGroup0[1], t);
				if (Norm2D(Pose0, tmp_tmp_0, Pose1, tmp_tmp_1) > EPSILON_2) {
					tmp_0 = tmp_tmp_0;
					break;
				}
				tmp_tmp_1 = ClosestPtPointLine(tmp_0, VecGroup1[0], VecGroup1[1], t);
				if (Norm2D(Pose0, tmp_tmp_0, Pose1, tmp_tmp_1) > EPSILON_2) {
					tmp_1 = tmp_tmp_1;
					break;
				}
				tmp_0 = tmp_tmp_0;
				tmp_1 = tmp_tmp_1;
				stop = Norm2D(Pose0, tmp_0, Pose1, tmp_1);
			}

			tmp = (tmp_0 + tmp_1) / 2;
			return true;

		}
		else if (SupportTag0 == 3 && SupportTag1 == 2) {
			while (stop > EPSILON_2) {

				GJK2D(VecGroup0, CDC0, PHS0, tmp_tmp_0, tmp_1);
				if (Norm2D(Pose0, tmp_tmp_0, Pose1, tmp_tmp_1) > EPSILON_2) {
					tmp_0 = tmp_tmp_0;
					break;
				}
				tmp_tmp_1 = ClosestPtPointLine(tmp_0, VecGroup1[0], VecGroup1[1], t);
				if (Norm2D(Pose0, tmp_tmp_0, Pose1, tmp_tmp_1) > EPSILON_2) {
					tmp_1 = tmp_tmp_1;
					break;
				}
				tmp_0 = tmp_tmp_0;
				tmp_1 = tmp_tmp_1;
				stop = Norm2D(Pose0, tmp_0, Pose1, tmp_1);
			}

			tmp = (tmp_0 + tmp_1) / 2;
			return true;
		}
		else if (SupportTag0 == 2 && SupportTag1 == 3) {
			while (stop > EPSILON_2) {

				GJK2D(VecGroup1, CDC1, PHS1, tmp_tmp_1, tmp_0);
				if (Norm2D(Pose0, tmp_tmp_0, Pose1, tmp_tmp_1) > EPSILON_2) {
					tmp_1 = tmp_tmp_1;
					break;
				}
				tmp_tmp_0 = ClosestPtPointLine(tmp_1, VecGroup0[0], VecGroup0[1], t);
				if (Norm2D(Pose0, tmp_tmp_0, Pose1, tmp_tmp_1) > EPSILON_2) {
					tmp_0 = tmp_tmp_0;
					break;
				}
				tmp_0 = tmp_tmp_0;
				tmp_1 = tmp_tmp_1;
				stop = Norm2D(Pose0, tmp_0, Pose1, tmp_1);
			}

			tmp = (tmp_0 + tmp_1) / 2;
			return true;

		}
		else {
			while (stop > EPSILON_2) {

				GJK2D(VecGroup0, CDC0, PHS0, tmp_tmp_0, tmp_1);
				if (Norm2D(Pose0, tmp_tmp_0, Pose1, tmp_tmp_1) > EPSILON_2) {
					tmp_0 = tmp_tmp_0;
					break;
				}
				GJK2D(VecGroup1, CDC1, PHS1, tmp_tmp_1, tmp_0);
				if (Norm2D(Pose0, tmp_tmp_0, Pose1, tmp_tmp_1) > EPSILON_2) {
					tmp_1 = tmp_tmp_1;
					break;
				}
				tmp_0 = tmp_tmp_0;
				tmp_1 = tmp_tmp_1;
				stop = Norm2D(Pose0, tmp_0, Pose1, tmp_1);
			}

			tmp = (tmp_0 + tmp_1) / 2;
			return true;
		}
	}

	// Euclidean Projection
	void PHContactEngine::EuclideanProjection(CONTACT SupportTag0, CONTACT SupportTag1, const Vec3f& z_0, Vec3f& tmp) {

		// 移動前の点
		Vec3d prev_point = OldCoP;
		// 移動後の点
		Vec3d current_point = z_0;

		// 切り替えながら判定するための変数
		int point_on_where = 0;
		double t = 0.0;

		// 停止条件
		double stop = Norm3D(prev_point, current_point);

		if (SupportTag0 == CONTACT::LINE && SupportTag1 == CONTACT::LINE) {

			while (stop > EPSILON_2) {
				prev_point = current_point;
				current_point = point_on_where == 0 ? ClosestPtPointLine(current_point, VecGroup1[0], VecGroup1[1], t) :
					ClosestPtPointLine(current_point, VecGroup0[0], VecGroup0[1], t);
				// current_pointがどこの物体に対して引かれたをずっとビットを反対させるため
				point_on_where = point_on_where ^ 1;
				stop = Norm3D(prev_point, current_point);
			}

			tmp = current_point;

		}
		else if (SupportTag0 == CONTACT::FACE && SupportTag1 == CONTACT::LINE) {

			Vec3d tmp_current_point;

			while (stop > EPSILON_2) {
				prev_point = current_point;
				if (point_on_where == 0) current_point = ClosestPtPointLine(current_point, VecGroup1[0], VecGroup1[1], t);
				else {
					GJK2D(VecGroup0, CDC0, PHS0, (Vec3f)tmp_current_point, current_point);
					current_point = tmp_current_point;
				}
				point_on_where = point_on_where ^ 1;
				stop = Norm3D(prev_point, current_point);
			}

			tmp = current_point;
		}
		else if (SupportTag0 == CONTACT::LINE && SupportTag1 == CONTACT::FACE) {

			Vec3d tmp_current_point;

			while (stop > EPSILON_2) {
				prev_point = current_point;
				if (point_on_where == 1) current_point = ClosestPtPointLine(current_point, VecGroup0[0], VecGroup0[1], t);
				else {
					GJK2D(VecGroup1, CDC1, PHS1, (Vec3f)tmp_current_point, current_point);
					current_point = tmp_current_point;
				}
				point_on_where = point_on_where ^ 1;
				stop = Norm3D(prev_point, current_point);
			}

			tmp = current_point;

		}
		else
		{
			while (stop > EPSILON_2) {

				Vec3d tmp_current_point;

				while (stop > EPSILON_2) {
					prev_point = current_point;
					if (point_on_where == 1) {
						GJK2D(VecGroup0, CDC0, PHS0, (Vec3f)tmp_current_point, current_point);
						current_point = tmp_current_point;
					}
					else {
						GJK2D(VecGroup1, CDC1, PHS1, (Vec3f)tmp_current_point, current_point);
						current_point = tmp_current_point;
					}
					point_on_where = point_on_where ^ 1;
					stop = Norm3D(prev_point, current_point);
				}

				tmp = current_point;
			}
		}
	}
}