/*
Convex物体に対しての接触力を計算するためのプログラミングとなります
*/

/*
made by KIM
*/

#ifndef PHCONTACTFORCYL_H
#define PHCONTACTFORCYL_H

#include <Physics/PHContactPoint.h>

namespace Spr {
	;
	class PHShapePairForLCP;

	class PHContactCylinder : public PHContactPoint {
	private:
		/// ZMP更新前
		Vec3d old_cop = Vec3d();
		/// ZMP更新後
		Vec3d new_cop = Vec3d();

		/// サポートタグによる点の候補
		std::vector<Vec3d> can_vec_0;
		std::vector<Vec3d> can_vec_1;
		
		/// サポートタグによる値
		int st_0, st_1;

		double  frictionMargin;   ///< 最大摩擦力と作用している摩擦力とのマージン
		double  copMargin;        ///< CoPから接触多角形境界までのマージン

		/// ZMP位置算出用
		Vec3d cont_point;
		PHSolid* a0;
		PHSolid* a1;
		CDConvex* c0;
		CDConvex* c1;
		Matrix3d cont_local;
		Quaterniond w2x;
		Posed pose_0, pose_1;

		// 3次元の線に一番近い点を探すものであり、出力として３次元の座標となる
		static Vec3d	ClosestPtPointLine(const Vec3d& old_ZMP, const Vec3d& tmp_1, const Vec3d& tmp_2, double& t);
		// ３次元の三角形に対して一番近い点を見つけ出し、出力として３次元の座標とする
		static Vec3d	ClosestPtPointTriangle(const Vec3d& old_ZMP, const Vec3d& tmp_1, const Vec3d& tmp_2, const Vec3d& tmp_3);
		// 3次元の線に一番近い点を探し、出力として、一番近い点から入力としていれた点までのベクトルを返す
		static Vec3d	NearVecOnLine(const Vec3d& old_ZMP, const Vec3d& tmp_1, const Vec3d& tmp_2);
		// 2DのConvex hullとある2Dのある点を入力とし、その点がConvex hullの中に入っているかどうかを判断しboolean値を返す
		static	bool	InOut2D(const std::vector<Vec3d>& tmp_vec, const Vec3d& pt);
		// 2DのConvex hullと、あるベクトルを入れることでサポートポイントの点の座標を取り出し、出力はその座標の番号となる
		static int		PointSupport(const std::vector<Vec3d>& tmp_vec, const Vec3d& normal);
		// ３D上の二つの点に対して2D上の座標に変換してからその距離を返す
		inline double	Norm2D(const Vec3d& tmp_0, const Vec3d& tmp_1);
		// GJKのバリエーション
		bool			GJKIteration(int st_0, int st_1, const Vec3f& z0, const Vec3f& z1, Vec3f& tmp);
		// どのように接触しているかを判断するものであり、Integer型としてその接触の形を分類する
		int				HowToContact(CDConvex* s0, CDConvex* s1);
		// 点Ptがc0とc1に属しているかどうかを判断し、Boolean値を返す
		inline bool		InOut(const Vec3d& pt, CDConvex* c0, CDConvex* c1);
		// ptがあるものに属しているかどうかを判断し、Boolean値を返す
		inline bool		IsInside(CDConvex* c, PHSolid* s, const Vec3d& pt);
		// 
		bool			IsInside2D(const std::vector<Vec3d>& c_vec, CDConvex* c, PHSolid* s, const Vec3d& pt);
		Vec3d			FindOnEdge(const Vec3d& old_ZMP, const std::vector<Vec3d>& c_vec);

	public:

		/// コンストラクタ
		PHContactCylinder() {}
		PHContactCylinder(const Matrix3d& local, PHShapePairForLCP* sp, Vec3d p, PHSolid* s0, PHSolid* s1);
		
		// ----- PHConstraintの機能をオーバーライド
		virtual bool Iterate();

		// ----- このクラスで実装する機能
		bool	FindClosestPoint(Vec3f& cp, const Vec3d& old_cop, int cts);
		void	Projection(SpatialVector& fnew, const int i, bool& updat);
		bool	GJK2D(const std::vector<Vec3d>& c_vec, CDConvex* c, PHSolid* s, Vec3f& tmp, const Vec3d& old_cop);

	};

}

#endif
