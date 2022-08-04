#ifndef CGOGN_SIMULATION_INTEGRATION_SCHEMES_H
#define CGOGN_SIMULATION_INTEGRATION_SCHEMES_H

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{
namespace simulation
{

template <typename Mesh>
struct explicit_euler
{
	/**
	 * x(t+dt) = x(t)+dt*v(t)
	 * v(t+dt) = v(t)+dt*Fext(x(t),v(t))
	 */
	template <typename T>
	using Attribute = typename mesh_traits<Mesh>::template Attribute<T>;
	using Vec3 = geometry::Vec3;

	template <typename FUNC, typename... Params, typename CELL>
	void integrate(const Mesh& m, const FUNC& force, double time_step, Params&... params)
	{
		force();
		integration_aux(m, time_step, params...);
	}

	template <typename FUNC, typename... Params, typename CELL>
	void integrate_aux(const Mesh& m, double time_step, std::vector<CELL>& particules, Attribute<Vec3>* forces_ext,
					   Attribute<Vec3>* speed, Attribute<Vec3>* position, Attribute<double>* masse, Params&... params)
	{
		integrate_aux(m, time_step, particules, forces_ext, speed, position, masse);
		integrate_aux<Params...>(m, time_step, params...);
	}

	template <typename CELL>
	void integrate_aux(const Mesh& m, double time_step, std::vector<CELL>& particules, Attribute<Vec3>* forces_ext,
					   Attribute<Vec3>* speed, Attribute<double>* masse, Attribute<Vec3>* position)
	{
		for (auto& c : particules)
		{
			value<Vec3>(m, position, c) = value<Vec3>(m, position, c) + time_step * value<Vec3>(m, speed, c);
			value<Vec3>(m, speed, c) = 0.995 * value<Vec3>(m, speed, c) +
									   time_step * value<Vec3>(m, forces_ext, c) / value<double>(m, masse, c);

			value<Vec3>(m, forces_ext, c) = Vec3(0, 0, 0);
		}
	}
};

template <typename Mesh>
struct symplectic_euler
{
	/**
	 * x(t+dt) = x(t)+dt*v(t+dt)
	 * v(t+dt) = v(t)+dt*Fext(x(t),v(t))
	 */
	template <typename T>
	using Attribute = typename mesh_traits<Mesh>::template Attribute<T>;
	using Vec3 = geometry::Vec3;

	template <typename FUNC, typename... Params, typename CELL>
	void integrate(const Mesh& m, const FUNC& force, double time_step, Params&... params)
	{
		force();
		integration_aux(m, time_step, params...);
	}

	template <typename FUNC, typename... Params, typename CELL>
	void integrate_aux(const Mesh& m, double time_step, std::vector<CELL>& particules, Attribute<Vec3>* forces_ext,
					   Attribute<Vec3>* speed, Attribute<Vec3>* position, Attribute<double>* masse, Params&... params)
	{
		integrate_aux(m, time_step, particules, forces_ext, speed, position, masse);
		integrate_aux<Params...>(m, time_step, params...);
	}

	template <typename CELL>
	void integrate_aux(const Mesh& m, double time_step, std::vector<CELL>& particules, Attribute<Vec3>* forces_ext,
					   Attribute<Vec3>* speed, Attribute<double>* masse, Attribute<Vec3>* position)
	{
		for (auto& c : particules)
		{
			value<Vec3>(m, speed, c) = 0.995 * value<Vec3>(m, speed, c) +
									   time_step * value<Vec3>(m, forces_ext, c) / value<double>(m, masse, c);
			value<Vec3>(m, position, c) = value<Vec3>(m, position, c) + time_step * value<Vec3>(m, speed, c);
			value<Vec3>(m, forces_ext, c) = Vec3(0, 0, 0);
		}
	}
};

template <typename Mesh>
struct Runge_Kutta4
{
	template <typename T>
	using Attribute = typename mesh_traits<Mesh>::template Attribute<T>;
	using Vec3 = geometry::Vec3;

	template <typename FUNC, typename... Params, typename CELL>
	void integrate(const Mesh& m, const FUNC& force, double time_step, Params&... params)
	{
		force();
		integration_part1(m, time_step, params...);
		force();
		integration_part2(m, time_step, params...);
		force();
		integration_part3(m, time_step, params...);
		force();
		integration_part4(m, time_step, params...);
	}

	template <typename FUNC, typename... Params, typename CELL>
	void integrate_part1(const Mesh& m, double time_step, std::vector<CELL>& particules, Attribute<Vec3>* forces_ext,
						 Attribute<Vec3>* speed, Attribute<Vec3>* position, Attribute<double>* masse, Params&... params)
	{
		integrate_part1(m, time_step, particules, forces_ext, speed, position, masse);
		integrate_part1<Params...>(m, time_step, params...);
	}

	template <typename CELL>
	void integrate_part1(const Mesh& m, double time_step, std::vector<CELL>& particules, Attribute<Vec3>* forces_ext,
						 Attribute<Vec3>* speed, Attribute<double>* masse, Attribute<Vec3>* position)
	{

		std::shared_ptr<Attribute<std::array<Vec3, 9>>> RK_coeff;
		RK_coeff = get_attribute<std::array<Vec3, 9>, CELL>(m, "Integration_RK_coeff");
		if (RK_coeff == nullptr)
			RK_coeff = add_attribute<std::array<Vec3, 9>, CELL>(m, "Integration_RK_coeff");
		for (auto& c : particules)
		{
			// k1
			value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[0] =
				(time_step * value<Vec3>(m, forces_ext, c) / value<double>(m, masse, c));
			// j1
			value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[1] =
				time_step * (0.995 * value<Vec3>(m, speed, c) + value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[0]);

			value<Vec3>(m, position, c) += value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[1] / 2.0f;

			value<Vec3>(m, this->forces_ext_.get(), c) = value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[8];
		}
	}

	template <typename FUNC, typename... Params, typename CELL>
	void integrate_part2(const Mesh& m, double time_step, std::vector<CELL>& particules, Attribute<Vec3>* forces_ext,
						 Attribute<Vec3>* speed, Attribute<Vec3>* position, Attribute<double>* masse, Params&... params)
	{
		integrate_part2(m, time_step, particules, forces_ext, speed, position, masse);
		integrate_part2<Params...>(m, time_step, params...);
	}

	template <typename CELL>
	void integrate_part2(const Mesh& m, double time_step, std::vector<CELL>& particules, Attribute<Vec3>* forces_ext,
						 Attribute<Vec3>* speed, Attribute<double>* masse, Attribute<Vec3>* position)
	{

		std::shared_ptr<Attribute<std::array<Vec3, 9>>> RK_coeff;
		RK_coeff = get_attribute<std::array<Vec3, 9>, CELL>(m, "Integration_RK_coeff");
		if (RK_coeff == nullptr)
			RK_coeff = add_attribute<std::array<Vec3, 9>, CELL>(m, "Integration_RK_coeff");
		for (auto& c : particules)
		{
			// k2
			value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[2] =
				(time_step * value<Vec3>(m, forces_ext, c) / value<double>(m, masse, c));
			// j2
			value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[3] =
				time_step *
				(0.995 * value<Vec3>(m, speed, c) + value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[2] / 2.0f);

			value<Vec3>(m, position, c) -= value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[1] / 2.0f;

			value<Vec3>(m, position, c) += value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[3] / 2.0f;

			value<Vec3>(m, forces_ext, c) = value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[8];
		}
	}

	template <typename FUNC, typename... Params, typename CELL>
	void integrate_part3(const Mesh& m, double time_step, std::vector<CELL>& particules, Attribute<Vec3>* forces_ext,
						 Attribute<Vec3>* speed, Attribute<Vec3>* position, Attribute<double>* masse, Params&... params)
	{
		integrate_part3(m, time_step, particules, forces_ext, speed, position, masse);
		integrate_part3<Params...>(m, time_step, params...);
	}

	template <typename CELL>
	void integrate_part3(const Mesh& m, double time_step, std::vector<CELL>& particules, Attribute<Vec3>* forces_ext,
						 Attribute<Vec3>* speed, Attribute<double>* masse, Attribute<Vec3>* position)
	{

		std::shared_ptr<Attribute<std::array<Vec3, 9>>> RK_coeff;
		RK_coeff = get_attribute<std::array<Vec3, 9>, CELL>(m, "Integration_RK_coeff");
		if (RK_coeff == nullptr)
			RK_coeff = add_attribute<std::array<Vec3, 9>, CELL>(m, "Integration_RK_coeff");
		for (auto& c : particules)
		{
			// k3
			value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[4] =
				(time_step * value<Vec3>(m, forces_ext, c) / value<double>(m, masse, c));
			// j3
			value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[5] =
				time_step *
				(0.995 * value<Vec3>(m, speed, c) + value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[4] / 2.0f);

			value<Vec3>(m, position, c) -= value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[3] / 2.0f;

			value<Vec3>(m, position, c) += value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[5];

			value<Vec3>(m, forces_ext, c) = value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[8];
		}
	}

	template <typename FUNC, typename... Params, typename CELL>
	void integrate_part4(const Mesh& m, double time_step, std::vector<CELL>& particules, Attribute<Vec3>* forces_ext,
						 Attribute<Vec3>* speed, Attribute<Vec3>* position, Attribute<double>* masse, Params&... params)
	{
		integrate_part4(m, time_step, particules, forces_ext, speed, position, masse);
		integrate_part4<Params...>(m, time_step, params...);
	}

	template <typename CELL>
	void integrate_part4(const Mesh& m, double time_step, std::vector<CELL>& particules, Attribute<Vec3>* forces_ext,
						 Attribute<Vec3>* speed, Attribute<double>* masse, Attribute<Vec3>* position)
	{

		std::shared_ptr<Attribute<std::array<Vec3, 9>>> RK_coeff;
		RK_coeff = get_attribute<std::array<Vec3, 9>, CELL>(m, "Integration_RK_coeff");
		if (RK_coeff == nullptr)
			RK_coeff = add_attribute<std::array<Vec3, 9>, CELL>(m, "Integration_RK_coeff");
		for (auto& c : particules)
		{
			// k4
			value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[6] =
				(time_step * value<Vec3>(m, forces_ext, c) / value<double>(m, masse, c));
			// j4
			value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[7] =
				time_step * (0.995 * value<Vec3>(m, speed, c) + value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[6]);

			value<Vec3>(m, position, c) -= value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[5];

			Vec3 k1 = value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[0];
			Vec3 k2 = value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[2];
			Vec3 k3 = value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[4];
			Vec3 k4 = value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[6];

			Vec3 j1 = value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[1];
			Vec3 j2 = value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[3];
			Vec3 j3 = value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[5];
			Vec3 j4 = value<std::array<Vec3, 9>>(m, RK_coeff.get(), c)[7];

			Vec3 diff_pos = 1. / 6. * (j1 + (2 * j2) + (2 * j3) + j4);

			Vec3 diff_speed = 1. / 6. * (k1 + 2 * k2 + 2 * k3 + k4);

			if (diff_pos.norm() < 1.0e-10)
				diff_pos = Vec3(0, 0, 0);

			if (diff_speed.norm() < 1.0e-10)
				diff_speed = Vec3(0, 0, 0);

			value<Vec3>(m, position, c) += diff_pos;

			value<Vec3>(m, speed, c) = 0.995 * value<Vec3>(m, speed, c) + diff_speed;
			value<Vec3>(m, forces_ext, c) = Vec3(0, 0, 0);
		}
	}
};

} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_INTEGRATION_SCHEMES_H
