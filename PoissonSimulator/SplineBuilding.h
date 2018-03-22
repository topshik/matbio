#pragma once
#include "stdafx.h"

using namespace std;
using namespace alglib;
using namespace boost;

#ifndef SPLINE_BUILDING_H
#define SPLINE_BUILDING_H

namespace spline_building {

	spline1dinterpolant Build_death_splines(string death_kernel_string, double death_cutoff, int death_spline_nodes) {

		double death_r;

		exprtk::symbol_table<double>	death_kernel_symbol_table;
		exprtk::expression<double>		death_kernel_expression;
		exprtk::parser<double>			death_kernel_parser;

		death_kernel_symbol_table.add_variable("r", death_r);
		death_kernel_symbol_table.add_constants();
		death_kernel_expression.register_symbol_table(death_kernel_symbol_table);
		death_kernel_parser.compile(death_kernel_string, death_kernel_expression);

		real_1d_array x_1d_array;
		real_1d_array y_1d_array;

		x_1d_array.setlength(death_spline_nodes + 3);
		y_1d_array.setlength(death_spline_nodes + 3);

		for (int i = 0; i < death_spline_nodes + 3; i++) {
			x_1d_array[i] = death_cutoff*i / (death_spline_nodes + 3 - 1);
			death_r = x_1d_array[i];
			if (i < death_spline_nodes) {
				y_1d_array[i] = death_kernel_expression.value();
			}
			else {
				y_1d_array[i] = 0;
			}
		}
		spline1dinterpolant death_kernel_spline;
		spline1dbuildmonotone(x_1d_array, y_1d_array, death_kernel_spline);
		return death_kernel_spline;
	}

	spline1dinterpolant Build_birth_splines(string birth_kernel_string, double birth_cutoff, int birth_spline_nodes) {
		double birth_r;

		exprtk::symbol_table<double>	birth_kernel_symbol_table;
		exprtk::expression<double>		birth_kernel_expression;
		exprtk::parser<double>			birth_kernel_parser;

		birth_kernel_symbol_table.add_variable("r", birth_r);
		birth_kernel_symbol_table.add_constants();
		birth_kernel_expression.register_symbol_table(birth_kernel_symbol_table);
		birth_kernel_parser.compile(birth_kernel_string, birth_kernel_expression);

		real_1d_array x_1d_array;
		real_1d_array y_1d_array;

		x_1d_array.setlength(birth_spline_nodes + 3);
		y_1d_array.setlength(birth_spline_nodes + 3);

		for (int i = 0; i < birth_spline_nodes + 3; i++) {
			x_1d_array[i] = birth_cutoff*i / (birth_spline_nodes + 3 - 1);
			birth_r = x_1d_array[i];
			if (i < birth_spline_nodes) {
				y_1d_array[i] = birth_kernel_expression.value();
			}
			else {
				y_1d_array[i] = 0;
			}
		}
		spline1dinterpolant birth_kernel_spline;
		spline1dbuildmonotone(x_1d_array, y_1d_array, birth_kernel_spline);

		real_1d_array x_quantile_1d_array;
		real_1d_array y_quantile_1d_array;

		x_quantile_1d_array.setlength(birth_spline_nodes);
		y_quantile_1d_array.setlength(birth_spline_nodes);

		double approx_const = spline1dintegrate(birth_kernel_spline, birth_cutoff);
		for (int i = 0; i < birth_spline_nodes; i++) {
			x_quantile_1d_array[i] = (double)i / (birth_spline_nodes - 1);
			y_quantile_1d_array[i] =
				math::tools::newton_raphson_iterate(
					[=](double y) {return make_tuple(
						spline1dintegrate(birth_kernel_spline, y) / approx_const - x_quantile_1d_array[i],
						spline1dcalc(birth_kernel_spline, y) / approx_const); },
					1e-10, 0.0, birth_cutoff, numeric_limits<double>::digits);
		}

		spline1dinterpolant birth_reverse_cdf_spline;
		spline1dbuildmonotone(x_quantile_1d_array, y_quantile_1d_array, birth_reverse_cdf_spline);
		return birth_reverse_cdf_spline;
	}

}
#endif 