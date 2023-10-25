/*
 * This file is part of Multem.
 * Copyright 2022 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * Multem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version of the License, or
 * (at your option) any later version.
 *
 * Multem is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Multem. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef PN_FACT_H
	#define PN_FACT_H

	#include <algorithm>
	#include <vector>

	#include "macros.h"
	#include "const_enum.h"

	namespace mt
	{
		/***************************************************************************************/
		/******************** factorization of a number in small prime number ******************/
		// factorization: 2^a × 3^b × 5^c × 7^d, where a>1 
		// https:// docs.nvidia.com/cuda/cufft/index.html

		struct PN_Fact
		{
			public:
				PN_Fact()
				{
					load_prime_numbers();
				}

				dt_int32 operator()(dt_int64 n, eDat_Sel_Typ dst = edst_closest)
				{
					auto p_idx = std::min_element(pnv.begin(), pnv.end(), [n](dt_int64 p0, dt_int64 pe){ return std::abs(n-p0) < std::abs(n-pe); });

					auto pn = static_cast<dt_int32>(*p_idx);

					switch(dst)
					{
						case edst_less_than:
						{
							if (pn>=n)
							{
								pn = static_cast<dt_int32>(*(p_idx-1));
							}
						}
						break;
						case edst_eless_than:
						{
							if (pn>n)
							{
								pn = static_cast<dt_int32>(*(p_idx-1));
							}
						}
						break;
						case edst_greater_than:
						{
							if (pn<=n)
							{
								pn = static_cast<dt_int32>(*(p_idx+1));
							}
						}
						break;
						case edst_egreater_than:
						{
							if (pn<n)
							{
								pn = static_cast<dt_int32>(*(p_idx+1));
							}
						}
						break;
					}
					return pn;
				}

			private:
				std::vector<dt_int64> pnv;

				void load_prime_numbers()
				{
					dt_int64 b_2 = 2;
					dt_int64 b_3 = 3;
					dt_int64 b_5 = 5;
					dt_int64 b_7 = 7;

					dt_int32 e_2_m = 16;
					dt_int32 e_3_m = 7;
					dt_int32 e_5_m = 5;
					dt_int32 e_7_m = 4;

					dt_int64 pn_0 = std::pow(b_2, 5);
					dt_int64 pn_e = static_cast<dt_int64>(std::pow(b_2, e_2_m));

					// reserve enough space
					pnv.reserve((e_2_m+1)*(e_3_m + 1)*(e_5_m + 1)*(e_7_m + 1));

					// fill in first numbers
					for(auto ie=1; ie<6; ie++)
					{
						auto p = static_cast<dt_int64>(std::pow(b_2, ie));
						pnv.push_back(p);
					}

					// calculate the other numbers
					for(auto ip7=0; ip7<=e_7_m; ip7++)
					{
						auto p7 = static_cast<dt_int64>(std::pow(b_7, ip7));

						for(auto ip5=0; ip5<=e_5_m; ip5++)
						{
							auto p5 = static_cast<dt_int64>(std::pow(b_5, ip5));

							for(auto ip3=0; ip3<=e_3_m; ip3++)
							{
								auto p3 = static_cast<dt_int64>(std::pow(b_3, ip3));

								for(auto ip2=1; ip2<=e_2_m; ip2++)
								{
									auto p2 = static_cast<dt_int64>(std::pow(b_2, ip2));

									auto p = p7*p5*p3*p2;

									if ((pn_0<p) && (p<=pn_e))
									{
										pnv.push_back(p);
									}
								}
							}
						}
					}
					std::sort(pnv.begin(), pnv.end());
					pnv.erase(std::unique(pnv.begin(), pnv.end()), pnv.end());
					pnv.shrink_to_fit();
				}
		};
	}

#endif