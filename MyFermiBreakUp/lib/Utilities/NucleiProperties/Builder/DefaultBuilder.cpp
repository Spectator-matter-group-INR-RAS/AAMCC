//
// Created by Artem Novikov on 20.05.2023.
//

#include "DefaultBuilder.h"

namespace properties {

void DefaultBuilder::BuildTable(FermiNucleiProperties::MassMap& data) const {
  data.emplace(NucleiData{1_m, 0_c}, 939.565);
  data.emplace(NucleiData{1_m, 1_c}, 938.272);
  data.emplace(NucleiData{2_m, 1_c}, 1875.61);
  data.emplace(NucleiData{3_m, 1_c}, 2808.92);
  data.emplace(NucleiData{3_m, 2_c}, 2808.39);
  data.emplace(NucleiData{3_m, 3_c}, 2821.62);
  data.emplace(NucleiData{4_m, 1_c}, 3750.09);
  data.emplace(NucleiData{4_m, 2_c}, 3727.38);
  data.emplace(NucleiData{4_m, 3_c}, 3749.77);
  data.emplace(NucleiData{5_m, 1_c}, 4689.85);
  data.emplace(NucleiData{5_m, 2_c}, 4667.68);
  data.emplace(NucleiData{5_m, 3_c}, 4667.62);
  data.emplace(NucleiData{5_m, 4_c}, 4692.57);
  data.emplace(NucleiData{6_m, 1_c}, 5630.33);
  data.emplace(NucleiData{6_m, 2_c}, 5605.53);
  data.emplace(NucleiData{6_m, 3_c}, 5601.52);
  data.emplace(NucleiData{6_m, 4_c}, 5605.3);
  data.emplace(NucleiData{6_m, 5_c}, 5633.73);
  data.emplace(NucleiData{7_m, 1_c}, 6569.08);
  data.emplace(NucleiData{7_m, 2_c}, 6545.51);
  data.emplace(NucleiData{7_m, 3_c}, 6533.83);
  data.emplace(NucleiData{7_m, 4_c}, 6534.18);
  data.emplace(NucleiData{7_m, 5_c}, 6545.58);
  data.emplace(NucleiData{8_m, 2_c}, 7482.54);
  data.emplace(NucleiData{8_m, 3_c}, 7471.37);
  data.emplace(NucleiData{8_m, 4_c}, 7454.85);
  data.emplace(NucleiData{8_m, 5_c}, 7472.32);
  data.emplace(NucleiData{8_m, 6_c}, 7483.95);
  data.emplace(NucleiData{9_m, 2_c}, 8423.36);
  data.emplace(NucleiData{9_m, 3_c}, 8406.87);
  data.emplace(NucleiData{9_m, 4_c}, 8392.75);
  data.emplace(NucleiData{9_m, 5_c}, 8393.31);
  data.emplace(NucleiData{9_m, 6_c}, 8409.29);
  data.emplace(NucleiData{10_m, 2_c}, 9363.09);
  data.emplace(NucleiData{10_m, 3_c}, 9346.46);
  data.emplace(NucleiData{10_m, 4_c}, 9325.5);
  data.emplace(NucleiData{10_m, 5_c}, 9324.44);
  data.emplace(NucleiData{10_m, 6_c}, 9327.57);
  data.emplace(NucleiData{10_m, 7_c}, 9350.16);
  data.emplace(NucleiData{11_m, 3_c}, 10285.6);
  data.emplace(NucleiData{11_m, 4_c}, 10264.6);
  data.emplace(NucleiData{11_m, 5_c}, 10252.5);
  data.emplace(NucleiData{11_m, 6_c}, 10254.0);
  data.emplace(NucleiData{11_m, 7_c}, 10267.2);
  data.emplace(NucleiData{12_m, 3_c}, 11225.3);
  data.emplace(NucleiData{12_m, 4_c}, 11201.0);
  data.emplace(NucleiData{12_m, 5_c}, 11188.7);
  data.emplace(NucleiData{12_m, 6_c}, 11174.9);
  data.emplace(NucleiData{12_m, 7_c}, 11191.7);
  data.emplace(NucleiData{12_m, 8_c}, 11205.8);
  data.emplace(NucleiData{13_m, 3_c}, 12166.2);
  data.emplace(NucleiData{13_m, 4_c}, 12141.0);
  data.emplace(NucleiData{13_m, 5_c}, 12123.4);
  data.emplace(NucleiData{13_m, 6_c}, 12109.5);
  data.emplace(NucleiData{13_m, 7_c}, 12111.2);
  data.emplace(NucleiData{13_m, 8_c}, 12128.5);
  data.emplace(NucleiData{14_m, 4_c}, 13078.8);
  data.emplace(NucleiData{14_m, 5_c}, 13062.0);
  data.emplace(NucleiData{14_m, 6_c}, 13040.9);
  data.emplace(NucleiData{14_m, 7_c}, 13040.2);
  data.emplace(NucleiData{14_m, 8_c}, 13044.8);
  data.emplace(NucleiData{14_m, 9_c}, 13068.3);
  data.emplace(NucleiData{15_m, 4_c}, 14020.1);
  data.emplace(NucleiData{15_m, 5_c}, 13998.8);
  data.emplace(NucleiData{15_m, 6_c}, 13979.2);
  data.emplace(NucleiData{15_m, 7_c}, 13968.9);
  data.emplace(NucleiData{15_m, 8_c}, 13971.2);
  data.emplace(NucleiData{15_m, 9_c}, 13984.6);
  data.emplace(NucleiData{16_m, 4_c}, 14959.3);
  data.emplace(NucleiData{16_m, 5_c}, 14938.5);
  data.emplace(NucleiData{16_m, 6_c}, 14914.5);
  data.emplace(NucleiData{16_m, 7_c}, 14906.0);
  data.emplace(NucleiData{16_m, 8_c}, 14895.1);
  data.emplace(NucleiData{16_m, 9_c}, 14910.0);
  data.emplace(NucleiData{16_m, 10_c}, 14922.8);
  data.emplace(NucleiData{17_m, 5_c}, 15876.6);
  data.emplace(NucleiData{17_m, 6_c}, 15853.4);
  data.emplace(NucleiData{17_m, 7_c}, 15839.7);
  data.emplace(NucleiData{17_m, 8_c}, 15830.5);
  data.emplace(NucleiData{17_m, 9_c}, 15832.8);
  data.emplace(NucleiData{17_m, 10_c}, 15846.8);
  data.emplace(NucleiData{18_m, 5_c}, 16816.2);
  data.emplace(NucleiData{18_m, 6_c}, 16788.7);
  data.emplace(NucleiData{18_m, 7_c}, 16776.4);
  data.emplace(NucleiData{18_m, 8_c}, 16762.0);
  data.emplace(NucleiData{18_m, 9_c}, 16763.2);
  data.emplace(NucleiData{18_m, 10_c}, 16767.1);
  data.emplace(NucleiData{18_m, 11_c}, 16786.3);
  data.emplace(NucleiData{19_m, 5_c}, 17754.6);
  data.emplace(NucleiData{19_m, 6_c}, 17727.7);
  data.emplace(NucleiData{19_m, 7_c}, 17710.7);
  data.emplace(NucleiData{19_m, 8_c}, 17697.6);
  data.emplace(NucleiData{19_m, 9_c}, 17692.3);
  data.emplace(NucleiData{19_m, 10_c}, 17695.0);
  data.emplace(NucleiData{19_m, 11_c}, 17705.7);
  data.emplace(NucleiData{19_m, 12_c}, 17724.1);
  data.emplace(NucleiData{20_m, 5_c}, 18694.5);
  data.emplace(NucleiData{20_m, 6_c}, 18664.4);
  data.emplace(NucleiData{20_m, 7_c}, 18648.1);
  data.emplace(NucleiData{20_m, 8_c}, 18629.6);
  data.emplace(NucleiData{20_m, 9_c}, 18625.3);
  data.emplace(NucleiData{20_m, 10_c}, 18617.7);
  data.emplace(NucleiData{20_m, 11_c}, 18631.1);
  data.emplace(NucleiData{20_m, 12_c}, 18641.3);
}

} // namespace properties
