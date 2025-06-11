//
// Created by Artem Novikov on 25.05.2023.
//

#include <CLHEP/Units/PhysicalConstants.h>

#include "DefaultBuilder.h"

#include "../Fragments/StableFermiFragment.h"
#include "../Fragments/B9FermiFragment.h"
#include "../Fragments/Be8FermiFragment.h"
#include "../Fragments/He5FermiFragment.h"
#include "../Fragments/Li5FermiFragment.h"

using CLHEP::MeV;

namespace pool {

void DefaultBuilder::Build(FermiFragmentPool::Map& map) const {
  static const StableFermiFragment Fragment00(  1_m, 0_c,  2,  0.00*MeV ); AddFragment(map, Fragment00);
  static const StableFermiFragment Fragment01(  1_m, 1_c,  2,  0.00*MeV ); AddFragment(map, Fragment01);
  static const StableFermiFragment Fragment02(  2_m, 1_c,  3,  0.00*MeV ); AddFragment(map, Fragment02);
  static const StableFermiFragment Fragment03(  3_m, 1_c,  2,  0.00*MeV ); AddFragment(map, Fragment03);
  static const StableFermiFragment Fragment04(  3_m, 2_c,  2,  0.00*MeV ); AddFragment(map, Fragment04);
  static const StableFermiFragment Fragment05(  4_m, 2_c,  1,  0.00*MeV ); AddFragment(map, Fragment05);
  static const He5FermiFragment    Fragment06(  5_m, 2_c,  4, 16.76*MeV ); AddFragment(map, Fragment06);  /// He5
  static const Li5FermiFragment    Fragment07(  5_m, 3_c,  4, 16.66*MeV ); AddFragment(map, Fragment07);  /// Li5
  static const StableFermiFragment Fragment08(  6_m, 2_c,  1,  0.00*MeV ); AddFragment(map, Fragment08);
  static const StableFermiFragment Fragment09(  6_m, 3_c,  3,  0.00*MeV ); AddFragment(map, Fragment09);

  static const StableFermiFragment Fragment10(  6_m, 3_c,  1,  3.56*MeV ); AddFragment(map, Fragment10);
  static const StableFermiFragment Fragment11(  7_m, 3_c,  4,  0.00*MeV ); AddFragment(map, Fragment11);
  static const StableFermiFragment Fragment12(  7_m, 3_c,  2,  0.48*MeV ); AddFragment(map, Fragment12);
  static const StableFermiFragment Fragment13(  7_m, 4_c,  4,  0.00*MeV ); AddFragment(map, Fragment13);
  static const StableFermiFragment Fragment14(  7_m, 4_c,  2,  0.43*MeV ); AddFragment(map, Fragment14);
  static const StableFermiFragment Fragment15(  8_m, 3_c,  5,  0.00*MeV ); AddFragment(map, Fragment15);
  static const StableFermiFragment Fragment16(  8_m, 3_c,  3,  0.98*MeV ); AddFragment(map, Fragment16);
  static const Be8FermiFragment    Fragment17(  8_m, 4_c,  1,  0.00*MeV ); AddFragment(map, Fragment17);  /// Be8
  static const StableFermiFragment Fragment18(  9_m, 4_c,  4,  0.00*MeV ); AddFragment(map, Fragment18);
  static const B9FermiFragment     Fragment19(  9_m, 5_c,  4,  0.00*MeV ); AddFragment(map, Fragment19);  /// B9

  static const StableFermiFragment Fragment20( 10_m, 4_c,  1,  0.00*MeV ); AddFragment(map, Fragment20);
  static const StableFermiFragment Fragment21( 10_m, 4_c,  5,  3.37*MeV ); AddFragment(map, Fragment21);
  static const StableFermiFragment Fragment22( 10_m, 4_c,  8,  5.96*MeV ); AddFragment(map, Fragment22);
  static const StableFermiFragment Fragment23( 10_m, 4_c,  1,  6.18*MeV ); AddFragment(map, Fragment23);
  static const StableFermiFragment Fragment24( 10_m, 4_c,  5,  6.26*MeV ); AddFragment(map, Fragment24);
  static const StableFermiFragment Fragment25( 10_m, 5_c,  7,  0.00*MeV ); AddFragment(map, Fragment25);
  static const StableFermiFragment Fragment26( 10_m, 5_c,  3,  0.72*MeV ); AddFragment(map, Fragment26);
  static const StableFermiFragment Fragment27( 10_m, 5_c,  1,  1.74*MeV ); AddFragment(map, Fragment27);
  static const StableFermiFragment Fragment28( 10_m, 5_c,  3,  2.15*MeV ); AddFragment(map, Fragment28);
  static const StableFermiFragment Fragment29( 10_m, 5_c,  5,  3.59*MeV ); AddFragment(map, Fragment29);

  static const StableFermiFragment Fragment30( 10_m, 6_c,  3,  0.00*MeV ); AddFragment(map, Fragment30);
  static const StableFermiFragment Fragment31( 10_m, 6_c,  5,  3.35*MeV ); AddFragment(map, Fragment31);
  static const StableFermiFragment Fragment32( 11_m, 5_c,  4,  0.00*MeV ); AddFragment(map, Fragment32);
  static const StableFermiFragment Fragment33( 11_m, 5_c,  2,  2.13*MeV ); AddFragment(map, Fragment33);
  static const StableFermiFragment Fragment34( 11_m, 5_c,  6,  4.44*MeV ); AddFragment(map, Fragment34);
  static const StableFermiFragment Fragment35( 11_m, 5_c,  4,  5.02*MeV ); AddFragment(map, Fragment35);
  static const StableFermiFragment Fragment36( 11_m, 5_c, 10,  6.76*MeV ); AddFragment(map, Fragment36);
  static const StableFermiFragment Fragment37( 11_m, 5_c,  6,  7.29*MeV ); AddFragment(map, Fragment37);
  static const StableFermiFragment Fragment38( 11_m, 5_c,  4,  7.98*MeV ); AddFragment(map, Fragment38);
  static const StableFermiFragment Fragment39( 11_m, 5_c,  6,  8.56*MeV ); AddFragment(map, Fragment39);

  static const StableFermiFragment Fragment40( 11_m, 6_c,  4,  0.00*MeV ); AddFragment(map, Fragment40);
  static const StableFermiFragment Fragment41( 11_m, 6_c,  2,  2.00*MeV ); AddFragment(map, Fragment41);
  static const StableFermiFragment Fragment42( 11_m, 6_c,  6,  4.32*MeV ); AddFragment(map, Fragment42);
  static const StableFermiFragment Fragment43( 11_m, 6_c,  4,  4.80*MeV ); AddFragment(map, Fragment43);
  static const StableFermiFragment Fragment44( 11_m, 6_c,  2,  6.34*MeV ); AddFragment(map, Fragment44);
  static const StableFermiFragment Fragment45( 11_m, 6_c,  8,  6.48*MeV ); AddFragment(map, Fragment45);
  static const StableFermiFragment Fragment46( 11_m, 6_c,  6,  6.90*MeV ); AddFragment(map, Fragment46);
  static const StableFermiFragment Fragment47( 11_m, 6_c,  4,  7.50*MeV ); AddFragment(map, Fragment47);
  static const StableFermiFragment Fragment48( 11_m, 6_c,  4,  8.10*MeV ); AddFragment(map, Fragment48);
  static const StableFermiFragment Fragment49( 11_m, 6_c,  6,  8.42*MeV ); AddFragment(map, Fragment49);

  static const StableFermiFragment Fragment50( 11_m, 6_c,  8,  8.66*MeV ); AddFragment(map, Fragment50);
  static const StableFermiFragment Fragment51( 12_m, 5_c,  3,  0.00*MeV ); AddFragment(map, Fragment51);
  static const StableFermiFragment Fragment52( 12_m, 5_c,  5,  0.95*MeV ); AddFragment(map, Fragment52);
  static const StableFermiFragment Fragment53( 12_m, 5_c,  5,  1.67*MeV ); AddFragment(map, Fragment53);
  static const StableFermiFragment Fragment54( 12_m, 5_c,  4,  2.65*MeV ); AddFragment(map, Fragment54);
  static const StableFermiFragment Fragment55( 12_m, 6_c,  1,  0.00*MeV ); AddFragment(map, Fragment55);
  static const StableFermiFragment Fragment56( 12_m, 6_c,  5,  4.44*MeV ); AddFragment(map, Fragment56);
  static const StableFermiFragment Fragment57( 13_m, 6_c,  2,  0.00*MeV ); AddFragment(map, Fragment57);
  static const StableFermiFragment Fragment58( 13_m, 6_c,  2,  3.09*MeV ); AddFragment(map, Fragment58);
  static const StableFermiFragment Fragment59( 13_m, 6_c,  4,  3.68*MeV ); AddFragment(map, Fragment59);

  static const StableFermiFragment Fragment60( 13_m, 6_c,  6,  3.85*MeV ); AddFragment(map, Fragment60);
  static const StableFermiFragment Fragment61( 13_m, 7_c,  2,  0.00*MeV ); AddFragment(map, Fragment61);
  static const StableFermiFragment Fragment62( 14_m, 6_c,  1,  0.00*MeV ); AddFragment(map, Fragment62);
  static const StableFermiFragment Fragment63( 14_m, 6_c,  3,  6.09*MeV ); AddFragment(map, Fragment63);
  static const StableFermiFragment Fragment64( 14_m, 6_c,  8,  6.69*MeV ); AddFragment(map, Fragment64);
  static const StableFermiFragment Fragment65( 14_m, 6_c,  6,  6.96*MeV ); AddFragment(map, Fragment65);
  static const StableFermiFragment Fragment66( 14_m, 6_c,  5,  7.34*MeV ); AddFragment(map, Fragment66);
  static const StableFermiFragment Fragment67( 14_m, 7_c,  3,  0.00*MeV ); AddFragment(map, Fragment67);
  static const StableFermiFragment Fragment68( 14_m, 7_c,  1,  2.31*MeV ); AddFragment(map, Fragment68);
  static const StableFermiFragment Fragment69( 14_m, 7_c,  3,  3.95*MeV ); AddFragment(map, Fragment69);

  static const StableFermiFragment Fragment70( 14_m, 7_c,  1,  4.92*MeV ); AddFragment(map, Fragment70);
  static const StableFermiFragment Fragment71( 14_m, 7_c,  5,  5.11*MeV ); AddFragment(map, Fragment71);
  static const StableFermiFragment Fragment72( 14_m, 7_c,  3,  5.69*MeV ); AddFragment(map, Fragment72);
  static const StableFermiFragment Fragment73( 14_m, 7_c,  7,  5.83*MeV ); AddFragment(map, Fragment73);
  static const StableFermiFragment Fragment74( 14_m, 7_c,  3,  6.20*MeV ); AddFragment(map, Fragment74);
  static const StableFermiFragment Fragment75( 14_m, 7_c,  7,  6.44*MeV ); AddFragment(map, Fragment75);
  static const StableFermiFragment Fragment76( 14_m, 7_c,  5,  7.03*MeV ); AddFragment(map, Fragment76);
  static const StableFermiFragment Fragment77( 15_m, 7_c,  2,  0.00*MeV ); AddFragment(map, Fragment77);
  static const StableFermiFragment Fragment78( 15_m, 7_c,  8,  5.28*MeV ); AddFragment(map, Fragment78);
  static const StableFermiFragment Fragment79( 15_m, 7_c,  4,  6.32*MeV ); AddFragment(map, Fragment79);

  static const StableFermiFragment Fragment80( 15_m, 7_c, 10,  7.22*MeV ); AddFragment(map, Fragment80);
  static const StableFermiFragment Fragment81( 15_m, 7_c,  8,  7.57*MeV ); AddFragment(map, Fragment81);
  static const StableFermiFragment Fragment82( 15_m, 7_c,  2,  8.31*MeV ); AddFragment(map, Fragment82);
  static const StableFermiFragment Fragment83( 15_m, 7_c,  4,  8.57*MeV ); AddFragment(map, Fragment83);
  static const StableFermiFragment Fragment84( 15_m, 7_c, 14,  9.15*MeV ); AddFragment(map, Fragment84);
  static const StableFermiFragment Fragment85( 15_m, 7_c, 14,  9.79*MeV ); AddFragment(map, Fragment85);
  static const StableFermiFragment Fragment86( 15_m, 7_c,  8, 10.00*MeV ); AddFragment(map, Fragment86);
  static const StableFermiFragment Fragment87( 15_m, 8_c,  2,  0.00*MeV ); AddFragment(map, Fragment87);
  static const StableFermiFragment Fragment88( 15_m, 8_c,  8,  5.22*MeV ); AddFragment(map, Fragment88);
  static const StableFermiFragment Fragment89( 15_m, 8_c,  4,  6.18*MeV ); AddFragment(map, Fragment89);

  static const StableFermiFragment Fragment90( 15_m, 8_c, 10,  6.83*MeV ); AddFragment(map, Fragment90);
  static const StableFermiFragment Fragment91( 15_m, 8_c,  8,  7.28*MeV ); AddFragment(map, Fragment91);
  static const StableFermiFragment Fragment92( 16_m, 7_c,  5,  0.00*MeV ); AddFragment(map, Fragment92);
  static const StableFermiFragment Fragment93( 16_m, 7_c,  1,  0.12*MeV ); AddFragment(map, Fragment93);
  static const StableFermiFragment Fragment94( 16_m, 7_c,  7,  0.30*MeV ); AddFragment(map, Fragment94);
  static const StableFermiFragment Fragment95( 16_m, 7_c,  3,  0.40*MeV ); AddFragment(map, Fragment95);
  static const StableFermiFragment Fragment96( 16_m, 8_c,  1,  0.00*MeV ); AddFragment(map, Fragment96);
  static const StableFermiFragment Fragment97( 16_m, 8_c,  8,  6.10*MeV ); AddFragment(map, Fragment97);
  static const StableFermiFragment Fragment98( 16_m, 8_c,  5,  6.92*MeV ); AddFragment(map, Fragment98);
  static const StableFermiFragment Fragment99( 16_m, 8_c,  3,  7.12*MeV ); AddFragment(map, Fragment99);
}

void DefaultBuilder::AddFragment(FermiFragmentPool::Map& map, const FermiFragment& fragment) const {
  map.emplace(NucleiData{fragment.GetMassNumber(), fragment.GetChargeNumber()}, &fragment);
}

} // namespace pool
