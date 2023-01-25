/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <LibintIntegrals/BasisSetHandler.h>
#include <LibintIntegrals/LibintIntegrals.h>
#include <Utils/Constants.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Settings.h>
#include <gmock/gmock.h>
#include <ctime>

using namespace Scine;
using namespace Integrals;

using namespace testing;

class TwoBodyIntsTest : public Test {};

TEST_F(TwoBodyIntsTest, Test2body) {
  std::stringstream h2("2\n\n"
                       "H 0 0 0\n"
                       "H 1.2 0 0");
  auto scineAtoms = Utils::XyzStreamHandler::read(h2);

  using IndexType = std::array<int, 4>;
  std::map<IndexType, double> pyScfEriMap;

  pyScfEriMap.insert({{0, 0, 0, 0}, 0.8971413572489109});
  pyScfEriMap.insert({{1, 0, 0, 0}, 0.5045902661946233});
  pyScfEriMap.insert({{1, 0, 1, 0}, 0.2999243474108791});
  pyScfEriMap.insert({{1, 1, 0, 0}, 0.507267480114002});
  pyScfEriMap.insert({{1, 1, 1, 0}, 0.3239341914533078});
  pyScfEriMap.insert({{1, 1, 1, 1}, 0.3940445135380343});
  pyScfEriMap.insert({{2, 0, 2, 0}, 0.16045646128458113});
  pyScfEriMap.insert({{2, 1, 2, 0}, 0.10166761452628235});
  pyScfEriMap.insert({{2, 1, 2, 1}, 0.06893333251070229});
  pyScfEriMap.insert({{2, 2, 0, 0}, 0.799671506518246});
  pyScfEriMap.insert({{2, 2, 1, 0}, 0.47249239596550685});
  pyScfEriMap.insert({{2, 2, 1, 1}, 0.4962129991093084});
  pyScfEriMap.insert({{2, 2, 2, 2}, 0.8242232905265838});
  pyScfEriMap.insert({{3, 0, 3, 0}, 0.16045646128458113});
  pyScfEriMap.insert({{3, 1, 3, 0}, 0.10166761452628235});
  pyScfEriMap.insert({{3, 1, 3, 1}, 0.06893333251070229});
  pyScfEriMap.insert({{3, 2, 3, 2}, 0.05046265044040281});
  pyScfEriMap.insert({{3, 3, 0, 0}, 0.799671506518246});
  pyScfEriMap.insert({{3, 3, 1, 0}, 0.47249239596550685});
  pyScfEriMap.insert({{3, 3, 1, 1}, 0.4962129991093084});
  pyScfEriMap.insert({{3, 3, 2, 2}, 0.7232979896457781});
  pyScfEriMap.insert({{3, 3, 3, 3}, 0.8242232905265838});
  pyScfEriMap.insert({{4, 0, 4, 0}, 0.16045646128458113});
  pyScfEriMap.insert({{4, 1, 4, 0}, 0.10166761452628235});
  pyScfEriMap.insert({{4, 1, 4, 1}, 0.06893333251070229});
  pyScfEriMap.insert({{4, 2, 4, 2}, 0.050462650440402795});
  pyScfEriMap.insert({{4, 3, 4, 3}, 0.050462650440402795});
  pyScfEriMap.insert({{4, 4, 0, 0}, 0.799671506518246});
  pyScfEriMap.insert({{4, 4, 1, 0}, 0.47249239596550685});
  pyScfEriMap.insert({{4, 4, 1, 1}, 0.4962129991093084});
  pyScfEriMap.insert({{4, 4, 2, 2}, 0.7232979896457781});
  pyScfEriMap.insert({{4, 4, 3, 3}, 0.7232979896457781});
  pyScfEriMap.insert({{4, 4, 4, 4}, 0.8242232905265837});
  pyScfEriMap.insert({{5, 0, 0, 0}, 0.17116648346941746});
  pyScfEriMap.insert({{5, 0, 1, 0}, 0.10421861103427668});
  pyScfEriMap.insert({{5, 0, 1, 1}, 0.11624632277323764});
  pyScfEriMap.insert({{5, 0, 2, 0}, 0.03550331171153253});
  pyScfEriMap.insert({{5, 0, 2, 1}, 0.0256185662050172});
  pyScfEriMap.insert({{5, 0, 2, 2}, 0.17288209217843925});
  pyScfEriMap.insert({{5, 0, 3, 3}, 0.1583823156007121});
  pyScfEriMap.insert({{5, 0, 4, 4}, 0.1583823156007121});
  pyScfEriMap.insert({{5, 0, 5, 0}, 0.0498204840433603});
  pyScfEriMap.insert({{5, 1, 0, 0}, 0.2032715434118313});
  pyScfEriMap.insert({{5, 1, 1, 0}, 0.1317237473992584});
  pyScfEriMap.insert({{5, 1, 1, 1}, 0.16319488597198312});
  pyScfEriMap.insert({{5, 1, 2, 0}, 0.04348506376975387});
  pyScfEriMap.insert({{5, 1, 2, 1}, 0.034455293464668156});
  pyScfEriMap.insert({{5, 1, 2, 2}, 0.21306349915127773});
  pyScfEriMap.insert({{5, 1, 3, 3}, 0.19430030777437562});
  pyScfEriMap.insert({{5, 1, 4, 4}, 0.19430030777437562});
  pyScfEriMap.insert({{5, 1, 5, 0}, 0.06852119159779566});
  pyScfEriMap.insert({{5, 1, 5, 1}, 0.10859350105043604});
  pyScfEriMap.insert({{5, 2, 0, 0}, 0.18584140868314752});
  pyScfEriMap.insert({{5, 2, 1, 0}, 0.11590150288979031});
  pyScfEriMap.insert({{5, 2, 1, 1}, 0.1316339018760785});
  pyScfEriMap.insert({{5, 2, 2, 0}, 0.05650626394988069});
  pyScfEriMap.insert({{5, 2, 2, 1}, 0.0404488140893028});
  pyScfEriMap.insert({{5, 2, 2, 2}, 0.20173315659836946});
  pyScfEriMap.insert({{5, 2, 3, 3}, 0.17197648182006686});
  pyScfEriMap.insert({{5, 2, 4, 4}, 0.17197648182006686});
  pyScfEriMap.insert({{5, 2, 5, 0}, 0.06232627778983818});
  pyScfEriMap.insert({{5, 2, 5, 1}, 0.08604708851796407});
  pyScfEriMap.insert({{5, 2, 5, 2}, 0.08172625576600219});
  pyScfEriMap.insert({{5, 3, 3, 0}, 0.02417719512688181});
  pyScfEriMap.insert({{5, 3, 3, 1}, 0.01646015169183101});
  pyScfEriMap.insert({{5, 3, 3, 2}, 0.008863018633008181});
  pyScfEriMap.insert({{5, 3, 5, 3}, 0.006287425808184565});
  pyScfEriMap.insert({{5, 4, 4, 0}, 0.02417719512688181});
  pyScfEriMap.insert({{5, 4, 4, 1}, 0.01646015169183101});
  pyScfEriMap.insert({{5, 4, 4, 2}, 0.008863018633008181});
  pyScfEriMap.insert({{5, 4, 5, 4}, 0.006287425808184565});
  pyScfEriMap.insert({{5, 5, 0, 0}, 0.4329465035895467});
  pyScfEriMap.insert({{5, 5, 1, 0}, 0.2886882511459753});
  pyScfEriMap.insert({{5, 5, 1, 1}, 0.37483587156644715});
  pyScfEriMap.insert({{5, 5, 2, 0}, 0.10723601876498609});
  pyScfEriMap.insert({{5, 5, 2, 1}, 0.08888209252263984});
  pyScfEriMap.insert({{5, 5, 2, 2}, 0.4733419945343873});
  pyScfEriMap.insert({{5, 5, 3, 3}, 0.4117361160007401});
  pyScfEriMap.insert({{5, 5, 4, 4}, 0.4117361160007401});
  pyScfEriMap.insert({{5, 5, 5, 0}, 0.17116648346941749});
  pyScfEriMap.insert({{5, 5, 5, 1}, 0.2957732012070589});
  pyScfEriMap.insert({{5, 5, 5, 2}, 0.21932521319374862});
  pyScfEriMap.insert({{5, 5, 5, 5}, 0.8971413572489109});
  pyScfEriMap.insert({{6, 0, 0, 0}, 0.29577320120705897});
  pyScfEriMap.insert({{6, 0, 1, 0}, 0.1770461343613643});
  pyScfEriMap.insert({{6, 0, 1, 1}, 0.19330002183834524});
  pyScfEriMap.insert({{6, 0, 2, 0}, 0.021841005637725123});
  pyScfEriMap.insert({{6, 0, 2, 1}, 0.015724532021340134});
  pyScfEriMap.insert({{6, 0, 2, 2}, 0.28015421843182475});
  pyScfEriMap.insert({{6, 0, 3, 3}, 0.27736521609567444});
  pyScfEriMap.insert({{6, 0, 4, 4}, 0.27736521609567444});
  pyScfEriMap.insert({{6, 0, 5, 0}, 0.06852119159779567});
  pyScfEriMap.insert({{6, 0, 5, 1}, 0.08931051119967236});
  pyScfEriMap.insert({{6, 0, 5, 2}, 0.07940304268397011});
  pyScfEriMap.insert({{6, 0, 5, 5}, 0.2032715434118313});
  pyScfEriMap.insert({{6, 0, 6, 0}, 0.10859350105043607});
  pyScfEriMap.insert({{6, 1, 0, 0}, 0.3408731587587354});
  pyScfEriMap.insert({{6, 1, 1, 0}, 0.21998609979440986});
  pyScfEriMap.insert({{6, 1, 1, 1}, 0.27361818189105797});
  pyScfEriMap.insert({{6, 1, 2, 0}, 0.030684222285006415});
  pyScfEriMap.insert({{6, 1, 2, 1}, 0.024673898670214037});
  pyScfEriMap.insert({{6, 1, 2, 2}, 0.33787120218693417});
  pyScfEriMap.insert({{6, 1, 3, 3}, 0.33336674981076697});
  pyScfEriMap.insert({{6, 1, 4, 4}, 0.33336674981076697});
  pyScfEriMap.insert({{6, 1, 5, 0}, 0.09194410995484553});
  pyScfEriMap.insert({{6, 1, 5, 1}, 0.13889611442500882});
  pyScfEriMap.insert({{6, 1, 5, 2}, 0.10843837784441013});
  pyScfEriMap.insert({{6, 1, 5, 5}, 0.3408731587587355});
  pyScfEriMap.insert({{6, 1, 6, 0}, 0.13889611442500885});
  pyScfEriMap.insert({{6, 1, 6, 1}, 0.21047319041990956});
  pyScfEriMap.insert({{6, 2, 0, 0}, 0.11280260876047243});
  pyScfEriMap.insert({{6, 2, 1, 0}, 0.07021769029992721});
  pyScfEriMap.insert({{6, 2, 1, 1}, 0.07963155244133215});
  pyScfEriMap.insert({{6, 2, 2, 0}, 0.06335655647293138});
  pyScfEriMap.insert({{6, 2, 2, 1}, 0.04356423204054653});
  pyScfEriMap.insert({{6, 2, 2, 2}, 0.11957743164287343});
  pyScfEriMap.insert({{6, 2, 3, 3}, 0.1055758621712681});
  pyScfEriMap.insert({{6, 2, 4, 4}, 0.10557586217126812});
  pyScfEriMap.insert({{6, 2, 5, 0}, 0.04331207102082376});
  pyScfEriMap.insert({{6, 2, 5, 1}, 0.058558253202364643});
  pyScfEriMap.insert({{6, 2, 5, 2}, 0.05790799586710096});
  pyScfEriMap.insert({{6, 2, 5, 5}, 0.14387853945185763});
  pyScfEriMap.insert({{6, 2, 6, 0}, 0.05213660701482233});
  pyScfEriMap.insert({{6, 2, 6, 1}, 0.07145738209918806});
  pyScfEriMap.insert({{6, 2, 6, 2}, 0.046324428069679104});
  pyScfEriMap.insert({{6, 3, 3, 0}, 0.057265343623830206});
  pyScfEriMap.insert({{6, 3, 3, 1}, 0.03902326337821374});
  pyScfEriMap.insert({{6, 3, 3, 2}, 0.006696202016495291});
  pyScfEriMap.insert({{6, 3, 5, 3}, 0.010835174709857428});
  pyScfEriMap.insert({{6, 3, 6, 3}, 0.023215675300076073});
  pyScfEriMap.insert({{6, 4, 4, 0}, 0.057265343623830206});
  pyScfEriMap.insert({{6, 4, 4, 1}, 0.03902326337821374});
  pyScfEriMap.insert({{6, 4, 4, 2}, 0.006696202016495291});
  pyScfEriMap.insert({{6, 4, 5, 4}, 0.010835174709857428});
  pyScfEriMap.insert({{6, 4, 6, 4}, 0.023215675300076073});
  pyScfEriMap.insert({{6, 5, 0, 0}, 0.2886882511459753});
  pyScfEriMap.insert({{6, 5, 1, 0}, 0.19121058701903482});
  pyScfEriMap.insert({{6, 5, 1, 1}, 0.24778172763834222});
  pyScfEriMap.insert({{6, 5, 2, 0}, 0.06319523145777861});
  pyScfEriMap.insert({{6, 5, 2, 1}, 0.05177114014655456});
  pyScfEriMap.insert({{6, 5, 2, 2}, 0.30693303357192364});
  pyScfEriMap.insert({{6, 5, 3, 3}, 0.2771783914805481});
  pyScfEriMap.insert({{6, 5, 4, 4}, 0.2771783914805481});
  pyScfEriMap.insert({{6, 5, 5, 0}, 0.10421861103427671});
  pyScfEriMap.insert({{6, 5, 5, 1}, 0.17704613436136424});
  pyScfEriMap.insert({{6, 5, 5, 2}, 0.13080009259229597});
  pyScfEriMap.insert({{6, 5, 5, 5}, 0.5045902661946233});
  pyScfEriMap.insert({{6, 5, 6, 0}, 0.1317237473992584});
  pyScfEriMap.insert({{6, 5, 6, 1}, 0.2199860997944098});
  pyScfEriMap.insert({{6, 5, 6, 2}, 0.0878862239014157});
  pyScfEriMap.insert({{6, 5, 6, 5}, 0.2999243474108791});
  pyScfEriMap.insert({{6, 6, 0, 0}, 0.3748358715664472});
  pyScfEriMap.insert({{6, 6, 1, 0}, 0.2477817276383423});
  pyScfEriMap.insert({{6, 6, 1, 1}, 0.325113388608948});
  pyScfEriMap.insert({{6, 6, 2, 0}, 0.054061154647772794});
  pyScfEriMap.insert({{6, 6, 2, 1}, 0.044635306462380404});
  pyScfEriMap.insert({{6, 6, 2, 2}, 0.3809038930508819});
  pyScfEriMap.insert({{6, 6, 3, 3}, 0.366770807444573});
  pyScfEriMap.insert({{6, 6, 4, 4}, 0.366770807444573});
  pyScfEriMap.insert({{6, 6, 5, 0}, 0.11624632277323765});
  pyScfEriMap.insert({{6, 6, 5, 1}, 0.19330002183834524});
  pyScfEriMap.insert({{6, 6, 5, 2}, 0.13995832055233776});
  pyScfEriMap.insert({{6, 6, 5, 5}, 0.507267480114002});
  pyScfEriMap.insert({{6, 6, 6, 0}, 0.16319488597198314});
  pyScfEriMap.insert({{6, 6, 6, 1}, 0.27361818189105797});
  pyScfEriMap.insert({{6, 6, 6, 2}, 0.09434662780029085});
  pyScfEriMap.insert({{6, 6, 6, 5}, 0.3239341914533078});
  pyScfEriMap.insert({{6, 6, 6, 6}, 0.3940445135380343});
  pyScfEriMap.insert({{7, 0, 0, 0}, -0.2193252131937487});
  pyScfEriMap.insert({{7, 0, 1, 0}, -0.130800092592296});
  pyScfEriMap.insert({{7, 0, 1, 1}, -0.1399583205523378});
  pyScfEriMap.insert({{7, 0, 2, 0}, -0.05333839291748142});
  pyScfEriMap.insert({{7, 0, 2, 1}, -0.03640432101031121});
  pyScfEriMap.insert({{7, 0, 2, 2}, -0.22369592285605533});
  pyScfEriMap.insert({{7, 0, 3, 3}, -0.19942022994504388});
  pyScfEriMap.insert({{7, 0, 4, 4}, -0.1994202299450439});
  pyScfEriMap.insert({{7, 0, 5, 0}, -0.06232627778983822});
  pyScfEriMap.insert({{7, 0, 5, 1}, -0.07940304268397018});
  pyScfEriMap.insert({{7, 0, 5, 2}, -0.07904627488155397});
  pyScfEriMap.insert({{7, 0, 5, 5}, -0.18584140868314758});
  pyScfEriMap.insert({{7, 0, 6, 0}, -0.08604708851796412});
  pyScfEriMap.insert({{7, 0, 6, 1}, -0.10843837784441014});
  pyScfEriMap.insert({{7, 0, 6, 2}, -0.056543106261907866});
  pyScfEriMap.insert({{7, 0, 6, 5}, -0.11590150288979038});
  pyScfEriMap.insert({{7, 0, 6, 6}, -0.13163390187607854});
  pyScfEriMap.insert({{7, 0, 7, 0}, 0.08172625576600223});
  pyScfEriMap.insert({{7, 1, 0, 0}, -0.1438785394518576});
  pyScfEriMap.insert({{7, 1, 1, 0}, -0.08788622390141568});
  pyScfEriMap.insert({{7, 1, 1, 1}, -0.09434662780029081});
  pyScfEriMap.insert({{7, 1, 2, 0}, -0.03839865406705271});
  pyScfEriMap.insert({{7, 1, 2, 1}, -0.027256477556459782});
  pyScfEriMap.insert({{7, 1, 2, 2}, -0.15188539030652465});
  pyScfEriMap.insert({{7, 1, 3, 3}, -0.13349755191025142});
  pyScfEriMap.insert({{7, 1, 4, 4}, -0.13349755191025142});
  pyScfEriMap.insert({{7, 1, 5, 0}, -0.04331207102082377});
  pyScfEriMap.insert({{7, 1, 5, 1}, -0.05213660701482234});
  pyScfEriMap.insert({{7, 1, 5, 2}, -0.05654310626190785});
  pyScfEriMap.insert({{7, 1, 5, 5}, -0.11280260876047236});
  pyScfEriMap.insert({{7, 1, 6, 0}, -0.058558253202364643});
  pyScfEriMap.insert({{7, 1, 6, 1}, -0.07145738209918802});
  pyScfEriMap.insert({{7, 1, 6, 2}, -0.04061198931369482});
  pyScfEriMap.insert({{7, 1, 6, 5}, -0.07021769029992721});
  pyScfEriMap.insert({{7, 1, 6, 6}, -0.07963155244133216});
  pyScfEriMap.insert({{7, 1, 7, 0}, 0.057907995867100964});
  pyScfEriMap.insert({{7, 1, 7, 1}, 0.0463244280696791});
  pyScfEriMap.insert({{7, 2, 0, 0}, -0.2853531269735247});
  pyScfEriMap.insert({{7, 2, 1, 0}, -0.17437822753147048});
  pyScfEriMap.insert({{7, 2, 1, 1}, -0.19104846580194246});
  pyScfEriMap.insert({{7, 2, 2, 0}, -0.08696646782687228});
  pyScfEriMap.insert({{7, 2, 2, 1}, -0.06017767155985875});
  pyScfEriMap.insert({{7, 2, 2, 2}, -0.30991963031425435});
  pyScfEriMap.insert({{7, 2, 3, 3}, -0.26055742543274163});
  pyScfEriMap.insert({{7, 2, 4, 4}, -0.26055742543274163});
  pyScfEriMap.insert({{7, 2, 5, 0}, -0.09074725632834978});
  pyScfEriMap.insert({{7, 2, 5, 1}, -0.11811848707816276});
  pyScfEriMap.insert({{7, 2, 5, 2}, -0.11930542817766498});
  pyScfEriMap.insert({{7, 2, 5, 5}, -0.2853531269735247});
  pyScfEriMap.insert({{7, 2, 6, 0}, -0.11811848707816276});
  pyScfEriMap.insert({{7, 2, 6, 1}, -0.15316832618151804});
  pyScfEriMap.insert({{7, 2, 6, 2}, -0.08561630600246517});
  pyScfEriMap.insert({{7, 2, 6, 5}, -0.1743782275314705});
  pyScfEriMap.insert({{7, 2, 6, 6}, -0.19104846580194246});
  pyScfEriMap.insert({{7, 2, 7, 0}, 0.11930542817766507});
  pyScfEriMap.insert({{7, 2, 7, 1}, 0.08561630600246517});
  pyScfEriMap.insert({{7, 2, 7, 2}, 0.1795410770511904});
  pyScfEriMap.insert({{7, 3, 3, 0}, -0.030697482625577373});
  pyScfEriMap.insert({{7, 3, 3, 1}, -0.020574279542097167});
  pyScfEriMap.insert({{7, 3, 3, 2}, -0.014217679958253654});
  pyScfEriMap.insert({{7, 3, 5, 3}, -0.008579518074314089});
  pyScfEriMap.insert({{7, 3, 6, 3}, -0.013924675467440325});
  pyScfEriMap.insert({{7, 3, 7, 3}, 0.012134961853727665});
  pyScfEriMap.insert({{7, 4, 4, 0}, -0.030697482625577373});
  pyScfEriMap.insert({{7, 4, 4, 1}, -0.020574279542097167});
  pyScfEriMap.insert({{7, 4, 4, 2}, -0.014217679958253654});
  pyScfEriMap.insert({{7, 4, 5, 4}, -0.008579518074314089});
  pyScfEriMap.insert({{7, 4, 6, 4}, -0.013924675467440325});
  pyScfEriMap.insert({{7, 4, 7, 4}, 0.012134961853727663});
  pyScfEriMap.insert({{7, 5, 0, 0}, -0.10723601876498609});
  pyScfEriMap.insert({{7, 5, 1, 0}, -0.06319523145777861});
  pyScfEriMap.insert({{7, 5, 1, 1}, -0.0540611546477728});
  pyScfEriMap.insert({{7, 5, 2, 0}, -0.044463638043451534});
  pyScfEriMap.insert({{7, 5, 2, 1}, -0.031157457925883862});
  pyScfEriMap.insert({{7, 5, 2, 2}, -0.12682174726047626});
  pyScfEriMap.insert({{7, 5, 3, 3}, -0.09474336352639352});
  pyScfEriMap.insert({{7, 5, 4, 4}, -0.09474336352639352});
  pyScfEriMap.insert({{7, 5, 5, 0}, -0.03550331171153254});
  pyScfEriMap.insert({{7, 5, 5, 1}, -0.02184100563772512});
  pyScfEriMap.insert({{7, 5, 5, 2}, -0.05333839291748143});
  pyScfEriMap.insert({{7, 5, 5, 5}, -5.719667761977114e-17});
  pyScfEriMap.insert({{7, 5, 6, 0}, -0.04348506376975387});
  pyScfEriMap.insert({{7, 5, 6, 1}, -0.03068422228500642});
  pyScfEriMap.insert({{7, 5, 6, 2}, -0.0383986540670527});
  pyScfEriMap.insert({{7, 5, 6, 5}, -2.4863871911576338e-17});
  pyScfEriMap.insert({{7, 5, 6, 6}, -4.665579778592396e-18});
  pyScfEriMap.insert({{7, 5, 7, 0}, 0.05650626394988069});
  pyScfEriMap.insert({{7, 5, 7, 1}, 0.06335655647293137});
  pyScfEriMap.insert({{7, 5, 7, 2}, 0.08696646782687226});
  pyScfEriMap.insert({{7, 5, 7, 5}, 0.16045646128458113});
  pyScfEriMap.insert({{7, 6, 0, 0}, -0.08888209252263983});
  pyScfEriMap.insert({{7, 6, 1, 0}, -0.05177114014655457});
  pyScfEriMap.insert({{7, 6, 1, 1}, -0.044635306462380404});
  pyScfEriMap.insert({{7, 6, 2, 0}, -0.03115745792588387});
  pyScfEriMap.insert({{7, 6, 2, 1}, -0.021263040256737385});
  pyScfEriMap.insert({{7, 6, 2, 2}, -0.09786582617251932});
  pyScfEriMap.insert({{7, 6, 3, 3}, -0.08022143347148329});
  pyScfEriMap.insert({{7, 6, 4, 4}, -0.0802214334714833});
  pyScfEriMap.insert({{7, 6, 5, 0}, -0.0256185662050172});
  pyScfEriMap.insert({{7, 6, 5, 1}, -0.015724532021340134});
  pyScfEriMap.insert({{7, 6, 5, 2}, -0.03640432101031119});
  pyScfEriMap.insert({{7, 6, 5, 5}, -2.6971877010665357e-17});
  pyScfEriMap.insert({{7, 6, 6, 0}, -0.034455293464668156});
  pyScfEriMap.insert({{7, 6, 6, 1}, -0.024673898670214037});
  pyScfEriMap.insert({{7, 6, 6, 2}, -0.027256477556459782});
  pyScfEriMap.insert({{7, 6, 6, 5}, -1.581349013992337e-17});
  pyScfEriMap.insert({{7, 6, 7, 0}, 0.0404488140893028});
  pyScfEriMap.insert({{7, 6, 7, 1}, 0.04356423204054653});
  pyScfEriMap.insert({{7, 6, 7, 2}, 0.06017767155985876});
  pyScfEriMap.insert({{7, 6, 7, 5}, 0.10166761452628235});
  pyScfEriMap.insert({{7, 6, 7, 6}, 0.06893333251070229});
  pyScfEriMap.insert({{7, 7, 0, 0}, 0.47334199453438736});
  pyScfEriMap.insert({{7, 7, 1, 0}, 0.3069330335719237});
  pyScfEriMap.insert({{7, 7, 1, 1}, 0.3809038930508819});
  pyScfEriMap.insert({{7, 7, 2, 0}, 0.12682174726047626});
  pyScfEriMap.insert({{7, 7, 2, 1}, 0.09786582617251932});
  pyScfEriMap.insert({{7, 7, 2, 2}, 0.520805378744601});
  pyScfEriMap.insert({{7, 7, 3, 3}, 0.44399727541714773});
  pyScfEriMap.insert({{7, 7, 4, 4}, 0.44399727541714773});
  pyScfEriMap.insert({{7, 7, 5, 0}, 0.1728820921784393});
  pyScfEriMap.insert({{7, 7, 5, 1}, 0.28015421843182475});
  pyScfEriMap.insert({{7, 7, 5, 2}, 0.2236959228560552});
  pyScfEriMap.insert({{7, 7, 5, 5}, 0.799671506518246});
  pyScfEriMap.insert({{7, 7, 6, 0}, 0.21306349915127773});
  pyScfEriMap.insert({{7, 7, 6, 1}, 0.3378712021869342});
  pyScfEriMap.insert({{7, 7, 6, 2}, 0.15188539030652465});
  pyScfEriMap.insert({{7, 7, 6, 5}, 0.47249239596550685});
  pyScfEriMap.insert({{7, 7, 6, 6}, 0.4962129991093084});
  pyScfEriMap.insert({{7, 7, 7, 0}, -0.20173315659836957});
  pyScfEriMap.insert({{7, 7, 7, 1}, -0.11957743164287343});
  pyScfEriMap.insert({{7, 7, 7, 2}, -0.3099196303142542});
  pyScfEriMap.insert({{7, 7, 7, 5}, -8.293002923021744e-18});
  pyScfEriMap.insert({{7, 7, 7, 7}, 0.8242232905265838});
  pyScfEriMap.insert({{8, 0, 3, 0}, 0.015340856254019312});
  pyScfEriMap.insert({{8, 0, 3, 1}, 0.011237827091036434});
  pyScfEriMap.insert({{8, 0, 3, 2}, 0.00826051746819102});
  pyScfEriMap.insert({{8, 0, 5, 3}, 0.005579310206420053});
  pyScfEriMap.insert({{8, 0, 6, 3}, 0.008097433829685353});
  pyScfEriMap.insert({{8, 0, 7, 3}, -0.00753987665610914});
  pyScfEriMap.insert({{8, 0, 8, 0}, 0.006287425808184568});
  pyScfEriMap.insert({{8, 1, 3, 0}, 0.01927705268139468});
  pyScfEriMap.insert({{8, 1, 3, 1}, 0.01534726334549713});
  pyScfEriMap.insert({{8, 1, 3, 2}, 0.010290634742205871});
  pyScfEriMap.insert({{8, 1, 5, 3}, 0.008097433829685353});
  pyScfEriMap.insert({{8, 1, 6, 3}, 0.011398425049425152});
  pyScfEriMap.insert({{8, 1, 7, 3}, -0.010214572627416846});
  pyScfEriMap.insert({{8, 1, 8, 0}, 0.01083517470985743});
  pyScfEriMap.insert({{8, 1, 8, 1}, 0.023215675300076073});
  pyScfEriMap.insert({{8, 2, 3, 0}, 0.01851527971343749});
  pyScfEriMap.insert({{8, 2, 3, 1}, 0.013716550331242344});
  pyScfEriMap.insert({{8, 2, 3, 2}, 0.011895063090150605});
  pyScfEriMap.insert({{8, 2, 5, 3}, 0.007539876656109137});
  pyScfEriMap.insert({{8, 2, 6, 3}, 0.010214572627416846});
  pyScfEriMap.insert({{8, 2, 7, 3}, -0.010485392974308172});
  pyScfEriMap.insert({{8, 2, 8, 0}, 0.008579518074314093});
  pyScfEriMap.insert({{8, 2, 8, 1}, 0.013924675467440325});
  pyScfEriMap.insert({{8, 2, 8, 2}, 0.012134961853727665});
  pyScfEriMap.insert({{8, 3, 0, 0}, 0.08115976802285974});
  pyScfEriMap.insert({{8, 3, 1, 0}, 0.050940587328480766});
  pyScfEriMap.insert({{8, 3, 1, 1}, 0.05831519499448904});
  pyScfEriMap.insert({{8, 3, 2, 0}, 0.01851527971343749});
  pyScfEriMap.insert({{8, 3, 2, 1}, 0.013716550331242344});
  pyScfEriMap.insert({{8, 3, 2, 2}, 0.08348336220560676});
  pyScfEriMap.insert({{8, 3, 3, 3}, 0.08110953974505454});
  pyScfEriMap.insert({{8, 3, 4, 4}, 0.07476204599198895});
  pyScfEriMap.insert({{8, 3, 5, 0}, 0.02449662658682637});
  pyScfEriMap.insert({{8, 3, 5, 1}, 0.03384472430561728});
  pyScfEriMap.insert({{8, 3, 5, 2}, 0.030707902575380325});
  pyScfEriMap.insert({{8, 3, 5, 5}, 0.08115976802285974});
  pyScfEriMap.insert({{8, 3, 6, 0}, 0.03384472430561729});
  pyScfEriMap.insert({{8, 3, 6, 1}, 0.046363806828598766});
  pyScfEriMap.insert({{8, 3, 6, 2}, 0.02203387412531316});
  pyScfEriMap.insert({{8, 3, 6, 5}, 0.05094058732848077});
  pyScfEriMap.insert({{8, 3, 6, 6}, 0.05831519499448904});
  pyScfEriMap.insert({{8, 3, 7, 0}, -0.030707902575380346});
  pyScfEriMap.insert({{8, 3, 7, 1}, -0.022033874125313158});
  pyScfEriMap.insert({{8, 3, 7, 2}, -0.04472897676758341});
  pyScfEriMap.insert({{8, 3, 7, 5}, -0.01851527971343749});
  pyScfEriMap.insert({{8, 3, 7, 6}, -0.013716550331242344});
  pyScfEriMap.insert({{8, 3, 7, 7}, 0.08348336220560676});
  pyScfEriMap.insert({{8, 3, 8, 3}, 0.013471479181925924});
  pyScfEriMap.insert({{8, 4, 4, 3}, 0.003173746876532797});
  pyScfEriMap.insert({{8, 4, 8, 4}, 0.0008247844397097458});
  pyScfEriMap.insert({{8, 5, 3, 0}, 0.02982284952032839});
  pyScfEriMap.insert({{8, 5, 3, 1}, 0.02471797181214046});
  pyScfEriMap.insert({{8, 5, 3, 2}, 0.017643987788957916});
  pyScfEriMap.insert({{8, 5, 5, 3}, 0.015340856254019309});
  pyScfEriMap.insert({{8, 5, 6, 3}, 0.019277052681394676});
  pyScfEriMap.insert({{8, 5, 7, 3}, -0.018515279713437485});
  pyScfEriMap.insert({{8, 5, 8, 0}, 0.024177195126881806});
  pyScfEriMap.insert({{8, 5, 8, 1}, 0.057265343623830206});
  pyScfEriMap.insert({{8, 5, 8, 2}, 0.030697482625577373});
  pyScfEriMap.insert({{8, 5, 8, 5}, 0.16045646128458113});
  pyScfEriMap.insert({{8, 6, 3, 0}, 0.024717971812140468});
  pyScfEriMap.insert({{8, 6, 3, 1}, 0.02029272928143766});
  pyScfEriMap.insert({{8, 6, 3, 2}, 0.013338197657181184});
  pyScfEriMap.insert({{8, 6, 5, 3}, 0.011237827091036429});
  pyScfEriMap.insert({{8, 6, 6, 3}, 0.01534726334549713});
  pyScfEriMap.insert({{8, 6, 7, 3}, -0.013716550331242345});
  pyScfEriMap.insert({{8, 6, 8, 0}, 0.016460151691831012});
  pyScfEriMap.insert({{8, 6, 8, 1}, 0.03902326337821374});
  pyScfEriMap.insert({{8, 6, 8, 2}, 0.020574279542097167});
  pyScfEriMap.insert({{8, 6, 8, 5}, 0.10166761452628235});
  pyScfEriMap.insert({{8, 6, 8, 6}, 0.06893333251070229});
  pyScfEriMap.insert({{8, 7, 3, 0}, -0.017643987788957916});
  pyScfEriMap.insert({{8, 7, 3, 1}, -0.013338197657181184});
  pyScfEriMap.insert({{8, 7, 3, 2}, -0.012602070194764439});
  pyScfEriMap.insert({{8, 7, 5, 3}, -0.008260517468191014});
  pyScfEriMap.insert({{8, 7, 6, 3}, -0.010290634742205871});
  pyScfEriMap.insert({{8, 7, 7, 3}, 0.011895063090150603});
  pyScfEriMap.insert({{8, 7, 8, 0}, -0.008863018633008183});
  pyScfEriMap.insert({{8, 7, 8, 1}, -0.0066962020164952896});
  pyScfEriMap.insert({{8, 7, 8, 2}, -0.014217679958253654});
  pyScfEriMap.insert({{8, 7, 8, 5}, -2.1359882077647224e-19});
  pyScfEriMap.insert({{8, 7, 8, 7}, 0.05046265044040281});
  pyScfEriMap.insert({{8, 8, 0, 0}, 0.4117361160007401});
  pyScfEriMap.insert({{8, 8, 1, 0}, 0.2771783914805482});
  pyScfEriMap.insert({{8, 8, 1, 1}, 0.3667708074445731});
  pyScfEriMap.insert({{8, 8, 2, 0}, 0.09474336352639351});
  pyScfEriMap.insert({{8, 8, 2, 1}, 0.08022143347148329});
  pyScfEriMap.insert({{8, 8, 2, 2}, 0.44399727541714773});
  pyScfEriMap.insert({{8, 8, 3, 3}, 0.40033533712023217});
  pyScfEriMap.insert({{8, 8, 4, 4}, 0.39197234092357564});
  pyScfEriMap.insert({{8, 8, 5, 0}, 0.15838231560071211});
  pyScfEriMap.insert({{8, 8, 5, 1}, 0.27736521609567444});
  pyScfEriMap.insert({{8, 8, 5, 2}, 0.19942022994504377});
  pyScfEriMap.insert({{8, 8, 5, 5}, 0.799671506518246});
  pyScfEriMap.insert({{8, 8, 6, 0}, 0.19430030777437562});
  pyScfEriMap.insert({{8, 8, 6, 1}, 0.33336674981076697});
  pyScfEriMap.insert({{8, 8, 6, 2}, 0.13349755191025142});
  pyScfEriMap.insert({{8, 8, 6, 5}, 0.47249239596550685});
  pyScfEriMap.insert({{8, 8, 6, 6}, 0.4962129991093084});
  pyScfEriMap.insert({{8, 8, 7, 0}, -0.17197648182006695});
  pyScfEriMap.insert({{8, 8, 7, 1}, -0.10557586217126812});
  pyScfEriMap.insert({{8, 8, 7, 2}, -0.26055742543274174});
  pyScfEriMap.insert({{8, 8, 7, 5}, -7.8658052814688e-18});
  pyScfEriMap.insert({{8, 8, 7, 7}, 0.7232979896457781});
  pyScfEriMap.insert({{8, 8, 8, 3}, 0.08110953974505454});
  pyScfEriMap.insert({{8, 8, 8, 8}, 0.8242232905265838});
  pyScfEriMap.insert({{9, 0, 4, 0}, 0.015340856254019312});
  pyScfEriMap.insert({{9, 0, 4, 1}, 0.011237827091036434});
  pyScfEriMap.insert({{9, 0, 4, 2}, 0.00826051746819102});
  pyScfEriMap.insert({{9, 0, 5, 4}, 0.005579310206420053});
  pyScfEriMap.insert({{9, 0, 6, 4}, 0.008097433829685353});
  pyScfEriMap.insert({{9, 0, 7, 4}, -0.007539876656109141});
  pyScfEriMap.insert({{9, 0, 9, 0}, 0.006287425808184568});
  pyScfEriMap.insert({{9, 1, 4, 0}, 0.01927705268139468});
  pyScfEriMap.insert({{9, 1, 4, 1}, 0.01534726334549713});
  pyScfEriMap.insert({{9, 1, 4, 2}, 0.010290634742205871});
  pyScfEriMap.insert({{9, 1, 5, 4}, 0.008097433829685353});
  pyScfEriMap.insert({{9, 1, 6, 4}, 0.011398425049425152});
  pyScfEriMap.insert({{9, 1, 7, 4}, -0.010214572627416846});
  pyScfEriMap.insert({{9, 1, 9, 0}, 0.01083517470985743});
  pyScfEriMap.insert({{9, 1, 9, 1}, 0.023215675300076073});
  pyScfEriMap.insert({{9, 2, 4, 0}, 0.01851527971343749});
  pyScfEriMap.insert({{9, 2, 4, 1}, 0.013716550331242344});
  pyScfEriMap.insert({{9, 2, 4, 2}, 0.011895063090150605});
  pyScfEriMap.insert({{9, 2, 5, 4}, 0.0075398766561091375});
  pyScfEriMap.insert({{9, 2, 6, 4}, 0.010214572627416846});
  pyScfEriMap.insert({{9, 2, 7, 4}, -0.010485392974308172});
  pyScfEriMap.insert({{9, 2, 9, 0}, 0.008579518074314093});
  pyScfEriMap.insert({{9, 2, 9, 1}, 0.013924675467440325});
  pyScfEriMap.insert({{9, 2, 9, 2}, 0.012134961853727663});
  pyScfEriMap.insert({{9, 3, 4, 3}, 0.003173746876532797});
  pyScfEriMap.insert({{9, 3, 8, 4}, 0.0008247844397097458});
  pyScfEriMap.insert({{9, 3, 9, 3}, 0.0008247844397097458});
  pyScfEriMap.insert({{9, 4, 0, 0}, 0.08115976802285974});
  pyScfEriMap.insert({{9, 4, 1, 0}, 0.050940587328480766});
  pyScfEriMap.insert({{9, 4, 1, 1}, 0.05831519499448904});
  pyScfEriMap.insert({{9, 4, 2, 0}, 0.01851527971343749});
  pyScfEriMap.insert({{9, 4, 2, 1}, 0.013716550331242344});
  pyScfEriMap.insert({{9, 4, 2, 2}, 0.08348336220560676});
  pyScfEriMap.insert({{9, 4, 3, 3}, 0.07476204599198895});
  pyScfEriMap.insert({{9, 4, 4, 4}, 0.08110953974505454});
  pyScfEriMap.insert({{9, 4, 5, 0}, 0.02449662658682637});
  pyScfEriMap.insert({{9, 4, 5, 1}, 0.03384472430561728});
  pyScfEriMap.insert({{9, 4, 5, 2}, 0.030707902575380325});
  pyScfEriMap.insert({{9, 4, 5, 5}, 0.08115976802285974});
  pyScfEriMap.insert({{9, 4, 6, 0}, 0.03384472430561729});
  pyScfEriMap.insert({{9, 4, 6, 1}, 0.046363806828598766});
  pyScfEriMap.insert({{9, 4, 6, 2}, 0.022033874125313158});
  pyScfEriMap.insert({{9, 4, 6, 5}, 0.05094058732848077});
  pyScfEriMap.insert({{9, 4, 6, 6}, 0.05831519499448904});
  pyScfEriMap.insert({{9, 4, 7, 0}, -0.030707902575380343});
  pyScfEriMap.insert({{9, 4, 7, 1}, -0.022033874125313158});
  pyScfEriMap.insert({{9, 4, 7, 2}, -0.04472897676758341});
  pyScfEriMap.insert({{9, 4, 7, 5}, -0.01851527971343749});
  pyScfEriMap.insert({{9, 4, 7, 6}, -0.013716550331242344});
  pyScfEriMap.insert({{9, 4, 7, 7}, 0.08348336220560676});
  pyScfEriMap.insert({{9, 4, 8, 3}, 0.011821910302506431});
  pyScfEriMap.insert({{9, 4, 8, 8}, 0.07476204599198895});
  pyScfEriMap.insert({{9, 4, 9, 4}, 0.013471479181925924});
  pyScfEriMap.insert({{9, 5, 4, 0}, 0.02982284952032839});
  pyScfEriMap.insert({{9, 5, 4, 1}, 0.02471797181214046});
  pyScfEriMap.insert({{9, 5, 4, 2}, 0.017643987788957916});
  pyScfEriMap.insert({{9, 5, 5, 4}, 0.015340856254019309});
  pyScfEriMap.insert({{9, 5, 6, 4}, 0.019277052681394676});
  pyScfEriMap.insert({{9, 5, 7, 4}, -0.018515279713437485});
  pyScfEriMap.insert({{9, 5, 9, 0}, 0.024177195126881806});
  pyScfEriMap.insert({{9, 5, 9, 1}, 0.057265343623830206});
  pyScfEriMap.insert({{9, 5, 9, 2}, 0.030697482625577373});
  pyScfEriMap.insert({{9, 5, 9, 5}, 0.16045646128458113});
  pyScfEriMap.insert({{9, 6, 4, 0}, 0.024717971812140468});
  pyScfEriMap.insert({{9, 6, 4, 1}, 0.02029272928143766});
  pyScfEriMap.insert({{9, 6, 4, 2}, 0.013338197657181186});
  pyScfEriMap.insert({{9, 6, 5, 4}, 0.011237827091036429});
  pyScfEriMap.insert({{9, 6, 6, 4}, 0.01534726334549713});
  pyScfEriMap.insert({{9, 6, 7, 4}, -0.013716550331242345});
  pyScfEriMap.insert({{9, 6, 9, 0}, 0.016460151691831012});
  pyScfEriMap.insert({{9, 6, 9, 1}, 0.03902326337821374});
  pyScfEriMap.insert({{9, 6, 9, 2}, 0.020574279542097167});
  pyScfEriMap.insert({{9, 6, 9, 5}, 0.10166761452628235});
  pyScfEriMap.insert({{9, 6, 9, 6}, 0.06893333251070229});
  pyScfEriMap.insert({{9, 7, 4, 0}, -0.017643987788957916});
  pyScfEriMap.insert({{9, 7, 4, 1}, -0.013338197657181186});
  pyScfEriMap.insert({{9, 7, 4, 2}, -0.012602070194764439});
  pyScfEriMap.insert({{9, 7, 5, 4}, -0.008260517468191014});
  pyScfEriMap.insert({{9, 7, 6, 4}, -0.010290634742205871});
  pyScfEriMap.insert({{9, 7, 7, 4}, 0.011895063090150603});
  pyScfEriMap.insert({{9, 7, 9, 0}, -0.008863018633008183});
  pyScfEriMap.insert({{9, 7, 9, 1}, -0.0066962020164952896});
  pyScfEriMap.insert({{9, 7, 9, 2}, -0.014217679958253654});
  pyScfEriMap.insert({{9, 7, 9, 5}, -2.1359882077647224e-19});
  pyScfEriMap.insert({{9, 7, 9, 7}, 0.050462650440402795});
  pyScfEriMap.insert({{9, 8, 4, 3}, 0.004181498098328303});
  pyScfEriMap.insert({{9, 8, 8, 4}, 0.003173746876532797});
  pyScfEriMap.insert({{9, 8, 9, 3}, 0.003173746876532797});
  pyScfEriMap.insert({{9, 8, 9, 8}, 0.050462650440402795});
  pyScfEriMap.insert({{9, 9, 0, 0}, 0.4117361160007401});
  pyScfEriMap.insert({{9, 9, 1, 0}, 0.2771783914805482});
  pyScfEriMap.insert({{9, 9, 1, 1}, 0.3667708074445731});
  pyScfEriMap.insert({{9, 9, 2, 0}, 0.09474336352639351});
  pyScfEriMap.insert({{9, 9, 2, 1}, 0.0802214334714833});
  pyScfEriMap.insert({{9, 9, 2, 2}, 0.44399727541714773});
  pyScfEriMap.insert({{9, 9, 3, 3}, 0.39197234092357564});
  pyScfEriMap.insert({{9, 9, 4, 4}, 0.4003353371202322});
  pyScfEriMap.insert({{9, 9, 5, 0}, 0.15838231560071211});
  pyScfEriMap.insert({{9, 9, 5, 1}, 0.27736521609567444});
  pyScfEriMap.insert({{9, 9, 5, 2}, 0.19942022994504377});
  pyScfEriMap.insert({{9, 9, 5, 5}, 0.799671506518246});
  pyScfEriMap.insert({{9, 9, 6, 0}, 0.19430030777437562});
  pyScfEriMap.insert({{9, 9, 6, 1}, 0.33336674981076697});
  pyScfEriMap.insert({{9, 9, 6, 2}, 0.13349755191025142});
  pyScfEriMap.insert({{9, 9, 6, 5}, 0.47249239596550685});
  pyScfEriMap.insert({{9, 9, 6, 6}, 0.4962129991093084});
  pyScfEriMap.insert({{9, 9, 7, 0}, -0.17197648182006695});
  pyScfEriMap.insert({{9, 9, 7, 1}, -0.10557586217126813});
  pyScfEriMap.insert({{9, 9, 7, 2}, -0.26055742543274174});
  pyScfEriMap.insert({{9, 9, 7, 5}, -7.8658052814688e-18});
  pyScfEriMap.insert({{9, 9, 7, 7}, 0.7232979896457781});
  pyScfEriMap.insert({{9, 9, 8, 3}, 0.07476204599198895});
  pyScfEriMap.insert({{9, 9, 8, 8}, 0.7232979896457781});
  pyScfEriMap.insert({{9, 9, 9, 4}, 0.08110953974505454});
  pyScfEriMap.insert({{9, 9, 9, 9}, 0.8242232905265837});

  LibintIntegrals eval;

  Utils::Settings& settings = eval.settings();

  settings.modifyBool("use_pure_spherical", false);

  auto basis = eval.initializeBasisSet("def2-svp", scineAtoms);

  Utils::Integrals::IntegralSpecifier specifier;
  specifier.op = Utils::Integrals::Operator::Coulomb;

  auto result_map = LibintIntegrals::evaluate(specifier, basis, basis);

  const auto& result = result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];

  int dim = static_cast<int>(basis.nbf());

  for (auto const& elem : pyScfEriMap) {
    EXPECT_THAT(result(elem.first.at(0) * dim + elem.first.at(1), elem.first.at(2) * dim + elem.first.at(3)),
                DoubleNear(elem.second, 1e-8));
  }
}

TEST_F(TwoBodyIntsTest, Test2body8foldSymmetry) {
  std::stringstream h2("5\n\n"
                       "H -2.4 0 0\n"
                       "H -1.2 0 0\n"
                       "H 0 0 0\n"
                       "H 1.2 0 0\n"
                       "H 2.4 0 0\n");

  auto scineAtoms = Utils::XyzStreamHandler::read(h2);

  LibintIntegrals eval;

  Utils::Settings& settings = eval.settings();

  settings.modifyBool("use_pure_spherical", false);

  auto basis = eval.initializeBasisSet("def2-svp", scineAtoms);
  Utils::Integrals::IntegralSpecifier specifier;
  specifier.op = Utils::Integrals::Operator::Coulomb;

  auto result_map = LibintIntegrals::evaluate(specifier, basis, basis);

  const auto& result = result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];

  int dim = static_cast<int>(basis.nbf());

  int ij;
  int ji;
  int kl;
  int lk;
  double res;
  // Loop structure for 8-fold symmetry
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j <= i; ++j) {
      for (int k = 0; k <= i; ++k) {
        const auto l_max = (i == k) ? j : k;
        for (int l = 0; l <= l_max; ++l) {
          ij = i * dim + j;
          kl = k * dim + l;
          ji = j * dim + i;
          lk = l * dim + k;
          res = result(ij, kl);
          EXPECT_THAT(result(ij, lk), DoubleNear(res, 1e-8));
          EXPECT_THAT(result(ji, kl), DoubleNear(res, 1e-8));
          EXPECT_THAT(result(ji, lk), DoubleNear(res, 1e-8));
          EXPECT_THAT(result(kl, ij), DoubleNear(res, 1e-8));
          EXPECT_THAT(result(kl, ji), DoubleNear(res, 1e-8));
          EXPECT_THAT(result(lk, ij), DoubleNear(res, 1e-8));
          EXPECT_THAT(result(lk, ji), DoubleNear(res, 1e-8));
        }
      }
    }
  }
}

TEST_F(TwoBodyIntsTest, Test2body4foldSymmetry) {
  std::stringstream h2_1("2\n\n"
                         "H 0 0 0\n"
                         "H 1.2 0 0");
  std::stringstream h2_2("2\n\n"
                         "H 0.5 0.5 0.5\n"
                         "H 1.7 0.5 0.5");

  auto scineAtoms1 = Utils::XyzStreamHandler::read(h2_1);

  auto scineAtoms2 = Utils::XyzStreamHandler::read(h2_2);

  LibintIntegrals eval;

  Utils::Settings& settings = eval.settings();

  settings.modifyBool("use_pure_spherical", false);

  auto basis1 = eval.initializeBasisSet("def2-svp", scineAtoms1);
  auto basis2 = eval.initializeBasisSet("def2-tzvp", scineAtoms2);

  Utils::Integrals::IntegralSpecifier specifier;
  specifier.op = Utils::Integrals::Operator::Coulomb;

  auto result_map = LibintIntegrals::evaluate(specifier, basis1, basis2);

  const auto& result = result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];

  int dim1 = static_cast<int>(basis1.nbf());
  int dim2 = static_cast<int>(basis2.nbf());

  int ij;
  int ji;
  int kl;
  int lk;
  double res;
  // Loop structure for 4-fold symmetry
  for (int i = 0; i < dim1; ++i) {
    for (int j = 0; j <= i; ++j) {
      for (int k = 0; k < dim2; ++k) {
        for (int l = 0; l <= k; ++l) {
          ij = i * dim1 + j;
          kl = k * dim2 + l;
          ji = j * dim1 + i;
          lk = l * dim2 + k;
          res = result(ij, kl);
          // std::cout << i << " " << j << " " << k << " " << l << " " << std::endl;
          EXPECT_THAT(result(ij, lk), DoubleNear(res, 1e-8));
          // std::cout << j << " " << i << " " << k << " " << l << " " << std::endl;
          EXPECT_THAT(result(ji, kl), DoubleNear(res, 1e-8));
          // std::cout << j << " " << i << " " << l << " " << k << " " << std::endl;
          EXPECT_THAT(result(ji, lk), DoubleNear(res, 1e-8));
        }
      }
    }
  }
}

TEST_F(TwoBodyIntsTest, Test2body2foldSymmetry) {
  std::stringstream h2_1("2\n\n"
                         "H 0 0 0\n"
                         "H 1.2 0 0");
  std::stringstream h2_2("2\n\n"
                         "H 0.5 0.5 0.5\n"
                         "H 1.7 0.5 0.5");

  auto scineAtoms1 = Utils::XyzStreamHandler::read(h2_1);

  auto scineAtoms2 = Utils::XyzStreamHandler::read(h2_2);

  LibintIntegrals eval;

  Utils::Settings& settings = eval.settings();

  settings.modifyBool("use_pure_spherical", false);

  auto basis1 = eval.initializeBasisSet("def2-svp", scineAtoms1);
  auto basis2 = eval.initializeBasisSet("def2-tzvp", scineAtoms2);

  Utils::Integrals::IntegralSpecifier specifier;
  specifier.op = Utils::Integrals::Operator::CoulombCOM;
  specifier.totalMass = 42;

  auto result_map = LibintIntegrals::evaluate(specifier, basis1, basis2);

  const auto& result = result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];

  int dim1 = static_cast<int>(basis1.nbf());
  int dim2 = static_cast<int>(basis2.nbf());

  int ij;
  int ji;
  int kl;
  int lk;
  double res;
  // Loop structure for 1-fold symmetry
  for (int i = 0; i < dim1; ++i) {
    for (int j = 0; j < i; ++j) {
      for (int k = 0; k < dim2; ++k) {
        for (int l = 0; l < k; ++l) {
          ij = i * dim1 + j;
          kl = k * dim2 + l;
          ji = j * dim1 + i;
          lk = l * dim2 + k;
          res = result(ij, kl);
          // std::cout << i << " " << j << " " << k << " " << l << " " << std::endl;
          EXPECT_THAT(result(ji, lk), DoubleNear(res, 1e-8));
        }
      }
    }
  }
}

TEST_F(TwoBodyIntsTest, TestCoulombDerivative) {
  // Reference

  std::stringstream h2("2\n\n"
                       "H 0 0 0\n"
                       "H 1.2 0 0");
  auto atoms = Utils::XyzStreamHandler::read(h2);

  // Central difference:
  std::stringstream h2_h0("2\n\n"
                          "H 0 0 0\n"
                          "H 1.2001 0 0");
  auto atoms_h0 = Utils::XyzStreamHandler::read(h2_h0);

  std::stringstream h2_h1("2\n\n"
                          "H 0 0 0\n"
                          "H 1.1999 0 0");
  auto atoms_h1 = Utils::XyzStreamHandler::read(h2_h1);

  LibintIntegrals eval;

  Utils::Settings& settings = eval.settings();

  settings.modifyBool("use_pure_spherical", false);

  std::string name = "def2-svp";

  auto basis = eval.initializeBasisSet(name, atoms);
  Utils::Integrals::IntegralSpecifier specifier;
  specifier.op = Utils::Integrals::Operator::Coulomb;
  specifier.derivOrder = 1;

  auto result_map = LibintIntegrals::evaluate(specifier, basis, basis);
  auto result = result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::x, 0}];
  result += result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::x, 1}];

  Utils::Integrals::IntegralSpecifier testSpecifier;
  testSpecifier.op = Utils::Integrals::Operator::Coulomb;
  std::vector<Eigen::MatrixXd> centralDifference;

  // Derivative in x-direction:
  auto tmpBasis = eval.initializeBasisSet(name, atoms_h0);
  auto tmp_map = LibintIntegrals::evaluate(testSpecifier, tmpBasis, tmpBasis);
  auto tmp_res = tmp_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
  centralDifference.push_back(tmp_res);
  tmpBasis = eval.initializeBasisSet(name, atoms_h1);
  tmp_map = LibintIntegrals::evaluate(testSpecifier, tmpBasis, tmpBasis);
  tmp_res = tmp_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
  centralDifference.push_back(tmp_res);

  Eigen::MatrixXd xDerivative =
      1 / (0.0002 * Utils::Constants::bohr_per_angstrom) * (centralDifference[0] - centralDifference[1]);

  auto shell2bf = basis.shell2bf();

  size_t i;
  size_t j;
  size_t k;
  size_t l;
  size_t ij;
  size_t kl;
  auto dim = basis.nbf();
  // loop over permutationally-unique set of shells
  for (size_t s1 = 0; s1 < basis.size(); ++s1) {
    for (size_t s2 = 0; s2 <= s1; ++s2) {
      for (size_t s3 = 0; s3 <= s1; ++s3) {
        const auto s4_max = (s1 == s3) ? s2 : s3;
        for (size_t s4 = 0; s4 <= s4_max; ++s4) {
          auto bf1 = shell2bf[s1];    // first basis function in first shell
          auto n1 = basis[s1].size(); // number of basis functions in first shell
          auto bf2 = shell2bf[s2];    // first basis function in 2nd shell
          auto n2 = basis[s2].size(); // number of basis functions in 2nd shell
          auto bf3 = shell2bf[s3];    // first basis function in 3rd shell
          auto n3 = basis[s3].size(); // number of basis functions in 3d shell
          auto bf4 = shell2bf[s4];    // first basis function in 4th shell
          auto n4 = basis[s4].size(); // number of basis functions in 4th shell
          for (size_t f1 = 0, f1234 = 0; f1 < n1; ++f1)
            for (size_t f2 = 0; f2 < n2; ++f2)
              for (size_t f3 = 0; f3 < n3; ++f3)
                for (size_t f4 = 0; f4 < n4; ++f4, ++f1234) {
                  i = bf1 + f1;
                  j = bf2 + f2;
                  k = bf3 + f3;
                  l = bf4 + f4;
                  ij = i * dim + j;
                  kl = k * dim + l;
                  if (basis[s1].getShift() == basis[s2].getShift() && basis[s1].getShift() == basis[s3].getShift() &&
                      basis[s1].getShift() == basis[s4].getShift()) {
                    auto tmp = result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::x, 0}](ij, kl);
                    tmp += result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::x, 1}](ij, kl);
                    tmp += result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::x, 2}](ij, kl);
                    tmp += result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::x, 3}](ij, kl);
                    EXPECT_THAT(xDerivative(ij, kl), DoubleNear(tmp, 1e-6));
                  }
                  else if (basis[s1].getShift() == basis[s2].getShift() &&
                           basis[s1].getShift() == basis[s3].getShift() && basis[s1].getShift() != basis[s4].getShift()) {
                    auto tmp = result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::x, 0}](ij, kl);
                    tmp += result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::x, 1}](ij, kl);
                    tmp += result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::x, 2}](ij, kl);
                    EXPECT_THAT(xDerivative(ij, kl), DoubleNear(tmp, 1e-6));
                  }
                  else if (basis[s1].getShift() == basis[s2].getShift() &&
                           basis[s1].getShift() != basis[s3].getShift() && basis[s1].getShift() != basis[s4].getShift()) {
                    auto tmp = result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::x, 0}](ij, kl);
                    tmp += result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::x, 1}](ij, kl);
                    EXPECT_THAT(xDerivative(ij, kl), DoubleNear(tmp, 1e-6));
                  }
                }
        }
      }
    }
  }
}
