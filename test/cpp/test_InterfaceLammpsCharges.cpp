#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE interface_lammps_charges

#include "InterfaceLammps.h"

#include <boost/test/unit_test.hpp>

#include <vector>

using namespace nnp;

namespace
{
class InterfaceLammpsChargesHarness : public InterfaceLammps
{
public:
    InterfaceLammpsChargesHarness()
    {
        numElements = 1;
        normalize = false;
        convCharge = 1.0;
        structure.numAtomsPerElement.assign(numElements, 0);
    }

    void prepare(NNPType type, std::vector<double> const& charges)
    {
        nnpType = type;
        structure.numAtoms = charges.size();
        structure.atoms.clear();
        structure.atoms.reserve(charges.size());
        structure.numAtomsPerElement.assign(numElements, 0);

        for (size_t i = 0; i < charges.size(); ++i)
        {
            Atom atom;
            atom.index = i;
            atom.element = 0;
            atom.tag = static_cast<int64_t>(i + 1);
            atom.charge = charges[i];
            structure.atoms.push_back(atom);
            structure.numAtomsPerElement[0]++;
        }
    }

    void setAtomIndex(size_t pos, size_t index)
    {
        structure.atoms.at(pos).index = index;
    }
};
} // namespace

BOOST_AUTO_TEST_CASE(get_charges_returns_values_for_q)
{
    InterfaceLammpsChargesHarness interface;
    std::vector<double> charges{0.42, -0.17, 0.05};
    interface.prepare(InterfaceLammps::NNPType::HDNNP_Q, charges);

    std::vector<double> out(charges.size(), 0.0);
    interface.getCharges(out.data());

    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), charges.begin(), charges.end());
}

BOOST_AUTO_TEST_CASE(get_charges_ignores_non_charge_types)
{
    InterfaceLammpsChargesHarness interface;
    std::vector<double> charges{1.0, -1.0};
    interface.prepare(InterfaceLammps::NNPType::HDNNP_2G, charges);

    std::vector<double> out(charges.size(), 5.0);
    interface.getCharges(out.data());

    BOOST_CHECK_CLOSE(out[0], 5.0, 1.0e-12);
    BOOST_CHECK_CLOSE(out[1], 5.0, 1.0e-12);
}

BOOST_AUTO_TEST_CASE(get_charges_respects_original_indices)
{
    InterfaceLammpsChargesHarness interface;
    std::vector<double> charges{0.3, -0.6};
    interface.prepare(InterfaceLammps::NNPType::HDNNP_Q, charges);
    interface.setAtomIndex(0, 3);
    interface.setAtomIndex(1, 1);

    std::vector<double> out(4, 0.0);
    interface.getCharges(out.data());

    BOOST_CHECK_SMALL(out[0], 1.0e-12);
    BOOST_CHECK_CLOSE(out[1], -0.6, 1.0e-12);
    BOOST_CHECK_SMALL(out[2], 1.0e-12);
    BOOST_CHECK_CLOSE(out[3], 0.3, 1.0e-12);
}
