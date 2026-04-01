from helena import enumerate_variables, enumerate_vectors, variable_unit_conversion
import pytest


@pytest.mark.parametrize("variable", ["AR2+", "something_not_in_header"])
def test_enumerate_variables(variable, long_header, snaptolshot):
    assert snaptolshot == enumerate_variables(variable, long_header)


@pytest.mark.parametrize("variable", ["AR2+", "E FLUX-R", "E FLUX-Z"])
def test_enumerate_vectors(variable, long_header, snaptolshot):
    assert snaptolshot == enumerate_vectors(variable, long_header, prefixes=None)


@pytest.mark.parametrize(
    "variable, units",
    [
        ("E FLUX-Z", "SI"),
        ("E FLUX-Z", "CGS"),
        ("VR-NEUTRAL", "SI"),
        ("VR-NEUTRAL", "CGS"),
        ("VR-ION+", "SI"),
        ("VR-ION+", "CGS"),
        ("ERADIAL", "SI"),
        ("EDRADIAL", "CGS"),
        ("RHO", "SI"),
        ("RHO", "CGS"),
        ("POW-ALL", "SI"),
        ("POW-ALL", "CGS"),
    ],
)
def test_variable_unit_conversion(variable, units, long_profile, snaptolshot):
    snaptolshot.assert_allclose(
        variable_unit_conversion(long_profile, variable, units, atomic_species=None)
    )
