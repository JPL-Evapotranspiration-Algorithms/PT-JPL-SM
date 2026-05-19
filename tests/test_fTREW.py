import numpy as np
import pytest

from PTJPLSM.fTREW import calculate_fTREW


def _legacy_inline_fTREW(PET_Wm2, canopy_height_meters, soil_moisture, field_capacity, wilting_point, canopy_buffer_sensitivity):
    stress_onset_weight = (1 / (1 + PET_Wm2)) - (canopy_buffer_sensitivity / (1 + canopy_height_meters))
    CHscalar = np.sqrt(canopy_height_meters)

    with np.errstate(divide='ignore', invalid='ignore'):
        WPCH = np.clip(np.where(CHscalar == 0, 0, wilting_point / CHscalar), 0, 1)
        CR = (1 - stress_onset_weight) * (field_capacity - WPCH) + WPCH
        fTREW = np.clip(1 - ((CR - soil_moisture) / (CR - WPCH)) ** CHscalar, 0, 1)

    return np.where(np.isnan(fTREW), 0, fTREW)


def test_calculate_fTREW_matches_legacy_inline_behavior():
    PET_Wm2 = np.array([0.5, 1.2, 2.1], dtype=float)
    canopy_height_meters = np.array([0.5, 2.0, 5.0], dtype=float)
    soil_moisture = np.array([0.10, 0.25, 0.35], dtype=float)
    field_capacity = np.array([0.30, 0.35, 0.40], dtype=float)
    wilting_point = np.array([0.05, 0.10, 0.15], dtype=float)

    expected = _legacy_inline_fTREW(
        PET_Wm2,
        canopy_height_meters,
        soil_moisture,
        field_capacity,
        wilting_point,
        canopy_buffer_sensitivity=0.1,
    )
    actual = calculate_fTREW(
        PET_Wm2=PET_Wm2,
        canopy_height_meters=canopy_height_meters,
        soil_moisture=soil_moisture,
        field_capacity=field_capacity,
        wilting_point=wilting_point,
        canopy_buffer_sensitivity=0.1,
    )

    np.testing.assert_allclose(actual, expected)
    assert np.all((actual >= 0) & (actual <= 1))


def test_calculate_fTREW_replaces_nan_with_zero():
    PET_Wm2 = np.array([0.0], dtype=float)
    canopy_height_meters = np.array([1.0], dtype=float)
    soil_moisture = np.array([0.2], dtype=float)
    field_capacity = np.array([0.2], dtype=float)
    wilting_point = np.array([0.2], dtype=float)

    f_trew = calculate_fTREW(
        PET_Wm2=PET_Wm2,
        canopy_height_meters=canopy_height_meters,
        soil_moisture=soil_moisture,
        field_capacity=field_capacity,
        wilting_point=wilting_point,
    )

    assert np.isfinite(f_trew[0])
    assert f_trew[0] == 0


def test_calculate_fTREW_validates_sensitivity_range():
    with pytest.raises(ValueError):
        calculate_fTREW(
            PET_Wm2=np.array([1.0]),
            canopy_height_meters=np.array([1.0]),
            soil_moisture=np.array([0.2]),
            field_capacity=np.array([0.3]),
            wilting_point=np.array([0.1]),
            canopy_buffer_sensitivity=1.5,
        )
