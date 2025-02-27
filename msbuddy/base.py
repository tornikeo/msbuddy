# ==============================================================================
# Copyright (C) 2024 Shipei Xing <philipxsp@hotmail.com>
#
# Licensed under the Apache License 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at https://github.com/Philipbear/msbuddy/blob/main/LICENSE
# ==============================================================================
"""
File: base.py
Author: Shipei Xing
Email: philipxsp@hotmail.com
GitHub: Philipbear
Description: base classes for msbuddy
"""

import logging
from typing import Union, Tuple, List

import numpy as np
from numba import njit

from msbuddy.utils import read_formula, form_arr_to_str

mass_i = 1.0033548  # mass of neutron
mass_e = 0.0005486  # mass of electron

logging.basicConfig(level=logging.WARNING)


class Formula:
    alphabet = ["C", "H", "Br", "Cl", "F", "I", "K", "N", "Na", "O", "P", "S"]
    mass_arr = np.array([12.000000, 1.007825, 78.918336, 34.968853, 18.998403, 126.904473, 38.963707, 14.003074,
                         22.989769, 15.994915, 30.973762, 31.972071], dtype=np.float32)

    # array: np.array(int)
    # charge: int
    # mass: float
    # dbe: float
    # isotope: int # M+0, M+1

    def __init__(self,
                 array: np.array,
                 charge: int,
                 mass: Union[float, None] = None,
                 isotope: int = 0):
        self.array = np.int16(array)
        self.charge = charge
        self.dbe = calc_formula_dbe(array)
        self.isotope = isotope

        # fill in mass directly from formula database, otherwise calculate
        if mass is None:
            self.mass = calc_formula_mass(np.float32(array), charge, isotope)
        else:
            self.mass = mass

    def __bool__(self):
        # if all elements are 0, then return False
        return bool(np.any(self.array))

    def __str__(self) -> str:
        if not self:
            return 'Null'
        return form_arr_to_str(self.array)


@njit
def calc_formula_dbe(arr):
    """
    calculate dbe of a formula
    :return: dbe
    """
    # ["C", "H", "Br", "Cl", "F", "I", "K", "N", "Na", "O", "P", "S"]
    # DBE: C + 1 - (H + F + Cl + Br + I + Na + K)/2 + (N + P)/2
    dbe = arr[0] + 1 - (arr[1] + arr[4] + arr[3] + arr[2] + arr[5] + arr[8] + arr[6]) / 2 + (arr[7] + arr[10]) / 2
    return dbe


@njit
def calc_formula_mass(array, charge, isotope):
    """
    calculate monoisotopic mass of a formula, charge adjusted
    :return: mass
    """
    ele_mass_arr = np.array([12.000000, 1.007825, 78.918336, 34.968853, 18.998403, 126.904473, 38.963707, 14.003074,
                             22.989769, 15.994915, 30.973762, 31.972071], dtype=np.float32)
    if charge == 0:
        mass = float(np.sum(array * ele_mass_arr) + isotope * mass_i)
    else:
        mass = float((np.sum(array * ele_mass_arr) - charge * mass_e + isotope * mass_i) / abs(charge))
    return mass


class Spectrum:
    def __init__(self, mz_array: np.array, int_array: np.array):
        """
        :param mz_array: np.array
        :param int_array: np.array
        """
        # convert to np array
        mz_array = np.array(mz_array, dtype=np.float32)
        int_array = np.array(int_array, dtype=np.float32)

        # mz_array and int_array must have the same length
        if len(mz_array) != len(int_array):
            raise ValueError("mz_array and int_array must have the same length")

        # sort by mz_array
        idx = np.argsort(mz_array)
        self.mz_array = np.float32(mz_array[idx])
        self.int_array = np.float32(int_array[idx])

    def __str__(self) -> str:
        return f"Spectrum(mz={self.mz_array}, int={self.int_array})"

    def __bool__(self):
        return self.mz_array.size != 0

    def __len__(self):
        return self.mz_array.size


class Adduct:
    """
    Adduct class
    :attr string: str
    :attr pos_mode: bool
    :attr charge: int
    :attr m: int
    :attr net_formula: Formula
    :attr loss_formula: Formula, the sum loss formula within the adduct
    """

    # for example, [M+H-H2O]+ has loss formula H2O; [M+Na]+ has no loss formula;
    # [M - 2H2O - Cl] - has loss formula H4O2Cl

    def __init__(self,
                 string: Union[str, None],
                 pos_mode: bool,
                 report_invalid: bool = False):
        """
        :param string: adduct string, e.g. [M+H]+
        :param pos_mode: True for positive mode, False for negative mode
        """
        if string is None:
            self.string = "[M+H]+" if pos_mode else "[M-H]-"
        else:
            self.string = string.replace(" ", "")

        self.pos_mode = pos_mode

        if self._check_common():
            pass
        elif self._check_valid_style():
            self._exchange()
            if self._check_valid_character():
                self._calc_m()
                self._calc_charge()
                self._calc_loss_and_net_formula()
            else:
                self._invalid(string, report_invalid)
        else:
            self._invalid(string, report_invalid)

    def __str__(self):
        return f"{self.string}; charge: {self.charge}; m: {self.m}"

    def _check_valid_style(self) -> bool:
        """
        check whether the adduct string is valid
        :return: True for valid, False for invalid
        """
        if self.string.count(']') == 1 and self.string.count('[') == 1 and 'M' in self.string:
            if ("+" in self.string[len(self.string) - 1] and self.pos_mode) \
                    or ("-" in self.string[len(self.string) - 1] and not self.pos_mode):
                return True
        return False

    # return true for valid characters in the string (except the common abbr.)
    def _check_valid_character(self) -> bool:
        valid_elements = {"-", "+", "C", "H", "Br", "Cl", "F", "I", "K", "N", "Na", "O", "P", "S"}

        m_index = self.string.index('M')
        right_index = self.string.index(']')
        segment = self.string[m_index + 1: right_index]

        i = 0
        while i < len(segment):
            # Check for two-character elements
            if i < len(segment) - 1 and segment[i:i + 2] in valid_elements:
                i += 2
            # Check for one-character elements and digits
            elif segment[i] in valid_elements or segment[i].isdigit():
                i += 1
            else:
                return False
        return True

    # return default value for invalid adduct
    def _invalid(self, string, report_invalid: bool = False):
        if report_invalid:
            raise ValueError("Adduct cannot be parsed: " + string)
        else:
            if self.pos_mode:
                logging.warning("Adduct cannot be parsed: " + string + ", set to [M+H]+")
                self.string = "[M+H]+"
                self.charge = +1
                self.loss_formula = None
                self.net_formula = Formula(array=[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], charge=0)
            else:
                logging.warning("Adduct cannot be parsed: " + string + ", set to [M-H]-")
                self.string = "[M-H]-"
                self.charge = -1
                self.loss_formula = Formula(array=[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], charge=0)
                self.net_formula = Formula(array=[0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], charge=0)
            self.m = 1

    def _check_common(self) -> bool:
        if self.pos_mode and self.string in ['[M+H]+', '[M+NH4]+', '[M+Na]+', '[M+K]+', '[M+H-H2O]+', '[M-H2O+H]+']:
            self.charge = +1
            self.m = 1
            if self.string == '[M+H]+':
                self.loss_formula = None
                self.net_formula = Formula(array=[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], charge=0)
            elif self.string == '[M+NH4]+':
                self.loss_formula = None
                self.net_formula = Formula(array=[0, 4, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], charge=0)
            elif self.string == '[M+Na]+':
                self.loss_formula = None
                self.net_formula = Formula(array=[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], charge=0)
            elif self.string == '[M+K]+':
                self.loss_formula = None
                self.net_formula = Formula(array=[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], charge=0)
            else:
                self.loss_formula = Formula(array=[0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], charge=0)
                self.net_formula = Formula(array=[0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0], charge=0)
            return True
        elif (not self.pos_mode) and self.string in ['[M-H]-', '[M+Cl]-', '[M+Br]-', '[M-H2O-H]-', '[M-H-H2O]-']:
            self.charge = -1
            self.m = 1
            if self.string == '[M-H]-':
                self.loss_formula = Formula(array=[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], charge=0)
                self.net_formula = Formula(array=[0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], charge=0)
            elif self.string == '[M+Cl]-':
                self.loss_formula = None
                self.net_formula = Formula(array=[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], charge=0)
            elif self.string == '[M+Br]-':
                self.loss_formula = None
                self.net_formula = Formula(array=[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], charge=0)
            else:
                self.loss_formula = Formula(array=[0, 3, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], charge=0)
                self.net_formula = Formula(array=[0, -3, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0], charge=0)
            return True
        return False

    def _calc_charge(self):
        if self.pos_mode:
            _substr = self.string.split(']')[1].replace('+', '')
            if _substr == '':
                self.charge = +1
            else:
                self.charge = int(_substr)
        else:
            _substr = self.string.split(']')[1].replace('-', '')
            if _substr == '':
                self.charge = -1
            else:
                self.charge = -int(_substr)

    def _calc_m(self):
        idx = self.string.index('M')
        _substr = self.string[idx - 1]
        if _substr.isnumeric():
            self.m = int(_substr)
        else:
            self.m = 1

    # exchange adduct name
    def _exchange(self):
        self.string = self.string.replace("ACN", "C2H3N")
        self.string = self.string.replace("IsoProp", "C3H8O")
        self.string = self.string.replace("DMSO", "C2H6OS")
        self.string = self.string.replace("FA", "HCOOH")
        self.string = self.string.replace("HFA", "HCOOH")
        self.string = self.string.replace("Hac", "CH3COOH")
        self.string = self.string.replace("HAc", "CH3COOH")
        self.string = self.string.replace("HAC", "CH3COOH")
        self.string = self.string.replace("Ac", "CH3COO")
        self.string = self.string.replace("ac", "CH3COO")
        self.string = self.string.replace("AC", "CH3COO")
        self.string = self.string.replace("TFA", "CF3COOH")
        self.string = self.string.replace("MeOH", "CH4O")

    def _calc_loss_and_net_formula(self):
        i = 1
        add_index = []
        loss_index = []
        add = ''
        loss = ''
        repeat_add = ''
        repeat_loss = ''

        # get the starting index of add and loss elements
        for element in self.string[1: len(self.string) - 2]:
            if element == "+":
                add_index.append(i)
            elif element == "-":
                loss_index.append(i)
            else:
                pass
            i += 1

        # get the losing and adding elements
        while len(add_index) != 0 and len(loss_index) != 0:
            if loss_index[0] < add_index[0]:
                if len(loss_index) == 1:
                    loss = loss + self.string[loss_index[0]:add_index[0]]
                else:
                    if loss_index[1] > add_index[0]:
                        loss = loss + self.string[loss_index[0]:add_index[0]]
                    else:
                        loss = loss + self.string[loss_index[0]:loss_index[1]]
                loss_index.remove(loss_index[0])
            else:
                if len(add_index) == 1:
                    add = add + self.string[add_index[0]:loss_index[0]]
                else:
                    if add_index[1] < loss_index[0]:
                        add = add + self.string[add_index[0]:add_index[1]]
                    else:
                        add = add + self.string[add_index[0]:loss_index[0]]
                add_index.remove(add_index[0])
        if len(add_index) == 0 and len(loss_index) == 0:
            self._invalid(self.string)
            return
        elif len(add_index) == 0:
            right_index = self.string.index(']')
            loss = loss + self.string[loss_index[0]: right_index]
        else:
            right_index = self.string.index(']')
            add = add + self.string[add_index[0]: right_index]

        # deal with situation like "2H2O"
        for index, element in enumerate(loss):
            if element == "-":
                if loss[len(loss) - 1] == '-' or loss[index + 1] == '-':
                    self._invalid(self.string)
                    return
                elif loss[index + 1].isnumeric():
                    repeat_element = ''
                    repeat_time = int(loss[index + 1])
                    for e in loss[index + 2:]:
                        if e != "-":
                            repeat_element = repeat_element + e
                        else:
                            break
                    if repeat_time > 2:
                        repeat_element = repeat_element * (repeat_time - 1)
                    list1 = list(loss)
                    list1[index + 1] = ' '
                    loss = (''.join(list1))
                    repeat_loss = repeat_loss + repeat_element
        if len(repeat_loss) != 0:
            loss = loss + '-' + repeat_loss
        loss = loss.replace(" ", "")

        for index2, element2 in enumerate(add):
            if element2 == '+':
                if add[len(add) - 1] == '+' or add[index2 + 1] == '+':
                    self._invalid(self.string)
                    return
                elif add[index2 + 1].isnumeric():
                    repeat_element2 = ''
                    repeat_time2 = int(add[index2 + 1])
                    for e2 in add[index2 + 2:]:
                        if e2 != "+":
                            repeat_element2 = repeat_element2 + e2
                        else:
                            break
                    if repeat_time2 > 2:
                        repeat_element2 = repeat_element2 * (repeat_time2 - 1)
                    list2 = list(add)
                    list2[index2 + 1] = ''
                    add = (''.join(list2))
                    repeat_add = repeat_add + repeat_element2
        if len(repeat_add) != 0:
            add = add + '+' + repeat_add
        add = add.replace(" ", "")

        loss_array = read_formula(loss.replace("-", ""))
        add_array = read_formula(add.replace("+", ""))
        net_array = np.subtract(add_array, loss_array)

        self.net_formula = Formula(array=net_array, charge=0)
        self.loss_formula = Formula(array=loss_array, charge=0)


def check_adduct(adduct_str: str) -> Tuple[bool, bool]:
    """
    given an adduct string, check whether it is valid, and determine its ion mode
    :param adduct_str: adduct string
    :return: (valid, pos_mode)
    """
    # check whether the adduct string is valid
    adduct_str = adduct_str.replace(" ", "")
    # check '[' and ']' and 'M'
    if adduct_str.count(']') != 1 or adduct_str.count('[') != 1 or 'M' not in adduct_str:
        return False, False
    # check sign after ']'
    if adduct_str[len(adduct_str) - 1] not in ['+', '-']:
        return False, False

    # ion mode
    pos_mode = True if adduct_str[len(adduct_str) - 1] == '+' else False
    return True, pos_mode


class ProcessedMS1:
    """
    Processed MS1 class
    :attr mz_tol: float
    :attr ppm: bool
    :attr idx_array: np.array, for retrieving raw indices of peaks, record raw indices of selected peaks e.g., [3,5,6]
    :attr mz_array: np.array
    :attr int_array: np.array
    """

    def __init__(self, mz: float, raw_spec: Spectrum, charge: int,
                 mz_tol: float, ppm: bool,
                 isotope_bin_mztol: float, max_isotope_cnt: int):
        if raw_spec:
            self.mz_tol = mz_tol
            self.ppm = ppm
            self._find_ms1_isotope(mz, raw_spec, charge, isotope_bin_mztol, max_isotope_cnt)
        else:
            self.idx_array = np.array([], dtype=int)
            self.mz_array = np.array([], dtype=np.float32)
            self.int_array = np.array([], dtype=np.float32)

    def __bool__(self):
        return self.mz_array.size != 0

    def __str__(self) -> str:
        return f"ProcessedMS1(idx={self.idx_array}, mz={self.mz_array}, int={self.int_array})"

    def __len__(self):
        return self.mz_array.size

    def _find_ms1_isotope(self, mz: float, raw_spec: Spectrum, charge: int,
                          isotope_bin_mztol: float, max_isotope_cnt: int):
        """
        find the MS1 isotope pattern from MS1 raw spectrum, record all peaks (could have multiple peaks in one bin)
        :param mz: precursor mz
        :param raw_spec: raw ms1 spectrum
        :param charge: precursor charge
        :param isotope_bin_mztol: isotope bin m/z tolerance
        :param max_isotope_cnt: maximum isotope count
        :return: idx array, mz array, int array
        """

        tmp_mz_diff = mz * self.mz_tol * 1e-6 if self.ppm else self.mz_tol

        # find the closest peak to the precursor mz
        m0_found, idx = _find_m0(mz, raw_spec.mz_array, tmp_mz_diff)
        if not m0_found:
            self.idx_array = np.array([])
            self.mz_array = np.array([], dtype=np.float32)
            self.int_array = np.array([], dtype=np.float32)
            return

        # fill in M0
        self.idx_array = np.array([idx])
        self.mz_array = np.array([raw_spec.mz_array[idx]])
        self.int_array = np.array([raw_spec.int_array[idx]])

        # find other isotope peaks
        idx_arr = _find_iso_peaks(mz, raw_spec.mz_array, charge, isotope_bin_mztol, max_isotope_cnt)

        self.idx_array = np.concatenate((self.idx_array, idx_arr))
        self.mz_array = np.concatenate((self.mz_array, raw_spec.mz_array[idx_arr]))
        self.int_array = np.concatenate((self.int_array, raw_spec.int_array[idx_arr]))


def _find_m0(mz: float, mz_array: np.array, mz_diff: float) -> Tuple[bool, int]:
    """
    find ms1 isotope M0 peak (for class ProcessedMS1)
    :param mz: float
    :param mz_array: np.array
    :param mz_diff: float
    :return: bool (M0 found) & index of the closest peak to the precursor mz
    """
    m0_found = False
    idx = 0
    tmp_mz_diff = mz_diff
    for k, m in enumerate(mz_array):
        if abs(m - mz) <= tmp_mz_diff:
            m0_found = True
            idx = k
            tmp_mz_diff = abs(m - mz)
    return m0_found, idx


def _find_iso_peaks(mz: float, mz_arr: np.array, charge: int,
                    isotope_bin_mztol: float, max_isotope_cnt: int) -> np.array:
    """
    find ms1 isotope M1, M2, ... peaks (for class ProcessedMS1)
    :param mz: float
    :param mz_arr: np.array
    :param charge: int
    :param isotope_bin_mztol: float
    :param max_isotope_cnt: int
    :return: idx array, mz array, int array
    """
    idx_array = np.array([], dtype=int)
    tmp_mz = mz
    for k in range(max_isotope_cnt - 1):
        tmp_mz += mass_i / abs(charge)
        found = False
        for j, m in enumerate(mz_arr):
            if abs(m - tmp_mz) <= isotope_bin_mztol:
                found = True
                idx_array = np.append(idx_array, j)
        if not found:
            break
    return idx_array


class ProcessedMS2:
    """
    Processed MS2 class
    :attr mz_tol: float
    :attr ppm: bool
    :attr idx_array: np.array, for retrieving raw indices of peaks, record raw indices of selected peaks e.g., [3,5,6]
    :attr mz_array: np.array
    :attr int_array: np.array
    """

    def __init__(self, mz: float, raw_spec: Spectrum,
                 mz_tol: float, ppm: bool,
                 rel_int_denoise_cutoff: float,
                 top_n_per_50_da: int):
        if raw_spec:
            self.mz_tol = mz_tol
            self.ppm = ppm
            self._preprocess(mz, raw_spec, rel_int_denoise_cutoff, top_n_per_50_da)
        else:
            self.idx_array = np.array([], dtype=int)
            self.mz_array = np.array([], dtype=np.float32)
            self.int_array = np.array([], dtype=np.float32)

    def __bool__(self):
        return self.mz_array.size != 0

    def __str__(self) -> str:
        return f"ProcessedMS2(idx={self.idx_array}, mz={self.mz_array}, int={self.int_array})"

    def __len__(self):
        return len(self.mz_array)

    def _preprocess(self, mz: float, raw_spec: Spectrum,
                    rel_int_denoise_cutoff: float, top_n_per_50_da: int):
        """
        preprocess MS2 spectrum, denoise (optional), de-precursor
        :param mz: precursor mz
        :param raw_spec: raw ms2 spectrum
        :param rel_int_denoise_cutoff: relative intensity cutoff
        :param top_n_per_50_da: max fragment count reserved for every 50 Da
        :return: fill self.idx_array, self.mz_array, self.int_array
        """

        # de-precursor
        if len(raw_spec.mz_array) > 0:
            self._deprecursor(mz, raw_spec)

        # denoise
        if len(self.mz_array) > 0:
            self._denoise(rel_int_denoise_cutoff)

        # keep top_n_per_50_da peaks in each 50 Da
        if len(self.mz_array) > 0:
            self._keep_top_n_per_50_da(top_n_per_50_da)

        # top n fragment
        if len(self.mz_array) > 0:
            top_n_frag = _calc_top_n_frag(mz)
            if len(self.mz_array) > top_n_frag:
                idx = np.argsort(self.int_array)
                reserved_idx = self.int_array >= self.int_array[idx[-top_n_frag]]
                self.idx_array = self.idx_array[reserved_idx]
                self.mz_array = self.mz_array[reserved_idx]
                self.int_array = self.int_array[reserved_idx]

    def _denoise(self, rel_int_denoise_cutoff: float):
        """
        denoise MS2 spectrum
        :param rel_int_denoise_cutoff: relative intensity cutoff
        :return: fill self.idx_array, self.mz_array, self.int_array
        """

        # sort by intensity, increasing
        idx = np.argsort(self.int_array)
        sorted_int = self.int_array[idx]

        final_int_threshold = rel_int_denoise_cutoff * sorted_int[-1]

        # refill self.idx_array, self.mz_array, self.int_array
        reserved_idx = self.int_array >= final_int_threshold
        self.idx_array = self.idx_array[reserved_idx]
        self.mz_array = self.mz_array[reserved_idx]
        self.int_array = self.int_array[reserved_idx]

    def _deprecursor(self, mz: float, raw_spec: Spectrum):
        """
        de-precursor
        :param mz: precursor mz
        :return: fill self.idx_array, self.mz_array, self.int_array
        """

        # remove peaks >= precursor mz - 1.5 Da
        mz_tol = mz * self.mz_tol * 1e-6 if self.ppm else self.mz_tol

        idx = raw_spec.mz_array < (mz - 1.5 - mz_tol)

        self.idx_array = np.arange(len(raw_spec.mz_array))[idx]
        self.mz_array = raw_spec.mz_array[idx]
        self.int_array = raw_spec.int_array[idx]

    def _keep_top_n_per_50_da(self, top_n_per_50_da: int):
        """
        keep top_n_per_50_da peaks in each 50 Da
        :param top_n_per_50_da: max fragment count reserved for every 50 Da
        :return: fill self.idx_array, self.mz_array, self.int_array
        """
        if top_n_per_50_da <= 1:
            return

        bool_arr = ms2_denoise(mz_arr=self.mz_array, int_arr=self.int_array, top_n=top_n_per_50_da)
        self.idx_array = self.idx_array[bool_arr]
        self.mz_array = self.mz_array[bool_arr]
        self.int_array = self.int_array[bool_arr]

    def normalize_intensity(self, method: str = 'sum'):
        """
        normalize intensity
        :param method: 'sum' or 'max'
        :return: fill self.int_array
        """
        if method == 'sum':
            self.int_array = self.int_array / np.sum(self.int_array)
        elif method == 'max':
            self.int_array = self.int_array / np.max(self.int_array)


@njit
def ms2_denoise(mz_arr, int_arr, top_n=6, per_mass_range=50):
    """
    keep top_n peaks in each per_mass_range
    """
    group_arr = mz_arr // per_mass_range
    bool_arr = np.array([True] * len(int_arr))

    for group in np.unique(group_arr):
        group_idx = np.where(group_arr == group)[0]
        if len(group_idx) > top_n:
            int_threshold = np.sort(int_arr[group_idx])[-top_n]
            bool_arr[group_idx] = int_arr[group_idx] >= int_threshold

    return bool_arr


def _calc_top_n_frag(pre_mz: float) -> int:
    """
    calculate top n frag No., a linear function of precursor m/z (for class ProcessedMS2)
    :param pre_mz: precursor m/z
    :param max_frag_reserved: max fragment count reserved
    :return: top n frag No. (int)
    """
    top_n = int(20 + 0.05 * pre_mz)
    return top_n


class MS2Explanation:
    """
    MS2Explanation class, used for storing MS2 explanation.
    """

    def __init__(self, idx_array: np.array,
                 explanation_list: List[Union[Formula, None]]):
        self.idx_array = idx_array  # indices of peaks in MS2 spectrum
        self.explanation_list = explanation_list  # List[Formula], isotope peaks are included

    def __str__(self):
        out_str = ""
        for i in range(len(self.idx_array)):
            out_str += "idx: {}, formula: {}\n".format(self.idx_array[i], str(self.explanation_list[i]))
        return out_str

    def __len__(self):
        return len(self.idx_array)

    def __bool__(self):
        return len(self) > 0


class CandidateFormula:
    """
    CandidateFormula is a class for storing a candidate formula. It's used in MetaFeature.candidate_formula_list.
    precursor formula in CandidateFormula is a neutral formula
    """

    def __init__(self, formula: Formula,
                 charged_formula: Union[Formula, None] = None,
                 mz_error: Union[float, None] = None,
                 exp_ms2_sum_int: Union[float, None] = None,
                 ms2_raw_explanation: Union[MS2Explanation, None] = None,
                 db_existed: bool = False,
                 optimal_formula: bool = False,
                 ms2_refined_explanation: Union[MS2Explanation, None] = None):
        self.formula = formula  # neutral formula
        self.charged_formula = charged_formula  # charged formula
        self.mz_error = mz_error  # mz error in ppm or Da, depending on config
        self.estimated_prob = None  # estimated probability (ML score, not normalized)
        self.normed_estimated_prob = None  # normalized estimated probability considering all candidate formulas
        self.estimated_fdr = None  # estimated FDR
        self.formula_feature_array = None  # formula feature array
        self.ms1_isotope_similarity = None
        self.exp_ms2_sum_int = exp_ms2_sum_int  # sum of MS2 intensity of peaks explained (during bottom-up search)
        self.ms2_raw_explanation = ms2_raw_explanation  # ms2 explanation during subformula assignment
        self.db_existed = db_existed  # whether this formula is in the formula database
        # self.optimal_formula = optimal_formula
        # self.ms2_refined_explanation = ms2_refined_explanation  # re-annotate frags using global optim.

    def __str__(self):
        cf_str = form_arr_to_str(self.formula.array) + \
                 "  est_prob: " + str(self.estimated_prob) + "  est_fdr: " + str(self.estimated_fdr)
        if self.ms2_raw_explanation:
            cf_str += " ms2_raw_exp: " + str(len(self.ms2_raw_explanation))

        return cf_str


class MetaFeature:
    """
    MetaFeature class, used for storing a metabolic feature.
    """

    def __init__(self,
                 identifier: Union[str, int],
                 mz: float,
                 charge: int,
                 rt: Union[float, None] = None,  # retention time in seconds
                 adduct: Union[str, None] = None,
                 ms1: Union[Spectrum, None] = None,
                 ms2: Union[Spectrum, None] = None):
        """
        Instantiate a new MetaFeature object.
        :param identifier: identifier (str or int)
        :param mz: precursor m/z (float)
        :param charge: precursor charge (int)
        :param rt: retention time in seconds (float)
        :param adduct: adduct string, e.g., "[M+H]+" (str)
        :param ms1: MS1 spectrum (Spectrum class)
        :param ms2: MS2 spectrum (Spectrum class)
        """
        if charge == 0:
            raise ValueError("Charge cannot be 0.")
        self.identifier = identifier
        self.mz = mz
        self.rt = rt

        pos_mode = charge > 0
        self.adduct = Adduct(adduct, pos_mode)  # type: Adduct

        self.ms1_raw = ms1
        self.ms1_processed = None  # type: ProcessedMS1 or None
        self.ms2_raw = ms2
        self.ms2_processed = None  # type: ProcessedMS2 or None
        self.candidate_formula_list = None  # Union[List[CandidateFormula], None]

    def __str__(self):
        mf_str = "mz " + str(self.mz) + "  adduct " + self.adduct.string
        if self.candidate_formula_list:
            mf_str += "  cf count: " + str(len(self.candidate_formula_list))
        else:
            mf_str += "  cf count: 0"
        return mf_str

    def data_preprocess(self, ppm: bool, ms1_tol: float, ms2_tol: float,
                        isotope_bin_mztol: float, max_isotope_cnt: int,
                        rel_int_denoise_cutoff: float,
                        top_n_per_50_da: int):
        """
        Data preprocessing.
        :param ppm: whether to use ppm as m/z tolerance
        :param ms1_tol: m/z tolerance for MS1, used for MS1 isotope pattern
        :param ms2_tol: m/z tolerance for MS2, used for MS2 data
        :param isotope_bin_mztol: m/z tolerance for isotope bin, used for MS1 isotope pattern
        :param max_isotope_cnt: maximum isotope count, used for MS1 isotope pattern
        :param rel_int_denoise_cutoff: relative intensity cutoff for MS2 denoise
        :param top_n_per_50_da: max fragment count reserved for every 50 Da
        :return: fill in ms1_processed and ms2_processed for each metaFeature
        """
        if self.ms1_raw:
            self.ms1_processed = ProcessedMS1(self.mz, self.ms1_raw,
                                              self.adduct.charge,
                                              ms1_tol, ppm, isotope_bin_mztol, max_isotope_cnt)
        if self.ms2_raw:
            self.ms2_processed = ProcessedMS2(self.mz, self.ms2_raw,
                                              ms2_tol, ppm, rel_int_denoise_cutoff, top_n_per_50_da)

    def summarize_result(self) -> dict:
        """
        Summarize the annotation result for a MetaFeature.
        :return: dict
        """
        result = {'identifier': self.identifier, 'mz': self.mz, 'rt': self.rt, 'adduct': self.adduct.string,
                  'formula_rank_1': None, 'estimated_fdr': None, 'formula_rank_2': None, 'formula_rank_3': None,
                  'formula_rank_4': None, 'formula_rank_5': None}

        if self.candidate_formula_list:
            result['formula_rank_1'] = form_arr_to_str(self.candidate_formula_list[0].formula.array)
            result['estimated_fdr'] = self.candidate_formula_list[0].estimated_fdr
            if len(self.candidate_formula_list) > 1:
                result['formula_rank_2'] = form_arr_to_str(self.candidate_formula_list[1].formula.array)
            if len(self.candidate_formula_list) > 2:
                result['formula_rank_3'] = form_arr_to_str(self.candidate_formula_list[2].formula.array)
            if len(self.candidate_formula_list) > 3:
                result['formula_rank_4'] = form_arr_to_str(self.candidate_formula_list[3].formula.array)
            if len(self.candidate_formula_list) > 4:
                result['formula_rank_5'] = form_arr_to_str(self.candidate_formula_list[4].formula.array)
        return result
