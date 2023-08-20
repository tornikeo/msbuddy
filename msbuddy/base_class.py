import sys
import math
import numpy as np
from typing import Union, Tuple, List
from chemparse import parse_formula


class Formula:
    alphabet = ["C", "H", "Br", "Cl", "F", "I", "K", "N", "Na", "O", "P", "S"]
    mass_arr = np.array([12.000000, 1.007825, 78.918336, 34.968853, 18.998403, 126.904473, 38.963707, 14.003074,
                         22.989769, 15.994915, 30.973762, 31.972071])

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
        self.array = array
        self.charge = charge
        self.__calc_formula_dbe()
        self.isotope = isotope

        # fill in mass directly from formula database, otherwise calculate
        if mass is None:
            self.__calc_formula_mass()
        else:
            self.mass = mass

    def __bool__(self):
        # if all elements are 0, then return False
        return bool(np.any(self.array))

    def __str__(self) -> str:
        if not self:
            return 'Null'

        return form_arr_to_str(self.array) + "  charge:" + str(self.charge) + "  mass:" + str(self.mass) + \
            "  dbe:" + str(self.dbe) + "  isotope:" + str(self.isotope)

    def __calc_formula_dbe(self):
        """
        calculate dbe of a formula
        :return: fill in dbe
        """
        # ["C", "H", "Br", "Cl", "F", "I", "K", "N", "Na", "O", "P", "S"]
        # DBE: C + 1 - (H + F + Cl + Br + I + Na + K)/2 + (N + P)/2
        arr = self.array
        dbe = arr[0] + 1 - (arr[1] + arr[4] + arr[3] + arr[2] + arr[5] + arr[8] + arr[6]) / 2 + (arr[7] + arr[10]) / 2

        self.dbe = dbe

    def __calc_formula_mass(self):
        """
        calculate monoisotopic mass of a formula, charge adjusted
        :return: fill in mass
        """
        if self.charge == 0:
            self.mass = float(np.sum(self.array * Formula.mass_arr) + self.isotope * 1.0033548)
        else:
            self.mass = float((np.sum(self.array * Formula.mass_arr) - self.charge * 0.00054858 +
                               self.isotope * 1.0033548) / abs(self.charge))


# global functions related to formula
def read_formula(form: str) -> np.array:
    """
    read neutral formula string and return a 12-dim array, return None if invalid
    :param form: formula string in neutral form, e.g., C6H5NO2
    :return: 12-dim array
    """
    # alphabet = ["C", "H", "Br", "Cl", "F", "I", "K", "N", "Na", "O", "P", "S"]

    # parse formula string
    try:
        parsed = parse_formula(form)
    except:
        return None

    # check whether contains elements not in alphabet
    for element in parsed.keys():
        if element not in Formula.alphabet:
            return None

    # convert to 12-dim array
    array = np.zeros(12, dtype=int)
    for i, element in enumerate(Formula.alphabet):
        if element in parsed.keys():
            array[i] = parsed[element]
    return array


def form_arr_to_str(form_arr: np.array) -> str:
    """
    convert 12-dim array to formula string
    :param form_arr: 12-dim array
    :return: formula string
    """

    def decode(element: str, cnt: int) -> str:
        if not cnt:
            return ''
        if cnt == 1:
            return element
        return element + str(cnt)

    formula_str = ''.join([decode(s, c) for s, c in zip(Formula.alphabet, form_arr)])
    return formula_str


class Spectrum:
    def __init__(self, mz_array: np.array, int_array: np.array):
        """
        :param mz_array: np.array
        :param int_array: np.array
        """

        # mz_array and int_array must have the same length
        if len(mz_array) != len(int_array):
            raise ValueError("mz_array and int_array must have the same length")

        # sort by mz_array
        idx = np.argsort(mz_array)
        self.mz_array = mz_array[idx]
        self.int_array = int_array[idx]

    def __str__(self) -> str:
        return f"Spectrum(mz={self.mz_array}, int={self.int_array})"

    def __bool__(self):
        return self.mz_array.size != 0

    def __len__(self):
        return self.mz_array.size


class Adduct:
    # string: str
    # pos_mode: bool
    # charge: int
    # m: int
    # net_formula: Formula
    # loss_formula: Formula  # the sum loss formula within the adduct
    # # for example, [M+H-H2O]+ has loss formula H2O; [M+Na]+ has no loss formula;
    # [M - 2H2O - Cl] - has loss formula H4O2Cl

    def __init__(self,
                 string: Union[str, None],
                 pos_mode: bool):
        """
        :param string: adduct string, e.g. [M+H]+
        :param pos_mode: True for positive mode, False for negative mode
        """
        if string is None:
            self.string = "[M+H]+" if pos_mode else "[M-H]-"
        else:
            self.string = string.replace(" ", "")

        self.pos_mode = pos_mode

        if self._check_valid_style():
            self._exchange()
            if self._check_common():
                pass
            else:
                if self._check_valid_character():
                    self._calc_m()
                    self._calc_charge()
                    self._calc_loss_and_net_formula()
                else:
                    self._invalid()
        else:
            self._invalid()

    def __str__(self):
        return f"{self.string} has charge {self.charge}, and m is {self.m}"

    # return true for valid input style (have M, [], -/+ at last, and at least one add/loss item)
    def _check_valid_style(self) -> bool:
        if self.string.count(']') == 1 and self.string.count('[') == 1 and 'M' in self.string:
            if ("+" in self.string[len(self.string) - 1] and self.pos_mode) \
                    or ("-" in self.string[len(self.string) - 1] and not self.pos_mode):
                return True
        return False

    # return true for valid characters in the string (except the common abbr.)
    def _check_valid_character(self) -> bool:
        valid_character = ["-", "+", "C", "H", "Br", "Cl", "F", "I", "K", "N", "Na", "O", "P", "S",
                           "1", "2", "3", "4", "5", "6", "7", "8", "9", "0"]
        m_index = self.string.index('M')
        right_index = self.string.index(']')
        return all([valid in valid_character for valid in self.string[m_index + 1: right_index]])

    # return default value for invalid adduct
    def _invalid(self):
        if self.pos_mode:
            sys.stdout.write("Invalid adduct for positive mode, set to [M+H]+\n")
            self.string = "[M+H]+"
            self.charge = +1
            self.loss_formula = None
            self.net_formula = Formula(array=[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], charge=0)
        else:
            sys.stdout.write("Invalid adduct for negative mode, set to [M-H]-\n")
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
            if self.string[len(self.string) - 2] == "]":
                self.charge = +1
            else:
                try:
                    self.charge = int(self.string[len(self.string) - 2])
                except ValueError:
                    self.charge = +1
        else:
            if self.string[len(self.string) - 2] == "]":
                self.charge = -1
            else:
                try:
                    self.charge = -int(self.string[len(self.string) - 2])
                except ValueError:
                    self.charge = -1

    def _calc_m(self):
        count = 0
        for element in self.string:
            if element != "M":
                count = count + 1
            else:
                break
        if self.string[count - 1].isnumeric():
            self.m = self.string[count - 1]
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
            i = i + 1

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
            self._invalid()
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
                    self._invalid()
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
                    self._invalid()
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


class ProcessedMS1:
    # MS1 preprocessing, find MS1 isotope pattern
    # mz_tol: float
    # ppm: bool
    # idx_array: np.array # for retrieving raw indices of peaks, record raw indices of selected peaks e.g., [3,5,6]
    # mz_array: np.array
    # int_array: np.array

    def __init__(self, mz: float, raw_spec: Spectrum, charge: int,
                 mz_tol: float, ppm: bool,
                 isotope_bin_mztol: float, max_isotope_cnt: int):
        if raw_spec:
            self.mz_tol = mz_tol
            self.ppm = ppm
            self._find_ms1_isotope(mz, raw_spec, charge, isotope_bin_mztol, max_isotope_cnt)
        else:
            self.idx_array = np.array([])
            self.mz_array = np.array([])
            self.int_array = np.array([])

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
        m0_found, idx = _func_a(mz, raw_spec.mz_array, tmp_mz_diff)
        if not m0_found:
            return

        # fill in M0
        self.idx_array = np.array([idx])
        self.mz_array = np.array([raw_spec.mz_array[idx]])
        self.int_array = np.array([raw_spec.int_array[idx]])

        # find other isotope peaks
        idx_arr = _func_b(mz, raw_spec.mz_array, charge, isotope_bin_mztol, max_isotope_cnt)

        self.idx_array = np.concatenate((self.idx_array, idx_arr))
        self.mz_array = np.concatenate((self.mz_array, raw_spec.mz_array[idx_arr]))
        self.int_array = np.concatenate((self.int_array, raw_spec.int_array[idx_arr]))


# @njit
def _func_a(mz: float, mz_array: np.array, mz_diff: float) -> Tuple[bool, int]:
    """
    helper function a for find ms1 isotope (for class ProcessedMS1)
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


# @njit
def _func_b(mz: float, mz_arr: np.array, charge: int,
            isotope_bin_mztol: float, max_isotope_cnt: int) -> np.array:
    """
    helper function b for find ms1 isotope (for class ProcessedMS1)
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
        tmp_mz += 1.003355 / abs(charge)
        found = False
        for j, m in enumerate(mz_arr):
            if abs(m - tmp_mz) <= isotope_bin_mztol:
                found = True
                idx_array = np.append(idx_array, j)
        if not found:
            break
    return idx_array


class ProcessedMS2:
    # MS2 preprocessing, de-precursor, denoise, top n fragments
    # mz_tol: float
    # ppm: bool
    # idx_array: np.array  # for retrieving raw indices of peaks
    # mz_array: np.array
    # int_array: np.array

    def __init__(self, mz: float, raw_spec: Spectrum,
                 mz_tol: float, ppm: bool,
                 denoise: bool,
                 rel_int_denoise: bool, rel_int_denoise_cutoff: float,
                 max_noise_frag_ratio: float, max_noise_rsd: float,
                 max_frag_reserved: int, use_all_frag: bool = False):
        if raw_spec:
            self.mz_tol = mz_tol
            self.ppm = ppm
            self._preprocess(mz, raw_spec, denoise, rel_int_denoise, rel_int_denoise_cutoff,
                             max_noise_frag_ratio, max_noise_rsd, max_frag_reserved, use_all_frag)
        else:
            self.idx_array = np.array([])
            self.mz_array = np.array([])
            self.int_array = np.array([])

    def __bool__(self):
        return self.mz_array.size != 0

    def __str__(self) -> str:
        return f"ProcessedMS2(idx={self.idx_array}, mz={self.mz_array}, int={self.int_array})"

    def __len__(self):
        return len(self.mz_array)

    def _preprocess(self, mz: float, raw_spec: Spectrum,
                    denoise: bool, rel_int_denoise: bool, rel_int_denoise_cutoff: float,
                    max_noise_frag_ratio: float, max_noise_rsd: float, max_frag_reserved: int,
                    use_all_frag: bool):
        """
        preprocess MS2 spectrum, denoise (optional), de-precursor
        :param mz: precursor mz
        :param raw_spec: raw ms2 spectrum
        :param denoise: whether to denoise
        :param max_noise_frag_ratio: max noise fragment ratio
        :param max_noise_rsd: max noise RSD
        :param rel_int_denoise: whether to use relative intensity to denoise
        :param rel_int_denoise_cutoff: relative intensity cutoff
        :param max_frag_reserved: max fragment count reserved
        :return: fill self.idx_array, self.mz_array, self.int_array
        """

        # de-precursor
        self._deprecursor(mz, raw_spec)

        # denoise, only perform when >= 10 peaks left
        if denoise and len(self.mz_array) >= 10:
            self._denoise(rel_int_denoise, rel_int_denoise_cutoff, max_noise_frag_ratio, max_noise_rsd)

        # top n frag
        if not use_all_frag:
            top_n_frag = _calc_top_n_frag(mz, max_frag_reserved)
            if top_n_frag < len(self.mz_array):
                idx = np.argsort(self.int_array)
                reserved_idx = self.int_array >= self.int_array[idx[-top_n_frag]]
                self.idx_array = self.idx_array[reserved_idx]
                self.mz_array = self.mz_array[reserved_idx]
                self.int_array = self.int_array[reserved_idx]

    def _denoise(self, rel_int_denoise: bool, rel_int_denoise_cutoff: float,
                 max_noise_frag_ratio: float, max_noise_rsd: float):
        """
        denoise MS2 spectrum
        :param rel_int_denoise: whether to use relative intensity to denoise
        :param rel_int_denoise_cutoff: relative intensity cutoff
        :param max_noise_frag_ratio: max noise fragment ratio
        :param max_noise_rsd: max noise RSD
        :return: fill self.idx_array, self.mz_array, self.int_array
        """

        # sort by intensity, increasing
        idx = np.argsort(self.int_array)
        # sorted_mz = self.mz_array[idx]
        sorted_int = self.int_array[idx]

        if rel_int_denoise:
            final_int_threshold = rel_int_denoise_cutoff * sorted_int[-1]
        else:
            # at least 3 peaks are used to estimate RSD, step 0.05, round to 0.05 to determine m_start
            m_start = round(3.0 / len(self.mz_array) * 20.0)
            m_end = math.floor(max_noise_frag_ratio / 0.05) + 1
            if m_end <= m_start:
                return

            final_int_threshold = 0.0
            max_int_threshold = sorted_int[round(max_noise_frag_ratio * len(self.mz_array)) - 1]
            for m in range(m_start, m_end):
                if round(m * 0.05 * len(self.mz_array)) < 3:
                    continue
                sub_sorted_int = sorted_int[0:round(m * 0.05 * len(self.mz_array))]
                noise_mean = np.mean(sub_sorted_int)
                noise_sd = np.std(sub_sorted_int, ddof=1)
                noise_rsd = noise_sd / noise_mean
                if noise_rsd <= max_noise_rsd:
                    tmp_int_threshold = noise_mean + noise_sd * 3.0
                    if tmp_int_threshold > max_int_threshold:
                        break
                    else:
                        final_int_threshold = tmp_int_threshold

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


def _calc_top_n_frag(pre_mz: float, max_frag_reserved: int) -> int:
    """
    calculate top n frag No., a linear function of precursor m/z (for class ProcessedMS2)
    :param pre_mz: precursor m/z
    :param max_frag_reserved: max fragment count reserved
    :return: top n frag No. (int)
    """
    # if pre_mz < 1000:
    #     top_n = int(80 - 0.06 * pre_mz)
    # else:
    #     top_n = 20

    if pre_mz < 1000:
        top_n = int(90 - 0.08 * pre_mz)
    else:
        top_n = int(30 - 0.02 * pre_mz)
    return min(top_n, max_frag_reserved)


class MS2Explanation:
    # class for storing MS2 explanation
    def __init__(self, idx_array: np.array,
                 explanation_array: List[Union[Formula, None]]):

        self.idx_array = idx_array  # raw MS2 peak index, fragment index
        self.explanation_array = explanation_array  # List[Formula], isotope peaks are included

    def __str__(self):
        out_str = ""
        for i in range(len(self.idx_array)):
            out_str += "idx: {}, formula: {}\n".format(self.idx_array[i], str(self.explanation_array[i]))
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
                 optimal_formula: bool = False,
                 ms1_isotope_similarity: Union[float, None] = None,
                 ms2_raw_explanation: Union[MS2Explanation, None] = None,
                 ms2_refined_explanation: Union[MS2Explanation, None] = None):

        self.formula = formula  # neutral formula
        self.ml_a_prob = None  # ml_a score for model A (formula feasibility)
        self.estimated_prob = None  # estimated probability (ml_b score for model B)
        self.normed_estimated_prob = None  # normalized estimated probability considering all candidate formulas
        self.estimated_fdr = None  # estimated FDR
        self.optimal_formula = optimal_formula
        self.ms1_isotope_similarity = ms1_isotope_similarity
        self.ms2_raw_explanation = ms2_raw_explanation  # ms2 explanation during precursor formula annotation
        self.ms2_refined_explanation = ms2_refined_explanation  # re-annotate frags using global optim.
        # self.ml_array = None # ml_a feature array for model B (formula candidate ranking)

    def __str__(self):
        return f'{self.formula.__str__()}' + "  ml_a_prob: " + str(self.ml_a_prob) + \
            "  estimated_prob: " + str(self.estimated_prob) + "  estimated_fdr: " + str(self.estimated_fdr)


class MetaFeature:
    def __init__(self,
                 identifier: Union[str, int],
                 mz: float,
                 charge: int,
                 rt: Union[float, None] = None,  # retention time in seconds
                 adduct: Union[str, None] = None,
                 ms1: Union[Spectrum, None] = None,
                 ms2: Union[Spectrum, None] = None):
        """
        MetaFeature class.
        :param identifier: identifier of the metabolic feature (str or int)
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
        return "Meta feature: mz " + str(self.mz) + "  adduct: " + str(self.adduct)

    def data_preprocess(self, ppm: bool, ms1_tol: float, ms2_tol: float,
                        isotope_bin_mztol: float, max_isotope_cnt: int,
                        ms2_denoise: bool,
                        rel_int_denoise: bool, rel_int_denoise_cutoff: float,
                        max_noise_frag_ratio: float, max_noise_rsd: float,
                        max_frag_reserved: int, use_all_frag: bool):
        """
        Data preprocessing.
        :param ppm: whether to use ppm as m/z tolerance
        :param ms1_tol: m/z tolerance for MS1, used for MS1 isotope pattern
        :param ms2_tol: m/z tolerance for MS2, used for MS2 data
        :param isotope_bin_mztol: m/z tolerance for isotope bin, used for MS1 isotope pattern
        :param max_isotope_cnt: maximum isotope count, used for MS1 isotope pattern
        :param ms2_denoise: whether to denoise MS2 spectrum
        :param rel_int_denoise: whether to use relative intensity for MS2 denoise
        :param rel_int_denoise_cutoff: relative intensity cutoff for MS2 denoise
        :param max_noise_frag_ratio: maximum noise fragment ratio, used for MS2 denoise
        :param max_noise_rsd: maximum noise RSD, used for MS2 denoise
        :param max_frag_reserved: max fragment number reserved, used for MS2 data
        :param use_all_frag: whether to use all fragments for annotation
        :return: fill in ms1_processed and ms2_processed for each metaFeature
        """
        if self.ms1_raw:
            self.ms1_processed = ProcessedMS1(self.mz, self.ms1_raw,
                                              self.adduct.charge,
                                              ms1_tol, ppm, isotope_bin_mztol, max_isotope_cnt)
        if self.ms2_raw:
            self.ms2_processed = ProcessedMS2(self.mz, self.ms2_raw,
                                              ms2_tol, ppm, ms2_denoise,
                                              rel_int_denoise, rel_int_denoise_cutoff,
                                              max_noise_frag_ratio,
                                              max_noise_rsd, max_frag_reserved, use_all_frag)

    def summarize_result(self) -> dict:
        """
        Summarize the annotation result for a MetaFeature.
        :return: dict
        """
        result = {'identifier': self.identifier, 'mz': self.mz, 'rt': self.rt, 'adduct': self.adduct.string,
                  'formula_rank_1': None, 'estimated_fdr': None}
        if self.candidate_formula_list:
            result['formula_rank_1'] = form_arr_to_str(self.candidate_formula_list[0].formula.array)
            result['estimated_fdr'] = self.candidate_formula_list[0].estimated_fdr
        return result


# test
if __name__ == '__main__':
    a = Adduct(string='[M+2H]2+', pos_mode=True)
    print(a.loss_formula)
    print(a.net_formula)

    # b = Adduct(string="[M-2ACN-3H+2DMSO]-", pos_mode=False)
    # print(b.loss_formula)
    # print(b.net_formula)

    c = Adduct(string="[M-H]-", pos_mode=True)  # invalid adduct test
    print(c)
