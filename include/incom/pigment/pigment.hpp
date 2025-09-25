//  MIT License
//
//  Copyright (c) 2025 Mcchal Lurie
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.


// pigment.hpp
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace incom {
namespace pigment {
namespace detail {
inline constexpr int    SAMPLE_SIZE = 38;
inline constexpr double GAMMA       = 2.4;
inline constexpr double EPS_MIN     = 1e-12;
} // namespace detail

using arr_dbl3  = std::array<double, 3>;
using arr_int3  = std::array<int, 3>;
using arr_dbl38 = std::array<double, detail::SAMPLE_SIZE>;

namespace detail {
// Definition of constant spectral data and conversion matrices
// Contains arrays for various spectra (White, Cyan, Magenta, Yellow, Red, Green, Blue).
namespace BASE_SPECTRA {
inline constexpr arr_dbl38 W = {
    1.00116072718764,  1.00116065159728,  1.00116031922747,  1.00115867270789,  1.00115259844552,  1.00113252528998,
    1.00108500663327,  1.00099687889453,  1.00086525152274,  1.0006962900094,   1.00050496114888,  1.00030808187992,
    1.00011966602013,  0.999952765968407, 0.999821836899297, 0.999738609557593, 0.999709551639612, 0.999731930210627,
    0.999799436346195, 0.999900330316671, 1.00002040652611,  1.00014478793658,  1.00025997903412,  1.00035579697089,
    1.00042753780269,  1.00047623344888,  1.00050720967508,  1.00052519156373,  1.00053509606896,  1.00054022097482,
    1.00054272816784,  1.00054389569087,  1.00054448212151,  1.00054476959992,  1.00054489887762,  1.00054496254689,
    1.00054498927058,  1.000544996993};
inline constexpr arr_dbl38 C = {
    0.970585001322962,  0.970592498143425,  0.970625348729891,  0.970786806119017,  0.971368673228248,
    0.973163230621252,  0.976740223158765,  0.981587605491377,  0.986280265652949,  0.989949147689134,
    0.99249270153842,   0.994145680405256,  0.995183975033212,  0.995756750110818,  0.99591281828671,
    0.995606157834528,  0.994597600961854,  0.99221571549237,   0.986236452783249,  0.967943337264541,
    0.891285004244943,  0.536202477862053,  0.154108119001878,  0.0574575093228929, 0.0315349873107007,
    0.0222633920086335, 0.0182022841492439, 0.016299055973264,  0.0153656239334613, 0.0149111568733976,
    0.0146954339898235, 0.0145964146717719, 0.0145470156699655, 0.0145228771899495, 0.0145120341118965,
    0.0145066940939832, 0.0145044507314479, 0.0145038009464639};
inline constexpr arr_dbl38 M = {
    0.990673557319988, 0.990671524961979,  0.990662582353421,  0.990618107644795,  0.99045148087871,
    0.989871081400204, 0.98828660875964,   0.984290692797504,  0.973934905625306,  0.941817838460145,
    0.817390326195156, 0.432472805065729,  0.13845397825887,   0.0537347216940033, 0.0292174996673231,
    0.021313651750859, 0.0201349530181136, 0.0241323096280662, 0.0372236145223627, 0.0760506552706601,
    0.205375471942399, 0.541268903460439,  0.815841685086486,  0.912817704123976,  0.946339830166962,
    0.959927696331991, 0.966260595230312,  0.969325970058424,  0.970854536721399,  0.971605066528128,
    0.971962769757392, 0.972127272274509,  0.972209417745812,  0.972249577678424,  0.972267621998742,
    0.97227650946215,  0.972280243306874,  0.97228132482656};
inline constexpr arr_dbl38 Y = {
    0.0210523371789306, 0.0210564627517414, 0.0210746178695038, 0.0211649058448753, 0.0215027957272504,
    0.0226738799041561, 0.0258235649693629, 0.0334879385639851, 0.0519069663740307, 0.100749014833473,
    0.239129899706847,  0.534804312272748,  0.79780757864303,   0.911449894067384,  0.953797963004507,
    0.971241615465429,  0.979303123807588,  0.983380119507575,  0.985461246567755,  0.986435046976605,
    0.986738250670141,  0.986617882445032,  0.986277776758643,  0.985860592444056,  0.98547492767621,
    0.985176934765558,  0.984971574014181,  0.984846303415712,  0.984775351811199,  0.984738066625265,
    0.984719648311765,  0.984711023391939,  0.984706683300676,  0.984704554393091,  0.98470359630937,
    0.984703124077552,  0.98470292561509,   0.984702868122795};
inline constexpr arr_dbl38 R = {
    0.0315605737777207, 0.0315520718330149, 0.0315148215513658, 0.0313318044982702, 0.0306729857725527,
    0.0286480476989607, 0.0246450407045709, 0.0192960753663651, 0.0142066612220556, 0.0102942608878609,
    0.0076191460521811, 0.005898041083542,  0.0048233247781713, 0.0042298748350633, 0.0040599171299341,
    0.0043533695594676, 0.0053434425970201, 0.0076917201010463, 0.0135969795736536, 0.0316975442661115,
    0.107861196355249,  0.463812603168704,  0.847055405272011,  0.943185409393918,  0.968862150696558,
    0.978030667473603,  0.982043643854306,  0.983923623718707,  0.984845484154382,  0.985294275814596,
    0.985507295219825,  0.985605071539837,  0.985653849933578,  0.985677685033883,  0.985688391806122,
    0.985693664690031,  0.985695879848205,  0.985696521463762};
inline constexpr arr_dbl38 G = {
    0.0095560747554212, 0.0095581580120851, 0.0095673245444588, 0.0096129126297349, 0.0097837090401843,
    0.010378622705871,  0.0120026452378567, 0.0160977721473922, 0.026706190223168,  0.0595555440185881,
    0.186039826532826,  0.570579820116159,  0.861467768400292,  0.945879089767658,  0.970465486474305,
    0.97841363028445,   0.979589031411224,  0.975533536908632,  0.962288755397813,  0.92312157451312,
    0.793434018943111,  0.459270135902429,  0.185574103666303,  0.0881774959955372, 0.05436302287667,
    0.0406288447060719, 0.034221520431697,  0.0311185790956966, 0.0295708898336134, 0.0288108739348928,
    0.0284486271324597, 0.0282820301724731, 0.0281988376490237, 0.0281581655342037, 0.0281398910216386,
    0.0281308901665811, 0.0281271086805816, 0.0281260133612096};
inline constexpr arr_dbl38 B = {
    0.979404752502014,  0.97940070684313,   0.979382903470261,  0.979294364945594,  0.97896301460857,
    0.977814466694043,  0.974724321133836,  0.967198482343973,  0.949079657530575,  0.900850128940977,
    0.76315044546224,   0.465922171649319,  0.201263280451005,  0.0877524413419623, 0.0457176793291679,
    0.0284706050521843, 0.020527176756985,  0.0165302792310211, 0.0145135107212858, 0.0136003508637687,
    0.0133604258769571, 0.013548894314568,  0.0139594356366992, 0.014443425575357,  0.0148854440621406,
    0.0152254296999746, 0.0154592848180209, 0.0156018026485961, 0.0156824871281936, 0.0157248764360615,
    0.0157458108784121, 0.0157556123350225, 0.0157605443964911, 0.0157629637515278, 0.0157640525629106,
    0.015764589232951,  0.0157648147772649, 0.0157648801149616};
} // namespace BASE_SPECTRA

// detail::CIE Color Matching Functions weighted by D65 Standard Illuminant
namespace CIE {
inline constexpr std::array<arr_dbl38, 3> CMF = {
    {{0.0000646919989576, 0.0002194098998132, 0.0011205743509343, 0.0037666134117111, 0.011880553603799,
      0.0232864424191771, 0.0345594181969747, 0.0372237901162006, 0.0324183761091486, 0.021233205609381,
      0.0104909907685421, 0.0032958375797931, 0.0005070351633801, 0.0009486742057141, 0.0062737180998318,
      0.0168646241897775, 0.028689649025981,  0.0426748124691731, 0.0562547481311377, 0.0694703972677158,
      0.0830531516998291, 0.0861260963002257, 0.0904661376847769, 0.0850038650591277, 0.0709066691074488,
      0.0506288916373645, 0.035473961885264,  0.0214682102597065, 0.0125164567619117, 0.0068045816390165,
      0.0034645657946526, 0.0014976097506959, 0.000769700480928,  0.0004073680581315, 0.0001690104031614,
      0.0000952245150365, 0.0000490309872958, 0.0000199961492222},
     {0.000001844289444,  0.0000062053235865, 0.0000310096046799, 0.0001047483849269, 0.0003536405299538,
      0.0009514714056444, 0.0022822631748318, 0.004207329043473,  0.0066887983719014, 0.0098883960193565,
      0.0152494514496311, 0.0214183109449723, 0.0334229301575068, 0.0513100134918512, 0.070402083939949,
      0.0878387072603517, 0.0942490536184085, 0.0979566702718931, 0.0941521856862608, 0.0867810237486753,
      0.0788565338632013, 0.0635267026203555, 0.05374141675682,   0.042646064357412,  0.0316173492792708,
      0.020885205921391,  0.0138601101360152, 0.0081026402038399, 0.004630102258803,  0.0024913800051319,
      0.0012593033677378, 0.000541646522168,  0.0002779528920067, 0.0001471080673854, 0.0000610327472927,
      0.0000343873229523, 0.0000177059860053, 0.000007220974913},
     {0.000305017147638,
      0.0010368066663574,
      0.0053131363323992,
      0.0179543925899536,
      0.0570775815345485,
      0.113651618936287,
      0.17335872618355,
      0.196206575558657,
      0.186082370706296,
      0.139950475383207,
      0.0891745294268649,
      0.0478962113517075,
      0.0281456253957952,
      0.0161376622950514,
      0.0077591019215214,
      0.0042961483736618,
      0.0020055092122156,
      0.0008614711098802,
      0.0003690387177652,
      0.0001914287288574,
      0.0001495555858975,
      0.0000923109285104,
      0.0000681349182337,
      0.0000288263655696,
      0.0000157671820553,
      0.0000039406041027,
      0.000001584012587,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0}}};
}

// Conversion matrices and constants for various color space transformations.
//
// see https://github.com/w3c/csswg-drafts/issues/5922}
// see https://github.com/color-js/color.js/blob/main/src/spaces/srgb-linear.js}
// see https://github.com/color-js/color.js/blob/main/src/spaces/oklab.js}
namespace CONVERSION {
inline constexpr std::array<arr_dbl3, 3> RGB_XYZ = {{{0.41239079926595934, 0.357584339383878, 0.1804807884018343},
                                                     {0.21263900587151027, 0.715168678767756, 0.07219231536073371},
                                                     {0.01933081871559182, 0.11919477979462598, 0.9505321522496607}}};
inline constexpr std::array<arr_dbl3, 3> XYZ_RGB = {{{3.2409699419045226, -1.537383177570094, -0.4986107602930034},
                                                     {-0.9692436362808796, 1.8759675015077202, 0.04155505740717559},
                                                     {0.05563007969699366, -0.20397695888897652, 1.0569715142428786}}};
inline constexpr std::array<arr_dbl3, 3> XYZ_LMS = {{{0.819022437996703, 0.3619062600528904, -0.1288737815209879},
                                                     {0.0329836539323885, 0.9292868615863434, 0.0361446663506424},
                                                     {0.0481771893596242, 0.2642395317527308, 0.6335478284694309}}};
inline constexpr std::array<arr_dbl3, 3> LMS_XYZ = {{{1.2268798758459243, -0.5578149944602171, 0.2813910456659647},
                                                     {-0.0405757452148008, 1.112286803280317, -0.0717110580655164},
                                                     {-0.0763729366746601, -0.4214933324022432, 1.5869240198367816}}};
inline constexpr std::array<arr_dbl3, 3> LMS_LAB = {{{0.210454268309314, 0.7936177747023054, -0.0040720430116193},
                                                     {1.9779985324311684, -2.4285922420485799, 0.450593709617411},
                                                     {0.0259040424655478, 0.7827717124575296, -0.8086757549230774}}};
inline constexpr std::array<arr_dbl3, 3> LAB_LMS = {{{1.0, 0.3963377773761749, 0.2158037573099136},
                                                     {1.0, -0.1055613458156586, -0.0638541728258133},
                                                     {1.0, -0.0894841775298119, -1.2914855480194092}}};
} // namespace CONVERSION

// Primary template: defaults to false
template <typename T>
struct is_std_array : std::false_type {};

// Partial specialization for std::array
template <typename T, std::size_t N>
struct is_std_array<std::array<T, N>> : std::true_type {};

// Helper variable template (C++14+)
template <typename T>
inline constexpr bool is_std_array_v = is_std_array<T>::value;


template <typename M, typename V>
requires std::ranges::range<M> && std::ranges::range<std::ranges::range_value_t<M>> && std::ranges::range<V> &&
         is_std_array_v<M>
inline constexpr auto mulMatVec(M const &matrix, V const &vec) -> std::array<double, std::tuple_size<M>{}> {
    std::size_t const                        rows = std::tuple_size<M>{};
    std::size_t const                        cols = vec.size();
    std::array<double, std::tuple_size<M>{}> out;
    for (std::size_t i = 0; i < rows; ++i) {
        out[i] = 0.0;
        for (std::size_t j = 0; j < cols; ++j) { out[i] += matrix[i][j] * vec[j]; }
    }
    return out;
}

template <typename M, typename V>
requires std::ranges::range<M> && std::ranges::range<std::ranges::range_value_t<M>> && std::ranges::range<V>
inline constexpr auto mulMatVec(M const &matrix, V const &vec) {
    std::size_t const   rows = matrix.size();
    std::size_t const   cols = vec.size();
    std::vector<double> out(rows, 0.0);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) { out[i] += matrix[i][j] * vec[j]; }
    }
    return out;
}

inline constexpr double clamp01(double x) {
    return (x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x));
}

inline constexpr double uncompand(double x) {
    return x > 0.04045 ? std::pow((x + 0.055) / 1.055, detail::GAMMA) : x / 12.92;
}

inline constexpr double compand(double x) {
    return x > 0.0031308 ? 1.055 * std::pow(x, 1.0 / detail::GAMMA) - 0.055 : x * 12.92;
}

inline constexpr bool inGamut(const arr_dbl3 &lRGB, double epsilon = 0.0) {
    for (double x : lRGB) {
        if (x < -epsilon || x > 1.0 + epsilon) { return false; }
    }
    return true;
}

inline constexpr double deltaEOK(const arr_dbl3 &OK1, const arr_dbl3 &OK2) {
    if (OK1.size() != 3 || OK2.size() != 3) { throw std::invalid_argument("deltaEOK expects vectors of length 3"); }
    double d0 = OK1[0] - OK2[0];
    double d1 = OK1[1] - OK2[1];
    double d2 = OK1[2] - OK2[2];
    return std::sqrt(d0 * d0 + d1 * d1 + d2 * d2);
}

// Computes the Kubelka–Munk absorption/scattering parameter KS for a given spectral reflectance R.
//
// In Kubelka–Munk theory, the KS function reflects the ratio that controls the conversion from spectral
// reflectance to an equivalent absorption/scattering coeffidetail::CIEnt. The formulation
// <code>(1 - R)² / (2 * R)</code> is a common approximation that assumes a diffusely scattering medium.
inline constexpr double KS_func(double R) {
    return (1.0 - R) * (1.0 - R) / (2.0 * R);
}

// Computes the Kubelka–Munk mixing coeffidetail::CIEnt KM from a given KS value.
//
// The KM function transforms the KS parameter into a measure that can be linearly mixed.
// This conversion is essential because the Kubelka–Munk model assumes that when pigments are
// mixed, the resulting reflectance is a function of the weighted combination of the pigment
// absorption and scattering properties. The formula used here:
//
// <pre>
// KM(KS) = 1 + KS - √(KS² + 2KS)
// </pre>
//
// provides the appropriate transformation for blending multiple pigment spectra.
inline constexpr double KM_func(double KS) {
    return 1.0 + KS - std::sqrt(std::pow(KS, 2) + (2.0 * KS));
}
} // namespace detail

// Converter is a uninstantiable 'static utility' class
class Converter {
public:
    static constexpr arr_dbl3 sRGB_to_lRGB(const arr_int3 &sRGB);
    static constexpr arr_int3 lRGB_to_sRGB(const arr_dbl3 &lRGB);
    static constexpr arr_dbl3 XYZ_to_lRGB(const arr_dbl3 &XYZ);
    static constexpr arr_dbl3 lRGB_to_XYZ(const arr_dbl3 &lRGB);
    static constexpr arr_dbl3 XYZ_to_OKLab(const arr_dbl3 &XYZ);
    static constexpr arr_dbl3 OKLab_to_XYZ(const arr_dbl3 &OKLab);
    static constexpr arr_dbl3 OKLab_to_OKLCh(const arr_dbl3 &OKLab);
    static constexpr arr_dbl3 OKLCh_to_OKLab(const arr_dbl3 &OKLCh);

    static constexpr arr_dbl38          parseSpectralReflectanceFromLRGB(const arr_dbl3 &lRGB);
    static constexpr std::array<int, 4> parseCssColor(const std::string &str);

    Converter()                             = delete; // no default constructor
    Converter(const Converter &)            = delete; // no copy
    Converter &operator=(const Converter &) = delete; // no assignment
    ~Converter()                            = delete; // no destructor
    void *operator new(std::size_t)         = delete; // no heap
    void *operator new[](std::size_t)       = delete; // no heap array
    void *operator new(std::size_t, void *) = delete; // blocked placement-new
};

class Color {
public:
    // Public data
    arr_int3  sRGB; // 0..255 (no alpha, ignoring alpha)
    arr_dbl3  lRGB; // linear RGB [0..1]
    arr_dbl3  XYZ;  // XYZ space
    arr_dbl38 R;    // spectral reflectance, size = SIZE

    // Tinting strength (default 1)
    double tintingStrength = 1.0;
    double _luminance      = 0.0;

private:
    // Lazy-computed caches
    mutable std::optional<arr_dbl3>  _OKLab;
    mutable std::optional<arr_dbl3>  _OKLCh;
    mutable std::optional<arr_dbl38> _KS;

public:
    // Constructors
    Color(const std::string &css) {
        auto c     = Converter::parseCssColor(css);
        sRGB       = {c[0], c[1], c[2]};
        lRGB       = Converter::sRGB_to_lRGB(sRGB);
        R          = Converter::parseSpectralReflectanceFromLRGB(lRGB);
        XYZ        = Converter::lRGB_to_XYZ(lRGB);
        _luminance = std::max(XYZ[1], detail::EPS_MIN);
    }
    Color(const arr_int3 &srgb_) {
        sRGB       = srgb_;
        lRGB       = Converter::sRGB_to_lRGB(sRGB);
        R          = Converter::parseSpectralReflectanceFromLRGB(lRGB);
        XYZ        = Converter::lRGB_to_XYZ(lRGB);
        _luminance = std::max(XYZ[1], detail::EPS_MIN);
    }
    Color(const arr_dbl3 &lrgb_) {
        sRGB       = Converter::lRGB_to_sRGB(lrgb_);
        lRGB       = lrgb_;
        R          = Converter::parseSpectralReflectanceFromLRGB(lRGB);
        XYZ        = Converter::lRGB_to_XYZ(lRGB);
        _luminance = std::max(XYZ[1], detail::EPS_MIN);
    }
    Color(const std::vector<double> &spectralR) {
        if ((int)spectralR.size() != detail::SAMPLE_SIZE) {
            throw std::invalid_argument("spectralR must have length SIZE");
        }
        for (size_t i = 0; auto &r_item : R) { r_item = spectralR[i++]; }
        XYZ        = detail::mulMatVec(detail::CIE::CMF, R);
        lRGB       = Converter::XYZ_to_lRGB(XYZ);
        sRGB       = Converter::lRGB_to_sRGB(lRGB);
        _luminance = std::max(XYZ[1], detail::EPS_MIN);
    }
    Color(const arr_dbl38 &spectralR) {
        if ((int)spectralR.size() != detail::SAMPLE_SIZE) {
            throw std::invalid_argument("spectralR must have length SIZE");
        }
        R          = spectralR;
        XYZ        = detail::mulMatVec(detail::CIE::CMF, R);
        lRGB       = Converter::XYZ_to_lRGB(XYZ);
        sRGB       = Converter::lRGB_to_sRGB(lRGB);
        _luminance = std::max(XYZ[1], detail::EPS_MIN);
    }

    // Getters
    const arr_dbl3 &OKLab() const {
        if (! _OKLab) { _OKLab = Converter::XYZ_to_OKLab(XYZ); }
        return *_OKLab;
    }
    const arr_dbl3 &OKLCh() const {
        if (! _OKLCh) { _OKLCh = Converter::OKLab_to_OKLCh(OKLab()); }
        return *_OKLCh;
    }
    const arr_dbl38 &KS() const {
        if (! _KS) {
            _KS           = arr_dbl38{};
            auto &_KS_ref = _KS.value();
            for (int i = 0; i < detail::SAMPLE_SIZE; ++i) { _KS_ref[i] = detail::KS_func(R[i]); }
        }
        return _KS.value();
    }

    double luminance() const { return _luminance; }

    // Methods
private:
    // gamutMap (binary search of chroma)
    static constexpr Color gamutMap(const Color &color, double jnd = 0.03, double eps = 0.0001);

public:
    // Mixes multiple colors using a model based on the Kubelka–Munk theory.
    //
    // This function implements a mixing algorithm that is inspired by the Kubelka–Munk theory,
    // which models how light interacts with diffusely scattering and absorbing layers (such as pigments or paints).
    // The approach is as follows:
    //
    // - For each wavelength band (with SIZE samples), compute a weighted average of the KS values.
    // - Weights are determined by the square of a factor that considers both the square-root of the color's
    //   luminance and its tinting strength multiplied by a user-specified factor.
    // - The resulting weighted KS average is then converted back using the KM function to obtain the
    //   mixed spectral reflectance.
    //
    // In effect, this method blends pigments based on their optical absorption and scattering properties,
    // providing a physically motivated approximation for pigment mixing as described by Kubelka–Munk.
    //
    // Each argument should be provided as an array of two elements: [Color, factor]. The factor determines
    // the influence of that particular color in the overall mix (simple weight)
    static constexpr Color mix(const std::vector<std::pair<Color, double>> &colors);

    // Generates a palette of colors transitioning between two colors.
    static constexpr std::vector<Color> palette(const Color &a, const Color &b, int sz);

    // Interpolates between multiple colors based on a parameter t.
    // Each additional argument should be an array with two elements: [Color, position].
    static constexpr Color gradient(double t, const std::vector<std::pair<Color, double>> &stops);

    constexpr Color       toGamut(const std::string &method = "map") const;
    constexpr std::string toString(const std::string &format = "hex", const std::string &method = "map") const;
};


// ###############################
// ### IMPLEMENTATION PART    ####
// ###############################

constexpr arr_dbl3 Converter::sRGB_to_lRGB(const arr_int3 &sRGB) {
    return {detail::uncompand(sRGB[0] / 255.0), detail::uncompand(sRGB[1] / 255.0), detail::uncompand(sRGB[2] / 255.0)};
}
constexpr arr_int3 Converter::lRGB_to_sRGB(const arr_dbl3 &lRGB) {
    return {int(std::round(detail::compand(lRGB[0]) * 255.0)), int(std::round(detail::compand(lRGB[1]) * 255.0)),
            int(std::round(detail::compand(lRGB[2]) * 255.0))};
}
constexpr arr_dbl3 Converter::XYZ_to_lRGB(const arr_dbl3 &XYZ) {
    return detail::mulMatVec(detail::CONVERSION::XYZ_RGB, XYZ);
}
constexpr arr_dbl3 Converter::lRGB_to_XYZ(const arr_dbl3 &lRGB) {
    return detail::mulMatVec(detail::CONVERSION::RGB_XYZ, lRGB);
}
constexpr arr_dbl3 Converter::XYZ_to_OKLab(const arr_dbl3 &XYZ) {
    // lms = mat * XYZ
    arr_dbl3 lms = detail::mulMatVec(detail::CONVERSION::XYZ_LMS, XYZ);
    for (double &v : lms) { v = std::cbrt(v); }
    arr_dbl3 lab = detail::mulMatVec(detail::CONVERSION::LMS_LAB, lms);
    return lab;
}
constexpr arr_dbl3 Converter::OKLab_to_XYZ(const arr_dbl3 &OKLab) {
    arr_dbl3 lms = detail::mulMatVec(detail::CONVERSION::LAB_LMS, OKLab);
    for (double &v : lms) { v = v * v * v; }
    arr_dbl3 XYZ = detail::mulMatVec(detail::CONVERSION::LMS_XYZ, lms);
    return XYZ;
}
constexpr arr_dbl3 Converter::OKLab_to_OKLCh(const arr_dbl3 &OKLab) {
    double h = std::atan2(OKLab[2], OKLab[1]) * 180.0 / M_PI;
    if (h < 0.0) { h += 360.0; }
    return {OKLab[0], std::sqrt(std::pow(OKLab[1], 2) + std::pow(OKLab[2], 2)), h};
}
constexpr arr_dbl3 Converter::OKLCh_to_OKLab(const arr_dbl3 &OKLCh) {
    return {OKLCh[0], OKLCh[1] * std::cos(OKLCh[2] * M_PI / 180.0), OKLCh[1] * std::sin(OKLCh[2] * M_PI / 180.0)};
}


constexpr arr_dbl38 Converter::parseSpectralReflectanceFromLRGB(const arr_dbl3 &lRGB) {
    // This is the inverse of lRGB_to_R
    // The code subtracts minimum w, c, m, y, r, g, b contributions using detail::BASE_SPECTRA.

    double const w   = std::min(std::min(lRGB[0], lRGB[1]), lRGB[2]);
    arr_dbl3     cpy = {lRGB[0] - w, lRGB[1] - w, lRGB[2] - w};

    double const c = std::min(cpy[1], cpy[2]);
    double const m = std::min(cpy[0], cpy[2]);
    double const y = std::min(cpy[0], cpy[1]);
    double const r = std::max(0.0, std::min(cpy[0] - cpy[2], cpy[0] - cpy[1]));
    double const g = std::max(0.0, std::min(cpy[1] - cpy[2], cpy[1] - cpy[0]));
    double const b = std::max(0.0, std::min(cpy[2] - cpy[1], cpy[2] - cpy[0]));

    arr_dbl38 R{};
    for (int i = 0; i < detail::SAMPLE_SIZE; ++i) {
        // ensure at least epsilon
        R[i] = std::max(detail::EPS_MIN, w * detail::BASE_SPECTRA::W[i] + c * detail::BASE_SPECTRA::C[i] +
                                             m * detail::BASE_SPECTRA::M[i] + y * detail::BASE_SPECTRA::Y[i] +
                                             r * detail::BASE_SPECTRA::R[i] + g * detail::BASE_SPECTRA::G[i] +
                                             b * detail::BASE_SPECTRA::B[i]);
    }
    return R;
}

// CSS color parsing (#hex or rgb(...))
constexpr std::array<int, 4> Converter::parseCssColor(const std::string &str) {
    if (str.empty()) { throw std::invalid_argument("Empty CSS color string"); }
    if (str[0] == '#') {
        std::string hex = str.substr(1);
        if (hex.size() == 3) {
            // e.g. "#abc" → "aabbcc"
            std::string h2;
            h2.reserve(6);
            for (char c : hex) {
                h2.push_back(c);
                h2.push_back(c);
            }
            hex = h2;
        }
        if (hex.size() != 6 && hex.size() != 8) { throw std::invalid_argument("Invalid hex color length"); }
        int r = std::stoi(hex.substr(0, 2), nullptr, 16);
        int g = std::stoi(hex.substr(2, 2), nullptr, 16);
        int b = std::stoi(hex.substr(4, 2), nullptr, 16);
        int a = 255;
        if (hex.size() == 8) { a = std::stoi(hex.substr(6, 2), nullptr, 16); }
        return {r, g, b, a};
    }
    else if (str.rfind("rgb", 0) == 0) {
        auto open  = str.find('(');
        auto close = str.find(')');
        if (open == std::string::npos || close == std::string::npos) {
            throw std::invalid_argument("Malformed rgb() string");
        }
        std::string              inner = str.substr(open + 1, close - open - 1);
        std::vector<std::string> parts;
        {
            std::istringstream iss(inner);
            std::string        part;
            while (std::getline(iss, part, ',')) { parts.push_back(part); }
        }
        if (parts.size() < 3) { throw std::invalid_argument("rgb() must have 3 components"); }
        int comps[4] = {0, 0, 0, 255};
        for (int i = 0; i < 3; ++i) {
            std::string s = parts[i];
            // trim whitespace
            while (! s.empty() && std::isspace((unsigned char)s.front())) { s.erase(s.begin()); }
            while (! s.empty() && std::isspace((unsigned char)s.back())) { s.pop_back(); }
            bool   isPercent = (s.back() == '%');
            double val       = std::stod(s);
            if (isPercent && i < 3) { comps[i] = int(std::round(val * 2.55)); }
            else { comps[i] = int(std::round(val)); }
        }
        return {comps[0], comps[1], comps[2], comps[3]};
    }
    else { throw std::invalid_argument("Unsupported CSS color format: " + str); }
}

// gamutMap (binary search of chroma)
constexpr Color Color::gamutMap(const Color &color, double jnd, double eps) {
    double L = color.OKLCh()[0];
    if (L >= 1.0) { return Color(arr_int3{255, 255, 255}); }
    if (L <= 0.0) { return Color(arr_int3{0, 0, 0}); }
    if (detail::inGamut(color.lRGB)) { return color; }

    double h    = color.OKLCh()[2];
    double lo   = 0.0;
    double hi   = color.OKLCh()[1];
    bool   loIn = true;

    // clipped = lRGB_to_OKLab( clamp each component )
    arr_dbl3 clippedOK = Converter::XYZ_to_OKLab(Converter::lRGB_to_XYZ(
        {detail::clamp01(color.lRGB[0]), detail::clamp01(color.lRGB[1]), detail::clamp01(color.lRGB[2])}));
    double   E         = detail::deltaEOK(clippedOK, Converter::XYZ_to_OKLab(Converter::lRGB_to_XYZ(color.lRGB)));
    if (E < jnd) {
        // short-circuit
        arr_dbl3 tmpXYZ  = Converter::OKLab_to_XYZ(clippedOK);
        arr_dbl3 tmplRGB = Converter::XYZ_to_lRGB(tmpXYZ);
        return Color(tmplRGB);
    }

    while (hi - lo > eps) {
        double mid = 0.5 * (lo + hi);
        auto   ok  = Converter::OKLCh_to_OKLab({L, mid, h});
        auto   XYZ = Converter::OKLab_to_XYZ(ok);
        auto   lr  = Converter::XYZ_to_lRGB(XYZ);

        if (loIn && detail::inGamut(lr)) { lo = mid; }
        else {
            auto clipped2OK = Converter::XYZ_to_OKLab(
                Converter::lRGB_to_XYZ({detail::clamp01(lr[0]), detail::clamp01(lr[1]), detail::clamp01(lr[2])}));
            double dist = detail::deltaEOK(clipped2OK, ok);
            if (dist < jnd) {
                if (jnd - dist < eps) { break; }
                else {
                    loIn = false;
                    lo   = mid;
                }
            }
            else { hi = mid; }
        }
    }

    auto ok  = Converter::OKLCh_to_OKLab({L, lo, h});
    auto xyz = Converter::OKLab_to_XYZ(ok);
    auto lr  = Converter::XYZ_to_lRGB(xyz);
    return Color(lr);
}

constexpr Color Color::mix(const std::vector<std::pair<Color, double>> &colors) {
    arr_dbl38 Rm{};
    for (int i = 0; i < detail::SAMPLE_SIZE; ++i) {
        double ksSum   = 0.0;
        double concSum = 0.0;
        for (auto const &[c, factor] : colors) {
            double const conc  = std::pow(factor, 2) * std::pow(c.tintingStrength, 2) * c.luminance();
            concSum           += conc;
            ksSum             += c.KS()[i] * conc;
        }
        Rm[i] = detail::KM_func(ksSum / concSum);
    }
    return Color(Rm);
}

// Generates a palette of colors transitioning between two colors.
constexpr std::vector<Color> Color::palette(const Color &a, const Color &b, int sz) {
    std::vector<Color> out;
    out.reserve(sz);
    for (int i = 0; i < sz; ++i) { out.push_back(mix({{a, double(sz - 1 - i)}, {b, double(i)}})); }
    return out;
}

// Interpolates between multiple colors based on a parameter t.
// Each additional argument should be an array with two elements: [Color, position].
constexpr Color Color::gradient(double t, const std::vector<std::pair<Color, double>> &stops) {
    std::pair<Color, double> a      = {Color(arr_int3{0, 0, 0}), 0.0};
    std::pair<Color, double> b      = a;
    bool                     foundA = false, foundB = false;
    for (auto const &pr : stops) {
        if (! foundA || pr.second > a.second) {
            if (pr.second <= t) {
                a      = pr;
                foundA = true;
            }
        }
        if (! foundB || pr.second < b.second) {
            if (pr.second >= t) {
                b      = pr;
                foundB = true;
            }
        }
    }
    if (! foundA) { return b.first; }
    if (! foundB) { return a.first; }
    if (a.second == b.second) { return a.first; }
    double f = (t - a.second) / (b.second - a.second);
    return mix({{a.first, 1.0 - f}, {b.first, f}});
}


constexpr Color Color::toGamut(const std::string &method) const {
    std::string m = method;
    std::transform(m.begin(), m.end(), m.begin(), ::tolower);
    if (m == "clip") {
        arr_int3 clipped;
        for (int i = 0; i < 3; ++i) { clipped[i] = int(std::round(std::clamp(sRGB[i], 0, 255))); }
        return Color(clipped);
    }
    else if (m == "map") { return gamutMap(*this); }
    else { throw std::invalid_argument("Unknown gamut mapping method: " + method); }
}

constexpr std::string Color::toString(const std::string &format, const std::string &method) const {
    std::vector<int> outRGB(3);
    if (! detail::inGamut(lRGB)) {
        std::string m = method;
        std::transform(m.begin(), m.end(), m.begin(), ::tolower);
        if (m == "clip") {
            for (int i = 0; i < 3; ++i) { outRGB[i] = int(std::round(std::clamp(sRGB[i], 0, 255))); }
        }
        else if (m == "map") {
            Color mapped = gamutMap(*this);
            outRGB       = {mapped.sRGB[0], mapped.sRGB[1], mapped.sRGB[2]};
        }
        else { throw std::invalid_argument("Unknown method in toString: " + method); }
    }
    else { outRGB = {sRGB[0], sRGB[1], sRGB[2]}; }

    if (format == "hex" || format == "HEX" || format == "Hex") {
        std::ostringstream oss;
        oss << "#";
        for (int c : outRGB) { oss << std::uppercase << std::hex << std::setw(2) << std::setfill('0') << (c & 0xFF); }
        return oss.str();
    }
    else if (format == "rgb" || format == "RGB") {
        std::ostringstream oss;
        oss << "rgb(" << outRGB[0] << "," << outRGB[1] << "," << outRGB[2] << ")";
        return oss.str();
    }
    else { throw std::invalid_argument("Unknown format in toString: " + format); }
}

} // namespace pigment
} // namespace incom

#ifndef INCOM_INCPIG_NAMESPACE_ALIAS
#define INCOM_INCPIG_NAMESPACE_ALIAS
namespace incpig = incom::pigment;
#endif