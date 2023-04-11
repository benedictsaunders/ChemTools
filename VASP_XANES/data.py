
from utils import *
class potentials:
    defaults = {
        "H": "H",
        "He": "He",
        "Li": "Li_sv",
        "B":"B",
        "C":"C",
        "N":"N",
        "O":"O",
        "F":"F",
        "Ne":"Ne",
        "Na":"Na_pv",
        "Mg":"Mg",
        "Al":"Al",
        "Si":"Si",
        "P":"P",
        "S":"S",
        "Cl":"Cl",
        "Ar":"Ar",
        "K":"K_sv",
        "Ca":"Ca_sv",
        "Sc":"Sc_sv",
        "Ti":"Ti_sv",
        "V":"V_sv",
        "Cr":"Cr_sv",
        "Mn":"Mn_pv",
        "Fe":"Fe",
        "Co":"Co",
        "Ni":"Ni",
        "Cu":"Cu",
        "Zn":"Zn",
        "Ga":"Ga_d",
        "Ge":"Ge_d",
        "As":"As",
        "Se":"Se",
        "Br":"Br",
        "Kr":"Kr",
        "Rb":"Rb_sv",
        "Sr":"Sr_sv",
        "Y":"Y_sv",
        "Zr":"Zr_sv",
        "Nb":"Nb_sv",
        "Mo":"Mo_sv",
        "Tc":"Tc_pv",
        "Ru":"Ru_pv",
        "Rh":"Rh_pv",
        "Pd":"Pd",
        "Ag":"Ag",
        "Cd":"Cd",
        "In":"In_d",
        "Sn":"Sn_d",
        "Sb":"Sb",
        "Te":"Te",
        "I":"I",
        "Xe":"Xe",
        "Cs":"Cs_sv",
        "Ba":"Ba_sv",
        "La":"La",
        "Ce":"Ce",
        "Pr":"Pr_3",
        "Nd":"Nd_3",
        "Pm":"Pm_3",
        "Sm":"Sm_3",
        "Eu":"Eu_2",
        "Gd":"Gd_3",
        "Tb":"Tb_3",
        "Dy":"Dy_3",
        "Ho":"Ho_3",
        "Er":"Er_3",
        "Tm":"Tm_3",
        "Yb":"Yb_2",
        "Lu":"Lu_3",
        "Hf":"Hf_pv",
        "Ta":"Ta_pv",
        "W":"W_sv",
        "Re":"Re",
        "Os":"Os",
        "Ir":"Ir",
        "Pt":"Pt",
        "Au":"Au",
        "Hg":"Hg",
        "Tl":"Tl_d",
        "Pb":"Pb_d",
        "Bi":"Bi_d",
        "Po":"Po_d",
        "At":"At",
        "Rn":"Rn",
        "Fr":"Fr_sv",
        "Ra":"Ra_sv",
        "Ac":"Ac",
        "Th":"Th",
        "Pa":"Pa",
        "U":"U",
        "Np":"Np",
        "Pu":"Pu",
        "Am":"Am",
        "Cm":"Cm",
    }

class XrayNotation:
    edge = AliasDict({
        "K":{
        "n":1, "l":0, "s":0.5, "j":0.5
        },
        "L1":{
        "n":2, "l":0, "s":0.5, "j":0.5
        },
        "L2":{
        "n":2, "l":1, "s":0.5, "j":0.5
        },
        "L3":{
        "n":2, "l":1, "s":0.5, "j":1.5
        }
    })
XrayNotation.edge.add_alias("K", "K1")
